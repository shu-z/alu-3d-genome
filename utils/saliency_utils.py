

import os
from itertools import combinations
import gzip
import json
import math
import re, ast
import sys

from basenji import dataset, seqnn, dna_io, layers
from Bio import SeqIO
from Bio import SeqIO, pairwise2
from Bio.Seq import Seq
from matplotlib.backends.backend_pdf import PdfPages
from pybedtools import BedTool
import matplotlib.pyplot as plt
import numpy as np 
import pandas as pd
import pyBigWig
import pysam
import tensorflow as tf #to calculate gradients

import akita_utils_forplotting as utils
import akita_utils_scoring as scoring_utils

# Paths - edit for your environment
SUPREMO_DIR = os.environ.get("SUPREMO_DIR", "/pollard/home/szhang20/akita_variant_scoring")
AKITA_DIR = os.environ.get("AKITA_DIR", "/pollard/home/szhang20/akita")
DATA_DIR = os.environ.get("ALU_DATA_DIR", "/pollard/data/projects/shzhang")
POLLARD_DATA = os.environ.get("POLLARD_DATA", "/pollard/data")

sys.path.append(f'{SUPREMO_DIR}/')  

########################################################################

#need to load in seqnn model for later 
repo_path=f'{SUPREMO_DIR}/'
model_file  = f'{repo_path}/Akita_model/model_best.h5'
params_file = f'{repo_path}/Akita_model/params.json'

with open(params_file) as params_open:
    params = json.load(params_open)
    params_model = params['model']
    params_train = params['train']

params_model['augment_shift']=0

seq_length = params_model['seq_length']
target_length = params_model['target_length']

seqnn_model = seqnn.SeqNN(params_model)

hic_diags = 2
tlen = (target_length-hic_diags) * (target_length-hic_diags+1) // 2

bin_size = seq_length//target_length

seqnn_model.restore(model_file)
print('Akita successfully loaded')

hic_params = params['model']['head_hic']
cropping = hic_params[5]['cropping']
target_length_cropped = target_length - 2 * cropping

half_patch_size = round(seq_length/2)

########################################################################

def compute_nucleotide_grads(seq, model=seqnn_model):
    """
    Compute nucleotide-level gradients for a sequence using the SeqNN model.
    
    Args:
        seq: string, DNA sequence (~1Mb)
        model: SeqNN/Basenji model
    Returns:
        grads: [seq_len, 4] array of gradients w.r.t one-hot input
    """

    # Convert sequence to one-hot float32
    seq_1hot_np = dna_io.dna_1hot(seq).astype(np.float32)
    seq_1hot_np = np.expand_dims(seq_1hot_np, 0)  # add batch dimension

    # Keras Input
    seq_input = tf.keras.Input(shape=(seqnn_model.seq_length, 4), dtype=tf.float32)

    # Pass through internal Keras model
    pred_tensor = model.model(seq_input)

    # Temporary model for gradient calculation
    temp_model = tf.keras.Model(inputs=seq_input, outputs=pred_tensor)

    # Wrap input as variable
    seq_var = tf.Variable(seq_1hot_np, dtype=tf.float32)

    # Compute gradient
    with tf.GradientTape() as tape:
        tape.watch(seq_var)
        pred = temp_model(seq_var)
        loss = tf.reduce_sum(pred)  # scalar to get gradients

    grads = tape.gradient(loss, seq_var)  # shape [1, seq_len, 4]
    grads = grads[0].numpy()  # remove batch dimension

    return grads  # [seq_len, 4]

plt.rcParams.update({
    'xtick.labelsize': 14,
    'ytick.labelsize': 14,
    'axes.labelsize': 16
})

def bin_and_smooth_saliency(REF_saliency, ALT1_saliency, 
                            len_ALT1_insert, 
                            insert_chrom, insert_pos, insert_Aluname,
                            bin_size=2048, smooth_window=5):
    """
    Plots saliency after binning (spatial averaging) and then applying a rolling window mean.

    The original 'window' (2048) is used here as the 'bin_size'.
    'smooth_window' is the size of the rolling mean applied to the binned data.
    """
    
    def rolling_window_mean(arr, window):
        """Compute rolling window mean with same-length output (centered)."""
        if window >= len(arr):
            return np.full_like(arr, np.nanmean(arr))
        kernel = np.ones(window) / window
        return np.convolve(arr, kernel, mode="same")

    def bin_data(arr, bin_size):
        """Bins the array by taking the mean of non-overlapping chunks."""
        L = len(arr)
        L_trimmed = L - (L % bin_size)
        arr_trimmed = arr[:L_trimmed]
        
        # Reshape and take the mean along the new axis (the bins)
        # The new shape is (number_of_bins, bin_size)
        binned_arr = arr_trimmed.reshape(-1, bin_size).mean(axis=1)
        return binned_arr

    ###############################################
    #align and crop REF with ALTs 

    min_len = min(len(REF_saliency),
                  len(ALT1_saliency) - len_ALT1_insert,
                  len(ALT2_saliency) - len_ALT2_insert)
    
    half_len_crop = min_len // 2
    
    # REF Cropping
    ref_mid = len(REF_saliency) // 2
    start_idx_ref = ref_mid - half_len_crop
    end_idx_ref = start_idx_ref + min_len 
    ref_aligned = REF_saliency[start_idx_ref : end_idx_ref]

    # ALT1 Cropping (using correct start index)
    L_ALT1_non_insert = len(ALT1_saliency) - len_ALT1_insert
    alt1_insert_start_idx = L_ALT1_non_insert // 2 
    alt1_left = ALT1_saliency[:alt1_insert_start_idx]
    alt1_right_start = alt1_insert_start_idx + len_ALT1_insert
    alt1_right = ALT1_saliency[alt1_right_start:]
    alt1_aligned = np.concatenate([
        alt1_left[-half_len_crop:], 
        alt1_right[:min_len - half_len_crop] 
    ])

    ################################################
    #bin and smooth the data 

    # 1. Bin the aligned data
    ref_binned = bin_data(ref_aligned, bin_size)
    alt1_binned = bin_data(alt1_aligned, bin_size)

    # 2. Smooth the binned data
    ref_sm = rolling_window_mean(ref_binned, smooth_window)
    alt1_sm = rolling_window_mean(alt1_binned, smooth_window)
    
    # X-axis for the BINNED data (distance from insert site in bins)
    L_binned = len(ref_sm)
    half_len_binned = L_binned // 2
    
    # The X-axis represents the center of each bin, scaled by bin_size
    # We use min_len for centering to include the potential leftover bases in the last bin
    x_binned = np.arange(-half_len_binned, L_binned - half_len_binned) * bin_size 
    x_no_bin = np.arange(-min_len//2, min_len//2)  

    #############################################
    #add in CTCF chip
    #min_len so window is same size as cropped sequences 
    start = insert_pos - min_len // 2
    end   = insert_pos + min_len // 2
    print(f'{insert_chrom}_{start}_{end}')

    # --- Load ChIP-seq bigWig ---
    bw = pyBigWig.open(f"{DATA_DIR}/alus/data/ENCFF883SWT_CTCF_FF_hg19.bigWig")
    chip_vals = np.array(bw.values(insert_chrom, start, end+1))
    bw.close()

    # Replace NaNs with 0
    chip_vals = np.nan_to_num(chip_vals)

    # Match x-axis (no binning here, could bin/smooth if needed)
    x_chip = np.arange(-min_len//2, min_len//2)

    ######################################
    # plot 
    
    col_HIV='#B33E52'
    col_Alu='#78B33E'
    col_CTCF='#7ABECC'
    
    fig, axes = plt.subplots(3, 1, figsize=(20, 8), sharex=True)
    fig.suptitle(f"Saliency around {insert_chrom}:{insert_pos}", fontsize=24)

    axes[0].plot(x_no_bin, ref_aligned, color="gray")
    axes[0].set_ylabel("REF")

    axes[1].plot(x_no_bin, alt1_aligned, color=col_HIV)
    axes[1].set_ylabel(f"HIV insert")

    # add in CTCF chip track 
    axes[2].plot(x_chip, chip_vals, color=col_CTCF)
    axes[2].set_ylabel("CTCF ChIP-seq")
    axes[2].set_xlabel("Distance from insert site (bp)")

    # vertical line at insert site
    for ax in axes:
        ax.axvline(0, color="black", linestyle="--", alpha=0.7)
        ax.set_xlim(x_binned.min(), x_binned.max())

    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.show()
    
