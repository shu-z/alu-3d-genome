#!/usr/bin/env python
# coding: utf-8

import os
from collections import Counter
import json
import math
import os, io, gzip 

from basenji import dataset, dna_io, layers
from basenji import seqnn_gpu as seqnn
from Bio.Seq import Seq
from cooltools.lib.numutils import (
     observed_over_expected, adaptive_coarsegrain, 
     interpolate_bad_singletons, interp_nan, set_diag
)
from scipy.stats import spearmanr
from statsmodels.stats.multitest import fdrcorrection
import numpy as np
import pandas as pd
import pysam
import tensorflow as tf

# Paths - edit for your environment
SUPREMO_DIR = os.environ.get("SUPREMO_DIR", "/pollard/home/szhang20/akita_variant_scoring")
AKITA_DIR = os.environ.get("AKITA_DIR", "/pollard/home/szhang20/akita")
POLLARD_DATA = os.environ.get("POLLARD_DATA", "/pollard/data")

# # # # # # # # # # # # # # # # # # 

#THIS ONE HAS ALL FUNCTIONS IN UTILS
#IE ONES FROM SCORES, SEQUENCES

#######

# Load required packages

#import cooltools
#from cooltools.lib.plotting import *

#from scores import get_DS_from_vector, get_scores, get_scores_BND, get_scores_SV
#from sequences import get_sequences, get_sequences_BND

half_patch_size = 2**19
MB = 2*half_patch_size
pixel_size = 2048
bins = 448
nt = ['A', 'T', 'C', 'G']

#katie's home dir has reference files 
#home_dir = '/pollard/home/ketringjoni/'
home_dir=f'{AKITA_DIR}/run_akita/refs/'
base_path = f"{AKITA_DIR}/"

hg38_lengths = pd.read_table(f'{home_dir}chrom_lengths', header = None, names = ['CHROM', 'chrom_max'])
#hg38_lengths = pd.read_table(f'{POLLARD_DATA}/projects/kgjoni/genome/hg38/chrom_lengths', header = None, names = ['CHROM', 'chrom_max'])

#gene_annot = pd.read_csv(f'{home_dir}genome/hg38/gene_annot')
#currently not being used
#gene_exp = pd.read_table(f'{home_dir}Akita/CBTN_collab/43062258-25f8-436d-980f-1cf1c9d44408.kallisto.abundance.tsv')
#centromere_coords = pd.read_table(f'{home_dir}genome/hg38/centromere_coords', header = None, names = ['CHROM', 'centromere'])

centromere_coords = pd.read_table(f'{home_dir}centromere_coords').rename(columns = {'chrom':'CHROM','start':'centromere'})[['CHROM', 'centromere']]

#centromere_coords = pd.read_table(f'{POLLARD_DATA}/projects/kgjoni/genome/hg38/centromere_coords').rename(columns = {'chrom':'CHROM','start':'centromere'})[['CHROM', 'centromere']]

fasta_file = f'{AKITA_DIR}/run_akita/refs/hg38.fa'
fasta_open = pysam.Fastafile(fasta_file)

# Note: if you want to get sequence starting at POS you should use this command on POS-1 since the genome starts at 0 with this command but 1 with how positions are annotated

# # # # # # # # # # # # # # # # # # 
# # # # # # Load model # # # # # #

#import tensor2tensor

if tf.__version__[0] == '1':
    tf.compat.v1.enable_eager_execution()

#limit each process's GPU VRAM
vram_lim = 7 #in GB
gpus = tf.config.list_physical_devices('GPU')
if gpus:
    tf.config.set_logical_device_configuration(
        gpus[0],
        [tf.config.LogicalDeviceConfiguration(memory_limit=vram_lim*1024)]  
    )

print(f'memory limit per process is {vram_lim}gb!')
    
#so it doesn't take up all the cpus 
#os.environ["OMP_NUM_THREADS"] = "1"
#tf.config.threading.set_intra_op_parallelism_threads(2)
#tf.config.threading.set_inter_op_parallelism_threads(2)

#idk exactly what this does but prevents the hd5 read error for loading model_file
os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"

#model_path = '{}/orig_files'.format(base_path)
model_path=f'{SUPREMO_DIR}/Akita_model'

params_file = model_path+'/params.json'
model_file  = model_path+'/model_best.h5'
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
print('model successfully loaded')

hic_params = params['model']['head_hic']
cropping = hic_params[5]['cropping']
target_length_cropped = target_length - 2 * cropping

# # # # # # # # # # # # # # # # # # 
# # # # # File processing # # # # #

def read_vcf(path, is_gzip=True):
    """
    Read vcf files and relabel their columns
    """
    
    if is_gzip:
        with gzip.open(path, 'rt') as f:
            lines = [l for l in f if not l.startswith('##')]
        
    else: #assumes ends in .vcf
        with open(path, 'r') as f:
            lines = [l for l in f if not l.startswith('##')]

    return pd.read_csv(
        io.StringIO(''.join(lines)),
        dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
               'QUAL': str, 'FILTER': str, 'INFO': str},
        sep='\t'
    ).rename(columns={'#CHROM': 'CHROM'})

# # # # # # # # # # # # # # # # # # 
# # Upstream helper functions # # #

def get_variant_position(CHR, POS, var_len, half_left, half_right):

    # Define variant position with respect to chromosome start and end

    # Get last coordinate of chromosome
    chrom_max = int(hg38_lengths[hg38_lengths.CHROM == CHR[3:]]['chrom_max']) 
    
    # Get centromere coordinate
    centromere = int(centromere_coords[centromere_coords.CHROM == CHR]['centromere'])

    # If variant too close to the beginning of chromosome
    if POS - half_left <= 0: 
        var_position = "chrom_start"
        print("Warning: Variant with start {}:{} not centered; too close to start of chromosome".format(CHR, POS))

    # If variant too close to left end of centromere
    elif POS + var_len - 1 + half_right > centromere and POS - half_left < centromere: 
        
        if POS > centromere:
            var_position = "chrom_centro_right"
            print("Warning: Variant with start {}:{} not centered; too close to right end of centromere".format(CHR, POS))

        elif POS <= centromere:
            var_position = "chrom_centro_left"
            print("Warning: Variant with start {}:{} not centered; too close to left end of centromere".format(CHR, POS))

    # If variant too close to the end of chromosome
    elif POS + var_len - 1 + half_right > chrom_max: 
        var_position = "chrom_end"
        print("Warning: Variant not centered; too close to end of chromosome")

    else:
        var_position = "chrom_mid"

    return var_position

# # # # # # # # # # # # # # # # # # 
# # # # # Running Akita # # # # # # 

def vector_from_seq(seq):
    
    # Get predicted matrix from ~1MB sequence
    
    seq_1hot = dna_io.dna_1hot(seq) 

    ensemble_shifts=[0,-5,5]
    sequence = tf.keras.Input(shape=(seqnn_model.seq_length, 4), name='sequence')
    sequences = [sequence]

    if len(ensemble_shifts) > 1:
        sequences = layers.EnsembleShift(ensemble_shifts)(sequences)
    sequences

    pred_targets = seqnn_model.predict(np.expand_dims(seq_1hot,0))[0,:,0]    

    return pred_targets

def vector_from_seq_with_grads(seq):
    # One-hot encode DNA sequence -> numpy
    seq_1hot_np = dna_io.dna_1hot(seq).astype(np.float32)  # shape [seq_len, 4]

    # Wrap in tf.Variable so gradients can flow
    seq_1hot = tf.Variable(seq_1hot_np[np.newaxis, ...], dtype=tf.float32)

    with tf.GradientTape() as tape:
        # Make sure to watch the input
        tape.watch(seq_1hot)

        # Forward pass (call model directly, not predict)
        pred_targets = seqnn_model(seq_1hot)  # [1, bins, 1]
        pred_targets = pred_targets[0, :, 0]                  # drop batch/channel

        # Define scalar loss so gradients can flow back
        loss = tf.reduce_sum(pred_targets)

    grads = tape.gradient(loss, seq_1hot)

    # Convert to numpy for convenience
    return pred_targets, grads

# # # # # # # # # # # # # # # # # # 
# # # # Processing output # # # # #

def from_upper_triu(vector_repr, matrix_len, num_diags):
    z = np.zeros((matrix_len,matrix_len))
    triu_tup = np.triu_indices(matrix_len,num_diags)
    z[triu_tup] = vector_repr
    for i in range(-num_diags+1,num_diags):
        set_diag(z, np.nan, i)
    return z + z.T

def mat_from_vector(pred_targets):
    mat = from_upper_triu(pred_targets,target_length_cropped,2)
    #mat = interp_all_nans(mat) 
    return mat

def upper_left(pred_matrix):
    # take upper left quarter of the matrix
    upper_left_qt = pred_matrix[:int(bins/2),:int(bins/2)]
    
    return upper_left_qt
    
def lower_right(pred_matrix):
    # take lower right quarter of the matrix
    lower_right_qt = pred_matrix[int(bins/2):,int(bins/2):]
    
    return lower_right_qt

#for plotting 

def plot_two_mats(REF_pred, ALT_pred):
    
    plt.figure(figsize=(12,10))
    target_index = 0
    vmin=-2; vmax=2

    # plot pred
    plt.subplot(121) 
    im = plt.matshow(REF_pred, fignum=False, cmap= 'RdBu_r', vmax=vmax, vmin=vmin)
    plt.title('Reference Prediction')

    # plot target 
    plt.subplot(122) 
    im = plt.matshow(ALT_pred, fignum=False, cmap= 'RdBu_r', vmax=vmax, vmin=vmin)
    plt.title('Alternate Prediction')

    plt.tight_layout()
    plt.show()

 ################## SCORES 

# # # # # # # # # # # # # # # # # # 
# # Getting disruption scores # # #

def get_DS_from_vector(vector1, vector2):

    #print('vector 1 shape', str(vector1.shape))
    #print(vector1.shape) 
    # Get disruption score between wt and variant vectors
    
    # A vecotr of sum of squared differences for each row of matrix 
    sub_vec = [x - y for x,y in zip(vector1, vector2)]
    
    # Get number of bins that are not nan (not including ones that overlap deletion)
    non_nan_values = np.count_nonzero(np.invert(np.isnan(sub_vec)))
    
    MSE = np.nansum([x**2 for x in sub_vec])/non_nan_values

    return MSE

def get_disruption_track(mat1, mat2):

    # A vector of mean of squared differences for each row of matrix 
    
    sub_mat = mat1 - mat2

    disruption_list = []
    for i in range(0,len(sub_mat)):
        
        # Get non-nan values in ith row
        non_nan_values = np.count_nonzero(~np.isnan(sub_mat[i]))
        
        # if there are values in ith row, get mean of squared differences
        if non_nan_values > 0:
            row_sum = np.nansum([x**2 for x in sub_mat[i]])
            row_avg = row_sum/non_nan_values
            
        # if the whole ith row is nan, MSE is nan for ith row
        else:
            row_avg = np.nan
        
        # save mean of squared differences for ith row
        disruption_list.append(row_avg)

    disruption_list = np.array(disruption_list) 
    
    return disruption_list

def get_correlation_track(wt_pred, del_pred_masked):

    # Get vector with correlation score between wt and masked deletion matrix for each row
    
    spearman_list = []
    for i in range(0,len(wt_pred)):
        
        non_nan_values = np.count_nonzero(~np.isnan(del_pred_masked[i]))
        
        # if there are values in ith row, get mean of squared differences
        if non_nan_values > 0:
            spearman_results = spearmanr(wt_pred[i], del_pred_masked[i], axis=0, nan_policy='omit') # ignores nans
            spearman_list.append(spearman_results)
        else:
            spearman_list.append([np.nan, np.nan])

    # Get nan indexes
    nan_indexes = np.unique(np.argwhere(np.isnan(spearman_list))[:,0])

    # Remove nans and correct for multiple testing with FDR

    # Remove nans
    spearman_array = np.array(spearman_list) 
    spearman_array_nonan = np.delete(spearman_array, nan_indexes, axis = 0)

    # FDR correct pvalues
    pvals = spearman_array_nonan[:,1]
    spearman_array_nonan[:,1] = fdrcorrection(pvals, alpha=0.05, method='indep', is_sorted=False)[1]

    # Make non-significant correlations 0
    for i in range(0,len(spearman_array_nonan)):
        if spearman_array_nonan[i,1] >= 0.05:
            spearman_array_nonan[i,0] = 0

    corr = spearman_array_nonan[:,0]

    # Re-insert nans
    for i in nan_indexes:
        corr = np.insert(corr, i, np.nan, axis=0)

    return corr

#########################################

def get_scores(CHR, POS, REF, ALT, SVTYPE, shift):
    print('getting scores')
    
#     try: # get_sequence will give an error if N composition > 5%
    REF_seq, ALT_seq = get_sequences(CHR, POS, REF, ALT, SVTYPE, shift)
    
    #this is where akita is actually rn
    #output vector size is (99681,)
    REF_vector, ALT_vector = vector_from_seq(REF_seq), vector_from_seq(ALT_seq)
    #print(REF_vector.shape)
    
    REF_pred1=0
    ALT_pred1=0
    #shouldn't need to mask if inversion
    if (len(REF) > pixel_size/2 or len(ALT) > pixel_size/2) and SVTYPE != "INV":
        #print('REF or ALT is larger than pixel size/2')
        #print('pixel size', pixel_size/2)
        #print('len REF: ', len(REF))
        #print('len ALT: ', len(ALT))

        REF_pred1, ALT_pred1 = mat_from_vector(REF_vector), mat_from_vector(ALT_vector)

        # mask matrices as long as not inversion?
        
        REF_pred, ALT_pred = mask_matrices(REF, ALT, REF_pred1, ALT_pred1, SVTYPE)

        # mask vectors
        REF_vector = REF_pred[np.triu_indices(len(REF_pred), 2)]
        ALT_vector = ALT_pred[np.triu_indices(len(ALT_pred), 2)]

    correlation, corr_pval = spearmanr(REF_vector, ALT_vector, nan_policy='omit')

    if corr_pval >= 0.05:
        correlation = 1 # this used to be 0, change to 1 on 4/26/22

    disruption_score = get_DS_from_vector(REF_vector, ALT_vector)

    #turn final REF and ALT vectors back into matrices 
    REF_mat, ALT_mat = mat_from_vector(REF_vector), mat_from_vector(ALT_vector)

    return REF_mat, ALT_mat, REF_pred1, ALT_pred1
    #print('non bnd shape')
    #print(REF_mat.shape)
    
    #get disruption and correlation tracks 
#     disruption_track=get_disruption_track(REF_mat, ALT_mat)
#     correlation_track=get_correlation_track(REF_mat, ALT_mat)
#     return disruption_score, correlation, disruption_track, correlation_track 

    #return REF_mat, ALT_mat
    
        #THISI SONLY FOR THE SCHAMETIAC

def get_scores_BND(REF_pred_L, REF_pred_R, ALT_pred):
   
    #print('shapes in get_scores_BND')
    #print(REF_pred_L.shape)
    #print(REF_pred_R.shape)
    #print(ALT_pred.shape)
    
    #inputs are all (448, 448)

    # Get REF and ALT vectors, excluding diagonal 
    indexes = np.triu_indices(bins/2, 2)
    
    #shape of below is (24753,)
    REF_UL = upper_left(REF_pred_L)[indexes]
    REF_LR = lower_right(REF_pred_R)[indexes]
    ALT_UL = upper_left(ALT_pred)[indexes]
    ALT_LR = lower_right(ALT_pred)[indexes]

    #print('ref_ul shape', str(REF_UL.shape))
    #print('ref_lr shape', str(REF_LR.shape))
    
    seq_REF = np.append(REF_UL, REF_LR)
    #print('seq ref shape' , str(seq_REF.shape))
    #shape is (49506,)? difers than non BND SVs?
    seq_ALT = np.append(ALT_UL, ALT_LR)
    
    # Get disruption score 
    mse = get_DS_from_vector(seq_REF, seq_ALT)
    
    # Get spearman correlation
    correlation, corr_pval = spearmanr(seq_REF, seq_ALT)
    
    if corr_pval >= 0.05:
        correlation = 0

    #turn final REF and ALT vectors back into matrices 
    #THIS DOENS'T WOKR FOR BND
    #REF_mat, ALT_mat = mat_from_vector(seq_REF), mat_from_vector(seq_ALT)

    #instead, just fill in matrices of np.nan bc BND doesn't include the interchromosomal tiles
    REF_mat, ALT_mat=np.empty((bins, bins,)), np.empty((bins, bins,))
    REF_mat[:], ALT_mat[:] = np.nan, np.nan

    #fill in upper lefts 
    REF_mat[:int(bins/2),:int(bins/2)]=upper_left(REF_pred_L)
    ALT_mat[:int(bins/2),:int(bins/2)]=upper_left(ALT_pred)
    
    #fill in lower rights 
    REF_mat[int(bins/2):,int(bins/2):]=lower_right(REF_pred_R)
    ALT_mat[int(bins/2):,int(bins/2):]=lower_right(ALT_pred)

#     #get disruption and correlation tracks 
#     disruption_track=get_disruption_track(REF_mat, ALT_mat)
#     correlation_track=get_correlation_track(REF_mat, ALT_mat)

     #return mse, correlation, disruption_track, correlation_track

    return REF_mat, ALT_mat

def get_scores_SV(CHR, POS, ALT, END, SVTYPE, shift):
    print('getting SCORES')
    
    # Get new REF and ALT alleles

    if SVTYPE in ["DEL", "DUP", "INV"]:
        
        if SVTYPE == "DEL":

            # Get REF and ALT allele sequences first
            REF = fasta_open.fetch(CHR, POS - 1, END).upper()
            #ALT is first nucleotide of REF you doofus
            ALT = REF[0]

        elif SVTYPE == "DUP":

            # Insert duplicated sequence before POS
            ALT = fasta_open.fetch(CHR, POS - 1, END).upper()
            REF = ALT[0]

        elif SVTYPE == "INV":

            REF = fasta_open.fetch(CHR, POS - 1, END).upper()
            ALT = REF[0] + str(Seq(REF[1:]).reverse_complement())

        if len(REF) > 2*half_patch_size or len(ALT) > 2*half_patch_size:
            # Variant larger than prediction window
            mse, correlation = np.nan, np.nan
            
        else:
            #only for schematic
            REF_mat, ALT_mat, REF_mat_before, ALT_mat_before = get_scores(CHR, POS, REF, ALT, SVTYPE, shift)
            
            #REF_mat, ALT_mat = get_scores(CHR, POS, REF, ALT, SVTYPE, shift)

    elif SVTYPE == "BND":

        REF_seq_L, REF_seq_R, ALT_seq = get_sequences_BND(CHR, POS, ALT, shift)

        REF_pred_L, REF_pred_R, ALT_pred = [mat_from_vector(vector) for vector in \
                                            [vector_from_seq(seq) for seq in [REF_seq_L, REF_seq_R, ALT_seq]]]

        # Get disruption score and correlation for this variant
        #mse, correlation, mse_track, corr_track = get_scores_BND(REF_pred_L, REF_pred_R, ALT_pred)
        
        #REF_mat, ALT_mat = get_scores_BND(REF_pred_L, REF_pred_R, ALT_pred)
        
        #this only for schematic
        REF_mat, ALT_mat, REF_mat_before, ALT_mat_before = get_scores_BND(REF_pred_L, REF_pred_R, ALT_pred)

    # haven't addressed big insertions (<INS> only have beginning and end of inserted sequence - ignore)

    #return REF_mat, ALT_mat

    #this only for schematic
    return REF_mat, ALT_mat, REF_mat_before, ALT_mat_before 

#####################################

def get_bin(x):
    
    # Get the bin number based on bp number in a sequence (ex: 2500th bp is in the 2nd bin)
    
    x_bin = math.ceil(x/pixel_size) - 32
    
    return x_bin

def mask_matrices(REF, ALT, REF_pred, ALT_pred, SVTYPE):

    # Mask reference and alternate predicted matrices, based on the type of variant, when they are centered in the sequence
    
    #shu note - don't have this function, using SVTYPE for now
    #variant_type = get_variant_type(REF, ALT)
    print('masking matrices')
    variant_type=SVTYPE
    print(variant_type)

    # Insertions: Mask REF, add nans if necessary, and mirror nans to ALT
    #if variant_type in ["Insertion", "Ins_sub"]:
    
    #i think this applies to dup but i should double check
    if variant_type in ["INS"]:

    # start with just the middle of the chromosome

        # Adjust reference sequence

        REF_len = len(REF)

        REF_half_left = math.ceil((MB - REF_len)/2) # if the REF allele is odd, shift right
        REF_half_right = math.floor((MB - REF_len)/2)

        # change REF allele to nans
        var_start = get_bin(REF_half_left - 1)
        var_end = get_bin(REF_half_left - 1 + REF_len)

        REF_pred_masked = REF_pred

        REF_pred_masked[var_start:var_end + 1, :] = np.nan
        REF_pred_masked[:, var_start:var_end + 1] = np.nan

        # Get start and end of ALT allele
        ALT_len = len(ALT)
        
        ALT_half_left = math.ceil((MB - ALT_len)/2) 
        ALT_half_right = math.floor((MB - ALT_len)/2)

        var_start_ALT = get_bin(ALT_half_left - 1)
        var_end_ALT = get_bin(ALT_half_right - 1 + ALT_len)

        # If the ALT allele falls on more bins than the REF allele, adjust ALT allele 
            # (add nan(s) to var and remove outside bin(s))
            # Otherwise don't mask
        
        if var_end_ALT - var_start_ALT > var_end - var_start:

            # Insert the rest of the nans corresponding to the ALT allele
            to_add = (var_end_ALT - var_start_ALT) - (var_end - var_start)

            for j in range(var_start, var_start + to_add): # range only includes the first variable 
                REF_pred_masked = np.insert(REF_pred_masked, j, np.nan, axis = 0)
                REF_pred_masked = np.insert(REF_pred_masked, j, np.nan, axis = 1)

            # Chop off the outside of the REF matrix 
            to_remove = (len(REF_pred_masked) - 448)/2

            REF_pred_masked = REF_pred_masked[math.floor(to_remove) : -math.ceil(to_remove), math.floor(to_remove) : -math.ceil(to_remove)]
            # remove less on the left bc that's where you put one less part of the variant with odd number of bp

            assert(len(REF_pred_masked) == 448)

        # Adjust alternate sequence

        # make all nans in REF_pred also nan in ALT_pred

        # change all nan
        REF_pred_novalues = REF_pred_masked.copy()

        REF_pred_novalues[np.invert(np.isnan(REF_pred_novalues))] = 0

        ALT_pred_masked = ALT_pred + REF_pred_novalues

        assert(len(ALT_pred_masked) == 448)

    #masking is different for DUP
    #DUP is separate because DUP is centered between the two identical/duplicated sequences 
    
    #THIS PART NOT CHANGED YET!!!!!!!!!!
    elif variant_type in ["DUP"]:

        # Adjust reference sequence
        
        #REF_len should be a single nucleotide 
        REF_len = len(REF)
        print('REF len')
        print(REF_len)

        h_alt=len(ALT)/2

        #half left to shift left an extra half of ALT!!! because not centered!!!!!!!!!
        #same goe sfor half right
        REF_half_left = math.ceil((MB - REF_len-h_alt)/2) # if the REF allele is odd, shift right
        REF_half_right = math.floor((MB - REF_len+h_alt)/2)

        # change REF allele to nans
        var_start = get_bin(REF_half_left - 1)
        var_end = get_bin(REF_half_left - 1 + h_alt)
        
        print('ref indices')
        print(var_start)
        print(var_end)

        REF_pred_masked = REF_pred

        REF_pred_masked[var_start:var_end + 1, :] = np.nan
        REF_pred_masked[:, var_start:var_end + 1] = np.nan

        # Get start and end of ALT allele

        #need to change this for dup!! since not centered!!!!!
        ALT_len = len(ALT)
        
        ALT_half_left = math.ceil((MB - ALT_len-h_alt)/2) 
        ALT_half_right = math.floor((MB - ALT_len+h_alt)/2)

        var_start_ALT = get_bin(ALT_half_left - 1)
        var_end_ALT = get_bin(ALT_half_left - 1 + h_alt)

        print('alt indices')
        print(var_start_ALT)
        print(var_end_ALT)

        # If the ALT allele falls on more bins than the REF allele, adjust ALT allele 
            # (add nan(s) to var and remove outside bin(s))
            # Otherwise don't mask
        
        if var_end_ALT - var_start_ALT > var_end - var_start:

            # Insert the rest of the nans corresponding to the ALT allele
            to_add = (var_end_ALT - var_start_ALT) - (var_end - var_start)

            for j in range(var_start, var_start + to_add): # range only includes the first variable 
                REF_pred_masked = np.insert(REF_pred_masked, j, np.nan, axis = 0)
                REF_pred_masked = np.insert(REF_pred_masked, j, np.nan, axis = 1)

            # Chop off the outside of the REF matrix 
            to_remove = (len(REF_pred_masked) - 448)/2

            REF_pred_masked = REF_pred_masked[math.floor(to_remove) : -math.ceil(to_remove), math.floor(to_remove) : -math.ceil(to_remove)]
            # remove less on the left bc that's where you put one less part of the variant with odd number of bp

            assert(len(REF_pred_masked) == 448)

        # Adjust alternate sequence

        # make all nans in REF_pred also nan in ALT_pred

        # change all nan
        REF_pred_novalues = REF_pred_masked.copy()

        REF_pred_novalues[np.invert(np.isnan(REF_pred_novalues))] = 0

        ALT_pred_masked = ALT_pred + REF_pred_novalues

        assert(len(ALT_pred_masked) == 448)

    # Deletions: Mask ALT, add nans if necessary, and mirror nans to REF
    elif variant_type in ["DEL"]:
        
        ALT_len = len(ALT)

        ALT_half_left = math.ceil((MB - ALT_len)/2) # if the ALT allele is odd, shift right
        ALT_half_right = math.floor((MB - ALT_len)/2)

        # change ALT allele to nans
        var_start = get_bin(ALT_half_left - 1)
        var_end = get_bin(ALT_half_left - 1 + ALT_len)

        ALT_pred_masked = ALT_pred

        ALT_pred_masked[var_start:var_end + 1, :] = np.nan
        ALT_pred_masked[:, var_start:var_end + 1] = np.nan

        # Get start and end of ALT allele
        REF_len = len(REF)
        
        REF_half_left = math.ceil((MB - REF_len)/2) 
        REF_half_right = math.floor((MB - REF_len)/2)

        var_start_REF = get_bin(REF_half_left - 1)
        var_end_REF = get_bin(REF_half_right - 1 + REF_len)

        # If the REF allele falls on more bins than the ALT allele, adjust REF allele 
            # (add nan(s) to var and remove outside bin(s))
            # Otherwise don't mask
        
        if var_end_REF - var_start_REF > var_end - var_start:

            # Insert the rest of the nans corresponding to the REF allele
            to_add = (var_end_REF - var_start_REF) - (var_end - var_start)

            for j in range(var_start, var_start + to_add): # range only includes the first variable 
                ALT_pred_masked = np.insert(ALT_pred_masked, j, np.nan, axis = 0)
                ALT_pred_masked = np.insert(ALT_pred_masked, j, np.nan, axis = 1)

            # Chop off the outside of the ALT matrix 
            to_remove = (len(ALT_pred_masked) - 448)/2

            ALT_pred_masked = ALT_pred_masked[math.floor(to_remove) : -math.ceil(to_remove), math.floor(to_remove) : -math.ceil(to_remove)]
            # remove less on the left bc that's where you put one less part of the variant with odd number of bp

            assert(len(ALT_pred_masked) == 448)

        # Adjust Reference sequence

        # make all nans in ALT_pred also nan in REF_pred

        # change all nan
        ALT_pred_novalues = ALT_pred_masked.copy()

        ALT_pred_novalues[np.invert(np.isnan(ALT_pred_novalues))] = 0

        REF_pred_masked = REF_pred + ALT_pred_novalues

        assert(len(REF_pred_masked) == 448)

    # SNPs or MNPs: Mask REF and mirror nans to ALT
    elif variant_type in ['SNP', 'MNP']:
        
        # start with just the middle of the chromosome

        # Adjust reference sequence

        REF_len = len(REF)

        REF_half_left = math.ceil((MB - REF_len)/2) # if the REF allele is odd, shift right
        REF_half_right = math.floor((MB - REF_len)/2)

        # change REF allele to nans
        var_start = get_bin(REF_half_left - 1)
        var_end = get_bin(REF_half_left - 1 + REF_len)

        REF_pred_masked = REF_pred

        REF_pred_masked[var_start:var_end + 1, :] = np.nan
        REF_pred_masked[:, var_start:var_end + 1] = np.nan

        # Adjust alternate sequence

        # make all nans in REF_pred also nan in ALT_pred

        # change all nan
        REF_pred_novalues = REF_pred_masked.copy()

        REF_pred_novalues[np.invert(np.isnan(REF_pred_novalues))] = 0

        ALT_pred_masked = ALT_pred + REF_pred_novalues

        assert(len(ALT_pred_masked) == 448)

    return REF_pred_masked, ALT_pred_masked

############## SEQUENCES

# # # # # # # # # # # # # # # # # # 
# # # Generating sequence # # # # #

#get sequence but this one also takes in a none option for shift
#this mostly applies to deletions since we aren't working with SNVs 
def get_sequences(CHR, POS, REF, ALT, SVTYPE, shift):
  
    # Get reference and alternate sequence from REF and ALT allele using reference genome 
    
    # Get reference sequence

    REF_len = len(REF)
    
    if shift =='none':
        REF_half_left = math.ceil((MB - REF_len)/2) # if the REF allele is odd, shift right
        REF_half_right = math.floor((MB - REF_len)/2)

    if shift == 'left':
        REF_half_left = math.ceil((MB - REF_len)/2) + 1 # if the REF allele is odd, shift right
        REF_half_right = math.floor((MB - REF_len)/2) - 1

    if shift == 'right':
        REF_half_left = math.ceil((MB - REF_len)/2) - 1 # if the REF allele is odd, shift right
        REF_half_right = math.floor((MB - REF_len)/2) + 1

    # Annotate whether variant is close to beginning or end of chromosome
    var_position = get_variant_position(CHR, POS, REF_len, REF_half_left, REF_half_right)
    #print('var_position', var_position)

    # Get last coordinate of chromosome
    chrom_max = int(hg38_lengths[hg38_lengths.CHROM == CHR[3:]]['chrom_max']) 

    # Get centromere coordinate
    centromere = int(centromere_coords[centromere_coords.CHROM == CHR]['centromere'])

    # Get start and end of reference sequence
    
    #most var positions should be chrom_mid
    if var_position == "chrom_mid":
        REF_start = POS - REF_half_left
        REF_stop = REF_start + MB 
        
        #if DUP, make sure whole DUP is centered -- this means shifting 
        #the REF to the left by half of the length of the DUP (center is between two dup sequences)
        if SVTYPE=='DUP':
            REF_start = math.floor(POS - REF_half_left - len(ALT)/2)
            REF_stop = REF_start + MB 
            #print('ref start ', REF_start)
            
            if REF_start<0:
                #print(REF_start)
                var_position=='chrom_start'
                REF_start=0
                REF_stop=MB

    elif var_position == "chrom_start": 
        REF_start = 0
        REF_stop = MB

    elif var_position == "chrom_centro_left": 
        REF_start = centromere - MB
        REF_stop = centromere

    elif var_position == "chrom_centro_right": 
        REF_start = centromere
        REF_stop = centromere + MB

    elif var_position == "chrom_end": 
        REF_start = chrom_max - MB
        REF_stop = chrom_max

    # Adjust variant position and sequences to the left and right of variant to reflect augmentation
    # for variants at edges: 
        # when you can't shift left, shift right twice for the 'left' sequence
        # when you can't shift right, shift left twice for the 'right' sequence
    if shift == 'left':
        if var_position in ["chrom_start", "chrom_centro_right"]:
            REF_start = REF_start + 2
            REF_stop = REF_stop + 2
            POS = POS - 2
        elif var_position in ["chrom_end", "chrom_centro_left"]: 
            REF_start = REF_start - 1
            REF_stop = REF_stop - 1
            POS = POS + 1
        else:
            POS = POS + 1

    if shift == 'right':
        if var_position in ["chrom_start", "chrom_centro_right"]: 
            REF_start = REF_start + 1
            REF_stop = REF_stop + 1
            POS = POS - 1
        elif var_position in ["chrom_end", "chrom_centro_left"]: 
            REF_start = REF_start - 2
            REF_stop = REF_stop - 2
            POS = POS + 2
        else:
            POS = POS - 1

    # Get reference sequence
    REF_seq = fasta_open.fetch(CHR, REF_start, REF_stop).upper()

    # Warn if Ns are more than 5% of sequence
    if Counter(REF_seq)['N']/MB*100 > 5:
        print("Error: N composition greater than 5%")
        raise ValueError('N composition greater than 5%')

    #check if this still works with shifts
    # Make sure that reference sequence matches given REF

    if var_position == "chrom_mid":
        if SVTYPE=='DUP':
            #assertion needs to be different 
            pass
        
        else:
            assert(REF_seq[(REF_half_left - 1) : (REF_half_left - 1 + REF_len)] == REF)

    elif var_position == "chrom_start": 
        assert(REF_seq[(POS - 1) : (POS - 1 + REF_len)] == REF) 

    elif var_position == "chrom_centro_right": 
        POS_adj = POS - centromere
        assert(REF_seq[(POS_adj - 1) : (POS_adj - 1 + REF_len)] == REF) 

    elif var_position in ["chrom_end", "chrom_centro_left"]: 
        assert(REF_seq[-(REF_stop - POS + 1) : -(REF_stop - POS + 1 - REF_len)] == REF)

    assert(len(REF_seq) == MB)

    # For SNPs, MNPs, Insertions, or Ins_subs: 
    #and duplications too 
    if len(REF) <= len(ALT):
    #but the above is also true for IMPRECISE dels, invs, for example?
    #no, REF and ALT are sequences/nts, not the stuff in the variant df

        # Create alternate sequence: change REF sequence at position from REF to ALT

        ALT_seq = REF_seq

        if var_position == "chrom_mid":
            
            #need to treat DUP separately to center whole duplication 
            #before, was treating as inversion 
            if SVTYPE=='DUP':

                #the REF starts left by an extra len(ALT)/2
                h_alt=len(ALT)/2
                
                #if decimal, round down for left side -- math.floor, round up for right side -- math.ceiling           
                #this is to offset the rounding from earlier 
                
                #get indices of where to break REF
                left_alt_idx=math.floor((REF_half_left - 1) - h_alt)
                right_alt_idx=math.ceil((REF_half_left - 1 + REF_len) - h_alt)
                
                #get final ALT sequence, should be centered between the dup sequences
                ALT_seq = ALT_seq[:left_alt_idx] + ALT + ALT_seq[right_alt_idx:]
                
            else:
                #this assumes is insertion, etc
                ALT_seq = ALT_seq[:(REF_half_left - 1)] + ALT + ALT_seq[(REF_half_left - 1 + REF_len):]

        elif var_position == "chrom_start": 
            ALT_seq = ALT_seq[:(POS - 1)] + ALT + ALT_seq[(POS - 1 + REF_len):]

        elif var_position == "chrom_centro_right": 
            POS_adj = POS - centromere
            ALT_seq = ALT_seq[:(POS_adj - 1)] + ALT + ALT_seq[(POS_adj - 1 + REF_len):]

        elif var_position == "chrom_centro_left": 
            ALT_seq = ALT_seq[:-(centromere - POS + 1)] + ALT + ALT_seq[-(centromere - POS + 1 - REF_len):]
            
        elif var_position == "chrom_end":
            ALT_seq = ALT_seq[:-(chrom_max - POS + 1)] + ALT + ALT_seq[-(chrom_max - POS + 1 - REF_len):]

        # Chop off ends of alternate sequence if it's longer 
        if len(ALT_seq) > len(REF_seq):
            to_remove = (len(ALT_seq) - len(REF_seq))/2

            if to_remove == 0.5:
                ALT_seq = ALT_seq[1:]
            else:
                ALT_seq = ALT_seq[math.ceil(to_remove) : -math.floor(to_remove)]

    # For Deletions of Del_subs
    elif len(REF) > len(ALT):

        del_len = len(REF) - len(ALT)

        to_add_left = math.ceil(del_len/2)
        to_add_right = math.floor(del_len/2)

        # Get start and end of reference sequence

        if var_position == "chrom_mid":
            ALT_start = REF_start - to_add_left
            ALT_stop = REF_stop + to_add_right

        elif var_position == "chrom_start": 
            ALT_start = 0
            ALT_stop = MB + del_len

        elif var_position == "chrom_centro_left": 
            ALT_start = centromere - MB - del_len
            ALT_stop = centromere

        elif var_position == "chrom_centro_right": 
            ALT_start = centromere
            ALT_stop = centromere + MB + del_len

        elif var_position == "chrom_end": 
            ALT_start = chrom_max - MB - del_len
            ALT_stop = chrom_max

        # Adjust variant position and sequences to the left and right of variant to reflect augmentation
        # for variants at edges: 
            # when you can't shift left, shift right twice for the 'left' sequence
            # when you can't shift right, shift left twice for the 'right' sequence
        if shift == 'left':
            if var_position in ["chrom_start", "chrom_centro_right"]:
                ALT_start = ALT_start + 2
                ALT_stop = ALT_stop + 2
            elif var_position in ["chrom_end", "chrom_centro_left"]: 
                ALT_start = ALT_start - 1
                ALT_stop = ALT_stop - 1

        if shift == 'right':
            if var_position in ["chrom_start", "chrom_centro_right"]: 
                ALT_start = ALT_start + 1
                ALT_stop = ALT_stop + 1
            elif var_position in ["chrom_end", "chrom_centro_left"]: 
                ALT_start = ALT_start - 2
                ALT_stop = ALT_stop - 2

        # Get alternate sequence
        ALT_seq = fasta_open.fetch(CHR, ALT_start, ALT_stop).upper()

        # Make sure that alternate sequence matches REF at POS

        if var_position == "chrom_mid":
            assert(ALT_seq[(REF_half_left - 1 + to_add_left) : (REF_half_left - 1 + to_add_left + REF_len)] == REF)

        elif var_position == "chrom_start": 
            assert(ALT_seq[(POS - 1) : (POS - 1 + REF_len)] == REF)

        elif var_position == "chrom_centro_right": 
            POS_adj = POS - centromere
            assert(ALT_seq[(POS_adj - 1) : (POS_adj - 1 + REF_len)] == REF)

        elif var_position == "chrom_centro_left": 
            assert(ALT_seq[-(centromere - POS + 1) : -(centromere - POS - REF_len + 1)] == REF)
            
        elif var_position == "chrom_end": 
            assert(ALT_seq[-(chrom_max - POS + 1) : -(chrom_max - POS - REF_len + 1)] == REF)

        # Change alternate sequence to match ALT at POS

        if var_position == "chrom_mid":
            # [:N] does not include N but [N:] includes N
            ALT_seq = ALT_seq[:(REF_half_left - 1 + to_add_left)] + ALT + ALT_seq[(REF_half_left - 1 + to_add_left + REF_len):] 

        elif var_position == "chrom_start": 
            ALT_seq = ALT_seq[:(POS - 1)] + ALT + ALT_seq[(POS - 1 + REF_len):]

        elif var_position == "chrom_centro_right": 
            POS_adj = POS - centromere
            ALT_seq = ALT_seq[:(POS_adj - 1)] + ALT + ALT_seq[(POS_adj - 1 + REF_len):]

        elif var_position == "chrom_centro_left": 
            ALT_seq = ALT_seq[:-(centromere - POS + 1)] + ALT + ALT_seq[-(centromere - POS - REF_len + 1):]
            
        elif var_position == "chrom_end": 
            ALT_seq = ALT_seq[:-(chrom_max - POS + 1)] + ALT + ALT_seq[-(chrom_max - POS - REF_len + 1):]

    assert(len(ALT_seq) == MB)

    #print(len(REF_seq), len(ALT_seq))
    return REF_seq, ALT_seq

#altered to take in shift
#i think this will throw an error if out of bounds -- need to alter like script above 
def get_sequences_BND(CHR, POS, ALT, shift):

    if '[' in ALT:

        if ALT[0] in nt:

            # t[p[

            CHR2 = ALT.split(':')[0].split('[')[1]
            POS2 = int(ALT.split('[')[1].split(':')[1])
            ALT_t = ALT.split('[')[0]
            
            if shift=='none':

                ALT_left = fasta_open.fetch(CHR, POS - half_patch_size, POS).upper() # don't inlcude POS
                ALT_right = fasta_open.fetch(CHR2, POS2 + 1, POS2 + 1 + half_patch_size).upper() 

                REF_for_left = fasta_open.fetch(CHR, POS - half_patch_size, POS + half_patch_size).upper()
                REF_for_right = fasta_open.fetch(CHR2, POS2 - half_patch_size, POS2 + half_patch_size).upper() 
            
            if shift=='left':
                ALT_left = fasta_open.fetch(CHR, POS - half_patch_size-1, POS).upper() # don't inlcude POS
                ALT_right = fasta_open.fetch(CHR2, POS2 + 1, POS2 + 1 + half_patch_size-1).upper() 

                REF_for_left = fasta_open.fetch(CHR, POS - half_patch_size-1, POS + half_patch_size-1).upper()
                REF_for_right = fasta_open.fetch(CHR2, POS2 - half_patch_size-1, POS2 + half_patch_size-1).upper() 
                
            if shift=='right':
                ALT_left = fasta_open.fetch(CHR, POS - half_patch_size+1, POS).upper() # don't inlcude POS
                ALT_right = fasta_open.fetch(CHR2, POS2 + 1, POS2 + 1 + half_patch_size+1).upper() 

                REF_for_left = fasta_open.fetch(CHR, POS - half_patch_size+1, POS + half_patch_size+1).upper()
                REF_for_right = fasta_open.fetch(CHR2, POS2 - half_patch_size+1, POS2 + half_patch_size+1).upper() 

        if ALT[0] not in nt:

            #  [p[t

            CHR2 = ALT.split(':')[0].split('[')[1]
            POS2 = int(ALT.split('[')[1].split(':')[1])
            ALT_t = ALT.split('[')[2]

            if shift=='none':
                ALT_left_revcomp = fasta_open.fetch(CHR2, POS2 + 1, POS2 + 1 + half_patch_size).upper() # don't include POS2
                ALT_left = str(Seq(ALT_left_revcomp).reverse_complement())
                ALT_right = fasta_open.fetch(CHR, POS + 1, POS + 1 + half_patch_size).upper()

                REF_for_left_revcomp = fasta_open.fetch(CHR2, POS2 - half_patch_size, POS2 + half_patch_size).upper() 
                REF_for_left = str(Seq(REF_for_left_revcomp).reverse_complement())
                REF_for_right = fasta_open.fetch(CHR, POS - half_patch_size, POS + half_patch_size).upper() 
                
            if shift=='left':
                ALT_left_revcomp = fasta_open.fetch(CHR2, POS2 + 1, POS2 + 1 + half_patch_size+1).upper() # don't include POS2
                ALT_left = str(Seq(ALT_left_revcomp).reverse_complement())
                ALT_right = fasta_open.fetch(CHR, POS + 1, POS + 1 + half_patch_size-1).upper()

                REF_for_left_revcomp = fasta_open.fetch(CHR2, POS2 - half_patch_size+1, POS2 + half_patch_size+1).upper() 
                REF_for_left = str(Seq(REF_for_left_revcomp).reverse_complement())
                REF_for_right = fasta_open.fetch(CHR, POS - half_patch_size-1, POS + half_patch_size-1).upper() 

            if shift=='right':
                ALT_left_revcomp = fasta_open.fetch(CHR2, POS2 + 1, POS2 + 1 + half_patch_size-1).upper() # don't include POS2
                ALT_left = str(Seq(ALT_left_revcomp).reverse_complement())
                ALT_right = fasta_open.fetch(CHR, POS + 1, POS + 1 + half_patch_size+1).upper()

                REF_for_left_revcomp = fasta_open.fetch(CHR2, POS2 - half_patch_size-1, POS2 + half_patch_size-1).upper() 
                REF_for_left = str(Seq(REF_for_left_revcomp).reverse_complement())
                REF_for_right = fasta_open.fetch(CHR, POS - half_patch_size+1, POS + half_patch_size+1).upper() 

    if ']' in ALT:

        if ALT[0] in nt:

            # t]p]

            CHR2 = ALT.split(':')[0].split(']')[1]
            POS2 = int(ALT.split(']')[1].split(':')[1])
            ALT_t = ALT.split(']')[0]
            
            if shift=='none':

                ALT_left = fasta_open.fetch(CHR, POS - half_patch_size, POS).upper() # don't include POS
                ALT_right_revcomp = fasta_open.fetch(CHR2, POS2 - half_patch_size, POS2).upper() # don't include POS2
                ALT_right = str(Seq(ALT_right_revcomp).reverse_complement())

                REF_for_left = fasta_open.fetch(CHR, POS - half_patch_size, POS + half_patch_size).upper()
                REF_for_right_revcomp = fasta_open.fetch(CHR2, POS2 - half_patch_size, POS2 + half_patch_size).upper()
                REF_for_right = str(Seq(REF_for_right_revcomp).reverse_complement())
                
            if shift=='left':

                ALT_left = fasta_open.fetch(CHR, POS - half_patch_size-1, POS).upper() # don't include POS
                ALT_right_revcomp = fasta_open.fetch(CHR2, POS2 - half_patch_size+1, POS2).upper() # don't include POS2
                ALT_right = str(Seq(ALT_right_revcomp).reverse_complement())

                REF_for_left = fasta_open.fetch(CHR, POS - half_patch_size-1, POS + half_patch_size-1).upper()
                REF_for_right_revcomp = fasta_open.fetch(CHR2, POS2 - half_patch_size+1, POS2 + half_patch_size+1).upper()
                REF_for_right = str(Seq(REF_for_right_revcomp).reverse_complement())

            if shift=='right':

                ALT_left = fasta_open.fetch(CHR, POS - half_patch_size+1, POS).upper() # don't include POS
                ALT_right_revcomp = fasta_open.fetch(CHR2, POS2 - half_patch_size-1, POS2).upper() # don't include POS2
                ALT_right = str(Seq(ALT_right_revcomp).reverse_complement())

                REF_for_left = fasta_open.fetch(CHR, POS - half_patch_size+1, POS + half_patch_size+1).upper()
                REF_for_right_revcomp = fasta_open.fetch(CHR2, POS2 - half_patch_size-1, POS2 + half_patch_size-1).upper()
                REF_for_right = str(Seq(REF_for_right_revcomp).reverse_complement())

        if ALT[0] not in nt:

            # ]p]t

            CHR2 = ALT.split(':')[0].split(']')[1]
            POS2 = int(ALT.split(']')[1].split(':')[1])
            ALT_t = ALT.split(']')[2]
            
            if shift=='none':
                ALT_left = fasta_open.fetch(CHR2, POS2 - half_patch_size, POS2).upper() # don't include POS2
                ALT_right = fasta_open.fetch(CHR, POS + 1, POS + 1 + half_patch_size).upper()

                REF_for_left = fasta_open.fetch(CHR2, POS2 - half_patch_size, POS2 + half_patch_size).upper() 
                REF_for_right = fasta_open.fetch(CHR, POS - half_patch_size, POS + half_patch_size).upper() 
                
            if shift=='left':
                ALT_left = fasta_open.fetch(CHR2, POS2 - half_patch_size-1, POS2).upper() # don't include POS2
                ALT_right = fasta_open.fetch(CHR, POS + 1, POS + 1 + half_patch_size-1).upper()

                REF_for_left = fasta_open.fetch(CHR2, POS2 - half_patch_size-1, POS2 + half_patch_size-1).upper() 
                REF_for_right = fasta_open.fetch(CHR, POS - half_patch_size-1, POS + half_patch_size-1).upper() 
                
            if shift=='right':
                ALT_left = fasta_open.fetch(CHR2, POS2 - half_patch_size+1, POS2).upper() # don't include POS2
                ALT_right = fasta_open.fetch(CHR, POS + 1, POS + 1 + half_patch_size+1).upper()

                REF_for_left = fasta_open.fetch(CHR2, POS2 - half_patch_size+1, POS2 + half_patch_size+1).upper() 
                REF_for_right = fasta_open.fetch(CHR, POS - half_patch_size+1, POS + half_patch_size+1).upper() 

    ALT_seq = ALT_left + ALT_t + ALT_right

    # chop off the sides if longer than ~1MB
    if len(ALT_seq) > MB:
        to_remove = (len(ALT_seq) - MB)/2

        if to_remove == 0.5:
            ALT_seq = ALT_seq[1:]
        else:
            ALT_seq = ALT_seq[math.ceil(to_remove) : -math.floor(to_remove)]

    return REF_for_left, REF_for_right, ALT_seq

