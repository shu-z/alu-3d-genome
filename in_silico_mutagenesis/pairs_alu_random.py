

import os
import collections 
import io
import math 
import sys
import sys, os, psutil, time, re, random 

from astropy.convolution import convolve
from astropy.convolution import Gaussian2DKernel
from cooltools.lib.numutils import interpolate_bad_singletons, set_diag, interp_nan
from cooltools.lib.numutils import observed_over_expected, adaptive_coarsegrain
from pybedtools import BedTool
from scipy.stats import pearsonr
import basenji
import cooler
import h5py
import intervaltree
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyBigWig
import pysam 

import akita_utils_forplotting as utils
import akita_utils_scoring as scoring_utils
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "utils"))
from io_helpers import read_fasta, save_figures_to_pdf

# Paths - edit for your environment
SUPREMO_DIR = os.environ.get("SUPREMO_DIR", "/pollard/home/szhang20/akita_variant_scoring")
AKITA_DIR = os.environ.get("AKITA_DIR", "/pollard/home/szhang20/akita")
SLURM_DIR = os.environ.get("SLURM_DIR", "/pollard/home/szhang20/slurm")
PROJECT_DIR = os.environ.get("ALU_PROJECT_DIR", "/pollard/home/szhang20/alu")
DATA_DIR = os.environ.get("ALU_DATA_DIR", "/pollard/data/projects/shzhang")
POLLARD_DATA = os.environ.get("POLLARD_DATA", "/pollard/data")


#this is the one that used to be in ~/akita_variant_scoring/akdemir_collab/bin/

MB=1048576

# hic imports

#plotting imports 
# import cooltools.lib.plotting
# from matplotlib.colors import LogNorm
# from mpl_toolkits.axes_grid1 import make_axes_locatable
# from matplotlib.ticker import EngFormatter
# from matplotlib.backends.backend_pdf import PdfPages
# import seaborn as sns

#other stuff 

            

##################################################################################

def fetch_delete_centered(
    chrom,
    alu_mid,
    deletions,
    hg38_fa,
    MB,
    chrom_length=None,
    debug=False
):
    """
    Fetch sequence centered on alu_mid unless alu_mid is deleted (then center on the deletion block
    that contains alu_mid). If deletions remove bases inside the desired MB window, expand fetch
    bounds to compensate from left/right so that after deletions you can return MB bases centered
    on the chosen center in the edited sequence. Pad with 'N' only if chromosome edges prevent
    obtaining MB bases.

    Args:
        chrom (str): chromosome name (same convention as hg38_fa.fetch)
        alu_mid (int): genomic coordinate of chosen midpoint (0-based or 1-based? see note)
        deletions (list of (start,end)): list of deletion intervals in genome coords (0-based, end-exclusive)
        hg38_fa: an object supporting .fetch(chrom, start, end) with 0-based half-open coords
        MB (int): desired final length
        chrom_length (int, optional): chromosome length to avoid fetching beyond end; if None
            function will try to fetch anyway and pad if needed
        debug (bool): print debugging information if True

    Returns:
        str: sequence of length MB (with 'N' padding if necessary)
    """

    half = MB // 2

    # 1) Find block that contains alu_mid (if any)
    center_block = None
    for ds, de in deletions:
        if ds <= alu_mid < de:
            center_block = (ds, de)
            break

    # chosen center coordinate (we use integer midpoint of block if center is deleted)
    if center_block is not None:
        center_coord = (center_block[0] + center_block[1]) // 2
    else:
        center_coord = alu_mid

    # 2) ideal window [ideal_start, ideal_end) around center_coord (genomic coords)
    ideal_start = center_coord - half
    ideal_end   = ideal_start + MB

    # 3) compute how many deleted bases fall inside the ideal window, split left/right
    left_comp = 0
    right_comp = 0
    # Also compute total deleted bases strictly left of center within the fetch region later
    # For now compute deletions overlapping the ideal window
    for ds, de in deletions:
        ov_s = max(ds, ideal_start)
        ov_e = min(de, ideal_end)
        if ov_s < ov_e:
            deleted_inside = ov_e - ov_s
            # allocate to left/right by whether overlap lies before/after center_coord
            if ov_e <= center_coord:
                left_comp += deleted_inside
            elif ov_s >= center_coord:
                right_comp += deleted_inside
            else:
                # split across center
                left_comp += center_coord - ov_s
                right_comp += ov_e - center_coord

    # 4) Expand fetch region by the computed compensation amounts
    fetch_start = ideal_start - left_comp
    fetch_end   = ideal_end + right_comp

    # Ensure coords are integers
    fetch_start = int(max(0, fetch_start))
    fetch_end = int(fetch_end)

    # If chrom_length provided, cap fetch_end and adjust fetch_start if necessary
    if chrom_length is not None:
        if fetch_end > chrom_length:
            fetch_end = chrom_length
            # if we hit chromosome end we might want more left to attempt MB bases,
            # but we will handle padding after deletions
    # 5) Fetch sequence
    seq = hg38_fa.fetch(chrom, fetch_start, fetch_end)
    seq_list = list(seq)  # mutable for deletion

    # 6) Apply deletions relative to fetch_start
    # Also compute how many bases deleted strictly to the left of center_coord within the fetched region
    deleted_before_center = 0
    for ds, de in deletions:
        rel_s = max(0, ds - fetch_start)
        rel_e = min(len(seq_list), de - fetch_start)
        if rel_s < rel_e:
            # remove bases
            for i in range(rel_s, rel_e):
                seq_list[i] = ""
            # compute deleted-before-center contribution:
            # overlap between deletion and [fetch_start, center_coord)
            ov_left = max(0, min(de, center_coord) - max(ds, fetch_start))
            deleted_before_center += ov_left

    seq_after_del = "".join(seq_list)

    # 7) Determine center index after deletions (index in seq_after_del)
    center_idx_original = center_coord - fetch_start
    center_idx_after = center_idx_original - deleted_before_center

    # 8) Compute final slice indices in the edited sequence
    final_start = center_idx_after - half
    final_end   = final_start + MB

    # If final indices are out of bounds, pad accordingly
    left_pad = right_pad = 0
    if final_start < 0:
        left_pad = -final_start
        final_start = 0
        final_end = MB  # shift right to maintain MB length
    if final_end > len(seq_after_del):
        right_pad = final_end - len(seq_after_del)
        final_end = len(seq_after_del)
        # final_start = final_end - MB  # Not needed; we'll extract and pad

    seq_final = seq_after_del[final_start:final_end]
    if left_pad:
        print('padding')
        seq_final = "N" * left_pad + seq_final
    if right_pad:
        print('padding')
        seq_final = seq_final + "N" * right_pad

    # Final safety: if still not MB (shouldn't be), pad/truncate
    if len(seq_final) < MB:
        print('padding')
        seq_final = seq_final + "N" * (MB - len(seq_final))
    elif len(seq_final) > MB:
        seq_final = seq_final[:MB]

    #print(seq_final[524288-10:524288+10])

    if debug:
        print("DEBUG fetch_delete_centered_compensate")
        print("chrom", chrom)
        print("alu_mid", alu_mid, "center_block", center_block, "center_coord", center_coord)
        print("ideal_start,end", ideal_start, ideal_end)
        print("left_comp,right_comp", left_comp, right_comp)
        print("fetch_start,fetch_end", fetch_start, fetch_end)
        print("fetched_len", len(seq))
        print("deleted_before_center", deleted_before_center)
        print("center_idx_original", center_idx_original, "center_idx_after", center_idx_after)
        print("final_start,final_end (in edited seq)", final_start, final_end)
        print("len seq_after_del", len(seq_after_del))
        print("len seq_final", len(seq_final))
        print("----------------------")

    assert len(seq_final) == MB, f"Final seq length {len(seq_final)} != {MB}"
    return seq_final

def subtract_intervals(allowed, forbidden):
    """
    allowed: list of (start, end)
    forbidden: list of (start, end)
    returns list of allowed intervals with forbidden removed
    """
    result = []
    for a_start, a_end in allowed:
        current = [(a_start, a_end)]
        for f_start, f_end in forbidden:
            new_current = []
            for c_start, c_end in current:
                # no overlap
                if f_end <= c_start or f_start >= c_end:
                    new_current.append((c_start, c_end))
                else:
                    # left remainder
                    if f_start > c_start:
                        new_current.append((c_start, f_start))
                    # right remainder
                    if f_end < c_end:
                        new_current.append((f_end, c_end))
            current = new_current
        result.extend(current)
    return result

##################################################################################
#run!


#res_dir=f'{DATA_DIR}/alus/results/20251125_aluhg38_all_scores/'



#results


def main():
    sys.path.append(f"{SUPREMO_DIR}/")
    nt = ['A', 'T', 'C', 'G']
    pd.set_option('display.max_columns', None)
    ModelSeq = collections.namedtuple('ModelSeq', ['chr', 'start', 'end', 'label'])
    chrom_list=[f'chr{i}' for i in range(1, 23)] + ['chrX']
    alu_dir=f'{PROJECT_DIR}/'
    hg38_fa=pysam.FastaFile(f'{SUPREMO_DIR}/data/hg38.fa')
    Alu=pd.read_csv(f'{alu_dir}data/Alu_hg38_repeatmasker_all.txt', sep='\t', index_col=0)
    res_dir=f'{SLURM_DIR}/alu_pair_dels/'
    alu_top=pd.read_csv(f'{res_dir}20251210_alutop1_1000sample.txt', sep='\t', index_col=0)
    results = []
    start_time=time.time()
    for idx, row in alu_top.iterrows():
        print(idx)

        try:
            chrom = row['CHROM']
            start = int(row['POS'])
            end   = int(row['END']) 
            Alu1_name = row['repName']

            # Alu1 midpoint
            alu1_mid = math.ceil((start + end) / 2)

            #REF_half_left = math.ceil((seq_length - REF_len)/2) - shift # if the REF allele is odd, shift right
            #REF_half_right = math.floor((seq_length - REF_len)/2) + shift

            # REF sequence
            REF_seq = fetch_delete_centered(chrom, alu1_mid, [], hg38_fa, MB)

            REF_seq_orig=hg38_fa.fetch(chrom, row['REF_start'], row['REF_stop'])

            #assert(REF_seq==REF_seq_orig)

            # Alu1 deletion
            Alu1_DEL_seq = fetch_delete_centered(chrom, alu1_mid, [(start, end)], hg38_fa, MB, debug=False)

            # Predictions for REF and Alu1
            REF_pred = utils.vector_from_seq(REF_seq)
            REF_pred_mat = utils.mat_from_vector(REF_pred)

            Alu1_pred = utils.vector_from_seq(Alu1_DEL_seq)
            Alu1_pred_mat = utils.mat_from_vector(Alu1_pred)

            MSE_Alu1 = scoring_utils.mse(REF_pred_mat, Alu1_pred_mat)
            Spearman_Alu1 = scoring_utils.spearman(REF_pred_mat, Alu1_pred_mat)

            #print(f'MSE_Alu1: {MSE_Alu1}')

            # Nearby Alus (this includes the center Alu)
            region_near = 200000 #how far away to look 
            nearby_all = Alu[
                (Alu['genoName'] == chrom) &
                (Alu['genoStart'] > start - region_near) &
                (Alu['genoEnd']   < end + region_near)
            ]

            # now set up sampling for random 300bp nearby (or really, the alu length)
            sample_len = (end-start)
            region_start = alu1_mid - region_near
            region_end   = alu1_mid + region_near

            # Forbidden intervals: center Alu + nearby Alus
            forbidden = [(start, end)]
            for _, r in nearby_all.iterrows():
                forbidden.append((int(r['genoStart']), int(r['genoEnd'])))

            # get all non-overlapping intervals
            allowed = [(region_start, region_end)]
            allowed_clean = subtract_intervals(allowed, forbidden)

            # Keep only intervals large enough
            allowed_clean = [
                (s, e) for s, e in allowed_clean if (e - s) >= sample_len
            ]

            if len(allowed_clean) == 0:
                print("No valid non-Alu region found")
                continue

            #sample 10 random 300bp (or matching Alu length) region not overlapping Alus
            for idx2 in range(0,10):
                try:

                    # Sample a 300 bp window
                    interval = allowed_clean[np.random.randint(len(allowed_clean))]
                    sample_start = np.random.randint(interval[0], interval[1] - sample_len + 1)
                    sample_end   = sample_start + sample_len
                    Random_name=f'{chrom}_{sample_start}_{sample_end}'

                    # Random deletion
                    Random_DEL_seq = fetch_delete_centered(chrom, alu1_mid, [(sample_start, sample_end)], hg38_fa, MB)

                    # Both deletions
                    Alu1_Random_DEL_seq = fetch_delete_centered(
                        chrom, alu1_mid, [(start, end), (sample_start, sample_end)], hg38_fa, MB
                    )

                    # Predictions
                    Random_pred = utils.vector_from_seq(Random_DEL_seq)
                    Random_pred_mat = utils.mat_from_vector(Random_pred)

                    Alu1_Random_pred = utils.vector_from_seq(Alu1_Random_DEL_seq)
                    Alu1_Random_pred_mat = utils.mat_from_vector(Alu1_Random_pred)

                    # Scoring
                    MSE_Random = scoring_utils.mse(REF_pred_mat, Random_pred_mat)
                    Spearman_Random = scoring_utils.spearman(REF_pred_mat, Random_pred_mat)

                    MSE_Alu1_Random = scoring_utils.mse(REF_pred_mat, Alu1_Random_pred_mat)
                    Spearman_Alu1_Random = scoring_utils.spearman(REF_pred_mat, Alu1_Random_pred_mat)

                    MSE_Alu1_v_Random = scoring_utils.mse(Alu1_pred_mat, Random_pred_mat)
                    Spearman_Alu1_v_Random = scoring_utils.spearman(Alu1_pred_mat, Random_pred_mat)

                    MSE_Alu1_v_Alu1_2 = scoring_utils.mse(Alu1_pred_mat, Alu1_Random_pred_mat)
                    Spearman_Alu1_v_Alu1_2 = scoring_utils.spearman(Alu1_pred_mat, Alu1_Random_pred_mat)
                    MSE_Random_v_Alu1_2 = scoring_utils.mse(Random_pred_mat, Alu1_Random_pred_mat)
                    Spearman_Random_v_Alu1_2 = scoring_utils.spearman(Random_pred_mat, Alu1_Random_pred_mat)

                    # Plotting
                    #map_list = [REF_pred_mat, Alu1_pred_mat, Random_pred_mat, Alu1_Random_pred_mat]
                    #query_seq = (chrom, start, end)
                    #scores_list = [MSE_Alu1, MSE_Random, MSE_Alu1_Random]
                    #figs = plot_maps(map_list, query_seq, scores_list, Alu1_name, Random_name, alu1_mid, (elem2_start, elem2_end))
                    #figures.append(figs)

                    # Store results
                    results.append({
                        'Alu1_idx': idx,
                        'Alu1_name': Alu1_name,
                        'Alu1_chrom': chrom,
                        'Alu1_start': start,
                        'Alu1_end': end,
                        'Random_idx': idx2,
                        'Random_name': Random_name,
                        'Random_start': sample_start,
                        'Random_end': sample_end,
                        'MSE_Alu1_DEL': MSE_Alu1,
                        'MSE_Random_DEL': MSE_Random,
                        'MSE_Alu1_Random_DEL': MSE_Alu1_Random,
                        'MSE_Alu1_v_Random': MSE_Alu1_v_Random,
                        'MSE_Alu1_v_Alu1_Random': MSE_Alu1_v_Alu1_2,
                        'MSE_Random_v_Alu1_Random': MSE_Random_v_Alu1_2,
                        'CORR_Alu1_DEL': Spearman_Alu1,
                        'CORR_Random_DEL': Spearman_Random,
                        'CORR_Alu1_Random_DEL': Spearman_Alu1_Random,
                        'CORR_Alu1_v_Random': Spearman_Alu1_v_Random,
                        'CORR_Alu1_v_Alu1_Random': Spearman_Alu1_v_Alu1_2,
                        'CORR_Random_v_Alu1_Random': Spearman_Random_v_Alu1_2,
                        'n_Alu_100kb': len(nearby_all)
                    })
                except:
                    print(f'error in Random idx {idx2}')

        except:
            print(f'error in {idx}')
    end_time=time.time()
    print('time', end_time-start_time)
    df_results = pd.DataFrame(results)
    df_results.to_csv(f'{res_dir}20251217_alutop1_1000sample_200kbwindow_randombp_DELpairs.txt', sep='\t')


if __name__ == '__main__':
    main()
