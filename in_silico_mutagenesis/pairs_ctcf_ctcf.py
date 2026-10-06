

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

##################################################################################
#run!


#look at CTCF pairs 
#add in ID column to use as name 


#pick a random 1000 CTCF sites 

#res_dir=f'{DATA_DIR}/alus/results/20251125_aluhg38_all_scores/'
#alu_top=pd.read_csv(f'{res_dir}20251210_alutop1_1000sample.txt', sep='\t', index_col=0)



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
    track_cols=['CHROM', 'POS', 'END', 'tf', 'rel_score', '-log10(pval)', 'strand']
    chrom_list=[f'chr{i}' for i in range(1, 23)] + ['chrX']
    CTCF_sites=pd.read_csv(f'{PROJECT_DIR}/data/MA0139.1.tsv', sep='\t', header=None, names=track_cols)
    CTCF_sites=CTCF_sites[CTCF_sites['CHROM'].isin(chrom_list)]
    CTCF_sites["CTCF_id"] = CTCF_sites[["CHROM", "POS", "END"]].astype(str).agg("_".join, axis=1)
    print(CTCF_sites.head())
    CTCF_to_test=CTCF_sites.sample(n=1000, random_state=729)
    res_dir=f'{SLURM_DIR}/alu_pair_dels/'
    results = []
    start_time=time.time()
    for idx, row in CTCF_to_test.iterrows():
        print(idx)

        try:
            chrom = row['CHROM']
            start = int(row['POS'])
            end   = int(row['END']) 
            CTCF1_name = row['CTCF_id']

            # CTCF1 midpoint
            CTCF1_mid = math.ceil((start + end) / 2)

            #REF_half_left = math.ceil((seq_length - REF_len)/2) - shift # if the REF allele is odd, shift right
            #REF_half_right = math.floor((seq_length - REF_len)/2) + shift

            # REF sequence
            REF_seq = fetch_delete_centered(chrom, CTCF1_mid, [], hg38_fa, MB)

            #REF_seq_orig=hg38_fa.fetch(chrom, row['REF_start'], row['REF_stop'])

            #assert(REF_seq==REF_seq_orig)

            # CTCF1 deletion
            CTCF1_DEL_seq = fetch_delete_centered(chrom, CTCF1_mid, [(start, end)], hg38_fa, MB, debug=False)

            # Predictions for REF and CTCF1
            REF_pred = utils.vector_from_seq(REF_seq)
            REF_pred_mat = utils.mat_from_vector(REF_pred)

            CTCF1_pred = utils.vector_from_seq(CTCF1_DEL_seq)
            CTCF1_pred_mat = utils.mat_from_vector(CTCF1_pred)

            MSE_CTCF1 = scoring_utils.mse(REF_pred_mat, CTCF1_pred_mat)
            Spearman_CTCF1 = scoring_utils.spearman(REF_pred_mat, CTCF1_pred_mat)

            #print(f'MSE_CTCF1: {MSE_CTCF1}')

            # Nearby CTCFs
            CTCF_near = 200000
            nearby_all = CTCF_sites[
                (CTCF_sites['CHROM'] == chrom) &
                (CTCF_sites['POS'] > start - CTCF_near) &
                (CTCF_sites['END']   < end + CTCF_near)
            ]

            #make sure to remove the center CTCF!!
            nearby_all = nearby_all[~((nearby_all["POS"] == start) & (nearby_all["END"] == end))]

            #print('len nearby', len(nearby_all))
            nearby = nearby_all.sample(n=min(10, len(nearby_all)), random_state=0)

            for idx2, elem2 in nearby.iterrows():
                try:
                    elem2_start = int(elem2['POS'])
                    elem2_end   = int(elem2['END'])
                    CTCF2_name   = elem2['CTCF_id']

                    # CTCF2 deletion
                    CTCF2_DEL_seq = fetch_delete_centered(chrom, CTCF1_mid, [(elem2_start, elem2_end)], hg38_fa, MB)

                    # Both deletions
                    CTCF1_CTCF2_DEL_seq = fetch_delete_centered(
                        chrom, CTCF1_mid, [(start, end), (elem2_start, elem2_end)], hg38_fa, MB
                    )

                    # Predictions
                    CTCF2_pred = utils.vector_from_seq(CTCF2_DEL_seq)
                    CTCF2_pred_mat = utils.mat_from_vector(CTCF2_pred)

                    CTCF1_CTCF2_pred = utils.vector_from_seq(CTCF1_CTCF2_DEL_seq)
                    CTCF1_CTCF2_pred_mat = utils.mat_from_vector(CTCF1_CTCF2_pred)

                    # Scoring
                    MSE_CTCF2 = scoring_utils.mse(REF_pred_mat, CTCF2_pred_mat)
                    Spearman_CTCF2 = scoring_utils.spearman(REF_pred_mat, CTCF2_pred_mat)

                    MSE_CTCF1_CTCF2 = scoring_utils.mse(REF_pred_mat, CTCF1_CTCF2_pred_mat)
                    Spearman_CTCF1_CTCF2 = scoring_utils.spearman(REF_pred_mat, CTCF1_CTCF2_pred_mat)

                    MSE_CTCF1_v_CTCF2 = scoring_utils.mse(CTCF1_pred_mat, CTCF2_pred_mat)
                    Spearman_CTCF1_v_CTCF2 = scoring_utils.spearman(CTCF1_pred_mat, CTCF2_pred_mat)

                    MSE_CTCF1_v_CTCF1_2 = scoring_utils.mse(CTCF1_pred_mat, CTCF1_CTCF2_pred_mat)
                    Spearman_CTCF1_v_CTCF1_2 = scoring_utils.spearman(CTCF1_pred_mat, CTCF1_CTCF2_pred_mat)
                    MSE_CTCF2_v_CTCF1_2 = scoring_utils.mse(CTCF2_pred_mat, CTCF1_CTCF2_pred_mat)
                    Spearman_CTCF2_v_CTCF1_2 = scoring_utils.spearman(CTCF2_pred_mat, CTCF1_CTCF2_pred_mat)

                    # Plotting
                    #map_list = [REF_pred_mat, CTCF1_pred_mat, CTCF2_pred_mat, CTCF1_CTCF2_pred_mat]
                    #query_seq = (chrom, start, end)
                    #scores_list = [MSE_CTCF1, MSE_CTCF2, MSE_CTCF1_CTCF2]
                    #figs = plot_maps(map_list, query_seq, scores_list, CTCF1_name, CTCF2_name, CTCF1_mid, (elem2_start, elem2_end))
                    #figures.append(figs)

                    # Store results
                    results.append({
                        'CTCF1_idx': idx,
                        'CTCF1_name': CTCF1_name,
                        'CTCF1_chrom': chrom,
                        'CTCF1_start': start,
                        'CTCF1_end': end,
                        'CTCF2_idx': idx2,
                        'CTCF2_name': CTCF2_name,
                        'CTCF2_start': elem2_start,
                        'CTCF2_end': elem2_end,
                        'MSE_CTCF1_DEL': MSE_CTCF1,
                        'MSE_CTCF2_DEL': MSE_CTCF2,
                        'MSE_CTCF1_CTCF2_DEL': MSE_CTCF1_CTCF2,
                        'MSE_CTCF1_v_CTCF2': MSE_CTCF1_v_CTCF2,
                        'MSE_CTCF1_v_CTCF1_2': MSE_CTCF1_v_CTCF1_2,
                        'MSE_CTCF2_v_CTCF1_2': MSE_CTCF2_v_CTCF1_2,
                        'CORR_CTCF1_DEL': Spearman_CTCF1,
                        'CORR_CTCF2_DEL': Spearman_CTCF2,
                        'CORR_CTCF1_CTCF2_DEL': Spearman_CTCF1_CTCF2,
                        'CORR_CTCF1_v_CTCF2': Spearman_CTCF1_v_CTCF2,
                        'CORR_CTCF1_v_CTCF1_2': Spearman_CTCF1_v_CTCF1_2,
                        'CORR_CTCF2_v_CTCF1_2': Spearman_CTCF2_v_CTCF1_2,
                        'n_CTCF_100kb': len(nearby_all)
                    })
                except:
                    print(f'error in CTCF2 idx {idx2}')

        except:
            print(f'error in {idx}')
    end_time=time.time()
    print('time', end_time-start_time)
    df_results = pd.DataFrame(results)
    df_results.to_csv(f'{res_dir}20251217_CTCF_1000sample_200kbwindow_DELpairs.txt', sep='\t')


if __name__ == '__main__':
    main()
