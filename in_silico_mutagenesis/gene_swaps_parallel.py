

import os
import math
import multiprocessing as mp
import sys

from Bio import pairwise2
from Bio.Seq import Seq
import numpy as np
import pandas as pd
import pysam

import akita_utils_forplotting as utils
import akita_utils_scoring as scoring_utils

# Paths - edit for your environment
SUPREMO_DIR = os.environ.get("SUPREMO_DIR", "/pollard/home/szhang20/akita_variant_scoring")
AKITA_DIR = os.environ.get("AKITA_DIR", "/pollard/home/szhang20/akita")
PROJECT_DIR = os.environ.get("ALU_PROJECT_DIR", "/pollard/home/szhang20/alu")

mp.set_start_method("spawn", force=True) 
sys.path.append(f'{SUPREMO_DIR}/')

MB = 1048576
HG38_PATH = f'{SUPREMO_DIR}/data/hg38.fa'

# Per-worker fasta handle (pysam.FastaFile is NOT fork-safe)
_hg38_fa = None

def init_worker():
    global _hg38_fa
    _hg38_fa = pysam.FastaFile(HG38_PATH)

def rc(seq):
    return str(Seq(seq).reverse_complement())

def align_alu_consensus(alu_new, alu_consensus):
    alignments = pairwise2.align.globalms(
        alu_new, alu_consensus,
        2, -1, -5, -0.5,
        one_alignment_only=True, penalize_end_gaps=False
    )
    alignment_score = alignments[0].score
    matches = sum(a == b for a, b in zip(alignments[0].seqA, alignments[0].seqB) if a != '-' and b != '-')
    length = alignments[0].end - alignments[0].start
    percent_identity = matches / length * 100
    gaps_in_alunew = alignments[0].seqA.count('-')
    gaps_in_alucons = alignments[0].seqB.count('-')
    return alignment_score, percent_identity, gaps_in_alunew, gaps_in_alucons

def insert_region(target_row, insert_row, compute_rc=False):
    global _hg38_fa
    try:
        chrom = target_row['CHROM']
        target_start = target_row['POS']
        target_end = target_row['END']
        target_strand=target_row['strand']
        
        shift = 0
        alu_len = target_end - target_start
        REF_half_left = math.ceil((MB - alu_len) / 2) - shift

        REF_start = target_start - REF_half_left
        REF_stop = REF_start + MB

        REF_seq = _hg38_fa.fetch(chrom, REF_start, REF_stop).upper()
        REF_pred = utils.vector_from_seq(REF_seq)

        to_add_left = math.ceil(alu_len / 2)
        to_add_right = math.floor(alu_len / 2)
        ALT_start = REF_start - to_add_left
        ALT_stop = REF_stop + to_add_right
        DEL_seq = (_hg38_fa.fetch(chrom, ALT_start, target_start).upper() +
                   _hg38_fa.fetch(chrom, target_end, ALT_stop).upper())
        DEL_pred = utils.vector_from_seq(DEL_seq)

        insert_chrom = insert_row['CHROM']
        insert_start = insert_row['POS']
        insert_end = insert_row['END']
        insert_strand = insert_row['strand']
        insert_alu_name = insert_row['repName']

        # make sure inserted alu in same orientation as original alu
        # also turn all to plus strand for alignment 
        insert_seq = _hg38_fa.fetch(insert_chrom, insert_start, insert_end).upper()
        target_seq = _hg38_fa.fetch(chrom, target_start, target_end).upper()

        # flip insert if strands differ, so it matches target orientation
        insert_seq_final = rc(insert_seq) if target_strand != insert_strand else insert_seq

        # align both sequences in '+' orientation
        target_seq_toalign = rc(target_seq) if target_strand == '-' else target_seq
        insert_seq_toalign = rc(insert_seq) if insert_strand == '-' else insert_seq

        alignment_score, percent_identity, gaps_in_aluorig, gaps_in_aluinsert = align_alu_consensus(
            target_seq_toalign.upper(), insert_seq_toalign.upper())

        len_insert = len(insert_seq_final)
        region_pad = (MB - len_insert) / 2
        left_start = int(np.floor(target_start - region_pad))
        right_end = int(np.floor(target_end + region_pad))

        left_seq_INS = _hg38_fa.fetch(chrom, left_start, target_start)
        right_seq_INS = _hg38_fa.fetch(chrom, target_end, right_end)

        ALT_insert_seq = (left_seq_INS + insert_seq_final + right_seq_INS).upper()
        ALT_insert_pred = utils.vector_from_seq(ALT_insert_seq)

        results = {
            'alu_insert_repName': insert_alu_name,
            'alu_insert_alignscore': alignment_score,
            'align_percentidentity': percent_identity,
            'gaps_in_aluorig': gaps_in_aluorig,
            'gaps_in_aluinsert': gaps_in_aluinsert,
            'insert_DEL_mse': scoring_utils.mse(DEL_pred, ALT_insert_pred),
            'insert_DEL_corr': scoring_utils.spearman(DEL_pred, ALT_insert_pred),
            'orig_DEL_mse': scoring_utils.mse(DEL_pred, REF_pred),
            'orig_DEL_corr': scoring_utils.spearman(DEL_pred, REF_pred),
            'orig_insert_mse': scoring_utils.mse(REF_pred, ALT_insert_pred),
            'orig_insert_corr': scoring_utils.spearman(REF_pred, ALT_insert_pred)
        }

        # Optional reverse-complement scoring
        if compute_rc:
            REF_seq_RC = rc(REF_seq)
            DEL_seq_RC = rc(DEL_seq)
            ALT_insert_seq_RC = rc(ALT_insert_seq)

            REF_pred_RC = utils.vector_from_seq(REF_seq_RC)
            DEL_pred_RC = utils.vector_from_seq(DEL_seq_RC)
            ALT_insert_pred_RC = utils.vector_from_seq(ALT_insert_seq_RC)

            results.update({
                'insert_DEL_mse_RC':  scoring_utils.mse(DEL_pred_RC, ALT_insert_pred_RC),
                'insert_DEL_corr_RC': scoring_utils.spearman(DEL_pred_RC, ALT_insert_pred_RC),
                'orig_DEL_mse_RC':    scoring_utils.mse(DEL_pred_RC, REF_pred_RC),
                'orig_DEL_corr_RC':   scoring_utils.spearman(DEL_pred_RC, REF_pred_RC),
                'orig_insert_mse_RC': scoring_utils.mse(REF_pred_RC, ALT_insert_pred_RC),
                'orig_insert_corr_RC':scoring_utils.spearman(REF_pred_RC, ALT_insert_pred_RC),
            })

        return results

    except Exception as e:
        print(f"Error: {e}")
        return None

def process_pair(args):
    target_idx, target_row, insert_idx, insert_row, orig_category, insert_catgory, compute_rc = args
    res = insert_region(target_row, insert_row, compute_rc=compute_rc)
    if res is None:
        return None
    keep_cols_target = ['CHROM', 'POS', 'END', 'strand','repName', 'mse_mean', 'corr_mean', 'orig_idx']
    keep_cols_insert=['mse_mean', 'orig_idx']
    row_dict = {col: target_row[col] for col in keep_cols_target}
    row_dict['alu_orig_category'] = orig_category
    row_dict['alu_insert_category'] = insert_catgory
    row_dict['alu_insert_orig_idx']=insert_row['orig_idx']
    row_dict['alu_insert_mse_mean']=insert_row['mse_mean']
    row_dict.update(res)
    return row_dict

def build_tasks(df_target, df_insert, orig_category, insert_catgory, compute_rc):
    tasks = []
    for t_idx, t_row in df_target.iterrows():
        for i_idx, i_row in df_insert.iterrows():
            tasks.append((t_idx, t_row.to_dict(), i_idx, i_row.to_dict(), orig_category, insert_catgory, compute_rc))
    return tasks

if __name__ == '__main__':
    res_dir = f'{PROJECT_DIR}/results/paper_results/'
    out_dir = f'{res_dir}20260511_gene_swaps/'
    COMPUTE_RC = True
    N_WORKERS = 1

    alu_annot_withtracks = pd.read_csv(
        f'{res_dir}20260316_merged_aluhg38_all_featureannot_NA_blacklist_filter_region_alu.txt',
        sep='\t', index_col=0
    )
    no_na_alu = alu_annot_withtracks.dropna(subset='mse_mean')
    p99 = no_na_alu["mse_mean"].quantile(0.99)
    p50 = no_na_alu["mse_mean"].quantile(0.50)
    print('p99, p50', p99, p50)

    # gene_list = ['ABR','AHNAK','ANKRD11','CFAP74','CRK',
    #              'DIAPH1','DNAH2','FAM193A','FRMD5','GNB1',
    #              'ITGAL','MECP2','NGEF','NXN','PITPNC1',
    #              'PTPRS','SPATS2','SPG7','SPIDR','SRPK2',
    #              'TAF7','TUBA1C','YWHAE','AC244517.5','ABCA17P','SMG1P5']

    gene_list = [
    'ABR', 'AC104532.1', 'AHNAK', 'ANKRD11', 'ATP1B2',
    'CAPN1', 'CFAP74', 'CPNE7', 'CRK', 'DIAPH1',
    'DNAH2', 'FAM193A', 'FRMD5', 'GNB1', 'HIST1H2BE',
    'ITGAL', 'MAP1B', 'MECP2', 'MRPS27', 'NAA38',
    'NGEF', 'NLRP8', 'PCDHA1', 'PCDHGA1', 'PCDHGA11',
    'PCDHGA2', 'PITPNC1', 'PTPRS', 'RASSF3', 'SART1',
    'SPATS2', 'SPG7', 'SPIDR', 'SRPK2', 'TAF7',
    'TMEM94', 'TP53', 'TRIM38', 'TUBA1C', 'WRAP53',
    'YWHAE', 'ZNF473', 'ZNF771'
]

    for gene_i in gene_list:
        print(f'running {gene_i}...')
        gene_df = alu_annot_withtracks[alu_annot_withtracks['single_gene_name'] == gene_i]
        low_score = gene_df[gene_df['mse_mean'] < p50]
        high_score = gene_df[gene_df['mse_mean'] > p99]

        tasks = []
        tasks += build_tasks(high_score, high_score, 'high', 'high', COMPUTE_RC)
        tasks += build_tasks(low_score, low_score, 'low', 'low', COMPUTE_RC)
        tasks += build_tasks(low_score, high_score, 'low', 'high', COMPUTE_RC)
        tasks += build_tasks(high_score, low_score, 'high', 'low', COMPUTE_RC)

        print(f'    {len(tasks)} tasks to run')

        with mp.Pool(processes=N_WORKERS, initializer=init_worker) as pool:
            results = pool.map(process_pair, tasks, chunksize=1)

        rows_out = [r for r in results if r is not None]
        cross_insert_df = pd.DataFrame(rows_out)
        cross_insert_df.to_csv(f'{out_dir}{gene_i}_geneswap.txt',sep='\t', index=False)
        print(f'    done, wrote {len(rows_out)} rows')