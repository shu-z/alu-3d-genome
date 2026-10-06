

import os
from itertools import combinations
import gzip
import math
import re, ast
import sys

from Bio import SeqIO
from Bio import SeqIO, pairwise2
from Bio.Seq import Seq
from pybedtools import BedTool
import numpy as np 
import pandas as pd
import pysam

import akita_utils_forplotting as utils
import akita_utils_scoring as scoring_utils
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "utils"))
from io_helpers import read_fasta

# Paths - edit for your environment
SUPREMO_DIR = os.environ.get("SUPREMO_DIR", "/pollard/home/szhang20/akita_variant_scoring")
AKITA_DIR = os.environ.get("AKITA_DIR", "/pollard/home/szhang20/akita")
SLURM_DIR = os.environ.get("SLURM_DIR", "/pollard/home/szhang20/slurm")
PROJECT_DIR = os.environ.get("ALU_PROJECT_DIR", "/pollard/home/szhang20/alu")
DATA_DIR = os.environ.get("ALU_DATA_DIR", "/pollard/data/projects/shzhang")
POLLARD_DATA = os.environ.get("POLLARD_DATA", "/pollard/data")

sys.path.append(f'{SUPREMO_DIR}/')  
MB=1048576
hg38_fa=pysam.FastaFile(f'{SUPREMO_DIR}/data/hg38.fa')


def align_alu_consensus(alu_new, alu_consensus, global_local='global'):
    """    
    Defaults for both global and local:
    Gap open penalty: -0.5
    Gap extension penalty: -0.1
    Mismatch penalty: 0 (matches contribute +1, and mismatches contribute 0)

    returns:
        alignment score: alignment score between the two sequences, scaled by the length of the alu_new sequence 
        identity_percentage: percent of alignment that has perfect match 
    """

    if global_local=='global':
        alignments = pairwise2.align.globalxx(alu_new, alu_consensus, one_alignment_only=True)
        #print(alignments)
    
        #scale alignment score by length of new alu sequence 
        #alignment_score=alignments[0].score/len(alu_new)

        #actually, let's just get the score itself and we can scale both later 
        alignments=pairwise2.align.globalms(alu_new, alu_consensus,
                                    2, -1,        # match, mismatch
                                    -5, -0.5,     # gap open, gap extend
                                one_alignment_only=True,penalize_end_gaps=False)

        alignment_score=alignments[0].score

        # get percent identity
        matches = sum(a == b for a, b in zip(alignments[0].seqA, alignments[0].seqB) if a != '-' and b != '-')
        length = alignments[0].end - alignments[0].start
        percent_identity = matches / length * 100

        # How many gaps?
        gaps_in_alunew = alignments[0].seqA.count('-')
        gaps_in_alucons = alignments[0].seqB.count('-')
  
        return alignment_score, percent_identity, gaps_in_alunew, gaps_in_alucons

# Function to insert one region into another and return metrics
def insert_region(target_row, insert_row):
    try:
        chrom = target_row['CHROM']
        target_start = target_row['POS']
        target_end = target_row['END']
        target_mid = np.floor((target_start + target_end)/2)

        #for now, shift is 0
        shift=0
        
        alu_len=target_row['END']-target_row['POS']
        REF_half_left = math.ceil((MB-alu_len)/2) - shift # if the REF allele is odd, shift right
        REF_half_right = math.floor((MB-alu_len)/2) + shift

        REF_start = target_row['POS'] - REF_half_left
        REF_stop = REF_start + MB 
    
        # REF sequence of target
        REF_seq = hg38_fa.fetch(chrom, REF_start, REF_stop).upper()
        REF_pred = utils.vector_from_seq(REF_seq)

        #ALT sequence (DEL)
        to_add_left = math.ceil(alu_len/2)
        to_add_right = math.floor(alu_len/2) 

        ALT_start = REF_start - to_add_left
        ALT_stop = REF_stop + to_add_right
        DEL_seq_left = hg38_fa.fetch(chrom, ALT_start, target_start).upper()
        DEL_seq_right = hg38_fa.fetch(chrom, target_end, ALT_stop).upper()
        DEL_seq=DEL_seq_left + DEL_seq_right
        DEL_pred = utils.vector_from_seq(DEL_seq)

        # Sequence of insertion
        insert_chrom = insert_row['CHROM']
        insert_start = insert_row['POS']
        insert_end = insert_row['END']
        insert_strand = insert_row['strand']
        insert_alu_name=insert_row['repName']
        
        insert_seq = hg38_fa.fetch(insert_chrom, insert_start, insert_end).upper()
        if insert_strand == '-':
            insert_seq = str(Seq(insert_seq).reverse_complement())
        
        # Alignment score
        alignment_score, percent_identity, gaps_in_aluorig, gaps_in_aluinsert = align_alu_consensus(
            hg38_fa.fetch(chrom, target_start, target_end).upper(),
            insert_seq.upper(),
            global_local='global'
        )
        
        # Padding for MB
        len_insert = len(insert_seq)
        region_pad = (MB - len_insert)/2
        left_start = np.floor(target_start - region_pad)
        right_end = np.floor(target_end + region_pad)
        left_seq_INS = hg38_fa.fetch(chrom, int(left_start), int(target_start))
        right_seq_INS = hg38_fa.fetch(chrom, int(target_end), int(right_end))
        
        ALT_insert_seq = (left_seq_INS + insert_seq + right_seq_INS).upper()
        ALT_insert_pred = utils.vector_from_seq(ALT_insert_seq)
        
        # Scoring
        MSE_insert_DEL = scoring_utils.mse(DEL_pred, ALT_insert_pred)
        Spearman_insert_DEL = scoring_utils.spearman(DEL_pred, ALT_insert_pred)
        MSE_orig_DEL = scoring_utils.mse(DEL_pred, REF_pred)
        Spearman_orig_DEL = scoring_utils.spearman(DEL_pred, REF_pred)
        MSE_orig_insert = scoring_utils.mse(REF_pred, ALT_insert_pred)
        Spearman_orig_insert = scoring_utils.spearman(REF_pred, ALT_insert_pred)

        return {
            'alu_insert_repName':insert_alu_name,
            'alu_insert_alignscore': alignment_score,
            'align_percentidentity': percent_identity,
            'gaps_in_aluorig' : gaps_in_aluorig,
            'gaps_in_aluinsert' : gaps_in_aluinsert,
            'insert_DEL_mse': MSE_insert_DEL,
            'insert_DEL_corr': Spearman_insert_DEL,
            'orig_DEL_mse': MSE_orig_DEL,
            'orig_DEL_corr': Spearman_orig_DEL,
            'orig_insert_mse': MSE_orig_insert,
            'orig_insert_corr': Spearman_orig_insert
        }
        
    except Exception as e:
        print(f"Error inserting {chrom}:{insert_start}-{insert_end} into {chrom}:{target_start}-{target_end}", e)
        return None

##################### read in necessary files 
#alu_dfam=pd.read_csv(f'{DATA_DIR}/alus/data/dfam/alu_dfam_consensus.txt', sep='\t', index_col=0)

#res_dir=f'{DATA_DIR}/alus/results/20251125_aluhg38_all_scores/'
res_dir=f'{SLURM_DIR}/'
out_dir=f'{PROJECT_DIR}/results/20260218_intron_swaps/'

alu_annot_withtracks=pd.read_csv(f'{res_dir}20260109_aluhg38_all_annot_GC_mapp_alu_repeats_genetrackoverlap.txt', sep='\t', index_col=0)
no_na_alu=alu_annot_withtracks.dropna(subset='mse_mean')

p99 = no_na_alu["mse_mean"].quantile(0.99)
p50 = no_na_alu["mse_mean"].quantile(0.50)

gene_list=['ABR','AHNAK','ANKRD11','CFAP74','CRK',
           'DIAPH1','DNAH2','FAM193A','FRMD5','GNB1',
           'ITGAL','MECP2','NGEF','NXN','PITPNC1',
           'PTPRS','SPATS2','SPG7','SPIDR','SRPK2',
           'TAF7', 'TUBA1C','YWHAE','AC244517.5','ABCA17P','SMG1P5']

for gene_i in gene_list:

    print(f'running {gene_i}...')

#look at PIT
    gene_df=alu_annot_withtracks[alu_annot_withtracks['single_gene_name']==gene_i]
    low_score=gene_df[gene_df['mse_mean']<p50]
    high_score=gene_df[gene_df['mse_mean']>p99]

    rows_out = []

    # print('    inserting high into low')
    # for idx_low, low_row in low_score.iterrows():
    #     print(f"{idx_low}")
    #     for idx_high, high_row in high_score.iterrows():
    #         res = insert_region(low_row, high_row)
    #         if res:
    #             row_dict = {
    #                 **low_row[['CHROM','POS','END','SVTYPE','SVLEN','strand','repName','mse_mean','corr_mean']].to_dict(),
    #                 'alu_category': 'low',          # indicate original row category
    #                 'alu_insert_idx': idx_high,
    #                 **res
    #             }
    #             rows_out.append(row_dict)
    
    # print('    inserting low into high')
    # for idx_high, high_row in high_score.iterrows():
    #     print(f"Processing high_score row {idx_high}")
    #     for idx_low, low_row in low_score.iterrows():
    #         res = insert_region(high_row, low_row)
    #         if res:
    #             row_dict = {
    #                 **high_row[['CHROM','POS','END','SVTYPE','SVLEN','strand','repName','mse_mean','corr_mean']].to_dict(),
    #                 'alu_category': 'high',         # indicate original row category
    #                 'alu_insert_idx': idx_low,
    #                 **res
    #             }
    #             rows_out.append(row_dict)

    print('    inserting high into high')
    for idx_high, high_row in high_score.iterrows():
        print(f"{idx_high}")
        for idx_high_insert, high_row_insert in high_score.iterrows():
            res = insert_region(high_row, high_row_insert)
            if res:
                row_dict = {
                    **high_row[['CHROM','POS','END','SVTYPE','SVLEN','strand','repName','mse_mean','corr_mean']].to_dict(),
                    'alu_category': 'high',          # indicate original row category
                    'alu_insert_idx': idx_high_insert,
                    **res
                }
                rows_out.append(row_dict)
    
    print('    inserting low into low')
    for idx_low, low_row in low_score.iterrows():
        print(f"Processing low_score row {idx_low}")
        for idx_low_insert, low_row_insert in low_score.iterrows():
            res = insert_region(low_row, low_row_insert)
            if res:
                row_dict = {
                    **low_row[['CHROM','POS','END','SVTYPE','SVLEN','strand','repName','mse_mean','corr_mean']].to_dict(),
                    'alu_category': 'low',         # indicate original row category
                    'alu_insert_idx': idx_low_insert,
                    **res
                }
                rows_out.append(row_dict)
    
    # Convert to DataFrame
    cross_insert_df = pd.DataFrame(rows_out)
    cross_insert_df.to_csv(f'{out_dir}{gene_i}_highlowmatch_insert.txt', sep='\t', index=False)

