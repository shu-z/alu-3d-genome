

import os
from functools import partial
from itertools import combinations
import gzip
import multiprocessing as mp
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
PROJECT_DIR = os.environ.get("ALU_PROJECT_DIR", "/pollard/home/szhang20/alu")
DATA_DIR = os.environ.get("ALU_DATA_DIR", "/pollard/data/projects/shzhang")
POLLARD_DATA = os.environ.get("POLLARD_DATA", "/pollard/data")

MB=1048576


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

##################### read in necessary files 

#alu_dfam=pd.read_csv(f'{DATA_DIR}/alus/data/dfam/alu_dfam_consensus.txt', sep='\t', index_col=0)

#read in previous split of high/low scoring alus 
#annotate for high/low scores 



def process_row(idx_row, hg38_fa_path, dfam_csv_path):
    """
    Single-row worker — opens its own file handles so it's fork-safe.
    Returns (idx, dict of {col: val}) or (idx, None) on error.
    """
    import pysam
    import pandas as pd
    import numpy as np
    from Bio.Seq import Seq
    import sys
    sys.path.append(f'{SUPREMO_DIR}/')
    import akita_utils_forplotting as utils
    import akita_utils_scoring as scoring_utils

    idx, row = idx_row
    hg38_fa  = pysam.FastaFile(hg38_fa_path)
    alu_dfam = pd.read_csv(dfam_csv_path, sep='\t', index_col=0)
    MB       = 1048576
    results  = {}

    try:
        chrom      = row['CHROM']
        alu_start  = row['POS']
        alu_end    = row['END']
        alu_strand = row['strand']
        alu_mid    = np.floor((alu_start + alu_end) / 2)
        alu_seq    = hg38_fa.fetch(chrom, int(alu_start), int(alu_end))
        revcomp    = (alu_strand == '-')

        REF_seq  = hg38_fa.fetch(chrom, int(np.floor(alu_mid - MB/2)), int(np.floor(alu_mid + MB/2))).upper()
        REF_pred = utils.vector_from_seq(REF_seq)

        REF_seq_RC=str(Seq(REF_seq).reverse_complement())
        REF_pred_RC=utils.vector_from_seq(REF_seq_RC)

        #also compare to deletion of the alu 
        left_seq_DEL = hg38_fa.fetch(chrom, int(np.floor(alu_start - MB/2)), alu_start)
        right_seq_DEL = hg38_fa.fetch(chrom, alu_end, int(np.floor(alu_end + MB/2)))

        DEL_seq=(left_seq_DEL+right_seq_DEL).upper()
        DEL_pred = utils.vector_from_seq(DEL_seq)

        DEL_seq_RC=str(Seq(DEL_seq).reverse_complement())
        DEL_pred_RC=utils.vector_from_seq(DEL_seq_RC)

        for _, dfam_row in alu_dfam.iterrows():
            try:
                alu_name     = dfam_row['AluName']
                consensus_seq = dfam_row['Consensus_seq']

                if revcomp:
                    #means on + strand, Alu is reversed. so we reverse. 
                    consensus_seq = str(Seq(consensus_seq).reverse_complement())

                alignment_score, percent_identity, gaps_in_alunew, gaps_in_alucons  = align_alu_consensus(
                                                        alu_seq.upper(), consensus_seq.upper(), global_local='global')
                len_alu      = len(consensus_seq)
                region_pad   = (MB - len_alu) / 2
                left_seq     = hg38_fa.fetch(chrom, int(np.floor(alu_start - region_pad)), int(alu_start))
                right_seq    = hg38_fa.fetch(chrom, int(alu_end), int(np.floor(alu_end + region_pad)))
                ALT_seq      = (left_seq + consensus_seq + right_seq).upper()
                ALT_pred     = utils.vector_from_seq(ALT_seq)

                ALT_seq_RC=str(Seq(ALT_seq).reverse_complement())
                ALT_pred_RC=utils.vector_from_seq(ALT_seq_RC)

                results[f'mse_{alu_name}']        = scoring_utils.mse(REF_pred, ALT_pred)
                results[f'corr_{alu_name}']       = scoring_utils.spearman(REF_pred, ALT_pred)
                results[f'mse_RC_{alu_name}']        = scoring_utils.mse(REF_pred_RC, ALT_pred_RC)
                results[f'corr_RC_{alu_name}']       = scoring_utils.spearman(REF_pred_RC, ALT_pred_RC)

                results[f'mse_DEL_{alu_name}']        = scoring_utils.mse(DEL_pred, ALT_pred)
                results[f'corr_DEL_{alu_name}']       = scoring_utils.spearman(DEL_pred, ALT_pred)
                results[f'mse_DEL_RC_{alu_name}']        = scoring_utils.mse(DEL_pred_RC, ALT_pred_RC)
                results[f'corr__DEL_RC_{alu_name}']       = scoring_utils.spearman(DEL_pred_RC, ALT_pred_RC)

                results[f'alignscore_{alu_name}']      = alignment_score
                results[f'percent_identity_{alu_name}']    = percent_identity
                results[f'gaps_alunew_{alu_name}']  = gaps_in_alunew
                results[f'gaps_alucons_{alu_name}'] = gaps_in_alucons

            except Exception as e:
                print(f'  [idx {idx}] ERROR on {dfam_row["AluName"]}: {e}')

    except Exception as e:
        print(f'[idx {idx}] ERROR: {e}')
        return idx, None

    print(f'[idx {idx}] done — {len(results)} scores')
    return idx, results

N_WORKERS = 9

if __name__ == '__main__':
    # spawn so TF/pysam initialise fresh in each child
    mp.set_start_method('spawn', force=True)

    rows = list(highlow_scores.iterrows())   # [(idx, row), ...]

    worker = partial(process_row,
                     hg38_fa_path=f'{SUPREMO_DIR}/data/hg38.fa',
                     dfam_csv_path=f'{PROJECT_DIR}/data/alu_dfam_consensus.txt')

    with mp.Pool(processes=N_WORKERS) as pool:
        for idx, results in pool.imap_unordered(worker, rows, chunksize=1):
            if results:
                for col, val in results.items():
                    highlow_scores.loc[idx, col] = val

    highlow_scores.to_csv(f'{PROJECT_DIR}/results/paper_results/20260414_highlow1000_insertconsensus_withDEL.txt', sep='\t')
    print('finished')


def main():
    sys.path.append(f'{SUPREMO_DIR}/')  
    res_dir=f'{PROJECT_DIR}/results/'
    hg38_fa=pysam.FastaFile(f'{SUPREMO_DIR}/data/hg38.fa')
    alu_dfam=pd.read_csv(f'{PROJECT_DIR}/data/alu_dfam_consensus.txt', sep='\t', index_col=0)
    alu_highsamp_all=pd.read_csv(f'{res_dir}paper_results/20260413_alu_sample/aluhigh_1000samp.txt', sep='\t', index_col=0)
    alu_lowsamp_all=pd.read_csv(f'{res_dir}paper_results/20260413_alu_sample/alulow_1000samp.txt', sep='\t', index_col=0)
    alu_highsamp=alu_highsamp_all[['CHROM', 'POS', 'END', 'strand', 'orig_idx', 'repName', 'mse_mean', 'corr_mean']]
    alu_lowsamp=alu_lowsamp_all[['CHROM', 'POS', 'END', 'strand', 'orig_idx', 'repName', 'mse_mean', 'corr_mean']]
    alu_highsamp['mse_high'] = 1
    alu_lowsamp['mse_high'] = 0 
    highlow_scores=pd.concat([alu_highsamp, alu_lowsamp], axis=0)


if __name__ == '__main__':
    main()
