

import os
from itertools import combinations
import gzip
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
        alignment_score=alignments[0].score/len(alu_new)
    
        #get just percentage of alignment  
        #identity_percentage=pairwise2.format_alignment(*alignments[0]).count('|') / len(alu_new)    

    if global_local=='local':
        alignments = pairwise2.align.localxx(alu_new, alu_consensus, one_alignment_only=True)
        #print(alignments)
    
        #scale alignment score by length of new alu sequence 
        alignment_score=alignments[0].score/len(alu_new)
    
        #get just percentage of alignment  
        #identity_percentage=pairwise2.format_alignment(*alignments[0]).count('|') / len(alu_new)     

    #alignment score and identity percentage end up being the same 
    return alignment_score

##################### read in necessary files 

#get top scores 


#let's look at the top 5% of either metric

#add in annotations if high in mse, corr, or both 


#RUN FOR EACH ROW 


def main():
    sys.path.append(f'{SUPREMO_DIR}/')  
    hg38_fa=pysam.FastaFile(f'{SUPREMO_DIR}/data/hg38.fa')
    alu_dfam=pd.read_csv(f'{DATA_DIR}/alus/data/dfam/alu_dfam_consensus.txt', sep='\t', index_col=0)
    alu_scores=pd.read_csv(f'{PROJECT_DIR}/results/202520513_hg38Alus_100sample_repName_DEL_combinedscores.txt', 
                           sep='\t', index_col=0)
    q90_corr = alu_scores['corr_mean'].quantile(0.10)
    q90_mse = alu_scores['mse_mean'].quantile(0.90)
    scores_top_mse = alu_scores[alu_scores['mse_mean'] >= q90_mse].copy()
    scores_top_corr = alu_scores[alu_scores['corr_mean'] <= q90_corr].copy()
    scores_top=pd.concat([scores_top_mse, scores_top_corr])
    scores_top=scores_top.drop_duplicates()
    scores_top = scores_top.reset_index()
    scores_top['mse_high']= np.where(scores_top['mse_mean'] >= q90_mse, 1, 0)
    scores_top['corr_high']= np.where(scores_top['corr_mean'] <= q90_corr, 1, 0)
    print('len scores top', len(scores_top))
    for idx, row in scores_top.iterrows():

        revcomp=False  
        print(idx)
        try:

            chrom, alu_start, alu_end, alu_len, alu_strand= row['CHROM'], row['POS'], row['END'], row['SVLEN'], row['strand']        
            alu_mid=np.floor((alu_start+alu_end)/2)
            alu_seq=hg38_fa.fetch(chrom, alu_start, alu_end)      

            #get REF seq and akita pred 
            REF_seq = hg38_fa.fetch(chrom, int(np.floor(alu_mid-MB/2)), int(np.floor(alu_mid+MB/2))).upper()
            REF_pred=utils.vector_from_seq(REF_seq)
            #REF_pred_mat=utils.mat_from_vector(REF_pred)

            ########################
            #run insertion once for Alu in each family 

            #need to check if Alu is on +/- strand! If -, replace with RC seq 
            if alu_strand=='-':
                revcomp=True

            for dfam_idx, row in alu_dfam.iterrows():
                try:
                    alu_name, consensus_seq=row['AluName'], row['Consensus_seq'] 

                    if revcomp:
                        consensus_seq = str(Seq(consensus_seq).reverse_complement())          

                    #get alignment
                    align_score=align_alu_consensus(alu_seq.upper(), consensus_seq.upper(), global_local='global')

                    #get length and padding 
                    len_alu=len(consensus_seq)
                    region_pad=(MB-len_alu)/2

                    #get ALT seq 
                    left_start=np.floor(alu_start-region_pad)
                    right_end=np.floor(alu_end+region_pad)
                    left_seq_INS = hg38_fa.fetch(chrom, int(left_start), int(alu_start))
                    right_seq_INS = hg38_fa.fetch(chrom, int(alu_end), int(right_end))
                    ALT_insert_seq= (left_seq_INS + consensus_seq + right_seq_INS).upper()

                    #make ALT prediction 
                    ALT_insert_pred=utils.vector_from_seq(ALT_insert_seq)
                    #ALT_insert_pred_mat=utils.mat_from_vector(ALT_insert_pred)

                    #get MSE and CORR
                    MSE_insert=scoring_utils.mse(REF_pred, ALT_insert_pred)
                    Spearman_insert=scoring_utils.spearman(REF_pred, ALT_insert_pred)

                    #add to original df 
                    scores_top.loc[idx, f'mse_{alu_name}']=MSE_insert   
                    scores_top.loc[idx, f'corr_{alu_name}']=Spearman_insert 
                    scores_top.loc[idx, f'alignscore_{alu_name}']=align_score   

                except:
                    print('error in inserting consensus Alu at index ', 'dfam_idx')

        except:
            print('error in Alu at index ', idx)
    scores_top.to_csv(f'{DATA_DIR}/alus/20251022_topalu_insertAlu.txt', sep='\t')
    print('finished')


if __name__ == '__main__':
    main()
