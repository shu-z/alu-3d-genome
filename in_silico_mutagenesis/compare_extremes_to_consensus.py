"""Converted from MSE_extremes_vs_Alu_consensus.ipynb; logic unchanged."""

import gzip
import os
import subprocess

from Bio import SeqIO
import matplotlib.pyplot as plt
import numpy as np 
import pandas as pd
import seaborn as sns

import alu_consensus_utils as utils 

# Paths - edit for your environment
PROJECT_DIR = os.environ.get("ALU_PROJECT_DIR", "/pollard/home/szhang20/alu")

# %load_ext autoreload
# %autoreload 2

#ok, so what do we want to do
#find the assigned consensus sequence for each Alu
#see how much it deviates from sequence 
#try calling features like promoter box, etc? 
#plot if there are significant differences between the top and bottom Alus
#KS test? difference of two distributions

#data information 
res_dir=f'{PROJECT_DIR}/results/'

top_bottom=pd.read_csv(f'{res_dir}20240129/DEL_mse_topbottom.txt', sep='\t', index_col=0)
alu_consensus=pd.read_csv(f'{res_dir}20240129/20240129_alu_consensus_seq_data.txt', sep='\t', index_col=0)

top_bottom

alu_consensus

def analyze_gc_seqs(var_index, alu_seq):

    #get indices of 5' A promoter box TGGCTCACGCC
    A_promoter_idx=alu_seq.find('TGGCTCACGCC')

    #get indices of 5' B promoter box 3’ B box  GNTCGAGAC
    B_promoter_idx=utils.find_substring_with_exception(alu_seq, 'GNTCGAGAC')

    #get indices of split middle thing AAAAATACAAAAAA
    mid_box_start=alu_seq.find('AAAAATACAAAAA')

    #get indices of polya tail
    #we are going to hope there is only one otherwise this will return many indices
    polyA_len = utils.count_polyA(alu_seq)
    polyA_idx= len(alu_seq)-polyA_len

    #get GC content of each half/arm
    total_GC=utils.calculate_gc(alu_seq, 0, len(alu_seq))
    RA_GC=-1
    LA_GC=-1
    RA_len=-1
    LA_len=-1
    
#     print('mid box start', mid_box_start)
#     print('polyA', polyA_idx)
#     print('alu len', len(alu_seq))
    if mid_box_start!= -1:
        LA_GC=utils.calculate_gc(alu_seq, 0, mid_box_start)
        LA_len=mid_box_start
        if polyA_idx!=-1 and polyA_idx >(mid_box_start+13):
            RA_GC=utils.calculate_gc(alu_seq, (mid_box_start+13), polyA_idx)
            RA_len=polyA_idx-(mid_box_start+13)

    AB_promoter_GC=-1
    Bpromoter_mid_GC=-1
    if A_promoter_idx!=-1 and B_promoter_idx!=-1:
        AB_promoter_GC=utils.calculate_gc(alu_seq, A_promoter_idx, B_promoter_idx)
        
    if B_promoter_idx!=-1 and mid_box_start!=-1:
        Bpromoter_mid_GC=utils.calculate_gc(alu_seq, B_promoter_idx, mid_box_start)

    return (var_index, len(alu_seq), A_promoter_idx, B_promoter_idx, mid_box_start, polyA_len, polyA_idx, 
             LA_len, RA_len, total_GC, LA_GC, RA_GC, AB_promoter_GC, Bpromoter_mid_GC)

top_bottom.head()

alu_consensus.head()

alu_gcanalyze_list=[]
for i in top_bottom.index:
    #print(i)
    alu = top_bottom.loc[i]

    #get alu sequence 
    alu_seq=utils.get_seq_hg38(alu.CHROM, alu.POS, alu.END)

    #pull out features of promoter box, etc 
    gcanalyze_df=analyze_gc_seqs(alu.var_index, alu_seq)
    alu_gcanalyze_list.append(gcanalyze_df)

    #find matching Alu from consensus
    consensus_match=alu_consensus[alu_consensus['Alu_name'] == alu.repName]
    
    if len(consensus_match)>0:
        #run alignment between consensus and alu
        identity_percentage, alignment_score= utils.align_alu_consensus(alu_seq, consensus_match.Alu_seq.values[0])
        top_bottom.loc[i, 'align_identity_percent']=identity_percentage
        top_bottom.loc[i, 'alignment_score']=alignment_score

    else:
        pass

print(pd.DataFrame(top_bottom.iloc[412]))

colnames=['var_index', 'Alu_len', 'A_promoter_idx', 'B_promoter_idx', 'mid_box_start', 'polyA_len', 'polyA_idx', 
         'LA_len', 'RA_len',  'total_GC','LA_GC',  'RA_GC', 'AB_promoter_GC', 'Bpromoter_mid_GC']
gc_analyze_df=pd.DataFrame(alu_gcanalyze_list, columns=colnames)

#combine these two dfs basedon the var index
consensus_analyzed_topbottom=pd.merge(top_bottom, gc_analyze_df, on='var_index')
consensus_analyzed_topbottom
consensus_analyzed_topbottom.replace(-1, np.nan, inplace=True)

consensus_analyzed_topbottom

consensus_analyzed_topbottom.columns

#just look at a few things quickly 

value_to_plot='align_identity_percent'

fig, (ax1, ax2) = plt.subplots(1, 2, sharex=True, gridspec_kw={'width_ratios': [4, 4]}, figsize=(12, 6))

sns.histplot(data=consensus_analyzed_topbottom, ax=ax1,
             x=value_to_plot, hue='score_category', 
             bins=50, element="step", stat="density", common_norm=False)
ax1.set_xlabel('Frequency')

sns.boxplot(data=consensus_analyzed_topbottom, ax=ax2, orient='v',
             y=value_to_plot, x='score_category')

# Display the histogram
plt.show()

consensus_analyzed_topbottom.to_csv(f'{PROJECT_DIR}/results/20240208/20240208_alu_topbottom_vs_consensus.txt', sep='\t')
