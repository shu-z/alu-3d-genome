"""Converted from 20260505_alu_GC_mutate.ipynb; logic unchanged."""

import os
import re
import subprocess

from pybedtools import BedTool
import numpy as np 
import pandas as pd 
import pybedtools
import pyranges as pr 
import pysam

# Paths - edit for your environment
PROJECT_DIR = os.environ.get("ALU_PROJECT_DIR", "/pollard/home/szhang20/alu")
DATA_DIR = os.environ.get("ALU_DATA_DIR", "/pollard/data/projects/shzhang")
POLLARD_DATA = os.environ.get("POLLARD_DATA", "/pollard/data")

#shuffle same bp 

#ok what do we want to do 
#take all high scoring Alus (1%) -- delete 300bp regions within the 1Mb window
    #maybe still sample these Alus? like pick 1000? 
#also grab 1000 random sequences from the not top disruptive Alus? 
#do all sequence augmentations? 
#maintain GC composition
#if in intron, stay in intron. avoid deleting another Alu

# ## read in results

data_dir=f'{DATA_DIR}/alus/data/'
resgc_dir=f'{PROJECT_DIR}/results/paper_results/20260504_GCmutate/'

os.listdir(resgc_dir)

score_df_list=[]
for window in ['1000', '10000']:

        alu_high=pd.read_csv(f'{resgc_dir}/alu_high_{window}window_10_GC_scores', sep='\t')
        print(alu_high.head())
        print(len(alu_high))

#         #merge on indices 
#         alu_score_i = pd.merge(df_alu, df_score, left_index=True, right_index=True, how='inner')

#         assert len(alu_score_i) == len(df_score)
#         assert len(alu_score_i) == len(df_alu)

#         score_df_list.append(alu_score_i)

#     random_score_df=pd.concat(score_df_list)
#     random_score_df=random_score_df.reset_index(drop=True)

#     #calculate mean, median, sd for scores per row 
#     mse_cols = [col for col in random_score_df.columns if col.startswith('mse')]
#     corr_cols = [col for col in random_score_df.columns if col.startswith('corr')]
#     df_mse = random_score_df[mse_cols]
#     df_corr = random_score_df[corr_cols]

#     # calcualte and add to alu scores 
#     random_score_df['mse_mean'] = df_mse.mean(axis=1)
#     random_score_df['mse_median'] = df_mse.median(axis=1)
#     random_score_df['mse_std'] = df_mse.std(axis=1)

#     random_score_df['corr_mean'] = df_corr.mean(axis=1)
#     random_score_df['corr_median'] = df_corr.median(axis=1)
#     random_score_df['corr_std'] = df_corr.std(axis=1)

#     #merge back to get original MSEs 
#     alu_scores_merge=alu_scores[['orig_idx', 'mse_mean']]
#     random_score_df=random_score_df.merge(alu_scores_merge,  on="original_idx", how='left')

#     random_score_df.to_csv(f'{res300_dir}alu_300bp_{window}_random_scores.txt', sep='\t')

random_score_df

random_score_df.to_csv(f'{out_dir}alu_300bp_random_scores.txt', sep='\t')

out_dir

out_dir

#repeat for 300bp random at smaller window
score_df_list=[]
for i in range(0, 10):
    df_score=pd.read_csv(f'{out_dir}output/alu_300bp_10kb_random_{i}_scores', sep='\t')
    df_alu=pd.read_csv(f'{out_dir}input/alu_300bp_10kb_random_{i}.txt', sep='\t')

    #merge on indices 
    alu_score_i = pd.merge(df_alu, df_score, left_index=True, right_index=True, how='inner')

    assert len(alu_score_i) == len(df_score)
    assert len(alu_score_i) == len(df_alu)

    score_df_list.append(alu_score_i)
    
random_score_df=pd.concat(score_df_list)
random_score_df=random_score_df.reset_index(drop=True)

random_score_df

#calculate mean, median, sd for scores per row 
mse_cols = [col for col in random_score_df.columns if col.startswith('mse')]
corr_cols = [col for col in random_score_df.columns if col.startswith('corr')]
df_mse = random_score_df[mse_cols]
df_corr = random_score_df[corr_cols]

# calcualte and add to alu scores 
random_score_df['mse_mean'] = df_mse.mean(axis=1)
random_score_df['mse_median'] = df_mse.median(axis=1)
random_score_df['mse_std'] = df_mse.std(axis=1)

random_score_df['corr_mean'] = df_corr.mean(axis=1)
random_score_df['corr_median'] = df_corr.median(axis=1)
random_score_df['corr_std'] = df_corr.std(axis=1)

#merge back to get original MSEs 
alu_scores_merge=alu_scores[['orig_idx', 'mse_mean']]
alu_scores_merge.columns=['original_idx', 'mse_mean_aluDEL']
random_score_df=random_score_df.merge(alu_scores_merge,  on="original_idx", how='left')

#merge back to get original MSEs 
random_score_df

random_score_df.to_csv(f'{out_dir}alu_300bp_10kb_random_scores.txt', sep='\t')

random_score_df

out_dir
