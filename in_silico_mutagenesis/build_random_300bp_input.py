"""Converted from Alu_300bp_random.ipynb; logic unchanged."""

import sys
import os
import re
import subprocess

from pybedtools import BedTool
import numpy as np 
import pandas as pd 
import pybedtools
import pyranges as pr 
import pysam
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "utils"))
from io_helpers import parse_attributes

# Paths - edit for your environment
SUPREMO_DIR = os.environ.get("SUPREMO_DIR", "/pollard/home/szhang20/akita_variant_scoring")
AKITA_DIR = os.environ.get("AKITA_DIR", "/pollard/home/szhang20/akita")
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

data_dir=f'{DATA_DIR}/alus/data/'
res_dir=f'{PROJECT_DIR}/results/'

#read in gencode and intron annotations 
#function expand info in gencode list 

pd.set_option('display.max_columns', None)

#read in gencode information 
gencode_annot_path=f'{data_dir}gencode/gencode.v29.annotation.gtf.gz'
gencode_df = pd.read_csv(gencode_annot_path, sep='\t', comment='#', header=None,compression='gzip', low_memory=False,
    names=["chrom", "source", "feature", "start", "end", "score", "strand", "score2", "info"])

#apply to gencode info column 
parsed = gencode_df['info'].apply(parse_attributes)
gencode_df = pd.concat([gencode_df, pd.DataFrame(parsed.tolist())], axis=1)

#make introns file for each gene by getting largest 'gene' feature
#then remove any place with exon, CDS, UTR

gencode_df_genes=gencode_df[gencode_df['feature']=='gene']
gencode_df_exon_CDS_UTR=gencode_df[gencode_df['feature'].isin(['exon', 'CDS', 'UTR'])]

gencode_df['feature'].value_counts()

#alu_scores_all=pd.read_csv(f'{res_dir}20260109_aluhg38_all_annot_GC_mapp_alu_repeats_genetrackoverlap.txt', sep='\t', index_col=0)
alu_scores_all=pd.read_csv(f'{res_dir}paper_results/20260316_merged_aluhg38_all_featureannot_NA_blacklist_filter_region_alu.txt', sep='\t', index_col=0)
alu_scores_all

#read in annotated/scored Alus 

#let's say we want to shuffle 5000 of the top 1% 
#the nanother 5000 randomly selected from the rest of the genome
#for each alu we shuffle 10x with regions nerably --> total 100,000 shuffles 

# alu_scores=alu_scores_all.dropna(subset='mse_mean')

# mse_1perc= alu_scores['mse_mean'].quantile(0.99)
# top1_mse = alu_scores[alu_scores['mse_mean'] >= mse_1perc]
# top1_mse_sample=top1_mse.sample(n=5000, replace=False, random_state=729)

# mse_50perc= alu_scores['mse_mean'].quantile(0.50)
# bottom50_mse=alu_scores[alu_scores['mse_mean'] <= mse_50perc]
# bottom50_mse_sample=bottom50_mse.sample(n=5000, replace=False, random_state=729)

# alu_to_shuffle=pd.concat([top1_mse_sample, bottom50_mse_sample], axis=0)

#read in previous split of high/low scoring alus 
alu_highsamp=pd.read_csv(f'{res_dir}paper_results/20260413_alu_sample/aluhigh_1000samp.txt', sep='\t', index_col=0)
alu_lowsamp=pd.read_csv(f'{res_dir}paper_results/20260413_alu_sample/alulow_1000samp.txt', sep='\t', index_col=0)
alu_to_shuffle=pd.concat([alu_highsamp, alu_lowsamp], axis=0)

alu_highsamp

alu_to_shuffle

#read in CTCFs 
track_cols=['chr', 'start', 'end', 'tf', 'rel_score', '-log10(pval)', 'strand']
chrom_list=[f'chr{i}' for i in range(1, 23)] + ['chrX']
CTCF_sites=pd.read_csv(f'{PROJECT_DIR}/data/MA0139.1.tsv', sep='\t', header=None, names=track_cols)
CTCF_sites=CTCF_sites[CTCF_sites['chr'].isin(chrom_list)]

#turn all into BEDs
col_keep=['chrom', 'start', 'end', 'feature', 'gene_type', 'gene_name']
gencode_genes_BED=pybedtools.BedTool.from_dataframe(gencode_df_genes[col_keep])
gencode_exons_BED=pybedtools.BedTool.from_dataframe(gencode_df[gencode_df['feature']=='exon'][col_keep])
gencode_exon_CDS_UTR_BED=pybedtools.BedTool.from_dataframe(gencode_df_exon_CDS_UTR[col_keep])

CTCF_BED=pybedtools.BedTool.from_dataframe(CTCF_sites[['chr', 'start', 'end']])   

alu_to_shuffle_BED=pybedtools.BedTool.from_dataframe(alu_to_shuffle[['CHROM', 'POS', 'END', 'orig_idx']])

alu_all_BED=pybedtools.BedTool.from_dataframe(alu_scores_all[['CHROM', 'POS', 'END', 'orig_idx']])

alu_overlapping_CTCF = alu_to_shuffle_BED.intersect(
    CTCF_BED, wa=True, f=0.5, u=True
)
alu_overlap_CTCF_df = alu_overlapping_CTCF.to_dataframe(names=['CHROM','POS','END','orig_idx'])

alu_overlap_CTCF_df

chrom_sizes=pd.read_csv(f'{SUPREMO_DIR}/data/chrom_lengths_hg38', sep='\t', names=['chrom','size'])
centromeres=pd.read_csv(f'{SUPREMO_DIR}/data/centromere_coords_hg38', sep='\t')
chrom_sizes['chrom']='chr'+chrom_sizes['chrom']

# Intersect: keep only alu regions that overlap genes
alu_overlapping_genes = alu_to_shuffle_BED.intersect(
    gencode_genes_BED, wa=True, f=0.5, u=True
)

# Intersect: keep only alu regions that do NOT overlap genes
alu_nonoverlapping_genes = alu_to_shuffle_BED.intersect(
    gencode_genes_BED, f=0.5, v=True
)

# Convert back to pandas dataframes
alu_overlap_df = alu_overlapping_genes.to_dataframe(names=['CHROM','POS','END','orig_idx'])
alu_nonoverlap_df = alu_nonoverlapping_genes.to_dataframe(names=['CHROM','POS','END','orig_idx'])

alu_nonoverlap_df

def get_gc_content(genome_fasta, chrom, start, end):
    """Calculate GC content for a genomic region using pysam."""
    seq = genome_fasta.fetch(chrom, start, end).upper()
    gc_count = seq.count('G') + seq.count('C')
    total = len(seq)
    return gc_count / total if total > 0 else 0

def fast_shuffle_list_logic(alu_df, chromsizes_df=None, must_overlap=None, must_not_overlap=None, 
                             n_shuffles=10, slop=100_000):
    
    def merge_to_pr(df_list):
        if not df_list: return None
        
        standardized_dfs = []
        for df in df_list:
            # Take ONLY the first 3 columns (chr, start, end) 
            subset = df.iloc[:, :3].copy()
            subset.columns = ['Chromosome', 'Start', 'End']
            standardized_dfs.append(subset)
        
        combined = pd.concat(standardized_dfs)
        return pr.PyRanges(combined).merge()

    def genome_pr_from_chromsizes(chromsizes_df):
        # chromsizes_df columns: Chromosome, End
        genome = chromsizes_df.copy()
        genome.columns = ['Chromosome', 'End']
        genome["Start"] = 0
        genome = genome[["Chromosome", "Start", "End"]]
        return pr.PyRanges(genome)

    # Create master 'Allowed' and 'Forbidden' maps
    mo_pr = merge_to_pr(must_overlap)
    mno_pr = merge_to_pr(must_not_overlap)

    # 2. Pre-calculate the "Valid Universe" 
    if mo_pr is not None:
        valid_universe = mo_pr
        if mno_pr is not None:
            valid_universe = valid_universe.subtract(mno_pr)
    elif mno_pr is not None:
        # If no must_overlap is provided, get valid from chromsizes 
        valid_universe = genome_pr_from_chromsizes(chromsizes_df)
        valid_universe = valid_universe.subtract(mno_pr)
    else:
        valid_universe = genome_pr_from_chromsizes(chromsizes_df)

    results = []

    # Process each gene
    for idx, row in alu_df.iterrows():
        print(idx)
        chrom = row['CHROM']
        orig_start = row['POS']
        orig_end = row['END']
        length = orig_end - orig_start
        
        # 3. Define the local slop window
        window_start = max(0, orig_start - slop)
        window_end = orig_end + slop
        window_pr = pr.from_dict({"Chromosome": [chrom], "Start": [window_start], "End": [window_end]})
        
        # 4. Find valid intervals in this slop window
        if valid_universe is not None:
            valid_intervals = valid_universe.intersect(window_pr).as_df()
        else:
            # If no constraints, the whole window is valid
            valid_intervals = window_pr.as_df()
        
        if valid_intervals.empty:
            continue
            
        # Adjust: The start position must allow the whole length to fit
        valid_intervals['End'] = valid_intervals['End'] - length
        valid_intervals = valid_intervals[valid_intervals['End'] > valid_intervals['Start']]
        
        if valid_intervals.empty:
            continue

        # 5. Fast Weighted Sampling
        valid_intervals['weights'] = valid_intervals['End'] - valid_intervals['Start']
        total_weight = valid_intervals['weights'].sum()
        
        # Sample n positions based on length of available segments
        chosen_indices = np.random.choice(valid_intervals.index, size=n_shuffles, 
                                          p=valid_intervals['weights'] / total_weight)
        
        for i, interval_idx in enumerate(chosen_indices):
            interval = valid_intervals.loc[interval_idx]
            new_start = np.random.randint(interval['Start'], interval['End'])
            
            results.append({
                'CHROM': chrom,
                'POS': new_start,
                'END': new_start + length,
                'original_idx': row['orig_idx'],
                'shuffle_num': i + 1
            })
            
    return pd.DataFrame(results)

# ### 10kb window

hg38_fa_path = f'{SUPREMO_DIR}/data/hg38.fa'
n_shuff=10
window=10_000

#all must not overlap other Alus, CTCF sites, 

#separate regions that do and don't overlap genes 
random_overlap_genes = fast_shuffle_list_logic(alu_overlap_df,
    must_overlap=[gencode_df_genes[col_keep]],  
    must_not_overlap=[alu_scores_all[['CHROM', 'POS', 'END', 'orig_idx']], CTCF_sites[['chr', 'start', 'end']], centromeres], 
    n_shuffles=n_shuff,slop=window)

#separate regions that do and don't overlap genes 
random_nooverlap_genes = fast_shuffle_list_logic(alu_nonoverlap_df,                                              
    chromsizes_df=chrom_sizes, must_overlap=None,  
    must_not_overlap=[gencode_df_genes[col_keep], alu_scores_all[['CHROM', 'POS', 'END', 'orig_idx']], 
                      CTCF_sites[['chr', 'start', 'end']], centromeres], 
    n_shuffles=n_shuff,slop=window)

#add in column showing where from and combine 
random_overlap_genes['shuffle_gene_overlap']=True
random_nooverlap_genes['shuffle_gene_overlap']=False

all_shuffled=pd.concat([random_overlap_genes, random_nooverlap_genes], axis=0, ignore_index=True)
all_shuffled['REF']='-'
all_shuffled['ALT']='-'
all_shuffled['SVTYPE']='DEL'
all_shuffled['SVLEN']=all_shuffled['END']-all_shuffled['POS']

col_first=['CHROM', 'POS', 'END', 'REF', 'ALT', 'SVTYPE', 'SVLEN', 'original_idx']
remaining_columns = [col for col in all_shuffled.columns if col not in col_first]
col_order = col_first + remaining_columns

all_shuffled=all_shuffled[col_order]
all_shuffled

out_dir=f'{PROJECT_DIR}/results/paper_results/300bp_random/10kb/'
os.makedirs(f'{out_dir}input/', exist_ok=True)
os.makedirs(f'{out_dir}output/', exist_ok=True)
os.makedirs(f'{out_dir}logs/', exist_ok=True)
chunk_size=1000

for i, start in enumerate(range(0, len(all_shuffled), chunk_size)):

    print(i, start)

    chunk = all_shuffled[start : start + chunk_size]
    chunk.to_csv(f"{out_dir}input/alu_300bp_10kb_random_{i}.txt", sep = '\t',index=False)

# ### repeat for 100kb

n_shuff=10
window=100_000

#all must not overlap other Alus, CTCF sites, 

#separate regions that do and don't overlap genes 
random_overlap_genes = fast_shuffle_list_logic(alu_overlap_df,
    must_overlap=[gencode_df_genes[col_keep]],  
    must_not_overlap=[alu_scores_all[['CHROM', 'POS', 'END', 'orig_idx']], CTCF_sites[['chr', 'start', 'end']], centromeres], 
    n_shuffles=n_shuff,slop=window)

#separate regions that do and don't overlap genes 
random_nooverlap_genes = fast_shuffle_list_logic(alu_nonoverlap_df,                                              
    chromsizes_df=chrom_sizes, must_overlap=None,  
    must_not_overlap=[gencode_df_genes[col_keep], alu_scores_all[['CHROM', 'POS', 'END', 'orig_idx']], 
                      CTCF_sites[['chr', 'start', 'end']], centromeres], 
    n_shuffles=n_shuff,slop=window)

#add in column showing where from and combine 
random_overlap_genes['shuffle_gene_overlap']=True
random_nooverlap_genes['shuffle_gene_overlap']=False

all_shuffled=pd.concat([random_overlap_genes, random_nooverlap_genes], axis=0, ignore_index=True)
all_shuffled['REF']='-'
all_shuffled['ALT']='-'
all_shuffled['SVTYPE']='DEL'
all_shuffled['SVLEN']=all_shuffled['END']-all_shuffled['POS']

col_first=['CHROM', 'POS', 'END', 'REF', 'ALT', 'SVTYPE', 'SVLEN', 'original_idx']
remaining_columns = [col for col in all_shuffled.columns if col not in col_first]
col_order = col_first + remaining_columns

all_shuffled=all_shuffled[col_order]
all_shuffled

out_dir=f'{PROJECT_DIR}/results/paper_results/300bp_random/100kb/'
os.makedirs(f'{out_dir}', exist_ok=True)
os.makedirs(f'{out_dir}input/', exist_ok=True)
os.makedirs(f'{out_dir}output/', exist_ok=True)
os.makedirs(f'{out_dir}logs/', exist_ok=True)
chunk_size=1000

for i, start in enumerate(range(0, len(all_shuffled), chunk_size)):

    print(i, start)

    chunk = all_shuffled[start : start + chunk_size]
    chunk.to_csv(f"{out_dir}input/alu_300bp_100kb_random_{i}.txt", sep = '\t',index=False)

os.listdir(f'{res300_dir}/{window}/')

# ## read in results

alu_highsamp=pd.read_csv(f'{res_dir}paper_results/20260413_alu_sample/aluhigh_1000samp.txt', sep='\t', index_col=0)
alu_lowsamp=pd.read_csv(f'{res_dir}paper_results/20260413_alu_sample/alulow_1000samp.txt', sep='\t', index_col=0)

res300_dir=f'{PROJECT_DIR}/results/paper_results/300bp_random/'
for window in ['10kb', '100kb']:

    score_df_list=[]
    for i in range(0, 20):
        df_score=pd.read_csv(f'{res300_dir}/{window}/output/alu_300bp_{window}_random_{i}_scores', sep='\t')
        df_alu=pd.read_csv(f'{res300_dir}/{window}/input/alu_300bp_{window}_random_{i}.txt', sep='\t')

        #merge on indices 
        alu_score_i = pd.merge(df_alu, df_score, left_index=True, right_index=True, how='inner')

        assert len(alu_score_i) == len(df_score)
        assert len(alu_score_i) == len(df_alu)

        score_df_list.append(alu_score_i)

    random_score_df=pd.concat(score_df_list)
    random_score_df=random_score_df.reset_index(drop=True)

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
    alu_scores_merge=alu_scores_all[['orig_idx', 'mse_mean']]
    alu_scores_merge.columns=['original_idx', 'mse_mean_aluDEL']
    random_score_df=random_score_df.merge(alu_scores_merge,  on="original_idx", how='left')

    #add disruption label back in 
    random_score_df['disruption'] = np.where(
        random_score_df['original_idx'].isin(alu_highsamp['orig_idx']), 'high',
        np.where(random_score_df['original_idx'].isin(alu_lowsamp['orig_idx']), 'low', np.nan)
    )

    random_score_df.to_csv(f'{res300_dir}alu_300bp_{window}_random_scores_combined.txt', sep='\t')

random_score_df['disruption'].value_counts()

len(np.unique(random_score_df['original_idx']))
