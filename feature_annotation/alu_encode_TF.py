"""Converted from alu_encode_TF.ipynb; logic unchanged."""

import os
import re
import subprocess

from Bio.Seq import Seq
from scipy.stats import mannwhitneyu
import matplotlib.pyplot as plt
import numpy as np 
import pandas as pd 
import pysam
import scipy.stats
import seaborn as sns

# Paths - edit for your environment
SUPREMO_DIR = os.environ.get("SUPREMO_DIR", "/pollard/home/szhang20/akita_variant_scoring")
AKITA_DIR = os.environ.get("AKITA_DIR", "/pollard/home/szhang20/akita")
PROJECT_DIR = os.environ.get("ALU_PROJECT_DIR", "/pollard/home/szhang20/alu")
DATA_DIR = os.environ.get("ALU_DATA_DIR", "/pollard/data/projects/shzhang")
POLLARD_DATA = os.environ.get("POLLARD_DATA", "/pollard/data")

#from plotnine import *

data_dir=f'{DATA_DIR}/alus/data/'

alu_res_dir=f'{PROJECT_DIR}/results/paper_results/'
alu_annot_withtracks_het=pd.read_csv(f'{alu_res_dir}20260316_merged_aluhg38_all_featureannot_NA_blacklist_filter_region_alu.txt', 
                                     sep='\t', index_col=0)

no_na_alu=alu_annot_withtracks_het.dropna(subset='mse_mean')
p99 = no_na_alu["mse_mean"].quantile(0.99)
p50 = no_na_alu["mse_mean"].quantile(0.50)

df_top1 = no_na_alu[no_na_alu["mse_mean"] >= p99]

df_top1

fasta = pysam.Fastafile(f'{SUPREMO_DIR}/data/hg38.fa')

def write_fasta(df, fasta, outpath):
    """df needs CHROM, POS, END, orig_idx columns"""
    with open(outpath, 'w') as f:
        for _, row in df.iterrows():
            seq = fasta.fetch(row['CHROM'], int(row['POS']), int(row['END']))
            f.write(f">{row['orig_idx']}\n{seq}\n")

# Define groups
top_thresh = no_na_alu['mse_mean'].quantile(0.99)
bot_thresh = no_na_alu['mse_mean'].quantile(0.50)

top = no_na_alu[no_na_alu['mse_mean'] >= top_thresh]
bot = no_na_alu[no_na_alu['mse_mean'] <= bot_thresh].sample(len(top), random_state=729)

write_fasta(top, fasta, f'{PROJECT_DIR}/results/meme_motif/high_mse_alus.fa')
write_fasta(bot, fasta, f'{PROJECT_DIR}/results/meme_motif/low_mse_alus.fa')

print(f"Written {len(top)} foreground and {len(bot)} background sequences")

fasta = pysam.Fastafile(f'{SUPREMO_DIR}/data/hg38.fa')

def write_fasta(df, fasta, outpath):
    """df needs CHROM, POS, END, orig_idx columns"""
    with open(outpath, 'w') as f:
        for _, row in df.iterrows():
            seq = fasta.fetch(row['CHROM'], int(row['POS']), int(row['END']))
            f.write(f">{row['orig_idx']}\n{seq}\n")

# Define groups
top_thresh = no_na_alu['mse_mean'].quantile(0.99)
bot_thresh = no_na_alu['mse_mean'].quantile(0.50)

top = no_na_alu[(no_na_alu['mse_mean'] >= top_thresh) & (no_na_alu['SVLEN'] >= 250)]
bot = no_na_alu[(no_na_alu['mse_mean'] <= bot_thresh) & (no_na_alu['SVLEN'] >= 250)].sample(len(top), random_state=729)

write_fasta(top, fasta, f'{PROJECT_DIR}/results/meme_motif/high_mse_alus_len250.fa')
write_fasta(bot, fasta, f'{PROJECT_DIR}/results/meme_motif/low_mse_alus_len250.fa')

print(f"Written {len(top)} foreground and {len(bot)} background sequences")

#split out again via alu family 

fasta = pysam.Fastafile(f'{SUPREMO_DIR}/data/hg38.fa')

def write_fasta(df, fasta, outpath):
    """df needs CHROM, POS, END, orig_idx, strand columns"""
    with open(outpath, 'w') as f:
        for _, row in df.iterrows():
            seq = fasta.fetch(row['CHROM'], int(row['POS']), int(row['END']))
            if row['strand'] == '-':
                seq = str(Seq(seq).reverse_complement())
            f.write(f">{row['orig_idx']}\n{seq}\n")

# Adjust this to whatever column holds the family/subfamily name
FAM_COL = 'repName'

for fam in ['AluJ', 'AluS', 'AluY']:
    sub = no_na_alu[no_na_alu[FAM_COL].str.startswith(fam, na=False)]

    top_thresh = sub['mse_mean'].quantile(0.99)
    bot_thresh = sub['mse_mean'].quantile(0.50)

    top = sub[(sub['mse_mean'] >= top_thresh) & (sub['SVLEN'] >= 250)]
    bot = sub[(sub['mse_mean'] <= bot_thresh) & (sub['SVLEN'] >= 250)].sample(len(top), random_state=729)

    write_fasta(top, fasta, f'{PROJECT_DIR}/results/meme_motif/high_mse_{fam}_len250.fa')
    write_fasta(bot, fasta, f'{PROJECT_DIR}/results/meme_motif/low_mse_{fam}_len250.fa')

    print(f"{fam}: {len(top)} foreground, {len(bot)} background")

plt.hist(top['SVLEN'], bins=30)
plt.show()

plt.hist(bot['SVLEN'], bins=30)
plt.show()

flank = 100
seq = fasta.fetch(row['CHROM'], max(0, int(row['POS']) - flank), int(row['END']) + flank)

# --- cell 13: did not parse; preserved as comments ---

# # Run STREME in discriminative mode
# streme \
#     --p high_mse_alus.fa \       # foreground (positives)
#     --n low_mse_alus.fa \        # background/control (negatives)
#     --dna \
#     --minw 6 \
#     --maxw 20 \
#     --thresh 0.05 \
#     --oc streme_output/          # output directory

# ### Parse ENCODE ChIP_seq
# #### pick ZNF, heterochromatin related TFs

#to get encode tsv information
#under 'Experiment Search', need to use 'List' as sort, and then download tsv 
chip_exp_metadata=pd.read_csv(f'{data_dir}encode/ChIP/encode_TF_ChIP_metadata.tsv', sep='\t', skiprows=1)

#marks i want

#constitutive heterochromatin 
#H3K9me3
#TRIM28 (recruited by KRABZNF)
#SETDB1 (recruited by KRABZNF)
#CBX proteins (part of HP1), binds to H3K9me3
#DNMT1

#Facultatitve heterochromatin 
#H3K9me2, H3K27me3
#EZH2, SUZ12, RNF2 #components of polycomb complex 
#REST, YY1

#BDP1

#anything that starts with ZNF

list_to_contain=['TRIM28', 'SETDB1', 'CBX', 'DNMT1', 'EZH2', 'SUZ12', 'RNF2', 'REST', 'YY1', 'ZNF', 'BDP1']
# Build a regex pattern that matches any of the substrings
pattern_tocontain = '|'.join(list_to_contain)

chip_exp_metadata['Biosample treatment'].value_counts()

#filter to those with replicates, and cell lines/primary cells  
filtered = chip_exp_metadata[
    chip_exp_metadata['Target of assay'].str.contains(pattern_tocontain) &
    chip_exp_metadata['Biological replicate'].str.contains(',') &
    (chip_exp_metadata['Biosample classification'] != 'tissue') & 
    chip_exp_metadata['Biosample treatment'].isna()]

filtered['Files']

# Get all accessions from the Files column (flatten the comma-separated lists)
valid_accessions = set(
    acc.strip().split('/')[-2]  # -2 instead of -1 due to trailing slash
    for val in filtered['Files'].dropna()
    for acc in val.split(',')
)

with open(f'{data_dir}encode/ChIP/encode_TF_ChIP.txt') as f:
    header = next(f)
    urls = [line.strip() for line in f if line.strip()]

filtered_urls = [
    url for url in urls
    if url.endswith('.bed.gz') and url.split('/')[-1].split('.')[0] in valid_accessions
]

with open(f'{data_dir}encode/ChIP/encode_TF_ChIP_bedfiltered_hetTF.txt', 'w') as f:
    f.write(header)
    f.write('\n'.join(filtered_urls))

len(filtered_urls)

chip_exp_metadata

chip_exp_metadata

# ### read in results

def rank_biserial_r(u_stat, n1, n2):
    """Effect size for Mann-Whitney U."""
    return 1 - (2 * u_stat) / (n1 * n2)

def compare_top_bottom(df, overlap_cols, metadata_df=None):
    """
    Compare top 1% vs bottom 50% by mse_mean across all overlap columns.
    
    Parameters:
        df: your alu_df_copy with mse_mean and overlap columns
        overlap_cols: list of overlap column names to test
        metadata_df: optional df with column metadata (e.g. TF name, cell type)
    """
    # Define groups
    top_thresh = df['mse_mean'].quantile(0.99)
    bot_thresh = df['mse_mean'].quantile(0.50)
    
    top = df[df['mse_mean'] >= top_thresh]
    bot = df[df['mse_mean'] <= bot_thresh]
        
    rows = []
    for col in overlap_cols:
        a = bot[col].dropna()
        b = top[col].dropna()
        
        if len(a) == 0 or len(b) == 0:
            continue
        
        u_stat, pval = mannwhitneyu(a, b, alternative='two-sided')
        rbc = rank_biserial_r(u_stat, len(a), len(b))
        
        # Parse column name to extract feature and window
        # e.g. ENCFF643LUR_overlap_10kb -> accession=ENCFF643LUR, window=10kb
        parts = col.split('_overlap_')
        accession = parts[0] if len(parts) == 2 else col
        window = parts[1] if len(parts) == 2 else 'unknown'
        
        row = {
            'column': col,
            'accession_samp': accession,
            'window': window,
            'mean_top': a.mean(),
            'mean_bot': b.mean(),
            'median_top': a.median(),
            'median_bot': b.median(),
            'mann_whitney_u': u_stat,
            'pval': pval,
            'rank_biserial_r': rbc,  # +1 = top always higher, -1 = top always lower
        }

        rows.append(row)

    results = pd.DataFrame(rows)
    
    # Multiple testing correction (Benjamini-Hochberg)
    from statsmodels.stats.multitest import multipletests
    _, pval_adj, _, _ = multipletests(results['pval'], method='fdr_bh')
    results['pval_adj'] = pval_adj
    
    results = results.sort_values('pval_adj')
    return results

test1=pd.read_csv(f'{PROJECT_DIR}/results/paper_results//aluhg38_feature_overlap_encode_ChIP_TF.csv')

test1_alu=test1.merge(no_na_alu[['orig_idx', 'mse_mean']], on='orig_idx', how='left')
test1_alu

metadata_df = chip_exp_metadata.set_index('Accession') 

overlap_cols = [c for c in test1_alu.columns if '_overlap_100kb' in c]

results_df = compare_top_bottom(test1_alu, overlap_cols, metadata_df=metadata_df)

file_to_meta = (
    chip_exp_metadata[['Accession', 'Target gene symbol', 'Biosample term name', 'Files']]
    .dropna(subset=['Files'])
    .assign(Files=lambda x: x['Files'].str.split(','))
    .explode('Files')
    .assign(file_acc=lambda x: x['Files'].str.strip().str.split('/').str[-2])
    .set_index('file_acc')
    [['Accession', 'Target gene symbol', 'Biosample term name']]
)

#merge results with metadata
results_df = results_df.merge(file_to_meta, left_on='accession_samp', right_index=True, how='left')

results_df.sort_values('rank_biserial_r', ascending=False)[0:25]

plt.hist(results_df['rank_biserial_r'], bins=50)
plt.show()

results_df.to_csv(f'{OUT_DIR}/encode_TF_overlap_mse_comparison.csv', index=False)
print(results_df.head(20))

#look at methylation overlap for different cell types

#read in methylation overlaps

sns.set_style("whitegrid")

feat_list = df_plot.columns[df_plot.columns.str.startswith('H1ESC_WGBS')]
subfamilies = ['Ancestral', 'AluJ', 'AluS', 'AluY']
groups      = ['bottom', 'top']

for feat in feat_list:
    print(feat)

    # ── all subfamilies side by side, one row per disruption group ────────────
    bins_global = np.linspace(df_plot[feat].min(), df_plot[feat].max(), 40)

    fig, axes = plt.subplots(
        len(groups), len(subfamilies),
        sharex=True, sharey=True,
        figsize=(4 * len(subfamilies), 2.5 * len(groups))
    )

    # normalise axes to always be 2D: axes[row, col]
    if len(groups) == 1 and len(subfamilies) == 1:
        axes = np.array([[axes]])
    elif len(groups) == 1:
        axes = axes[np.newaxis, :]
    elif len(subfamilies) == 1:
        axes = axes[:, np.newaxis]

    for col, subfamily in enumerate(subfamilies):
        df_sub = df_plot[df_plot['subfamily'] == subfamily]
        axes[0, col].set_title(subfamily, fontsize=10)

        for row, g in enumerate(groups):
            ax     = axes[row, col]
            subset = df_sub[df_sub['disruption'] == g]
            ax.hist(subset[feat], bins=bins_global, density=True, alpha=0.8)

            # y-axis label only on leftmost column
            if col == 0:
                ax.set_ylabel(g, fontsize=8)

    # shared x label
    fig.text(0.5, 0.01, 'Methylation Count', ha='center', fontsize=10)
    plt.suptitle(f'{feat}', y=1.01, fontsize=13)
    plt.tight_layout()
    plt.show()
