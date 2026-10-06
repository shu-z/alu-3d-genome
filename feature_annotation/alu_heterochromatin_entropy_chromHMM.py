"""Converted from alu_heterochromatin_entropy_chromHMM.ipynb; logic unchanged."""

from functools import reduce
import glob
import os
import random
import re
import subprocess

from pybedtools import BedTool, cleanup
from scipy.stats import mannwhitneyu, chi2_contingency
import matplotlib.pyplot as plt
import numpy as np 
import pandas as pd 
import pybedtools
import scipy.stats
import seaborn as sns

# Paths - edit for your environment
PROJECT_DIR = os.environ.get("ALU_PROJECT_DIR", "/pollard/home/szhang20/alu")
DATA_DIR = os.environ.get("ALU_DATA_DIR", "/pollard/data/projects/shzhang")
POLLARD_DATA = os.environ.get("POLLARD_DATA", "/pollard/data")

data_dir=f'{DATA_DIR}/alus/data/'
h3k27me3_dir=f'{data_dir}encode/H3K27me3/'

colnames=['chrom', 'start', 'end', 'name', 'score', 'strand', 'signal', 'pval', 'qval', 'peak_summit']
bed_files=[f for f in os.listdir(h3k27me3_dir) if f.endswith('.bed.gz')]

colnames=['chrom', 'start', 'end', 'name', 'score', 'strand', 'signal', 'pval', 'qval', 'peak_summit']

test1=pd.read_csv(f'{h3k27me3_dir}ENCFF998IEU.bed.gz', sep='\t', header=None, names=colnames)
test1['SVLEN']=test1['end']-test1['start']

plt.hist(test1['signal'], bins=30)

plt.hist(test1[test1['SVLEN']<5000]['SVLEN'], bins=30)

plt.scatter(test1['SVLEN'], test1['signal'])

#h3k27me3 file metadata 
metadata=pd.read_csv(f'{h3k27me3_dir}H3K27me3_ChIP_primary_cellline_differentiated_metadata.tsv', sep='\t', skiprows=1)
metadata[metadata['Files'].str.contains('ENCFF998IEU')]

#cols to keep from metadata 
cols_to_keep = ['Accession', 'Assay name', 'Assay title','Biosample classification', 'Target', 'Target of assay',
       'Target gene symbol', 'Biosample summary', 'Biosample term name','Description', 'Lab', 'Project', 'Status', 'Files',
       'Biosample accession', 'Biological replicate',
       'Technical replicate', 'Organism', 'Life stage',
       'Biosample age', 'Biosample treatment']  

rows = []
for fname in bed_files:
    prefix = fname.replace('.bed.gz', '') 
    
    # Find rows in metadata where the Files column contains this prefix
    mask = metadata['Files'].str.contains(prefix, na=False)
    matched = metadata[mask]
    
    if not matched.empty:
        row = matched[cols_to_keep].iloc[0].to_dict()
    else:
        row = {col: None for col in cols_to_keep}
    
    row = {'filename': fname, **row}
    rows.append(row)

result_df = pd.DataFrame(rows, columns=['filename'] + cols_to_keep)

result_df

result_df['Biosample term name'].values

#read in alus with overlap scores 
RES_DIR = f'{DATA_DIR}/alus/results/20251125_aluhg38_all_scores/'
OUT_DIR = f'{PROJECT_DIR}/results/20260208_feature_overlap_all/'

alu_annot = pd.read_csv(f'{RES_DIR}20260123_aluhg38_all_annot_GC_mapp_alu_repeats_genetrackoverlap.txt', sep='\t', index_col=0)

print(alu_annot.columns.values)

os.listdir(OUT_DIR)

alu_annot = pd.read_csv(f'{RES_DIR}20260123_aluhg38_all_annot_GC_mapp_alu_repeats_genetrackoverlap.txt', sep='\t', index_col=0)

bw_avg=pd.read_csv(f'{OUT_DIR}aluhg38_feature_averageBW.csv')
GC_avg=pd.read_csv(f'{OUT_DIR}aluhg38_feature_GC.csv')
feat_overlap=pd.read_csv(f'{OUT_DIR}aluhg38_feature_overlap_bed.csv')
feat_count=pd.read_csv(f'{OUT_DIR}aluhg38_feature_count_bed.csv')
repeat_1mb=pd.read_csv(f'{OUT_DIR}aluhg38_feature_overlap_repeatsonly.csv')
repeat_other=pd.read_csv(f'{OUT_DIR}aluhg38_feature_overlap_repeatsonly_short.csv')
phylop_1mb=pd.read_csv(f'{OUT_DIR}aluhg38_phyloP_countBW_1Mb.csv')

#merge original alu with all annotations 
alu_annot_subset=alu_annot[['CHROM', 'POS', 'END', 'repName', 'orig_idx', 'mse_mean', 'corr_mean', 
                            'REF_start', 'REF_stop', 'var_rel_pos_REF', 
                            'single_gene_name', 'single_gene_type', 'single_gene_overlap_bp',
                           'swScore', 'milliDel', 'id','#bin', 'milliIns', 'milliDiv', 'SVLEN']]
df_list=[alu_annot_subset, GC_avg, bw_avg, feat_count, feat_overlap, repeat_other, repeat_1mb, phylop_1mb]

merged = reduce(lambda l, r: pd.merge(l, r, on=['CHROM', 'POS', 'END', 'orig_idx'], how="outer"), df_list)

#add in disruption
merged_nona=merged[~merged['mse_mean'].isna()]
q99=np.quantile(merged_nona['mse_mean'], 0.99)
q50=np.quantile(merged_nona['mse_mean'], 0.50)
merged_nona['disruption']=np.where(merged_nona['mse_mean'] >=q99, 'top', 
                            np.where(merged_nona['mse_mean'] <=q50, 'bottom', 'middle'))

#merge original alu with all annotations 
alu_annot_subset=alu_annot[['CHROM', 'POS', 'END', 'repName', 'orig_idx', 'mse_mean', 'corr_mean', 
                            'REF_start', 'REF_stop', 'var_rel_pos_REF', 
                            'single_gene_name', 'single_gene_type', 'single_gene_overlap_bp',
                           'swScore', 'milliDel', 'id','#bin', 'milliIns', 'milliDiv', 'SVLEN']]
df_list=[alu_annot_subset, GC_avg, bw_avg, feat_count, feat_overlap, repeat_other, repeat_1mb, phylop_1mb]

merged = reduce(lambda l, r: pd.merge(l, r, on=['CHROM', 'POS', 'END', 'orig_idx'], how="outer"), df_list)

repeat_other

#add in disruption
merged_nona=merged[~merged['mse_mean'].isna()]
q99=np.quantile(merged_nona['mse_mean'], 0.99)
q50=np.quantile(merged_nona['mse_mean'], 0.50)
merged_nona['disruption']=np.where(merged_nona['mse_mean'] >=q99, 'top', 
                            np.where(merged_nona['mse_mean'] <=q50, 'bottom', 'middle'))

#what do we want to do 
#we want to be able to separate high/low facultative heterochromatin marks
#also want to separate high/low disruption. look at distributions for both

#want to shw that high h3k27me3 and high disruption 
#want to assess if these are also cell-type specific heterochromatin regions 
#could lead to hypothesis that alus here help with cell-type specific folding 

#start with overlaps for 10 cell lines first
#do we care about all windows still? could simplify which overlaps we run

#what about overall intensity of peakiness? -- this might be more relevant for CTCF
    #could we categorize into low, medium, high peak? quantile of peak value per sample?

high_mse=merged_nona[merged_nona['disruption']=='top']
low_mse=merged_nona[merged_nona['disruption']=='bottom']

high_mse

plt.hist(high_mse['HFF_H3K27me3_overlap_1Mb'])

cutoff=0.10
high_mse_high_h3k27me3=high_mse[high_mse['HFF_H3K27me3_overlap_1Mb']>cutoff]
high_mse_low_h3k27me3=high_mse[high_mse['HFF_H3K27me3_overlap_1Mb']<cutoff]
order = sorted(high_mse_high_h3k27me3['repName'].dropna().unique())

prop = (high_mse_high_h3k27me3['repName'].value_counts(normalize=True).sort_index().reset_index())
prop.columns = ['column_name', 'proportion']

plt.figure(figsize=(12, 5))
sns.barplot(data=prop, x='column_name', y='proportion')
plt.xticks(rotation=45)
plt.ylabel("Proportion")
plt.ylim(0, 0.15)
plt.tight_layout()
plt.show()

prop = (high_mse_low_h3k27me3['repName'].value_counts(normalize=True).sort_index().reset_index())
prop.columns = ['column_name', 'proportion']

plt.figure(figsize=(12, 5))
sns.barplot(data=prop, x='column_name', y='proportion')
plt.xticks(rotation=45)
plt.ylabel("Proportion")
plt.ylim(0, 0.15)
plt.tight_layout()
plt.show()

def analyze_region(region_row, bed_files, bed_dir, n_bins=200):
    """
    For a single region, compute:
    - coverage fraction per cell type
    - summary stats across cell types
    - binned coverage profile per cell type
    - spatial consistency metrics
    """
    chrom, start, end = region_row['CHROM'], int(region_row['REF_start_x']), int(region_row['REF_stop_x'])
    region_len = end - start
    region_bt  = pybedtools.BedTool(f"{chrom}\t{start}\t{end}", from_string=True)
    
    # Bin the region into n_bins equal windows
    bins = pybedtools.BedTool.window_maker(region_bt, b=region_bt, n=n_bins)
    
    coverage_fracs = {}   # total coverage fraction per cell
    bin_profiles   = {}   # per-bin coverage profile per cell
    
    for fname in bed_files:
        cell = fname.replace('.bed.gz', '')
        bed  = pybedtools.BedTool(os.path.join(bed_dir, fname))
        
        # Clip peaks to this region only
        clipped = bed.intersect(region_bt, u=True)
        
        if len(clipped) == 0:
            coverage_fracs[cell] = 0.0
            bin_profiles[cell]   = np.zeros(n_bins)
            continue
        
        # Total coverage fraction for the whole region
        cov = region_bt.coverage(clipped)
        df  = cov.to_dataframe(names=['chr','start','end','count','bp','len','frac'])
        coverage_fracs[cell] = df['frac'].values[0]
        
        # Per-bin coverage profile
        bin_cov = bins.coverage(clipped)
        bin_df  = bin_cov.to_dataframe(names=['chr','start','end','count','bp','len','frac'])
        bin_profiles[cell] = bin_df['frac'].values
    
    #summarize across cell types 
    frac_vec = np.array(list(coverage_fracs.values()))
    mean_cov = frac_vec.mean()
    cv       = frac_vec.std() / mean_cov if mean_cov > 0 else 0.0
    
    # Entropy over coverage fractions
    def coverage_entropy(vals, bins=10):
        if vals.sum() == 0:
            return 0.0
        counts, _ = np.histogram(vals, bins=bins, range=(0, 1))
        probs = counts / counts.sum()
        probs = probs[probs > 0]
        return -np.sum(probs * np.log2(probs))
    
    ent = coverage_entropy(frac_vec)
    
    # 'spatial' metrics 
    profile_matrix = np.vstack(list(bin_profiles.values()))  # shape: (n_cells, n_bins)
    
    # Mean profile across cell types — where in the region is methylation on average
    mean_profile = profile_matrix.mean(axis=0)
    print(mean_profile)
    
    # Per-bin std — high std in a bin means cell types disagree about that sub-region
    bin_std = profile_matrix.std(axis=0)
    
    # Spatial entropy — entropy of the mean profile (is methylation spread or focal?)
    # High = spread evenly across the region, Low = concentrated in one sub-region
    spatial_entropy = coverage_entropy(mean_profile)
    
    # Pairwise correlation between cell type profiles
    # High mean correlation = cell types methylate the same sub-regions
    # Low = they methylate different parts
    if profile_matrix.shape[0] > 1:
        corr_matrix = np.corrcoef(profile_matrix)
        # Take upper triangle (pairwise correlations, excluding self)
        upper = corr_matrix[np.triu_indices_from(corr_matrix, k=1)]
        mean_pairwise_corr = np.nanmean(upper)
    else:
        mean_pairwise_corr = np.nan
    
    return {
        'chr': chrom, 'start': start, 'end': end,
        'orig_idx': row['orig_idx'], 'mse_mean' : row['mse_mean'],
        'HFF_H3K27me3_overlap_1Mb': row['HFF_H3K27me3_overlap_1Mb'],
        # Per cell coverage
        **{f'frac_{cell}': frac for cell, frac in coverage_fracs.items()},
        # Summary
        'mean_coverage':       mean_cov,
        'cv_coverage':         cv,           # CV is coeffiecient of variation, high = variable total coverage across cells
        'entropy_coverage':    ent,           # high = heterogeneous coverage levels
        # Spatial
        'spatial_entropy':     spatial_entropy,    # high = spread across region, low = focal
        'mean_pairwise_corr':  mean_pairwise_corr, # high = same location, low = different parts
        'mean_bin_std':        bin_std.mean(),      # average disagreement per bin
        # Profiles (for plotting)
        '_mean_profile':       mean_profile,
        '_bin_std':            bin_std,
        '_profile_matrix':     profile_matrix,
    }

query_df_unique

#filter to small subset 
query_filelist=['K562', 'HCT116', 'GM12878', 'H1', 'foreskin fibroblast', 'neuron', 'cardiac muscle cell']
query_df=result_df[(result_df['Biosample term name'].isin(query_filelist)) & (result_df['Biological replicate']!='1')]
#filter out if same Biosample term name
query_df_unique=query_df.drop_duplicates(subset=['Biosample term name'])
query_beds=query_df_unique['filename'].tolist()

high_mse_samp=high_mse.sample(n=1000, random_state=729)
low_mse_samp=low_mse.sample(n=1000, random_state=729)

# Run across all regions
results_test = []
for idx, row in high_mse_samp[0:2].iterrows():
    print(idx)
    r = analyze_region(row, query_beds, h3k27me3_dir)
    results.append(r)

# Separate scalar results from array results for the dataframe
scalar_cols = [k for k in results[0].keys() if not k.startswith('_')]
results_df_test  = pd.DataFrame([{k: r[k] for k in scalar_cols} for r in results_test])

# Run across all regions
results = []
for idx, row in high_mse_samp.iterrows():
    print(idx)
    r = analyze_region(row, query_beds, h3k27me3_dir)
    results.append(r)

# Separate scalar results from array results for the dataframe
scalar_cols = [k for k in results[0].keys() if not k.startswith('_')]
results_df  = pd.DataFrame([{k: r[k] for k in scalar_cols} for r in results])

results_df.to_csv(f'{PROJECT_DIR}/results/20260226_h3k27me3_entropy_highmse.txt', sep='\t')

# Run across all regions
results = []
for idx, row in low_mse_samp.iterrows():
    print(idx)
    r = analyze_region(row, query_beds, h3k27me3_dir)
    results.append(r)

# Separate scalar results from array results for the dataframe
scalar_cols = [k for k in results[0].keys() if not k.startswith('_')]
results_df  = pd.DataFrame([{k: r[k] for k in scalar_cols} for r in results])

results_df.to_csv(f'{PROJECT_DIR}/results/20260226_h3k27me3_entropy_lowmse.txt', sep='\t')

high_annot=pd.read_csv(f'{PROJECT_DIR}/results/20260226_h3k27me3_entropy_highmse.txt', sep='\t', index_col=0)
low_annot=pd.read_csv(f'{PROJECT_DIR}/results/20260226_h3k27me3_entropy_lowmse.txt', sep='\t', index_col=0)
high_annot['disruption']='high'
low_annot['disruption']='low'

combined_annot=pd.concat([high_annot, low_annot], axis=0)

combined_annot

combined_annot.columns[6:]

sns.set_style("whitegrid")        

for feat in combined_annot.columns[6:]:
    print(feat)
    groups = combined_annot["disruption"].unique()
    bins = np.linspace(combined_annot[feat].min(), combined_annot[feat].max(), 50)
    #bins = np.linspace(0, 1, 50)
    
    fig, axes = plt.subplots(len(groups), 1, sharex=True, sharey=True, figsize=(6,7))
    
    for ax, g in zip(axes, groups):
        subset = combined_annot[combined_annot["disruption"] == g]
        ax.hist(subset[feat], bins=bins, density=True)
        ax.set_title(g)
    
    plt.xlabel("value")
    plt.suptitle(feat)
    plt.tight_layout()
    plt.show()

alu_dfam=pd.read_csv(f'{DATA_DIR}/alus/data/dfam/alu_dfam_consensus.txt', sep='\t', index_col=0)
alu_dfam

merged.sort_values(by='SVLEN', ascending=False)['SVLEN']

merged.iloc[968130].values

plt.hist(merged[merged['SVLEN']>200]['SVLEN'], bins=80)
plt.show()

merged.columns.values

merged_full=merged_nona[merged_nona['SVLEN']>275]

#plt.scatter(merged['mse_mean'], merged['HFF_H3K27me3_overlap_10kb'])
sns.set_style("whitegrid")        # options: white, dark, whitegrid, darkgrid, ticks

#df_plot=merged_nona
df_plot=merged_full

feat_list=['swScore', 'milliDel', 'milliIns', 'milliDiv',
    'promoter_count_1Mb', 'HFF_CTCF_count_1Mb', 'gene_PC_overlap_1Mb',
            'mapp36_1Mb', 'phyloP100_highfrac_1Mb', 'phyloP100_lowfrac_1Mb', 'phastCon100_1Mb', 'phastCon100_alu',
           'H1ESC_WGBS_count_alu', 'H1ESC_WGBS_count_1Mb', 
           'HFF_H3K9me3_overlap_alu', 'HFF_H3K9me3_overlap_1Mb',
          'HFF_H3K27me3_overlap_alu', 'HFF_H3K27me3_overlap_1Mb']
for feat in feat_list:
    print(feat)
    groups = df_plot["disruption"].unique()
    bins = np.linspace(df_plot[feat].min(), df_plot[feat].max(), 50)
    #bins = np.linspace(0, 1, 50)
    
    fig, axes = plt.subplots(len(groups), 1, sharex=True, sharey=True, figsize=(6,7))
    
    for ax, g in zip(axes, groups):
        subset = df_plot[df_plot["disruption"] == g]
        ax.hist(subset[feat], bins=bins, density=True)
        ax.set_title(g)
    
    plt.xlabel("value")
    plt.suptitle(feat)
    plt.tight_layout()
    plt.show()

# ## CTCF shared vs. non-shared peaks

data_dir=f'{DATA_DIR}/alus/data/'
CTCF_dir=f'{data_dir}encode/CTCF/'

colnames=['chrom', 'start', 'end', 'name', 'score', 'strand', 'signal', 'pval', 'qval', 'peak_summit']
bed_files=[f for f in os.listdir(CTCF_dir) if f.endswith('.bed.gz')]

#peak_len_mean=[]
#peak_len_median=[]
for i in random.sample(list(range(0,len(bed_files))), 15):
    #print(i)

    test1=pd.read_csv(f'{CTCF_dir}{bed_files[i]}', sep='\t', names=colnames)
    test1['peak_len']=test1['end']-test1['start']

    #peak_len_mean.append(np.mean(test1['peak_len']))
    #peak_len_median.append(np.median(test1['peak_len']))

    cols_of_interest = ['peak_len', 'signal', 'qval']
    
    fig, axes = plt.subplots(1, len(cols_of_interest), figsize=(4 * len(cols_of_interest), 4))
    
    for ax, col in zip(axes, cols_of_interest):
        ax.hist(test1[col], bins=30)
        ax.set_xlabel(col)
        ax.set_ylabel('Count')
        ax.set_title(f'{col} Distribution')
    
    plt.tight_layout()
    plt.show()

plt.hist(peak_len_mean, bins=30)
plt.show()

plt.hist(peak_len_median, bins=30)
plt.show()

test1

#params 

MERGE_DIST = 50        # merge peaks within this many bp
ENTROPY_THRESH = 0.5     # below = low entropy
COVERAGE_THRESH = 0.3     # above = constitutive, below = consistently absent

#encode bed files 
colnames=['chrom', 'start', 'end', 'name', 'score', 'strand', 'signal', 'pval', 'qval', 'peak_summit']
SIGNAL_COL = 'signal'

#CTCF file metadata 
metadata=pd.read_csv(f'{CTCF_dir}encode_CTCF_ChIP_primary_cellline_differentiated.tsv', sep='\t', skiprows=1)

#cols to keep from metadata 
cols_to_keep = ['Accession', 'Assay name', 'Assay title','Biosample classification', 'Target', 'Target of assay',
       'Target gene symbol', 'Biosample summary', 'Biosample term name','Description', 'Lab', 'Project', 'Status', 'Files',
       'Biosample accession', 'Biological replicate',
       'Technical replicate', 'Organism', 'Life stage',
       'Biosample age', 'Biosample treatment']  

rows = []
for fname in bed_files:
    prefix = fname.replace('.bed.gz', '') 
    
    # Find rows in metadata where the Files column contains this prefix
    mask = metadata['Files'].str.contains(prefix, na=False)
    matched = metadata[mask]
    
    if not matched.empty:
        row = matched[cols_to_keep].iloc[0].to_dict()
    else:
        row = {col: None for col in cols_to_keep}
    
    row = {'filename': fname, **row}
    rows.append(row)

result_df = pd.DataFrame(rows, columns=['filename'] + cols_to_keep)

#filter result_df
result_df=result_df[result_df['Biological replicate']!='1']
#query_filelist=['K562', 'HCT116', 'GM12878', 'H1', 'foreskin fibroblast', 'neuron', 'cardiac muscle cell']
#query_df=result_df[(result_df['Biosample term name'].isin(query_filelist)) & (result_df['Biological replicate']!='1')]
#filter out if same Biosample term name
query_df=result_df
query_df_unique=query_df.drop_duplicates(subset=['Biosample term name'])
query_beds=query_df_unique['filename'].tolist()
query_names=query_df_unique['Biosample term name'].tolist()

query_names

query_df_unique['Biosample term name'].unique()

MAX_PEAK_SIZE = 1000
MIN_SIGNAL    = 10
MIN_QVAL      = 3       # -log10(FDR) > 5, i.e. FDR < 0.00001
HALF_WIDTH    = 250     # summit-centred window half-width
MERGE_DIST    = -100

all_beds      = [pybedtools.BedTool(os.path.join(CTCF_dir, f)) for f in query_beds]
filtered_beds = []

for bed in all_beds:
    filtered = bed.filter(
        lambda r: len(r)             <= MAX_PEAK_SIZE and
                  float(r.fields[6]) >= MIN_SIGNAL    and   # signalValue
                  float(r.fields[8]) >= MIN_QVAL            # qValue
    ).saveas()

    # re-centre each peak on its summit and fix to HALF_WIDTH on each side
    # narrowPeak field 9 is the summit offset from peak start
    def summit_center(r):
        summit   = int(r.start) + int(r.fields[9])
        r.start  = max(0, summit - HALF_WIDTH)
        r.end    = summit + HALF_WIDTH
        return r

    #filtered = filtered.each(summit_center).saveas()

    if len(filtered) < len(bed):
        print(f'before: {len(bed)}  after: {len(filtered)}')

    filtered_beds.append(filtered)

print("Building merged peak universe...")
cat_bed = pybedtools.BedTool.cat(*filtered_beds, postmerge=False).sort()
universe = cat_bed.merge(d=MERGE_DIST, c=1, o="count")
print(f"Universe: {len(universe)} merged peak regions")

universe_df=universe.to_dataframe()
universe_df['region_len']=universe_df['end']-universe_df['start']

#also get total number of cell line scontributing 
for name, filtered_bed in zip(query_names, filtered_beds):
    intersect = universe.intersect(filtered_bed, c=True, F=0.4)
    universe_df[name] = [int(f.fields[-1]) > 0 for f in intersect]

universe_df['n_cell_lines'] = universe_df[query_names].sum(axis=1)
universe_df['is_shared']    = universe_df['n_cell_lines'] == len(query_names)
universe_df['is_unique']    = universe_df['n_cell_lines'] == 1

print(universe_df['n_cell_lines'].value_counts().sort_index())
print(f"\nshared across all: {universe_df['is_shared'].sum()}")
print(f"unique to one:     {universe_df['is_unique'].sum()}")

#add in categorical labels 
n = len(query_names)
conditions = [
    universe_df['n_cell_lines'] >= 0.80 * n,
    universe_df['n_cell_lines'] >= 0.50 * n,
    universe_df['n_cell_lines'] >= 0.10 * n,
    universe_df['n_cell_lines'] >= 0.02 * n,
    universe_df['n_cell_lines'] == 1,
]
labels = [
    'constitutive',
    'broadly_shared',
    'lineage_restricted',
    'cell_type_specific',
    'singleton',
]

universe_df['peak_class'] = np.select(conditions, labels, default='not_found')
universe_df

universe_df.to_csv(f'{CTCF_dir}CTCF_shared_peaks.txt', sep='\t')

CTCF_dir

plt.hist(universe_df['n_cell_lines'], bins=30)
plt.xlabel('# samples containing merged CTCF peak')
plt.ylabel('count')

plt.hist(universe_df['region_len'], bins=30)
plt.xlabel('# samples containing merged CTCF peak')
plt.ylabel('count')

np.percentile(universe_df['n_cell_lines'], [25, 50, 75])

#read these back in

universe_df.to_csv(f'{CTCF_dir}CTCF_shared_peaks.txt', sep='\t')

#read in alus with overlap scores 
RES_DIR = f'{DATA_DIR}/alus/results/20251125_aluhg38_all_scores/'
OUT_DIR = f'{PROJECT_DIR}/results/20260208_feature_overlap_all/'

alu_annot = pd.read_csv(f'{RES_DIR}20260123_aluhg38_all_annot_GC_mapp_alu_repeats_genetrackoverlap.txt', sep='\t', index_col=0)

os.listdir(OUT_DIR)

CTCF_shared=pd.read_csv(f'{OUT_DIR}aluhg38_feature_count_CTCF_sharedpeaks.csv')
segdup_selfchain=pd.read_csv(f'{OUT_DIR}aluhg38_feature_overlap_segdup_selfchain.csv')
phyloP_1Mb=pd.read_csv(f'{OUT_DIR}aluhg38_phyloP_countBW_1Mb.csv')
phyloP=pd.read_csv(f'{OUT_DIR}aluhg38_phyloP_countBW_short_0_2.csv')

#merge original alu with all annotations 
alu_annot_subset=alu_annot[['CHROM', 'POS', 'END', 'repName', 'orig_idx', 'mse_mean', 'corr_mean', 
                            'REF_start', 'REF_stop', 'var_rel_pos_REF', 
                            'single_gene_name', 'single_gene_type', 'single_gene_overlap_bp',
                           'swScore', 'milliDel', 'id','#bin', 'milliIns', 'milliDiv', 'SVLEN']]
df_list=[alu_annot_subset, CTCF_shared, segdup_selfchain, phyloP, phyloP_1Mb]

merged = reduce(lambda l, r: pd.merge(l, r, on=['CHROM', 'POS', 'END', 'orig_idx'], how="outer"), df_list)

#add in disruption
merged_nona=merged[~merged['mse_mean'].isna()]
q99=np.quantile(merged_nona['mse_mean'], 0.99)
q50=np.quantile(merged_nona['mse_mean'], 0.50)
merged_nona['disruption']=np.where(merged_nona['mse_mean'] >=q99, 'top', 
                            np.where(merged_nona['mse_mean'] <=q50, 'bottom', 'middle'))

merged_nona

#plt.scatter(merged['mse_mean'], merged['HFF_H3K27me3_overlap_10kb'])
sns.set_style("whitegrid")        # options: white, dark, whitegrid, darkgrid, ticks

df_plot=merged_nona
#df_plot=merged_full

feat_list=merged_nona.columns[merged_nona.columns.str.startswith('H1ESC_WGBS')]

for feat in feat_list:
    print(feat)
    groups = df_plot["disruption"].unique()
    bins = np.linspace(df_plot[feat].min(), df_plot[feat].max(), 40)
    #bins = np.linspace(0, 1, 50)
    
    fig, axes = plt.subplots(len(groups), 1, sharex=True, sharey=True, figsize=(3,7))
    
    for ax, g in zip(axes, groups):
        subset = df_plot[df_plot["disruption"] == g]
        ax.hist(subset[feat], bins=bins, density=True)
        ax.set_title(g)
    
    plt.xlabel("value")
    plt.suptitle(feat)
    plt.tight_layout()
    plt.show()

df_plot = merged_nona

# get subfamily from first 4 chars of repName
df_plot = df_plot.copy()
df_plot['subfamily'] = df_plot['repName'].str[:4]

#replace 
SUBFAMILY_RENAME = {
    'FAM': 'Ancestral',
    'FLAM': 'Ancestral',
    'FRAM': 'Ancestral',
}

df_plot['subfamily'] = df_plot['repName'].str[:4].replace(SUBFAMILY_RENAME)

df_plot['subfamily'].value_counts()

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

# # chromHMM states

chromHMM_dir=f'{DATA_DIR}/alus/data/annotations/chromhmm/'
state_dir=f'{chromHMM_dir}mnemonics/'
pd.set_option('display.max_columns', None)

metadata=pd.read_csv(f'{chromHMM_dir}Roadmap_metadata.tsv', sep='\t')
metadata_keep=metadata[metadata['Manual Use Train (Core)']==1]
metadata_keep

test1=pd.read_csv(f'{state_dir}E129_15_coreMarks_hg38lift_mnemonics.bed.gz', sep='\t', 
                  header=None, names=['chrom', 'start', 'end', 'state_name'])
test1['region_len']=test1['end']-test1['start']
test1

##### config 
STATE_DIR  = state_dir      
OUTPUT_DIR = f"{DATA_DIR}/alus/data/annotations/chromhmm/state_universe/"  # <-- change as needed

os.makedirs(OUTPUT_DIR, exist_ok=True)

COL_NAMES  = ['chrom', 'start', 'end', 'state_name']
FILE_GLOB  = os.path.join(STATE_DIR, "*_coreMarks_hg38lift_mnemonics.bed.gz")
all_files   = sorted(glob.glob(FILE_GLOB))

files_want=metadata_keep['Epigenome ID (EID)'].tolist()

#read all sample files 
files_want=metadata_keep['Epigenome ID (EID)'].tolist()
all_files   = [f for f in all_files if os.path.basename(f).split("_")[0] in files_want]

sample_dfs  = {}   # { sample_id : DataFrame }

print(f"Found {len(all_files)} sample files")
for fpath in all_files:
    fname     = os.path.basename(fpath)
    sample_id = fname.split("_")[0]          # e.g. "E129"
    df        = pd.read_csv(fpath, sep='\t', header=None, names=COL_NAMES)
    df['region_len']  = df['end'] - df['start']
    df['sample_id']   = sample_id
    sample_dfs[sample_id] = df
    print(f"  Loaded {sample_id}: {len(df):,} rows")

all_samples = sorted(sample_dfs.keys())      # query_names equivalent
print(f"\nSamples: {all_samples}")

#get unique states across all files
all_states = sorted(
    set().union(*[set(df['state_name'].unique()) for df in sample_dfs.values()])
)
print(f"States found: {all_states}\n")

all_states

OUTPUT_DIR

#process each state 
#for state in all_states:
for state in ['8_ZNF/Rpts','9_Het']:
    print(f"\n{'='*60}")
    print(f"Processing state: {state}")

    # collect per-sample BedTools for this state (only samples that have it)
    state_beds   = []   # pybedtools.BedTool objects
    state_names  = []   # matching sample IDs

    tmp_files = []      # track temp files to clean up

    for sample_id, df in sample_dfs.items():
        state_df = df[df['state_name'] == state][['chrom', 'start', 'end']].copy()
        if state_df.empty:
            continue

        # write temp bed for this sample+state
        tmp_path = os.path.join(OUTPUT_DIR, f"_tmp_{sample_id}_{state}.bed")
        tmp_path = tmp_path.replace(' ', '_')
        tmp_path = tmp_path.replace('/', '_')

        state_df.to_csv(tmp_path, sep='\t', header=False, index=False)
        tmp_files.append(tmp_path)

        state_beds.append(pybedtools.BedTool(tmp_path))
        state_names.append(sample_id)

    if len(state_beds) == 0:
        print(f"  No samples found for state {state}, skipping.")
        continue

    # build universe from state-specific beds
    cat_bed  = pybedtools.BedTool.cat(*state_beds, postmerge=False).sort()

    # first pass: merge with d=0 to get universe, then compute median length
    universe_raw = cat_bed.merge(d=0, c=1, o="count")
    universe_df  = universe_raw.to_dataframe()
    universe_df['region_len'] = universe_df['end'] - universe_df['start']

    median_len = universe_df['region_len'].median()
    MERGE_DIST = -int(median_len / 2)
    print(f"  Median region len (d=0 universe): {median_len:.1f} bp  →  MERGE_DIST: {MERGE_DIST}")

    # second pass: rebuild universe with state-specific merge distance
    universe = cat_bed.merge(d=MERGE_DIST, c=1, o="count")
    print(f"  Universe: {len(universe)} merged peak regions")

    universe_df = universe.to_dataframe()
    universe_df['region_len'] = universe_df['end'] - universe_df['start']

    # intersect each sample against the universe
    for name, bed in zip(state_names, state_beds):
        intersect = universe.intersect(bed, c=True, F=0.4)
        universe_df[name] = [int(f.fields[-1]) > 0 for f in intersect]

    # summary columns
    universe_df['n_cell_lines'] = universe_df[state_names].sum(axis=1)
    universe_df['is_shared']    = universe_df['n_cell_lines'] == len(state_names)
    universe_df['is_unique']    = universe_df['n_cell_lines'] == 1

    print(universe_df['n_cell_lines'].value_counts().sort_index())
    print(f"  Shared across all: {universe_df['is_shared'].sum()}")
    print(f"  Unique to one:     {universe_df['is_unique'].sum()}")

    # classify peaks
    universe_df['frac'] = universe_df['n_cell_lines'] / len(state_names)
    
    conditions = [
        universe_df['n_cell_lines'] == 1,                                       
        universe_df['frac'] <= 0.10,                                          
        universe_df['frac'] <= 0.50,                                      
        universe_df['frac'] <= 0.80,                                        
        universe_df['frac'] <= 1.00,                                          
    ]

    labels = [
        'singleton',          # exactly 1 sample
        'cell_type_specific', # >1 sample but ≤25%
        'lineage_restricted', # 25–50%
        'broadly_shared',     # 50–75%
        'constitutive',       # 75–100%
    ]

    #np.select fills out first one that matches that criteria, and onwards
    universe_df['peak_class'] = np.select(conditions, labels, default='not_found')
    universe_df['state_name'] = state
    universe_df['n_samples_total'] = len(state_names)

    # save
    safe_state = str(state).replace('/', '_').replace(' ', '_')
    out_path   = os.path.join(OUTPUT_DIR, f"universe_{safe_state}.csv")
    universe_df.to_csv(out_path, index=False)
    print(f"  Saved → {out_path}")

    # clean up temp bed files
    for tmp in tmp_files:
        os.remove(tmp)
    pybedtools.cleanup()

print("\nDone.")

test1['state_name'].value_counts()

#add in categorical labels 
n = len(query_names)
conditions = [
    universe_df['n_cell_lines'] >= 0.80 * n,
    universe_df['n_cell_lines'] >= 0.50 * n,
    universe_df['n_cell_lines'] >= 0.10 * n,
    universe_df['n_cell_lines'] >= 0.02 * n,
    universe_df['n_cell_lines'] == 1,
]
labels = [
    'constitutive',
    'broadly_shared',
    'lineage_restricted',
    'cell_type_specific',
    'singleton',
]

universe_df['peak_class'] = np.select(conditions, labels, default='not_found')
universe_df

# ## read in results

#read in merged
merged_path=f'{PROJECT_DIR}/results/paper_results/20260316_merged_aluhg38_all_featureannot_NA_blacklist_filter_region_alu.txt'
merged=pd.read_csv(merged_path, sep='\t')

#add in disruption
merged_nona=merged[~merged['mse_mean'].isna()]
q99=np.quantile(merged_nona['mse_mean'], 0.99)
q50=np.quantile(merged_nona['mse_mean'], 0.50)
merged_nona['disruption']=np.where(merged_nona['mse_mean'] >=q99, 'top', 
                            np.where(merged_nona['mse_mean'] <=q50, 'bottom', 'middle'))

#read in states
chromHMM_shared=pd.read_csv(f'{chromHMM_dir}overlap/aluhg38_feature_count_10_TssBiv.csv')

chromHMM_shared

merged_nona.columns.values

alu_annot_subset=merged_nona[['CHROM', 'POS', 'END', 'repName', 'orig_idx', 'mse_mean', 'corr_mean', 'disruption',
                                'single_gene_name', 'single_gene_type', 'single_gene_overlap_bp', 
                              'phyloP100_highfrac_1Mb', 'phyloP100_lowfrac_1Mb', 
                              'H1ESC_WGBS_count_1Mb', 'H1ESC_WGBS_count_100kb', 'H1ESC_WGBS_count_10kb', 'H1ESC_WGBS_count_1kb',
                                'HFF_H3K9me3_overlap_1kb', 'HFF_H3K9me3_overlap_10kb',
                               'HFF_H3K9me3_overlap_100kb', 'HFF_H3K9me3_overlap_1Mb',
                               'HFF_H3K27me3_overlap_alu', 'HFF_H3K27me3_overlap_1kb',
                               'HFF_H3K27me3_overlap_10kb', 'HFF_H3K27me3_overlap_100kb']]

alu_annot_subset

plt.hist(alu_annot_subset['phyloP100_highfrac_1Mb'], 100)
plt.show()

for state_file in os.listdir(f'{chromHMM_dir}overlap/'):
    print(state_file)

    if (state_file=='aluhg38_feature_count_all_states.csv') or (state_file=='count'):
        continue

    if (state_file!='aluhg38_feature_count_8_ZNF_Rpts.csv'):
        continue

    state_name = state_file.split('aluhg38_feature_count_')[1].replace('.csv', '')

    chromHMM_shared=pd.read_csv(f'{chromHMM_dir}overlap/{state_file}')

    #merge original alu with all annotations 

    df_list=[alu_annot_subset, chromHMM_shared]
    
    merged_state = reduce(lambda l, r: pd.merge(l, r, on=['CHROM', 'POS', 'END', 'orig_idx'], how="left"), df_list)

    #plot distributions 
    sns.set_style("whitegrid")        # options: white, dark, whitegrid, darkgrid, ticks
    
    df_plot=merged_state
    
    feat_list=df_plot.columns[df_plot.columns.str.startswith(state_name)]
    feat_list_1mb=feat_list[feat_list.str.endswith('1Mb')]
    for feat in feat_list:
        print(feat)
        #groups = df_plot["disruption"].unique()
        groups=['top', 'bottom']
        bins = np.linspace(df_plot[feat].min(), df_plot[feat].max(), 30)
        #bins = np.linspace(0, 1, 50)
        
        fig, axes = plt.subplots(len(groups), 1, sharex=True, sharey=True, figsize=(3,5))
        
        for ax, g in zip(axes, groups):
            subset = df_plot[df_plot["disruption"] == g]
            ax.hist(subset[feat], bins=bins, density=True)
            ax.set_title(g)
        
        plt.xlabel("value")
        plt.suptitle(feat)
        plt.tight_layout()
        plt.show()

# ### plot chromHMM states vs phyloP categories

chromHMM_dir

for state_file in os.listdir(f'{chromHMM_dir}overlap/'):
    print(state_file)
    if (state_file=='aluhg38_feature_count_all_states.csv') or (state_file=='count'):
        continue

    state_name = state_file.split('aluhg38_feature_count_')[1].replace('.csv', '')
    chromHMM_shared=pd.read_csv(f'{chromHMM_dir}overlap/{state_file}')

    df_list=[alu_annot_subset, chromHMM_shared]
    merged_state = reduce(lambda l, r: pd.merge(l, r, on=['CHROM', 'POS', 'END', 'orig_idx'], how="left"), df_list)
    
    sns.set_style("whitegrid")
    df_plot = merged_state
    
    feat_list = df_plot.columns[df_plot.columns.str.startswith(state_name)]
    feat_list_1mb=feat_list[feat_list.str.endswith('1Mb')]

    phylop_col='phyloP100_highfrac_1Mb'
    q10=np.quantile(alu_annot_subset[phylop_col], 0.10)
    q90=np.quantile(alu_annot_subset[phylop_col], 0.90)
    df_plot['phylop_cat'] = np.where(df_plot[phylop_col] >= q90, 'high', 
                                     np.where(df_plot[phylop_col] <= q10, 'low', 'mid'))

    phylop_vals = ['high', 'low']   # phyloP_high = 0/1 or False/True
    groups = ['top', 'bottom']

    for feat in feat_list_1mb:
        print(feat)
        bins = np.linspace(df_plot[feat].min(), df_plot[feat].max(), 30)

        # rows = disruption groups, cols = phyloP_high
        fig, axes = plt.subplots(
            len(groups), len(phylop_vals),
            sharex=True, sharey=True,
            figsize=(6, 3)
        )

        for col_idx, phylop in enumerate(phylop_vals):
            phylop_label = f'phyloP_high={(phylop)}'

            df_phylop = df_plot[df_plot['phylop_cat'] == phylop]

            for row_idx, g in enumerate(groups):
                ax = axes[row_idx, col_idx]
                subset = df_phylop[df_phylop['disruption'] == g]
                ax.hist(subset[feat], bins=bins, density=True)

                # row labels (disruption) on left-most column only
                if col_idx == 0:
                    ax.set_ylabel(g)
                # column labels (phyloP) on top row only
                if row_idx == 0:
                    ax.set_title(phylop_label)

        plt.suptitle(feat)
        fig.text(0.5, 0.01, 'value', ha='center')
        plt.tight_layout()
        plt.show()

os.listdir(f'{chromHMM_dir}overlap/')

state_file_interest=['aluhg38_feature_count_9_Het.csv',
                     'aluhg38_feature_count_7_Enh.csv',
                     'aluhg38_feature_count_6_EnhG.csv',
                    'aluhg38_feature_count_10_TssBiv.csv',
                    'aluhg38_feature_count_11_BivFlnk.csv',
                    'aluhg38_feature_count_12_EnhBiv.csv',
                     'aluhg38_feature_count_13_ReprPC.csv',
                     'aluhg38_feature_count_14_ReprPCWk.csv']

# ### plot chrom state with methylation

#for state_file in os.listdir(f'{chromHMM_dir}overlap/'):
for state_file in state_file_interest:
    print(state_file)
    if (state_file=='aluhg38_feature_count_all_states.csv') or (state_file=='count'):
        continue

    state_name = state_file.split('aluhg38_feature_count_')[1].replace('.csv', '')
    chromHMM_shared=pd.read_csv(f'{chromHMM_dir}overlap/{state_file}')

    df_list=[alu_annot_subset, chromHMM_shared]
    merged_state = reduce(lambda l, r: pd.merge(l, r, on=['CHROM', 'POS', 'END', 'orig_idx'], how="left"), df_list)
    
    sns.set_style("whitegrid")
    df_plot = merged_state
    
    feat_list = df_plot.columns[df_plot.columns.str.startswith(state_name)]
    feat_list_1mb=feat_list[feat_list.str.endswith('1Mb')]

    #HFF_H3K27me3_overlap_1kb, H1ESC_WGBS_count_10kb, HFF_H3K9me3_overlap_1kb
    phylop_col='HFF_H3K27me3_overlap_100kb'
    q10=np.quantile(alu_annot_subset[phylop_col], 0.10)
    q90=np.quantile(alu_annot_subset[phylop_col], 0.50)
    df_plot['methyl_cat'] = np.where(df_plot[phylop_col] >= q90, 'high', 
                                     np.where(df_plot[phylop_col] <= q10, 'low', 'mid'))

    phylop_vals = ['high', 'low']   # phyloP_high = 0/1 or False/True
    groups = ['top', 'bottom']

    for feat in feat_list_1mb:
        print(feat)
        bins = np.linspace(df_plot[feat].min(), df_plot[feat].max(), 30)

        # rows = disruption groups, cols = phyloP_high
        fig, axes = plt.subplots(
            len(groups), len(phylop_vals),
            sharex=True, sharey=True,
            figsize=(6, 3)
        )

        for col_idx, phylop in enumerate(phylop_vals):
            phylop_label = f'{phylop_col}={(phylop)}'

            df_phylop = df_plot[df_plot['methyl_cat'] == phylop]

            for row_idx, g in enumerate(groups):
                ax = axes[row_idx, col_idx]
                subset = df_phylop[df_phylop['disruption'] == g]
                ax.hist(subset[feat], bins=bins, density=True)

                # row labels (disruption) on left-most column only
                if col_idx == 0:
                    ax.set_ylabel(g)
                # column labels (phyloP) on top row only
                if row_idx == 0:
                    ax.set_title(phylop_label)

        plt.suptitle(feat)
        fig.text(0.5, 0.01, 'value', ha='center')
        plt.tight_layout()
        plt.show()

#for state_file in os.listdir(f'{chromHMM_dir}overlap/'):
for state_file in state_file_interest:
    print(state_file)
    if (state_file=='aluhg38_feature_count_all_states.csv') or (state_file=='count'):
        continue
    state_name = state_file.split('aluhg38_feature_count_')[1].replace('.csv', '')
    chromHMM_shared=pd.read_csv(f'{chromHMM_dir}overlap/{state_file}')
    
    df_list=[alu_annot_subset, chromHMM_shared]
    merged_state = reduce(lambda l, r: pd.merge(l, r, on=['CHROM', 'POS', 'END', 'orig_idx'], how="left"), df_list)
    
    sns.set_style("whitegrid")
    df_plot = merged_state
    
    feat_list = df_plot.columns[df_plot.columns.str.startswith(state_name)]
    feat_list_1mb=feat_list[feat_list.str.endswith('100kb')]
    methylation_col='H1ESC_WGBS_count_100kb'  # adjust to your actual methylation column name
    groups = ['top', 'bottom']
    for feat in feat_list_1mb:
        print(feat)
        fig, axes = plt.subplots(
            len(groups), 1,
            sharex=True, sharey=True,
            figsize=(5, 6)
        )
        # shared extents so both panels are directly comparable
        sub_all = df_plot.dropna(subset=[feat, methylation_col])
        extent = (
            sub_all[feat].min(), sub_all[feat].max(),
            sub_all[methylation_col].min(), sub_all[methylation_col].max()
        )
        for row_idx, g in enumerate(groups):
            ax = axes[row_idx]
            subset = df_plot[df_plot['disruption'] == g].dropna(subset=[feat, methylation_col])
            if len(subset) > 0:
                hb = ax.hexbin(
                    subset[feat], subset[methylation_col],
                    gridsize=30, cmap='viridis', mincnt=1, extent=extent
                )
                fig.colorbar(hb, ax=ax, label='count')
            ax.set_ylabel(f'{g}\nmethylation')
        axes[-1].set_xlabel(feat)
        plt.suptitle(feat)
        plt.tight_layout()
        plt.show()

#read in methylation from different cell types -- see if this changes at all

methyl_overlap=pd.read_csv(f'{PROJECT_DIR}/results/paper_results/aluhg38_WGBS_count.csv')
methyl_overlap

alu_annot_methyl=alu_annot_subset.merge(methyl_overlap, on=['CHROM', 'POS', 'END', 'orig_idx'])

# Pull all ENC methylation columns ending with '1Mb'
enc_1mb_cols = [c for c in df_plot.columns if c.endswith('1Mb') and ('ENC' in c)]

corr_matrix = alu_annot_methyl[enc_1mb_cols].corr()

# Visualize
fig, ax = plt.subplots(figsize=(max(6, 0.5 * len(enc_1mb_cols)), max(5, 0.5 * len(enc_1mb_cols))))
sns.heatmap(corr_matrix,annot=True, fmt='.2f',
    cmap='coolwarm', center=0,vmin=-1, vmax=1,
    square=True,cbar_kws={'label': 'Pearson r'},ax=ax,)
plt.tight_layout()
plt.show()

#for state_file in os.listdir(f'{chromHMM_dir}overlap/')[0:2]:
for state_file in state_file_interest:

    print(state_file)
    if (state_file == 'aluhg38_feature_count_all_states.csv') or (state_file == 'count'):
        continue
    state_name = state_file.split('aluhg38_feature_count_')[1].replace('.csv', '')
    chromHMM_shared = pd.read_csv(f'{chromHMM_dir}overlap/{state_file}')
    
    df_list = [alu_annot_methyl, chromHMM_shared]
    merged_state = reduce(lambda l, r: pd.merge(l, r, on=['CHROM', 'POS', 'END', 'orig_idx'], how="left"), df_list)
    
    sns.set_style("whitegrid")
    df_plot = merged_state
    
    # chromHMM features for this state at 100kb
    feat_list = df_plot.columns[df_plot.columns.str.startswith(state_name)]
    feat_list_1mb = feat_list[feat_list.str.endswith('1Mb')]
    
    # methylation cols: end with '100kb' AND contain 'WGBS_ENC'
    methylation_cols = [c for c in df_plot.columns if c.endswith('1Mb') and ('ENC' in c)]
    
    groups = ['top', 'bottom']
    
    # ---- 1. Side-by-side histogram comparison across all methylation cols ----
    if len(methylation_cols) > 0:
        fig, axes = plt.subplots(
            1, len(methylation_cols),
            figsize=(4 * len(methylation_cols), 4),
            sharey=False,
        )

        for ax, methylation_col in zip(axes, methylation_cols):
            sub_all = df_plot.dropna(subset=[methylation_col])
            if len(sub_all) == 0:
                continue
            bins = np.linspace(sub_all[methylation_col].min(), sub_all[methylation_col].max(), 40)
            
            for g in groups:
                subset = df_plot[df_plot['disruption'] == g].dropna(subset=[methylation_col])
                if len(subset) > 0:
                    ax.hist(
                        subset[methylation_col],
                        bins=bins,
                        alpha=0.5,
                        density=True,
                        label=f'{g} (n={len(subset)})',
                        edgecolor='black',
                        linewidth=0.3,
                    )
            ax.set_xlabel(methylation_col, fontsize=8)
            ax.set_ylabel('density')
            ax.legend(fontsize=8)
            ax.set_title(methylation_col, fontsize=9)
        
        plt.suptitle(f'{state_name} — methylation distributions', y=1.02)
        plt.tight_layout()
        plt.show()
    
    # ---- 2. Hexbin: chromHMM feature vs methylation, top/bottom rows ----
    for feat in feat_list_1mb:
        for methylation_col in methylation_cols:
            print(f'{feat}  vs  {methylation_col}')
            
            fig, axes = plt.subplots(1, len(groups), sharex=True, sharey=True,figsize=(8, 4),)
            
            sub_all = df_plot.dropna(subset=[feat, methylation_col])
            if len(sub_all) == 0:
                plt.close(fig)
                continue
            extent = (
                sub_all[feat].min(), sub_all[feat].max(),
                sub_all[methylation_col].min(), sub_all[methylation_col].max(),
            )
            
            for row_idx, g in enumerate(groups):
                ax = axes[row_idx]
                subset = df_plot[df_plot['disruption'] == g].dropna(subset=[feat, methylation_col])
                if len(subset) > 0:
                    hb = ax.hexbin(
                        subset[feat], subset[methylation_col],
                        gridsize=30, cmap='viridis', mincnt=1, extent=extent,
                    )
                    fig.colorbar(hb, ax=ax, label='count')
                ax.set_ylabel(f'{g}\n{methylation_col}', fontsize=8)
            
            axes[-1].set_xlabel(feat)
            plt.suptitle(f'{feat}  vs  {methylation_col}', fontsize=10)
            plt.tight_layout()
            plt.show()
