"""Converted from 20250519_Alu_sequencestructure_windowmutate.ipynb; logic unchanged."""

from collections import defaultdict
import math
import os
import subprocess

from Bio import pairwise2
from Bio.Seq import Seq
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd 
import pysam
import seaborn as sns

# Paths - edit for your environment
SUPREMO_DIR = os.environ.get("SUPREMO_DIR", "/pollard/home/szhang20/akita_variant_scoring")
AKITA_DIR = os.environ.get("AKITA_DIR", "/pollard/home/szhang20/akita")
PROJECT_DIR = os.environ.get("ALU_PROJECT_DIR", "/pollard/home/szhang20/alu")
DATA_DIR = os.environ.get("ALU_DATA_DIR", "/pollard/data/projects/shzhang")
POLLARD_DATA = os.environ.get("POLLARD_DATA", "/pollard/data")

res_dir=f'{PROJECT_DIR}/results/'

# res_dir=f'{PROJECT_DIR}/results/'
# DEL_MSE=pd.read_csv(f'{res_dir}202520513_hg38Alus_100sample_repName_DEL_MSE_topbottom.txt', sep='\t', index_col=0)
# DEL_top=DEL_MSE[DEL_MSE.disruption=='top']
# DEL_top

#read in alu_scores and combine 
#res_dir_all=f'{DATA_DIR}/alus/results/20251125_aluhg38_all_scores/'
#alu_all=pd.read_csv(f'{res_dir_all}20260109_aluhg38_all_annot_GC_mapp_alu_repeats_geneoverlap.txt', sep='\t', index_col=0)

alu_all=pd.read_csv(f'{res_dir}paper_results/20260316_merged_aluhg38_all_featureannot_NA_blacklist_filter_region_alu.txt', sep='\t', index_col=0)
alu_all

#subset to top 1% MSE, drop na
alu_all_noNA = alu_all.dropna(subset=['mse_mean'])
top1_quantile= alu_all_noNA['mse_mean'].quantile(0.99)
top1_mse = alu_all_noNA[alu_all_noNA['mse_mean'] >= top1_quantile]
top1_mse

alu_consensus=pd.read_csv(f'{res_dir}20240129/20240129_alu_consensus_seq_data.txt', sep='\t', index_col=0)
alu_consensus.head()

alu_top100_repName= (alu_all_noNA.sort_values("mse_mean", ascending=False).groupby("repName", group_keys=False).head(100))

window=30
stride=10

# Initialize output list
output = []

# Process each region
for idx, row in alu_top100_repName.iterrows():
    chrom = row["CHROM"]
    start = row["POS"]
    end = row["END"]
    orig_idx=row["orig_idx"]
    
    for window_start in range(start, end - window + 1, stride):
        window_end = window_start + window
        output.append({
            "CHROM": chrom,
            "POS": window_start,
            "REF": '-',
            "ALT": '-',
            "END": window_end,
            "SVTYPE": 'DEL',
            "SVLEN": window_end-window_start,
            "orig_idx": orig_idx
        })

# Convert to DataFrame
windows_df = pd.DataFrame(output)

windows_df

#split for run

out_dir=f'{PROJECT_DIR}/results/paper_results/20260413_alu_rollingwindow_top100/'
os.mkdir(out_dir)
os.mkdir(f'{out_dir}input/')
os.mkdir(f'{out_dir}output/')
os.mkdir(f'{out_dir}logs/')

chunk_size=10000

for i, start in enumerate(range(0, len(windows_df), chunk_size)):
#for i, start in enumerate(range(0, 200, chunk_size)):
    print(i, start)

    chunk = windows_df[start : start + chunk_size]
    chunk.to_csv(f"{out_dir}input/alu_top100_repName_rolling_{i}.txt", sep = '\t',index=False)

windows_df.to_csv(f'{out_dir}/top100_repName_window30_shift10.txt', sep='\t', index=False, header=True)

# #### repeat with 5100 sampled Alus

alu5100=pd.read_csv(f'{res_dir}paper_results/20260413_alu_sample/alu5100_annot.txt', sep='\t', index_col=0)
alu5100

window=30
stride=10

# Initialize output list
output = []

# Process each region
for idx, row in alu5100.iterrows():
    chrom = row["CHROM"]
    start = row["POS"]
    end = row["END"]
    orig_idx=row["orig_idx"]
    
    for window_start in range(start, end - window + 1, stride):
        window_end = window_start + window
        output.append({
            "CHROM": chrom,
            "POS": window_start,
            "REF": '-',
            "ALT": '-',
            "END": window_end,
            "SVTYPE": 'DEL',
            "SVLEN": window_end-window_start,
            "orig_idx": orig_idx
        })

# Convert to DataFrame
windows_df = pd.DataFrame(output)

#split for run

out_dir=f'{PROJECT_DIR}/results/paper_results/20260413_alu_rollingwindow_5100/'
os.mkdir(out_dir)
os.mkdir(f'{out_dir}input/')
os.mkdir(f'{out_dir}output/')
os.mkdir(f'{out_dir}logs/')

chunk_size=10000

for i, start in enumerate(range(0, len(windows_df), chunk_size)):
#for i, start in enumerate(range(0, 200, chunk_size)):
    print(i, start)

    chunk = windows_df[start : start + chunk_size]
    chunk.to_csv(f"{out_dir}input/alu_5100_repName_rolling_{i}.txt", sep = '\t',index=False)

windows_df.to_csv(f'{out_dir}alu_5100_repName_window30_shift10.txt', sep='\t', index=False, header=True)

# ### results

#read the windows back in 
score_df_list=[]
for i in range(0, 12):
    df_score=pd.read_csv(f'{out_dir}output/alu_top100_repName_rolling_{i}_scores', sep='\t')
    df_alu=pd.read_csv(f'{out_dir}input/alu_top100_repName_rolling_{i}.txt', sep='\t')

    #merge on indices 
    alu_score_i = pd.merge(df_alu, df_score, left_index=True, right_index=True, how='inner')

    assert len(alu_score_i) == len(df_score)
    assert len(alu_score_i) == len(df_alu)

    score_df_list.append(alu_score_i)
    
window_score_df=pd.concat(score_df_list)
window_score_df=window_score_df.reset_index(drop=True)

#calculate mean, median, sd for scores per row 
mse_cols = [col for col in window_score_df.columns if col.startswith('mse')]
corr_cols = [col for col in window_score_df.columns if col.startswith('corr')]
df_mse = window_score_df[mse_cols]
df_corr = window_score_df[corr_cols]

# calcualte and add to alu scores 
window_score_df['mse_mean'] = df_mse.mean(axis=1)
window_score_df['mse_median'] = df_mse.median(axis=1)
window_score_df['mse_std'] = df_mse.std(axis=1)

window_score_df['corr_mean'] = df_corr.mean(axis=1)
window_score_df['corr_median'] = df_corr.median(axis=1)
window_score_df['corr_std'] = df_corr.std(axis=1)

window_score_df=window_score_df.merge(alu_top100_repName[['orig_idx', 'strand']], how='left')
window_score_df

window_score_df.columns

#this is old interation
# scores=pd.read_csv(f'{res_dir}20250524_rolling_window/top_DEL_MSE_window30_shift10_scores', sep='\t')
# windows_df=pd.read_csv(f'{res_dir}20250524_rolling_window/top_DEL_MSE_window30_shift10.txt', sep='\t')
# windows_df['row_index']=windows_df.index
# windows_scores=windows_df.merge(scores, left_on='row_index', right_on='var_index', how='left')
# #windows_scores=windows_df.merge(DEL_top, left_on='original_index', right_on='row_idx')
# windows_scores

ref_fasta = pysam.FastaFile(f"{SUPREMO_DIR}/data/hg38.fa")

top1_quantile

#subset_high_MSE=alu_top100_repName[alu_top100_repName['mse_mean']>top1_quantile]
subset_high_MSE=alu_top100_repName

subset_high_MSE=subset_high_MSE[subset_high_MSE['strand']=='+']
subset_high_MSE

mapped_scores = []
windows_scores=window_score_df

for idx, row in subset_high_MSE.iterrows():
    chrom = row["CHROM"]
    start = int(row["POS"])
    end = int(row["END"])
    subfam = row["repName"]
    original_idx=row['orig_idx']
    
    # Pull sequence from hg38 ref 
    element_seq = ref_fasta.fetch(chrom, start, end).upper().replace("N", "")
    
    #check if consensus seq exists 
    consensus_row = alu_consensus[alu_consensus['Alu_name'] == subfam]

    if consensus_row.empty:
        print(f"Warning: no consensus found for subfamily {subfam}")
        continue

    consensus_seq = consensus_row['Alu_seq'].values[0]

    if consensus_seq is None or not element_seq:
        continue

    # Align element to its subfamily consensus
    alignments = pairwise2.align.globalms(
        consensus_seq, element_seq,
        2, -1, -5, -0.5,
        one_alignment_only=True
    )

    if not alignments:
        continue

    consensus_aln, element_aln, *_ = alignments[0]

    # Map: element-relative position → consensus position
    e_pos = c_pos = 0
    element_to_consensus = {}
    for e_char, c_char in zip(element_aln, consensus_aln):
        if e_char != '-':
            if c_char != '-':
                element_to_consensus[e_pos] = c_pos
            e_pos += 1
        if c_char != '-':
            c_pos += 1

    # Get windows for this element
    element_windows = windows_scores[windows_scores["orig_idx"] == original_idx]

    for _, win in element_windows.iterrows():
        win_start = int(win["POS"]) - start
        win_end = int(win["END"]) - start

        consensus_start = element_to_consensus.get(win_start)
        consensus_end = element_to_consensus.get(win_end - 1)  # inclusive end

        if consensus_start is not None and consensus_end is not None:
            mapped_scores.append({
                "subfamily": subfam,
                "element_index": idx,
                "consensus_start": consensus_start,
                "consensus_end": consensus_end,
                "score": win["mse_mean"]
            })
            
mapped_df = pd.DataFrame(mapped_scores)

mapped_df

coverage_dict={}
for subfam, seq in zip(alu_consensus["Alu_name"], alu_consensus["Alu_seq"]):
    length = len(seq)
    coverage_dict[subfam] = {
        "sum": np.zeros(length),
        "count": np.zeros(length)
    }

for _, row in mapped_df.iterrows():
    subfam = row["subfamily"]
    start = int(row["consensus_start"])
    end = int(row["consensus_end"])
    score = row["score"]

    if subfam not in coverage_dict:
        continue
    
    consensus_len = len(coverage_dict[subfam]["sum"])
    if start >= consensus_len:
        continue

    end = min(end, consensus_len)
    coverage_dict[subfam]["sum"][start:end] += score
    coverage_dict[subfam]["count"][start:end] += 1

coverage_profiles = []

for subfam, arrays in coverage_dict.items():
    length = len(arrays["sum"])
    avg = np.divide(
        arrays["sum"], arrays["count"],
        out=np.zeros(length),
        where=arrays["count"] > 0
    )

    df = pd.DataFrame({
        "position": np.arange(length),
        "score": avg,
        "subfamily": subfam
    })
    coverage_profiles.append(df)

plot_df = pd.concat(coverage_profiles, ignore_index=True)

plot_df

# Get list of unique subfamilies
subfamilies = sorted(plot_df["subfamily"].unique())
n_subfams = len(subfamilies)

# Set up subplot grid
n_cols = 4  
n_rows = math.ceil(n_subfams / n_cols)

fig, axes = plt.subplots(n_rows, n_cols, figsize=(5 * n_cols, 3.5 * n_rows), sharey=True)
axes = axes.flatten()

# Plot each subfamily separately
for i, subfam in enumerate(subfamilies):
    ax = axes[i]
    sub_df = plot_df[plot_df["subfamily"] == subfam]

    sns.lineplot(data=sub_df, x="position", y="score", ax=ax)
    ax.set_title(subfam)
    ax.set_xlabel("Consensus position")
    ax.set_ylabel("Score")

# Turn off unused axes
for j in range(i + 1, len(axes)):
    axes[j].axis("off")

plt.tight_layout()
plt.savefig(f"{PROJECT_DIR}/Alu_rolling_window.pdf", dpi=300, bbox_inches="tight", transparent=True)
plt.suptitle("Consensus-Aligned Coverage by Alu Subfamily", fontsize=16, y=1.02)
plt.show()

len(alu_consensus[alu_consensus['Alu_name']=='AluYh9']['Alu_seq'].values[0])

len(alu_consensus[alu_consensus['Alu_name']=='AluYm1']['Alu_seq'].values[0])

alu_consensus[alu_consensus['Alu_name']=='AluYh9']['Alu_seq'].values

alu_consensus[alu_consensus['Alu_name']=='AluYm1']['Alu_seq'].values

plt.figure(figsize=(12, 6))
sns.lineplot(data=plot_df, x="position", y="score", hue="subfamily")
plt.title("Consensus-Aligned Coverage Track")
plt.xlabel("Position along consensus")
plt.ylabel("Average window score")
#plt.tight_layout()
plt.show()

def build_alignment_map(seq1, seq2):
    """Returns a mapping from seq1 indices to seq2 indices via global alignment."""
    alignment = pairwise2.align.globalms(seq1, seq2, 2, -1, -0.5, -0.1)[0]
    aligned1, aligned2 = alignment.seqA, alignment.seqB
    idx1, idx2 = 0, 0
    mapping = {}
    for i in range(len(aligned1)):
        if aligned1[i] != '-' and aligned2[i] != '-':
            mapping[idx1] = idx2
        if aligned1[i] != '-':
            idx1 += 1
        if aligned2[i] != '-':
            idx2 += 1
    return mapping

# Example setup
k = 30  # or any window size
stride = 1  # you can change this for faster scan
master_len = len(master_consensus)
master_disruption = [[] for _ in range(master_len)]

# Precompute sub -> master mappings
sub_to_master_map = {
    subfam: build_alignment_map(sub_con, master_consensus)
    for subfam, sub_con in sub_consensuses.items()
}

for seq, subfam in zip(sequences, subfamily_labels):
    sub_con = sub_consensuses[subfam]
    seq_to_sub_map = build_alignment_map(seq, sub_con)
    sub_to_master = sub_to_master_map[subfam]

    # Slide window of size k
    for i in range(0, len(seq) - k + 1, stride):
        mutated_seq = mutate_window(seq, i, k)
        disruption = model_score(seq) - model_score(mutated_seq)

        # Map window positions to master
        master_positions = []
        for j in range(i, i + k):
            if j in seq_to_sub_map:
                sub_pos = seq_to_sub_map[j]
                if sub_pos in sub_to_master:
                    master_pos = sub_to_master[sub_pos]
                    if 0 <= master_pos < master_len:
                        master_positions.append(master_pos)

        for pos in master_positions:
            master_disruption[pos].append(disruption)

# Aggregate (mean)
mean_disruption = [np.mean(scores) if scores else 0 for scores in master_disruption]

# Plot
plt.plot(mean_disruption)
plt.title(f"Mean Disruption Score (Window Size = {k}) across Master Consensus")
plt.xlabel("Master Consensus Position (bp)")
plt.ylabel("Mean Disruption")
plt.tight_layout()
plt.show()
