"""Converted from 20260810_alu_random_gene_overlap.ipynb; logic unchanged."""

import sys
import os
import re
import subprocess
import time

from pybedtools import BedTool
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pybedtools
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "utils"))
from io_helpers import parse_attributes

# Paths - edit for your environment
SUPREMO_DIR = os.environ.get("SUPREMO_DIR", "/pollard/home/szhang20/akita_variant_scoring")
AKITA_DIR = os.environ.get("AKITA_DIR", "/pollard/home/szhang20/akita")
PROJECT_DIR = os.environ.get("ALU_PROJECT_DIR", "/pollard/home/szhang20/alu")
DATA_DIR = os.environ.get("ALU_DATA_DIR", "/pollard/data/projects/shzhang")
POLLARD_DATA = os.environ.get("POLLARD_DATA", "/pollard/data")

# # Random expectation for Alu-gene overlap
# 
# Question: if the scored Alus were distributed randomly throughout the genome
# (instead of at their real, evolutionarily-selected insertion sites), what
# proportion would we expect to fall into a gene?
# 
# Approach:
# 1. Load the scored/blacklist-filtered Alu set (`20260316_merged_aluhg38_all_featureannot_NA_blacklist_filter_region_alu.txt`)
#    and compute the **observed** proportion overlapping a gene body.
# 2. Shuffle those same Alu intervals (same chromosome-wide pool, same lengths, same count)
#    to random genomic positions, excluding ENCODE blacklist regions, using `bedtools shuffle`.
# 3. Repeat many times to build a **null distribution** of "proportion in a gene" under
#    random placement, and compare it to the observed value.
# 
# "Falls into a gene" is defined the same way the `single_gene_name` column in the
# scored file was originally annotated: an Alu is assigned to a gene if >=50% of the
# Alu's length is covered by a single GENCODE `gene` feature (any `gene_type`, not just
# protein-coding).

# ## Paths

data_dir = f'{DATA_DIR}/alus/data/'
res_dir = f'{PROJECT_DIR}/results/'

ALU_SCORED_PATH = f'{res_dir}paper_results/20260316_merged_aluhg38_all_featureannot_NA_blacklist_filter_region_alu.txt'
GENCODE_PATH = f'{data_dir}gencode/gencode.v29.annotation.gtf.gz'
BLACKLIST_PATH = f'{data_dir}encBlacklist.bed'
CHROM_LENGTHS_PATH = f'{SUPREMO_DIR}/data/chrom_lengths_hg38'

OUT_DIR = f'{res_dir}paper_results/20260810_random_gene_overlap/'
os.makedirs(OUT_DIR, exist_ok=True)

N_SHUFFLES = 100
SEED = 729
OVERLAP_FRAC = 0.5  # an Alu counts as "in a gene" if >=50% of it is covered by a gene feature

pybedtools.helpers.set_tempdir(OUT_DIR)

# ## Load scored Alus

alu_df = pd.read_csv(ALU_SCORED_PATH, sep='\t',
                      usecols=['CHROM', 'POS', 'END', 'orig_idx', 'single_gene_name'])
alu_df['POS'] = alu_df['POS'].astype(int)
alu_df['END'] = alu_df['END'].astype(int)

n_alu = len(alu_df)
print(f'{n_alu:,} scored, blacklist-filtered Alus loaded')

alu_bed = BedTool.from_dataframe(alu_df[['CHROM', 'POS', 'END', 'orig_idx']])

# ## Gene annotations (GENCODE v29, all gene types)


gencode_df = pd.read_csv(GENCODE_PATH, sep='\t', comment='#', header=None, compression='gzip', low_memory=False,
                          names=['chrom', 'source', 'feature', 'start', 'end', 'score', 'strand', 'score2', 'info'])

gencode_df_genes = gencode_df[gencode_df['feature'] == 'gene'].copy()
parsed = gencode_df_genes['info'].apply(parse_attributes)
gencode_df_genes = pd.concat([gencode_df_genes.reset_index(drop=True), pd.DataFrame(parsed.tolist())], axis=1)

print(f'{len(gencode_df_genes):,} GENCODE gene features')

gene_bed = BedTool.from_dataframe(gencode_df_genes[['chrom', 'start', 'end', 'gene_name', 'gene_type']])

# ## Blacklist regions + genome file (for shuffling)

blacklist_bed = BedTool(BLACKLIST_PATH)
print(f'{blacklist_bed.count():,} ENCODE blacklist regions (excluded from random placement)')

chrom_sizes = pd.read_csv(CHROM_LENGTHS_PATH, sep='\t', header=None, names=['chrom', 'size'])
chrom_sizes['chrom'] = 'chr' + chrom_sizes['chrom'].astype(str)

genome_file = f'{OUT_DIR}hg38_genome.txt'
chrom_sizes.to_csv(genome_file, sep='\t', header=False, index=False)

# ## Observed proportion of Alus falling in a gene
# 
# (Recomputed directly with the same >=50%-overlap rule used to originally annotate
# `single_gene_name`, as a sanity check against that column.)

def prop_in_gene(bed, n_total):
    hits = bed.intersect(gene_bed, u=True, f=OVERLAP_FRAC)
    n_hit = hits.count()
    return n_hit / n_total

observed_prop = prop_in_gene(alu_bed, n_alu)
observed_prop_from_column = alu_df['single_gene_name'].notna().mean()

print(f'observed proportion of Alus in a gene (recomputed):     {observed_prop:.4f}')
print(f'observed proportion of Alus in a gene (single_gene_name): {observed_prop_from_column:.4f}')

# ## Null distribution: shuffle Alus randomly across the genome
# 
# Each Alu is moved to a uniformly random position anywhere in the genome (not
# restricted to its original chromosome), keeping its length fixed, and rejecting
# placements that fall in a blacklist region. Repeated `N_SHUFFLES` times.
# 
# ~a few seconds per shuffle at this scale (~1.1M intervals), so the full run takes
# several minutes.

def shuffle_once_and_score(seed):
    """Shuffle alu_bed genome-wide (excl. blacklist) and return the proportion overlapping a gene.

    Deletes only the two temp files this call creates, so `alu_bed` and `gene_bed`
    (needed by every other iteration) are never touched -- `pybedtools.cleanup()`
    deletes *all* temp files from the session and would otherwise invalidate them.
    """
    shuffled_bed = alu_bed.shuffle(g=genome_file, excl=BLACKLIST_PATH, chrom=False, seed=seed)
    hits = shuffled_bed.intersect(gene_bed, u=True, f=OVERLAP_FRAC)
    prop = hits.count() / n_alu

    for fn in (shuffled_bed.fn, hits.fn):
        if os.path.exists(fn):
            os.remove(fn)

    return prop

null_props = []

t0 = time.time()
for i in range(N_SHUFFLES):
    null_props.append(shuffle_once_and_score(SEED + i))
    if (i + 1) % 10 == 0:
        elapsed = time.time() - t0
        print(f'{i + 1}/{N_SHUFFLES} shuffles done ({elapsed:.0f}s elapsed)')

null_props = np.array(null_props)
pybedtools.cleanup()
print(f'\ndone in {time.time() - t0:.0f}s')

# ## Compare observed vs. random expectation

null_mean = null_props.mean()
null_sd = null_props.std(ddof=1)
ci_lo, ci_hi = np.percentile(null_props, [2.5, 97.5])

# empirical two-sided p-value: how often does a random shuffle deviate from the
# null mean at least as much as the observed value does
null_p = (np.abs(null_props - null_mean) >= np.abs(observed_prop - null_mean)).mean()
null_p = max(null_p, 1 / (N_SHUFFLES + 1))  # floor at the resolution of the permutation test
z = (observed_prop - null_mean) / null_sd

summary = pd.Series({
    'n_alu': n_alu,
    'n_shuffles': N_SHUFFLES,
    'observed_proportion_in_gene': observed_prop,
    'random_expected_proportion_in_gene': null_mean,
    'random_sd': null_sd,
    'random_95CI_lo': ci_lo,
    'random_95CI_hi': ci_hi,
    'enrichment_ratio': observed_prop / null_mean,
    'z_score': z,
    'empirical_p': null_p,
})
summary

pd.DataFrame({'shuffle_proportion_in_gene': null_props}).to_csv(f'{OUT_DIR}null_distribution_proportion_in_gene.txt', sep='\t', index=False)
summary.to_csv(f'{OUT_DIR}summary.txt', sep='\t', header=False)
print(f'wrote results to {OUT_DIR}')

# ## Plot

fig, ax = plt.subplots(figsize=(6, 4), dpi=150)

ax.hist(null_props, bins=20, color='#6E7B91', edgecolor='white', linewidth=0.5,
        label=f'random placement\n(n={N_SHUFFLES} shuffles)')
ax.axvline(observed_prop, color='#C2410C', linewidth=2, zorder=5)
ax.text(observed_prop, ax.get_ylim()[1] * 0.97, f'  observed = {observed_prop:.3f}',
        color='#C2410C', va='top', ha='left', fontsize=10, fontweight='bold')

ax.set_xlabel('proportion of Alus overlapping a gene (\u226550% of Alu length)')
ax.set_ylabel('number of shuffles')
ax.set_title('Observed Alu-gene overlap vs. random genomic placement')
ax.spines[['top', 'right']].set_visible(False)
ax.legend(frameon=False, loc='upper left' if observed_prop > null_mean else 'upper right')

fig.tight_layout()
fig.savefig(f'{OUT_DIR}observed_vs_random_gene_overlap.png', dpi=300, bbox_inches='tight')
fig.show()
