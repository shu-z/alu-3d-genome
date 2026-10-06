"""Converted from plot_regions_tracks.ipynb; logic unchanged."""

import os
import subprocess
import sys, os, subprocess

import matplotlib
import numpy as np
import pandas as pd
import pygenometracks.tracksClass

# Paths - edit for your environment
PROJECT_DIR = os.environ.get("ALU_PROJECT_DIR", "/pollard/home/szhang20/alu")
DATA_DIR = os.environ.get("ALU_DATA_DIR", "/pollard/data/projects/shzhang")
POLLARD_DATA = os.environ.get("POLLARD_DATA", "/pollard/data")

# # Batch pygenometracks region plots
# 
# Given a list of `(chrom, start, end)` regions, generate the same MSE / genes / CTCF /
# H3K27me3 / H3K9me3 / phyloP100way / segdup / self-chain / RepeatMasker track plot for
# each one.
# 
# Edit the `regions` list near the bottom and re-run.

# %load_ext autoreload
# %autoreload 2

os.chdir("~/akita_variant_scoring/")

pd.set_option('display.max_columns', None)

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype']  = 42

# ## Paths / config
# 
# Same source files as the original notebook — edit here if they move.

res_dir = f'{PROJECT_DIR}/results/paper_results/'

# Alu annotation table (source of the MSE track + region subsetting)
alu = pd.read_csv(
    f'{res_dir}20260316_merged_aluhg38_all_featureannot_NA_blacklist_filter_region_alu.txt',
    sep='\t', index_col=0)

# Static annotation source files
SEGDUP_GZ    = f'{DATA_DIR}/alus/data/annotations/genomicSuperDups.txt.gz'
CHAINSELF_GZ = f'{DATA_DIR}/alus/data/annotations/chainSelf.txt.gz'
RMSK_GZ      = f'{DATA_DIR}/alus/data/repeats/rmsk.txt.gz'

# BED12 gene models (with exon/CDS structure) built from GENCODE v29 -- see
# ~/alu/bin/plot_maps/build_gene_bed12.py. One representative transcript per gene
# (appris_principal > basic > longest), matching the gene set in the older flat
# genes_oneper.gencode.v31.fixed.bed (genes not found in v29 fall back to a single block).
GENES_BED    = f'{DATA_DIR}/alus/data/annotations/genes_oneper_bed12.gencode.v29.bed'
CTCF_BW      = f'{DATA_DIR}/alus/data/encode/bw/ENCFF209TQB_CTCF.bigWig'
H3K27ME3_BW  = f'{DATA_DIR}/alus/data/encode/bw/HFF_H3K27me3_ENCFF854VTY.bigWig'
H3K9ME3_BW   = f'{DATA_DIR}/alus/data/encode/bw/HFF_H3K9me3_ENCFF542CZT.bigWig'
PHYLOP_BW    = f'{POLLARD_DATA}/wynton/consortia/goldenPath/hg38/phyloP100way/hg38.phyloP100way.bw'

# Output location for the finished PDFs
OUTDIR  = f'{PROJECT_DIR}/figs/paper_figs/region_tracks/'
# Scratch location for intermediate bed/bedgraph files (one subfolder per region)
WORKDIR = f'{PROJECT_DIR}/figs/paper_figs/region_tracks/_work/'

os.makedirs(OUTDIR, exist_ok=True)
os.makedirs(WORKDIR, exist_ok=True)

# Genome-wide MSE quantiles used to color-split the MSE track (same thresholds as original: p99 / median)
p99_mse = alu['mse_mean'].quantile(0.99)
p50_mse = alu['mse_mean'].quantile(0.50)

# ## Track-building helpers

gap_frac = 0.15  # fraction of bin width trimmed as a gap between bars, for the bar-style bed tracks

def bedgraph_to_bar_bed(frac_bedgraph, outfile, gap_frac=gap_frac):
    with open(frac_bedgraph) as fin, open(outfile, 'w') as fout:
        for line in fin:
            chrom, start, end, frac = line.strip().split('\t')
            start, end = int(start), int(end)
            gap = int((end - start) * gap_frac / 2)
            fout.write(f"{chrom}\t{start + gap}\t{end - gap}\t.\t{frac}\t.\n")

def make_binned_annotation_track(chrom, windows_bed, source_gz, chrom_col, start_col, end_col, name, workdir):
    """awk out one chrom's intervals from a UCSC-style gzipped table, take bedtools coverage
    against the region's fixed-size windows, and convert to a bar-plot bed."""
    frac_bedgraph = os.path.join(workdir, f'{name}_frac.bedgraph')
    bar_bed       = os.path.join(workdir, f'{name}_bar.bed')

    cmd = (
        f"zcat {source_gz} "
        f"| awk '${chrom_col}==\"{chrom}\" {{print ${chrom_col}\"\\t\"${start_col}\"\\t\"${end_col}}}' "
        f"| bedtools coverage -a {windows_bed} -b stdin "
        f"| awk '{{print $1\"\\t\"$2\"\\t\"$3\"\\t\"$7}}' > {frac_bedgraph}"
    )
    subprocess.run(cmd, shell=True, check=True)
    bedgraph_to_bar_bed(frac_bedgraph, bar_bed)
    return bar_bed

def build_tracks_ini(mse_mid, mse_low, mse_high, mse_all, segdup_bar, chainself_bar, rmsk_bar):
    return f"""
[mse mid]
file = {mse_mid}
title = MSE
height = 5
file_type = bedgraph
type = points:1
color = #e5e5e5
min_value = 0
max_value = 0.05
show_data_range = false

[mse low]
file = {mse_low}
file_type = bedgraph
type = points:1
color = #6AA1C8
min_value = 0
max_value = 0.05
overlay_previous = share-y
show_data_range = false

[mse high]
file = {mse_high}
file_type = bedgraph
type = points:1
color = #C15C6D
min_value = 0
max_value = 0.05
overlay_previous = share-y
show_data_range = false

[mse mean]
file = {mse_all}
file_type = bedgraph
type = line:1.5
color = #996600
summary_method = mean
binsize = 5000
min_value = 0
max_value = 0.05
overlay_previous = share-y
show_data_range = false

[spacer]
height = 0.3

[genes]
file = {GENES_BED}
title = Genes
height = 6
file_type = bed
color = #cccccc
labels = true
show_data_range = false
line_width = 0.5

[spacer]
height = 0.3

[ctcf]
file = {CTCF_BW}
title = CTCF
height = 1.5
file_type = bigwig
color = #CCAA7A
display = collapsed
show_data_range = false

[spacer]
height = 0.3

[h3k9me3]
file = {H3K27ME3_BW}
title = H3K27me3
height = 1.5
file_type = bigwig
color = #A55194
display = collapsed
show_data_range = false

[spacer]
height = 0.3

[h3k27me3]
file = {H3K9ME3_BW}
title = H3K9me3
height = 1.5
file_type = bigwig
color = #9C9EDE
display = collapsed
show_data_range = false

[spacer]
height = 0.3

[phylop100]
file = {PHYLOP_BW}
title = phyloP100way
height = 1.5
file_type = bigwig
color = #8CA252
display = collapsed
show_data_range = false

[spacer]
height = 0.3

[segdup]
file = {segdup_bar}
title = Segmental Dups
height = 1
file_type = bed
color = Greys
min_value = 0
max_value = 1
border_color = none
display = collapsed
labels = false
show_data_range = false

[chainself]
file = {chainself_bar}
title = Self Chain
height = 1
file_type = bed
color = Greys
min_value = 0
max_value = 1
border_color = none
display = collapsed
labels = false
show_data_range = false

[rmsk]
file = {rmsk_bar}
title = RepeatMasker
height = 1
file_type = bed
color = Greys
min_value = 0
max_value = 1
border_color = none
display = collapsed
labels = false
show_data_range = false

[spacer]
height = 0.3

[x-axis]
fontsize = 10

[spacer]
height = 0.3
"""

# ## Main entry point: one region in, one PDF out

def plot_region_tracks(chrom, start, end, binsize=5000, label=None):
    """Reproduce the Step-1/1b/2/3 pipeline from the original notebook for one region.
    Returns the path to the saved PDF."""
    tag     = label if label else f'{chrom}_{start}_{end}'
    workdir = os.path.join(WORKDIR, tag)
    os.makedirs(workdir, exist_ok=True)

    # Step 1b: binned annotation tracks (segdup / self-chain / RepeatMasker)
    windows_bed = os.path.join(workdir, 'windows.bed')
    with open(windows_bed, 'w') as f:
        for b in range(start, end, binsize):
            f.write(f"{chrom}\t{b}\t{min(b + binsize, end)}\n")

    segdup_bar    = make_binned_annotation_track(chrom, windows_bed, SEGDUP_GZ,    2, 3, 4, 'segdup',    workdir)
    chainself_bar = make_binned_annotation_track(chrom, windows_bed, CHAINSELF_GZ, 3, 5, 6, 'chainself', workdir)
    rmsk_bar      = make_binned_annotation_track(chrom, windows_bed, RMSK_GZ,      6, 7, 8, 'rmsk',      workdir)

    # Step 1: MSE bedgraphs split by genome-wide percentile
    alu_plot = alu[(alu['CHROM'] == chrom) & (alu['POS'] > start) & (alu['POS'] < end)].sort_values('POS').copy()
    alu_plot['chrom'] = chrom
    alu_plot['start'] = alu_plot['POS'].astype(int)
    alu_plot['end']   = alu_plot['start'] + 1
    cols = ['chrom', 'start', 'end', 'mse_mean']

    mse_all  = os.path.join(workdir, 'mse_all.bedgraph')
    mse_mid  = os.path.join(workdir, 'mse_mid.bedgraph')
    mse_low  = os.path.join(workdir, 'mse_low.bedgraph')
    mse_high = os.path.join(workdir, 'mse_high.bedgraph')

    alu_plot[cols].to_csv(mse_all, sep='\t', header=False, index=False)
    alu_plot[alu_plot['mse_mean'] >= p50_mse][cols].to_csv(mse_mid,  sep='\t', header=False, index=False)
    alu_plot[alu_plot['mse_mean'] <  p50_mse][cols].to_csv(mse_low,  sep='\t', header=False, index=False)
    alu_plot[alu_plot['mse_mean'] >  p99_mse][cols].to_csv(mse_high, sep='\t', header=False, index=False)

    # Step 2: tracks.ini
    tracks_ini = os.path.join(workdir, 'tracks.ini')
    with open(tracks_ini, 'w') as f:
        f.write(build_tracks_ini(mse_mid, mse_low, mse_high, mse_all, segdup_bar, chainself_bar, rmsk_bar))

    # Step 3: plot
    trp = pygenometracks.tracksClass.PlotTracks(
        tracks_ini,
        fig_width=15,
        fig_height=4,
        dpi=300,
        plot_regions=[(chrom, start, end)]
    )
    out_pdf = os.path.join(OUTDIR, f'region_tracks_{tag}_w15h4.pdf')
    trp.plot(out_pdf, chrom, start, end)
    return out_pdf

# ## Regions to plot
# 
# Edit this list — each entry is `(chrom, start, end)` or `(chrom, start, end, label)`.

regions = [
    ('chr5', 71000001, 72000000),
    ('chr5', 141000001, 142000000),
    ('chr16', 89000001, 90000000),
    ('chr6', 26000001, 27000000),
    ('chr9', 136000001, 137000000),
    ('chr16', 2000001, 3000000),
    ('chr16', 30000001, 31000000),
    ('chr19', 58000001, 58617616),
    ('chrX', 153000001, 154000000),
    ('chr12', 49000001, 50000000),
    ('chr6', 27000001, 28000000),
    ('chr17', 1000001, 2000000),
    ('chr19', 8000001, 9000000),
    ('chr4', 184000001, 185000000),
    ('chr6', 25000001, 26000000),
    ('chr9', 137000001, 138000000),
    ('chr17', 7000001, 8000000),
    ('chr2', 96000001, 97000000),
    ('chrX', 154000001, 155000000),
    ('chr11', 65000001, 66000000),
    ('chr5', 72000001, 73000000),
    ('chr5', 140000001, 141000000),
    ('chr16', 3000001, 4000000),
    ('chr14', 21000001, 22000000),
    ('chr4', 2000001, 3000000),
    ('chr19', 55000001, 56000000),
]

saved_pdfs = []
for region in regions:
    chrom, start, end = region[0], region[1], region[2]
    label = region[3] if len(region) > 3 else None
    try:
        out_pdf = plot_region_tracks(chrom, start, end, label=label)
        print(f'saved {chrom}:{start}-{end} -> {out_pdf}')
        saved_pdfs.append(out_pdf)
    except Exception as e:
        print(f'error plotting {chrom}:{start}-{end}: {e}')

saved_pdfs
