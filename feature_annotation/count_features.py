#!/usr/bin/env python3
"""Count genomic features overlapping Alus and their flanking windows."""

import argparse
from multiprocessing import Pool
from pathlib import Path
from pybedtools import BedTool, cleanup
import glob
import math
import numpy as np
import os
import pandas as pd
import pyBigWig
import pysam
import re
import time


SUPREMO_DIR = os.environ.get("SUPREMO_DIR", "/pollard/home/szhang20/akita_variant_scoring")
AKITA_DIR = os.environ.get("AKITA_DIR", "/pollard/home/szhang20/akita")
PROJECT_DIR = os.environ.get("ALU_PROJECT_DIR", "/pollard/home/szhang20/alu")
PROJECTS_DIR = os.environ.get("ALU_DATA_DIR", "/pollard/data/projects/shzhang")
POLLARD_DATA = os.environ.get("POLLARD_DATA", "/pollard/data")
DATA_DIR = f'{PROJECTS_DIR}/alus/data/'
RES_DIR  = f'{PROJECTS_DIR}/alus/results/20251125_aluhg38_all_scores/'
OUT_DIR=f'{PROJECT_DIR}/results/20260208_feature_overlap_all/'
CHROM_LENGTHS_PATH = '~/akita_variant_scoring/data/chrom_lengths_hg38'
FASTA_PATH         = f'{SUPREMO_DIR}/data/hg38.fa'
PROMOTER_PATH=f'{DATA_DIR}annotations/gencode.v48.promoter_windows.gff3.gz'
WGBS_PATH=f'{DATA_DIR}encode/WGBS/H1_WGBS_ENCFF434CNG.bed.gz'
CTCF_PATH=f'{DATA_DIR}encode/H1_CTCF_ENCFF093VEE.bed.gz'
SEQ_LENGTH   = 1048576
HALF_SEQ     = SEQ_LENGTH // 2
WINDOW_SIZES = [1000, 10000, 100000]
UCSC_TRACKS=f'{POLLARD_DATA}/wynton/consortia/goldenPath/hg38/'
CENTROMERE_PATH = '~/akita_variant_scoring/data/centromere_coords_hg38'
fasta_open = pysam.Fastafile(FASTA_PATH)
map_bw_36=f'{DATA_DIR}annotations/k36.Umap.MultiTrackMappability.bw'
map_bw_100=f'{DATA_DIR}annotations/k100.Umap.MultiTrackMappability.bw'
phyloP30=f'{UCSC_TRACKS}phyloP30way/hg38.phyloP30way.bw'
phyloP100=f'{UCSC_TRACKS}phyloP100way/hg38.phyloP100way.bw'
phastCon30=f'{UCSC_TRACKS}phastCons30way/hg38.phastCons30way.bw'
phastCon100=f'{UCSC_TRACKS}phastCons100way/hg38.phastCons100way.bw'
OUT_DIR__CHROMHMM  = f'{PROJECTS_DIR}/alus/data/annotations/chromhmm/overlap/'
UNIVERSE_DIR = f'{PROJECTS_DIR}/alus/data/annotations/chromhmm/state_universe/'
CHROM_LIST   = ['chr' + str(i) for i in range(1, 23)] + ['chrX', 'chrY']
PEAK_CLASSES = [
    'constitutive',
    'broadly_shared',
    'lineage_restricted',
    'cell_type_specific',
    'singleton',
]
OUT_DIR__CTCF  = f'{PROJECT_DIR}/results/20260208_feature_overlap_all/'
PROMOTER_PATH__CTCF = f'{DATA_DIR}annotations/gencode.v48.promoter_windows.gff3.gz'
WGBS_PATH__CTCF     = f'{DATA_DIR}encode/H1_WGBS_ENCFF434CNG.bed.gz'
CTCF_UNIVERSE_PATH = f'{DATA_DIR}encode/CTCF/CTCF_shared_peaks.txt'
RES_DIR__WGBS = f'{PROJECT_DIR}/results/paper_results/'
OUT_DIR__WGBS = RES_DIR__WGBS
WGBS_DIR = f'{DATA_DIR}encode/WGBS/'  # directory with WGBS bed.gz files


def parse_info_gff3(info_str):
    """
    Parse a GFF3-style info/attributes column:
    key=value;key=value;...
    """
    info = {}
    for field in info_str.strip().split(';'):
        if not field:
            continue
        if '=' in field:
            key, value = field.split('=', 1)
            info[key] = value
    return info


def prepare_alu_windows(alu_df: pd.DataFrame, chrom_lengths: pd.DataFrame, size: int) -> pd.DataFrame:
    df = alu_df.copy()
    df = df.merge(chrom_lengths, on='CHROM', how='left')
    df['alu_mid'] = (df['POS'] + df['END']) // 2
    df['start'] = (df['alu_mid'] - size//2).clip(lower=0)
    df['end'] = df[['alu_mid','chrom_length']].min(axis=1) + (size//2)
    return df[['CHROM','start','end','orig_idx']]


def load_chrom_lengths(path: str) -> pd.DataFrame:
    chrom_lengths = pd.read_table(path, header=None, names=['CHROM', 'chrom_max'])
    chrom_lengths_merge = chrom_lengths.copy()
    chrom_lengths_merge.columns = ['CHROM', 'chrom_length']
    chrom_lengths_merge['CHROM'] = 'chr' + chrom_lengths_merge['CHROM'].astype(str)
    return chrom_lengths_merge


def calculate_gc(chrom, start, end, row_idx):
    try:
        seq = fasta_open.fetch(chrom, start, end).upper()
        if len(seq) == 0:
            return float('nan')
        gc_count = seq.count('G') + seq.count('C')
        return gc_count / len(seq)
    except Exception as e:
        print(f'GC error row {row_idx}: {e}')
        return float('nan')


def average_region_bw(chrom, grab_start, grab_end, map_bw, idx):
    try:
        # stats returns a list (one value per bin). Here we ask for 1 bin = mean over region
        m = map_bw.stats(chrom, int(grab_start), int(grab_end), type="mean")[0]
        if m is None:
            return float('nan')
        return m
        
    except (ValueError, RuntimeError) as e:
        print(f'BigWig error row {idx}: {e}')
        return float('nan')


def count_region_bw(chrom, grab_start, grab_end, map_bw, idx):
    try:
        # Get raw per-base values across the region
        vals = map_bw.values(chrom, int(grab_start), int(grab_end))
        if vals is None:
            return float('nan'), float('nan')
        
        vals = np.array(vals)
        region_len   = len(vals)
        vals = vals[~np.isnan(vals)]  # remove NaN bases (gaps in coverage)
        
        # using phyloP > 2 for conserved, <0 for not conserved
        count_above =  int(np.sum(vals >=  2))
        count_below =  int(np.sum(vals <= -2))
        frac_above   = count_above / region_len
        frac_below   = count_below / region_len
        
        return frac_above, frac_below
    
    except Exception as e:
        print(f"Error at index {idx}: {e}")
        return float('nan'), float('nan')  


def prepare_alu_windows__bigwig(df, chrom_lengths_merge, size):
    df_region = df.copy()
    df_region = df_region.merge(chrom_lengths_merge, on='CHROM', how='left')
    df_region['alu_mid'] = ((df_region['POS'] + df_region['END']) // 2).astype(int)
    df_region['region_start'] = (df_region['alu_mid'] - size // 2).clip(lower=0)
    df_region['region_end'] = df_region['alu_mid'] + size // 2
    df_region['region_end'] = df_region[['region_end', 'chrom_length']].min(axis=1)
    return df_region


def count_region_bw__bigwig_mp(chrom, grab_start, grab_end, map_bw, idx):
    try:
        vals = map_bw.values(chrom, int(grab_start), int(grab_end), numpy=True)  # return numpy array directly
        if vals is None:
            return float('nan'), float('nan')
        
        region_len = len(vals)
        mask = ~np.isnan(vals)
        
        # avoid creating filtered array — operate on mask directly
        frac_above = np.sum((vals > 2) & mask) / region_len
        frac_below = np.sum((vals < -2) & mask) / region_len
        
        return frac_above, frac_below
    
    except Exception as e:
        print(f"Error at index {idx}: {e}")
        return float('nan'), float('nan')  


def prepare_alu_windows__bigwig_mp(df, chrom_lengths_merge, size):
    df_region = df.copy()
    df_region = df_region.merge(chrom_lengths_merge, on='CHROM', how='left')
    df_region['alu_mid'] = ((df_region['POS'] + df_region['END']) // 2).astype(int)
    df_region['region_start'] = (df_region['alu_mid'] - size // 2).clip(lower=0)
    df_region['region_end'] = df_region['alu_mid'] + size // 2
    df_region['region_end'] = df_region[['region_end', 'chrom_length']].min(axis=1)
    return df_region


def init_worker(bw_path):
    global bw_handle
    bw_handle = pyBigWig.open(bw_path)


def count_region_worker(args):
    chrom, start, end, idx = args
    vals = bw_handle.values(chrom, int(start), int(end), numpy=True)
    if vals is None:
        return float('nan'), float('nan')
    region_len = len(vals)
    mask = ~np.isnan(vals)
    return np.sum((vals >= 2) & mask) / region_len, np.sum((vals <= -2) & mask) / region_len


def prepare_alu_windows__chromhmm(alu_df: pd.DataFrame,
                        chrom_lengths: pd.DataFrame,
                        size: int) -> pd.DataFrame:
    """
    Build fixed-size windows centered on each Alu midpoint.
    Returns DataFrame with columns: CHROM, region_start, region_end, orig_idx
    """
    df = alu_df[['CHROM', 'POS', 'END', 'orig_idx']].copy()
    df = df.merge(chrom_lengths, on='CHROM', how='left')
    df['alu_mid']      = (df['POS'] + df['END']) // 2
    df['region_start'] = (df['alu_mid'] - size // 2).clip(lower=0)
    df['region_end']   = (df['alu_mid'] + size // 2).clip(upper=df['chrom_length'])
    return df[['CHROM', 'region_start', 'region_end', 'orig_idx']]


def compute_fraction_overlap_bychrom(region_df: pd.DataFrame,
                                     feature_df: pd.DataFrame,
                                     colname: str) -> pd.DataFrame:
    """
    For each row in region_df, compute the fraction of the region covered by
    features in feature_df (already a DataFrame with columns chrom/start/end).

    region_df must have columns: CHROM, region_start, region_end, orig_idx
    feature_df must have columns: chrom, start, end  (standard BED)

    Returns DataFrame with columns: orig_idx, colname
    """
    df = region_df.dropna(subset=['region_start', 'region_end']).copy()
    df['region_start'] = df['region_start'].astype(int)
    df['region_end']   = df['region_end'].astype(int)

    overlap_totals = {}

    for chrom in CHROM_LIST:
        df_chrom = df[df['CHROM'] == chrom]
        if df_chrom.empty:
            continue

        feat_chrom = feature_df[feature_df['chrom'] == chrom]
        if feat_chrom.empty:
            continue

        print(f'      {chrom} ({len(df_chrom)} regions, {len(feat_chrom)} features)')

        # build BedTools from DataFrames — columns must be in BED order
        alu_bed_chrom  = BedTool.from_dataframe(
            df_chrom[['CHROM', 'region_start', 'region_end', 'orig_idx']])
        feat_bed_chrom = BedTool.from_dataframe(
            feat_chrom[['chrom', 'start', 'end']]).sort().merge()

        inter = alu_bed_chrom.intersect(feat_bed_chrom, wo=True)

        for entry in inter:
            try:
                alu_id      = int(entry[3])
                overlap_len = int(entry[-1])
            except (ValueError, IndexError):
                print(f'      WARNING: unexpected entry format: {entry}')
                continue
            overlap_totals[alu_id] = overlap_totals.get(alu_id, 0) + overlap_len

        cleanup()

    # compute fraction: overlap / region length
    df[colname] = df.apply(
        lambda row: overlap_totals.get(row['orig_idx'], 0)
                    / (row['region_end'] - row['region_start'])
        if (row['region_end'] - row['region_start']) > 0 else 0.0,
        axis=1
    )

    return df[['orig_idx', colname]]


def peak_class_to_df(presence: pd.DataFrame, peak_class: str) -> pd.DataFrame:
    """Subset universe presence DataFrame to one peak class, return as DataFrame."""
    return presence[presence['peak_class'] == peak_class][['chrom', 'start', 'end']].copy()


def _run_all_windows(feat_name: str,
                     feat_df: pd.DataFrame,
                     alu_df_state: pd.DataFrame,
                     alu_annot: pd.DataFrame,
                     chrom_lengths_merge: pd.DataFrame) -> pd.DataFrame:
    """
    Run fraction overlap at Alu body, 1kb/10kb/100kb windows, and 1Mb window
    for one feature DataFrame, merging results back into alu_df_state.
    """

    ##### alu body 
    print('    alu body')
    alu_body_region = alu_annot[['CHROM', 'POS', 'END', 'orig_idx']].copy()
    # rename to the expected region_start / region_end interface
    alu_body_region = alu_body_region.rename(columns={'POS': 'region_start',
                                                       'END': 'region_end'})
    alu_df_state = alu_df_state.merge(
        compute_fraction_overlap_bychrom(
            alu_body_region, feat_df, f'{feat_name}_frac_alu'),
        on='orig_idx', how='left')

    ##### windows
    for size in WINDOW_SIZES:
        print(f'    {size // 1000}kb window')
        df_region = prepare_alu_windows__chromhmm(alu_annot, chrom_lengths_merge, size)
        alu_df_state = alu_df_state.merge(
            compute_fraction_overlap_bychrom(
                df_region, feat_df, f'{feat_name}_frac_{size // 1000}kb'),
            on='orig_idx', how='left')

    ##### 1mb window 
    print('    1Mb window')
    df_1mb = alu_annot[['CHROM', 'REF_start', 'REF_stop', 'orig_idx']].copy()
    df_1mb['REF_start'] = df_1mb['REF_start'].astype(int)
    df_1mb['REF_stop']  = df_1mb['REF_stop'].astype(int)
    df_1mb = df_1mb.rename(columns={'REF_start': 'region_start',
                                     'REF_stop':  'region_end'})
    alu_df_state = alu_df_state.merge(
        compute_fraction_overlap_bychrom(
            df_1mb, feat_df, f'{feat_name}_frac_1Mb'),
        on='orig_idx', how='left')

    return alu_df_state


def parse_info_gff3__ctcf(info_str):
    info = {}
    for field in info_str.strip().split(';'):
        if not field:
            continue
        if '=' in field:
            key, value = field.split('=', 1)
            info[key] = value
    return info


def prepare_alu_windows__ctcf(alu_df: pd.DataFrame,
                        chrom_lengths: pd.DataFrame,
                        size: int) -> pd.DataFrame:
    df = alu_df.copy()
    df = df.merge(chrom_lengths, on='CHROM', how='left')
    df['alu_mid'] = (df['POS'] + df['END']) // 2
    df['start']   = (df['alu_mid'] - size // 2).clip(lower=0)
    df['end']     = df[['alu_mid', 'chrom_length']].min(axis=1) + (size // 2)
    return df[['CHROM', 'start', 'end', 'orig_idx']]


def count_overlaps(query_bed: BedTool,
                   feature_bed: BedTool,
                   query_names: list,
                   count_col: str) -> pd.DataFrame:
    """Intersect query with feature, return just orig_idx + count column."""
    result = query_bed.intersect(feature_bed, c=True, F=0.5).to_dataframe(
        names=query_names + [count_col])
    return result[['orig_idx', count_col]]


def peak_class_to_bed(presence: pd.DataFrame, peak_class: str) -> BedTool:
    """Subset the universe presence DataFrame to one peak class and return a BedTool."""
    subset = presence[presence['peak_class'] == peak_class][['chrom', 'start', 'end']]
    return BedTool.from_dataframe(subset)


def _run_all_windows__ctcf(feat_name, feat_BED, alu_df_copy, alu_annot, chrom_lengths_merge):
    """
    Run overlap counts at Alu body, 1kb/10kb/100kb windows, and 1Mb window
    for one feature BED, and merge results back into alu_df_copy.
    """
    # Alu body overlap
    print('  alu body')
    alu_body_BED = BedTool.from_dataframe(
        alu_annot[['CHROM', 'POS', 'END', 'orig_idx']])
    alu_df_copy = alu_df_copy.merge(
        count_overlaps(alu_body_BED, feat_BED,
                       ['chrom', 'start', 'end', 'orig_idx'],
                       f'{feat_name}_count_alu'),
        on='orig_idx')

    # fixed-size windows
    for size in WINDOW_SIZES:
        print(f'  {size//1000}kb window')
        df_region     = prepare_alu_windows__ctcf(alu_df_copy, chrom_lengths_merge, size)
        df_region_BED = BedTool.from_dataframe(df_region)
        alu_df_copy   = alu_df_copy.merge(
            count_overlaps(df_region_BED, feat_BED,
                           ['chrom', 'start', 'end', 'orig_idx'],
                           f'{feat_name}_count_{size//1000}kb'),
            on='orig_idx')

    # 1 Mb window (Akita input window)
    print('  1Mb window')
    df_1mb = alu_annot[['CHROM', 'REF_start', 'REF_stop', 'orig_idx']].copy()
    df_1mb['REF_start'] = df_1mb['REF_start'].astype(int)
    df_1mb['REF_stop']  = df_1mb['REF_stop'].astype(int)
    df_1mb_BED          = BedTool.from_dataframe(df_1mb)
    alu_df_copy         = alu_df_copy.merge(
        count_overlaps(df_1mb_BED, feat_BED,
                       ['chrom', 'start', 'end', 'orig_idx'],
                       f'{feat_name}_count_1Mb'),
        on='orig_idx')

    return alu_df_copy


def prepare_alu_windows__wgbs(alu_df: pd.DataFrame, chrom_lengths: pd.DataFrame, size: int) -> pd.DataFrame:
    df = alu_df.copy()
    df = df.merge(chrom_lengths, on='CHROM', how='left')
    df['alu_mid'] = (df['POS'] + df['END']) // 2
    df['start'] = (df['alu_mid'] - size // 2).clip(lower=0)
    df['end'] = df[['alu_mid', 'chrom_length']].min(axis=1) + (size // 2)
    return df[['CHROM', 'start', 'end', 'orig_idx']]


def process_wgbs(wgbs_path: Path, alu_annot: pd.DataFrame, alu_df_copy: pd.DataFrame, 
                 chrom_lengths_merge: pd.DataFrame) -> pd.DataFrame:
    
    feat_name = wgbs_path.stem.replace('.bed', '')  # e.g. H1_WGBS_ENCFF434CNG
    print(f'\nProcessing {feat_name}')

    # Load and filter WGBS
    WGBS = pd.read_csv(wgbs_path, sep='\t', header=None)
    WGBS.columns = ['chrom', 'start', 'end', 'name', 'score', 'strand', 'thickStart', 'thickEnd',
                    'RGB', 'coverage', 'percent_reads_methylated', 'ref_geno', 'samp_geno', 'qual']
    WGBS_high_methyl = WGBS[(WGBS['coverage'] > 5) & (WGBS['percent_reads_methylated'] >= 80)]
    feat_BED = BedTool.from_dataframe(WGBS_high_methyl[['chrom', 'start', 'end']])

    result = alu_df_copy[['orig_idx']].copy()

    # Alu-level overlap
    print('  alu window')
    df_window = alu_annot[['CHROM', 'POS', 'END', 'orig_idx']]
    df_window_BED = BedTool.from_dataframe(df_window)
    overlap_counts = df_window_BED.intersect(feat_BED, c=True, F=0.5).to_dataframe(
        names=['chrom', 'start', 'end', 'orig_idx', f'{feat_name}_count_alu'])
    result = result.merge(overlap_counts[['orig_idx', f'{feat_name}_count_alu']], on='orig_idx')

    # Window overlaps
    for size in WINDOW_SIZES:
        print(f'  {size // 1000}kb window')
        df_region = prepare_alu_windows__wgbs(alu_df_copy, chrom_lengths_merge, size)
        df_region_BED = BedTool.from_dataframe(df_region)
        colname = f'{feat_name}_count_{size // 1000}kb'
        overlap_counts = df_region_BED.intersect(feat_BED, c=True, F=0.5).to_dataframe(
            names=['chrom', 'start', 'end', 'orig_idx', colname])
        result = result.merge(overlap_counts[['orig_idx', colname]], on='orig_idx')

    # 1Mb window
    print('  1Mb window')
    df_window = alu_annot[['CHROM', 'REF_start', 'REF_stop', 'orig_idx']].copy()
    df_window['REF_start'] = df_window['REF_start'].astype(int)
    df_window['REF_stop'] = df_window['REF_stop'].astype(int)
    df_window_BED = BedTool.from_dataframe(df_window)
    overlap_counts = df_window_BED.intersect(feat_BED, c=True, F=0.5).to_dataframe(
        names=['chrom', 'start', 'end', 'orig_idx', f'{feat_name}_count_1Mb'])
    result = result.merge(overlap_counts[['orig_idx', f'{feat_name}_count_1Mb']], on='orig_idx')

    cleanup()
    return result


def add_elements_args(parser):
        parser.add_argument("--chrom-lengths", default=CHROM_LENGTHS_PATH,
                            help="chromosome lengths table")
        parser.add_argument("--ctcf", default=CTCF_PATH,
                            help="CTCF peak bed")
        parser.add_argument("--fasta", default=FASTA_PATH,
                            help="reference FASTA")
        parser.add_argument("--out-dir", default=OUT_DIR,
                            help="directory for output")
        parser.add_argument("--promoter", default=PROMOTER_PATH,
                            help="promoter annotation gff3")
        parser.add_argument("--res-dir", default=RES_DIR,
                            help="directory holding the scored/annotated Alu table")
        parser.add_argument("--wgbs", default=WGBS_PATH,
                            help="WGBS bed")
        parser.add_argument("--window-sizes", default=WINDOW_SIZES, type=int, nargs="+",
                            help="window sizes in bp")


def run_elements(args):
    chrom_lengths_path = args.chrom_lengths
    ctcf_path = args.ctcf
    fasta_path = args.fasta
    out_dir = args.out_dir
    promoter_path = args.promoter
    res_dir = args.res_dir
    wgbs_path = args.wgbs
    window_sizes = args.window_sizes
    chrom_lengths = pd.read_table(chrom_lengths_path, header=None, names=['CHROM','chrom_max'])
    chrom_lengths_merge = chrom_lengths.copy()
    chrom_lengths_merge.columns = ['CHROM','chrom_length']
    chrom_lengths_merge['CHROM'] = 'chr' + chrom_lengths_merge['CHROM'].astype(str)
    fasta_open = pysam.Fastafile(fasta_path)
    promoter_df = pd.read_csv(promoter_path, sep='\t', comment='#', header=None,compression='gzip', low_memory=False,
    names=["chrom", "source", "feature", "start", "end", "score", "strand", "score2", "info"])
    parse_promoter = promoter_df['info'].apply(parse_info_gff3)
    promoter_df = pd.concat([promoter_df, pd.DataFrame(parse_promoter.tolist())], axis=1)
    WGBS=pd.read_csv(wgbs_path, sep='\t', header=None)
    WGBS.columns=['chrom', 'start', 'end', 'name', 'score', 'strand', 'thickStart', 'thickEnd', 
              'RGB', 'coverage', 'percent_reads_methylated', 'ref_geno', 'samp_geno', 'qual']
    WGBS_high_methyl=WGBS[(WGBS['coverage']>5) & (WGBS['percent_reads_methylated']>=80)]
    alu_annot = pd.read_csv(f'{res_dir}20260123_aluhg38_all_annot_GC_mapp_alu_repeats_genetrackoverlap.txt', sep='\t', index_col=0)
    alu_df_copy=alu_annot[['CHROM', 'POS', 'END','orig_idx']]
    feature_dict = {
   #     'alu' : alu_annot[['CHROM', 'POS', 'END','orig_idx']], 
   #     'promoter' : promoter_df[['chrom', 'start', 'end']],
        'HFF_CTCF' : pd.read_csv(ctcf_path, sep='\t', header=None),
    #    'H1ESC_WGBS' : WGBS_high_methyl
    }
    for feat, feat_df in feature_dict.items():
        print(f'processing feature: {feat}')
        feat_BED = BedTool.from_dataframe(feat_df)

        #alu overlap
        print('alu window')
        df_window= alu_annot[['CHROM','POS','END', 'orig_idx']]
        df_window_BED = BedTool.from_dataframe(df_window)

        #half of the other element must be covered to be considered 
        overlap_counts = df_window_BED.intersect(feat_BED, c=True, F=0.5).to_dataframe(
            names=['chrom','start','end', 'orig_idx', f'{feat}_count_alu'])
        #merge counts back
        alu_df_copy=alu_df_copy.merge(overlap_counts[['orig_idx', f'{feat}_count_alu']])

        for size in window_sizes:
            print(f'  Window size: {size}')
            
            df_region = prepare_alu_windows(alu_df_copy, chrom_lengths_merge, size)
            colname = f'{feat}_count_{size//1000}kb'
            df_region_BED = BedTool.from_dataframe(df_region)

            overlap_counts = df_region_BED.intersect(feat_BED, c=True, F=0.5).to_dataframe(
                names=['chrom','start','end','orig_idx', f'{feat}_count_{size//1000}kb'])
            alu_df_copy = alu_df_copy.merge(overlap_counts[['orig_idx', f'{feat}_count_{size//1000}kb']])

        #window overlap
        print('1Mb window')
        df_window= alu_annot[['CHROM','REF_start','REF_stop', 'orig_idx']]
        df_window['REF_start']=df_window['REF_start'].astype(int)
        df_window['REF_stop']=df_window['REF_stop'].astype(int)

        df_window_BED = BedTool.from_dataframe(df_window)
        
        #half of the other element must be covered to be considered 
        overlap_counts = df_window_BED.intersect(feat_BED, c=True, F=0.5).to_dataframe(
            names=['chrom','start','end', 'orig_idx', f'{feat}_count_1Mb'])
        #merge counts back
        alu_df_copy=alu_df_copy.merge(overlap_counts[['orig_idx', f'{feat}_count_1Mb']])
        
        cleanup() #cleans up pybedtools temp files 
    alu_df_copy.to_csv(f'{out_dir}/aluhg38_feature_count_bed_H1_CTCF.csv', index=False)
    print("Done!")


def add_bigwig_args(parser):
        parser.add_argument("--chrom-lengths", default=CHROM_LENGTHS_PATH,
                            help="chromosome lengths table")
        parser.add_argument("--out-dir", default=OUT_DIR,
                            help="directory for output")
        parser.add_argument("--res-dir", default=RES_DIR,
                            help="directory holding the scored/annotated Alu table")
        parser.add_argument("--window-sizes", default=WINDOW_SIZES, type=int, nargs="+",
                            help="window sizes in bp")


def run_bigwig(args):
    chrom_lengths_path = args.chrom_lengths
    out_dir = args.out_dir
    res_dir = args.res_dir
    window_sizes = args.window_sizes
    chrom_lengths = pd.read_table(chrom_lengths_path, header=None, names=['CHROM','chrom_max'])
    chrom_lengths_merge = chrom_lengths.copy()
    chrom_lengths_merge.columns = ['CHROM','chrom_length']
    chrom_lengths_merge['CHROM'] = 'chr' + chrom_lengths_merge['CHROM'].astype(str)
    alu_annot = pd.read_csv(f'{res_dir}20260123_aluhg38_all_annot_GC_mapp_alu_repeats_genetrackoverlap.txt', sep='\t', index_col=0)
    alu_df_copy = alu_annot[['CHROM', 'POS', 'END','orig_idx','REF_start','REF_stop']]
    feature_dict={
        #'mapp36' : map_bw_36,
        #'mapp100' : map_bw_100, 
        #'phyloP30' : phyloP30,
        'phyloP100' : phyloP100, 
        #'phastCon30' : phastCon30,
        #'phastCon100' : phastCon100
    }
    for feat, feat_path in feature_dict.items():
        print(f'processing feature: {feat}')
        feat_bw =pyBigWig.open(feat_path)
        print(f'  alu window')
        
        alu_df_copy[f'{feat}_highfrac_alu'], alu_df_copy[f'{feat}_lowfrac_alu'] = zip(*alu_df_copy.apply(
            lambda row: count_region_bw(row['CHROM'], row['POS'], row['END'], feat_bw, row.name), axis=1))

        for size in window_sizes:
            print(f'  Window size: {size}')
            df_region = prepare_alu_windows__bigwig(alu_df_copy, chrom_lengths_merge, size)
            #merge df_region with alu_df_copy
            df_region[f'{feat}_{size//1000}kb'] = df_region.apply(
            lambda row: average_region_bw(row['CHROM'], row['region_start'], row['region_end'], feat_bw, row.name), axis=1)

            df_region[f'{feat}_highfrac_{size//1000}kb'], df_region[f'{feat}_lowfrac_{size//1000}kb'] = zip(*df_region.apply(
            lambda row: count_region_bw(row['CHROM'], row['region_start'], row['region_end'], feat_bw, row.name), axis=1))
            
            alu_df_copy=alu_df_copy.merge(df_region[['orig_idx', f'{feat}_highfrac_{size//1000}kb', f'{feat}_lowfrac_{size//1000}kb']])

        # print(f'  1Mb window')
        
        # alu_df_copy[f'{feat}_highfrac_1Mb'], alu_df_copy[f'{feat}_lowfrac_1Mb'] = zip(*alu_df_copy.apply(
        #     lambda row: count_region_bw(row['CHROM'], row['REF_start'], row['REF_stop'], feat_bw, row.name), axis=1))
    
            #save intermediate df
            print('saving...')
            alu_df_copy.to_csv(f'{out_dir}/aluhg38_phyloP_countBW_short_-2_2.csv', index=False)
    print(alu_df_copy.head())
    alu_df_copy.to_csv(f'{out_dir}/aluhg38_phyloP_countBW_short_-2_2.csv', index=False)
    print("Done!")


def add_bigwig_mp_args(parser):
        parser.add_argument("--chrom-lengths", default=CHROM_LENGTHS_PATH,
                            help="chromosome lengths table")
        parser.add_argument("--out-dir", default=OUT_DIR,
                            help="directory for output")
        parser.add_argument("--res-dir", default=RES_DIR,
                            help="directory holding the scored/annotated Alu table")


def run_bigwig_mp(args):
    chrom_lengths_path = args.chrom_lengths
    out_dir = args.out_dir
    res_dir = args.res_dir
    chrom_lengths = pd.read_table(chrom_lengths_path, header=None, names=['CHROM','chrom_max'])
    chrom_lengths_merge = chrom_lengths.copy()
    chrom_lengths_merge.columns = ['CHROM','chrom_length']
    chrom_lengths_merge['CHROM'] = 'chr' + chrom_lengths_merge['CHROM'].astype(str)
    alu_annot = pd.read_csv(f'{res_dir}20260123_aluhg38_all_annot_GC_mapp_alu_repeats_genetrackoverlap.txt', sep='\t', index_col=0)
    alu_df_copy = alu_annot[['CHROM', 'POS', 'END','orig_idx','REF_start','REF_stop']]
    feature_dict={
        #'mapp36' : map_bw_36,
        #'mapp100' : map_bw_100, 
        #'phyloP30' : phyloP30,
        'phyloP100' : phyloP100, 
        #'phastCon30' : phastCon30,
        #'phastCon100' : phastCon100
    }
    for feat, feat_path in feature_dict.items():
        print(f'processing feature: {feat}')
        #feat_bw =pyBigWig.open(feat_path)
        # print(f'  alu window')
        
        # alu_df_copy[f'{feat}_highfrac_alu'], alu_df_copy[f'{feat}_lowfrac_alu'] = zip(*alu_df_copy.apply(
        #     lambda row: count_region_bw__bigwig_mp(row['CHROM'], row['POS'], row['END'], feat_bw, row.name), axis=1))

        # for size in WINDOW_SIZES:
        #     print(f'  Window size: {size}')
        #     df_region = prepare_alu_windows__bigwig_mp(alu_df_copy, chrom_lengths_merge, size)
        #     #merge df_region with alu_df_copy
        #     df_region[f'{feat}_{size//1000}kb'] = df_region.apply(
        #     lambda row: average_region_bw(row['CHROM'], row['region_start'], row['region_end'], feat_bw, row.name), axis=1)

        #     df_region[f'{feat}_highfrac_{size//1000}kb'], df_region[f'{feat}_lowfrac_{size//1000}kb'] = zip(*df_region.apply(
        #     lambda row: count_region_bw__bigwig_mp(row['CHROM'], row['region_start'], row['region_end'], feat_bw, row.name), axis=1))
            
        #     alu_df_copy=alu_df_copy.merge(df_region[['orig_idx', f'{feat}_highfrac_{size//1000}kb', f'{feat}_lowfrac_{size//1000}kb']])

        args = [
            (row['CHROM'], row['REF_start'], row['REF_stop'], row.name)
            for _, row in alu_df_copy.iterrows()
        ]

        with Pool(processes=18, initializer=init_worker, initargs=(feat_path,)) as pool:
            results = pool.map(count_region_worker, args)

        frac_above, frac_below = zip(*results)
        alu_df_copy[f'{feat}_highfrac_1Mb'] = list(frac_above)
        alu_df_copy[f'{feat}_lowfrac_1Mb']  = list(frac_below)
        
        # alu_df_copy[f'{feat}_highfrac_1Mb'], alu_df_copy[f'{feat}_lowfrac_1Mb'] = zip(*alu_df_copy.apply(
        #     lambda row: count_region_bw__bigwig_mp(row['CHROM'], row['REF_start'], row['REF_stop'], feat_bw, row.name), axis=1))

        #save intermediate df
        alu_df_copy.to_csv(f'{out_dir}/aluhg38_phyloP_countBW_1Mb_-2_2.csv', index=False)
    print(alu_df_copy.head())
    alu_df_copy.to_csv(f'{out_dir}/aluhg38_phyloP_countBW_1Mb_-2_2.csv', index=False)
    print("Done!")


def add_chromhmm_args(parser):
        parser.add_argument("--chrom-lengths", default=CHROM_LENGTHS_PATH,
                            help="chromosome lengths table")
        parser.add_argument("--out-dir", default=OUT_DIR__CHROMHMM,
                            help="directory for output")
        parser.add_argument("--res-dir", default=RES_DIR,
                            help="directory holding the scored/annotated Alu table")
        parser.add_argument("--universe-dir", default=UNIVERSE_DIR,
                            help="directory of chromHMM state universe beds")


def run_chromhmm(args):
    chrom_lengths_path = args.chrom_lengths
    out_dir = args.out_dir
    res_dir = args.res_dir
    universe_dir = args.universe_dir
    os.makedirs(out_dir, exist_ok=True)
    chrom_lengths_merge = pd.read_table(
        chrom_lengths_path, header=None, names=['CHROM', 'chrom_length'])
    chrom_lengths_merge['CHROM'] = 'chr' + chrom_lengths_merge['CHROM'].astype(str)
    alu_annot = pd.read_csv(
        f'{res_dir}20260123_aluhg38_all_annot_GC_mapp_alu_repeats_genetrackoverlap.txt',
        sep='\t', index_col=0)
    print(f'Loaded {len(alu_annot):,} Alu annotations')
    universe_files = sorted(glob.glob(os.path.join(universe_dir, 'universe_*.csv')))
    universe_files = [f for f in universe_files if '_8_' in f or '_9_' in f]
    print(f'Found {len(universe_files)} state universe files\n')
    all_states_results = []   # accumulate per-state DataFrames for wide merge
    for universe_path in universe_files:

        fname      = os.path.basename(universe_path)
        state_name = fname.replace('universe_', '').replace('.csv', '')
        print(f'\n{"="*60}')
        print(f'State: {state_name}')

        presence = pd.read_csv(universe_path)
        # presence must have columns: chrom, start, end, peak_class
        print('  Peak class counts:')
        print(presence['peak_class'].value_counts().to_string())

        # fresh per-state accumulator — starts with core Alu identity columns only
        alu_df_state = alu_annot[['CHROM', 'POS', 'END', 'orig_idx']].copy()

        # inner loop: peak classes 
        for peak_class in PEAK_CLASSES:
            feat_name = f'{state_name}_{peak_class}'
            n_peaks   = (presence['peak_class'] == peak_class).sum()
            print(f'\n  [{peak_class}]  {n_peaks} peaks')

            if n_peaks == 0:
                print('    no peaks, skipping')
                continue

            # pass feature as a plain DataFrame (not BedTool) to avoid
            # the BedTool[df-style indexing] bug
            feat_df      = peak_class_to_df(presence, peak_class)
            alu_df_state = _run_all_windows(
                feat_name, feat_df, alu_df_state, alu_annot, chrom_lengths_merge)
            cleanup()

        # tag and save per-state file
        alu_df_state['state_name'] = state_name
        per_state_path = os.path.join(out_dir, f'aluhg38_feature_count_{state_name}.csv')
        alu_df_state.to_csv(per_state_path, index=False)
        print(f'\n  Saved per-state → {per_state_path}')

        all_states_results.append(alu_df_state)
    print('\nBuilding combined wide output...')
    id_cols  = ['CHROM', 'POS', 'END', 'orig_idx']
    combined = all_states_results[0].drop(columns='state_name')
    for df in all_states_results[1:]:
        state_cols = id_cols + [c for c in df.columns
                                if c not in id_cols and c != 'state_name']
        combined = combined.merge(df[state_cols], on=id_cols, how='outer')
    combined_path = os.path.join(out_dir, 'aluhg38_feature_count_all_states.csv')
    combined.to_csv(combined_path, index=False)
    print(f'Saved combined → {combined_path}')
    print(f'Combined shape: {combined.shape}')


def add_ctcf_args(parser):
        parser.add_argument("--chrom-lengths", default=CHROM_LENGTHS_PATH,
                            help="chromosome lengths table")
        parser.add_argument("--ctcf-universe", default=CTCF_UNIVERSE_PATH,
                            help="CTCF universe bed")
        parser.add_argument("--fasta", default=FASTA_PATH,
                            help="reference FASTA")
        parser.add_argument("--out-dir", default=OUT_DIR__CTCF,
                            help="directory for output")
        parser.add_argument("--res-dir", default=RES_DIR,
                            help="directory holding the scored/annotated Alu table")


def run_ctcf(args):
    chrom_lengths_path = args.chrom_lengths
    ctcf_universe_path = args.ctcf_universe
    fasta_path = args.fasta
    out_dir = args.out_dir
    res_dir = args.res_dir
    chrom_lengths = pd.read_table(
        chrom_lengths_path, header=None, names=['CHROM', 'chrom_max'])
    chrom_lengths_merge = chrom_lengths.copy()
    chrom_lengths_merge.columns = ['CHROM', 'chrom_length']
    chrom_lengths_merge['CHROM'] = 'chr' + chrom_lengths_merge['CHROM'].astype(str)
    fasta_open = pysam.Fastafile(fasta_path)
    presence = pd.read_csv(ctcf_universe_path, sep='\t')
    print('CTCF universe peak class counts:')
    print(presence['peak_class'].value_counts())
    alu_annot  = pd.read_csv(
        f'{res_dir}20260123_aluhg38_all_annot_GC_mapp_alu_repeats_genetrackoverlap.txt',
        sep='\t', index_col=0)
    alu_df_copy = alu_annot[['CHROM', 'POS', 'END', 'orig_idx']].copy()
    for peak_class in PEAK_CLASSES:
        feat_name = f'CTCF_{peak_class}'
        n_peaks   = (presence['peak_class'] == peak_class).sum()
        print(f'\nprocessing {feat_name}  ({n_peaks} peaks)')

        if n_peaks == 0:
            print(f'  no peaks in class {peak_class}, skipping')
            continue

        feat_BED  = peak_class_to_bed(presence, peak_class)
        alu_df_copy = _run_all_windows__ctcf(
            feat_name, feat_BED, alu_df_copy, alu_annot, chrom_lengths_merge)
        cleanup()
    out_path = f'{out_dir}/aluhg38_feature_count_CTCF_sharedpeaks.csv'
    alu_df_copy.to_csv(out_path, index=False)
    print(f'\nDone! Saved to {out_path}')


def add_wgbs_args(parser):
        parser.add_argument("--chrom-lengths", default=CHROM_LENGTHS_PATH,
                            help="chromosome lengths table")
        parser.add_argument("--out-dir", default=OUT_DIR__WGBS,
                            help="directory for output")
        parser.add_argument("--res-dir", default=RES_DIR__WGBS,
                            help="directory holding the scored/annotated Alu table")
        parser.add_argument("--wgbs-dir", default=WGBS_DIR,
                            help="directory of WGBS beds")


def run_wgbs(args):
    chrom_lengths_path = args.chrom_lengths
    out_dir = args.out_dir
    res_dir = args.res_dir
    wgbs_dir = args.wgbs_dir
    chrom_lengths_merge = pd.read_table(chrom_lengths_path, header=None, names=['CHROM', 'chrom_length'])
    chrom_lengths_merge['CHROM'] = 'chr' + chrom_lengths_merge['CHROM'].astype(str)
    alu_annot = pd.read_csv(f'{res_dir}20260316_merged_aluhg38_all_featureannot_NA_blacklist_filter_region_alu.txt',
                            sep='\t', index_col=0)
    alu_df_copy = alu_annot[['CHROM', 'POS', 'END', 'orig_idx']].copy()
    wgbs_files = sorted(Path(wgbs_dir).glob('*.bed.gz'))
    print(f'Found {len(wgbs_files)} WGBS files')
    for wgbs_path in wgbs_files:
        result = process_wgbs(wgbs_path, alu_annot, alu_df_copy, chrom_lengths_merge)
        alu_df_copy = alu_df_copy.merge(result, on='orig_idx')

        # Save intermediate
        alu_df_copy.to_csv(f'{out_dir}/aluhg38_WGBS_count.csv', index=False)
    alu_df_copy.to_csv(f'{out_dir}/aluhg38_WGBS_count.csv', index=False)
    print('Done!')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    p_elements = sub.add_parser("elements", help="count promoter, CTCF and WGBS features overlapping each Alu and window")
    add_elements_args(p_elements)
    p_elements.set_defaults(func=run_elements)
    p_bigwig = sub.add_parser("bigwig", help="fraction of bases above/below a bigWig threshold (phyloP, phastCons, mappability)")
    add_bigwig_args(p_bigwig)
    p_bigwig.set_defaults(func=run_bigwig)
    p_bigwig_mp = sub.add_parser("bigwig_mp", help="same as bigwig, parallelised over a process pool for the 1Mb window")
    add_bigwig_mp_args(p_bigwig_mp)
    p_bigwig_mp.set_defaults(func=run_bigwig_mp)
    p_chromhmm = sub.add_parser("chromhmm", help="chromHMM state coverage per Alu and window, from the state universe")
    add_chromhmm_args(p_chromhmm)
    p_chromhmm.set_defaults(func=run_chromhmm)
    p_ctcf = sub.add_parser("ctcf", help="CTCF peak-class counts per Alu and window, from the CTCF universe")
    add_ctcf_args(p_ctcf)
    p_ctcf.set_defaults(func=run_ctcf)
    p_wgbs = sub.add_parser("wgbs", help="methylated CpG counts per Alu and window across all WGBS beds")
    add_wgbs_args(p_wgbs)
    p_wgbs.set_defaults(func=run_wgbs)
    args = parser.parse_args()
    args.func(args)


if __name__ == '__main__':
    main()
