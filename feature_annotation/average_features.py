#!/usr/bin/env python3
"""Average GC content and bigWig signal over Alus and their flanking windows."""

import os
import argparse
from pybedtools import BedTool, cleanup
import math
import numpy as np
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
RES_DIR = f'{PROJECTS_DIR}/alus/results/20251125_aluhg38_all_scores/'
DATA_DIR = f'{PROJECTS_DIR}/alus/data/'
UCSC_TRACKS=f'{POLLARD_DATA}/wynton/consortia/goldenPath/hg38/'
OUT_DIR=f'{PROJECT_DIR}/results/20260208_feature_overlap_all/'
CHROM_LENGTHS_PATH = '~/akita_variant_scoring/data/chrom_lengths_hg38'
CENTROMERE_PATH = '~/akita_variant_scoring/data/centromere_coords_hg38'
FASTA_PATH = f'{SUPREMO_DIR}/data/hg38.fa'
fasta_open = pysam.Fastafile(FASTA_PATH)
map_bw_36=f'{DATA_DIR}annotations/k36.Umap.MultiTrackMappability.bw'
map_bw_100=f'{DATA_DIR}annotations/k100.Umap.MultiTrackMappability.bw'
phyloP30=f'{UCSC_TRACKS}phyloP30way/hg38.phyloP30way.bw'
phyloP100=f'{UCSC_TRACKS}phyloP100way/hg38.phyloP100way.bw'
phastCon30=f'{UCSC_TRACKS}phastCons30way/hg38.phastCons30way.bw'
phastCon100=f'{UCSC_TRACKS}phastCons100way/hg38.phastCons100way.bw'
SEQ_LENGTH = 1048576
HALF_SEQ = SEQ_LENGTH // 2
WINDOW_SIZES = [1000, 10000, 100000]


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


def average_region_bw(chrom, start, end, map_bw, idx):
    """
    Get the average BigWig value in [start, end).
    """
    try:
        vals = map_bw.values(chrom, start, end, numpy=True)
        vals = vals[vals != None]
        if len(vals) == 0:
            return float('nan')
        return np.nanmean(vals)
    except (ValueError, RuntimeError) as e:
        print(f'BigWig error row {idx}: {e}')
        return float('nan')


def prepare_alu_windows(df, chrom_lengths_merge, size):
    df_region = df.copy()
    df_region = df_region.merge(chrom_lengths_merge, on='CHROM', how='left')
    df_region['alu_mid'] = ((df_region['POS'] + df_region['END']) // 2).astype(int)
    df_region['region_start'] = (df_region['alu_mid'] - size // 2).clip(lower=0)
    df_region['region_end'] = df_region['alu_mid'] + size // 2
    df_region['region_end'] = df_region[['region_end', 'chrom_length']].min(axis=1)
    return df_region


def average_region_bw__bigwig(chrom, grab_start, grab_end, map_bw, idx):
    try:
        # stats returns a list (one value per bin). Here we ask for 1 bin = mean over region
        m = map_bw.stats(chrom, int(grab_start), int(grab_end), type="mean")[0]
        if m is None:
            return float('nan')
        return m
        
    except (ValueError, RuntimeError) as e:
        print(f'BigWig error row {idx}: {e}')
        return float('nan')


def add_gc_args(parser):
        parser.add_argument("--chrom-lengths", default=CHROM_LENGTHS_PATH,
                            help="chromosome lengths table")
        parser.add_argument("--out-dir", default=OUT_DIR,
                            help="directory for output")
        parser.add_argument("--res-dir", default=RES_DIR,
                            help="directory holding the scored/annotated Alu table")
        parser.add_argument("--window-sizes", default=WINDOW_SIZES, type=int, nargs="+",
                            help="window sizes in bp")


def run_gc(args):
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
    print('calculating GC')
    alu_df_copy['GC_alu'] = alu_df_copy.apply(
        lambda row: calculate_gc(row['CHROM'], row['POS'], row['END'], row.name), axis=1)
    for size in window_sizes:
        print(f'  Window size: {size}')
        df_region = prepare_alu_windows(alu_df_copy, chrom_lengths_merge, size)
        #merge df_region with alu_df_copy
        df_region[f'GC_{size//1000}kb'] = df_region.apply(
        lambda row: calculate_gc(row['CHROM'], row['region_start'], row['region_end'], row.name), axis=1)
        alu_df_copy=alu_df_copy.merge(df_region[['orig_idx', f'GC_{size//1000}kb']])
    print(f'  1Mb window')
    alu_df_copy['GC_1Mb'] = alu_df_copy.apply(
        lambda row: calculate_gc(row['CHROM'], row['REF_start'], row['REF_stop'], row.name), axis=1)
    alu_df_copy.to_csv(f'{out_dir}/aluhg38_feature_GC.csv', index=False)
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
        'mapp36' : map_bw_36,
        'mapp100' : map_bw_100, 
        'phyloP30' : phyloP30,
        'phyloP100' : phyloP100, 
        'phastCon30' : phastCon30,
        'phastCon100' : phastCon100
    }
    for feat, feat_path in feature_dict.items():
        print(f'processing feature: {feat}')
        feat_bw =pyBigWig.open(feat_path)
        print(f'  alu window')
        alu_df_copy[f'{feat}_alu'] = alu_df_copy.apply(
            lambda row: average_region_bw__bigwig(row['CHROM'], row['POS'], row['END'], feat_bw, row.name), axis=1)

        for size in window_sizes:
            print(f'  Window size: {size}')
            df_region = prepare_alu_windows(alu_df_copy, chrom_lengths_merge, size)
            #merge df_region with alu_df_copy
            df_region[f'{feat}_{size//1000}kb'] = df_region.apply(
            lambda row: average_region_bw__bigwig(row['CHROM'], row['region_start'], row['region_end'], feat_bw, row.name), axis=1)
            alu_df_copy=alu_df_copy.merge(df_region[['orig_idx', f'{feat}_{size//1000}kb']])

        print(f'  1Mb window')
        alu_df_copy[f'{feat}_1Mb'] = alu_df_copy.apply(
            lambda row: average_region_bw__bigwig(row['CHROM'], row['REF_start'], row['REF_stop'], feat_bw, row.name), axis=1)
    
        #save intermediate df
        alu_df_copy.to_csv(f'{out_dir}/aluhg38_feature_averageBW.csv', index=False)
    print(alu_df_copy.head())
    alu_df_copy.to_csv(f'{out_dir}/aluhg38_feature_averageBW.csv', index=False)
    print("Done!")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    p_gc = sub.add_parser("gc", help="mean GC content of each Alu and its flanking windows")
    add_gc_args(p_gc)
    p_gc.set_defaults(func=run_gc)
    p_bigwig = sub.add_parser("bigwig", help="mean bigWig signal over each Alu and its flanking windows")
    add_bigwig_args(p_bigwig)
    p_bigwig.set_defaults(func=run_bigwig)
    args = parser.parse_args()
    args.func(args)


if __name__ == '__main__':
    main()
