#!/usr/bin/env python3
"""Fractional overlap of bed tracks with Alus and their flanking windows."""

import sys
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from pybedtools import BedTool, cleanup
import math
import numpy as np
import os
import pandas as pd
import pysam
import re
import time
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "utils"))
from io_helpers import parse_attributes


SUPREMO_DIR = os.environ.get("SUPREMO_DIR", "/pollard/home/szhang20/akita_variant_scoring")
AKITA_DIR = os.environ.get("AKITA_DIR", "/pollard/home/szhang20/akita")
PROJECT_DIR = os.environ.get("ALU_PROJECT_DIR", "/pollard/home/szhang20/alu")
PROJECTS_DIR = os.environ.get("ALU_DATA_DIR", "/pollard/data/projects/shzhang")
POLLARD_DATA = os.environ.get("POLLARD_DATA", "/pollard/data")
DATA_DIR = f'{PROJECTS_DIR}/alus/data/'
RES_DIR = f'{PROJECTS_DIR}/alus/results/20251125_aluhg38_all_scores/'
OUT_DIR=f'{PROJECT_DIR}/results/20260208_feature_overlap_all/'
CHROM_LENGTHS_PATH = '~/akita_variant_scoring/data/chrom_lengths_hg38'
FASTA_PATH = f'{SUPREMO_DIR}/data/hg38.fa'
GTF_PATH = f'{DATA_DIR}gencode/gencode.v29.annotation.gtf.gz'
REPEAT_PATH = f'{DATA_DIR}repeats/hg38.repeatmasker.txt'
H3K9_PATH = f'{DATA_DIR}encode/H1_H3K9me3_ENCFF310CZV.bed.gz'
H3K27_PATH = f'{DATA_DIR}encode/H1_H3K27me3_ENCFF254ACI.bed.gz'
SEQ_LENGTH = 1048576
HALF_SEQ = SEQ_LENGTH // 2
WINDOW_SIZES = [1000, 10000, 100000]
RES_DIR__CHIP = f'{PROJECT_DIR}/results/paper_results/'
OUT_DIR__CHIP = RES_DIR__CHIP
ENCODE_DIR = f'{PROJECTS_DIR}/alus/data/encode/ChIP/'  # directory with .bed.gz files
CHROM_LENGTHS_PATH__CHIP = f'{SUPREMO_DIR}/data/chrom_lengths_hg38'
WINDOW_SIZES__CHIP = [10000, 100000]  # 10kb, 100kb (1Mb handled separately)
N_WORKERS = 24  # adjust to your cluster




def read_gencode(gtf_path: str) -> pd.DataFrame:
    df = pd.read_csv(gtf_path, sep='\t', comment='#', header=None, compression='gzip', low_memory=False,
                     names=["chrom", "source", "feature", "start", "end", "score", "strand", "score2", "info"])
    df_attrs = df['info'].apply(parse_attributes)
    df = pd.concat([df, pd.DataFrame(df_attrs.tolist())], axis=1)
    df['gene_len'] = df['end'] - df['start']
    return df


def compute_fraction_overlap(alu_df: pd.DataFrame, feature_BED: BedTool, colname: str) -> pd.DataFrame:
    df = alu_df.copy().dropna(subset=['region_start','region_end'])
    df['region_start'] = df['region_start'].astype(int)
    df['region_end'] = df['region_end'].astype(int)

    alu_bed = BedTool.from_dataframe(df)
    feature_BED = feature_BED.sort().merge()
    inter = alu_bed.intersect(feature_BED, wo=True)

    overlap_totals = {}
    for entry in inter:
        #print(entry)
        #print('alu_id', entry[3])
        alu_id = int(entry[3])

        #print('alu_id', alu_id)

        overlap_len = int(entry[-1])
        overlap_totals[alu_id] = overlap_totals.get(alu_id, 0) + overlap_len
        # except:
        #     print(f"ERROR!!!!!! {entry}")

    df[colname] = df.apply(
        lambda row: overlap_totals.get(row['orig_idx'],0)/(row['region_end']-row['region_start']) 
        if row['region_end']-row['region_start']>0 else 0.0, axis=1
    )
    return df[['orig_idx', colname]]


def compute_fraction_overlap_bychrom(alu_df: pd.DataFrame, feature_BED: BedTool, colname: str) -> pd.DataFrame:
    # Prepare ALU dataframe
    df = alu_df.copy().dropna(subset=['region_start','region_end'])
    df['region_start'] = df['region_start'].astype(int)
    df['region_end'] = df['region_end'].astype(int)
    df = df.copy().dropna(subset=['region_start','region_end'])

    # Sort and merge features once
    #feature_BED = feature_BED.sort().merge()

    # Dictionary to store overlap per ALU
    overlap_totals = {}

    # Iterate chromosome by chromosome
    chrom_list=['chr' + str(i) for i in range(1,23)] + ['chrX', 'chrY']
    #chrom_list=['chr9']
    #for chrom in df['chrom'].unique():
    for chrom in chrom_list:
        print(f'    {chrom} ... ')
        df_chrom = df[df['CHROM'] == chrom]
        if df_chrom.empty:
            continue

        alu_bed_chrom = BedTool.from_dataframe(df_chrom)
        chrom_feat=feature_BED[feature_BED['genoName']==chrom]

        feature_chrom=BedTool.from_dataframe(chrom_feat)
        feature_chrom = feature_chrom.sort().merge()

        #feature_chrom = feature_BED.filter(lambda x: x.chrom == chrom).saveas()
        inter = alu_bed_chrom.intersect(feature_chrom, wo=True)

        for entry in inter:
            if len(entry) == 0:
                print('entry is empty ', entry)
                continue

            overlap_str = entry[-1]
            alu_id_str=entry[3]

            if overlap_str is None or overlap_str == "":
                print("Empty overlap field:", entry)
                continue

            try:
                overlap_len = int(overlap_str)
            except ValueError:
                print("Non-integer overlap field:", overlap_str, "in entry:", entry)
                continue

            try:
                alu_id = int(alu_id_str)
            except ValueError:
                print("Non-integer overlap field:", alu_id_str, "in entry:", entry)
                continue

            alu_id = int(entry[3])
            overlap_len = int(entry[-1])
            overlap_totals[alu_id] = overlap_totals.get(alu_id, 0) + overlap_len
        
        cleanup() #cleans up pybedtools temp files 

    # Compute fraction overlap
    df[colname] = df.apply(
        lambda row: overlap_totals.get(row['orig_idx'], 0) / (row['region_end'] - row['region_start'])
        if row['region_end'] - row['region_start'] > 0 else 0.0,
        axis=1
    )

    return df[['orig_idx', colname]]


def prepare_alu_windows(alu_df: pd.DataFrame, chrom_lengths: pd.DataFrame, size: int) -> pd.DataFrame:
    df = alu_df.copy()
    df = df.merge(chrom_lengths, on='CHROM', how='left')
    df['alu_mid'] = (df['POS'] + df['END']) // 2
    df['region_start'] = (df['alu_mid'] - size//2).clip(lower=0)
    df['region_end'] = df[['alu_mid','chrom_length']].min(axis=1) + (size//2)
    return df[['CHROM','region_start','region_end','orig_idx']]


def compute_fraction_overlap_bychrom__chip(alu_df: pd.DataFrame, feature_df: pd.DataFrame, colname: str) -> pd.DataFrame:
    df = alu_df.copy().dropna(subset=['region_start', 'region_end'])
    df['region_start'] = df['region_start'].astype(int)
    df['region_end'] = df['region_end'].astype(int)

    overlap_totals = {}
    chrom_list = ['chr' + str(i) for i in range(1, 23)] + ['chrX', 'chrY']

    for chrom in chrom_list:
        df_chrom = df[df['CHROM'] == chrom]
        if df_chrom.empty:
            continue

        chrom_feat = feature_df[feature_df.iloc[:, 0] == chrom]
        if chrom_feat.empty:
            continue

        alu_bed_chrom = BedTool.from_dataframe(df_chrom)
        feature_chrom = BedTool.from_dataframe(chrom_feat).sort().merge()
        inter = alu_bed_chrom.intersect(feature_chrom, wo=True)

        for entry in inter:
            if len(entry) == 0:
                continue
            try:
                alu_id = int(entry[3])
                overlap_len = int(entry[-1])
            except ValueError:
                continue
            overlap_totals[alu_id] = overlap_totals.get(alu_id, 0) + overlap_len

        cleanup()

    df[colname] = df.apply(
        lambda row: overlap_totals.get(row['orig_idx'], 0) / (row['region_end'] - row['region_start'])
        if row['region_end'] - row['region_start'] > 0 else 0.0,
        axis=1
    )
    return df[['orig_idx', colname]]


def prepare_alu_windows__chip(alu_df: pd.DataFrame, chrom_lengths: pd.DataFrame, size: int) -> pd.DataFrame:
    df = alu_df.copy()
    df = df.merge(chrom_lengths, on='CHROM', how='left')
    df['alu_mid'] = (df['POS'] + df['END']) // 2
    df['region_start'] = (df['alu_mid'] - size // 2).clip(lower=0)
    df['region_end'] = df[['alu_mid', 'chrom_length']].min(axis=1) + (size // 2)
    return df[['CHROM', 'region_start', 'region_end', 'orig_idx']]


def process_one_file(bed_path: str, alu_score_df: pd.DataFrame, chrom_lengths_merge: pd.DataFrame) -> pd.DataFrame:
    """Process a single ENCODE bed.gz file — all windows. Returns merged df with orig_idx."""
    bed_path = Path(bed_path)
    feat_name = bed_path.stem.replace('.bed', '')  # e.g. ENCFF643LUR
    print(f'[{feat_name}] Loading...')

    try:
        feat_df = pd.read_csv(bed_path, sep='\t', header=None, compression='gzip')
    except Exception as e:
        print(f'[{feat_name}] Failed to load: {e}')
        return None

    result_df = alu_score_df[['orig_idx']].copy()

    # 10kb and 100kb windows
    for size in WINDOW_SIZES__CHIP:
        colname = f'{feat_name}_overlap_{size // 1000}kb'
        print(f'[{feat_name}] Window {size // 1000}kb...')
        df_region = prepare_alu_windows__chip(alu_score_df, chrom_lengths_merge, size)
        df_feat = compute_fraction_overlap_bychrom__chip(df_region, feat_df, colname)
        result_df = result_df.merge(df_feat, on='orig_idx')

    # 1Mb window
    colname_1mb = f'{feat_name}_overlap_1Mb'
    print(f'[{feat_name}] Window 1Mb...')
    df_window = alu_score_df[['CHROM', 'REF_start', 'REF_stop', 'orig_idx']].copy()
    df_window.columns = ['CHROM', 'region_start', 'region_end', 'orig_idx']
    df_window = df_window.merge(chrom_lengths_merge, on='CHROM', how='left')
    df_window['region_start'] = df_window['region_start'].clip(lower=0)
    df_window['region_end'] = df_window[['region_end', 'chrom_length']].min(axis=1)
    df_feat_1mb = compute_fraction_overlap_bychrom__chip(df_window, feat_df, colname_1mb)
    result_df = result_df.merge(df_feat_1mb, on='orig_idx')

    print(f'[{feat_name}] Done.')
    return result_df


def add_tracks_args(parser):
        parser.add_argument("--chrom-lengths", default=CHROM_LENGTHS_PATH,
                            help="chromosome lengths table")
        parser.add_argument("--fasta", default=FASTA_PATH,
                            help="reference FASTA")
        parser.add_argument("--h3k27", default=H3K27_PATH,
                            help="H3K27me3 peak bed")
        parser.add_argument("--h3k9", default=H3K9_PATH,
                            help="H3K9me3 peak bed")
        parser.add_argument("--out-dir", default=OUT_DIR,
                            help="directory for output")
        parser.add_argument("--repeats", default=REPEAT_PATH,
                            help="repeatmasker table")
        parser.add_argument("--res-dir", default=RES_DIR,
                            help="directory holding the scored/annotated Alu table")
        parser.add_argument("--window-sizes", default=WINDOW_SIZES, type=int, nargs="+",
                            help="window sizes in bp")


def run_tracks(args):
    chrom_lengths_path = args.chrom_lengths
    fasta_path = args.fasta
    h3k27_path = args.h3k27
    h3k9_path = args.h3k9
    out_dir = args.out_dir
    repeat_path = args.repeats
    res_dir = args.res_dir
    window_sizes = args.window_sizes
    chrom_lengths = pd.read_table(chrom_lengths_path, header=None, names=['CHROM','chrom_max'])
    chrom_lengths_merge = chrom_lengths.copy()
    chrom_lengths_merge.columns = ['CHROM','chrom_length']
    chrom_lengths_merge['CHROM'] = 'chr' + chrom_lengths_merge['CHROM'].astype(str)
    fasta_open = pysam.Fastafile(fasta_path)
    repeat_all = pd.read_csv(repeat_path, sep='\t')
    feature_dict={}
    feature_dict['HFF_H3K9me3'] = pd.read_csv(h3k9_path, sep='\t', header=None)
    feature_dict['HFF_H3K27me3'] = pd.read_csv(h3k27_path, sep='\t', header=None)
    alu_annot = pd.read_csv(f'{res_dir}20260123_aluhg38_all_annot_GC_mapp_alu_repeats_genetrackoverlap.txt',
                             sep='\t', index_col=0)
    alu_score_df = alu_annot[['CHROM', 'POS', 'END','orig_idx','REF_start','REF_stop']]
    alu_df_copy=alu_annot[['CHROM', 'POS', 'END','orig_idx']]
    for feat, feat_df in feature_dict.items():
        print(f'processing feature: {feat}')
        feat_BED = BedTool.from_dataframe(feat_df)

        # #alu overlap
        print('alu window')
        df_window= alu_score_df[['CHROM','POS','END','orig_idx']]
        df_window.columns=['CHROM', 'region_start', 'region_end', 'orig_idx']
        #df_feat_window = compute_fraction_overlap_bychrom(df_window, feat_df, f'{feat}_overlap_alu')

        #if using non chrom one 
        df_feat_window = compute_fraction_overlap(df_window, feat_BED, f'{feat}_overlap_alu')

        alu_df_copy = alu_df_copy.merge(df_feat_window, on='orig_idx')

        for size in window_sizes:
            print(f'  Window size: {size}')
            df_region = prepare_alu_windows(alu_score_df, chrom_lengths_merge, size)
            colname = f'{feat}_overlap_{size//1000}kb'
            #df_feat = compute_fraction_overlap_bychrom(df_region, feat_df, colname)
            df_feat = compute_fraction_overlap(df_region, feat_BED, colname)

            alu_df_copy = alu_df_copy.merge(df_feat, on='orig_idx')

        #window overlap
        print('1Mb window')
        df_window = alu_score_df[['CHROM','REF_start','REF_stop','orig_idx']]
        df_window.columns=['CHROM', 'region_start', 'region_end', 'orig_idx']
        #add in clipping to see if this fixes things 
        df_window = df_window.merge(chrom_lengths_merge, on='CHROM', how='left')
        df_window['region_start'] = (df_window['region_start']).clip(lower=0)
        df_window['region_end'] = df_window[['region_end','chrom_length']].min(axis=1) 

        #df_feat_window = compute_fraction_overlap_bychrom(df_window, feat_df, f'{feat}_overlap_1Mb')
        df_feat_window = compute_fraction_overlap(df_window, feat_BED, f'{feat}_overlap_1Mb')

        alu_df_copy = alu_df_copy.merge(df_feat_window, on='orig_idx')
        cleanup() #cleans up pybedtools temp files 
        
        #save intermediate copy
        alu_df_copy.to_csv(f'{out_dir}/aluhg38_feature_overlap_H1_histone.csv', index=False)
    alu_df_copy.to_csv(f'{out_dir}/aluhg38_feature_overlap_H1_histone.csv', index=False)
    print("Done!")


def add_chip_args(parser):
        parser.add_argument("--chrom-lengths", default=CHROM_LENGTHS_PATH__CHIP,
                            help="chromosome lengths table")
        parser.add_argument("--encode-dir", default=ENCODE_DIR,
                            help="directory of ENCODE ChIP beds")
        parser.add_argument("--workers", default=N_WORKERS, type=int,
                            help="parallel worker processes")
        parser.add_argument("--out-dir", default=OUT_DIR__CHIP,
                            help="directory for output")
        parser.add_argument("--res-dir", default=RES_DIR__CHIP,
                            help="directory holding the scored/annotated Alu table")


def run_chip(args):
    chrom_lengths_path = args.chrom_lengths
    encode_dir = args.encode_dir
    n_workers = args.workers
    out_dir = args.out_dir
    res_dir = args.res_dir
    chrom_lengths_merge = pd.read_table(chrom_lengths_path, header=None, names=['CHROM', 'chrom_length'])
    chrom_lengths_merge['CHROM'] = 'chr' + chrom_lengths_merge['CHROM'].astype(str)
    alu_annot = pd.read_csv(f'{res_dir}20260316_merged_aluhg38_all_featureannot_NA_blacklist_filter_region_alu.txt',
        sep='\t', index_col=0)
    alu_score_df = alu_annot[['CHROM', 'POS', 'END', 'orig_idx', 'REF_start', 'REF_stop']]
    alu_df_copy = alu_annot[['CHROM', 'POS', 'END', 'orig_idx']].copy()
    encode_files = sorted(Path(encode_dir).glob('ENCF*.bed.gz'))
    print(f'Found {len(encode_files)} ENCODE bed.gz files.')
    results = []
    with ProcessPoolExecutor(max_workers=n_workers) as executor:
        futures = {
            executor.submit(process_one_file, str(f), alu_score_df, chrom_lengths_merge): f.name
            for f in encode_files
        }
        for future in as_completed(futures):
            fname = futures[future]
            try:
                res = future.result()
                if res is not None:
                    results.append(res)
                    print(f'Completed: {fname}')
            except Exception as e:
                print(f'ERROR processing {fname}: {e}')
    for res in results:
        alu_df_copy = alu_df_copy.merge(res, on='orig_idx', how='left')
    out_path = f'{out_dir}/aluhg38_feature_overlap_encode_ChIP_TF.csv'
    alu_df_copy.to_csv(out_path, index=False)
    print(f'Saved to {out_path}')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    p_tracks = sub.add_parser("tracks", help="fractional overlap with histone, repeat and gencode bed tracks")
    add_tracks_args(p_tracks)
    p_tracks.set_defaults(func=run_tracks)
    p_chip = sub.add_parser("chip", help="fractional overlap with every ENCODE ChIP bed in a directory, in parallel")
    add_chip_args(p_chip)
    p_chip.set_defaults(func=run_chip)
    args = parser.parse_args()
    args.func(args)


if __name__ == '__main__':
    main()
