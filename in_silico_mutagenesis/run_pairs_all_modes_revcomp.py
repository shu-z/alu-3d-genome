"""
run_del_pairs.py  —  unified deletion-pair scoring script (fwd + RC integrated)

Runs one or more modes in parallel, each scoring pairs of genomic deletions
and measuring their effect on Akita predictions.

Every sequence is predicted twice — once forward, once reverse-complemented —
and every MSE/CORR metric in the output appears in both _fwd and _rc forms.

Modes
─────
  alu_alu      : Alu1 (top-scoring) × nearby Alu2
  alu_highalu  : Alu1 × nearby high-scoring Alu2 (mse_mean > 0.0286)
  alu_random   : Alu1 × random non-Alu region (size-matched)
  alu_ctcf     : Alu1 × nearby CTCF site
  ctcf_ctcf    : CTCF1 (random 1000) × nearby CTCF2

Usage
─────
Run all modes in parallel (default):
    python run_del_pairs.py

Run specific modes:
    python run_del_pairs.py --modes alu_alu alu_ctcf

Run a single mode (no subprocess overhead):
    python run_del_pairs.py --modes alu_alu
"""

from concurrent.futures import ProcessPoolExecutor
import argparse
import math
import os
import sys
import time

import numpy as np
import pandas as pd
import pysam

# Paths - edit for your environment
SUPREMO_DIR = os.environ.get("SUPREMO_DIR", "/pollard/home/szhang20/akita_variant_scoring")
AKITA_DIR = os.environ.get("AKITA_DIR", "/pollard/home/szhang20/akita")
SLURM_DIR = os.environ.get("SLURM_DIR", "/pollard/home/szhang20/slurm")
PROJECT_DIR = os.environ.get("ALU_PROJECT_DIR", "/pollard/home/szhang20/alu")

MB = 1048576
n_samp=50

def revcomp(seq: str) -> str:
    comp = str.maketrans('ACGTacgtNn', 'TGCAtgcaNn')
    return seq.translate(comp)[::-1]

def fetch_delete_centered(chrom, alu_mid, deletions, hg38_fa, MB,
                          chrom_length=None, debug=False):
    half = MB // 2
    center_block = None
    for ds, de in deletions:
        if ds <= alu_mid < de:
            center_block = (ds, de)
            break
    center_coord = (center_block[0] + center_block[1]) // 2 if center_block else alu_mid

    ideal_start = center_coord - half
    ideal_end   = ideal_start + MB

    left_comp = right_comp = 0
    for ds, de in deletions:
        ov_s, ov_e = max(ds, ideal_start), min(de, ideal_end)
        if ov_s < ov_e:
            deleted_inside = ov_e - ov_s
            if ov_e <= center_coord:
                left_comp += deleted_inside
            elif ov_s >= center_coord:
                right_comp += deleted_inside
            else:
                left_comp  += center_coord - ov_s
                right_comp += ov_e - center_coord

    fetch_start = int(max(0, ideal_start - left_comp))
    fetch_end   = int(ideal_end + right_comp)
    if chrom_length and fetch_end > chrom_length:
        fetch_end = chrom_length

    seq      = hg38_fa.fetch(chrom, fetch_start, fetch_end)
    seq_list = list(seq)

    deleted_before_center = 0
    for ds, de in deletions:
        rel_s = max(0, ds - fetch_start)
        rel_e = min(len(seq_list), de - fetch_start)
        if rel_s < rel_e:
            for i in range(rel_s, rel_e):
                seq_list[i] = ""
            ov_left = max(0, min(de, center_coord) - max(ds, fetch_start))
            deleted_before_center += ov_left

    seq_after_del    = "".join(seq_list)
    center_idx_after = (center_coord - fetch_start) - deleted_before_center
    final_start      = center_idx_after - half
    final_end        = final_start + MB

    left_pad = right_pad = 0
    if final_start < 0:
        left_pad    = -final_start
        final_start = 0
        final_end   = MB
    if final_end > len(seq_after_del):
        right_pad = final_end - len(seq_after_del)
        final_end = len(seq_after_del)

    seq_final = seq_after_del[final_start:final_end]
    if left_pad:
        seq_final = "N" * left_pad + seq_final
    if right_pad:
        seq_final = seq_final + "N" * right_pad
    if len(seq_final) < MB:
        seq_final = seq_final + "N" * (MB - len(seq_final))
    elif len(seq_final) > MB:
        seq_final = seq_final[:MB]

    assert len(seq_final) == MB, f"Final seq length {len(seq_final)} != {MB}"
    return seq_final

def subtract_intervals(allowed, forbidden):
    result = []
    for a_start, a_end in allowed:
        current = [(a_start, a_end)]
        for f_start, f_end in forbidden:
            new_current = []
            for c_start, c_end in current:
                if f_end <= c_start or f_start >= c_end:
                    new_current.append((c_start, c_end))
                else:
                    if f_start > c_start:
                        new_current.append((c_start, f_start))
                    if f_end < c_end:
                        new_current.append((f_end, c_end))
            current = new_current
        result.extend(current)
    return result

def predict_both(seq, utils):
    """
    Predict contact map for `seq` in both forward and RC orientations.

    Returns
    -------
    fwd_mat : np.ndarray  — prediction on `seq` as-is
    rc_mat  : np.ndarray  — prediction on revcomp(seq)
    """
    fwd_vec = utils.vector_from_seq(seq)
    rc_vec  = utils.vector_from_seq(revcomp(seq))
    fwd_mat = utils.mat_from_vector(fwd_vec)
    rc_mat  = utils.mat_from_vector(rc_vec)
    return fwd_mat, rc_mat

def additivity_scores(REF_mat, ALT1_mat, ALT2_mat, ALT12_mat):
    """
    Test whether the double deletion effect is additive.

    Additivity means:  delta12 ≈ delta1 + delta2
    where delta = ALT - REF  (signed, per-bin difference)

    NaN positions (e.g. lower triangle, missing bins) are excluded
    from all summary statistics, consistent with the SuPreMo scoring approach.

    Returns
    -------
    mse_nonadd   : MSE(delta12, delta1+delta2); near 0 = additive
    mean_inter   : mean(delta12 - (delta1+delta2)); >0 synergistic, <0 buffering
    mean_delta1  : mean(ALT1 - REF)
    mean_delta2  : mean(ALT2 - REF)
    mean_delta12 : mean(ALT12 - REF)
    """
    delta1      = ALT1_mat  - REF_mat
    delta2      = ALT2_mat  - REF_mat
    delta12     = ALT12_mat - REF_mat
    interaction = delta12 - (delta1 + delta2)

    nan_mask = (
        np.isnan(REF_mat)   |
        np.isnan(ALT1_mat)  |
        np.isnan(ALT2_mat)  |
        np.isnan(ALT12_mat)
    )
    valid = ~nan_mask.reshape(-1)

    inter_flat   = interaction.reshape(-1)[valid]
    delta1_flat  = delta1.reshape(-1)[valid]
    delta2_flat  = delta2.reshape(-1)[valid]
    delta12_flat = delta12.reshape(-1)[valid]

    return (
        float(np.mean(inter_flat ** 2)),
        float(np.mean(inter_flat)),
        float(np.mean(delta1_flat)),
        float(np.mean(delta2_flat)),
        float(np.mean(delta12_flat)),
    )

def pair_metrics(REF_mat, ALT1_mat, ALT2_mat, ALT12_mat, scoring_utils,
                 names=('Alu1', 'Alu2')):
    """
    Compute all 6 pairwise MSE + 6 pairwise CORR + 5 additivity metrics
    for one orientation. Returned as a flat dict, keyed by the role names
    so the same function works for alu_alu, alu_random, alu_ctcf, ctcf_ctcf.

    `names` = (name1, name2), e.g. ('Alu1', 'Alu2') or ('Alu1', 'Random').
    """
    n1, n2 = names
    n12    = f'{n1}_{n2}'

    out = {
        f'MSE_{n1}_DEL':          scoring_utils.mse(REF_mat,  ALT1_mat),
        f'MSE_{n2}_DEL':          scoring_utils.mse(REF_mat,  ALT2_mat),
        f'MSE_{n12}_DEL':         scoring_utils.mse(REF_mat,  ALT12_mat),
        f'MSE_{n1}_v_{n2}':       scoring_utils.mse(ALT1_mat, ALT2_mat),
        f'MSE_{n1}_v_{n12}':      scoring_utils.mse(ALT1_mat, ALT12_mat),
        f'MSE_{n2}_v_{n12}':      scoring_utils.mse(ALT2_mat, ALT12_mat),
        f'CORR_{n1}_DEL':         scoring_utils.spearman(REF_mat,  ALT1_mat),
        f'CORR_{n2}_DEL':         scoring_utils.spearman(REF_mat,  ALT2_mat),
        f'CORR_{n12}_DEL':        scoring_utils.spearman(REF_mat,  ALT12_mat),
        f'CORR_{n1}_v_{n2}':      scoring_utils.spearman(ALT1_mat, ALT2_mat),
        f'CORR_{n1}_v_{n12}':     scoring_utils.spearman(ALT1_mat, ALT12_mat),
        f'CORR_{n2}_v_{n12}':     scoring_utils.spearman(ALT2_mat, ALT12_mat),
    }

    mse_nonadd, mean_inter, mean_d1, mean_d2, mean_d12 = additivity_scores(
        REF_mat, ALT1_mat, ALT2_mat, ALT12_mat)
    out.update({
        'MSE_nonadditive':       mse_nonadd,
        'mean_interaction':      mean_inter,
        f'mean_delta_{n1}':      mean_d1,
        f'mean_delta_{n2}':      mean_d2,
        f'mean_delta_{n12}':     mean_d12,
    })
    return out

def both_orient_metrics(REF_fwd, REF_rc, ALT1_fwd, ALT1_rc,
                        ALT2_fwd, ALT2_rc, ALT12_fwd, ALT12_rc,
                        scoring_utils, names=('Alu1', 'Alu2')):
    """Compute pair_metrics twice (fwd + rc) and merge with suffixes."""
    fwd = pair_metrics(REF_fwd, ALT1_fwd, ALT2_fwd, ALT12_fwd, scoring_utils, names)
    rc  = pair_metrics(REF_rc,  ALT1_rc,  ALT2_rc,  ALT12_rc,  scoring_utils, names)
    merged = {f'{k}_fwd': v for k, v in fwd.items()}
    merged.update({f'{k}_rc': v for k, v in rc.items()})
    return merged

################################################

def _run_alu_alu(cfg):
    sys.path.append(cfg['repo_path'])
    import akita_utils_forplotting as utils
    import akita_utils_scoring as scoring_utils

    hg38_fa = pysam.FastaFile(cfg['hg38_fa'])
    Alu     = pd.read_csv(cfg['alu_csv'],     sep='\t', index_col=0)
    alu_top = pd.read_csv(cfg['alu_top_csv'], sep='\t', index_col=0)
    results = []

    for idx, row in alu_top.iterrows():
        print(f'[alu_alu] {idx}')
        try:
            chrom, start, end = row['CHROM'], int(row['POS']), int(row['END'])
            Alu1_name = row['repName']
            alu1_mid  = math.ceil((start + end) / 2)

            REF_seq      = fetch_delete_centered(chrom, alu1_mid, [], hg38_fa, MB)
            Alu1_DEL_seq = fetch_delete_centered(chrom, alu1_mid, [(start, end)], hg38_fa, MB)

            REF_fwd,  REF_rc  = predict_both(REF_seq,      utils)
            Alu1_fwd, Alu1_rc = predict_both(Alu1_DEL_seq, utils)

            nearby_all = Alu[
                (Alu['CHROM'] == chrom) &
                (Alu['POS'] > start - cfg['nearby_window']) &
                (Alu['END']  < end   + cfg['nearby_window'])
            ]
            nearby_all = nearby_all[~((nearby_all['POS'] == start) & (nearby_all['END'] == end))]
            nearby = nearby_all.sample(n=min(n_samp, len(nearby_all)), random_state=0)

            for idx2, elem2 in nearby.iterrows():
                try:
                    e2s, e2e = int(elem2['POS']), int(elem2['END'])
                    Alu2_DEL_seq      = fetch_delete_centered(chrom, alu1_mid, [(e2s, e2e)], hg38_fa, MB)
                    Alu1_Alu2_DEL_seq = fetch_delete_centered(chrom, alu1_mid, [(start, end), (e2s, e2e)], hg38_fa, MB)

                    Alu2_fwd,      Alu2_rc      = predict_both(Alu2_DEL_seq,      utils)
                    Alu1_Alu2_fwd, Alu1_Alu2_rc = predict_both(Alu1_Alu2_DEL_seq, utils)

                    metrics = both_orient_metrics(
                        REF_fwd, REF_rc, Alu1_fwd, Alu1_rc,
                        Alu2_fwd, Alu2_rc, Alu1_Alu2_fwd, Alu1_Alu2_rc,
                        scoring_utils, names=('Alu1', 'Alu2'))

                    results.append({
                        'Alu1_idx': idx, 'Alu1_name': Alu1_name,
                        'Alu1_chrom': chrom, 'Alu1_start': start, 'Alu1_end': end,
                        'Alu2_idx': idx2, 'Alu2_name': elem2['repName'],
                        'Alu2_start': e2s, 'Alu2_end': e2e,
                        'n_Alu_nearby': len(nearby_all),
                        'nearby_window': cfg['nearby_window'],
                        **metrics,
                    })
                except Exception as e:
                    print(f'[alu_alu] error Alu2 {idx2}: {e}')
        except Exception as e:
            print(f'[alu_alu] error {idx}: {e}')

    pd.DataFrame(results).to_csv(cfg['out_alu_alu'], sep='\t')
    print(f'[alu_alu] done → {cfg["out_alu_alu"]}')

def _run_alu_highalu(cfg):
    sys.path.append(cfg['repo_path'])
    import akita_utils_forplotting as utils
    import akita_utils_scoring as scoring_utils

    hg38_fa = pysam.FastaFile(cfg['hg38_fa'])
    Alu_all = pd.read_csv(cfg['alu_csv'],     sep='\t', index_col=0)
    alu_top = pd.read_csv(cfg['alu_top_csv'], sep='\t', index_col=0)
    results = []

    # filter all Alu to high-scoring Alus
    Alu = Alu_all[Alu_all['mse_mean'] > 0.0286]

    for idx, row in alu_top.iterrows():
        print(f'[alu_highalu] {idx}')
        try:
            chrom, start, end = row['CHROM'], int(row['POS']), int(row['END'])
            Alu1_name = row['repName']
            alu1_mid  = math.ceil((start + end) / 2)

            REF_seq      = fetch_delete_centered(chrom, alu1_mid, [], hg38_fa, MB)
            Alu1_DEL_seq = fetch_delete_centered(chrom, alu1_mid, [(start, end)], hg38_fa, MB)

            REF_fwd,  REF_rc  = predict_both(REF_seq,      utils)
            Alu1_fwd, Alu1_rc = predict_both(Alu1_DEL_seq, utils)

            nearby_all = Alu[
                (Alu['CHROM'] == chrom) &
                (Alu['POS']  > start - cfg['nearby_window']) &
                (Alu['END']  < end   + cfg['nearby_window'])
            ]
            nearby_all = nearby_all[~((nearby_all['POS'] == start) & (nearby_all['END'] == end))]
            nearby = nearby_all.sample(n=min(n_samp, len(nearby_all)), random_state=0)

            for idx2, elem2 in nearby.iterrows():
                try:
                    e2s, e2e = int(elem2['POS']), int(elem2['END'])
                    Alu2_DEL_seq      = fetch_delete_centered(chrom, alu1_mid, [(e2s, e2e)], hg38_fa, MB)
                    Alu1_Alu2_DEL_seq = fetch_delete_centered(chrom, alu1_mid, [(start, end), (e2s, e2e)], hg38_fa, MB)

                    Alu2_fwd,      Alu2_rc      = predict_both(Alu2_DEL_seq,      utils)
                    Alu1_Alu2_fwd, Alu1_Alu2_rc = predict_both(Alu1_Alu2_DEL_seq, utils)

                    metrics = both_orient_metrics(
                        REF_fwd, REF_rc, Alu1_fwd, Alu1_rc,
                        Alu2_fwd, Alu2_rc, Alu1_Alu2_fwd, Alu1_Alu2_rc,
                        scoring_utils, names=('Alu1', 'Alu2'))

                    results.append({
                        'Alu1_idx': idx, 'Alu1_name': Alu1_name,
                        'Alu1_chrom': chrom, 'Alu1_start': start, 'Alu1_end': end,
                        'Alu2_idx': idx2, 'Alu2_name': elem2['repName'],
                        'Alu2_start': e2s, 'Alu2_end': e2e,
                        'n_Alu_nearby': len(nearby_all),
                        'nearby_window': cfg['nearby_window'],
                        **metrics,
                    })
                except Exception as e:
                    print(f'[alu_highalu] error Alu2 {idx2}: {e}')
        except Exception as e:
            print(f'[alu_highalu] error {idx}: {e}')

    pd.DataFrame(results).to_csv(cfg['out_alu_highalu'], sep='\t')
    print(f'[alu_highalu] done → {cfg["out_alu_highalu"]}')

def _run_alu_random(cfg):
    sys.path.append(cfg['repo_path'])
    import akita_utils_forplotting as utils
    import akita_utils_scoring as scoring_utils

    hg38_fa = pysam.FastaFile(cfg['hg38_fa'])
    Alu     = pd.read_csv(cfg['alu_csv'],     sep='\t', index_col=0)
    alu_top = pd.read_csv(cfg['alu_top_csv'], sep='\t', index_col=0)
    results = []

    for idx, row in alu_top.iterrows():
        print(f'[alu_random] {idx}')
        try:
            chrom, start, end = row['CHROM'], int(row['POS']), int(row['END'])
            Alu1_name = row['repName']
            alu1_mid  = math.ceil((start + end) / 2)

            REF_seq      = fetch_delete_centered(chrom, alu1_mid, [], hg38_fa, MB)
            Alu1_DEL_seq = fetch_delete_centered(chrom, alu1_mid, [(start, end)], hg38_fa, MB)

            REF_fwd,  REF_rc  = predict_both(REF_seq,      utils)
            Alu1_fwd, Alu1_rc = predict_both(Alu1_DEL_seq, utils)

            region_near = cfg['nearby_window']
            nearby_all  = Alu[
                (Alu['CHROM'] == chrom) &
                (Alu['POS'] > start - region_near) &
                (Alu['END']  < end   + region_near)
            ]
            sample_len    = end - start
            forbidden     = [(start, end)] + [(int(r['POS']), int(r['END'])) for _, r in nearby_all.iterrows()]
            allowed_clean = [
                (s, e) for s, e in subtract_intervals([(alu1_mid - region_near, alu1_mid + region_near)], forbidden)
                if (e - s) >= sample_len
            ]
            if not allowed_clean:
                print(f'[alu_random] no valid region for {idx}')
                continue

            for idx2 in range(n_samp):
                try:
                    interval     = allowed_clean[np.random.randint(len(allowed_clean))]
                    sample_start = np.random.randint(interval[0], interval[1] - sample_len + 1)
                    sample_end   = sample_start + sample_len

                    Rand_DEL_seq      = fetch_delete_centered(chrom, alu1_mid, [(sample_start, sample_end)], hg38_fa, MB)
                    Alu1_Rand_DEL_seq = fetch_delete_centered(chrom, alu1_mid, [(start, end), (sample_start, sample_end)], hg38_fa, MB)

                    Rand_fwd,      Rand_rc      = predict_both(Rand_DEL_seq,      utils)
                    Alu1_Rand_fwd, Alu1_Rand_rc = predict_both(Alu1_Rand_DEL_seq, utils)

                    metrics = both_orient_metrics(
                        REF_fwd, REF_rc, Alu1_fwd, Alu1_rc,
                        Rand_fwd, Rand_rc, Alu1_Rand_fwd, Alu1_Rand_rc,
                        scoring_utils, names=('Alu1', 'Random'))

                    results.append({
                        'Alu1_idx': idx, 'Alu1_name': Alu1_name,
                        'Alu1_chrom': chrom, 'Alu1_start': start, 'Alu1_end': end,
                        'Random_idx': idx2,
                        'Random_name': f'{chrom}_{sample_start}_{sample_end}',
                        'Random_start': sample_start, 'Random_end': sample_end,
                        'n_Alu_nearby': len(nearby_all),
                        'nearby_window': cfg['nearby_window'],
                        **metrics,
                    })
                except Exception as e:
                    print(f'[alu_random] error random {idx2}: {e}')
        except Exception as e:
            print(f'[alu_random] error {idx}: {e}')

    pd.DataFrame(results).to_csv(cfg['out_alu_random'], sep='\t')
    print(f'[alu_random] done → {cfg["out_alu_random"]}')

def _run_alu_ctcf(cfg):
    sys.path.append(cfg['repo_path'])
    import akita_utils_forplotting as utils
    import akita_utils_scoring as scoring_utils

    hg38_fa    = pysam.FastaFile(cfg['hg38_fa'])
    alu_top    = pd.read_csv(cfg['alu_top_csv'], sep='\t', index_col=0)
    CTCF_sites = _load_ctcf(cfg)
    results    = []

    for idx, row in alu_top.iterrows():
        print(f'[alu_ctcf] {idx}')
        try:
            chrom, start, end = row['CHROM'], int(row['POS']), int(row['END'])
            Alu1_name = row['repName']
            alu1_mid  = math.ceil((start + end) / 2)

            REF_seq      = fetch_delete_centered(chrom, alu1_mid, [], hg38_fa, MB)
            Alu1_DEL_seq = fetch_delete_centered(chrom, alu1_mid, [(start, end)], hg38_fa, MB)

            REF_fwd,  REF_rc  = predict_both(REF_seq,      utils)
            Alu1_fwd, Alu1_rc = predict_both(Alu1_DEL_seq, utils)

            nearby_all = CTCF_sites[
                (CTCF_sites['CHROM'] == chrom) &
                (CTCF_sites['POS'] > start - cfg['nearby_window']) &
                (CTCF_sites['END'] < end   + cfg['nearby_window'])
            ]
            nearby_all = nearby_all[~((nearby_all['POS'] < end) & (nearby_all['END'] > start))]
            nearby = nearby_all.sample(n=min(n_samp, len(nearby_all)), random_state=0)

            for idx2, elem2 in nearby.iterrows():
                try:
                    e2s, e2e = int(elem2['POS']), int(elem2['END'])
                    CTCF2_DEL_seq      = fetch_delete_centered(chrom, alu1_mid, [(e2s, e2e)], hg38_fa, MB)
                    Alu1_CTCF2_DEL_seq = fetch_delete_centered(chrom, alu1_mid, [(start, end), (e2s, e2e)], hg38_fa, MB)

                    CTCF2_fwd,      CTCF2_rc      = predict_both(CTCF2_DEL_seq,      utils)
                    Alu1_CTCF2_fwd, Alu1_CTCF2_rc = predict_both(Alu1_CTCF2_DEL_seq, utils)

                    metrics = both_orient_metrics(
                        REF_fwd, REF_rc, Alu1_fwd, Alu1_rc,
                        CTCF2_fwd, CTCF2_rc, Alu1_CTCF2_fwd, Alu1_CTCF2_rc,
                        scoring_utils, names=('Alu1', 'CTCF2'))

                    results.append({
                        'Alu1_idx': idx, 'Alu1_name': Alu1_name,
                        'Alu1_chrom': chrom, 'Alu1_start': start, 'Alu1_end': end,
                        'CTCF2_idx': idx2, 'CTCF2_name': elem2['CTCF_id'],
                        'CTCF2_start': e2s, 'CTCF2_end': e2e,
                        'n_CTCF_nearby': len(nearby_all),
                        'nearby_window': cfg['nearby_window'],
                        **metrics,
                    })
                except Exception as e:
                    print(f'[alu_ctcf] error CTCF2 {idx2}: {e}')
        except Exception as e:
            print(f'[alu_ctcf] error {idx}: {e}')

    pd.DataFrame(results).to_csv(cfg['out_alu_ctcf'], sep='\t')
    print(f'[alu_ctcf] done → {cfg["out_alu_ctcf"]}')

def _run_ctcf_ctcf(cfg):
    sys.path.append(cfg['repo_path'])
    import akita_utils_forplotting as utils
    import akita_utils_scoring as scoring_utils

    hg38_fa    = pysam.FastaFile(cfg['hg38_fa'])
    CTCF_sites = _load_ctcf(cfg)
    results    = []

    CTCF_to_test = CTCF_sites.sample(n=1000, random_state=729)

    for idx, row in CTCF_to_test.iterrows():
        print(f'[ctcf_ctcf] {idx}')
        try:
            chrom, start, end = row['CHROM'], int(row['POS']), int(row['END'])
            CTCF1_name = row['CTCF_id']
            ctcf1_mid  = math.ceil((start + end) / 2)

            REF_seq       = fetch_delete_centered(chrom, ctcf1_mid, [], hg38_fa, MB)
            CTCF1_DEL_seq = fetch_delete_centered(chrom, ctcf1_mid, [(start, end)], hg38_fa, MB)

            REF_fwd,   REF_rc   = predict_both(REF_seq,       utils)
            CTCF1_fwd, CTCF1_rc = predict_both(CTCF1_DEL_seq, utils)

            nearby_all = CTCF_sites[
                (CTCF_sites['CHROM'] == chrom) &
                (CTCF_sites['POS'] > start - cfg['nearby_window']) &
                (CTCF_sites['END'] < end   + cfg['nearby_window'])
            ]
            nearby_all = nearby_all[~((nearby_all['POS'] == start) & (nearby_all['END'] == end))]
            nearby = nearby_all.sample(n=min(n_samp, len(nearby_all)), random_state=0)

            for idx2, elem2 in nearby.iterrows():
                try:
                    e2s, e2e = int(elem2['POS']), int(elem2['END'])
                    CTCF2_DEL_seq       = fetch_delete_centered(chrom, ctcf1_mid, [(e2s, e2e)], hg38_fa, MB)
                    CTCF1_CTCF2_DEL_seq = fetch_delete_centered(chrom, ctcf1_mid, [(start, end), (e2s, e2e)], hg38_fa, MB)

                    CTCF2_fwd,       CTCF2_rc       = predict_both(CTCF2_DEL_seq,       utils)
                    CTCF1_CTCF2_fwd, CTCF1_CTCF2_rc = predict_both(CTCF1_CTCF2_DEL_seq, utils)

                    metrics = both_orient_metrics(
                        REF_fwd, REF_rc, CTCF1_fwd, CTCF1_rc,
                        CTCF2_fwd, CTCF2_rc, CTCF1_CTCF2_fwd, CTCF1_CTCF2_rc,
                        scoring_utils, names=('CTCF1', 'CTCF2'))

                    results.append({
                        'CTCF1_idx': idx, 'CTCF1_name': CTCF1_name,
                        'CTCF1_chrom': chrom, 'CTCF1_start': start, 'CTCF1_end': end,
                        'CTCF2_idx': idx2, 'CTCF2_name': elem2['CTCF_id'],
                        'CTCF2_start': e2s, 'CTCF2_end': e2e,
                        'n_CTCF_nearby': len(nearby_all),
                        'nearby_window': cfg['nearby_window'],
                        **metrics,
                    })
                except Exception as e:
                    print(f'[ctcf_ctcf] error CTCF2 {idx2}: {e}')
        except Exception as e:
            print(f'[ctcf_ctcf] error {idx}: {e}')

    pd.DataFrame(results).to_csv(cfg['out_ctcf_ctcf'], sep='\t')
    print(f'[ctcf_ctcf] done → {cfg["out_ctcf_ctcf"]}')

###### load ctcf 

def _load_ctcf(cfg):
    track_cols = ['CHROM', 'POS', 'END', 'tf', 'rel_score', '-log10(pval)', 'strand']
    chrom_list = [f'chr{i}' for i in range(1, 23)] + ['chrX']
    df = pd.read_csv(cfg['ctcf_csv'], sep='\t', header=None, names=track_cols)
    df = df[df['CHROM'].isin(chrom_list)]
    df['CTCF_id'] = df[['CHROM', 'POS', 'END']].astype(str).agg('_'.join, axis=1)
    return df

################################################

MODE_FN = {
    'alu_alu':     _run_alu_alu,
    'alu_highalu': _run_alu_highalu,
    'alu_random':  _run_alu_random,
    'alu_ctcf':    _run_alu_ctcf,
    #'ctcf_ctcf':   _run_ctcf_ctcf,
}

ALL_MODES = list(MODE_FN.keys())

def _run_mode(args):
    """Unpacks (mode, cfg) and dispatches — used by ProcessPoolExecutor."""
    mode, cfg = args
    t0 = time.time()
    MODE_FN[mode](cfg)
    print(f'[{mode}] elapsed {time.time()-t0:.0f}s')

def _parse():
    p = argparse.ArgumentParser(description='Unified deletion-pair scorer (fwd + RC)')
    p.add_argument('--modes', nargs='+', default=ALL_MODES, choices=ALL_MODES,
                   help='Which modes to run (default: all in parallel)')
    p.add_argument('--n_workers', type=int, default=None,
                   help='Max parallel processes (default: number of modes)')
    p.add_argument('--nearby_window', type=int, default=200000,
                   help='Search radius (bp) for nearby elements (default: 200000)')
    # path overrides
    p.add_argument('--repo_path',   default=f'{SUPREMO_DIR}/')
    p.add_argument('--hg38_fa',     default=f'{SUPREMO_DIR}/data/hg38.fa')
    p.add_argument('--alu_csv',     default=f'{PROJECT_DIR}/results/paper_results/20260316_merged_aluhg38_all_featureannot_NA_blacklist_filter_region_alu.txt')
    p.add_argument('--ctcf_csv',    default=f'{PROJECT_DIR}/data/MA0139.1.tsv')
    p.add_argument('--alu_top_csv', default=f'{PROJECT_DIR}/results/paper_results/aluhigh_1000samp_supremo.txt')
    p.add_argument('--res_dir',     default=f'{SLURM_DIR}/alu_pair_dels/')
    p.add_argument('--name',        default='alu_high')

    return p.parse_args()

if __name__ == '__main__':
    args = _parse()
    res  = args.res_dir.rstrip('/') + '/'
    os.makedirs(res, exist_ok=True)

    cfg = {
        'repo_path':     args.repo_path,
        'hg38_fa':       args.hg38_fa,
        'alu_csv':       args.alu_csv,
        'ctcf_csv':      args.ctcf_csv,
        'alu_top_csv':   args.alu_top_csv,
        'nearby_window': args.nearby_window,
        'out_alu_alu':     f'{res}{args.name}_alu_alu_DELpairs_fwdrc.txt',
        'out_alu_highalu': f'{res}{args.name}_alu_highalu_DELpairs_fwdrc.txt',
        'out_alu_random':  f'{res}{args.name}_alu_random_DELpairs_fwdrc.txt',
        'out_alu_ctcf':    f'{res}{args.name}_alu_ctcf_DELpairs_fwdrc.txt',
        'out_ctcf_ctcf':   f'{res}{args.name}_ctcf_ctcf_DELpairs_fwdrc.txt',
    }

    n_workers = args.n_workers or len(args.modes)
    pairs     = [(mode, cfg) for mode in args.modes]

    t0 = time.time()
    if len(pairs) == 1:
        _run_mode(pairs[0])
    else:
        with ProcessPoolExecutor(max_workers=n_workers) as ex:
            list(ex.map(_run_mode, pairs))

    print(f'\nAll done in {time.time()-t0:.0f}s')
