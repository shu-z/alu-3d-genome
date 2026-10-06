"""Converted from 20260303_alu_insertion_gradients.ipynb; logic unchanged."""

from __future__ import annotations

import math, pysam
import os
import subprocess
import sys

from cooltools.lib.plotting import gridspec_inches
import h5py
import matplotlib as mpl
import matplotlib.backends.backend_pdf as pdf_backend
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyBigWig

import akita_utils_scoring    as _su
import get_Akita_scores_utils as _akita
import plotting_utils         as _pu
import scoring                as _sc
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "utils"))
from io_helpers import save_figures_to_pdf

# Paths - edit for your environment
SUPREMO_DIR = os.environ.get("SUPREMO_DIR", "/pollard/home/szhang20/akita_variant_scoring")
AKITA_DIR = os.environ.get("AKITA_DIR", "/pollard/home/szhang20/akita")
PROJECT_DIR = os.environ.get("ALU_PROJECT_DIR", "/pollard/home/szhang20/alu")
DATA_DIR = os.environ.get("ALU_DATA_DIR", "/pollard/data/projects/shzhang")
POLLARD_DATA = os.environ.get("POLLARD_DATA", "/pollard/data")

# %load_ext autoreload
# %autoreload 2

os.chdir("~/akita_variant_scoring/")


#what we want to do, cohesively

#first, take a sampling of top... 1000 Alus in top 1%? 
#filter to full length Alus 
#insert every single consensus sequence into them? -- ignore FAM, FRAM, FLAM
#also insert other top alus?.. -- do we want to select based on some criteria?
    #repeatmasker annotations like millidiv, ins etc don't differ between alus
#want to get MSE diff between REF (Alu DEL) and each INS
#want to get MSE diff between each INS vs orig Alu
#want to get saliency track for REF (Alu DEL), orig Alu'

#will need to parallelize

#to quantify changes in saliency -- get REF-ALT, pick bins that have an absolute change? or fold change?

#finally, overlap bins with largest changes with our universal CTCF set 
    #questions to answer -- is saliency changing more at more cell-type specific CTCFs? 
    #how often is saliency changing at CTCF sites vs other regions? and what are the other regions?

#read in alus with overlap scores 
res_dir = f'{PROJECT_DIR}/results/paper_results/'
ins_grad_dir = f'{PROJECT_DIR}/bin/ins_gradients/'

scores_top = pd.read_csv(f'{res_dir}20260414_highlow1000_insertconsensus_withDEL.txt',sep='\t', index_col=0)
scores_top

#indices we want 
query_indices_no_PCDH = [
     330025, 971709, 308144, 307931, 894177, 308140, 775656, 380293, 989927, 973796, 979684,
     808325, 1110696, 1509, 594428, 278777, 592546, 1065376, 308135,
     330106, 971465, 537341, 577482, 946283, 615308, 904830, 308141, 784414,
    779116, 881832, 278761, 1084551, 924259, 971880, 308160, 1019933, 971750, 888842,
    575044, 65537, 379896, 256966, 876512, 139045, 923647, 1142121, 921284, 393618,
    923888, 91291, 24945, 65443, 953723, 1061595, 308125, 1071679, 594387, 1098967,
    1110842, 946861, 60696, 893629, 1093129, 294719, 946539, 1104964, 989951, 65977,
    946044, 1100941, 375619, 971602, 779207, 779133, 957880, 1058710, 960685,
    1132753, 1078157, 1149088, 66016, 946509, 184472, 615340, 808894, 1055987, 666378,
]

scores_query=scores_top[scores_top['orig_idx'].isin(query_indices_no_PCDH)]
scores_query

scores_query.to_csv(f'{ins_grad_dir}20260519_topMSE/top_ins_noPDCH.txt', sep='\t')

# ### Results and plot!

MM_PER_IN = 25.4
FIG_W_MM  = 120
FIG_W_IN  = FIG_W_MM / MM_PER_IN          # ≈ 4.72"
ROW_H_IN  = 0.55                          # per-track height; tweak to taste
 
# ggplot2-ish palette (still ggplot for colors, but classic for theme)
COL_REF   = '#3B6FB6'   # blue
COL_ALT   = '#C0504D'   # brick
COL_INS   = '#7F7F7F'   # grey for inserted Alu region
COL_CTCF  = '#4E8B57'   # green
COL_DELTA = '#8064A2'
 
# CTCF peak-class colors (Brewer Set2-ish)
PEAK_CLASS_COLORS = {
    'cell_type_specific' : '#B33E52',
    'lineage_restricted' : '#653EB3',
    'broadly_shared'     : '#0F8299',
    'constitutive'       : '#54990F',  
    'singleton'          : '#999999'
}
 
PLOT_RC = {
    'font.family'        : 'sans-serif',
    'font.sans-serif'    : ['Helvetica', 'Arial', 'DejaVu Sans'],
    'font.size'          : 7,
    'axes.titlesize'     : 8,
    'axes.labelsize'     : 7,
    'xtick.labelsize'    : 6,
    'ytick.labelsize'    : 6,
    'legend.fontsize'    : 6,
    # theme_classic: white background, no grid, only bottom/left spines
    'axes.linewidth'     : 0.6,
    'axes.edgecolor'     : 'black',
    'axes.facecolor'     : 'white',
    'figure.facecolor'   : 'white',
    'axes.grid'          : False,
    'axes.spines.top'    : False,
    'axes.spines.right'  : False,
    'axes.spines.left'   : True,
    'axes.spines.bottom' : True,
    'xtick.direction'    : 'out',
    'ytick.direction'    : 'out',
    'xtick.major.size'   : 2.5,
    'ytick.major.size'   : 2.5,
    'xtick.major.width'  : 0.5,
    'ytick.major.width'  : 0.5,
    'lines.linewidth'    : 0.9,
    'figure.dpi'         : 150,
    'savefig.dpi'        : 600,
    'savefig.bbox'       : 'tight',
    'pdf.fonttype'       : 42,
    'ps.fonttype'        : 42,
    'svg.fonttype'       : 'none',
}

mpl.rcParams.update(PLOT_RC)

BIN_SIZE  = 256
CROP_BINS = 32   # Akita edge crop = 32 * 2048 nt = 65536 bp each side
OUT_BINS   = 448

####### helper functions

def rolling_mean(arr: np.ndarray, window: int) -> np.ndarray:
    """Centred rolling mean, same-length output. NaN-safe."""
    if window >= len(arr):
        return np.full_like(arr, np.nanmean(arr))
    return np.convolve(arr, np.ones(window) / window, mode='same')

BIN_SIZE  = 512
CROP_BINS = 32   # Akita edge crop = 32 * 2048 nt = 65536 bp each side

def nt_to_saliency_binned(grads: np.ndarray,
                           bin_size: int = BIN_SIZE,
                           crop_bins: int = CROP_BINS) -> np.ndarray:
    """
    [seq_len, 4] float16 → binned saliency, with Akita edge crop applied.
    seq_len=2^20, bin_size=512 → 2048 raw bins → crop 32 each side → 1984 bins.
    """
    sal     = np.sum(np.abs(grads.astype(np.float32)), axis=1)
    L       = len(sal) - (len(sal) % bin_size)
    binned  = sal[:L].reshape(-1, bin_size).mean(axis=1)
    cropped = binned[crop_bins : len(binned) - crop_bins]
    print(f'  nt_to_saliency_binned: nt={len(sal)} → {len(binned)} bins '
          f'→ {len(cropped)} after crop')
    return cropped   # [1984]

def align_and_crop(REF_grads: np.ndarray, ALT_grads: np.ndarray,consensus_len_bp: int,
                   bin_size: int = BIN_SIZE, crop_bins: int = CROP_BINS, akita_bin_size: int = 2048
                   ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    REF and ALT are both 2^20 nt, same length.
    ALT contains a consensus_len insertion at its center — so ALT's flanks
    are each compressed by consensus_len//2 bp relative to REF.

    Comparable region (after Akita crop):
      - ALT : excise the center consensus_len nt, keep the flanks
      - REF : crop consensus_len//2 nt from each end to match ALT's flanks

    All of this happens AFTER the 32-Akita-bin (32*2048 nt) edge crop.

    Returns
    -------
    ref_out    : [n_comparable_bins]  REF saliency over comparable region
    alt_out    : [n_comparable_bins]  ALT saliency, insertion excised
    centers_bp : [n_comparable_bins]  bp relative to insert_pos
    alt_insert : [ins_bins]           ALT saliency within the inserted sequence
    """
    akita_crop_nt = crop_bins * akita_bin_size   # 65536 nt each side

    #### flatten to nt saliecny 
    ref_sal = np.sum(np.abs(REF_grads.astype(np.float32)), axis=1)   # [2^20]
    alt_sal = np.sum(np.abs(ALT_grads.astype(np.float32)), axis=1)   # [2^20]
    seq_len = len(ref_sal)   # 1048576

    ### apply akita crop 
    ref_cropped = ref_sal[akita_crop_nt : seq_len - akita_crop_nt]   # [917504 nt]
    alt_cropped = alt_sal[akita_crop_nt : seq_len - akita_crop_nt]   # [917504 nt]
    cropped_len = len(ref_cropped)   # seq_len - 2*akita_crop_nt

    # ── 3. locate insertion in cropped ALT
    ins_half  = consensus_len_bp // 2
    alt_mid   = cropped_len // 2
    ins_start = alt_mid - ins_half
    ins_end   = ins_start + consensus_len_bp   # use full consensus_len, not 2*ins_half

    alt_insert_nt = alt_cropped[ins_start : ins_end]

    # ── 4. excise insertion from ALT; crop REF by exact same amount ───────
    alt_comparable = np.concatenate([alt_cropped[:ins_start], alt_cropped[ins_end:]])

    actual_ins_len = ins_end - ins_start   # == consensus_len_bp, but explicit
    crop_each_side = actual_ins_len // 2
    # if odd, take the extra nt from the right crop so both are same length
    crop_left  = crop_each_side
    crop_right = actual_ins_len - crop_each_side
    ref_comparable = ref_cropped[crop_left : cropped_len - crop_right]

    if len(ref_comparable) != len(alt_comparable):
        raise ValueError(f'Comparable length mismatch: '
                         f'ref={len(ref_comparable)}, alt={len(alt_comparable)}, '
                         f'consensus_len={consensus_len_bp}')

    # ── 5. bin both to bin_size 
    def _bin(arr):
        L = len(arr) - (len(arr) % bin_size)
        #return arr[:L].reshape(-1, bin_size).mean(axis=1)
        return arr[:L].reshape(-1, bin_size).max(axis=1)

    ref_out       = _bin(ref_comparable)
    alt_out       = _bin(alt_comparable)
    alt_insert    = _bin(alt_insert_nt)

    # ── 6. coordinate axis: bp relative to insert_pos 
    # comparable region spans cropped_len - consensus_len_bp nt,
    # centered on insert_pos
    half_comparable_nt = len(ref_comparable) // 2
    centers_bp = (np.arange(len(ref_out)) + 0.5) * bin_size - half_comparable_nt

    print(f'  align_and_crop: seq={seq_len}, akita_crop={akita_crop_nt}, '
          f'ins=[{ins_start}:{ins_end}] ({consensus_len_bp} nt), '
          f'comparable={len(ref_comparable)} nt → {len(ref_out)} bins')

    assert len(ref_out) == len(alt_out), \
        f'bin mismatch: ref={len(ref_out)}, alt={len(alt_out)}'

    return ref_out, alt_out, centers_bp, alt_insert

def bin_ctcf(chrom: str, insert_pos: int,
             consensus_len_bp: int,
             bin_size: int = BIN_SIZE,
             crop_bins: int = CROP_BINS,
             akita_bin_size: int = 2048) -> tuple[np.ndarray, np.ndarray]:
    """
    Fetch CTCF over the same comparable window as align_and_crop:
      full window  = 2^20 nt centered on insert_pos
      akita crop   = 32 * 2048 nt each side
      then drop ins_half from each end to match REF comparable region
    """
    seq_len       = 1 << 20
    akita_crop_nt = crop_bins * akita_bin_size    # 65536
    ins_half      = consensus_len_bp // 2

    # comparable region in genomic coords
    comp_half = seq_len // 2 - akita_crop_nt - ins_half
    start     = insert_pos - comp_half
    end       = insert_pos + comp_half

    try:
        with pyBigWig.open(CTCF_BW) as bw:
            vals = np.nan_to_num(np.array(bw.values(chrom, start, end)))
    except Exception as e:
        print(f'  WARNING: CTCF load failed {chrom}:{start}-{end}: {e}')
        vals = np.zeros(end - start)

    L      = len(vals) - (len(vals) % bin_size)
    binned = vals[:L].reshape(-1, bin_size).mean(axis=1)

    half_comparable_nt = comp_half
    centers_bp = (np.arange(len(binned)) + 0.5) * bin_size - half_comparable_nt

    print(f'  bin_ctcf: {chrom}:{start}-{end} ({end-start} nt) → {len(binned)} bins')
    return binned, centers_bp

SAL_THRESHOLD = 10
def apply_saliency_threshold(ref_out: np.ndarray,
                              alt_out: np.ndarray,
                              threshold: float = SAL_THRESHOLD
                              ) -> tuple[np.ndarray, np.ndarray]:
    """
    Zero out bins where max(REF, ALT) < threshold.
    Operates elementwise — keeps array lengths unchanged.
    """
    mask = np.maximum(ref_out, alt_out) >= threshold
    return ref_out * mask, alt_out * mask

# def plot_variant(orig_idx, chrom, insert_pos, REF_sal, insertions_data, ctcf_vals, centers_bp):
#     n_ins  = len(insertions_data)
#     #n_rows = 1 + n_ins * 2 + 1
#     n_rows = 1 + n_ins + 1

#     fig, axes = plt.subplots(n_rows, 1, figsize=(20, 2.5 * n_rows), sharex=True)
#     fig.suptitle(f'{orig_idx}  {chrom}:{insert_pos}', fontsize=16)

#     sal_axes   = []
#     delta_axes = []
#     ax_idx = 0

#     ######## REF
#     axes[ax_idx].plot(centers_bp, REF_sal, color=COL_REF)
#     axes[ax_idx].set_ylabel('REF saliency')
#     axes[ax_idx].axvline(0, color='black', linestyle='--', alpha=0.5)
#     sal_axes.append(axes[ax_idx])
#     ax_idx += 1

#     ######## ALT + delta rows 
#     for ins in insertions_data:
#         name       = ins['alu_name']
#         alt_binned = ins['alt_binned']
#         delta      = alt_binned - REF_sal

#         axes[ax_idx].plot(centers_bp, alt_binned, color=COL_ALT)
#         axes[ax_idx].set_ylabel(f'{name}\nsaliency')
#         axes[ax_idx].axvline(0, color='black', linestyle='--', alpha=0.5)
#         sal_axes.append(axes[ax_idx])
#         ax_idx += 1

# #         axes[ax_idx].plot(centers_bp, delta, color=COL_DELTA)
# #         axes[ax_idx].axhline(0, color='black', linewidth=0.5, linestyle='--')
# #         axes[ax_idx].axvline(0, color='black', linestyle='--', alpha=0.5)
# #         axes[ax_idx].set_ylabel(f'{name}\ndelta')
# #         delta_axes.append(axes[ax_idx])
# #         ax_idx += 1

#     ######## CTCF 
#     axes[ax_idx].plot(centers_bp, ctcf_vals, color=COL_CTCF)
#     axes[ax_idx].set_ylabel('CTCF')
#     axes[ax_idx].set_xlabel('Distance from insertion site (bp)')
#     axes[ax_idx].axvline(0, color='black', linestyle='--', alpha=0.5)

#     ########link y-axes within each group 
#     for ax in sal_axes[1:]:
#         sal_axes[0].get_shared_y_axes().join(sal_axes[0], ax)
#         ax.autoscale()

# #     for ax in delta_axes[1:]:
# #         delta_axes[0].get_shared_y_axes().join(delta_axes[0], ax)
# #         ax.autoscale()

#     plt.tight_layout(rect=[0, 0.03, 1, 0.97])
#     return fig

def plot_var(orig_idx, alu_names=None):
    with h5py.File(f'{H5_dir}variant_{orig_idx}.h5', 'r') as var_hdf, \
        pdf_backend.PdfPages(f'{fig_dir}variant_{orig_idx}.pdf') as pdf:

        row=scores_query[scores_query['orig_idx']==orig_idx]
        chrom=row['CHROM'].values[0]
        alu_start=row['POS'].values[0]
        alu_end=row['END'].values[0]
        insert_pos= (alu_start + alu_end) // 2

        var_keys=var_hdf.keys()

        REF_grads  = var_hdf['REF_saliency_nt'][()]   # raw [seq_len, 4] — binned inside align_and_crop

        ins_keys = [k for k in var_keys if k.startswith('insertion_')]
        if alu_names:
            ins_keys = [f'insertion_{n}' for n in alu_names if f'insertion_{n}' in var_keys]

        insertions_data = []
        centers_bp      = None
        ref_out_final   = None                         

        for ins_key in ins_keys:
            alu_name      = ins_key.replace('insertion_', '')
            ins_grp       = var_hdf[ins_key]
            align_score   = float(ins_grp.attrs['align_score'])
            consensus_len = int(ins_grp.attrs['consensus_len'])
            ALT_grads     = ins_grp['ALT_saliency_nt'][()]

            ref_out, alt_out, centers_bp, alt_insert = align_and_crop(REF_grads, ALT_grads, consensus_len)
            ref_out, alt_out = apply_saliency_threshold(ref_out, alt_out)   # ← add this

            if ref_out_final is None:
                ref_out_final = ref_out  

            insertions_data.append(dict(
                alu_name    = alu_name,
                alt_binned  = alt_out,
                alt_insert  = alt_insert,
                align_score = align_score,
            ))

        # remove the separate REF_binned line entirely, use ref_out_final instead:
        # REF_binned = nt_to_saliency_binned(REF_grads)   ← delete this

        ctcf_binned, ctcf_centers = bin_ctcf(chrom, insert_pos, consensus_len)

#         assert np.allclose(centers_bp, ctcf_centers), \
#             f'saliency vs CTCF center mismatch!\n{centers_bp[:3]}…\n{ctcf_centers[:3]}…'

        fig = plot_variant(orig_idx, chrom, insert_pos,
                           ref_out_final,           # ← not REF_binned
                           insertions_data,
                           ctcf_binned, centers_bp)
        pdf.savefig(fig)
        plt.close(fig)

def splice_insertion(comparable: np.ndarray,
                     centers_bp: np.ndarray,
                     insert_payload: np.ndarray | None,
                     consensus_len_bp: int,
                     bin_size: int = 512
                     ) -> tuple[np.ndarray, np.ndarray]:
    """
    Split `comparable` at the center (where the insertion was excised) and
    splice `insert_payload` into the gap (NaN if payload is None).
 
    Alus are typically shorter than `bin_size` (≈300 bp vs 512 bp), so the
    insertion is rendered as at least one bin spanning [-ins_half, +ins_half]
    using the mean of whatever payload is provided. This keeps the x-axis
    geometry honest while still showing saliency *at* the Alu.
    """
    split = len(comparable) // 2
    left, right = comparable[:split], comparable[split:]
 
    n_ins_bins = max(1, consensus_len_bp // bin_size)
 
    if insert_payload is None or len(insert_payload) == 0:
        payload = np.full(n_ins_bins, np.nan)
    else:
        payload = np.asarray(insert_payload, dtype=float)
        if len(payload) == n_ins_bins:
            pass
        elif n_ins_bins == 1:
            # insertion shorter than a bin: collapse to one value
            payload = np.array([np.nanmean(payload)])
        elif len(payload) > n_ins_bins:
            payload = payload[:n_ins_bins]
        else:
            payload = np.concatenate(
                [payload, np.full(n_ins_bins - len(payload), np.nan)])
 
    full = np.concatenate([left, payload, right])
 
    # rebuild x-axis: left half keeps its negative coords, insertion spans
    # [-ins_half, +ins_half] bp, right half keeps its positive coords.
    ins_half = consensus_len_bp / 2.0
    left_centers  = centers_bp[:split] - ins_half
    right_centers = centers_bp[split:] + ins_half
    # insertion bin centers spread evenly across [-ins_half, +ins_half]
    ins_centers = np.linspace(-ins_half, ins_half, n_ins_bins, endpoint=False) \
                  + (consensus_len_bp / n_ins_bins) / 2.0
    full_centers = np.concatenate([left_centers, ins_centers, right_centers])
 
    return full, full_centers
 
def insertion_mean_saliency(ALT_grads: np.ndarray,
                            consensus_len_bp: int,
                            akita_crop_nt: int = 32 * 2048
                            ) -> float:
    """
    Mean per-nt saliency over the inserted region (computed at nt resolution,
    not binned). Used to populate the single insertion bin when the Alu is
    shorter than `bin_size`.
    """
    sal = np.sum(np.abs(ALT_grads.astype(np.float32)), axis=1)
    cropped = sal[akita_crop_nt : len(sal) - akita_crop_nt]
    cropped_len = len(cropped)
    ins_half = consensus_len_bp // 2
    mid = cropped_len // 2
    ins = cropped[mid - ins_half : mid - ins_half + consensus_len_bp]
    return float(np.mean(ins)) if len(ins) else float('nan')

def filter_peaks_to_window(peaks_df,
                           chrom: str,
                           insert_pos: int,
                           window_half_bp: int):
    """
    Return peaks overlapping [insert_pos - window_half, insert_pos + window_half]
    on `chrom`, with start/end converted to bp relative to insert_pos.
    """
    win_start = insert_pos - window_half_bp
    win_end   = insert_pos + window_half_bp
    sub = peaks_df[(peaks_df['chrom'] == chrom)
                   & (peaks_df['end']   >= win_start)
                   & (peaks_df['start'] <= win_end)].copy()
    sub['rel_start'] = sub['start'] - insert_pos
    sub['rel_end']   = sub['end']   - insert_pos
    # clip to window so rectangles don't overshoot the axis
    sub['rel_start'] = sub['rel_start'].clip(lower=-window_half_bp)
    sub['rel_end']   = sub['rel_end'].clip(upper= window_half_bp)
    return sub

# ─────────────────────────────────────────────────────────────────────────────
def plot_variant(orig_idx, chrom, insert_pos,
                 ref_full, insertions_data, ctcf_full, centers_full,
                 consensus_len_bp, peaks_df=None):
    n_ins  = len(insertions_data)
    has_peaks = peaks_df is not None and len(peaks_df) > 0
    n_rows = 1 + n_ins + 1 + (1 if has_peaks else 0)
 
    height_ratios = [1.0] * (1 + n_ins + 1)
    if has_peaks:
        height_ratios.append(0.35)
 
    fig_h = ROW_H_IN * sum(height_ratios) + 0.7
    fig, axes = plt.subplots(n_rows, 1,
                             figsize=(FIG_W_IN, fig_h),
                             sharex=True,
                             constrained_layout=True,
                             gridspec_kw={'height_ratios': height_ratios})
    if n_rows == 1:
        axes = [axes]
 
    fig.suptitle(f'{orig_idx}   {chrom}:{insert_pos:,}',
                 fontsize=8, y=1.00)
 
    ins_half = consensus_len_bp / 2.0
 
    def _shade_insertion(ax):
        # Use the larger of (consensus_len, min_visible_bp) so the highlight
        # is visible even at full-Mb x-axis scales. Alus are ~300 bp; on a
        # ±500 kb axis that's sub-pixel, so widen for visibility but keep
        # the true insertion edges marked with thin dashed lines.
        x_range = np.nanmax(centers_full) - np.nanmin(centers_full)
        min_visible = max(consensus_len_bp, x_range * 0.005)   # ≥0.5% of axis
        ax.axvspan(-min_visible/2, min_visible/2,
                   facecolor='#FFD27F', alpha=0.45, linewidth=0, zorder=0)
        # true insertion boundaries
        ax.axvline(-ins_half, color='#B26A00', linestyle='-',
                   linewidth=0.4, alpha=0.7)
        ax.axvline( ins_half, color='#B26A00', linestyle='-',
                   linewidth=0.4, alpha=0.7)
        ax.axvline(0, color='#555555', linestyle='--',
                   linewidth=0.5, alpha=0.5)
 
    sal_data = []   # keep refs to (ax, ydata) for clean y-sharing later
 
    # REF row
    ax = axes[0]
    _shade_insertion(ax)
    ax.plot(centers_full, ref_full, color=COL_REF)
    ax.set_ylabel('REF', rotation=0, ha='right', va='center', labelpad=12)
    sal_data.append((ax, ref_full))
 
    # ALT rows
    for i, ins in enumerate(insertions_data, start=1):
        ax = axes[i]
        _shade_insertion(ax)
        ax.plot(centers_full, ins['alt_full'], color=COL_ALT)
        label = f"{ins['alu_name']}\nalign={ins['align_score']:.2f}"
        ax.set_ylabel(label, rotation=0, ha='right', va='center', labelpad=12)
        sal_data.append((ax, ins['alt_full']))
 
    # CTCF continuous row
    ax_ctcf = axes[1 + n_ins]
    _shade_insertion(ax_ctcf)
    ax_ctcf.plot(centers_full, ctcf_full, color=COL_CTCF)
    ax_ctcf.set_ylabel('CTCF', rotation=0, ha='right', va='center', labelpad=12)
 
  # CTCF peaks row
    ax_bottom = ax_ctcf
    if has_peaks:
        ax_pk = axes[-1]
        _shade_insertion(ax_pk)
        # widen peaks that would render sub-pixel at this zoom
        x_range = np.nanmax(centers_full) - np.nanmin(centers_full)
        min_peak_w = x_range * 0.003   # ≥0.3% of axis width
 
        classes_present = []
        for _, pk in peaks_df.iterrows():
            cls = pk['peak_class']
            if cls not in classes_present:
                classes_present.append(cls)
            color = PEAK_CLASS_COLORS.get(cls, '#999999')
            start, end = pk['rel_start'], pk['rel_end']
            if end - start < min_peak_w:
                mid = (start + end) / 2
                start, end = mid - min_peak_w/2, mid + min_peak_w/2
            ax_pk.axvspan(start, end,
                          ymin=0.05, ymax=0.95,
                          facecolor=color, edgecolor='none',
                          alpha=0.9, zorder=2)
        ax_pk.set_ylim(0, 1)
        ax_pk.set_yticks([])
        ax_pk.set_ylabel('peaks', rotation=0, ha='right', va='center', labelpad=12)
        ax_pk.spines['left'].set_visible(False)
 
        handles = [mpl.patches.Patch(facecolor=PEAK_CLASS_COLORS.get(c, '#999999'),
                                     edgecolor='none', label=c)
                   for c in classes_present]
        # legend above the peaks strip (between CTCF and peaks rows) so it
        # never collides with the x-axis label
        ax_pk.legend(handles=handles, loc='lower center',
                     bbox_to_anchor=(0.5, 1.05),
                     ncol=min(len(handles), 3),
                     frameon=False, handlelength=1.2, handleheight=0.8,
                     columnspacing=1.2, borderaxespad=0.0)
        ax_bottom = ax_pk
        
    # share y across saliency tracks so REF/ALT are comparable
    if len(sal_data) > 1:
        ymax = max(np.nanmax(y) for _, y in sal_data)
        ymin = min(np.nanmin(y) for _, y in sal_data)
        pad  = (ymax - ymin) * 0.05 if ymax > ymin else 1.0
        for ax, _ in sal_data:
            ax.set_ylim(ymin - pad*0.1, ymax + pad)
 
    # tidy x ticks — use kb labels for readability
    def _kb(x, _):
        return f'{x/1000:g}'
    ax_bottom.xaxis.set_major_formatter(mpl.ticker.FuncFormatter(_kb))
    ax_bottom.set_xlabel('Distance from insertion site (kb)')
 
    return fig

# ─────────────────────────────────────────────────────────────────────────────
def plot_variant(orig_idx, chrom, insert_pos,
                 ref_full, insertions_data, ctcf_full, centers_full,
                 consensus_len_bp, peaks_df=None):
    n_ins  = len(insertions_data)
    has_peaks = peaks_df is not None and len(peaks_df) > 0
    n_rows = 1 + n_ins + 1 + (1 if has_peaks else 0)
 
    height_ratios = [1.0] * (1 + n_ins + 1)
    if has_peaks:
        height_ratios.append(0.35)
 
    fig_h = ROW_H_IN * sum(height_ratios) + 0.7
    fig, axes = plt.subplots(n_rows, 1,
                             figsize=(FIG_W_IN, fig_h),
                             sharex=True,
                             constrained_layout=True,
                             gridspec_kw={'height_ratios': height_ratios})
    if n_rows == 1:
        axes = [axes]
 
    fig.suptitle(f'{orig_idx}   {chrom}:{insert_pos:,}',
                 fontsize=8, y=1.00)
 
    ins_half = consensus_len_bp / 2.0
 
    def _shade_insertion(ax):
        # Use the larger of (consensus_len, min_visible_bp) so the highlight
        # is visible even at full-Mb x-axis scales. Alus are ~300 bp; on a
        # ±500 kb axis that's sub-pixel, so widen for visibility but keep
        # the true insertion edges marked with thin dashed lines.
        x_range = np.nanmax(centers_full) - np.nanmin(centers_full)
        min_visible = max(consensus_len_bp, x_range * 0.005)   # ≥0.5% of axis
        ax.axvspan(-min_visible/2, min_visible/2,
                   facecolor='#FFD27F', alpha=0.45, linewidth=0, zorder=0)
        # true insertion boundaries
        ax.axvline(-ins_half, color='#B26A00', linestyle='-',
                   linewidth=0.4, alpha=0.7)
        ax.axvline( ins_half, color='#B26A00', linestyle='-',
                   linewidth=0.4, alpha=0.7)
        ax.axvline(0, color='#555555', linestyle='--',
                   linewidth=0.5, alpha=0.5)
 
    sal_data = []   # keep refs to (ax, ydata) for clean y-sharing later
 
    # REF row
    ax = axes[0]
    _shade_insertion(ax)
    ax.plot(centers_full, ref_full, color=COL_REF)
    ax.set_ylabel('REF', rotation=0, ha='right', va='center', labelpad=12)
    sal_data.append((ax, ref_full))
 
    # ALT rows
    for i, ins in enumerate(insertions_data, start=1):
        ax = axes[i]
        _shade_insertion(ax)
        ax.plot(centers_full, ins['alt_full'], color=COL_ALT)
        label = f"{ins['alu_name']}\nalign={ins['align_score']:.2f}"
        ax.set_ylabel(label, rotation=0, ha='right', va='center', labelpad=12)
        sal_data.append((ax, ins['alt_full']))
 
    # CTCF continuous row
    ax_ctcf = axes[1 + n_ins]
    _shade_insertion(ax_ctcf)
    ax_ctcf.plot(centers_full, ctcf_full, color=COL_CTCF)
    ax_ctcf.set_ylabel('CTCF', rotation=0, ha='right', va='center', labelpad=12)
 
  # CTCF peaks row
    ax_bottom = ax_ctcf
    if has_peaks:
        ax_pk = axes[-1]
        _shade_insertion(ax_pk)
        # widen peaks that would render sub-pixel at this zoom
        x_range = np.nanmax(centers_full) - np.nanmin(centers_full)
        min_peak_w = x_range * 0.003   # ≥0.3% of axis width
 
        classes_present = []
        for _, pk in peaks_df.iterrows():
            cls = pk['peak_class']
            if cls not in classes_present:
                classes_present.append(cls)
            color = PEAK_CLASS_COLORS.get(cls, '#999999')
            start, end = pk['rel_start'], pk['rel_end']
            if end - start < min_peak_w:
                mid = (start + end) / 2
                start, end = mid - min_peak_w/2, mid + min_peak_w/2
            ax_pk.axvspan(start, end,
                          ymin=0.05, ymax=0.95,
                          facecolor=color, edgecolor='none',
                          alpha=0.9, zorder=2)
        ax_pk.set_ylim(0, 1)
        ax_pk.set_yticks([])
        ax_pk.set_ylabel('peaks', rotation=0, ha='right', va='center', labelpad=12)
        ax_pk.spines['left'].set_visible(False)
 
        handles = [mpl.patches.Patch(facecolor=PEAK_CLASS_COLORS.get(c, '#999999'),
                                     edgecolor='none', label=c)
                   for c in classes_present]
        # legend above the peaks strip (between CTCF and peaks rows) so it
        # never collides with the x-axis label
        ax_pk.legend(handles=handles, loc='lower center',
                     bbox_to_anchor=(0.5, 1.05),
                     ncol=min(len(handles), 3),
                     frameon=False, handlelength=1.2, handleheight=0.8,
                     columnspacing=1.2, borderaxespad=0.0)
        ax_bottom = ax_pk
        
    # share y across saliency tracks so REF/ALT are comparable
    if len(sal_data) > 1:
        ymax = max(np.nanmax(y) for _, y in sal_data)
        ymin = min(np.nanmin(y) for _, y in sal_data)
        pad  = (ymax - ymin) * 0.05 if ymax > ymin else 1.0
        for ax, _ in sal_data:
            ax.set_ylim(ymin - pad*0.1, ymax + pad)
 
    # tidy x ticks — use kb labels for readability
    def _kb(x, _):
        return f'{x/1000:g}'
    ax_bottom.xaxis.set_major_formatter(mpl.ticker.FuncFormatter(_kb))
    ax_bottom.set_xlabel('Distance from insertion site (kb)')
 
    return fig

def plot_var(orig_idx, alu_names=None, peaks_df=None):
    with plt.rc_context(PLOT_RC), \
         h5py.File(f'{H5_dir}variant_{orig_idx}.h5', 'r') as var_hdf, \
         pdf_backend.PdfPages(f'{fig_dir}variant_{orig_idx}.pdf') as pdf:
 
        row       = scores_query[scores_query['orig_idx'] == orig_idx]
        chrom     = row['CHROM'].values[0]
        alu_start = row['POS'].values[0]
        alu_end   = row['END'].values[0]
        insert_pos = (alu_start + alu_end) // 2
 
        REF_grads = var_hdf['REF_saliency_nt'][()]
        var_keys  = list(var_hdf.keys())
        ins_keys  = [k for k in var_keys if k.startswith('insertion_')]
        if alu_names:
            ins_keys = [f'insertion_{n}' for n in alu_names
                        if f'insertion_{n}' in var_keys]
 
        insertions_data    = []
        centers_full       = None
        ref_full_final     = None
        consensus_len_used = None
 
        for ins_key in ins_keys:
            alu_name      = ins_key.replace('insertion_', '')
            ins_grp       = var_hdf[ins_key]
            align_score   = float(ins_grp.attrs['align_score'])
            consensus_len = int(ins_grp.attrs['consensus_len'])
            ALT_grads     = ins_grp['ALT_saliency_nt'][()]
 
            ref_out, alt_out, centers_bp, alt_insert = align_and_crop(
                REF_grads, ALT_grads, consensus_len)
            ref_out, alt_out = apply_saliency_threshold(ref_out, alt_out)
 
            # If the Alu is shorter than bin_size, alt_insert will be empty.
            # Compute a single mean-saliency value at nt resolution instead.
            if len(alt_insert) == 0:
                alt_insert = np.array(
                    [insertion_mean_saliency(ALT_grads, consensus_len)])
 
            # splice insertion back in so we can see saliency at the Alu
            alt_full, c_full = splice_insertion(
                alt_out, centers_bp, alt_insert, consensus_len)
            ref_full, _      = splice_insertion(
                ref_out, centers_bp, None, consensus_len)   # NaN over insertion
 
            if ref_full_final is None:
                ref_full_final     = ref_full
                centers_full       = c_full
                consensus_len_used = consensus_len
 
            insertions_data.append(dict(
                alu_name    = alu_name,
                alt_full    = alt_full,
                align_score = align_score,
            ))
 
        # CTCF over the same comparable window, then pad to match
        ctcf_binned, ctcf_centers = bin_ctcf(chrom, insert_pos, consensus_len_used)
        ctcf_full, _ = splice_insertion(
            ctcf_binned, ctcf_centers, None, consensus_len_used)
        if len(ctcf_full) != len(centers_full):
            print(f'  WARN: ctcf len={len(ctcf_full)} vs centers={len(centers_full)}')
 
        # filter peaks to visible window
        peaks_in_window = None
        if peaks_df is not None:
            window_half = int(np.nanmax(np.abs(centers_full)))
            peaks_in_window = filter_peaks_to_window(
                peaks_df, chrom, insert_pos, window_half)
            
        fig = plot_variant(orig_idx, chrom, insert_pos,
                           ref_full_final, insertions_data,
                           ctcf_full, centers_full,
                           consensus_len_used,
                           peaks_df=peaks_in_window)
        pdf.savefig(fig)
        plt.close(fig)
 
    print(f'Saved → {fig_dir}variant_{orig_idx}.pdf')

# ─── Pair-deletion contact maps  (analogous to plot_alu_consensus_inserts) ────
# Run after the plot_var cell; re-imports are safe.

sys.path.append(f'{SUPREMO_DIR}/')
sys.path.insert(0, f'{SUPREMO_DIR}/scripts/')

_hg38 = pysam.FastaFile(f'{SUPREMO_DIR}/data/hg38.fa')
_MB   = 1_048_576
_BIN  = _akita.bin_size               # 2048
_TLC  = _akita.target_length_cropped  # 448

_cmap_pred = mcolors.LinearSegmentedColormap.from_list(
    'rdbu_dp',
    ['#1a1f6b','#3a62b0','#7a9fd4','#c4cfe8','#f7f7f7',
     '#f2b49a','#d4765f','#a8303f','#5c1020'], N=256)
_cmap_diff = mcolors.LinearSegmentedColormap.from_list(
    'prgn_dp',
    ['#3b2255','#6b3d7a','#9b6aaa','#c4a3c8','#e8dff0','#f7f7f7',
     '#ccddb8','#8aaf72','#4e8a4e','#2a5c35','#0f3320'], N=256)
_cmap_inter = mcolors.LinearSegmentedColormap.from_list(
    'puor_inter',
    ['#2d004b','#542788','#8073ac','#b2abd2','#d8daeb','#f7f7f7',
     '#fee0b6','#fdb863','#e08214','#b35806','#7f3b08'], N=256)

def _revcomp(seq):
    t = str.maketrans('ACGTacgtNn', 'TGCAtgcaNn')
    return seq.translate(t)[::-1]

def _fetch_del(chrom, center, deletions, fa=_hg38, mb=_MB):
    """1 MB sequence centred on `center` with `deletions` excised (mirrors run_del_pairs.py)."""
    half = mb // 2
    cb = next(((ds, de) for ds, de in deletions if ds <= center < de), None)
    cc = (cb[0] + cb[1]) // 2 if cb else center

    i0, i1 = cc - half, cc - half + mb
    lc = rc_ = 0
    for ds, de in deletions:
        os_, oe = max(ds, i0), min(de, i1)
        if os_ < oe:
            d = oe - os_
            if   oe <= cc:  lc  += d
            elif os_ >= cc: rc_ += d
            else:           lc  += cc - os_; rc_ += oe - cc

    fs = int(max(0, i0 - lc))
    fe = int(i1 + rc_)
    seq = list(fa.fetch(chrom, fs, fe))

    dbc = 0
    for ds, de in deletions:
        rs, re = max(0, ds - fs), min(len(seq), de - fs)
        if rs < re:
            for i in range(rs, re): seq[i] = ''
            dbc += max(0, min(de, cc) - max(ds, fs))

    s  = ''.join(seq)
    ci = (cc - fs) - dbc
    a, b = ci - half, ci - half + mb
    lp = rp = 0
    if a < 0:       lp = -a;          a = 0; b = mb
    if b > len(s):  rp = b - len(s);  b = len(s)
    out = s[a:b]
    if lp: out = 'N' * lp + out
    if rp: out = out + 'N' * rp
    return (out + 'N' * (mb - len(out)))[:mb]

def _additivity_scores(ref_m, alt1_m, alt2_m, alt12_m):
    """
    Additivity decomposition for one orientation.
    interaction = delta12 - (delta1 + delta2);  >0 synergistic, <0 buffering.
    Returns dict with mse_nonadd, mean_inter, mean_delta_alu1/ctcf2/both.
    """
    d1   = alt1_m  - ref_m
    d2   = alt2_m  - ref_m
    d12  = alt12_m - ref_m
    inter = d12 - (d1 + d2)
    valid = ~(np.isnan(ref_m) | np.isnan(alt1_m) | np.isnan(alt2_m) | np.isnan(alt12_m)).reshape(-1)
    iv   = inter.reshape(-1)[valid]
    return {
        'mse_nonadd':       float(np.mean(iv ** 2)),
        'mean_inter':       float(np.mean(iv)),
        'mean_delta_alu1':  float(np.mean(d1.reshape(-1)[valid])),
        'mean_delta_ctcf2': float(np.mean(d2.reshape(-1)[valid])),
        'mean_delta_both':  float(np.mean(d12.reshape(-1)[valid])),
    }

def get_ctcf_at_saliency_peaks(orig_idx, peaks_df, sal_threshold=SAL_THRESHOLD,
                                saliency_source='REF', window_bp=200_000):
    """
    Return CTCF peaks whose bin has saliency >= sal_threshold, using the same
    pipeline as plot_var (align_and_crop → apply_saliency_threshold) so the
    result matches exactly what is plotted.

    Parameters
    ----------
    saliency_source : 'REF' — use ref_out (REF saliency over comparable region)
                      str   — Alu name (e.g. 'AluYk12'); use alt_out
    """
    row        = scores_query[scores_query['orig_idx'] == orig_idx].iloc[0]
    chrom      = row['CHROM']
    insert_pos = (int(row['POS']) + int(row['END'])) // 2

    with h5py.File(f'{H5_dir}variant_{orig_idx}.h5', 'r') as hf:
        REF_grads = hf['REF_saliency_nt'][()]

        if saliency_source == 'REF':
            # need any insertion for consensus_len / ALT_grads to run align_and_crop
            ins_keys = [k for k in hf.keys() if k.startswith('insertion_')]
            if not ins_keys:
                raise ValueError('No insertions in h5; cannot run align_and_crop')
            ins_grp = hf[ins_keys[0]]
        else:
            key = f'insertion_{saliency_source}'
            if key not in hf:
                raise KeyError(f'{key!r} not in h5. Available: {list(hf.keys())}'
                               f'  (pass saliency_source="REF" or an Alu name)')
            ins_grp = hf[key]

        consensus_len = int(ins_grp.attrs['consensus_len'])
        ALT_grads     = ins_grp['ALT_saliency_nt'][()]

    ref_out, alt_out, centers_bp, _ = align_and_crop(REF_grads, ALT_grads, consensus_len)
    ref_out, alt_out = apply_saliency_threshold(ref_out, alt_out)
    sal = ref_out if saliency_source == 'REF' else alt_out

    nearby = peaks_df[
        (peaks_df['chrom'] == chrom) &
        (peaks_df['end']   >= insert_pos - window_bp) &
        (peaks_df['start'] <= insert_pos + window_bp)
    ].copy()

    if len(nearby) == 0:
        return nearby

    def _sal(r):
        mid_rel = (r['start'] + r['end']) / 2 - insert_pos
        idx = int(np.argmin(np.abs(centers_bp - mid_rel)))
        # reject if nearest bin centre is more than one bin away (outside plotted range)
        if np.abs(centers_bp[idx] - mid_rel) > BIN_SIZE:
            return 0.0
        return float(sal[idx])

    nearby['sal'] = nearby.apply(_sal, axis=1)
    return nearby[nearby['sal'] >= sal_threshold].drop(columns='sal').copy()

def build_del_pair_maps(orig_idx, ctcf2_start, ctcf2_end, ctcf2_chrom=None):
    """
    Build REF, Alu1_DEL, CTCF2_DEL, Alu1+CTCF2_DEL contact maps centred on the
    primary Alu.  Predictions run in both fwd and RC; fwd maps are returned for
    plotting; interaction / MSE / Spearman scores for both orientations go in meta.

    Parameters
    ----------
    orig_idx    : primary Alu's orig_idx in scores_query
    ctcf2_start : genomic start of the CTCF site to delete
    ctcf2_end   : genomic end of the CTCF site to delete
    ctcf2_chrom : chromosome (must match Alu; inferred if None)

    Returns
    -------
    REF_fwd, Alu1_DEL_fwd, CTCF2_DEL_fwd, Alu1_CTCF2_DEL_fwd : 448×448 ndarray
    meta : dict — coords, fwd+rc MSE/Spearman/additivity scores, map geometry
    """
    row   = scores_query[scores_query['orig_idx'] == orig_idx].iloc[0]
    chrom = row['CHROM']
    a1s   = int(row['POS'])
    a1e   = int(row['END'])
    ctr   = (a1s + a1e) // 2
    c2s, c2e = int(ctcf2_start), int(ctcf2_end)

    if ctcf2_chrom is not None and ctcf2_chrom != chrom:
        raise ValueError(f'CTCF2 chrom {ctcf2_chrom!r} != Alu1 chrom {chrom!r}')

    seqs = [
        _fetch_del(chrom, ctr, []),
        _fetch_del(chrom, ctr, [(a1s, a1e)]),
        _fetch_del(chrom, ctr, [(c2s, c2e)]),
        _fetch_del(chrom, ctr, [(a1s, a1e), (c2s, c2e)]),
    ]

    def _both(seq):
        fv = _akita.vector_from_seq(seq)
        rv = _akita.vector_from_seq(_revcomp(seq))
        return fv, rv, _akita.map_from_vector(fv), _akita.map_from_vector(rv)

    (ref_fv,  ref_rv,  REF_fwd,  REF_rc),  \
    (a1_fv,   a1_rv,   A1_fwd,   A1_rc),   \
    (c2_fv,   c2_rv,   C2_fwd,   C2_rc),   \
    (a1c2_fv, a1c2_rv, A1C2_fwd, A1C2_rc) = [_both(s) for s in seqs]

    add_fwd = _additivity_scores(REF_fwd, A1_fwd, C2_fwd, A1C2_fwd)
    add_rc  = _additivity_scores(REF_rc,  A1_rc,  C2_rc,  A1C2_rc)

    win_start = ctr - _MB // 2
    def _bin(bp): return _akita.get_bin(bp - win_start)

    meta = dict(
        orig_idx        = orig_idx,
        chrom           = chrom,
        alu1_start      = a1s,    alu1_end  = a1e,
        alu1_name       = row.get('repName', 'Alu1'),
        ctcf2_start     = c2s,    ctcf2_end = c2e,
        ctcf2_id        = f'{chrom}_{c2s}_{c2e}',
        alu1_lines      = [_bin(a1s), _bin(a1e)],
        ctcf2_lines     = [_bin(c2s), _bin(c2e)],
        rel_pos_map     = _bin(a1s),
        map_start_coord = win_start + 32 * _BIN,
        alu1_len_bins   = math.ceil((a1e - a1s) / _BIN),
        # ── fwd metrics ──
        mse_alu1_fwd    = _su.mse(ref_fv, a1_fv),
        mse_ctcf2_fwd   = _su.mse(ref_fv, c2_fv),
        mse_both_fwd    = _su.mse(ref_fv, a1c2_fv),
        corr_alu1_fwd   = _su.spearman(ref_fv, a1_fv),
        corr_ctcf2_fwd  = _su.spearman(ref_fv, c2_fv),
        corr_both_fwd   = _su.spearman(ref_fv, a1c2_fv),
        **{f'{k}_fwd': v for k, v in add_fwd.items()},
        # ── rc metrics ──
        mse_alu1_rc     = _su.mse(ref_rv, a1_rv),
        mse_ctcf2_rc    = _su.mse(ref_rv, c2_rv),
        mse_both_rc     = _su.mse(ref_rv, a1c2_rv),
        corr_alu1_rc    = _su.spearman(ref_rv, a1_rv),
        corr_ctcf2_rc   = _su.spearman(ref_rv, c2_rv),
        corr_both_rc    = _su.spearman(ref_rv, a1c2_rv),
        **{f'{k}_rc': v for k, v in add_rc.items()},
    )
    return REF_fwd, A1_fwd, C2_fwd, A1C2_fwd, meta

def plot_del_pair(REF_map, Alu1_DEL_map, CTCF2_DEL_map, Alu1_CTCF2_DEL_map, meta,
                  cmap_pred=_cmap_pred, cmap_diff=_cmap_diff, cmap_inter=_cmap_inter,
                  scale=1, gene_rows=3, show_insulation=True):
    """
    Triangle-heatmap plot for an Alu × CTCF pair deletion, matching the
    layout of plot_consensus_inserts.

    Maps shown are fwd-strand predictions.  MSE / Spearman panel labels show fwd
    values; the interaction panel title reports both fwd and rc additivity scores.

    Panels (top → bottom):
      REF | Alu1 DEL | CTCF2 DEL | Alu1+CTCF2 DEL
      Alu1_DEL−REF | CTCF2_DEL−REF | Alu1+CTCF2_DEL−REF
      Interaction (non-additive)
      [ Insulation track ]  [ Gene track ]
    """
    m     = meta
    lines = m['alu1_lines'] + m['ctcf2_lines']

    interact = ( (Alu1_CTCF2_DEL_map - REF_map)
               - (Alu1_DEL_map       - REF_map)
               - (CTCF2_DEL_map      - REF_map) )

    panels = [
        (REF_map,
         f"REF  ({m['alu1_name']} + CTCF2 present)",
         cmap_pred, -2, 2),
        (Alu1_DEL_map,
         f"Alu1 DEL   MSE={m['mse_alu1_fwd']:.4f}  ρ={m['corr_alu1_fwd']:.3f}",
         cmap_pred, -2, 2),
        (CTCF2_DEL_map,
         f"CTCF2 DEL ({m['ctcf2_id']})   MSE={m['mse_ctcf2_fwd']:.4f}  ρ={m['corr_ctcf2_fwd']:.3f}",
         cmap_pred, -2, 2),
        (Alu1_CTCF2_DEL_map,
         f"Alu1+CTCF2 DEL   MSE={m['mse_both_fwd']:.4f}  ρ={m['corr_both_fwd']:.3f}",
         cmap_pred, -2, 2),
        (Alu1_DEL_map - REF_map,
         "Alu1_DEL − REF",
         cmap_diff, -1, 1),
        (CTCF2_DEL_map - REF_map,
         "CTCF2_DEL − REF",
         cmap_diff, -1, 1),
        (Alu1_CTCF2_DEL_map - REF_map,
         "Alu1+CTCF2_DEL − REF",
         cmap_diff, -1, 1),
        (interact,
         f"Interaction (non-additive)"
         f"   mse_nonadd: fwd={m['mse_nonadd_fwd']:.4f}  rc={m['mse_nonadd_rc']:.4f}"
         f"   mean_inter: fwd={m['mean_inter_fwd']:.3f}  rc={m['mean_inter_rc']:.3f}",
         cmap_inter, -0.5, 0.5),
    ]

    # gene track
    try:
        genes_in_map = _pu.get_genes_in_map(
            m['chrom'], m['map_start_coord'],
            rel_pos_map=m['rel_pos_map'],
            SVTYPE='DEL', SVLEN=m['alu1_end'] - m['alu1_start'])
        has_genes = isinstance(genes_in_map, pd.DataFrame) and len(genes_in_map) > 0
    except Exception:
        has_genes = False

    if show_insulation:
        ins_ref   = _sc.insulation_track(REF_map)
        ins_alu1  = _sc.insulation_track(Alu1_DEL_map)
        ins_ctcf2 = _sc.insulation_track(CTCF2_DEL_map)
        ins_both  = _sc.insulation_track(Alu1_CTCF2_DEL_map)

    lw  = 0.5 * scale
    pw  = 2   * scale
    pw1 = 0.8 * scale
    gh  = pw1 * 2 * gene_rows / 3

    row_heights = []
    for _ in panels:
        row_heights += [pw, .1]
    row_heights = row_heights[:-1]
    if show_insulation:
        row_heights += [.1, pw1]
    if has_genes:
        row_heights += [.15, gh]

    fig, gs = gridspec_inches([pw * 2], row_heights)

    for ri, (mat, title, cm, vn, vx) in enumerate(panels):
        ax = plt.subplot(gs[ri * 2, 0])
        _pu.pcolormesh_45deg(plt, mat, lines, lw, cmap=cm, vmax=vx, vmin=vn)
        ax.annotate(title,
                    xy=(0, 0.9), xycoords='axes fraction',
                    xytext=(0, 4 * scale), textcoords='offset points',
                    fontsize=10 * scale, ha='left', va='bottom',
                    annotation_clip=False)

    if show_insulation:
        ins_row = len(panels) * 2
        ax_ins  = plt.subplot(gs[ins_row, 0])
        x = np.arange(_TLC)
        ax_ins.plot(x, ins_ref,   color='#3a62b0', lw=lw, label='REF')
        ax_ins.plot(x, ins_alu1,  color='#888888', lw=lw, label='Alu1 DEL',       ls='--')
        ax_ins.plot(x, ins_ctcf2, color='#a8303f', lw=lw, label='CTCF2 DEL',      ls='-.')
        ax_ins.plot(x, ins_both,  color='#2a5c35', lw=lw, label='Alu1+CTCF2 DEL', ls=':')
        ax_ins.set_xlim([0, _TLC])
        ax_ins.yaxis.tick_right()
        ax_ins.legend(fontsize=9 * scale, loc='upper right', ncol=2)
        ax_ins.set_ylabel('Insulation', fontsize=7 * scale, rotation=90)
        for ln in lines:
            ax_ins.axvline(x=ln, color='gray', linestyle='dashed', lw=lw)
        ax_ins.set_xticks([]); ax_ins.set_yticks([])

    if has_genes:
        gene_row = len(panels) * 2 + (2 if show_insulation else 0)
        ax_g = fig.add_subplot(gs[gene_row, 0])
        gene_track = list(genes_in_map[['Start', 'width']].to_records(index=False))
        gbh = 15 * gene_rows / 3
        bar_locs = []
        for lvl in range(gene_rows):
            loc = (lvl + 1) * (gbh / gene_rows) - 0.5
            ax_g.broken_barh(gene_track[lvl:][::gene_rows], (loc, 0.5 * scale),
                             facecolors='tab:blue')
            bar_locs.append(loc)
        bar_locs = (bar_locs * math.ceil(len(genes_in_map) / gene_rows))[:len(genes_in_map)]
        gr = genes_in_map.reset_index(drop=True)
        for i in range(len(gr)):
            ax_g.annotate(gr.loc[i, 'Gene'],
                          (gr.loc[i, 'Start'], bar_locs[i]),
                          rotation=45, ha='right', va='top',
                          annotation_clip=False, fontsize=10 * scale)
        for ln in lines:
            ax_g.axvline(x=ln, color='gray', linestyle='dashed', lw=lw)
        ax_g.set_ylabel('Genes', rotation=90)
        ax_g.set_ylim([0, gbh])
        ax_g.set_xlim([0, _TLC])
        ax_g.axis('off')

    fig.suptitle(
        f"idx {m['orig_idx']}  {m['chrom']}:{m['alu1_start']}-{m['alu1_end']}"
        f"  |  {m['alu1_name']} × CTCF2 {m['ctcf2_id']}"
        f"  |  mean_inter  fwd={m['mean_inter_fwd']:.3f}  rc={m['mean_inter_rc']:.3f}",
        fontsize=10 * scale, y=1.02)

    plt.show()
    return fig

#filepaths
H5_dir    = f'{PROJECT_DIR}/bin/ins_gradients/20260519_topMSE/h5_out/'
CTCF_BW    = f'{DATA_DIR}/alus/data/encode/bw/ENCFF209TQB_CTCF.bigWig'
fig_dir    = f'{PROJECT_DIR}/bin/ins_gradients/20260519_topMSE/pdf_out/'
BIN_SIZE   = 2048
SMOOTH_WIN = 5       # rolling mean window (in bins)

#also read in set of universal CTCFs and such 
CTCF_shared=pd.read_csv(f'{DATA_DIR}/alus/data/encode/CTCF/CTCF_shared_peaks.txt', sep='\t', index_col=0)
CTCF_shared_short=CTCF_shared[['chrom', 'start', 'end', 'n_cell_lines', 'is_shared', 'is_unique', 'peak_class']]
CTCF_shared_short

scores_query=pd.read_csv(f'{ins_grad_dir}20260519_topMSE/top_ins_noPDCH.txt', sep='\t', index_col=0)
scores_query

scores_query[scores_query['orig_idx']==330025]

indices= {
    330025: ['AluSc5', 'AluSg7', 'AluYk12'],
    308144: ['AluJo', 'AluYa5', 'AluYb8', 'AluSz6', 'AluSx1'],
    307931: ['AluSq2', 'AluSp', 'AluSg7', 'AluYa5', 'AluYb8'], 
    808325: ['AluSc', 'AluSg', 'AluYa8', 'AluYk12']
}

for idx in indices:
    alu_plot=indices[idx]
    plot_var(orig_idx=idx, alu_names=alu_plot, peaks_df=CTCF_shared_short)

scores_query

#pair alu ctcf del 
#308144
#330025

query_idx=308144

ctcf_candidates = get_ctcf_at_saliency_peaks(
    orig_idx      = query_idx,
    peaks_df      = CTCF_shared_short,
    sal_threshold = 4000,
    saliency_source    = 'REF',
    window_bp     = 600_000,
)
ctcf_candidates[['chrom','start','end','peak_class']]
ctcf_candidates

fig_dir_ctcf=f'{PROJECT_DIR}/bin/ins_gradients/20260519_topMSE/ctcf_pdf_out/'

with plt.rc_context(PLOT_RC), \
     pdf_backend.PdfPages(f'{fig_dir_ctcf}variant_{query_idx}_ctcf_puor.pdf') as pdf:

    for i in range(len(ctcf_candidates)):

        ctcf_row = ctcf_candidates.iloc[i]

        REF, A1, C2, A1C2, meta = build_del_pair_maps(
            orig_idx    = query_idx,
            ctcf2_start = ctcf_row['start'],
            ctcf2_end   = ctcf_row['end'],
        )
        print(f"mean_inter  fwd={meta['mean_inter_fwd']:.3f}  rc={meta['mean_inter_rc']:.3f}")
        print(f"mse_nonadd  fwd={meta['mse_nonadd_fwd']:.4f}  rc={meta['mse_nonadd_rc']:.4f}")
        fig = plot_del_pair(REF, A1, C2, A1C2, meta, scale=1.2, show_insulation=True, gene_rows=3)

        pdf.savefig(fig)
        plt.close(fig)

SAL_THRESHOLD

# ### Mutate based on alignment differences

res_dir=f'{PROJECT_DIR}/results/20250717_HIV/'
alu_insert_scores=pd.read_csv(f'{res_dir}20250724_HIVtopinserts_consensusAluinserts_scores.txt',
                              sep='\t', index_col=0)
pd.set_option('display.max_columns', None)
alu_insert_scores

#insert_pos_list=['58321293', '24233452', '53588197', '52780910', '58130752', '58028190']
insert_pos_list=['24233452']

#look at 24233452 -- AluYc vs AluJo(more similar to REF)

#look at alignment differences between Alus -- individually test the nucleotides that are different between them 
#between Alu seqs, do global alignment 

for pos in insert_pos_list:
    print('running insert position', pos)
    alu_insert_row=alu_insert_scores[alu_insert_scores['insert_position']==float(pos)]
    consensus_alu_names=ast.literal_eval(alu_insert_row['mse_subfam_names'].iloc[0])
    print(consensus_alu_names)

    for i in (combinations(consensus_alu_names,2)):
        print(i)
        alu_seq_a=alu_dfam_filtered[alu_dfam_filtered['AluName']==i[0]]['Consensus_seq'].iloc[0]
        alu_seq_b=alu_dfam_filtered[alu_dfam_filtered['AluName']==i[1]]['Consensus_seq'].iloc[0]

        alignments=pairwise2.align.globalxx(alu_seq_a, alu_seq_b, one_alignment_only=True)
        aln1, aln2, score, start, end = alignments[0]
        print(f'score {i}', score)
        
        print('alu1: ', aln1[50:100])
        print('alu2: ', aln2[50:100])
        
        diffs = []
        pos1 = pos2 = 0  # 1-based positions in original sequences

        for b1, b2 in zip(aln1, aln2):
            if b1 != '-':
                pos1 += 1
            if b2 != '-':
                pos2 += 1

            if b1 == b2:
                continue

            if b1 == '-':
                diffs.append({'type': 'INS', 'seq1_pos': pos1, 'seq2_pos': pos2, 'base_seq1': '-', 'base_seq2': b2})
            elif b2 == '-':
                diffs.append({'type': 'DEL', 'seq1_pos': pos1, 'seq2_pos': pos2, 'base_seq1': b1, 'base_seq2': '-'})
            else:
                diffs.append({'type': 'SUB', 'seq1_pos': pos1, 'seq2_pos': pos2, 'base_seq1': b1, 'base_seq2': b2})

        align_df=pd.DataFrame(diffs)
        print(align_df)
        break

#look at alignment differences between Alus -- individually test the nucleotides that are different between them 
#between Alu seqs, do global alignment 
#format so the first one is more like the REF 

def alu_pair_align(alu_tup):
    print('running insert position', alu_tup[0])
    alu_insert_row=alu_insert_scores[alu_insert_scores['insert_position']==float(alu_tup[0])]
    consensus_alu_names=ast.literal_eval(alu_insert_row['mse_subfam_names'].iloc[0])
    #print(consensus_alu_names)

    alu_seq_a=alu_dfam_filtered[alu_dfam_filtered['AluName']==alu_tup[1]]['Consensus_seq'].iloc[0]
    alu_seq_b=alu_dfam_filtered[alu_dfam_filtered['AluName']==alu_tup[2]]['Consensus_seq'].iloc[0]

    alignments=pairwise2.align.globalxx(alu_seq_a, alu_seq_b, one_alignment_only=True)
    aln1, aln2, score, start, end = alignments[0]
    print(f'score {alu_tup[1]}, {alu_tup[2]}', score)

    print('alu1: ', aln1[50:100])
    print('alu2: ', aln2[50:100])

    diffs = []
    pos1 = pos2 = -1  # 1-based positions in original sequences

    for b1, b2 in zip(aln1, aln2):
        if b1 != '-':
            pos1 += 1
        if b2 != '-':
            pos2 += 1

        if b1 == b2:
            continue

        if b1 == '-':
            diffs.append(['INS', pos1, pos2,  '-', b2])
        elif b2 == '-':
            diffs.append(['DEL', pos1,  pos2,  b1,  '-'])
        else:
            diffs.append(['SUB', pos1, pos2, b1, b2])

    align_df=pd.DataFrame(diffs)
    align_df.columns=['type', 'seq1_pos', 'seq2_pos', 'base_seq1', 'base_seq2']
                                
    return(align_df, aln1, aln2)

def apply_edits_single(seq, edits_df):
    seq_list = [b for b in seq if b != '-']
    shift = 0  # track net shift due to insertions/deletions
    
    # sort by position ascending
    for pos, group in edits_df.groupby("seq1_pos", sort=True):
        idx = pos + shift  # convert to 0-based with shift
        
        # separate deletions and insertions
        dels = group[group["type"] == "DEL"]
        ins = group[group["type"] == "INS"]
        
        # handle deletions: remove one base per DEL
        if len(dels) > 0:
            seq_list[idx : idx + len(dels)] = []  # remove the bases
            shift -= len(dels)
        
        # handle insertions: insert concatenated bases at current index
        if len(ins) > 0:
            ins_bases = "".join(ins["base_seq2"].tolist())
            seq_list[idx:idx] = list(ins_bases)
            shift += len(ins_bases)

    edited_seq = "".join(seq_list)
    #print('50 to 100', edited_seq[50:100])
    final_seq=edited_seq.replace("-", "")
    print(final_seq[50:100])
    
    return final_seq

edited_seq_list=[]

alu_pair_test=[('24233452', 'AluJo', 'AluYc')]
align_df, aln1, aln2=alu_pair_align(alu_pair_test[0])

print(align_df[0:15])

for idx in np.unique(align_df['seq1_pos']):
    print(idx)
    edited_seq=apply_edits_single(aln1, align_df[align_df['seq1_pos']==idx])
                     
    edited_seq_list.append((idx, edited_seq))

#look get predictions for each individual bp altered sequence 
chrom='chr19'
insert_pos=24233452
REF_seq = hg19_fa.fetch(chrom, int(np.floor(insert_pos-MB/2)), int(np.floor(insert_pos+MB/2)))
REF_pred = utils.vector_from_seq(REF_seq)
REF_pred_mat = utils.mat_from_vector(REF_pred)

Alu_name_list=['AluJo', 'AluYc']
Alu_map_list=[]

for Alu in Alu_name_list:

    consensus_seq_row = alu_dfam_filtered[alu_dfam_filtered['AluName']==Alu]
    insert_seq = consensus_seq_row['Consensus_seq'].iloc[0]

    len_alu = len(insert_seq)
    region_pad = (MB - len_alu)/2
    left_start = np.floor(insert_pos - region_pad)
    right_end = np.floor(insert_pos + region_pad)
    left_seq_INS = hg19_fa.fetch(chrom, int(left_start), int(insert_pos))
    right_seq_INS = hg19_fa.fetch(chrom, int(insert_pos), int(right_end))
    ALT_seq = left_seq_INS + insert_seq + right_seq_INS
    
    ALT_pred = utils.vector_from_seq(ALT_seq)
    ALT_pred_mat = utils.mat_from_vector(ALT_pred)
    
    Alu_map_list.append(ALT_pred_mat)

alt_map_list=[]
alt_idx_list=[]

for tup in edited_seq_list:
    
    print(f'running insert at {tup[0]}')

    insert_seq=tup[1]
    len_alu = len(insert_seq)
    region_pad = (MB - len_alu)/2
    left_start = np.floor(insert_pos - region_pad)
    right_end = np.floor(insert_pos + region_pad)
    left_seq_INS = hg19_fa.fetch(chrom, int(left_start), int(insert_pos))
    right_seq_INS = hg19_fa.fetch(chrom, int(insert_pos), int(right_end))
    ALT_seq = left_seq_INS + insert_seq + right_seq_INS
    
    ALT_pred = utils.vector_from_seq(ALT_seq)
    ALT_pred_mat = utils.mat_from_vector(ALT_pred)
    
    alt_map_list.append(ALT_pred_mat)
    alt_idx_list.append(tup[0])

plot_maps_edits(REF_pred_mat, Alu_map_list, Alu_name_list, alt_map_list, (chrom, insert_pos), alt_idx_list)

edited_seq_list=[]

alu_pair_test=[('58130752', 'AluSz6', 'AluYc')]
align_df, aln1, aln2=alu_pair_align(alu_pair_test[0])

print(align_df[0:5])

for idx in np.unique(align_df['seq1_pos']):
    print(idx)
    edited_seq=apply_edits_single(aln1, align_df[align_df['seq1_pos']==idx])
                     
    edited_seq_list.append((idx, edited_seq))

#look get predictions for each individual bp altered sequence 
chrom='chr19'
insert_pos=58130752
REF_seq = hg19_fa.fetch(chrom, int(np.floor(insert_pos-MB/2)), int(np.floor(insert_pos+MB/2)))
REF_pred = utils.vector_from_seq(REF_seq)
REF_pred_mat = utils.mat_from_vector(REF_pred)

Alu_name_list=['AluSz6', 'AluYc']
Alu_map_list=[]

for Alu in Alu_name_list:

    consensus_seq_row = alu_dfam_filtered[alu_dfam_filtered['AluName']==Alu]
    insert_seq = consensus_seq_row['Consensus_seq'].iloc[0]

    len_alu = len(insert_seq)
    region_pad = (MB - len_alu)/2
    left_start = np.floor(insert_pos - region_pad)
    right_end = np.floor(insert_pos + region_pad)
    left_seq_INS = hg19_fa.fetch(chrom, int(left_start), int(insert_pos))
    right_seq_INS = hg19_fa.fetch(chrom, int(insert_pos), int(right_end))
    ALT_seq = left_seq_INS + insert_seq + right_seq_INS
    
    ALT_pred = utils.vector_from_seq(ALT_seq)
    ALT_pred_mat = utils.mat_from_vector(ALT_pred)
    
    Alu_map_list.append(ALT_pred_mat)

alt_map_list=[]
alt_idx_list=[]

for tup in edited_seq_list:
    
    print(f'running insert at {tup[0]}')

    insert_seq=tup[1]
    len_alu = len(insert_seq)
    region_pad = (MB - len_alu)/2
    left_start = np.floor(insert_pos - region_pad)
    right_end = np.floor(insert_pos + region_pad)
    left_seq_INS = hg19_fa.fetch(chrom, int(left_start), int(insert_pos))
    right_seq_INS = hg19_fa.fetch(chrom, int(insert_pos), int(right_end))
    ALT_seq = left_seq_INS + insert_seq + right_seq_INS
    
    ALT_pred = utils.vector_from_seq(ALT_seq)
    ALT_pred_mat = utils.mat_from_vector(ALT_pred)
    
    alt_map_list.append(ALT_pred_mat)
    alt_idx_list.append(tup[0])

plot_maps_edits(REF_pred_mat, Alu_map_list, Alu_name_list, alt_map_list, (chrom, insert_pos), alt_idx_list)

#plot maps with 1 nuc changes 
def plot_maps_edits(REF_mat, Alu_mats, Alu_name_list, edit_maps, query_seq, names_list, max_cols=4):

    #get number of rows and cols 
    n_maps = len(edit_maps)+len(Alu_mats)+1
    ncols = min(max_cols, n_maps)
    nrows = math.ceil(n_maps / ncols)

    fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(5*ncols, 5*nrows))
    if n_maps == 1:
        axes = [axes]
    else:
        axes = axes.flatten()

    vmin, vmax = -2, 2

    # plot the reference first
    im = axes[0].matshow(REF_mat, cmap='RdBu_r', vmin=vmin, vmax=vmax)
    axes[0].set_title("REF", y=1.1)
    plt.colorbar(im, ax=axes[0], fraction=0.046, pad=0.04, ticks=[-2, -1, 0, 1, 2])
    
    #plot the ALT maps we want to see
    for i, (m, name) in enumerate(zip(Alu_mats, Alu_name_list), start=1):
        im = axes[i].matshow(m, cmap='RdBu_r', vmin=vmin, vmax=vmax)
        axes[i].set_title(f"{name}", y=1.1)    

    # plot shuffled/variant maps
    for i, (m, name) in enumerate(zip(edit_maps, names_list), start=len(Alu_mats)+1):
        im = axes[i].matshow(m, cmap='RdBu_r', vmin=vmin, vmax=vmax)
        #axes[i].set_title(f"{name} Insertion\nMSE = {score:.4f}", y=1.1)
        axes[i].set_title(f"{name}", y=1.1)
        #plt.colorbar(im, ax=axes[i], fraction=0.046, pad=0.04, ticks=[-2, -1, 0, 1, 2])

    # hide unused subplots if any
    for j in range(n_maps, len(axes)):
        axes[j].axis('off')

    fig.suptitle(
        f"Insertion site: {query_seq[0]}:{query_seq[1]}",
        fontsize=16, y=0.92
    )

    plt.tight_layout(rect=[0, 0, 1, 0.98])
    plt.subplots_adjust(top=0.9, wspace=0.2, hspace=0.2)

    #return fig

res_dir=f'{PROJECT_DIR}/results/20250717_HIV/'
alu_insert_scores=pd.read_csv(f'{res_dir}20250724_HIVtopinserts_consensusAluinserts_scores.txt',
                              sep='\t', index_col=0)
pd.set_option('display.max_columns', None)
alu_insert_scores
