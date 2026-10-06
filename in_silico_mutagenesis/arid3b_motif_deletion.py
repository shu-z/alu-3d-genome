"""
ARID3B motif-deletion examples in AluY elements.

For a handful of AluY x ARID3B-motif overlaps (same variant table / TF track
used in delete_motifs_arid3b.ipynb and delete_motifs_in_regions_fixed.ipynb),
draws one PDF page per example with:

  1. a schematic of the Alu showing exactly where the ARID3B motif(s) sit
     within it (position + width, in bp from the Alu start), and
  2. REF vs. motif-deleted per-nucleotide saliency, both zoomed on the Alu
     and over the full ~1 Mb Akita input window.

No HDF5 -- gradients are computed on the fly, same approach as
delete_motifs_arid3b.ipynb.
"""

import os
import json
import sys
import warnings

from basenji import dna_io
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pysam
import tensorflow as tf

# Paths - edit for your environment
SUPREMO_DIR = os.environ.get("SUPREMO_DIR", "/pollard/home/szhang20/akita_variant_scoring")
AKITA_DIR = os.environ.get("AKITA_DIR", "/pollard/home/szhang20/akita")
PROJECT_DIR = os.environ.get("ALU_PROJECT_DIR", "/pollard/home/szhang20/alu")
POLLARD_DATA = os.environ.get("POLLARD_DATA", "/pollard/data")

# ---- inputs -----------------------------------------------------------
repo_path   = f'{SUPREMO_DIR}/'
model_file  = f'{repo_path}Akita_model/model_best.h5'
params_file = f'{repo_path}Akita_model/params.json'
hg38_path   = f'{repo_path}data/hg38.fa'

alu_dir      = f'{PROJECT_DIR}/'
tf_track_dir = f'{POLLARD_DATA}/tf_binding/tfbs_motifs/JASPAR/UCSC_tracks/2024/hg38/'

TF_name       = 'ARID3B'
TF_track_id   = 'MA0601.2.tsv.gz'
MOTIF_MIN_SCORE = 850          # keep only strong ARID3B sites, as in the notebooks
ALU_SUBFAMILY = 'AluY'         # only AluYs, per request
N_EXAMPLES    = None           # None = plot every AluY/ARID3B overlap found
OUT_PDF       = f'{alu_dir}figs/20260810_arid3b_motif_deletion_examples_AluY.pdf'

MB    = 1_048_576              # Akita input length
SHIFT = 0                      # window shift in bp
FLANK = 150                    # bp of flanking sequence to show around the Alu

COL_REF, COL_DEL, COL_MOTIF = '#3B6FB6', '#C0504D', '#F2C14E'

matplotlib.rcParams.update({
    'pdf.fonttype': 42, 'ps.fonttype': 42, 'svg.fonttype': 'none',
    'font.size': 7, 'savefig.dpi': 600, 'figure.dpi': 150,
})
PLOT_RC = {
    'font.family'        : 'sans-serif',
    'font.sans-serif'    : ['Helvetica', 'Arial', 'DejaVu Sans'],
    'font.size'          : 7,
    'axes.titlesize'     : 8,
    'axes.labelsize'     : 7,
    'xtick.labelsize'    : 6,
    'ytick.labelsize'    : 6,
    'legend.fontsize'    : 6,
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
}

# ---- model / gradients -------------------------------------------------

def load_model_params(params_file):
    with open(params_file) as f:
        params = json.load(f)
    pm = params['model']
    pm['augment_shift'] = 0
    cropping  = pm['head_hic'][5]['cropping']
    tlen      = pm['target_length']
    tlen_crop = tlen - 2 * cropping
    bin_size  = pm['seq_length'] // tlen
    print('bin_size', bin_size)
    return pm, tlen_crop, bin_size

def load_akita(params_model, repo_path, model_file):
    sys.path.append(repo_path)
    from basenji import seqnn_gpu as seqnn
    print('using seqnn_gpu!')
    model = seqnn.SeqNN(params_model)
    model.restore(model_file)
    return model

def compute_nucleotide_grads(seq, model):
    """Per-nucleotide gradient of summed model output w.r.t. one-hot input. [seq_len, 4]"""
    seq_1hot_np = dna_io.dna_1hot(seq).astype(np.float32)
    seq_1hot_np = np.expand_dims(seq_1hot_np, 0)

    seq_input = tf.keras.Input(shape=(model.seq_length, 4), dtype=tf.float32)
    pred_tensor = model.model(seq_input)
    temp_model = tf.keras.Model(inputs=seq_input, outputs=pred_tensor)

    seq_var = tf.Variable(seq_1hot_np, dtype=tf.float32)
    with tf.GradientTape() as tape:
        tape.watch(seq_var)
        pred = temp_model(seq_var)
        loss = tf.reduce_sum(pred)
    grads = tape.gradient(loss, seq_var)
    return grads[0].numpy()

def saliency_from_grads(grads):
    return np.sum(np.abs(grads.astype(np.float32)), axis=1)

# ---- motif-deletion sequence construction ------------------------------

def _build_chrom_overlap_index(motif_df):
    """
    Precompute, per chromosome, motif starts sorted ascending plus a running
    max of ends, so 'does [pos, end) overlap any motif site' is a couple of
    binary searches instead of a full bedtools intersect (avoids depending
    on the bedtools CLI binary, which isn't on PATH in every conda env here).
    """
    index = {}
    for chrom, g in motif_df.groupby('chr'):
        starts = g['start'].to_numpy()
        ends   = g['end'].to_numpy()
        order  = np.argsort(starts)
        starts = starts[order]
        cummax_ends = np.maximum.accumulate(ends[order])
        index[chrom] = (starts, cummax_ends)
    return index

def _has_overlap(chrom_index, chrom, pos, end):
    """True if half-open interval [pos, end) overlaps any indexed motif site."""
    if chrom not in chrom_index:
        return False
    starts, cummax_ends = chrom_index[chrom]
    i = np.searchsorted(starts, end, side='left')
    return bool(i > 0 and cummax_ends[i - 1] > pos)

def _merge_spans(spans):
    """Merge overlapping or adjacent [start, end) spans."""
    spans = sorted((int(s), int(e)) for s, e in spans if int(e) > int(s))
    merged = []
    for s, e in spans:
        if merged and s <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], e))
        else:
            merged.append((s, e))
    return merged

def build_ref_seq(hg38_fa, chrom, alu_start, alu_end):
    """Fetch an MB window centered on the Alu midpoint."""
    alu_mid   = (alu_start + alu_end) // 2
    win_start = alu_mid - MB // 2 + SHIFT
    REF_seq   = hg38_fa.fetch(chrom, win_start, win_start + MB).upper()
    if len(REF_seq) != MB:
        raise ValueError(f'REF length {len(REF_seq)} != MB {MB}')
    return REF_seq, win_start

def _excise_and_repad(hg38_fa, chrom, REF_seq, win_start, merged_local):
    """
    Excise merged seq-local motif spans from REF_seq and repad from genome
    flanks to restore length MB. Returns (DEL_seq, del_gpos), where
    del_gpos[i] is the genomic coordinate of DEL_seq[i] (motif bases are
    simply absent, so they show up as a gap once mapped back to coordinates).
    """
    L       = len(REF_seq)
    win_end = win_start + L
    ref_gpos = np.arange(win_start, win_end, dtype=np.int64)

    in_motif = np.zeros(L, dtype=bool)
    for s, e in merged_local:
        in_motif[s:e] = True

    non_motif_seq  = ''.join(c for c, m in zip(REF_seq, in_motif) if not m)
    non_motif_gpos = ref_gpos[~in_motif]

    total_del = int(in_motif.sum())
    mid       = L // 2
    del_left  = sum(min(e, mid) - s for s, e in merged_local if s < mid)
    del_right = total_del - del_left

    left_flank  = hg38_fa.fetch(chrom, win_start - del_left, win_start).upper()
    right_flank = hg38_fa.fetch(chrom, win_end, win_end + del_right).upper()

    DEL_seq  = left_flank + non_motif_seq + right_flank
    del_gpos = np.concatenate([
        np.arange(win_start - del_left, win_start,           dtype=np.int64),
        non_motif_gpos,
        np.arange(win_end,              win_end + del_right, dtype=np.int64),
    ])

    if len(DEL_seq) != MB:
        raise ValueError(f'DEL length {len(DEL_seq)} != MB {MB}')
    return DEL_seq, del_gpos

def ref_vs_motif_deleted_saliency(hg38_fa, seqnn_model, chrom, alu_start, alu_end,
                                   motif_spans_genomic):
    """Per-nucleotide saliency for REF and motif-deleted sequences."""
    REF_seq, win_start = build_ref_seq(hg38_fa, chrom, alu_start, alu_end)
    ref_gpos = np.arange(win_start, win_start + MB, dtype=np.int64)

    local_spans  = [(s - win_start, e - win_start) for s, e in motif_spans_genomic]
    merged_local = _merge_spans(local_spans)

    DEL_seq, del_gpos = _excise_and_repad(hg38_fa, chrom, REF_seq, win_start, merged_local)

    ref_sal = saliency_from_grads(compute_nucleotide_grads(REF_seq, seqnn_model))
    del_sal = saliency_from_grads(compute_nucleotide_grads(DEL_seq, seqnn_model))
    return ref_gpos, ref_sal, del_gpos, del_sal

def align_tracks(ref_gpos, ref_sal, del_gpos, del_sal, region_start, region_end):
    x         = np.arange(region_start, region_end, dtype=np.int64)
    ref_track = np.full(x.shape, np.nan, dtype=np.float32)
    del_track = np.full(x.shape, np.nan, dtype=np.float32)
    rm = (ref_gpos >= region_start) & (ref_gpos < region_end)
    ref_track[ref_gpos[rm] - region_start] = ref_sal[rm]
    dm = (del_gpos >= region_start) & (del_gpos < region_end)
    del_track[del_gpos[dm] - region_start] = del_sal[dm]
    return x, ref_track, del_track

def _smooth(a, w):
    if w <= 1:
        return a
    k   = np.ones(w) / w
    num = np.convolve(np.nan_to_num(a, nan=0.0), k, mode='same')
    den = np.convolve((~np.isnan(a)).astype(float), k, mode='same')
    out = num / np.where(den == 0, np.nan, den)
    out[np.isnan(a)] = np.nan
    return out

# ---- plotting ------------------------------------------------------------

def draw_motif_location_schematic(ax, alu_start, alu_end, alu_strand,
                                   motif_spans_merged, motif_name, xlim):
    """
    Draws where the TF motif is actually being deleted within the Alu:
    a bar for the Alu body plus one highlighted block per (merged) motif
    span, positioned and sized in bp relative to the Alu start (x=0).
    """
    alu_len = alu_end - alu_start
    ax.add_patch(plt.Rectangle((0, 0.25), alu_len, 0.5,
                                facecolor='0.88', edgecolor='black',
                                linewidth=0.6, zorder=2))

    for s, e in motif_spans_merged:
        s_local = max(s, alu_start) - alu_start
        e_local = min(e, alu_end) - alu_start
        width   = e_local - s_local
        ax.add_patch(plt.Rectangle((s_local, 0.12), width, 0.76,
                                    facecolor=COL_MOTIF, edgecolor='#8a6d1f',
                                    linewidth=0.6, zorder=3))
        mid = (s_local + e_local) / 2
        ax.annotate(f'{motif_name}\n{width} bp @ +{s_local}',
                    (mid, 1.0), ha='center', va='bottom', fontsize=5, zorder=4)

    arrow = '5′ → 3′' if alu_strand == '+' else '3′ ← 5′'
    ax.text(alu_len / 2, 0.5, arrow, ha='center', va='center', fontsize=5,
            color='0.3', zorder=1)

    ax.set_xlim(*xlim)
    ax.set_ylim(0, 1.35)
    ax.set_yticks([])
    ax.set_ylabel('Alu +\nmotif', rotation=0, ha='right', va='center')
    ax.margins(x=0)
    for sp in ('top', 'right', 'left'):
        ax.spines[sp].set_visible(False)

def _draw_saliency_panels(axes, panels, xr, alu_start, alu_end, motif_spans_merged, has_flank):
    for ax, (track, col, lab) in zip(axes, panels):
        ax.fill_between(xr, 0, track, step='mid', color=col, linewidth=0)
        if has_flank:
            p0 = ax.axvspan(0, alu_end - alu_start, color='0.92', zorder=0)
            p0.set_clip_on(False)
        for s, e in motif_spans_merged:
            p = ax.axvspan(s - alu_start, e - alu_start, color=COL_MOTIF,
                            alpha=0.45, zorder=1, linewidth=0)
            p.set_clip_on(False)
        ax.set_ylabel(lab, rotation=0, ha='right', va='center')
        ax.margins(x=0)
        ax.set_ylim(0, 200)
        for sp in ('top', 'right'):
            ax.spines[sp].set_visible(False)

def plot_arid3b_deletion_example(hg38_fa, seqnn_model, chrom, alu_start, alu_end,
                                  alu_strand, orig_idx, repName, motif_spans,
                                  motif_name='ARID3B', flank=FLANK, smooth=3):
    """
    One figure per Alu: motif-location schematic on top, then REF vs.
    motif-deleted saliency zoomed on the Alu itself (+/- flank bp of
    context), no full ~1 Mb window.
    """
    motif_spans = [(int(s), int(e)) for s, e in motif_spans]
    if not motif_spans:
        print(f'  [skip] {chrom}:{alu_start:,}-{alu_end:,}: no motif spans')
        return None

    n_sites_raw    = len(motif_spans)
    merged_genomic = _merge_spans(motif_spans)
    n_sites_merged = len(merged_genomic)

    ref_gpos, ref_sal, del_gpos, del_sal = ref_vs_motif_deleted_saliency(
        hg38_fa, seqnn_model, chrom, alu_start, alu_end, motif_spans)

    x_alu, rt_alu, dt_alu = align_tracks(
        ref_gpos, ref_sal, del_gpos, del_sal, alu_start - flank, alu_end + flank)

    rt_alu, dt_alu = _smooth(rt_alu, smooth), _smooth(dt_alu, smooth)

    xr_alu  = x_alu - alu_start
    has_flank_alu = (x_alu.min() < alu_start) or (x_alu.max() > alu_end - 1)

    panels_alu = [(rt_alu, COL_REF, 'REF'), (dt_alu, COL_DEL, f'REF\n[{motif_name} del]')]

    n_rows = 1 + len(panels_alu)
    fig, axes = plt.subplots(
        n_rows, 1, figsize=(240 / 25.4, 0.7 + 0.95 * (n_rows - 1)),
        constrained_layout=True,
        gridspec_kw={'height_ratios': [0.5] + [1] * (n_rows - 1)},
    )
    axes = np.atleast_1d(axes)

    draw_motif_location_schematic(
        axes[0], alu_start, alu_end, alu_strand, merged_genomic, motif_name,
        xlim=(xr_alu.min(), xr_alu.max()))

    _draw_saliency_panels(axes[1:], panels_alu, xr_alu,
                           alu_start, alu_end, merged_genomic, has_flank_alu)

    if alu_strand == '-':
        for ax in axes:
            ax.invert_xaxis()

    axes[-1].set_xlabel(f'{chrom}:{alu_start:,}-{alu_end:,} ({alu_strand})')
    fig.suptitle(
        f'{repName}  orig_idx={orig_idx}  {chrom}:{alu_start:,}-{alu_end:,}  '
        f'{motif_name} sites: {n_sites_raw} raw / {n_sites_merged} merged',
        fontsize=8)
    return fig

# ---- main ------------------------------------------------------------

if __name__ == '__main__':
    hg38_fa = pysam.FastaFile(hg38_path)
    params_model, tlen_crop, bin_size = load_model_params(params_file)
    seqnn_model = load_akita(params_model, repo_path, model_file)

    # Alu variant table + ARID3B JASPAR track, same source as
    # delete_motifs_arid3b.ipynb / delete_motifs_in_regions_fixed.ipynb
    df_high = pd.read_csv(
        f'{alu_dir}results/paper_results/20260413_alu_sample/aluhigh_1000samp.txt',
        sep='\t', index_col=0)

    track_cols = ['chr', 'start', 'end', 'tf', 'rel_score', '-log10(pval)', 'strand']
    track_df = pd.read_csv(f'{tf_track_dir}{TF_track_id}', sep='\t', header=None,
                            skiprows=1, names=track_cols)
    track_df['start'] = track_df['start'].astype(int)
    track_df['end']   = track_df['end'].astype(int)
    track_df_query = track_df[track_df['rel_score'] >= MOTIF_MIN_SCORE]
    print(f'{TF_name} track loaded: {len(track_df):,} sites, '
          f'{len(track_df_query):,} >= score {MOTIF_MIN_SCORE}')

    motif_index = _build_chrom_overlap_index(track_df_query)
    overlap_mask = [
        _has_overlap(motif_index, row.CHROM, row.POS, row.END)
        for row in df_high.itertuples()
    ]
    df_overlap = df_high[overlap_mask]

    df_overlap_AluY = df_overlap[df_overlap['repName'].str.contains(ALU_SUBFAMILY)]
    print(f'{len(df_overlap_AluY)} {ALU_SUBFAMILY} elements overlap a strong {TF_name} motif')

    n_written = 0
    with warnings.catch_warnings(), plt.rc_context(PLOT_RC), PdfPages(OUT_PDF) as pdf:
        warnings.simplefilter('ignore')

        for _, row in df_overlap_AluY.iterrows():
            if N_EXAMPLES is not None and n_written >= N_EXAMPLES:
                break

            CHROM    = row['CHROM']
            POS      = int(row['POS'])
            END      = int(row['END'])
            strand   = row['strand']
            orig_idx = row['orig_idx']
            repName  = row['repName']

            sites = track_df_query[
                (track_df_query['chr']   == CHROM) &
                (track_df_query['start'] >= POS) &
                (track_df_query['end']   <= END)
            ]
            if sites.empty:
                continue
            motif_spans = list(zip(sites['start'].astype(int), sites['end'].astype(int)))

            print(f'Processing {repName} orig_idx={orig_idx}  {CHROM}:{POS}-{END}')
            fig = plot_arid3b_deletion_example(
                hg38_fa, seqnn_model, CHROM, POS, END, strand, orig_idx, repName,
                motif_spans, motif_name=TF_name)

            if fig is not None:
                pdf.savefig(fig, bbox_inches='tight')
                plt.close(fig)
                n_written += 1

    print(f'Saved {n_written} examples -> {OUT_PDF}')
