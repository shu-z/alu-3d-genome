"""
Whole ~1 Mb window REF vs. ARID3B-motif-deleted saliency for AluY elements,
binned at Akita's native 2048 bp resolution, with HFF CTCF ChIP-seq peaks
(ENCFF294RSZ) plotted underneath.

Companion to 20260810_arid3b_motif_deletion_examples.py, which zooms in on
the Alu itself -- this one shows the same REF/motif-deleted saliency
comparison but across the full model input window, binned down to bin
resolution, to see how the deletion's effect propagates relative to CTCF
sites in the window. Reuses that script's Akita-loading, sequence-excision,
and saliency machinery directly (imported as a module) rather than
duplicating it.
"""

import os
import importlib.util
import sys
import warnings

from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import pandas as pd
import pysam

# Paths - edit for your environment
PROJECT_DIR = os.environ.get("ALU_PROJECT_DIR", "/pollard/home/szhang20/alu")
DATA_DIR = os.environ.get("ALU_DATA_DIR", "/pollard/data/projects/shzhang")
POLLARD_DATA = os.environ.get("POLLARD_DATA", "/pollard/data")

_CORE_PATH = f'{PROJECT_DIR}/bin/20260810_arid3b_motif_deletion_examples.py'
_spec = importlib.util.spec_from_file_location('arid3b_core', _CORE_PATH)
core = importlib.util.module_from_spec(_spec)
sys.modules['arid3b_core'] = core
_spec.loader.exec_module(core)

CTCF_BED   = f'{DATA_DIR}/alus/data/encode/HFF_CTCF_ENCFF294RSZ.bed.gz'
BIN_SIZE   = 2048              # Akita's native bin resolution
N_EXAMPLES = None              # None = every AluY/ARID3B overlap found
OUT_PDF    = f'{core.alu_dir}figs/20260813_arid3b_window_binned_saliency_ctcf_AluY.pdf'

COL_CTCF   = '#4E8B57'
COL_ALUDEL = '#8064A2'   # matches ins_gradients notebook's COL_DELTA
CROP_BINS  = 32   # Akita edge crop: predictions unreliable within crop_bins*bin_size
                  # of the window edges (matches 20260303_alu_insertion_gradients.ipynb)

# figure sizing, matching ins_gradients/20260303_alu_insertion_gradients.ipynb (cell 14)
MM_PER_IN = 25.4
FIG_W_MM  = 120
FIG_W_IN  = FIG_W_MM / MM_PER_IN
ROW_H_IN  = 0.55

def bin_track(track, bin_size):
    """
    Max-pool a per-nucleotide track into contiguous bin_size-bp bins.
    Max (not mean) so saliency spikes survive binning instead of being
    averaged away -- matches the approach in
    ins_gradients/20260303_alu_insertion_gradients.ipynb (align_and_crop's
    _bin helper), which is why that notebook's saliency plots reach y-values
    in the thousands instead of the tens/hundreds a mean-bin gives.
    """
    n_bins  = len(track) // bin_size
    trimmed = track[:n_bins * bin_size].reshape(n_bins, bin_size)
    return np.nanmax(trimmed, axis=1)

def load_ctcf_peaks(path):
    """ENCODE narrowPeak: chrom start end name score strand signalValue pValue qValue summit."""
    cols = ['chrom', 'start', 'end', 'name', 'score', 'strand',
            'signalValue', 'pValue', 'qValue', 'summit']
    return pd.read_csv(path, sep='\t', header=None, names=cols)

def ctcf_peaks_in_window(ctcf_df, chrom, win_start, win_end):
    return ctcf_df[(ctcf_df['chrom'] == chrom) &
                    (ctcf_df['end'] > win_start) & (ctcf_df['start'] < win_end)]

def plot_window_binned_saliency_with_ctcf(hg38_fa, seqnn_model, ctcf_df,
                                           chrom, alu_start, alu_end, alu_strand,
                                           orig_idx, repName, motif_spans,
                                           motif_name='ARID3B', bin_size=BIN_SIZE):
    motif_spans = [(int(s), int(e)) for s, e in motif_spans]
    if not motif_spans:
        print(f'  [skip] {chrom}:{alu_start:,}-{alu_end:,}: no motif spans')
        return None

    n_sites_raw    = len(motif_spans)
    merged_genomic = core._merge_spans(motif_spans)
    n_sites_merged = len(merged_genomic)

    # REF once, then two deletions off the same REF sequence: the ARID3B
    # motif only, and the whole Alu ("AluDEL") -- avoids recomputing REF
    # gradients twice.
    REF_seq, win_start = core.build_ref_seq(hg38_fa, chrom, alu_start, alu_end)
    win_end  = win_start + core.MB
    ref_gpos = np.arange(win_start, win_end, dtype=np.int64)
    ref_sal  = core.saliency_from_grads(core.compute_nucleotide_grads(REF_seq, seqnn_model))

    local_motif_spans = [(s - win_start, e - win_start) for s, e in motif_spans]
    DEL_motif_seq, del_motif_gpos = core._excise_and_repad(
        hg38_fa, chrom, REF_seq, win_start, core._merge_spans(local_motif_spans))
    del_motif_sal = core.saliency_from_grads(
        core.compute_nucleotide_grads(DEL_motif_seq, seqnn_model))

    local_alu_span = [(alu_start - win_start, alu_end - win_start)]
    DEL_alu_seq, del_alu_gpos = core._excise_and_repad(
        hg38_fa, chrom, REF_seq, win_start, local_alu_span)
    del_alu_sal = core.saliency_from_grads(
        core.compute_nucleotide_grads(DEL_alu_seq, seqnn_model))

    # drop the outer CROP_BINS*bin_size on each side: Akita predictions
    # aren't reliable that close to the input edges
    crop_nt    = CROP_BINS * bin_size
    crop_start = win_start + crop_nt
    crop_end   = win_end - crop_nt

    _, rt_full, _ = core.align_tracks(ref_gpos, ref_sal, ref_gpos, ref_sal, crop_start, crop_end)
    _, mt_full, _ = core.align_tracks(del_motif_gpos, del_motif_sal, del_motif_gpos, del_motif_sal,
                                       crop_start, crop_end)
    _, at_full, _ = core.align_tracks(del_alu_gpos, del_alu_sal, del_alu_gpos, del_alu_sal,
                                       crop_start, crop_end)

    rt_binned = bin_track(rt_full, bin_size)
    mt_binned = bin_track(mt_full, bin_size)
    at_binned = bin_track(at_full, bin_size)
    n_bins    = len(rt_binned)

    bin_edges   = crop_start + np.arange(n_bins + 1) * bin_size
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    xr_bins     = bin_centers - alu_start

    ctcf_sites = ctcf_peaks_in_window(ctcf_df, chrom, crop_start, crop_end)

    height_ratios = [1, 1, 1, 0.6]
    fig_h = ROW_H_IN * sum(height_ratios) + 0.7
    fig, axes = plt.subplots(
        4, 1, figsize=(FIG_W_IN, fig_h), constrained_layout=True,
        gridspec_kw={'height_ratios': height_ratios}, sharex=True)

    sal_panels = [
        (axes[0], rt_binned, core.COL_REF, 'REF'),
        (axes[1], mt_binned, core.COL_DEL, f'REF\n[{motif_name} del]'),
        (axes[2], at_binned, COL_ALUDEL, 'REF\n[Alu del]'),
    ]
    for ax, track, col, lab in sal_panels:
        ax.bar(xr_bins, track, width=bin_size, color=col, linewidth=0)
        ax.axvspan(0, alu_end - alu_start, color='0.85', zorder=0)
        for s, e in merged_genomic:
            ax.axvspan(s - alu_start, e - alu_start, color=core.COL_MOTIF,
                       alpha=0.6, zorder=2, linewidth=0)
        ax.set_ylabel(lab, rotation=0, ha='right', va='center')
        ax.margins(x=0)
        for sp in ('top', 'right'):
            ax.spines[sp].set_visible(False)

    # cap/share the y-axis across the three saliency panels so REF, motif-del,
    # and Alu-del are on the same scale (matches plot_variant's sal_data
    # y-sharing in 20260303_alu_insertion_gradients.ipynb)
    all_sal = np.concatenate([rt_binned, mt_binned, at_binned])
    ymin, ymax = np.nanmin(all_sal), np.nanmax(all_sal)
    pad = (ymax - ymin) * 0.05 if ymax > ymin else 1.0
    for ax, _, _, _ in sal_panels:
        ax.set_ylim(ymin - pad * 0.1, ymax + pad)

    ax_ctcf = axes[3]
    for _, pk in ctcf_sites.iterrows():
        mid   = (pk['start'] + pk['end']) / 2 - alu_start
        width = max(pk['end'] - pk['start'], bin_size * 0.3)
        ax_ctcf.bar(mid, pk['signalValue'], width=width, color=COL_CTCF, linewidth=0)
    ax_ctcf.axvspan(0, alu_end - alu_start, color='0.85', zorder=0)
    ax_ctcf.set_ylabel('HFF\nCTCF', rotation=0, ha='right', va='center')
    ax_ctcf.margins(x=0)
    for sp in ('top', 'right'):
        ax_ctcf.spines[sp].set_visible(False)

    if alu_strand == '-':
        for ax in axes:
            ax.invert_xaxis()

    axes[-1].xaxis.set_major_formatter(mticker.FuncFormatter(lambda x, _: f'{x / 1000:g}'))
    axes[-1].set_xlabel(
        f'{chrom}:{alu_start:,}-{alu_end:,} ({alu_strand}), distance from Alu start (kb)')
    fig.suptitle(
        f'{repName}  orig_idx={orig_idx}  {chrom}:{alu_start:,}-{alu_end:,}  '
        f'{motif_name} sites: {n_sites_raw} raw / {n_sites_merged} merged  '
        f'(max-binned {bin_size} bp, edge-cropped {crop_nt:,} bp)', fontsize=8)
    return fig

if __name__ == '__main__':
    hg38_fa = pysam.FastaFile(core.hg38_path)
    params_model, tlen_crop, bin_size_model = core.load_model_params(core.params_file)
    seqnn_model = core.load_akita(params_model, core.repo_path, core.model_file)
    assert bin_size_model == BIN_SIZE, f'Akita bin size {bin_size_model} != {BIN_SIZE}'

    df_high = pd.read_csv(
        f'{core.alu_dir}results/paper_results/20260413_alu_sample/aluhigh_1000samp.txt',
        sep='\t', index_col=0)

    track_cols = ['chr', 'start', 'end', 'tf', 'rel_score', '-log10(pval)', 'strand']
    track_df = pd.read_csv(f'{core.tf_track_dir}{core.TF_track_id}', sep='\t', header=None,
                            skiprows=1, names=track_cols)
    track_df['start'] = track_df['start'].astype(int)
    track_df['end']   = track_df['end'].astype(int)
    track_df_query = track_df[track_df['rel_score'] >= core.MOTIF_MIN_SCORE]
    print(f'{core.TF_name} track loaded: {len(track_df):,} sites, '
          f'{len(track_df_query):,} >= score {core.MOTIF_MIN_SCORE}')

    motif_index = core._build_chrom_overlap_index(track_df_query)
    overlap_mask = [core._has_overlap(motif_index, r.CHROM, r.POS, r.END)
                     for r in df_high.itertuples()]
    df_overlap = df_high[overlap_mask]
    df_overlap_AluY = df_overlap[df_overlap['repName'].str.contains(core.ALU_SUBFAMILY)]
    print(f'{len(df_overlap_AluY)} {core.ALU_SUBFAMILY} elements overlap a strong {core.TF_name} motif')

    ctcf_df = load_ctcf_peaks(CTCF_BED)
    print(f'HFF CTCF peaks loaded: {len(ctcf_df):,}')

    n_written = 0
    with warnings.catch_warnings(), plt.rc_context(core.PLOT_RC), PdfPages(OUT_PDF) as pdf:
        warnings.simplefilter('ignore')

        for _, row in df_overlap_AluY.iterrows():
            if N_EXAMPLES is not None and n_written >= N_EXAMPLES:
                break

            CHROM, POS, END = row['CHROM'], int(row['POS']), int(row['END'])
            strand, orig_idx, repName = row['strand'], row['orig_idx'], row['repName']

            sites = track_df_query[
                (track_df_query['chr']   == CHROM) &
                (track_df_query['start'] >= POS) &
                (track_df_query['end']   <= END)
            ]
            if sites.empty:
                continue
            motif_spans = list(zip(sites['start'].astype(int), sites['end'].astype(int)))

            print(f'Processing {repName} orig_idx={orig_idx}  {CHROM}:{POS}-{END}')
            fig = plot_window_binned_saliency_with_ctcf(
                hg38_fa, seqnn_model, ctcf_df, CHROM, POS, END, strand, orig_idx, repName,
                motif_spans, motif_name=core.TF_name)

            if fig is not None:
                pdf.savefig(fig, bbox_inches='tight')
                plt.close(fig)
                n_written += 1

    print(f'Saved {n_written} examples -> {OUT_PDF}')
