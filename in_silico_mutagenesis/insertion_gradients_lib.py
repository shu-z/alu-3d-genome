"""
core.py  —  data loading, GPU setup, Akita model, and scoring functions.

All functions are pure (no global state) and accept explicit arguments,
making them easy to test individually or reuse in notebooks.
"""

import json
import math
import subprocess
import sys

from Bio import pairwise2
from Bio.Seq import Seq
import numpy as np
import pandas as pd
import pysam

from config import Config

#### prep sequences 

def rc_seq(seq: str, cfg: Config) -> str:

    if cfg.reverse_complement:
        seq = str(Seq(seq).reverse_complement())
    return seq

def _fetch_window(hg38_fa, chrom: str, center: int, half: int, shift: int) -> str:
    start = int(np.floor(center - half + shift))
    end   = int(np.floor(center + half + shift))
    return hg38_fa.fetch(chrom, start, end).upper()

def get_ref_alt_sequences(hg38_fa, chrom: str, alu_start: int, alu_end: int,
                          cfg: Config,
                          mode: str = 'INS',
                          consensus_seq: str = None) -> tuple[str, str]:
    """
    Generate REF and ALT sequences matching SuPreMo's get_sequences logic.

    mode='REF' : REF sequence only (ALT == REF, no variant)
    mode='INS' : Insert consensus_seq at alu_start within the REF window (requires consensus_seq)
    mode='DEL' : Delete the Alu (alu_start..alu_end) from the REF window

    Returns (REF_seq, ALT_seq)
    """
    if mode == 'INS' and consensus_seq is None:
        raise ValueError("consensus_seq is required for mode='INS'")

    half    = cfg.MB // 2
    alu_mid = int(np.floor((alu_start + alu_end) / 2))

    # REF window centered on alu_mid + shift
    win_start = int(np.floor(alu_mid - half + cfg.shift))
    win_end   = int(np.floor(alu_mid + half + cfg.shift))
    REF_seq   = hg38_fa.fetch(chrom, win_start, win_end).upper()

    if mode == 'REF':
        ALT_seq = REF_seq

    elif mode == 'INS':
        # Insert consensus at alu_start; crop symmetrically back to cfg.MB
        var_rel_pos = alu_start - win_start
        ALT_seq     = REF_seq[:var_rel_pos] + consensus_seq + REF_seq[var_rel_pos:]

        ins_len   = len(consensus_seq)
        to_remove = ins_len / 2
        if to_remove == 0.5:
            ALT_seq = ALT_seq[1:]
        else:
            ALT_seq = ALT_seq[math.ceil(to_remove) : -math.floor(to_remove)]

    elif mode == 'DEL':
        del_len      = alu_end - alu_start
        to_add_left  = math.ceil(del_len / 2)
        to_add_right = math.floor(del_len / 2)

        ALT_start = win_start - to_add_left
        ALT_stop  = win_end   + to_add_right
        ALT_full  = hg38_fa.fetch(chrom, ALT_start, ALT_stop).upper()

        var_rel_pos_ALT = alu_start - ALT_start
        ALT_seq = ALT_full[:var_rel_pos_ALT] + ALT_full[var_rel_pos_ALT + del_len:]

        if len(ALT_seq) != cfg.MB:
            raise ValueError(f'DEL ALT sequence length {len(ALT_seq)} != {cfg.MB}')

    else:
        raise ValueError(f"mode must be 'REF', 'INS', or 'DEL', got {mode!r}")

    if cfg.revcomp:
        REF_seq = str(Seq(REF_seq).reverse_complement())
        ALT_seq = str(Seq(ALT_seq).reverse_complement())

    return REF_seq, ALT_seq

def load_model_params(cfg: Config) -> tuple:
    """
    Parse Akita params.json.
    Returns (params_model, target_length_cropped, bin_size_bp).
    """
    with open(cfg.params_file) as f:
        params = json.load(f)
    pm = params['model']
    pm['augment_shift'] = 0
    cropping   = pm['head_hic'][5]['cropping']
    tlen       = pm['target_length']
    tlen_crop  = tlen - 2 * cropping
    bin_size   = pm['seq_length'] // tlen

    print('bin_size', bin_size)
    return pm, tlen_crop, bin_size

############################## set up gpu 

def configure_gpu(cfg: Config, worker_id: int) -> None:
    """
    Cap TF GPU memory for this process to cfg.gpu_mem_mb (or 1/n_workers of
    total if not set).  Must be called before any TF graph is built.
    """
    import tensorflow as tf

    gpus = tf.config.list_physical_devices('GPU')
    if not gpus:
        print(f'[worker {worker_id}] No GPU — using CPU')
        return

    if cfg.gpu_mem_mb is not None:
        mem_mb = cfg.gpu_mem_mb
    else:
        try:
            out    = subprocess.check_output(
                ['nvidia-smi', '--query-gpu=memory.total', '--format=csv,noheader,nounits'],
                text=True)
            total  = int(out.strip().split('\n')[0])
        except Exception:
            total  = 16_000
        mem_mb = int(total / cfg.n_workers)

    try:
        tf.config.set_logical_device_configuration(
            gpus[0], [tf.config.LogicalDeviceConfiguration(memory_limit=mem_mb)])
        print(f'[worker {worker_id}] GPU limit: {mem_mb} MB')
    except RuntimeError:
        tf.config.experimental.set_memory_growth(gpus[0], True)
        print(f'[worker {worker_id}] Fell back to memory-growth mode')

def load_akita(params_model: dict, cfg: Config):
    """Load and restore the Akita SeqNN model."""
    sys.path.append(cfg.repo_path)
    from basenji import seqnn_gpu as seqnn
    print('using seqnn_gpu!')
    model = seqnn.SeqNN(params_model)
    model.restore(cfg.model_file)
    return model

### gradient and scoring 

def compute_nucleotide_grads(seq: str, seqnn_model) -> np.ndarray:
    """
    Compute per-nucleotide gradients for a ~1 Mb DNA sequence via GradientTape.

    Returns
    -------
    grads : np.ndarray [seq_len, 4] float32
        Collapse to saliency with np.sum(np.abs(grads), axis=1).
    """
    import tensorflow as tf
    from basenji import dna_io

    seq_1hot = np.expand_dims(dna_io.dna_1hot(seq).astype(np.float32), 0)
    seq_input   = tf.keras.Input(shape=(seqnn_model.seq_length, 4), dtype=tf.float32)
    temp_model  = tf.keras.Model(inputs=seq_input,
                                  outputs=seqnn_model.model(seq_input))
    seq_var = tf.Variable(seq_1hot, dtype=tf.float32)
    with tf.GradientTape() as tape:
        tape.watch(seq_var)
        loss = tf.reduce_sum(temp_model(seq_var))
    return tape.gradient(loss, seq_var)[0].numpy()   # [seq_len, 4]

def bin_saliency(grads: np.ndarray, bin_size: int, target_length_cropped: int) -> np.ndarray:
    """
    Sum |gradients| across bases, then aggregate into Akita bins.
    Skips the 32-bin crop on each edge to match Akita's output window.

    Returns np.ndarray [target_length_cropped] float32.
    """
    sal = np.sum(np.abs(grads), axis=1)   # [seq_len]
    binned = np.empty(target_length_cropped, dtype=np.float32)
    for i in range(target_length_cropped):
        vals = sal[(i + 32) * bin_size : (i + 33) * bin_size]
        vals = vals[~np.isnan(vals)]

        #replace with mean vs sum in case there are nans 
        binned[i] = np.mean(vals) if len(vals) > 0 else np.nan
    return binned

def align_alu_consensus(alu_seq: str, consensus_seq: str) -> float:
    """Global pairwise alignment score, scaled by len(alu_seq)."""
    alns = pairwise2.align.globalxx(alu_seq, consensus_seq, one_alignment_only=True)
    return alns[0].score / len(alu_seq)

#### variant computation 

def compute_variant(worker_id, idx, row, cfg, hg38_fa, alu_dfam,
                    seqnn_model, utils, scoring_utils,
                    tlen_crop, bin_size) -> dict:
    """
    Compute REF gradients and all Alu insertions for one variant.

    Returns a dict ready to be written to HDF5 + the scores DataFrame updates.
    Keys:
        'idx'          : variant index
        'chrom', 'alu_start', 'alu_end', 'alu_strand'
        'REF_pred'     : np.ndarray
        'REF_grads'    : np.ndarray [seq_len, 4]
        'REF_sal_bin'  : np.ndarray [tlen_crop]
        'insertions'   : list of dicts, one per Alu family (see below)
        'scores'       : dict of {col_name: value} for the scores DataFrame

    Each insertion dict contains:
        'alu_name', 'align_score', 'mse', 'spearman',
        'ALT_pred', 'ALT_grads', 'ALT_sal_bin', 'delta_sal_bin'
    """
    tag = f'[worker {worker_id}]'

    chrom      = row['CHROM']
    alu_start  = int(row['POS'])
    alu_end    = int(row['END'])
    alu_strand = row['strand']
    orig_idx   = row['orig_idx']
    alu_mid    = int(np.floor((alu_start + alu_end) / 2))
    alu_seq    = hg38_fa.fetch(chrom, alu_start, alu_end)
    revcomp    = (alu_strand == '-')

    # REF (no DEL)
    REF_seq, DEL_seq=get_ref_alt_sequences(hg38_fa, chrom, alu_start, alu_end, cfg, mode='DEL') 

    REF_pred    = utils.vector_from_seq(REF_seq)
    REF_grads   = compute_nucleotide_grads(REF_seq, seqnn_model)
    REF_sal_bin = bin_saliency(REF_grads, bin_size, tlen_crop)

    DEL_pred    = utils.vector_from_seq(DEL_seq)
    DEL_grads   = compute_nucleotide_grads(DEL_seq, seqnn_model)
    DEL_sal_bin = bin_saliency(DEL_grads, bin_size, tlen_crop)

    mse_REF_DEL           = scoring_utils.mse(REF_pred, DEL_pred)
    spearman_REF_DEL      = scoring_utils.spearman(REF_pred, DEL_pred)
    
    insertions = []
    scores     = {}

    for _, dfam_row in alu_dfam.iterrows():

        try:

            alu_name      = dfam_row['AluName']
            consensus_seq = dfam_row['Consensus_seq']

            if revcomp:
                #means on + strand, Alu is reversed. so we reverse. 
                consensus_seq = str(Seq(consensus_seq).reverse_complement())

            REF_seq, ALT_seq=get_ref_alt_sequences(hg38_fa, chrom, alu_start, alu_end, cfg, 
                                                   mode='INS', consensus_seq= consensus_seq)
            align_score = align_alu_consensus(alu_seq.upper(), consensus_seq.upper())

            ALT_pred = utils.vector_from_seq(ALT_seq)
            mse_REF_INS = scoring_utils.mse(REF_pred, ALT_pred)
            spearman_REF_INS = scoring_utils.spearman(REF_pred, ALT_pred)

            mse_DEL_INS = scoring_utils.mse(DEL_pred, ALT_pred)
            spearman_DEL_INS = scoring_utils.spearman(DEL_pred, ALT_pred)

            ALT_grads = compute_nucleotide_grads(ALT_seq, seqnn_model)
            ALT_sal_bin = bin_saliency(ALT_grads, bin_size, tlen_crop)
            delta_sal_bin = ALT_sal_bin - REF_sal_bin

            insertions.append(dict(
                alu_name=alu_name, align_score=align_score, consensus_len=len(consensus_seq),
                mse_REFDEL=mse_REF_DEL, corr_REFDEL=spearman_REF_DEL,
                mse_REFINS=mse_REF_INS, corr_REFINS=spearman_REF_INS,
                ALT_pred=ALT_pred, ALT_grads=ALT_grads,
                ALT_sal_bin=ALT_sal_bin, delta_sal_bin=delta_sal_bin,
            ))
            
            scores[f'mse_REFDEL_{alu_name}'] = mse_REF_DEL
            scores[f'corr_REFDEL_{alu_name}'] = spearman_REF_DEL
            scores[f'mse_REFINS_{alu_name}'] = mse_REF_INS
            scores[f'corr_REFINS_{alu_name}'] = spearman_REF_INS
            scores[f'mse_DELINS_{alu_name}'] = mse_DEL_INS
            scores[f'corr_DELINS_{alu_name}'] = spearman_DEL_INS
            scores[f'alignscore_{alu_name}']  = align_score

        except Exception as e:
            print(f'  {tag} ERROR inserting {alu_name} @ variant {orig_idx}: {e}')

    return dict(orig_idx=orig_idx, chrom=chrom, alu_start=alu_start, alu_end=alu_end,
                alu_strand=alu_strand, shift=cfg.shift,
                revcomp=cfg.revcomp,
                REF_pred=REF_pred, REF_grads=REF_grads,
                REF_sal_bin=REF_sal_bin, insertions=insertions, scores=scores)