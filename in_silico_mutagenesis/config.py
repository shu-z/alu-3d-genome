"""
config.py  —  all paths and tunable constants in one place.
Edit this file; nothing else should need to change for a new run.
"""

import os
from dataclasses import dataclass

# Paths - edit for your environment
SUPREMO_DIR = os.environ.get("SUPREMO_DIR", "/pollard/home/szhang20/akita_variant_scoring")
AKITA_DIR = os.environ.get("AKITA_DIR", "/pollard/home/szhang20/akita")
PROJECT_DIR = os.environ.get("ALU_PROJECT_DIR", "/pollard/home/szhang20/alu")
DATA_DIR = os.environ.get("ALU_DATA_DIR", "/pollard/data/projects/shzhang")
POLLARD_DATA = os.environ.get("POLLARD_DATA", "/pollard/data")

@dataclass
class Config:
    # parallelism
    n_workers: int = 2
    # Hard GPU memory limit per worker (MB). None → auto: total_gpu_mb / n_workers.
    # Reduce manually (e.g. gpu_mem_mb=2500) if workers hit OOM.
    gpu_mem_mb: int = 10000

    # inputs 
    repo_path:   str = f'{SUPREMO_DIR}/'
    model_file:  str = f'{SUPREMO_DIR}/Akita_model/model_best.h5'
    params_file: str = f'{SUPREMO_DIR}/Akita_model/params.json'
    hg38_fa:     str = f'{SUPREMO_DIR}/data/hg38.fa'
    dfam_csv:    str = f'{DATA_DIR}/alus/data/dfam/alu_dfam_consensus.txt'
    scores_csv:  str = (f'{PROJECT_DIR}/bin/ins_gradients/20260519_topMSE/top_ins_noPDCH.txt')

    # outputs 
    out_dir=f'{PROJECT_DIR}/bin/ins_gradients/20260519_topMSE/'
    out_h5_dir:     str = f'{out_dir}h5_out/'
    out_scores: str = f'{out_dir}top_ins_noPDCH_scores.txt'

    # sequence augmentation options 

    revcomp: bool = False

    # shift: slide the 1 Mb window by this many bp before prediction.
    #   Positive → shift right (fetch further downstream),
    #   negative → shift left (fetch further upstream).
    #   Both REF and ALT windows are shifted by the same amount so the
    #   insertion stays centred relative to the window.
    shift: int = 0

    # akita model 
    MB:                    int = 1_048_576
    hdf5_compression:      str = 'gzip'
    hdf5_compression_opts: int = 4