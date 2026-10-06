
import argparse
import multiprocessing as mp
import os
import signal
import sys

import h5py
import numpy as np
import pandas as pd
import pysam

from insertion_gradients_lib import load_model_params, load_akita, configure_gpu, compute_variant
from config import Config

def _write_variant_h5(result: dict, out_dir: str, hdf_kwargs: dict) -> str:
    """Write one variant's result dict to its own HDF5 file. Returns the path."""

    orig_idx = result['orig_idx']
    path = os.path.join(out_dir, f'variant_{orig_idx}.h5')

    with h5py.File(path, 'w') as hdf:
        hdf.attrs['orig_idx']   = orig_idx
        hdf.attrs['chrom']      = result['chrom']
        hdf.attrs['alu_start']  = result['alu_start']
        hdf.attrs['alu_end']    = result['alu_end']
        hdf.attrs['alu_strand'] = str(result['alu_strand'])
        hdf.attrs['shift']      = result['shift']
        hdf.attrs['revcomp']    = result['revcomp']

        hdf.create_dataset('REF_pred',         data=result['REF_pred'].astype(np.float32),  **hdf_kwargs)
        hdf.create_dataset('REF_saliency_nt',  data=result['REF_grads'].astype(np.float16), **hdf_kwargs)
        hdf.create_dataset('REF_saliency_bin', data=result['REF_sal_bin'],                  **hdf_kwargs)

        for ins in result['insertions']:
            ig = hdf.require_group(f'insertion_{ins["alu_name"]}')
            ig.attrs['align_score']   = float(ins['align_score'])
            ig.attrs['consensus_len'] = float(ins['consensus_len'])
            ig.attrs['mse_REFDEL']    = float(ins['mse_REFDEL'])
            ig.attrs['corr_REFDEL']   = float(ins['corr_REFDEL'])
            ig.attrs['mse_REF_INS']   = float(ins['mse_REFINS'])
            ig.attrs['corr_REFINS']   = float(ins['corr_REFINS'])
            ig.create_dataset('ALT_pred',           data=ins['ALT_pred'].astype(np.float32),  **hdf_kwargs)
            ig.create_dataset('ALT_saliency_nt',    data=ins['ALT_grads'].astype(np.float16), **hdf_kwargs)
            ig.create_dataset('ALT_saliency_bin',   data=ins['ALT_sal_bin'],                  **hdf_kwargs)
            ig.create_dataset('delta_saliency_bin', data=ins['delta_sal_bin'],                **hdf_kwargs)

    return path

def _worker(worker_id: int, cfg: Config, result_queue: mp.Queue) -> None:
    """
    Spawned in a fresh process. Processes every cfg.n_workers-th variant
    (round-robin), writes each one to its own HDF5 file, and sends the
    score row back to the parent via result_queue.
    Sends None when finished so the parent knows this worker is done.
    """

    # Exit cleanly on Ctrl-C / SIGTERM from parent. Don't ignore — that
    # leaves orphaned children when the parent dies.
    def _bail(signum, frame):
        print(f'[worker {worker_id}] received signal {signum}, exiting')
        sys.exit(0)
    signal.signal(signal.SIGINT,  _bail)
    signal.signal(signal.SIGTERM, _bail)

    # GPU must be configured before any TF import
    configure_gpu(cfg, worker_id)

    sys.path.append(cfg.repo_path)
    import akita_utils_forplotting as utils
    import akita_utils_scoring as scoring_utils

    params_model, tlen_crop, bin_size = load_model_params(cfg)
    seqnn_model = load_akita(params_model, cfg)
    print(f'[worker {worker_id}] Akita loaded')

    hg38_fa  = pysam.FastaFile(cfg.hg38_fa)
    alu_dfam = pd.read_csv(cfg.dfam_csv, sep='\t', index_col=0)

    # prefiltered set
    scores_top = pd.read_csv(cfg.scores_csv, sep='\t', index_col=0)

    hdf_kwargs = dict(compression=cfg.hdf5_compression,
                      compression_opts=cfg.hdf5_compression_opts)

    chunk = scores_top.iloc[worker_id::cfg.n_workers]
    print(f'[worker {worker_id}] {len(chunk)} variants assigned')

    for idx, row in chunk.iterrows():
        print(f'[worker {worker_id}] variant {idx}')
        try:
            result = compute_variant(
                worker_id, idx, row, cfg,
                hg38_fa, alu_dfam, seqnn_model,
                utils, scoring_utils, tlen_crop, bin_size,
            )
            path = _write_variant_h5(result, cfg.out_h5_dir, hdf_kwargs)
            print(f'[worker {worker_id}] wrote {path}')
            # Send only the score row back — the heavy arrays are on disk.
            result_queue.put({'idx': idx, 'scores': result['scores']})
        except Exception as e:
            print(f'[worker {worker_id}] ERROR on variant {idx}: {e}')

    result_queue.put(None)   # signal this worker is done
    print(f'[worker {worker_id}] finished')

def run(cfg: Config, worker_id: int | None = None) -> None:
    """
    Launch workers and collect per-variant score rows.

    If worker_id is given, run only that one worker in the current process —
    useful for debugging or SLURM arrays.
    """
    scores_top = pd.read_csv(cfg.scores_csv, sep='\t', index_col=0)

    os.makedirs(cfg.out_h5_dir, exist_ok=True)

    hdf_kwargs = dict(compression=cfg.hdf5_compression,
                      compression_opts=cfg.hdf5_compression_opts)

    if worker_id is not None:
        # ── single-worker debug mode
        configure_gpu(cfg, worker_id)
        sys.path.append(cfg.repo_path)
        import akita_utils_forplotting as utils
        import akita_utils_scoring as scoring_utils

        params_model, tlen_crop, bin_size = load_model_params(cfg)
        seqnn_model = load_akita(params_model, cfg)
        hg38_fa     = pysam.FastaFile(cfg.hg38_fa)
        alu_dfam    = pd.read_csv(cfg.dfam_csv, sep='\t', index_col=0)

        chunk = scores_top.iloc[worker_id::cfg.n_workers]
        for idx, row in chunk.iterrows():
            result = compute_variant(
                worker_id, idx, row, cfg,
                hg38_fa, alu_dfam, seqnn_model,
                utils, scoring_utils, tlen_crop, bin_size,
            )
            _write_variant_h5(result, cfg.out_h5_dir, hdf_kwargs)
            for col, val in result['scores'].items():
                scores_top.loc[idx, col] = val

    else:
        # ── full parallel mode
        # 'spawn' ensures TF initialises fresh in each child so that the
        # GPU memory cap takes effect before any graph is built.
        mp.set_start_method('spawn', force=True)
        result_queue = mp.Queue()

        workers = [
            mp.Process(target=_worker, args=(wid, cfg, result_queue),
                       daemon=True)
            for wid in range(cfg.n_workers)
        ]
        for w in workers:
            w.start()

        def _shutdown():
            """Terminate every still-living worker, then kill stragglers."""
            for w in workers:
                if w.is_alive():
                    w.terminate()
            for w in workers:
                w.join(timeout=5)
                if w.is_alive():
                    w.kill()
                    w.join()

        try:
            done_count = 0
            while done_count < cfg.n_workers:
                # Timeout lets us notice if every worker has died without
                # sending its sentinel (e.g. segfault, OOM kill).
                try:
                    item = result_queue.get(timeout=1.0)
                except mp.queues.Empty:
                    if not any(w.is_alive() for w in workers):
                        print('All workers exited without finishing.')
                        break
                    continue

                if item is None:
                    done_count += 1
                else:
                    idx = item['idx']
                    for col, val in item['scores'].items():
                        scores_top.loc[idx, col] = val
        except KeyboardInterrupt:
            print('\nInterrupted — shutting down workers...')
            _shutdown()
            raise
        except Exception:
            _shutdown()
            raise
        else:
            for w in workers:
                w.join()

    print('writing scores...')
    scores_top.to_csv(cfg.out_scores, sep='\t')
    print(f'\nDone.\n  HDF5 dir → {cfg.out_h5_dir}\n  scores   → {cfg.out_scores}')

####### CLI

def _parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description='Alu insertion gradient scorer')
    p.add_argument('--n_workers',  type=int,   default=None)
    p.add_argument('--gpu_mem_mb', type=int,   default=None,
                   help='Hard GPU memory cap per worker in MB. '
                        'Default: total_gpu_mb / n_workers.')
    p.add_argument('--worker_id',  type=int,   default=None,
                   help='Debug: run only this worker in the current process.')
    p.add_argument('--shift', type=int, default=None,
                   help='Slide the 1 Mb window by this many bp before prediction. '
                        'Positive = right (downstream), negative = left (upstream).')
    p.add_argument('--revcomp', action='store_true', default=None,
                   help='RC every sequence (REF and ALT) before passing to Akita.')
    p.add_argument('--out_h5_dir', type=str,   default=None,
                   help='Directory to write per-variant HDF5 files into.')
    p.add_argument('--out_scores', type=str,   default=None)
    return p.parse_args()

if __name__ == '__main__':
    args = _parse_args()
    overrides = {k: v for k, v in vars(args).items()
                 if v is not None and k != 'worker_id'}
    cfg = Config(**overrides)
    run(cfg, worker_id=args.worker_id)


def main():
    """
    run.py  —  parallel orchestration with per-variant HDF5 output.

    Each variant is written to its own file: <out_dir>/variant_<orig_idx>.h5
    Workers write directly to disk (no queue), so there is no central writer
    bottleneck. The parent only spawns workers and collects per-variant score
    rows for the final TSV.

    Usage
    ─────
    Full parallel run (default 6 workers):
        python run.py

    Fewer workers or a GPU memory override:
        python run.py --n_workers 4 --gpu_mem_mb 2500

    Single-worker mode for debugging:
        python run.py --worker_id 0 --n_workers 6
    """


if __name__ == '__main__':
    main()
