

"""Run SuPreMo over a sharded input set, queueing jobs on one GPU.

The same queue drives every sharded scoring run in this directory; the shard
naming and output location are the only things that change between them. The
four runs used in the manuscript were:

    # 300 bp random deletions, 10 kb and 100 kb spacing
    run_supremo_sharded.py -d $RES/300bp_random/10kb/  -p alu_300bp_10kb_random  -n 20
    run_supremo_sharded.py -d $RES/300bp_random/100kb/ -p alu_300bp_100kb_random -n 20

    # rolling-window mutagenesis, top-100-per-subfamily and the broader 5100 set
    run_supremo_sharded.py -d $RES/20260413_alu_rollingwindow_top100/ \
                           -p alu_top100_repName_rolling -n 12
    run_supremo_sharded.py -d $RES/20260413_alu_rollingwindow_5100/ \
                           -p alu_5100_repName_rolling -n 12

where RES is <project>/results/paper_results.
"""
import argparse
import csv
import os
import subprocess
import time

import psutil  

# Paths - edit for your environment
SUPREMO_DIR = os.environ.get("SUPREMO_DIR", "/pollard/home/szhang20/akita_variant_scoring")
AKITA_DIR = os.environ.get("AKITA_DIR", "/pollard/home/szhang20/akita")
PROJECT_DIR = os.environ.get("ALU_PROJECT_DIR", "/pollard/home/szhang20/alu")

########################################################
GPU_ID = int(os.environ.get("CUDA_VISIBLE_DEVICES", 0))
MAX_UTIL = 80
MAX_MEM_MB = 78000
MAX_JOBS = 11
CHECK_INTERVAL = 120
SCORE_DIR = f"{PROJECT_DIR}/results/paper_results/300bp_random/10kb/"
LOG_DIR = f"{SCORE_DIR}logs/"
INPUT_PREFIX = "alu_300bp_10kb_random"
N_SHARDS = 20
FA_PATH = f"{SUPREMO_DIR}/data/hg38.fa"
SUPREMO_PATH = f"{SUPREMO_DIR}/scripts/SuPreMo.py"


########################################################

def get_gpu_status(gpu_id=0):
    cmd = [
        "nvidia-smi",
        "--query-gpu=utilization.gpu,memory.used",
        "--format=csv,noheader,nounits",
    ]
    out = subprocess.check_output(cmd).decode("utf-8").strip()
    util, mem = map(int, out.splitlines()[gpu_id].split(", "))
    return util, mem

def gpu_ready(util, mem):
    return util < MAX_UTIL and mem < MAX_MEM_MB

def get_process_mem_mb(pid):
    """Return RSS (MB) of pid + all its children."""
    try:
        proc = psutil.Process(pid)
        rss = proc.memory_info().rss
        for child in proc.children(recursive=True):
            try:
                rss += child.memory_info().rss
            except:
                pass
        return rss / (1024**2)
    except psutil.NoSuchProcess:
        return 0

########################################################

def start_job(script_args, script_id):
    script_name = os.path.basename(script_args[2])
    log_path = os.path.join(LOG_DIR, f"{script_name}.log")
    #csv_path = os.path.join(LOG_DIR, f"{script_name}_gpu.csv")

    f_log = open(log_path, "w")
    #f_csv = open(csv_path, "w", newline="")
    #csv_writer = csv.writer(f_csv)

    # csv_writer.writerow([
    #     "timestamp", "gpu_util_pct", "gpu_mem_mb",
    #     "proc_mem_mb"
    # ])

    p = subprocess.Popen(
        script_args,
        stdout=f_log,
        stderr=subprocess.STDOUT,
        env={**os.environ, "CUDA_VISIBLE_DEVICES": str(GPU_ID)},
    )

    return {
        "process": p,
        "stdout": f_log,
        #"csv_writer": csv_writer,
        #"csv_file": f_csv,
        "script_id": script_id,
        "script": " ".join(script_args),
        "start_time": time.time(),

        # stats
        "util_sum": 0,
        "mem_sum": 0,
        "proc_mem_sum": 0,
        "count": 0,
        "util_max": 0,
        "mem_max": 0,
        "proc_mem_max": 0,
    }

########################################################



########################################################
# Wait for remaining jobs
########################################################


########################################################


def main():
    global GPU_ID, LOG_DIR, MAX_UTIL, MAX_MEM_MB  # rebound below, read by helpers
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--score-dir", "-d", default=SCORE_DIR,
                        help="directory holding input/ and output/ for this run")
    parser.add_argument("--input-prefix", "-p", default=INPUT_PREFIX,
                        help="shard basename; '_<n>.txt' is appended")
    parser.add_argument("--shards", "-n", type=int, default=N_SHARDS,
                        help="number of input shards to score")
    parser.add_argument("--log-dir", default=None,
                        help="directory for job logs (default: <score-dir>/logs/)")
    parser.add_argument("--fasta", "-f", default=FA_PATH, help="reference FASTA")
    parser.add_argument("--supremo", "-s", default=SUPREMO_PATH, help="path to SuPreMo.py")
    parser.add_argument("--gpu-id", "-g", type=int, default=GPU_ID, help="GPU index to run on")
    parser.add_argument("--max-jobs", type=int, default=MAX_JOBS, help="maximum concurrent jobs")
    parser.add_argument("--max-util", type=float, default=MAX_UTIL,
                        help="hold off launching above this GPU utilisation percent")
    parser.add_argument("--max-mem", type=int, default=MAX_MEM_MB,
                        help="hold off launching above this GPU memory use, in MB")
    parser.add_argument("--check-interval", type=int, default=CHECK_INTERVAL,
                        help="seconds between GPU checks")
    args = parser.parse_args()

    GPU_ID = args.gpu_id
    MAX_UTIL = args.max_util
    MAX_MEM_MB = args.max_mem
    score_dir = args.score_dir
    LOG_DIR = args.log_dir or f"{score_dir}logs/"
    max_jobs, check_interval = args.max_jobs, args.check_interval

    os.makedirs(LOG_DIR, exist_ok=True)
    running = []
    finished = []
    for script_id in range(args.shards):
        filepath = f"{score_dir}input/{args.input_prefix}_{script_id}.txt"
        filename = f"{args.input_prefix}_{script_id}"

        script_args = [
            "python",
            args.supremo,
            filepath,
            "--file", filename,
            "--dir", f"{score_dir}output/",
            "--get_Akita_scores",
            "--genome", "hg38",
            "--fa", args.fasta,
            "--shift_by", "0", "1", "-1",
            "--revcomp", "add_revcomp"
        ]

        # Wait for slot
        while True:
            gpu_util, gpu_mem = get_gpu_status(GPU_ID)

            new_running = []
            for job in running:
                p = job["process"]
                proc_mem = get_process_mem_mb(p.pid)

                # update stats
                job["util_sum"] += gpu_util
                job["mem_sum"] += gpu_mem
                job["proc_mem_sum"] += proc_mem
                job["count"] += 1
                job["util_max"] = max(job["util_max"], gpu_util)
                job["mem_max"] = max(job["mem_max"], gpu_mem)
                job["proc_mem_max"] = max(job["proc_mem_max"], proc_mem)

                # log
                ts = time.time()
                # job["stdout"].write(
                #     f"{ts:.2f}\tGPU util: {gpu_util}%\t"
                #     f"GPU mem: {gpu_mem}MB\t"
                #     f"PROC mem: {proc_mem:.1f}MB\n"
                # )
                # job["stdout"].flush()

                # job["csv_writer"].writerow([ts, gpu_util, gpu_mem, proc_mem])
                # job["csv_file"].flush()

                if p.poll() is None:
                    new_running.append(job)
                else:
                    # finished
                    duration = time.time() - job["start_time"]
                    avg_util = job["util_sum"] / job["count"]
                    avg_mem = job["mem_sum"] / job["count"]
                    avg_proc_mem = job["proc_mem_sum"] / job["count"]

                    job["stdout"].write(
                        f"Finished in {duration:.2f}s | "
                        f"GPU util avg={avg_util:.1f}%, max={job['util_max']}% | "
                        f"GPU mem avg={avg_mem:.1f}MB, max={job['mem_max']}MB | "
                        f"PROC mem avg={avg_proc_mem:.1f}MB, "
                        f"max={job['proc_mem_max']:.1f}MB\n"
                    )
                    job["stdout"].close()
                    #job["csv_file"].close()

                    finished.append((job["script"], duration))

            running = new_running

            if len(running) < max_jobs and gpu_ready(gpu_util, gpu_mem):
                break

            time.sleep(check_interval)

        job = start_job(script_args, script_id)
        print(f"Started job: {filename}")
        running.append(job)
    while running:
        time.sleep(check_interval)
        gpu_util, gpu_mem = get_gpu_status(GPU_ID)

        new_running = []
        for job in running:
            p = job["process"]
            proc_mem = get_process_mem_mb(p.pid)

            job["util_sum"] += gpu_util
            job["mem_sum"] += gpu_mem
            job["proc_mem_sum"] += proc_mem
            job["count"] += 1
            job["util_max"] = max(job["util_max"], gpu_util)
            job["mem_max"] = max(job["mem_max"], gpu_mem)
            job["proc_mem_max"] = max(job["proc_mem_max"], proc_mem)

            ts = time.time()
            # job["stdout"].write(
            #     f"{ts:.2f}\tGPU util: {gpu_util}%\t"
            #     f"GPU mem: {gpu_mem}MB\tPROC mem: {proc_mem:.1f}MB\n"
            # )
            # job["stdout"].flush()
            # job["csv_writer"].writerow([ts, gpu_util, gpu_mem, proc_mem])
            # job["csv_file"].flush()

            if p.poll() is None:
                new_running.append(job)
            else:
                duration = time.time() - job["start_time"]
                avg_util = job["util_sum"] / job["count"]
                avg_mem = job["mem_sum"] / job["count"]
                avg_proc_mem = job["proc_mem_sum"] / job["count"]

                job["stdout"].write(
                    f"Finished in {duration:.2f}s | "
                    f"GPU util avg={avg_util:.1f}%, max={job['util_max']}% | "
                    f"GPU mem avg={avg_mem:.1f}MB, max={job['mem_max']}MB | "
                    f"PROC mem avg={avg_proc_mem:.1f}MB, "
                    f"max={job['proc_mem_max']:.1f}MB\n"
                )
                job["stdout"].close()
                # job["csv_file"].close()
                finished.append((job["script"], duration))

        running = new_running
    with open(os.path.join(LOG_DIR, "timing.log"), "w") as f:
        for name, t in finished:
            f.write(f"{name}: {t:.2f} sec\n")
    print("\nAll jobs finished!")


if __name__ == '__main__':
    main()
