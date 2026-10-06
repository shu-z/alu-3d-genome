"""Sequence-composition mutagenesis of sampled Alus, scored with SuPreMo.

Covers all three mutation types in one pass. Each is a SuPreMo flag paired with
a "<up_flank> <down_flank> <percent>" argument, built by make_command() below:

    GC   --gc_mutate    mutates G/C bases, GC_AT_PERCENTS
    AT   --at_mutate    mutates A/T bases, GC_AT_PERCENTS
    CpG  --CpG_mutate   adds CpGs as a percentage of the current count,
                        CPG_PERCENTS (larger, since CpG counts are small)

A fourth job per input, --shuffle_mutate, shuffles the Alu body as a control.

Pick types with --mutations, levels with --gc-at-percent / --cpg-percent, and how
much flanking sequence is mutated with --windows (0 means the Alu body only).
"""
import argparse
from multiprocessing import Pool
import os
import subprocess

# Paths - edit for your environment
SUPREMO_DIR = os.environ.get("SUPREMO_DIR", "/pollard/home/szhang20/akita_variant_scoring")
AKITA_DIR = os.environ.get("AKITA_DIR", "/pollard/home/szhang20/akita")
PROJECT_DIR = os.environ.get("ALU_PROJECT_DIR", "/pollard/home/szhang20/alu")

alu_dir  = f'{PROJECT_DIR}/'
res_dir  = f'{alu_dir}results/paper_results/'
out_dir  = f'{res_dir}20260504_GCmutate/'
fa_path  = f'{SUPREMO_DIR}/data/hg38.fa'
supremo  = f'{SUPREMO_DIR}/scripts/SuPreMo_GCmutate.py'

file_dict = {
    'alu_high': f'{res_dir}20260413_alu_sample/aluhigh_1000samp_supremo.txt',
    'alu_low':  f'{res_dir}20260413_alu_sample/alulow_1000samp_supremo.txt',
}

WINDOWS         = [0, 300, 1000, 10000, 100000]           # 0 = mutate variant body only
GC_AT_PERCENTS  = [10, 20, 30, 40, 50]
CPG_PERCENTS    = [50, 100]

def make_command(args, file_path, out_prefix, flag, mut_arg):
    """Build a single SuPreMo invocation."""
    return [
        "python", args.supremo,
        file_path,
        "--file", out_prefix, "--dir", args.out_dir,
        "--get_Akita_scores",
        "--shift", "0", "--revcomp", "add_revcomp",
        "--genome", "hg38", "--fa", args.fasta,
        flag, mut_arg,
    ]

def directions_for_window(window):
    """
    Return a list of (up_flank, down_flank, mutate_variant_flag, suffix) tuples
    describing what to mutate for this window.

    window == 0     -> mutate the variant body only, no flanks
    window == 300   -> mutate upstream-only and downstream-only (two jobs)
    other windows   -> mutate both flanks symmetrically (one job)
    """
    if window == 0:
        # '1' is the mutate_variant flag; up/down flanks both 0
        return [(0, 0, 1, 'variant')]
    if window == 300:
        return [
            (window, 0, 0, f'{window}windowup'),
            (0, window, 0, f'{window}windowdown'),
        ]
    return [(window, window, 0, f'{window}window')]

def build_commands(args):
    jobs = []
    kinds = {'GC': '--gc_mutate', 'AT': '--at_mutate'}
    for file_name, file_path in zip(args.names, args.inputs):

        # GC and AT mutation jobs
        for window in args.windows:
            for mut_percent in args.gc_at_percent:
                for up, down, mut_var, suffix in directions_for_window(window):
                    mut_opt = f'{up} {down} {mut_percent} {mut_var}'
                    for kind, flag in [(k, v) for k, v in kinds.items() if k in args.mutations]:
                        out_prefix = f'{file_name}_{suffix}_{mut_percent}_{kind}'
                        label = f"{file_name} window={suffix} mut={mut_percent} {kind}"
                        jobs.append((label, make_command(args, file_path, out_prefix, flag, mut_opt)))

        # CpG-boost jobs (different percent grid)
        for window in (args.windows if 'CpG' in args.mutations else []):
            for mut_percent in args.cpg_percent:
                for up, down, mut_var, suffix in directions_for_window(window):
                    mut_opt = f'{up} {down} {mut_percent} {mut_var}'
                    out_prefix = f'{file_name}_{suffix}_{mut_percent}_CpG'
                    label = f"{file_name} window={suffix} mut={mut_percent} CpG"
                    jobs.append((label, make_command(args, file_path, out_prefix, '--CpG_mutate', mut_opt)))

        # Shuffle: one job per file, no window concept
        if 'shuffle' in args.mutations:
            out_prefix = f'{file_name}_alu_shuff'
            label = f"{file_name} mut=shuff"
            jobs.append((label, make_command(args, file_path, out_prefix, '--shuffle_mutate', '0 0 1')))

    return jobs

def run_one(job):
    label, command = job
    try:
        subprocess.run(command, check=True)
        return (label, True, None)
    except subprocess.CalledProcessError as e:
        return (label, False, f"exit {e.returncode}")
    except Exception as e:
        return (label, False, repr(e))

N_WORKERS = 9


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--inputs", "-i", nargs="+", default=list(file_dict.values()),
                        help="SuPreMo input tables to mutate")
    parser.add_argument("--names", "-n", nargs="+", default=list(file_dict),
                        help="short name per input, used as the output prefix")
    parser.add_argument("--out-dir", "-o", default=out_dir, help="directory for SuPreMo output")
    parser.add_argument("--fasta", "-f", default=fa_path, help="reference FASTA")
    parser.add_argument("--supremo", "-s", default=supremo,
                        help="path to SuPreMo_GCmutate.py")
    parser.add_argument("--mutations", "-m", nargs="+", default=['GC', 'AT', 'CpG', 'shuffle'],
                        choices=['GC', 'AT', 'CpG', 'shuffle'], help="which mutation types to run")
    parser.add_argument("--windows", "-w", type=int, nargs="+", default=WINDOWS,
                        help="flanking window sizes in bp; 0 mutates the Alu body only")
    parser.add_argument("--gc-at-percent", type=int, nargs="+", default=GC_AT_PERCENTS,
                        help="GC and AT mutation levels, as percentages")
    parser.add_argument("--cpg-percent", type=int, nargs="+", default=CPG_PERCENTS,
                        help="CpG boost levels, as a percentage of the current CpG count")
    parser.add_argument("--workers", type=int, default=N_WORKERS, help="parallel worker processes")
    parser.add_argument("--dry-run", action="store_true",
                        help="list the jobs that would run, then exit")
    args = parser.parse_args()

    if len(args.names) != len(args.inputs):
        parser.error("--names must have one entry per --inputs")

    jobs = build_commands(args)
    if args.dry_run:
        for label, command in jobs:
            print(label)
        print(f"{len(jobs)} jobs")
        return

    os.makedirs(args.out_dir, exist_ok=True)
    print(f"Running {len(jobs)} jobs across {args.workers} workers")
    with Pool(processes=args.workers) as pool:
        for label, ok, err in pool.imap_unordered(run_one, jobs):
            status = "OK " if ok else "FAIL"
            print(f"[{status}] {label}" + (f" \u2014 {err}" if err else ""))


if __name__ == "__main__":
    main()