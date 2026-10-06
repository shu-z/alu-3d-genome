#!/usr/bin/env python3
"""Run the pairwise-deletion scripts together, one process each."""
import argparse
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

HERE = Path(__file__).resolve().parent
SCRIPTS = ['pairs_alu_random.py', 'pairs_ctcf_ctcf.py', 'pairs_alu_ctcf.py']


def run(script):
    subprocess.run([sys.executable, str(HERE / script)], check=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scripts", "-s", nargs="+", default=SCRIPTS,
                        help="scripts to run, relative to this directory")
    parser.add_argument("--workers", "-w", type=int, default=len(SCRIPTS),
                        help="how many to run at once")
    args = parser.parse_args()

    with ProcessPoolExecutor(max_workers=args.workers) as executor:
        list(executor.map(run, args.scripts))


if __name__ == '__main__':
    main()
