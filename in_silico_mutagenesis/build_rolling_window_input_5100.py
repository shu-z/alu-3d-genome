
import os

import pandas as pd

# Paths - edit for your environment
PROJECT_DIR = os.environ.get("ALU_PROJECT_DIR", "/pollard/home/szhang20/alu")


def main():
    """
    Rebuild the rolling-window input shards for the 5100-Alu sample.

    Fixes a copy-paste bug in 20250519_Alu_sequencestructure_windowmutate.ipynb
    (cell 12): it was supposed to build rolling windows over `alu5100` (a random
    sample of 100 Alus per repName, alu.groupby('repName').sample(n=100,
    random_state=729), saved as alu5100_annot.txt in alu_resample_final.ipynb)
    but instead re-used `alu_top100_repName` (the top 100 *by mse_mean* per
    repName) left over from the cell above it -- so the "5100" rolling-window
    run was actually just a duplicate of the "top100" run on identical input.

    Mirrors cells 12-13 of that notebook, with alu5100 substituted in.
    """
    res_dir = f'{PROJECT_DIR}/results/paper_results/'
    alu5100 = pd.read_csv(f'{res_dir}20260413_alu_sample/alu5100_annot.txt', sep='\t', index_col=0)
    print(f'alu5100: {len(alu5100)} rows, {alu5100["repName"].nunique()} repNames')
    window = 30
    stride = 10
    output = []
    for idx, row in alu5100.iterrows():
        chrom = row["CHROM"]
        start = row["POS"]
        end = row["END"]
        orig_idx = row["orig_idx"]

        for window_start in range(start, end - window + 1, stride):
            window_end = window_start + window
            output.append({
                "CHROM": chrom,
                "POS": window_start,
                "REF": '-',
                "ALT": '-',
                "END": window_end,
                "SVTYPE": 'DEL',
                "SVLEN": window_end - window_start,
                "orig_idx": orig_idx,
            })
    windows_df = pd.DataFrame(output)
    print(f'windows_df: {len(windows_df)} rows')
    out_dir = f'{res_dir}20260413_alu_rollingwindow_5100/'
    os.makedirs(out_dir, exist_ok=False)
    os.makedirs(f'{out_dir}input/')
    os.makedirs(f'{out_dir}output/')
    os.makedirs(f'{out_dir}logs/')
    chunk_size = 10000
    n_chunks = 0
    for i, start in enumerate(range(0, len(windows_df), chunk_size)):
        print(i, start)
        chunk = windows_df[start:start + chunk_size]
        chunk.to_csv(f"{out_dir}input/alu_5100_repName_rolling_{i}.txt", sep='\t', index=False)
        n_chunks += 1
    windows_df.to_csv(f'{out_dir}alu_5100_repName_window30_shift10.txt', sep='\t', index=False, header=True)
    print(f'Wrote {n_chunks} shards -> {out_dir}input/')


if __name__ == '__main__':
    main()
