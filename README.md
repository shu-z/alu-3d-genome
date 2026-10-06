# *Alu* Elements as Modulators of 3D Genome Folding

Code for a genome-wide *in silico* deletion screen of ~1.1 million *Alus*, scored with a
sequence-to-structure deep learning model, plus the annotation, perturbation and figure
code used to characterise the high-scoring subset.

Preprint: [biorxiv.org/content/10.64898/2026.09.16.752217v1](https://www.biorxiv.org/content/10.64898/2026.09.16.752217v1)

## Installation

Scoring uses [SuPreMo](https://github.com/pollard-lab/SuPreMo), wrapping the
[Akita](https://github.com/calico/basenji/tree/master/manuscripts/akita) contact-map model;
some analyses use [Enformer](https://github.com/deepmind/deepmind-research/tree/master/enformer).
Clone SuPreMo/Akita and Basenji separately and obtain the model weights, `hg38.fa` and
`CHM13.fa`.

Python 3.9 with `tensorflow` 2.8, `numpy` 1.21, `pandas` 1.4, `scipy` 1.7, `scikit-learn`
1.6, `matplotlib` 3.5, `h5py` 3.14, `cooltools` 0.5, `cooler` 0.8, `pysam` 0.19,
`biopython` 1.80, `pyBigWig` 0.3, `pybedtools`, `pyranges` 0.1, `shap` 0.48, `astropy` 5.0.
GPU scoring needs a CUDA runtime.

Each script has a `Paths` block at the top. Every root reads an environment variable and
falls back to the path it was run from:

```bash
export ALU_PROJECT_DIR=/path/to/alu                 # inputs, results, figures
export ALU_DATA_DIR=/path/to/data/projects/<user>   # bulk annotation data
export SUPREMO_DIR=/path/to/akita_variant_scoring   # SuPreMo and the Akita model
export AKITA_DIR=/path/to/akita                     # Akita checkout
export POLLARD_DATA=/path/to/data                   # shared mount (UCSC tracks)
export SLURM_DIR=/path/to/slurm                     # job logs
```

Shared helpers are in [utils/](utils). `akita_utils_scoring` and `akita_utils_forplotting`
are vendored from the `akita_variant_scoring` fork, not the public SuPreMo release; scripts
resolve them from `SUPREMO_DIR` at run time.

## Instructions

Scripts default to the values used for the manuscript. Pass `--help` for options.

### Worked example: scoring 10 *Alus*

[example/](example) scores ten hg38 *Alus* under all six augmentations (shifts of 0, +1,
-1 bp, each forward and reverse complemented).

```bash
bash example/run_supremo_example.sh
```

Scores and maps for the ten are checked in under
[example/expected_output/](example/expected_output), so
[combine_scores_and_plot.ipynb](example/combine_scores_and_plot.ipynb) runs without the
model. It collapses the six augmentations to a mean and SD per *Alu* and plots the scores
and contact maps.

### Genomic context annotation

Annotate scored *Alus* with chromatin and sequence features at the *Alu* and in 1 kb,
10 kb, 100 kb and 1 Mb windows.

```bash
python feature_annotation/count_features.py elements --window-sizes 1000 10000 100000
python feature_annotation/count_features.py chromhmm
python feature_annotation/count_features.py ctcf
python feature_annotation/count_features.py wgbs
python feature_annotation/count_features.py bigwig       # phyloP / phastCons / mappability
python feature_annotation/count_features.py bigwig_mp    # same, parallelised for the 1Mb window

python feature_annotation/overlap_features.py tracks
python feature_annotation/overlap_features.py chip --encode-dir /path/to/encode/ChIP/ --workers 8

python feature_annotation/average_features.py gc
python feature_annotation/average_features.py bigwig

python feature_annotation/20250206_alu_feature_elasticnet.py   # elastic net and logistic regression
bash feature_annotation/run_enformer_alu.sh                    # Enformer tracks
```

### Targeted *in silico* perturbation

**Sequence composition.** GC, AT and CpG mutation plus a shuffled control. Defaults
reproduce the 146 jobs used for the manuscript.

```bash
python in_silico_mutagenesis/mutate_composition.py                    # all 146 jobs
python in_silico_mutagenesis/mutate_composition.py --dry-run          # list, do not run
python in_silico_mutagenesis/mutate_composition.py -m GC --windows 0  # one type, Alu body only
python in_silico_mutagenesis/build_shuffled_controls.py
```

**Rolling-window mutagenesis.** Slide a mutated window across each *Alu*.

```bash
python in_silico_mutagenesis/build_rolling_window_input.py
python in_silico_mutagenesis/build_rolling_window_input_5100.py
python in_silico_mutagenesis/run_supremo_sharded.py \
    -d $RES/20260413_alu_rollingwindow_top100/ -p alu_top100_repName_rolling -n 12
```

**Random deletion controls.** Random 300 bp deletions at matched distances.

```bash
python in_silico_mutagenesis/build_random_300bp_input.py
python in_silico_mutagenesis/run_supremo_sharded.py \
    -d $RES/300bp_random/10kb/ -p alu_300bp_10kb_random -n 20
python in_silico_mutagenesis/run_supremo_sharded.py \
    -d $RES/300bp_random/100kb/ -p alu_300bp_100kb_random -n 20
```

**Transcription factor motifs.** Motif deletion, chiefly ARID3B, and binned saliency.

```bash
python in_silico_mutagenesis/arid3b_motif_deletion.py
python in_silico_mutagenesis/arid3b_binned_saliency.py
```

**Gene swaps.** Move *Alus* into and out of gene bodies and introns.

```bash
python in_silico_mutagenesis/gene_swaps_parallel.py
python in_silico_mutagenesis/gene_intron_swaps.py
python in_silico_mutagenesis/gene_overlap_null.py
```

**Consensus insertions and *Alu*–CTCF pairs.**

```bash
python in_silico_mutagenesis/insert_consensus_alu.py
python in_silico_mutagenesis/insert_consensus_alu_parallel.py
python in_silico_mutagenesis/run_insertion_gradients.py

python in_silico_mutagenesis/run_pairs_all_modes.py
python in_silico_mutagenesis/pairs_alu_ctcf.py
python in_silico_mutagenesis/pairs_ctcf_ctcf.py
python in_silico_mutagenesis/pairs_alu_random.py
```

[blank_canvas_ctcf.ipynb](in_silico_mutagenesis/blank_canvas_ctcf.ipynb) builds random
sequence at a chosen GC content with CTCF motifs inserted in either orientation.

### Figures

Scripts in [figures/](figures) are named for the panel they draw.

```bash
Rscript figures/fig1_all_alu_scores.R        # genome-wide score distributions
Rscript figures/fig1_karyotype.R             # disruption across chromosomes
Rscript figures/fig2_1KG_SV.R                # polymorphic Alus
Rscript figures/fig3_elasticnet.R            # feature model coefficients
Rscript figures/fig4_gene_swaps.R            # gene swaps
Rscript figures/fig5_gc_mutagenesis.R        # GC/AT/CpG mutagenesis
Rscript figures/fig5_random_300bp.R          # random deletion controls
Rscript figures/fig5_rolling_window.R        # rolling-window profiles
Rscript figures/fig6_consensus_insertions.R  # consensus-Alu insertions

Rscript figures/phylo_tree.R                 # subfamily phylogeny, not a numbered panel
python  figures/general_regions_tracks.py    # chromatin state tracks, used across figures
```

[fig1_individual_maps.ipynb](figures/fig1_individual_maps.ipynb) plots contact maps for the
top and bottom disruption deciles. `fig2_1KG_SV.R` reads `1KG_T2T_scored_SV.txt`, produced
by the 1000 Genomes scripts, which are not included here.

## Data

Inputs are not included and must be obtained separately:

- `hg38.fa`, `CHM13.fa` and the Akita model weights, from the SuPreMo/Akita release
- RepeatMasker *Alu* annotations for hg38 and CHM13
- GENCODE gene annotation (v29 and v48)
- ENCODE ChIP-seq, WGBS and CTCF peaks for H1 and HFFc6
- chromHMM state calls and the state universe
- the 1000 Genomes long-read structural variant callset

## Contact

shu.zhang@gladstone.ucsf.edu
