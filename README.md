# Alu Elements as Modulators of 3D Genome Architecture

Alu elements are among the most abundant mobile DNA in the human genome, but how an
individual insertion shapes nearby chromatin interactions is largely unknown. This
repository contains the code for a genome-wide *in silico* deletion screen of ~1.1 million
Alus, scored with a sequence-to-structure deep learning model, together with the
annotation, perturbation and figure code used to characterise the high-scoring subset.

- [Installation](#installation)
- [Instructions](#instructions)
  - [Worked example: scoring 10 Alus](#worked-example-scoring-10-alus)
  - [Genomic context annotation](#genomic-context-annotation)
  - [Targeted in silico perturbation](#targeted-in-silico-perturbation)
  - [Figures](#figures)
- [Data](#data)
- [Contact](#contact)

## Installation

Scoring is done with [SuPreMo](https://github.com/pollard-lab/SuPreMo), which wraps the
[Akita](https://github.com/calico/basenji/tree/master/manuscripts/akita) contact-map model;
a smaller set of analyses uses [Enformer](https://github.com/deepmind/deepmind-research/tree/master/enformer).
Clone SuPreMo/Akita and Basenji separately and obtain the Akita model weights and `hg38.fa`
(and `CHM13.fa` for the T2T and 1000 Genomes runs) before running anything here.

Analyses were run under a conda environment on Python 3.9.0, with `tensorflow` 2.8.0,
`numpy` 1.21.5, `pandas` 1.4.4, `scipy` 1.7.3, `scikit-learn` 1.6.1, `matplotlib` 3.5.1,
`h5py` 3.14.0, `cooltools` 0.5.1, `cooler` 0.8.11, `pysam` 0.19.0, `biopython` 1.80,
`pyBigWig` 0.3.18, `pybedtools`, `pyranges` 0.1.4, `shap` 0.48.0 and `astropy` 5.0.4.
GPU scoring additionally needs a working CUDA runtime for TensorFlow.

Every script carries a short `Paths` block at the top holding the machine-specific roots it
needs. Each root reads an environment variable and falls back to the path it was run from,
so you can either export the variables once or edit the block. Each script defines only the
roots it uses, so no one file contains all of these.

```bash
export ALU_PROJECT_DIR=/path/to/alu                    # inputs, results, figures
export ALU_DATA_DIR=/path/to/data/projects/<user>      # bulk annotation data
export SUPREMO_DIR=/path/to/akita_variant_scoring      # SuPreMo and the Akita model
export AKITA_DIR=/path/to/akita                        # Akita checkout
export POLLARD_DATA=/path/to/data                      # shared mount (UCSC tracks)
export SLURM_DIR=/path/to/slurm                        # job logs
```

The same variables work in the shell scripts and, via `Sys.getenv`, in the R figure
scripts. In code the block looks like this:

```python
PROJECT_DIR = os.environ.get("ALU_PROJECT_DIR", "/pollard/home/szhang20/alu")
SUPREMO_DIR = os.environ.get("SUPREMO_DIR", "/pollard/home/szhang20/akita_variant_scoring")
```

`ALU_DATA_DIR` covers both `DATA_DIR` and `PROJECTS_DIR`, which hold the same root under
two names (the second exists only in `feature_annotation/`, where the scripts already use
`DATA_DIR` for a subdirectory of it).

One wrinkle: in `feature_annotation/` the bulk data root is called `PROJECTS_DIR` rather
than `DATA_DIR`, because those scripts already define their own `DATA_DIR` for a
subdirectory of it. Same value, different name, to keep the two from colliding.

The R scripts in `figures/` carry the same kind of block, with `PROJECT_DIR` pointing at a
local copy of the results and figure directories.

Shared helpers live in [utils/](utils): the contact-map scoring functions, the Akita
plotting utilities, and `io_helpers.py` for the few small functions the analysis scripts
share (FASTA reading, GTF attribute parsing, multi-page PDF writing). Scripts reach them by
appending `utils/` to `sys.path`, so they work from any working directory.

Two caveats. `utils/akita_utils_forplotting.py` builds the Akita model when it is imported,
so any script importing it pays that cost up front. And `akita_utils_scoring` and
`akita_utils_forplotting` are vendored from the `akita_variant_scoring` fork rather than the
public SuPreMo release, which does not contain them; the scripts that use them resolve them
from `SUPREMO_DIR` at run time, so the copies here are for reference.

## Instructions

Scripts that take arguments default to the values used for the manuscript, so each can be
run bare to reproduce what was done; pass `--help` to see the options.

### Worked example: scoring 10 Alus

[example/](example) is a self-contained version of the screen, small enough to run in a
couple of minutes. It scores ten real hg38 Alus spanning the AluY, AluS and AluJ families
as deletions, under all six sequence augmentations: shifts of 0, +1 and -1 bp, each forward
and reverse complemented.

```bash
bash example/run_supremo_example.sh
```

That writes `example/output/sample_10alus_scores` (one column per augmentation) and
`sample_10alus_maps.npy` (the predicted reference and alternate contact maps). A copy of
both is checked in under [example/expected_output/](example/expected_output), so the
notebook runs without the model present.

[combine_scores_and_plot.ipynb](example/combine_scores_and_plot.ipynb) picks the six
augmentation columns up by prefix, collapses them to a mean and standard deviation per Alu,
and plots three things: disruption and correlation with their spread across augmentations,
the reference/alternate/difference maps for the most disrupted Alu, and the difference maps
for all ten on a shared scale. The spread is worth looking at — across these ten Alus the
standard deviation over augmentations is a median 37% of the mean MSE, against a 3.3-fold
range between the Alus themselves.

The full screen applied the same command to ~1.1M Alus, sharded across GPU jobs.

### Genomic context annotation

Annotate the scored Alus with chromatin and sequence features at the Alu itself and in
1 kb, 10 kb, 100 kb and 1 Mb windows around it, then model disruption from those features.

Three scripts cover this, each with one subcommand per feature type; run any of them
with `--help` to list the subcommands.

```bash
# counts of features overlapping each Alu and its windows
python feature_annotation/count_features.py elements --window-sizes 1000 10000 100000
python feature_annotation/count_features.py chromhmm
python feature_annotation/count_features.py ctcf
python feature_annotation/count_features.py wgbs
python feature_annotation/count_features.py bigwig       # phyloP / phastCons / mappability
python feature_annotation/count_features.py bigwig_mp    # same, parallelised for the 1Mb window

# fractional overlap with bed tracks
python feature_annotation/overlap_features.py tracks
python feature_annotation/overlap_features.py chip --encode-dir /path/to/encode/ChIP/ --workers 8

# mean GC content and bigWig signal in the same windows
python feature_annotation/average_features.py gc
python feature_annotation/average_features.py bigwig

# elastic net and logistic regression on the annotated table
python feature_annotation/20250206_alu_feature_elasticnet.py
```

Enformer tracks over the high- and low-disruption sets are produced by
[run_enformer_alu.sh](feature_annotation/run_enformer_alu.sh).

### Targeted in silico perturbation

All deliberate sequence changes live in [in_silico_mutagenesis/](in_silico_mutagenesis).

**Sequence composition.** Mutate GC, AT and CpG content within the Alu body and in flanking
windows, at a series of levels, with a shuffled-sequence control. All four run from one
script: `--mutations` picks the types, `--windows` how much flanking sequence goes with
them, and `--dry-run` lists the jobs without running them. The defaults reproduce the 146
jobs used for the manuscript.

```bash
python in_silico_mutagenesis/mutate_composition.py                    # all 146 jobs
python in_silico_mutagenesis/mutate_composition.py --dry-run          # list, do not run
python in_silico_mutagenesis/mutate_composition.py -m GC --windows 0  # one type, Alu body only

python in_silico_mutagenesis/build_shuffled_controls.py               # shuffled input prep
```

**Rolling-window mutagenesis.** Slide a mutated window across each Alu to localise which
part of the element carries the disruption signal.

```bash
python in_silico_mutagenesis/build_rolling_window_input.py
python in_silico_mutagenesis/build_rolling_window_input_5100.py
python in_silico_mutagenesis/run_supremo_sharded.py \
    -d $RES/20260413_alu_rollingwindow_top100/ -p alu_top100_repName_rolling -n 12
```

**Random deletion controls.** Delete random 300 bp segments at matched distances to
separate Alu identity from the act of deleting sequence. The same sharded runner scores
both spacings.

```bash
python in_silico_mutagenesis/build_random_300bp_input.py
python in_silico_mutagenesis/run_supremo_sharded.py \
    -d $RES/300bp_random/10kb/ -p alu_300bp_10kb_random -n 20
python in_silico_mutagenesis/run_supremo_sharded.py \
    -d $RES/300bp_random/100kb/ -p alu_300bp_100kb_random -n 20
```

**Transcription factor motifs.** Delete TF motifs inside Alus, chiefly ARID3B, and compute
binned saliency across the element.

```bash
python in_silico_mutagenesis/arid3b_motif_deletion.py
python in_silico_mutagenesis/arid3b_binned_saliency.py
```

**Gene swaps.** Move Alus into and out of gene bodies and introns to test how much of the
effect is set by genomic context.

```bash
python in_silico_mutagenesis/gene_swaps_parallel.py
python in_silico_mutagenesis/gene_intron_swaps.py
python in_silico_mutagenesis/gene_overlap_null.py
```

**Consensus insertions and Alu–CTCF pairs.** Insert consensus Alus at chosen positions,
measure insertion gradients, and score pairs of Alus and CTCF motifs together.

```bash
python in_silico_mutagenesis/insert_consensus_alu.py
python in_silico_mutagenesis/insert_consensus_alu_parallel.py
python in_silico_mutagenesis/run_insertion_gradients.py

python in_silico_mutagenesis/run_pairs_all_modes.py
python in_silico_mutagenesis/pairs_alu_ctcf.py
python in_silico_mutagenesis/pairs_ctcf_ctcf.py
python in_silico_mutagenesis/pairs_alu_random.py
```

[blank_canvas_ctcf.ipynb](in_silico_mutagenesis/blank_canvas_ctcf.ipynb) is kept as a
notebook: it generates random sequence at a chosen GC content, parses the CTCF JASPAR
matrix, and inserts motifs in either orientation, and is the worked example for the
synthetic-sequence panels.

### Figures

Scripts in [figures/](figures) are named for the manuscript figure they draw, and read the
tables produced above.

```bash
Rscript figures/fig1_all_alu_scores.R        # genome-wide score distributions
Rscript figures/fig1_karyotype.R             # disruption across chromosomes
Rscript figures/fig2_1KG_SV.R                # polymorphic Alus (see note below)
Rscript figures/fig3_elasticnet.R            # feature model coefficients
Rscript figures/fig4_gene_swaps.R            # gene swaps
Rscript figures/fig5_gc_mutagenesis.R        # GC/AT/CpG mutagenesis
Rscript figures/fig5_random_300bp.R          # random deletion controls
Rscript figures/fig5_rolling_window.R        # rolling-window profiles
Rscript figures/fig6_consensus_insertions.R  # consensus-Alu insertions

Rscript figures/phylo_tree.R                 # Alu subfamily phylogeny, not a numbered panel
python  figures/general_regions_tracks.py    # chromatin state tracks, used across figures
```

[fig1_individual_maps.ipynb](figures/fig1_individual_maps.ipynb) is kept as a notebook: it
selects the top and bottom disruption deciles genome-wide and plots their contact maps, and
is the worked example for the individual-locus panels.

The R scripts carry their own `Paths` block, reading `ALU_PROJECT_DIR` and falling back to
a local copy of the results and figure directories.

`fig2_1KG_SV.R` reads `1KG_T2T_scored_SV.txt`, which was produced by the 1000 Genomes
scoring and annotation scripts; those are not included in this repository.

## Data

Inputs are not included in this repository and must be obtained separately:

- `hg38.fa` and `CHM13.fa`, plus the Akita model weights, from the SuPreMo/Akita release
- RepeatMasker Alu annotations for hg38 and CHM13
- GENCODE gene annotation (v29 and v48 are both used)
- ENCODE ChIP-seq, WGBS and CTCF peak files for H1 and HFFc6
- chromHMM state calls and the chromHMM state universe
- the 1000 Genomes long-read structural variant callset

The paths to each are set in the `Paths` block at the top of the scripts that read them.

## Contact

Questions about the code can be directed to shu.zhang@gladstone.ucsf.edu.
