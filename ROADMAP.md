# SCEPTRE roadmap

Last reviewed: October 2026.

The order matters. Work on a lower priority waits unless it unblocks a higher
one. Topics under "Later" start only once a proposal for them is accepted.

## 1. GPU acceleration of association testing

**Goal.** Make full-screen permutation and conditional randomization test (CRT)
analyses run in minutes on a GPU instead of hours on a CPU cluster.

**Why.** Perturb-seq screens keep growing. Earlier GPU prototypes in pysceptre
ran a full-screen CRT analysis of the Gasperini et al. (2019) dataset in about
ten minutes, against about fifteen hours on 32 CPU cores. That speedup needs to hold on the maintained code
and be backed by a reproducible comparison.

**Done when.** On the current pysceptre release, full-screen permutation and CRT
analyses of a standard benchmark dataset run on a GPU, results agree with the
CPU implementation according to the [parity definition](specs/), and the timing
against a stated CPU baseline is recorded in the benchmark table.

**Where.** pysceptre.

**Status.** In progress.

## 2. Native data-model foundations

**Goal.** Build both packages on their ecosystem's standard data structures:
`SingleCellExperiment` in R, AnnData/MuData in Python.

**Why.** Standard containers let SCEPTRE fit into existing Bioconductor and
scverse workflows without conversion, and give both packages a clean, shared
model of the data.

**Done when.**

- R: the analysis object is built on `SingleCellExperiment`, disk-backed counts
  use DelayedArray/HDF5Array, data can be imported from an existing
  `SingleCellExperiment`, and objects saved by the released version still load.
- Python: an AnnData/MuData entry point, including disk-backed `.h5ad`/Zarr
  data, with explicit choice of the raw-count layer.
- Shared: a written mapping between the R and Python data models, and a
  versioned result-table schema used by both.

**Where.** sceptre (R), pysceptre, and specs in this repository.

**Status.** An R design exists and will land as a series of reviewed pull
requests after the Bioconductor review concludes. The Python entry point is not
started.

## Supporting work

These serve the two priorities above.

- **Bioconductor review.** The R package was
  [submitted](https://github.com/Bioconductor/BiocContributions/issues/176) in
  August 2026. During review, the submitted version receives only bug fixes and
  responses to reviewer comments.
- **R/Python parity.** A written definition of when the two implementations
  agree, plus small shared test fixtures that each package's own tests check.
- **Benchmarks.** One script per benchmark case and a single results table.

## Later

Each needs an accepted proposal before work starts.

- Improved gRNA assignment, including an overhaul of the mixture-based method.
- Per-pair power or minimum-detectable-effect estimates, building on
  [PerturbPlan](https://katsevich-lab.github.io/perturbplan/).
- Hardening the Nextflow pipeline: containers, tests, array jobs and feature
  parity with the packages.
- Additional association tests in R, such as the Wilcoxon test.
- Seurat interoperability.
- Distributed execution from within the R package.
- A shared native core used by both the R and Python packages.

## Reviewing this roadmap

The maintainers review the roadmap at least quarterly, and whenever evidence or
capacity changes. Changes to the ranked priorities are made by pull request to
this file.
