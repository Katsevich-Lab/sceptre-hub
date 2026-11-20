# SCEPTRE hub

SCEPTRE provides statistically rigorous methods and software for analyzing
single-cell CRISPR screens. It is developed as one coordinated effort across an
R package, a Python package and a Nextflow pipeline. This repository is where
that coordination happens: shared direction and decisions that span more than
one repository.

## Components

| Component | Role | Code and documentation |
| --- | --- | --- |
| sceptre (R) | Bioconductor package for interactive analysis | [Code](https://github.com/Katsevich-Lab/sceptre) · [Documentation](https://katsevich-lab.github.io/sceptre/) |
| pysceptre (Python) | Python implementation, including GPU-accelerated analysis | [Code](https://github.com/broadinstitute/pysceptre) · [Documentation](https://broadinstitute.github.io/pysceptre/) |
| sceptre-pipeline (Nextflow) | Distributed execution on clusters and the cloud | [Code](https://github.com/timothy-barry/sceptre-pipeline) |
| Book | Methods and worked analyses | [Hands-On Single-Cell CRISPR Screen Analysis](https://timothy-barry.github.io/sceptre-book/) |

## What lives here and what lives in each repository

**The rule:** a change that has to happen in more than one language is decided
here. A change that only affects how one language does something lives in that
language's repository.

This repository holds:

- [ROADMAP.md](ROADMAP.md): the development phases and priorities.
- [perturb-seq-data-format/](perturb-seq-data-format/): a draft of an on-disk
  format for Perturb-seq data, the same file for R and Python, meant for any
  software.
- [sceptre-analysis-format/](sceptre-analysis-format/): a draft of the format
  for sceptre's analyses of such a file.
- [site/](site/): the SCEPTRE site at
  <https://katsevich-lab.github.io/sceptre-hub/>, built from the roadmap and
  the two format folders.
- [LICENSE](LICENSE): the MIT licence, which covers everything here.

Each component repository holds its code, tests, API documentation, vignettes,
changelog, CI, releases and implementation issues.

## Status

Phase 1 of the [roadmap](ROADMAP.md) is in progress. The two formats are drafts
(version 0.1) and have not been reviewed.
