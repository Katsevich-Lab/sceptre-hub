---
title: SCEPTRE
toc: false
---

SCEPTRE provides statistically rigorous methods and software for analyzing
single-cell CRISPR screens. It is developed as one coordinated effort across an
R package, a Python package and a Nextflow pipeline.

## Software

::: {.grid}
::: {.g-col-12 .g-col-md-6 .format-card .border .rounded-3 .shadow-sm .p-4}
### sceptre

The R package, on Bioconductor, for interactive analysis.

[Code](https://github.com/Katsevich-Lab/sceptre) · [Documentation](https://katsevich-lab.github.io/sceptre/)
:::

::: {.g-col-12 .g-col-md-6 .format-card .border .rounded-3 .shadow-sm .p-4}
### pysceptre

The Python package, including GPU-accelerated analysis.

[Code](https://github.com/broadinstitute/pysceptre) · [Documentation](https://broadinstitute.github.io/pysceptre/)
:::

::: {.g-col-12 .g-col-md-6 .format-card .border .rounded-3 .shadow-sm .p-4}
### sceptre-pipeline

The Nextflow pipeline, for distributed analysis on clusters and the cloud.

[Code](https://github.com/timothy-barry/sceptre-pipeline)
:::

::: {.g-col-12 .g-col-md-6 .format-card .border .rounded-3 .shadow-sm .p-4}
### Book

*Hands-On Single-Cell CRISPR Screen Analysis*: the methods and worked analyses.

[Read](https://timothy-barry.github.io/sceptre-book/)
:::
:::

## For developers

- [File formats](file-formats/index.md): the Perturb-seq data file format and the sceptre analysis format.
- [Roadmap](roadmap.qmd): the development phases.
