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

## File formats

Two formats, defined separately. The second builds on the first.

::: {.grid}
::: {.g-col-12 .g-col-md-6 .format-card .border .rounded-3 .shadow-sm .p-4}
### Perturb-seq data file format

[![](https://img.shields.io/github/actions/workflow/status/Katsevich-Lab/sceptre-hub/perturb-seq-data-format.yml?branch=main&label={{version perturb-seq-data-format}}){fig-alt="Version and test status"}](https://github.com/Katsevich-Lab/sceptre-hub/actions/workflows/perturb-seq-data-format.yml)

The counts of a Perturb-seq dataset, and what is known about its cells, features and gRNAs. Independent of any analysis tool.

[Overview](perturb-seq-data-format/overview.qmd) · [Specification](perturb-seq-data-format/SPEC.md) · [Rationale](perturb-seq-data-format/RATIONALE.md)
:::

::: {.g-col-12 .g-col-md-6 .format-card .border .rounded-3 .shadow-sm .p-4}
### sceptre analysis format

[![](https://img.shields.io/github/actions/workflow/status/Katsevich-Lab/sceptre-hub/sceptre-analysis-format.yml?branch=main&label={{version sceptre-analysis-format}}){fig-alt="Version and test status"}](https://github.com/Katsevich-Lab/sceptre-hub/actions/workflows/sceptre-analysis-format.yml)

One analysis of a data file by sceptre: settings, gRNA assignments, quality control and results, stage by stage.

[Overview](sceptre-analysis-format/overview.qmd) · [Specification](sceptre-analysis-format/SPEC.md) · [Rationale](sceptre-analysis-format/RATIONALE.md)
:::
:::

## Roadmap

Development proceeds in four phases. The first, agreeing on the file formats,
is in progress. [The roadmap](roadmap.qmd) gives the phases and priorities.
