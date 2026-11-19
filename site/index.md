---
title: Formats for Perturb-seq data and sceptre analyses
toc: false
---

Two formats, defined separately. The second builds on the first.

| | Perturb-seq data file format | sceptre analysis format |
|---|---|---|
| **What it stores** | The counts of a Perturb-seq dataset, and what is known about its cells, features and gRNAs. Independent of any analysis tool | One analysis of a data file by sceptre: settings, gRNA assignments, quality control and results, stage by stage |
| **Start here** | [Overview](perturb-seq-data-format/overview.qmd) | [Overview](sceptre-analysis-format/overview.qmd) |
| **The rules** | [Specification](perturb-seq-data-format/SPEC.md) | [Specification](sceptre-analysis-format/SPEC.md) |
| **The reasons** | [Rationale](perturb-seq-data-format/RATIONALE.md) | [Rationale](sceptre-analysis-format/RATIONALE.md) |

Every rule has a number, such as S8.3 or T5.2, and a link of its own: a
mention of a rule anywhere on this site leads to it. Rules starting with S
belong to the data file format; rules starting with A or T belong to the
analysis format.

These pages are built from the Markdown files in the
[repository](https://github.com/Katsevich-Lab/sceptre-hub), which are the
source. The repository also holds the validators, the example files and the
test suites.
