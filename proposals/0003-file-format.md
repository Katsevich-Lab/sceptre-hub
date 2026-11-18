# 0003: A file format for Perturb-seq data and analyses

- **Status:** Draft
- **Author:** Eugene Katsevich (drafted with Claude)
- **Date:** 2026-11-18
- **Affects:** all

## Motivation

The R package, the Python package and the pipeline each store data and results
in their own way, and no two of them read the same files. An agreed format has
to come before the code is reorganised around it, because files outlive the
code that wrote them and every implementation is built against the format.

Two requirements drive the design: the same files on disk in R and in Python,
and data kept separate from analysis, because one dataset has several analyses.

## Proposal

The format is specified in [`file-format/`](../file-format/):

- [`SPEC.md`](../file-format/SPEC.md): the rules, with the exact names and types;
- [`RATIONALE.md`](../file-format/RATIONALE.md): the reasons for each choice;
- [`REVIEW.md`](../file-format/REVIEW.md): one page on where review is most needed;
- [`spec/`](../file-format/spec/): the machine-readable definitions.

In short: a dataset is one MuData `.h5mu` file, written once. An analysis is a
folder beside it or elsewhere, with one folder per stage holding Parquet tables
and a JSON record, and a cache of fitted models. Every file carries a version,
`MAJOR.MINOR`; major version 0 is a draft, frozen as 1.0 at the first release
that reads and writes it.

This proposal replaces [`contracts/file-formats.md`](../contracts/file-formats.md),
which is stale. That file is not changed here; it is to be replaced when this
proposal is decided.

## Alternatives considered

`RATIONALE.md` gives the alternatives for each choice. For the format as a
whole: keeping one file that holds data and results together, as the IGVF
pipeline does (section 9 of the rationale compares the two); and stepping
through intermediate formats, which would mean building readers and writers
more than once.

## Validation

- A prototype reader, writer and validator in R and in Python
  (`file-format/format/`), which write the same files and read each other's.
- Conformance files, valid and deliberately broken, on which the two
  validators agree (`file-format/conformance/`).
- One run of the data file on a real screen of 207,324 cells
  (`file-format/scale-test/`).

Not yet done: an analysis at scale, and review by anyone but the author.

## Implementation

Not started. The packages are unchanged; the prototype lives outside them.
Each package would adopt the format in its own phase, with a validator that
reads the hub's `spec/` files.

## Decision

Open.
