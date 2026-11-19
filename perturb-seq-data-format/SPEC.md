# The Perturb-seq data file format: specification

**Version 0.1 (draft), November 2026.**

**Status.** This is a first draft. It was written in about a week, by one
author, and implemented twice by that same author, in R and in Python. It has
been tested on datasets of 500 and 1,000 cells and on one real screen of
207,324 cells. Nobody else has reviewed it. It looks more finished than it is.

**What it is.** One MuData `.h5mu` file per screen: counts, cells, features
and gRNAs. The format is a set of names and restrictions on top of MuData and
AnnData, not a new container. It is meant for any software that reads or
writes Perturb-seq data, and nothing in it is specific to one tool. One
consumer is the [sceptre analysis format](../sceptre-analysis-format/), which
stores analyses of a data file and names the file by its identifier.

**Scope.** This document fixes how a Perturb-seq dataset is stored: what the
file contains, what its parts are called and in what form. It covers the
counts, the cells, the features and the gRNAs. It holds nothing an analysis
produces (S12). It does not cover raw sequencing output, alignment, or how
data gets into this format.

Two limits of this version:

- **One target per vector.** Several gRNAs delivered together against the same
  element are covered (S8.6). A vector whose gRNAs have different targets, a
  combinatorial screen, is not, and a reader refuses such a file (S8.7).
- **Target coordinates are genomic positions.** A screen that knocks down RNA
  names its transcript in `grna_target` and may leave the coordinates out. The
  optional `perturbation_type` (S10.5) records which kind of screen it was;
  there is no separate entry for the target molecule.

This document states the rules. It gives no reasons: for each rule, the
section of [`RATIONALE.md`](RATIONALE.md) with the same letter and number
(S3, S8, ...) explains the choice, its alternatives and the evidence.
`overview.html` is the introduction, and says what review is asked for.

## Versions and changes

- **The version** has two parts, `MAJOR.MINOR`, and is in every data file as
  `perturb_seq_format` (S10).
  - A **minor** increase adds something an older reader can safely ignore: an
    optional column, modality or metadata entry. Files that do not use the
    addition are unchanged.
  - A **major** increase is any change after which an older reader could
    misread a file: a rule removed or changed, or an addition that changes the
    meaning of what was already there.
- **What a reader does.** A reader MUST either read a file correctly or detect
  that it cannot. It reads any file whose major version it knows, ignoring what
  it does not know, and MUST refuse a file with another major version.
- **Major version 0 is the draft.** Until the format is frozen as 1.0, its
  rules may change without a major increase and no migration is promised. It
  is frozen separately from the sceptre analysis format, which has a version
  of its own.
- **How the specification changes.** Nothing enters this document until the R
  library and the Python library both implement it and a file in
  `conformance/` exercises it.

| Version | What it added |
|---|---|
| 0.1 | Everything in this document |

There is no dataset version. A dataset that is imported again gets a new `dataset_id`
(S10.4), and that is the protection an analysis needs.

## How to read this

**Three kinds of statement.** Every numbered statement is one of these, and
its table says which:

| Kind | What it is | How it is tested |
|---|---|---|
| **file** | True of every conforming file | By the validator, `validate_screen()`, in R and in Python |
| **behaviour** | What a reader or a writer does | By the conformance run and the checks in `verify/` |
| **guidance** | Advice or a definition. Not binding, and not tested | — |

**Key words.** In file and behaviour rules, **MUST**, **MUST NOT**, **SHOULD**
and **MAY** have their usual meanings. A rule without one of these words
states a fact about the format with the force of MUST. Rule wording is not
final; the sorting into kinds is.

**The Check column** says what tests a rule today:

- **V**: the validators check it, and a file in `conformance/` exercises it;
- **v**: the validators check it, with no conformance file;
- **tests**: a check in `verify/` or the conformance run exercises it, named in brackets;
- **none**: nothing tests it yet.

Of the 61 statements, 40 are file rules, 14 behaviour rules and 7
guidance. The validators check all 40 file rules. Tests or the conformance
run cover 9 of the behaviour rules; 5 have no test: S6.4, S8.4, S9.3, S12.1, S12.2.

**The machine-readable part** is [`spec/screen.json`](spec/screen.json): the
structure of a data file as data: required groups, number types, encodings,
allowed column kinds, required columns and metadata. Both validators read it.
HDF5 has no schema language, so this file is this project's own description.
The rules that need computing (the order of cells within a feature,
uniqueness, the same cells in every modality) are code in each language.

**Precedence.** Where `spec/screen.json` and the prose of this document
disagree, `spec/screen.json` is the rule and the prose is the mistake.

**Conforming.** A file conforms when it breaks no file rule, and otherwise it
does not. It is what a conforming writer produces and what the validators
pass. Writers are strict and readers are lenient: a reader accepts some files
that do not conform, where a behaviour rule says so (S4.4).

## Terms

| Term | Meaning |
|---|---|
| **Cell** | One row of every count matrix; one observation. Cells have identifiers |
| **Modality** | One kind of measurement on the cells, with its own count matrix, cell table and feature table |
| **Feature** | One column of a modality's count matrix: a gene, a gRNA, a peak, a protein. Features have identifiers |
| **gRNA** | A feature of the `grna` modality |
| **Target** | What a gRNA is aimed at, named by `grna_target`. Several gRNAs may share a target. A gRNA whose target is `non-targeting` is a negative control |
| **Response** | A feature of the modality an analysis tests. Which modality that is, is the analysis's choice |
| **Writer** | Software that produces a data file |
| **Reader** | Any software that uses the contents of a data file: a package, a script, a person's own code. A data file is a standard `.h5mu` and is opened with the standard libraries; the rules for readers bind whatever software then uses it |
| **Conforming** | Of a file: breaking no file rule. Of a reader or writer: following every behaviour rule |

## References

The format uses these encodings as they are, and is defined against the
versions it was tested with:

| What | Version |
|---|---|
| AnnData on-disk encodings, as found in the files: `dataframe`, `array`, `string-array`, `categorical` (each 0.2.0); `csc_matrix`, `dict` (each 0.1.0) | as written by anndata 0.13.4 (Python) and MuData 1.16.0 (R) |
| MuData on-disk layout, `encoding-version` 0.1.0 | as written by mudata 0.4.1 (Python) and MuData 1.16.0 (R) |
| HDF5 | as written by h5py 3.16.0 and rhdf5 2.56.0 |

A later version of one of these libraries that writes a file differently does
not change this document.

---

## 1. Layout

```
my_screen.h5mu
├── mod/
│   ├── grna/                         the gRNA modality (S2)
│   │   ├── X/                        counts: data, indices, indptr (S3, S4)
│   │   ├── obs/                      this modality's cell table (S5 to S7)
│   │   └── var/                      the gRNA table (S8)
│   └── rna/                          one modality per thing measured (S2)
│       ├── X/
│       ├── obs/
│       └── var/                      its feature table (S7)
├── obs/                              the shared cell table (S5)
└── uns/                              file metadata (S10)
```

A dataset is one data file. The other groups of a MuData file are as MuData
defines them, and this format says nothing about them. Analyses of a data file
are stored apart from it, and any number of them may refer to one file.

---

## 2. The data file

### S1. Container

| Rule | Kind | | Check |
|---|---|---|---|
| S1.1 | file | The data file MUST be an HDF5 file in the MuData layout, with the root attribute `encoding-type` equal to `MuData`. | V |
| S1.2 | guidance | Every element MUST use an AnnData or MuData encoding as those projects specify it. This format defines names and restrictions, not encodings. | — |

### S2. Modalities

| Rule | Kind | | Check |
|---|---|---|---|
| S2.1 | file | The file MUST have a modality named `grna` and at least one other modality. | V |
| S2.2 | file | A modality name MUST start with a letter and contain only letters, digits and underscores. | V |
| S2.3 | file | Every modality, and the shared cell table, MUST hold the same cells in the same order. | V |
| S2.4 | guidance | A modality other than `grna` SHOULD be named for what was measured. `rna` is the name for RNA. | — |
| S2.5 | behaviour | Which modality an analysis tests is recorded with the analysis (in the sceptre analysis format, in `analysis.json`). A reader MUST find it there, never by a fixed name. A package MAY choose it for the user when the file has exactly one modality besides `grna`. | tests (the conformance file with an extra modality, copied by both libraries) |
| S2.6 | behaviour | A reader MUST ignore modalities that the analysis does not name. A tool that reads a data file and writes it again MUST keep them. | V (a conformance file is copied and keeps its modalities) |
| S2.7 | file | Every modality of a file follows the rules of S3 to S7, whether or not an analysis uses it, so that a file is conforming or not on its own and any modality can be a response later. | V (a valid conformance case) |

### S3. Counts

Rules S3 and S4 apply to the count matrix `X` of every modality.

| Rule | Kind | | Check |
|---|---|---|---|
| S3.1 | file | `X` MUST be an AnnData `csc_matrix` of shape (cells, features): compressed by feature. | V |
| S3.2 | file | `indptr` MUST have one entry per feature plus one, start at 0, never decrease, and end at the number of stored values; `data` and `indices` MUST have that many entries. | V |
| S3.3 | file | Within a feature, cell positions MUST be within range and strictly increasing (no repeats). | V |
| S3.4 | file | Every stored value MUST be a positive integer: `X` holds raw counts, and zeros are not stored. | V |

### S4. Number types of the counts

**The index type** of a modality is one signed integer type for both of its
index arrays. It is 32-bit when the number of stored values, the number of
cells and the number of features are each at most 2,147,483,647, and 64-bit
otherwise. This is the type scipy chooses for a sparse matrix of that size, so
the two arrays always have the same width, and a matrix built by scipy already
has it. A file whose index arrays are wider than its size requires, or differ
from each other, does not conform; a reader still reads it (S4.4).

| Rule | Kind | | Check |
|---|---|---|---|
| S4.1 | file | `X/data` MUST be stored as 32-bit signed integers. | V |
| S4.2 | file | `X/indices` MUST be stored in the index type. | V |
| S4.3 | file | `X/indptr` MUST be stored in the index type. | V |
| S4.4 | behaviour | A reader MUST accept any integer or floating-point type for `X/data`, `X/indices` and `X/indptr`, provided the values are whole numbers, and read the same counts as from a conforming file. Such a file still does not conform. | tests (conformance files with other number types and mixed index widths, read with `mudata` and `MuData`) |

### S5. Cell covariates

| Rule | Kind | | Check |
|---|---|---|---|
| S5.1 | file | A column of the shared cell table MUST NOT be named with a modality's name followed by an underscore. | V |
| S5.2 | guidance | Facts about the cell (batch, donor, anything the user supplies) MUST be in the shared cell table. Measurements of one modality MUST be in that modality's cell table, without a prefix. | — |
| S5.3 | behaviour | The covariate table of a model is the shared columns together with the columns of `grna` and of the response modality, each modality's columns prefixed with the modality's name and an underscore: `rna_n_umis`, `grna_n_umis`. | tests in R (the R writer splits the table and the checks join it again: every R copy in the two `verify/` folders). In Python: none (testing it needs a package) |
| S5.4 | file | The name of a column of a cell table MUST start with a letter and contain only letters, digits, underscores and dots. | V |

### S6. Computed columns

| Rule | Kind | | Check |
|---|---|---|---|
| S6.1 | file | Every modality's cell table MUST have `n_umis` (the cell's total count) and `n_nonzero` (its number of features with a nonzero count). A writer computes them when they are absent. | V |
| S6.2 | file | `n_umis`, `n_nonzero` and `p_mito` MUST be numbers. | V |
| S6.3 | guidance | A modality that has mitochondrial genes MAY have `p_mito`: the proportion, from 0 to 1, of the cell's counts from them. | — |
| S6.4 | behaviour | No other value that follows from stored values may be required by a reader. | none |

### S7. Column kinds

"Cell tables" are the shared cell table and each modality's `obs`. "Feature
tables" are each modality's `var`.

| Rule | Kind | | Check |
|---|---|---|---|
| S7.1 | file | Every column of a cell table MUST be one of: 64-bit float (`array`); 32-bit signed integer (`array`); `categorical`. | V |
| S7.2 | file | Every column of a feature table MUST be one of those three, or text in the `string-array` encoding. | V |
| S7.3 | file | A categorical MUST have unique text categories and codes from −1 (missing) to the number of categories minus one. | V |
| S7.4 | file | An integer column MUST NOT hold a missing value. A writer stores an integer column that has missing values as 64-bit floats. | V |
| S7.5 | behaviour | A writer MUST store a text or true/false covariate as a categorical; true/false has the categories `FALSE`, `TRUE`, in that order. | tests (the encoding matrix; the types files) |
| S7.6 | behaviour | A categorical's stored category order is its level order, and the first category is the baseline of a model. | tests (the encoding matrix; the types files: levels not alphabetical) |
| S7.7 | file | Missing values are NaN in a float column and code −1 in a categorical. A text column of a feature table that has missing values is stored as a categorical. | V (through S7.3, S7.4) |
| S7.8 | file | Identifiers and text use `string-array`, never `nullable-string-array`. | V (through S7.1, S7.2) |

### S8. The gRNA table

The gRNA table is the feature table of `grna`: one row per gRNA.

| Rule | Kind | | Check |
|---|---|---|---|
| S8.1 | file | It MUST have a column `grna_target`, text, with no missing or empty values. | V |
| S8.2 | file | It MAY have `target_chr` (text) and `target_start`, `target_end` (numbers): the position of the gRNA's target. | V |
| S8.3 | behaviour | A non-targeting gRNA has `grna_target` equal to `non-targeting`. | tests (control groups) |
| S8.4 | behaviour | A gRNA has one target. A writer MUST refuse a gRNA listed under two targets. | none |
| S8.5 | file | `grna_group` MUST NOT be stored. It follows from `grna_target` and the analysis's integration strategy. | V |
| S8.6 | file | The gRNA table MAY have a text column `vector_id`, naming the vector that delivers each gRNA. If present it has no missing values: every gRNA is on exactly one vector, and a vector with one gRNA is one row. | V |
| S8.7 | file | All gRNAs with the same `vector_id` MUST have the same `grna_target`: a vector has one target. | V |
| S8.8 | file | A `vector_id` MUST NOT equal the identifier of a gRNA on another vector. A vector with one gRNA may be named after it. | V |
| S8.9 | behaviour | A reader MUST refuse a file that breaks S8.7 or S8.8. A writer MUST NOT sum counts over the gRNAs of a vector: counts are stored per gRNA, and combining them is a step of an analysis. | tests (`check_screen()` refuses the conformance files, in both languages) |
| S8.10 | file | Positions are **1-based with both ends included**: `target_start` is the first base of the target and `target_end` the last. Both are whole numbers, `target_start` is at least 1, and `target_end` is not before `target_start`. A single position has the two equal. | V |
| S8.11 | file | A gRNA has all three of `target_chr`, `target_start` and `target_end`, or none of them. A non-targeting gRNA has none. | V |
| S8.12 | guidance | This is the convention of GTF and GFF files, Ensembl, and Bioconductor's `GRanges`. An interval from a BED file (0-based, end excluded) is converted at import by adding 1 to its start. The format has one convention and no entry that says which is in use. Any position column added later follows it. | — |
| S8.13 | guidance | Chromosome names are text, and the format does not choose between `chr1` and `1`. They have to match the naming of whatever positions they are compared with. | — |
| S8.14 | file | A file in which any gRNA has a position MUST have `genome_build` (S10.5). A position cannot be compared with another without knowing the build of both, and the build cannot be worked out from the file. | V |

### S9. Identifiers

| Rule | Kind | | Check |
|---|---|---|---|
| S9.1 | file | Cell identifiers MUST be unique. | V |
| S9.2 | file | Feature identifiers MUST be unique across all modalities of the file. | V |
| S9.3 | behaviour | Identifiers are UTF-8 text. For cells without barcodes, a writer assigns `"1"`, `"2"`, ... and stores them. | none |
| S9.4 | file | An identifier MUST NOT contain a line feed. | V |

### S10. File metadata

| Rule | Kind | | Check |
|---|---|---|---|
| S10.1 | file | `uns` MUST have the scalar text entries `perturb_seq_format` (a version, `MAJOR.MINOR`) and `dataset_id` (not empty). | V |
| S10.2 | file | The entries this section names are the format's. `uns` MAY hold other entries, of other tools; they are not the format's and a reader ignores them. | V (a valid conformance case) |
| S10.3 | file | `perturb_seq_format` is the version of the data file format, `0.1` for this document. A reader MUST refuse a major version other than the one it knows. | V |
| S10.4 | behaviour | `dataset_id` identifies the dataset. A writer assigns a new random `dataset_id` when it imports data, and keeps it when it copies a data file. | tests (`verify/writers` of the sceptre analysis format's folder: a copied data file keeps its `dataset_id`). A new one on import: none (testing it needs a package) |
| S10.5 | file | `uns` MAY have these scalar text entries; `genome_build` is required where the gRNA table has positions (S8.14). None identifies the dataset (`dataset_id` does), and no reader's behaviour depends on them: `title`, a title for people; `source`, a URL, DOI or accession saying where the data came from; `perturbation_type`; `organism`, as a scientific name such as `Homo sapiens`; `genome_build`, such as `GRCh38`, the build of every coordinate in the file. | V |
| S10.6 | file | `perturbation_type`, if present, MUST be one of `CRISPRko` (knockout), `CRISPRi` (interference), `CRISPRa` (activation), `CRISPRbe` (base editing), `CRISPRkd` (knockdown of RNA). | V |
| S10.7 | guidance | A writer SHOULD record `perturbation_type` and `organism` when it knows them, and `genome_build` even in a file without positions. They are known at import, cannot be worked out from the file, and are lost for good if not written down. | — |

### S11. Compression

| Rule | Kind | | Check |
|---|---|---|---|
| S11.1 | file | The count arrays MAY be uncompressed or gzip-compressed (with or without the shuffle filter). No other filter is allowed. | v |

### S12. What the data file does not hold

| Rule | Kind | | Check |
|---|---|---|---|
| S12.1 | behaviour | The data file MUST NOT hold anything an analysis produces: assignments, which cells were kept, gRNA groups, pairs, or results. | none |
| S12.2 | behaviour | A data file is written once. An analysis MUST NOT modify it. | none |

---

## 3. Types

| Value | In the data file |
|---|---|
| Counts | 32-bit signed integer (S4) |
| Index arrays of the counts | the index type: 32-bit signed integer, or 64-bit past the limit (S4) |
| Other integers | 32-bit signed integer |
| Real numbers | 64-bit float |
| True/false | categorical `FALSE`, `TRUE` (S7.5) |
| Text | UTF-8, `string-array` |
| Missing | NaN; categorical code −1 |

---

## 4. Conformance

- **Files.** `conformance/` holds data files that conform and data files that
  each break one rule or a few. `conformance/cases.json` lists them with the
  rules a validator must report for each, and, for the files of S4.4, the
  conforming file a reader must get the same counts from.
  `verify/make_conformance.py` builds them from `example/`.
- **Reading.** A data file is read with the standard libraries:
  `mudata.read_h5mu` in Python, `MuData::readH5MU` in R. The format has no
  reader of its own. `check_screen()` in `format/` is the small check for the
  refusals of S8.9 and S10.3.
- **Validators.** `validate_screen(file)` exists with the same name in
  `format/perturb_seq_data.R` and `format/perturb_seq_data.py`. It reads the
  file directly, builds no package object, needs no other software, and
  returns every problem it finds as rows of `rule`, `where` and `message`. A
  file conforms when there are no rows.
- **No canonical writer.** The format is defined by this document and `spec/`.
  Any code whose files pass the validators and read back correctly in the
  standard libraries is a conforming writer. `write_screen()` and
  `fmt_write_screen()` in `format/` are temporary writers used by the checks.
- **The verification suite.** `verify/run.sh` runs every check of the format.
  It needs no package that analyses the data.
- **Agreement.** `verify/check_conformance.R` and `verify/check_conformance.py` run the two
  validators over every case. Both MUST report exactly the expected rules,
  and the standard reader of each language MUST read every file of S4.4 to
  the same counts as the conforming file it was made from.
- **What the validators do not check.** They do not recompute anything from
  the counts (S6), or check the meaning of a value, only its form.

## 5. What the two libraries support

"Tested" means a check in `verify/` or the conformance set passes in that
language, on the small examples. "At scale" means the one run on 207,324 cells
reported in `RATIONALE.md`, Appendix A.

| Feature | R | Python |
|---|---|---|
| Write and read a data file | tested, and at scale | tested, and at scale |
| A file written by the other language | tested | tested |
| A response modality with another name; modalities no analysis uses | tested | tested |
| Validators, on all conformance cases | tested | tested |
| Read counts and index arrays stored in another number type (S4.4) | tested | tested |
| gzip-compressed counts | validated, and read in a measurement (`RATIONALE.md`, Appendix A) | validated; read not tested |
| Reading counts from disk, one feature at a time | at scale, through MuData's backed reader and through rhdf5 | at scale, through h5py |
| Counts past the 32-bit index limit | not tested at that size (the 64-bit types are tested with the limit lowered) | not tested at that size (the 64-bit types are tested with the limit lowered) |
