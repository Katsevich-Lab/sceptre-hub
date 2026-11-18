# The sceptre file format: specification

**Format 0.1 (draft), 18 November 2026.**

**Status.** This is a first draft. It was written in about a week, by one
author, and implemented twice by that same author, in R and in Python. It has
been tested on datasets of 500 and 1,000 cells. Nobody else has reviewed it,
and it has not been run on a large dataset. It looks more finished than it is.

**Scope.** This document fixes how a Perturb-seq dataset and sceptre's
analyses of it are stored: which files exist, what they are called, what they
contain and in what form. It covers the data (counts, cells, features, gRNAs)
and the outputs of each stage of an analysis. It does not define what the
numbers mean: what a p-value tests, how a fold change is estimated, or what an
assignment method does are defined by sceptre's method documentation. That is
a choice, and reviewers may want it otherwise. It also does not cover raw
sequencing output, alignment, or how data gets into this format.

This document states the rules. It gives no reasons: for each rule, the
section of [`RATIONALE.md`](RATIONALE.md) with the same letter and number
(S3, A1, T5, ...) explains the choice, its alternatives and the evidence.
[`REVIEW.md`](REVIEW.md) is a one-page guide for reviewers.

## Versions and changes

- **The version** has two parts, `MAJOR.MINOR`, and is stamped in every data
  file (S10) and every analysis (A1).
  - A **minor** increase adds something an older reader can safely ignore: an
    optional column, table, file or metadata entry. Files that do not use the
    addition are unchanged.
  - A **major** increase is any change after which an older reader could
    misread a file: a rule removed or changed, or an addition that changes the
    meaning of what was already there.
- **What a reader does.** A reader MUST either read a file correctly or detect
  that it cannot. It reads any file whose major version it knows, ignoring what
  it does not know, and MUST refuse a file with another major version.
- **Major version 0 is the draft.** Until sceptre's first release that reads
  and writes these files, the rules may change without a major increase and no
  migration is promised. That release freezes the format as 1.0.
- **How the specification changes.** Nothing enters this document until the R
  library and the Python library both implement it and a file in
  `conformance/` exercises it.

| Version | What it added |
|---|---|
| 0.1 | Everything in this document |

## How to read this

**Three kinds of statement.** Every numbered statement is one of these, and
its table says which:

| Kind | What it is | How it is tested |
|---|---|---|
| **file** | True of every conforming file | By the validators, `validate_screen()` and `validate_analysis()`, in R and in Python |
| **behaviour** | What a reader or a writer does | By the round-trip checks in the two notebooks and by the conformance run |
| **guidance** | Advice or a definition. Not binding, and not tested | — |

**Key words.** In file and behaviour rules, **MUST**, **MUST NOT**, **SHOULD**
and **MAY** have their usual meanings. A rule without one of these words
states a fact about the format with the force of MUST. Rule wording is not
final; the sorting into kinds is.

**The Check column** says what tests a rule today:

- **V**: the validators check it, and a file in `conformance/` exercises it;
- **v**: the validators check it, with no conformance file;
- **tests**: a check in the notebooks exercises it, named in brackets;
- **none**: nothing tests it yet.

Of the 104 statements, 72 are file rules, 26 behaviour rules and 6
guidance. The validators check 69 of the file rules, and tests cover 3
more. Tests cover 17 of the behaviour rules, and the validators 2; 7
have no test: S6.4, S8.4, S9.3, S12.1, S12.2, A1.4, A6.2.

**The machine-readable parts** are in [`spec/`](spec/), and both validators and
both prototype libraries read them:

| File | What it defines |
|---|---|
| `spec/screen.json` | The structure of a data file: required groups, number types, encodings, allowed column kinds, required columns and metadata (section 2) |
| `spec/tables.json` | Every table of an analysis: its stage, key, columns and types (rules T0 to T6) |
| `spec/analysis.schema.json` | `analysis.json` (rule A1.1), as a JSON Schema |
| `spec/stage.schema.json` | The six stage records (rule A2.1), as a JSON Schema |

HDF5 has no schema language, so `spec/screen.json` is this project's own
description. The rules that need computing (the order of cells within a
feature, uniqueness, the same cells in every modality) are code in each
language.

**Precedence.** Where a file of `spec/` and the prose of this document
disagree, `spec/` is the rule and the prose is the mistake.

**Two levels.**

- A file is **conforming** when it breaks no file rule. It is what a conforming
  writer produces and what the validators pass.
- A file is **readable** when it is not conforming and a reader MUST still
  accept it. Today this covers one thing: counts stored in another integer or
  floating-point type (S4.4), as files from other tools have them.

The validators report, for a file, which of the three it is: conforming,
readable, or not readable.

## Terms

| Term | Meaning |
|---|---|
| **Cell** | One row of every count matrix; one observation. Cells have identifiers |
| **Modality** | One kind of measurement on the cells, with its own count matrix, cell table and feature table |
| **Feature** | One column of a modality's count matrix: a gene, a gRNA, a peak, a protein. Features have identifiers |
| **gRNA** | A feature of the `grna` modality |
| **Target** | What a gRNA is aimed at, named by `grna_target`. Several gRNAs may share a target. A gRNA whose target is `non-targeting` is a negative control |
| **Response** | A feature of the modality an analysis tests. Which modality that is, is the analysis's choice |
| **Pair** | A target (or, under some strategies, a gRNA) and a response, tested for association |
| **Stage** | One of the six steps of an analysis: `parameters`, `assignment`, `qc`, `calibration`, `power`, `discovery` |
| **Writer** | Software that produces a data file or the files of a stage |
| **Reader** | Software that takes a data file, or a data file and an analysis folder, and uses their contents |
| **Conforming** | Of a file: breaking no file rule. Of a reader or writer: following every behaviour rule |

## References

The format uses these encodings as they are, and is defined against the
versions it was tested with:

| What | Version |
|---|---|
| AnnData on-disk encodings, as found in the files: `dataframe`, `array`, `string-array`, `categorical` (each 0.2.0); `csc_matrix`, `dict` (each 0.1.0) | as written by anndata 0.13.4 (Python) and MuData 1.16.0 (R) |
| MuData on-disk layout, `encoding-version` 0.1.0 | as written by mudata 0.4.1 (Python) and MuData 1.16.0 (R) |
| HDF5 | as written by h5py 3.16.0 and rhdf5 2.56.0 |
| Apache Parquet | as written by pyarrow 26.0.0 and nanoparquet 0.5.2 |
| JSON Schema | draft 7 |

A later version of one of these libraries that writes a file differently does
not change this document.

---

## 1. Layout

```
my_screen.h5mu                        the data file (section 2)

my_analysis/                          one analysis (sections 3 and 4)
├── analysis.json
├── stages/<stage>.json               one per completed stage
├── pairs/discovery.parquet           pairs/positive_control.parquet
├── assignments.parquet               assignment_flagged_cells.parquet (optional)
├── cells.parquet
├── pair_qc/discovery.parquet         pair_qc/positive_control.parquet
├── negative_control_groups.parquet
├── results/<stage>.parquet           results/<stage>_resamples.parquet (optional)
└── fits/<stage>.parquet              (optional)
```

A dataset is one data file. An analysis is one folder. Any number of analyses
may refer to one data file, and an analysis may be stored anywhere.

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
| S2.5 | behaviour | A reader MUST find the modality an analysis tests through `analysis.json` (A1), never by a fixed name. A package MAY choose it for the user when the file has exactly one modality besides `grna`. | tests (the `atac` variant) |
| S2.6 | behaviour | A reader MUST ignore modalities that the analysis does not name. | tests (the `hashing` variant) |
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

| Rule | Kind | | Check |
|---|---|---|---|
| S4.1 | file | `X/data` MUST be stored as 32-bit signed integers. | V |
| S4.2 | file | `X/indices` MUST be stored as 32-bit signed integers. | V |
| S4.3 | file | `X/indptr` MUST be stored as 64-bit signed integers. | V |
| S4.4 | behaviour | A reader MUST accept any integer or floating-point type for these arrays, provided the values are whole numbers. A file that breaks only S4.1 to S4.3 in this way is readable. | tests (the readable conformance files) |

### S5. Cell covariates

| Rule | Kind | | Check |
|---|---|---|---|
| S5.1 | file | A column of the shared cell table MUST NOT be named with a modality's name followed by an underscore. | V |
| S5.2 | guidance | Facts about the cell (batch, donor, anything the user supplies) MUST be in the shared cell table. Measurements of one modality MUST be in that modality's cell table, without a prefix. | — |
| S5.3 | behaviour | The covariate table of a model is the shared columns together with the columns of `grna` and of the response modality, each modality's columns prefixed with the modality's name and an underscore: `rna_n_umis`, `grna_n_umis`. | tests (round trip) |
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
| S7.5 | behaviour | A writer MUST store a text or true/false covariate as a categorical; true/false has the categories `FALSE`, `TRUE`, in that order. | tests (the types file) |
| S7.6 | behaviour | A categorical's stored category order is its level order, and the first category is the baseline of a model. | tests (the types file: levels not alphabetical) |
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
| S10.1 | file | `uns/sceptre` MUST be an AnnData dictionary with the scalar text entries `format` (a version, `MAJOR.MINOR`), `moi` and `id` (not empty). | V |
| S10.2 | file | `moi` MUST be `high` or `low`. | V |
| S10.3 | file | `format` is the data file's own format version, `0.1` for this document. A reader MUST refuse a major version other than the one it knows. | V |
| S10.4 | behaviour | `id` identifies the dataset. A writer assigns a new random `id` when it imports data, and keeps the `id` when it copies a data file. | tests (subset and re-imported files) |

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

## 3. The analysis folder

### A1. `analysis.json`

| Rule | Kind | | Check |
|---|---|---|---|
| A1.1 | file | `analysis.json` MUST exist and conform to `spec/analysis.schema.json`. It has `format` (the analysis's own format version, `0.1` for this document; a reader MUST refuse a major version other than the one it knows), `response_modality` (the name of the modality the analysis tests; never `grna`), and `screen`. | V |
| A1.2 | file | `screen` describes the data file: its `id` (S10.4); `n_cells`, `n_responses`, `n_grnas`; and `cell_ids_sha256`, `response_ids_sha256`, `grna_ids_sha256`. "Responses" are the features of the response modality. A reader MUST refuse a data file for which any of these differs, or which has no modality of that name. | V |
| A1.3 | file | The SHA-256 of a list of identifiers is that of their UTF-8 bytes, in file order, joined by a line feed, with no line feed at the end. (S9.4 keeps this unambiguous.) | V (through A1.2) |
| A1.4 | behaviour | `screen.path` MAY give the data file's path relative to the folder. It is a hint; a reader MUST NOT rely on it. | none |

### A2. Stage records

The stages are `parameters`, `assignment`, `qc`, then `calibration`, `power`
and `discovery`. The last three do not depend on one another.

| Rule | Kind | | Check |
|---|---|---|---|
| A2.1 | file | A file in `stages/` MUST be named `<stage>.json` for one of the six stages, conform to `spec/stage.schema.json`, and have `stage` equal to its file name. | V |
| A2.2 | file | A stage's record MUST NOT exist unless the records of all earlier stages exist. | V |
| A2.3 | behaviour | A stage has run if and only if its record exists. A writer MUST write the record after the stage's tables. | tests (a stage is cleared) |
| A2.4 | file | A record holds the stage's `settings` (what the user chose, with defaults filled in), `provenance` (`package`, `version`, `language`: the implementation's language, which is `R` or `Python` for the two libraries here and may be any other, and `written_at`: UTC, `YYYY-MM-DDTHH:MM:SSZ`), and `files` (A5). It MAY hold `summaries`, which a reader MUST NOT need. | V (through A2.1) |

### A3. Settings and JSON

| Rule | Kind | | Check |
|---|---|---|---|
| A3.1 | file | The settings of each stage, their names and their allowed values are those of `spec/stage.schema.json`. | V (through A2.1) |
| A3.2 | file | A formula names covariates as S5.3 does (`rna_n_umis`). Setting names keep role words (`response_n_umis_range`, `response_fit_method`). | tests (both readers rebuild the design from the formula) |
| A3.3 | file | JSON files are UTF-8. A missing value is `null`. NaN and infinity MUST NOT appear. | V (a record that is not valid JSON fails A1.1 or A2.1) |
| A3.4 | behaviour | A writer MUST write a real number with 17 significant digits, so that the number read is the number written. | tests (the threshold reproduces every call) |
| A3.5 | behaviour | The resampling budget `B3` in force is that of the latest stage record that carries one. | tests (round trip) |

### A4. Rerunning

| Rule | Kind | | Check |
|---|---|---|---|
| A4.1 | file | The folder MUST NOT hold a table of this specification that no stage record lists. | V |
| A4.2 | behaviour | To rerun a stage, a writer MUST first delete the records and files of that stage and of every later stage, records first. | tests (a stage is cleared) |
| A4.3 | behaviour | A reader MUST ignore files it does not know. | V (a valid conformance case) |

### A5. Files of a stage

| Rule | Kind | | Check |
|---|---|---|---|
| A5.1 | file | Every file a record lists in `files` MUST exist and be readable. | V |
| A5.2 | file | The `rows` a record gives for a file MUST equal the file's number of rows (summed over the parts of a folder). | V |
| A5.3 | file | A record MUST NOT list a table that `spec/tables.json` assigns to another stage. | V |
| A5.4 | file | A record MUST list every table that `spec/tables.json` marks as required for its stage, even when the table has no rows. | V |
| A5.5 | guidance | A `files` entry MAY carry a `sha256` of the file. | — |

### A6. Writing

| Rule | Kind | | Check |
|---|---|---|---|
| A6.1 | guidance | A writer SHOULD write each file under a temporary name and rename it into place. | — |
| A6.2 | behaviour | No stage modifies a file that an earlier stage wrote. | none |

---

## 4. The tables

`spec/tables.json` gives, for each table: the stage that writes it, whether it
is required, its key, and its columns with their types. This section states
the rules that apply to all tables and the meaning of each table.

### T0. Every table

| Rule | Kind | | Check |
|---|---|---|---|
| T0.1 | file | A table MUST have every required column of its definition, each of the type given there, with no missing values except where the definition says a column is nullable. | V |
| T0.2 | file | The rows of each file MUST be sorted by the table's key, comparing identifiers as UTF-8 bytes. | V |
| T0.3 | file | A row group MUST NOT have more than 1,000,000 rows. | v |
| T0.4 | file | No two rows of a table may have the same key. | V |
| T0.5 | file | Cells, gRNAs, responses and targets are named by identifier, and every identifier MUST exist in the data file. | V |
| T0.6 | file | A table is one `.parquet` file, or a folder of that name holding part files with the same columns. Each part is sorted; parts need not be sorted across each other. | V (a valid conformance case) |
| T0.7 | file | An identifier column MAY be plain text or dictionary-encoded. | V (a valid conformance case) |
| T0.8 | behaviour | A reader MUST ignore columns it does not know. In the tables whose definition says `extra_columns: kept`, a reader keeps them. | V (a valid conformance case) |
| T0.9 | file | Tables are Parquet, compressed with snappy, zstd or gzip, or not at all. Types are those of section 5. | types: V (through T0.1); compression: none |

`grna_id` belongs to the key of the pair-QC and results tables when it is
present. T3.3 and T5.6 say when that is, so a table's key follows from the
analysis's strategy and the table's stage; `spec/tables.json` records it as
`present_under`.

### T1. `assignments.parquet` and `assignment_flagged_cells.parquet`

| Rule | Kind | | Check |
|---|---|---|---|
| T1.1 | file | `assignments.parquet` has one row per assigned (gRNA, cell): `grna_id`, `cell_id`. | V (through T0) |
| T1.2 | file | `assignment_flagged_cells.parquet` is optional. It lists, by `cell_id`, the cells that the assignment method judges not to carry exactly one gRNA beyond what the assignment table shows. | V (through T0) |
| T1.3 | behaviour | The cells without exactly one gRNA are those with zero or several assignment rows, together with the flagged cells. | tests (low MOI) |

### T2. `cells.parquet`

| Rule | Kind | | Check |
|---|---|---|---|
| T2.1 | file | One row per cell of the data file: `cell_id`, and `in_use`, true for the cells that passed QC. | V |

### T3. `pairs/` and `pair_qc/`

| Rule | Kind | | Check |
|---|---|---|---|
| T3.1 | file | `pairs/discovery.parquet` and `pairs/positive_control.parquet` hold the pairs to test as the user gave them: `grna_target`, `response_id`. Row order carries no meaning. | V (through T0) |
| T3.2 | file | `pair_qc/discovery.parquet` and `pair_qc/positive_control.parquet` hold the verdict of pairwise QC: the key columns and `pass_qc`. | V (through T0) |
| T3.3 | file | The pair-QC tables MUST have a `grna_id` column under the `singleton` and `bonferroni` strategies, where each gRNA is tested on its own, and MUST NOT have one under `union`. | V |

### T4. `negative_control_groups.parquet`

| Rule | Kind | | Check |
|---|---|---|---|
| T4.1 | file | One row per member of a negative-control group: `grna_group` (the group's name) and `grna_id`. | V (through T0) |
| T4.2 | behaviour | A reader MUST take a group's members from this table, never from the group's name. | tests (control groups) |
| T4.3 | file | A group's name appears in the calibration results in the column that names the tested unit: `grna_target` under the `union` strategy, `grna_id` otherwise. | tests (round trip under each strategy) |

### T5. `results/`

| Rule | Kind | | Check |
|---|---|---|---|
| T5.1 | file | A results table has the key columns, `pass_qc`, `p_value`, and at least one of `fold_change` and `log_2_fold_change`. `se_fold_change` and diagnostic columns MAY follow. | V |
| T5.2 | behaviour | A reader that finds no `log_2_fold_change` MUST compute it as `log(fold_change) / log(2)` with natural logarithms, in that form. A writer MUST NOT omit the column unless it equals that expression for every row. | tests (every example row) |
| T5.3 | file | `significant` MUST NOT be stored. The record of the `calibration` and `discovery` stages has `significance_threshold`: the largest p-value called significant, or `null` if none was. A row is significant when its `p_value` is at most the threshold. | V |
| T5.4 | file | Rows are sorted by key (T0.2), not by p-value. | V (through T0.2) |
| T5.5 | file | `results/<stage>_resamples.parquet` is optional: the key columns and `z_null`, a list of the pair's resampled statistics. | V (through T0) |
| T5.6 | file | A results table and its resamples table MUST have a `grna_id` column under the `singleton` strategy and, for the calibration stage, under `bonferroni`. They MUST NOT have one otherwise: not under `union`, and not for the power and discovery stages under `bonferroni`, whose results are per target. | V |

### T6. `fits/`

| Rule | Kind | | Check |
|---|---|---|---|
| T6.1 | file | `fits/<stage>.parquet` is listed by a stage if and only if that stage's record has a `fits` entry, which names the file and gives the fits' `identity`: `screen_id`, `formula`, `cells_in_use_sha256`, `response_fit_method`, `control_group`. | V |
| T6.2 | file | One row per response: `response_id`, `coefficients` (one per term of the design, in the design's order), `theta`. | V |
| T6.3 | behaviour | The fits are a cache. A reader MAY ignore them. A reader MUST NOT reuse a fit unless every part of its `identity` matches what the reader would use itself, and the fitting method is one that the reader's language shares with the language that made the fit (the record's `provenance.language`). | tests (a changed fit record) |

### T7. The design matrix

| Rule | Kind | | Check |
|---|---|---|---|
| T7.1 | file | The formula in the parameters record MUST be additive, each term one of: a covariate name (S5.4), `log(x)`, `log(x + 1)`, or `0` or `1` for the intercept. | V |
| T7.2 | file | The design matrix is not stored. The parameters record has `design.terms` (the names of the matrix's columns, with covariates named as in S5.3) and `design.sha256`. | V (through A2.1) |
| T7.3 | file | `design.sha256` is the SHA-256 of the matrix's values as little-endian 64-bit floats, column by column. | tests (both readers recompute the hash) |
| T7.4 | behaviour | A reader that fits models MUST rebuild the matrix from the formula and the covariates and refuse the analysis if the term names or the hash differ. | tests (both readers) |

---

## 5. Types

| Value | In the data file | In a table |
|---|---|---|
| Counts | 32-bit signed integer (S4) | — |
| Other integers | 32-bit signed integer | 32-bit signed integer |
| Real numbers | 64-bit float | 64-bit float |
| True/false | categorical `FALSE`, `TRUE` (S7.5) | Parquet boolean |
| Text | UTF-8, `string-array` | UTF-8, plain or dictionary-encoded |
| A list of numbers | — | list of 64-bit floats |
| Missing | NaN; categorical code −1 | null |

---

## 6. Conformance

- **Files.** `conformance/` holds files that conform, files that are readable
  and files that each break one rule. `conformance/cases.json` lists them with
  the rules a validator must report for each and the level it must conclude.
  `make_conformance.py` builds them from `example/`.
- **Validators.** `validate_screen(file)` and `validate_analysis(directory,
  screen_file)` exist with the same names in `format/sceptre_format.R` and
  `format/sceptre_format.py`. They read files directly, build no package
  object, and return every problem they find as rows of `rule`, `where`,
  `message` and `readable`. `conformance_level()` turns the rows into
  "conforming", "readable" or "not readable".
- **Agreement.** `check_conformance.R` and `check_conformance.py` run the two
  validators over every case. Both MUST report exactly the expected rules and
  level, and each reader MUST read every readable file to the same counts as
  the conforming file it was made from.
- **What the validators do not check.** They do not rebuild the design matrix
  (T7.4), recompute anything from the counts, or check the meaning of a value,
  only its form.

## 7. What the two libraries support

"Tested" means a check in the notebooks or the conformance set passes in that
language, on the small examples. Nothing here has been run on a large dataset.

| Feature | R | Python |
|---|---|---|
| Write and read a data file | tested | tested |
| A file written by the other language | tested | tested |
| High MOI; low MOI with the "maximum" assignment method and flagged cells | tested | tested |
| The `union`, `singleton` and `bonferroni` strategies | tested | tested |
| Resampled statistics as a list column | tested | tested |
| A table as a folder of part files | tested | tested |
| A response modality with another name; modalities no analysis uses | tested | tested |
| Reuse of saved fits, and refusal when the method differs | tested | tested |
| Continue an analysis the other language stopped | tested | tested |
| Validators, on all conformance cases | tested | tested |
| Read counts stored in another number type (readable files) | tested | tested |
| A covariate name with a dot in a formula | tested | tested through an alias: pysceptre's own formula parser does not take dots |
| gzip-compressed counts | validated, and read in an earlier measurement | validated; read not tested |
| Disk-backed reads | not tested | not tested |
| Writing a stage while it runs (the libraries write a finished analysis) | not implemented | not implemented |
