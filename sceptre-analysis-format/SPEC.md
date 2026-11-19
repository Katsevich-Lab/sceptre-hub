# The sceptre analysis format: specification

**Version 0.1 (draft), November 2026.**

**Status.** This is a first draft. It was written in about a week, by one
author, and implemented twice by that same author, in R and in Python. It has
been tested on analyses of datasets of 500 and 1,000 cells. Its tables have
been written and read at real size with made-up content, on Linux; no real
analysis has, and nobody else has reviewed it. It looks more finished than it
is.

**What it builds on.** An analysis is an analysis of one data file in the
[Perturb-seq data file format](../perturb-seq-data-format/SPEC.md). This document does not repeat that
format: rules numbered S (S7.6, S10.4, ...) are rules of [its
specification](../perturb-seq-data-format/SPEC.md), and the terms cell, modality, feature, gRNA and target
are defined there. The dependency runs one way. A data file knows nothing
about analyses; an analysis names its data file by identifier (A1).

**Scope.** This document fixes how sceptre's analyses of a Perturb-seq
dataset are stored: which files exist, what they are called, what they
contain and in what form. It covers the outputs of each stage of an analysis.
It does not define what the numbers mean: what a p-value tests, how a fold
change is estimated, or what an assignment method does are defined by
sceptre's method documentation. That is a choice, and reviewers may want it
otherwise.

This document states the rules. It gives no reasons: for each rule, the
section of [`RATIONALE.md`](RATIONALE.md) with the same letter and number
(A1, T5, ...) explains the choice, its alternatives and the evidence.
`overview.html` is the introduction, and says what review is asked for.

## Versions and changes

- **The version** has two parts, `MAJOR.MINOR`, and is in every analysis as
  `format` in `analysis.json` (A1.1). It is separate from the version of the
  data file format, and the two change independently. An analysis records
  which version of the data file format its data file had (A1.2).
  - A **minor** increase adds something an older reader can safely ignore: an
    optional column, table, file or setting. Files that do not use the
    addition are unchanged.
  - A **major** increase is any change after which an older reader could
    misread a file: a rule removed or changed, or an addition that changes the
    meaning of what was already there.
- **What a reader does.** A reader MUST either read a file correctly or detect
  that it cannot. It reads any analysis whose major version it knows, ignoring
  what it does not know, and MUST refuse an analysis with another major
  version.
- **Major version 0 is the draft.** Until the format is frozen as 1.0, its
  rules may change without a major increase and no migration is promised. It
  is frozen at sceptre's first release that reads and writes it, separately
  from the data file format.
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
| **file** | True of every conforming analysis | By the validator, `validate_analysis()`, in R and in Python |
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

Of the 66 statements, 44 are file rules, 18 behaviour rules and 4
guidance. The validators check 42 of the file rules (T0.9 in part), a test covers 1 more
(T4.3), and 1 has no test (A3.2). Tests or the conformance run
cover 11 of the behaviour rules (A1.8 in part, T7.4 in R only); 7 have no
test: A1.5, A1.9, A3.5, A6.2, T1.3, T7.3, T7.5. Seven of the untested rules
(A1.5, A1.9, A3.2, A3.5, T1.3, T7.3, T7.5) can only be tested with a package
that fits models or builds its own object.

**The machine-readable parts** are in [`spec/`](spec/), and both validators and
both prototype libraries read them:

| File | What it defines |
|---|---|
| `spec/tables.json` | Every table of an analysis: its stage, key, columns and types (rules T0 to T6) |
| `spec/analysis.schema.json` | `analysis.json` (rule A1.1), as a JSON Schema |
| `spec/stage.schema.json` | The six stage records (rule A2.1), as a JSON Schema |
| `spec/fits.schema.json` | The identity each part of the fit cache carries: how its fits were made (rule T6.1), as a JSON Schema |

**Precedence.** Where a file of `spec/` and the prose of this document
disagree, `spec/` is the rule and the prose is the mistake.

**Conforming.** An analysis conforms when it breaks no file rule, and
otherwise it does not. It is what a conforming writer produces and what the
validators pass.

**An analysis and its data file.** Checking an analysis in full needs its
data file, because some rules are about the two together: that the link
matches the file (A1.2), that the response modality is in it (A1.1), and that
every cell, gRNA, response and target named in a table exists in it (T0.5,
T2.1). Everything else can be checked from the folder alone: the records and
their schemas, which stages have run, the files each stage lists, and each
table's columns, types, keys and order. `validate_analysis()` takes the data
file as an optional second argument and checks what it can with what it is
given.

## Terms

Cell, modality, feature, gRNA and target are defined in the [data file
format's specification](../perturb-seq-data-format/SPEC.md#terms).

| Term | Meaning |
|---|---|
| **Response** | A feature of the modality an analysis tests. Which modality that is, is the analysis's choice |
| **Pair** | A target (or, under some strategies, a gRNA) and a response, tested for association |
| **Stage** | One of the six steps of an analysis: `parameters`, `assignment`, `qc`, `calibration`, `power`, `discovery` |
| **Writer** | Software that produces the files of a stage |
| **Reader** | Any software that uses the contents of an analysis folder: a package, a script, a person's own code |
| **Conforming** | Of an analysis: breaking no file rule. Of a reader or writer: following every behaviour rule |

## References

The format uses these encodings as they are, and is defined against the
versions it was tested with:

| What | Version |
|---|---|
| The Perturb-seq data file format | 0.1 |
| Apache Parquet | as written by pyarrow 26.0.0 and nanoparquet 0.5.2 |
| JSON Schema | draft 7 |

A later version of one of these libraries that writes a file differently does
not change this document.

---

## 1. Layout

```
my_analysis/                          one analysis
├── analysis.json
├── parameters/                       one folder per completed stage:
│   ├── record.json                   the stage's record, written last
│   ├── discovery_pairs.parquet
│   └── positive_control_pairs.parquet
├── assignment/
│   ├── record.json
│   ├── assignments.parquet
│   └── flagged_cells.parquet         (optional)
├── qc/
│   ├── record.json
│   ├── cells.parquet
│   ├── discovery_pair_qc.parquet
│   └── positive_control_pair_qc.parquet
├── calibration/
│   ├── record.json
│   ├── negative_control_groups.parquet
│   ├── results.parquet
│   └── resamples.parquet             (optional)
├── power/
│   ├── record.json
│   ├── results.parquet
│   └── resamples.parquet             (optional)
├── discovery/
│   ├── record.json
│   ├── results.parquet
│   └── resamples.parquet             (optional)
└── fits/                             NOT a stage: a cache, with rules of its own (T6)
    └── response/
        └── part-….parquet            one part file per job that added fits;
                                      each says inside how its fits were made
```

An analysis is one folder. Any number of analyses may refer to one data file,
and an analysis may be stored anywhere.

---

## 2. The analysis folder

### A1. `analysis.json`

| Rule | Kind | | Check |
|---|---|---|---|
| A1.1 | file | `analysis.json` MUST exist and conform to `spec/analysis.schema.json`. It has `format` (the analysis's own format version, `0.1` for this document; a reader MUST refuse a major version other than the one it knows), `response_modality` (the name of the modality the analysis tests; never `grna`), and `screen`. | V |
| A1.2 | file | `screen` describes the data file: `id`, which is the data file's `dataset_id` (S10.4); `format`, which is the data file's `perturb_seq_format` (S10.3) when the analysis was created; `n_cells`, `n_responses`, `n_grnas`; and `cell_ids_sha256`, `response_ids_sha256`, `grna_ids_sha256`. "Responses" are the features of the response modality. A reader MUST refuse a data file for which `id`, a count or a hash differs, or which has no modality of that name. `format` is recorded for people and tools and is not part of the identity a reader checks. | V |
| A1.3 | file | The SHA-256 of a list of identifiers is that of their UTF-8 bytes, in file order, joined by a line feed, with no line feed at the end. (S9.4 keeps this unambiguous.) | V (through A1.2) |
| A1.4 | file | `screen.path` is optional. If present it MUST be a relative POSIX path or an `http` or `https` URL. It MUST NOT be an absolute local path. | V |
| A1.5 | behaviour | `screen.path` is a hint. A reader takes the data file's location from its caller and MUST NOT depend on the hint. A writer records a relative path only when the data file is in the analysis folder's parent or grandparent folder or below one of them, and otherwise records nothing. | none (testing it needs a package) |
| A1.6 | file | `screen.title` is optional: the data file's `title` (S10.5), copied so that a person can tell which dataset an analysis belongs to. It is not part of the identity a reader checks. | V (through A1.1) |
| A1.7 | guidance | Where a dataset lives is for each package and each user's own configuration to resolve, for instance by a search path or a table from dataset `id` to location. Nothing specific to one machine belongs in `analysis.json`. A package that follows the hint by itself SHOULD do so only for locations the user has configured, since a folder from elsewhere could point a reader at any file or URL. | — |
| A1.8 | behaviour | Which modality an analysis tests is recorded in `analysis.json` (A1.1); the data file does not mark one. A reader MUST find it there, never by a fixed name. A package MAY choose it for the user when the data file has exactly one modality besides `grna`. | tests (the validators take the name from `analysis.json`: the conformance folder that names a modality its data file lacks). The package's choice: none |
| A1.9 | behaviour | A reader MUST ignore every modality of the data file other than `grna` and the one the analysis names. | none (testing it needs a package) |

### A2. Stage records

The stages are `parameters`, `assignment`, `qc`, then `calibration`, `power`
and `discovery`. The last three do not depend on one another.

| Rule | Kind | | Check |
|---|---|---|---|
| A2.1 | file | A stage's record is the file `record.json` in the folder named for the stage. It MUST conform to `spec/stage.schema.json` and have `stage` equal to the folder's name. | V |
| A2.2 | file | A stage's record MUST NOT exist unless the records of all earlier stages exist. | V |
| A2.3 | behaviour | A stage has run if and only if its record exists. A writer MUST write the record after the stage's tables, all in the stage's folder. A writer that finds a stage's folder without a record MUST delete the folder before writing the stage. | tests (`verify/writers`: a stage is written again over an unfinished attempt) |
| A2.4 | file | A record holds the stage's `settings` (what the user chose, with defaults filled in), `provenance` (`package`, `version`, `language`: the implementation's language, which is `R` or `Python` for the two libraries here and may be any other, and `written_at`: UTC, `YYYY-MM-DDTHH:MM:SSZ`), and `files` (A5). It MAY hold `summaries`, which a reader MUST NOT need. | V (through A2.1) |

### A3. Settings and JSON

| Rule | Kind | | Check |
|---|---|---|---|
| A3.1 | file | The settings of each stage, their names and their allowed values are those of `spec/stage.schema.json`. | V (through A2.1) |
| A3.2 | file | A formula names covariates as T7.4 does (`rna_n_umis`). Setting names keep role words (`response_n_umis_range`, `response_fit_method`). | none (testing it needs a package) |
| A3.3 | file | JSON files are UTF-8. A missing value is `null`. NaN and infinity MUST NOT appear. | V (a record that is not valid JSON fails A1.1 or A2.1) |
| A3.4 | behaviour | A writer MUST write a real number with 17 significant digits, so that the number read is the number written. | tests (the threshold reproduces every call) |
| A3.5 | behaviour | The resampling budget `B3` in force is that of the latest stage record that carries one. | none (testing it needs a package) |
| A3.6 | file | The parameters record's `moi` setting is `high` or `low`: the preset from which a package chose its defaults, such as the control group and what is done with cells carrying several gRNAs. It is a choice of the analysis, and the data file does not record it. | V (through A2.1) |

### A4. Rerunning

| Rule | Kind | | Check |
|---|---|---|---|
| A4.1 | file | A stage's folder MUST NOT hold one of that stage's tables that its record does not list. A folder holding tables and no record is a stage that did not finish. | V |
| A4.2 | behaviour | To rerun a stage, a writer MUST first delete the folders of that stage and of every later stage, each record before the rest of its folder. If the `qc` stage or an earlier one is among them, the writer SHOULD delete the fit cache too (T6.6). | tests (`verify/writers`: a retried stage) |
| A4.3 | behaviour | A reader MUST ignore files it does not know. | V (a valid conformance case) |

### A5. Files of a stage

| Rule | Kind | | Check |
|---|---|---|---|
| A5.1 | file | Every file a record lists in `files` MUST exist in the record's folder and be readable. | V |
| A5.2 | file | The `rows` a record gives for a file MUST equal the file's number of rows (summed over the parts of a folder). | V |
| A5.3 | file | A `files` path names a table in the record's own folder: a file name ending in `.parquet`, with no folder part. A stage cannot list a file outside its folder, and a reader MUST refuse a record that tries to, or whose `stage` is not its folder's name. | V (through A2.1; `read_analysis()` refuses the conformance files, in both languages) |
| A5.4 | file | A record MUST list every table that `spec/tables.json` marks as required for its stage, even when the table has no rows. | V |
| A5.5 | guidance | A `files` entry MAY carry a `sha256` of the file. | — |

### A6. Writing

| Rule | Kind | | Check |
|---|---|---|---|
| A6.1 | guidance | A writer SHOULD write each file under a temporary name and rename it into place. | — |
| A6.2 | behaviour | No stage modifies a file that an earlier stage wrote. | none |

---

## 3. The tables

`spec/tables.json` gives, for each table: its path, the stage that writes it
(which is the folder it is in), whether it is required, its key, and its
columns with their types. This section states
the rules that apply to all tables and the meaning of each table.

### T0. Every table

| Rule | Kind | | Check |
|---|---|---|---|
| T0.1 | file | A table MUST have every required column of its definition, each of the type given there, with no missing values except where the definition says a column is nullable. | V |
| T0.2 | file | The rows of each file MUST be sorted by the table's key, comparing identifiers as UTF-8 bytes. | V |
| T0.3 | file | A row group MUST NOT have more than 1,000,000 rows. | v |
| T0.4 | file | No two rows of a table may have the same key. The fit cache is the one exception: two of its part files may hold the same response (T6.4), though one file may not. | V |
| T0.5 | file | Cells, gRNAs, responses and targets are named by identifier, and every identifier MUST exist in the data file. | V |
| T0.6 | file | A table is one `.parquet` file, or a folder of that name holding part files with the same columns. Each part is sorted; parts need not be sorted across each other. | V (a valid conformance case) |
| T0.7 | file | An identifier column MAY be plain text or dictionary-encoded. | V (a valid conformance case) |
| T0.8 | behaviour | A reader MUST ignore columns it does not know. In the tables whose definition says `extra_columns: kept`, a reader keeps them. | V (a valid conformance case) |
| T0.9 | file | Tables are Parquet, compressed with snappy, zstd or gzip, or not at all. Types are those of section 4. | types: V (through T0.1); compression: none |

`grna_id` belongs to the key of the pair-QC and results tables when it is
present. T3.3 and T5.6 say when that is, so a table's key follows from the
analysis's strategy and the table's stage; `spec/tables.json` records it as
`present_under`.

### T1. `assignment/assignments.parquet` and `assignment/flagged_cells.parquet`

| Rule | Kind | | Check |
|---|---|---|---|
| T1.1 | file | `assignment/assignments.parquet` has one row per assigned (gRNA, cell): `grna_id`, `cell_id`. | V (through T0) |
| T1.2 | file | `assignment/flagged_cells.parquet` is optional. It lists, by `cell_id`, the cells that the assignment method judges not to carry exactly one gRNA beyond what the assignment table shows. | V (through T0) |
| T1.3 | behaviour | The cells without exactly one gRNA are those with zero or several assignment rows, together with the flagged cells. | none (testing it needs a package) |

### T2. `qc/cells.parquet`

| Rule | Kind | | Check |
|---|---|---|---|
| T2.1 | file | One row per cell of the data file: `cell_id`, and `in_use`, true for the cells that passed QC. | V |

### T3. The pairs to test, and their QC

| Rule | Kind | | Check |
|---|---|---|---|
| T3.1 | file | `parameters/discovery_pairs.parquet` and `parameters/positive_control_pairs.parquet` hold the pairs to test as the user gave them: `grna_target`, `response_id`. Row order carries no meaning. | V (through T0) |
| T3.2 | file | `qc/discovery_pair_qc.parquet` and `qc/positive_control_pair_qc.parquet` hold the verdict of pairwise QC: the key columns and `pass_qc`. | V (through T0) |
| T3.3 | file | The pair-QC tables MUST have a `grna_id` column under the `singleton` and `bonferroni` strategies, where each gRNA is tested on its own, and MUST NOT have one under `union`. | V |

### T4. `calibration/negative_control_groups.parquet`

| Rule | Kind | | Check |
|---|---|---|---|
| T4.1 | file | One row per member of a negative-control group: `grna_group` (the group's name) and `grna_id`. | V (through T0) |
| T4.2 | behaviour | A reader MUST take a group's members from this table, never from the group's name. | tests (control groups) |
| T4.3 | file | A group's name appears in the calibration results in the column that names the tested unit: `grna_target` under the `union` strategy, `grna_id` otherwise. | tests (derived values, under each strategy) |

### T5. Results: `<stage>/results.parquet`

| Rule | Kind | | Check |
|---|---|---|---|
| T5.1 | file | A results table has the key columns, `pass_qc`, `p_value` and `log_2_fold_change`. `fold_change`, `se_fold_change` and diagnostic columns MAY follow. | V (through T0.1) |
| T5.2 | file | In a table that has `fold_change`, a row's `log_2_fold_change` MUST equal `log(fold_change) / log(2)` to a relative difference of at most 10^-9, wherever neither is missing. A reader takes the stored value and does not recompute it. | V |
| T5.3 | file | `significant` MUST NOT be stored. The record of the `calibration` and `discovery` stages has `significance_threshold`: the largest p-value called significant, or `null` if none was. A row is significant when its `p_value` is at most the threshold. | V |
| T5.4 | file | Rows are sorted by key (T0.2), not by p-value. | V (through T0.2) |
| T5.5 | file | `<stage>/resamples.parquet` is optional: the key columns and `z_null`, a list of the pair's resampled statistics. | V (through T0) |
| T5.6 | file | A results table and its resamples table MUST have a `grna_id` column under the `singleton` strategy and, for the calibration stage, under `bonferroni`. They MUST NOT have one otherwise: not under `union`, and not for the power and discovery stages under `bonferroni`, whose results are per target. | V |

### T6. The fit cache: `fits/response/`

**The fit cache follows different rules from the stage folders.** It holds
fitted models for the responses, which the three testing stages share. No
stage owns it: it has no stage record, no stage lists its files, it fills as
jobs add to it, and it is cleared at a different moment from any testing
stage. Because it never finishes, it cannot have a record written last, so
each part describes itself. A fit depends on the data, the formula, the fitting method, the control
group and the cells that passed QC, all of which are fixed once the parameters
and QC stages have run.

| Rule | Kind | | Check |
|---|---|---|---|
| T6.1 | file | `fits/response/`, if it exists, MUST hold part files (`*.parquet`) and nothing else. Each part MUST carry the **identity** of the fits it holds: JSON text under the key `sceptre_fit_identity` in the Parquet file's own key-value metadata, conforming to `spec/fits.schema.json`: `screen_id`, `formula`, `cells_in_use_sha256` (the hash of A1.3 over the identifiers of the cells in use), `response_fit_method`, `control_group`. Other keys in that metadata are ignored. | V |
| T6.2 | file | The fits are the rows of the part files: one row per response, with `response_id`, `coefficients` and `theta`. In a part whose identity is that of the analysis, `coefficients` has one value per term of the design, in the design's order. | V |
| T6.3 | file | A fit cache MUST NOT exist before the parameters and QC stages have run. A part whose identity is not that of the analysis (another formula, method, control group, dataset or set of cells) is allowed: it is unused, and the validators do not report it. | V |
| T6.4 | behaviour | Each writer adds a part file of its own, carrying the identity of the fits in it, and MUST NOT modify another's. No file of the cache is shared between writers. Two parts MAY hold the same response; a reader takes either. | V (a valid conformance case with two parts); tests |
| T6.5 | behaviour | The fits are a cache. A reader MAY ignore them, and deleting them changes no result beyond numerical tolerance. A reader MUST NOT take a fit from a part unless every field of that part's identity equals what the reader would use itself and the fitting method is one its language shares. It MUST ignore any other part, and a part with no identity. Identities are compared as parsed values, not as text, and two formulas are equal when they are the same once white space is removed. | tests (parts with another method, another set of cells, no identity, and a formula spaced differently) |
| T6.6 | behaviour | A writer SHOULD delete the cache when the `qc` stage or an earlier stage is rerun, and MUST keep it when a testing stage is rerun. Deleting is housekeeping: no result depends on it (T6.5), and it keeps a plain read of the folder from mixing old parts with current ones. | V (a valid conformance case); tests |
| T6.7 | guidance | `fits/grna/` is reserved for fitted models of the gRNAs, which this version does not define. It would work the same way, with identity fields of its own. A table of another kind of fit is a folder of its own, never a column in this one. | — |

### T7. The design matrix

| Rule | Kind | | Check |
|---|---|---|---|
| T7.1 | file | The formula in the parameters record MUST be additive, each term one of: a covariate name (T7.4), `log(x)`, `log(x + 1)`, or `0` or `1` for the intercept. | V |
| T7.2 | file | The design matrix is not stored. The parameters record has `design.terms` (the names of the matrix's columns, with covariates named as in T7.4). | V (through A2.1) |
| T7.3 | behaviour | A reader that fits models MUST rebuild the matrix from the formula and the covariates and refuse the analysis if the term names differ. | none (testing it needs a package) |
| T7.4 | behaviour | The covariate table of a model is the columns of the data file's shared cell table together with the columns of the cell tables of `grna` and of the response modality, each modality's columns prefixed with the modality's name and an underscore: `rna_n_umis`, `grna_n_umis`. The data file's column names (S5.1, S5.3) keep these names distinct and usable in a formula. | tests in R (the data file format's R writer splits such a table and the checks join it again: every R copy in `verify/` here and in that format's `verify/`). In Python: none (testing it needs a package) |
| T7.5 | behaviour | The baseline of a categorical covariate in the design matrix is its first stored category (S7.6). | none (testing it needs a package) |

---

## 4. Types

| Value | In a table | In a JSON record |
|---|---|---|
| Integers | 32-bit signed integer | a JSON number |
| Real numbers | 64-bit float | a JSON number, written with 17 significant digits (A3.4) |
| True/false | Parquet boolean | `true`, `false` |
| Text | UTF-8, plain or dictionary-encoded | UTF-8 (A3.3) |
| A list of numbers | list of 64-bit floats | a JSON array |
| Missing | null | `null`; NaN and infinity do not appear (A3.3) |

---

## 5. Conformance

- **Files.** `conformance/` holds analysis folders that conform and analysis
  folders that each break one rule or a few. `conformance/cases.json` lists
  them with the rules a validator must report for each, and the data file each
  is an analysis of, which is under `../perturb-seq-data-format/`.
  `verify/make_conformance.py` builds them from `example/`.
- **Validators.** `validate_analysis(directory, screen_file)` exists with the
  same name in `format/sceptre_analysis.R` and `format/sceptre_analysis.py`,
  which load the data file format's code. It reads files directly, builds no
  package object, and returns every problem it finds as rows of `rule`,
  `where` and `message`. An analysis conforms when there are no rows. Whether
  the data file itself conforms is the data file format's validator's
  question.
- **No canonical writer.** The format is defined by this document and `spec/`.
  Any code whose files pass the validators and read back correctly in the
  standard libraries is a conforming writer. The writers in `format/` are
  temporary.
- **The verification suite.** `verify/run.sh` runs every check of the format.
  It needs neither sceptre nor pysceptre. `verify/run.sh format` runs the
  checks that use no writer, and `verify/run.sh writers` the checks of the
  two writers in `format/`.
- **Agreement.** `verify/check_conformance.R` and `verify/check_conformance.py` run the two
  validators over every case. Both MUST report exactly the expected rules.
- **What the validators do not check.** They do not rebuild the design matrix
  (T7.3), recompute anything from the counts, or check the meaning of a value,
  only its form.

## 6. What the two libraries support

"Tested" means a check in `verify/` or the conformance set passes in that
language, on the small examples, or that the one check with both packages
passed it (`RATIONALE.md`, Appendix B).

| Feature | R | Python |
|---|---|---|
| Write and read an analysis | tested | tested |
| An analysis written by the other language | tested | tested |
| High MOI; low MOI with the "maximum" assignment method and flagged cells | tested | tested |
| The `union`, `singleton` and `bonferroni` strategies | tested | tested |
| Resampled statistics as a list column | tested | tested |
| A table as a folder of part files | tested | tested |
| An analysis of a response modality with another name | tested | tested |
| Reuse of saved fits, and refusal when the method differs | tested | tested |
| Continue an analysis the other language stopped | tested | tested |
| Validators, on all conformance cases | tested | tested |
| A covariate name with a dot in a formula | tested | tested through an alias: pysceptre's own formula parser does not take dots |
| The tables of an analysis at real size | measured with made-up content (`scale-test/`); no real analysis | measured with made-up content (`scale-test/`); no real analysis |
| Writing a stage while it runs (the libraries write a finished analysis) | not implemented | not implemented |
