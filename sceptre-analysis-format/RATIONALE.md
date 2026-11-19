# The sceptre analysis format: rationale

[`SPEC.md`](SPEC.md) states the rules of the sceptre analysis format, numbered,
with no reasons. This document gives the reasons: for each choice, the
alternatives, the evidence, what it locks in, and how it compares with the
IGVF pipeline's format. Sections here and there share their letters and
numbers (A1, T5, ...), so rule `T5.2` of the specification is explained in
section T5 below.

An analysis is an analysis of one data file in the Perturb-seq data file
format. [That format's rationale](../perturb-seq-data-format/RATIONALE.md) holds what the two share, and this
document refers to it instead of repeating it: the
sections S1 to S12 that the S numbers here cite, how the evidence is tagged
and where the IGVF and code facts come from, and [how a format is defined
and checked](../perturb-seq-data-format/RATIONALE.md#7-how-the-format-is-defined-and-checked).

**Status:** draft (format 0.1), November 2026. Every choice below is
implemented in the prototype beside this file. None has been reviewed;
`overview.html` is the introduction, and says what review is asked for. Until
sceptre's first release that reads and writes these files, the format may
change and no migration code is promised. At that release it is frozen as
format 1.0. [Section 7](#7-fixed-now-free-to-grow) says which parts would be
expensive to change after that and which can grow.

Each choice is labelled **Decided** (the maintainer has ruled on it) or
**Adopted** (the author's recommendation, taken up so that the draft is
complete, and more open to challenge).

**How to read this.** Section 1 gives the whole format in a page. Sections 2
to 4 defend each choice in a fixed form: the decision, the alternatives, why,
and what it locks into files. Where this document and `SPEC.md` disagree,
`SPEC.md` is the rule and the disagreement is a mistake to fix. Evidence is
tagged *(verified)*, *(measured)*, *(inferred)* or *(not independently
checked)*, as the data file format's rationale defines.

---

## The principles

Every choice below follows from these. Where a choice needed a trade-off, the
earlier principle wins. Several are the [data file format's
principles](../perturb-seq-data-format/RATIONALE.md#2-the-eight-principles), applied to an analysis folder.

| # | Principle | Why | What follows |
|---|---|---|---|
| 1 | **The same files from R and Python.** | An analysis is started in one language and continued or read in the other, and what a file means should not depend on which wrote it. | Tables are Parquet, with the few column types both languages read and write correctly (T0.1). A real number in a record is written with 17 significant digits, so the number read is the number written (A3.4). Nothing depends on arithmetic giving the same bits on two platforms: a stored value is read, never rebuilt (`log_2_fold_change`, T5), and the design matrix has no hash (T7). Saved fits cross languages only where the fitting method is the same (T6.5). |
| 2 | **An analysis is kept apart from its data.** | One dataset is analysed many times, with different settings. The data is large and fixed at import; an analysis is small and gets redone. | An analysis is a folder of its own. It names its data file by identifier and by hashes of the identifier lists (A1), never modifies it, can be stored anywhere, and can be read for its results without it. Every analysis choice lives here: which modality is the response (A1), the low or high MOI preset (A3.6), the assignments (T1), the cells in use (T2) and the grouping of gRNAs (T4). |
| 3 | **Every file is written once.** | A file edited in place can be left half-written, cannot be written safely by parallel jobs, and no longer shows which step produced what. | No stage edits a file an earlier stage wrote: stages add files (A6.2). Rerunning a stage deletes its files and those of later stages (A4.2). Parallel jobs write part files. A stage's record is written last, so an unfinished stage can be told from a finished one (A2). `significant` is not stored, because it depends on the whole table (T5.3). |
| 4 | **Identifiers, never positions.** | A row number means something only for one file in one order, and R counts from 1 where Python counts from 0. | Tables name cells, gRNAs and responses by identifier, and every identifier exists in the data file (T0.5). Tables survive subsetting, reordering and being written in pieces. |
| 5 | **Existing conventions before invented ones.** | An existing format comes with readers, tools and users. An invented encoding needs code in two languages and has no other reader. | Parquet and JSON are used as specified, and the records are defined by JSON Schemas. How a part of the fit cache was made is stored in Parquet's own key-value metadata (T6.1). The data is a file of the Perturb-seq data file format, not something redefined here. |
| 6 | **Readers ignore what they do not know.** | The format will grow, and files outlive the code that wrote them. Older software has to keep working on newer files where it can, and fail cleanly where it cannot. | Columns and files the format does not define are ignored, or kept where a table says so (T0.8). They can be added with a minor version. Software refuses an analysis of a major version it does not know (A1.1). |
| 7 | **What follows from other stored values is not stored.** | Two copies of one fact can disagree, and the file then does not say which is right. | There is no `significant` column (T5.3) and no count of gRNAs per cell (T2). The exception is a value that cannot be rebuilt to the same number on every platform: `log_2_fold_change` is stored beside `fold_change`, and the validators check that the two agree (T5.2). |
| 8 | **Files are harder to change than code.** | Code is corrected in the next release. A file already written stays as it is for years. | How today's code reads data is evidence when the pattern is sound. Where it is not, the format is chosen on its merits and the code change is listed in [section 5](#5-code-changes-this-implies): the members of a control group are a table, not text inside the group's name (T4). |

---

## 1. The format at a glance

```
my_analysis/                          ONE ANALYSIS. Any number per dataset, stored anywhere.
├── analysis.json                     which dataset, which modality is the response; format version
├── parameters/                       ONE FOLDER PER COMPLETED STAGE
│   ├── record.json                   the stage's settings and files; written last
│   ├── discovery_pairs.parquet       the pairs to test
│   └── positive_control_pairs.parquet
├── assignment/
│   ├── record.json
│   ├── assignments.parquet           one row per assigned (gRNA, cell)
│   └── flagged_cells.parquet         optional: cells the method flags
├── qc/
│   ├── record.json
│   ├── cells.parquet                 one row per cell: kept or not
│   ├── discovery_pair_qc.parquet     which pairs pass pairwise QC
│   └── positive_control_pair_qc.parquet
├── calibration/
│   ├── record.json
│   ├── negative_control_groups.parquet   which gRNAs form each control group
│   ├── results.parquet               one row per tested pair
│   └── resamples.parquet             optional: resampled statistics
├── power/                            record.json, results.parquet, resamples.parquet
├── discovery/                        record.json, results.parquet, resamples.parquet
└── fits/                             A CACHE, NOT A STAGE (T6)
    └── response/                     one part file per job, each saying inside how its fits were made
```

Any `.parquet` table may instead be a folder of part files with the same name.

The data file and each analysis carry separate format versions, because the
data file is read on its own and the two can change independently.

---

## 2. Types in tables and records

One rule per kind of value, the same in both languages. The types of the data
file are in [its rationale](../perturb-seq-data-format/RATIONALE.md#3-number-and-text-types).

| Value | Type in the file | Why |
|---|---|---|
| Integers in tables | 32-bit signed integer | R's integer is 32-bit. *(measured)* R's Parquet reader returns 64-bit integers as doubles. No value in the format needs more than 32 bits. |
| Real numbers | 64-bit float (double) | p-values need the range. Never 32-bit. |
| True/false in Parquet tables | Parquet boolean | Round-trips both ways, with missing values. *(measured)* |
| Text | UTF-8 | In Parquet, plain or dictionary-encoded. |
| A list of numbers per row | Parquet list of doubles | *(measured)* Round-trips both ways between the two Parquet libraries. Used for fit coefficients and resampled statistics. |
| Missing values | Parquet: null. JSON: `null`. | No NaN or infinity in JSON (A3). |
| Sort order | byte order of the UTF-8 text | *(measured)* R's radix sort and Python's byte sort agree; R's default sort follows the locale and does not. |

**IGVF.** Effect sizes in
result tables may be rounded to 32-bit floats, with p-values left as doubles
(`bin/compact_result_dtypes.py`). This format keeps every real number a double:
one rule, and no rounding between what a package computed and what it reads back.

---

## 3. The analysis folder

### A1. Linking an analysis to its dataset and its response modality — *Adopted*

- **Decision.** `analysis.json` records:
  - `format`, the analysis's format version;
  - `response_modality`, the name of the modality this analysis tests (S2). It is always written, even when the file has only one candidate;
  - `screen`: the data file's `dataset_id` (S10), as `id`, the data file's format version, its numbers of cells, responses and gRNAs, a SHA-256 of each ordered identifier list, and a relative path as a hint.

  A reader refuses a dataset whose `id` or hashes differ, or that has no modality of that name.
- **The response modality is found here, never by a fixed name** *(Decided)*. The data file names its modalities for what was measured and marks none as the one to test (S2). So a reader takes the response modality from `analysis.json`, and ignores the modalities the analysis does not name: a file may carry a hashing modality, for instance.
- **In the packages.** Both packages call the response modality `response` in memory today, and users write `response_n_umis` in formulas. The packages' own readers and writers translate between that and the file's names, so no package behaviour changes *(measured: every round trip in [Appendix B](#appendix-b-what-the-prototype-implements) passed)*.
- **Not in this format.** Testing two modalities in one analysis. Under this naming it is two analyses of one file.
- **Why the response modality is here and not in the parameters record.** The link's response count and hash are those of the response modality, so the two belong together, and nothing in an analysis survives a change of response modality.
- **Alternatives for the link.**
  - A path only: breaks when either side is moved.
  - Size and checksum of the data file (the fork's newest branch, *verified, `io/analysis_files.py` on `sceptre3-dsfix`*): breaks if the file is recompressed, and a checksum of a 100 GB file is slow.
  - Requiring the analysis to live inside the dataset's folder: impossible for a shared, read-only dataset.
- **Why two checks.** The `id` survives moving and copying and changes when the dataset is re-imported. The identifier hashes catch a file that was subset by a generic tool that kept the `id`. *(measured)* Hashing 10 million cell identifiers takes 0.3 s in Python and 2.4 s in R.
- **Not caught.** Counts edited in place with the identifiers and `id` unchanged.
- **The data file's format version is recorded, not checked.** `screen.format` is the data file's `perturb_seq_format` when the analysis was created. The two formats have separate versions, so an analysis says which version of the other it was written against. It is for people and tools; the identity a reader checks is the `id` and the hashes.
- **Defaults are for packages.** A package may choose the response modality for the user when the file has one candidate. The file is always explicit.
- **Convention.** Analyses live in `analyses/<name>/` beside the data file unless there is a reason not to.
- **Analyses far from their data** *(Decided; details Adopted)*. A data file may sit read-only on a server while its analyses live on a laptop or elsewhere on the server. The link allows this: it is the `id` and the hashes, a reader is told where the data file is and never follows `screen.path`, the data file knows nothing about its analyses, and an analysis folder can be read for its results with no data file at all. Four rules support it:
  - *The hint is a relative path or a URL, never an absolute local path.* An absolute path names a user and a server's layout, and means nothing on another machine. The schema enforces the form and the validators report it (A1.4).
  - *A writer records a path only when the data file is near*: in the analysis folder's parent or grandparent, or below one of them. Farther away, the relative path is a chain of `..` that names directories and helps nobody, so the writer records nothing.
  - *A dataset may have a `title`*, for people, stored in the data file and copied into the link. An error can then say which dataset an analysis was run on instead of quoting a random `id`. The name is not part of the identity check.
  - *A dataset may have a `source`*: a URL, DOI or accession, such as a GEO accession for a published screen.
- **No dataset version** *(Decided)*. A version would be bookkeeping for every writer. The `id` changes when a dataset is imported again, which is the protection an analysis needs.
- **Finding the data is package behaviour, not format.** A package can resolve a dataset `id` to a location through a search path or a small table in the user's own configuration. Nothing specific to one machine goes in `analysis.json`. If a package ever follows the hint by itself, it should do so only for locations the user has configured: a folder from an untrusted source could otherwise point a reader at any file or URL. None of this is built.
- **Precedents** *(not independently checked)*. The common pattern is to identify data by a name or a hash, keep its location as a separate, replaceable mapping, allow only relative paths or URLs, and carry a label for people.
  - *Frictionless Data Package*: a file reference is a relative POSIX path or an HTTP(S) URL; absolute paths are disallowed, as a security risk and as rarely useful.
  - *BIDS*: references name a dataset, not a location; a separate table maps each name to a local or remote address; source datasets are recorded by URL, DOI and version.
  - *DVC*: data is identified by a content hash with a relative path; the physical location is a named remote configured outside the pointer file.
- **IGVF.** Its checkpoint manifest links a continuation to its parent by a SHA-256 of the whole parent file, plus cell and feature counts, the pipeline commit and a hash of the gene annotation file (`bin/preinference_checkpoint.py`). That catches an in-place edit of the counts, which this format's link does not. It costs a full read of the file and breaks if the file is recompressed. An optional whole-file checksum could be added here for those who want that guarantee ([section 7](#7-fixed-now-free-to-grow)).

### A2. One record per stage, written last — *Adopted*

- **Decision.** Each completed stage has a folder named for it, holding its tables and `record.json`, which is written after the tables. The record holds the settings that stage used, small summaries, provenance (package, version, language, time), and the list of files the stage wrote with their row counts.
- **Alternative.** One `analysis.json` holding everything, rewritten by each stage.
- **Why.**
  - It keeps principle 3: no file is ever rewritten. The single record is the one file every stage would modify.
  - sceptre's settings arrive stage by stage (`assign_grnas(method = ...)`, `run_qc(...)`) *(verified, both codebases)*, so they have a natural home.
  - The record is a commit marker. A stage whose record is missing did not finish, whatever files exist.
  - The calibration check, power check and discovery analysis do not depend on one another's results *(verified, same rank in both codebases)*. A pipeline could run them at once, and they would race to rewrite a single record.
- **Cost.** Up to six small files instead of one. Reading all settings means reading several.
- **A folder per stage** *(Decided)*. Each stage's record and tables are together in a folder named for the stage. The alternative is records in one folder (`stages/`) and tables grouped by kind (`pairs/`, `pair_qc/`, `results/`, `fits/`), with stage ownership recorded only inside the records.
  - *Grouping by kind mixes two organising ideas*, and uses two vocabularies for one three-way split: pair sets named by pair type and results by stage, so `pairs/positive_control` goes with `results/power`.
  - *Grouping by kind needs rules only to police the layout.* A record could list another stage's table, and clearing a stage means finding its files through its record. With a folder per stage, a record names files in its own folder and nothing else (the schema enforces it), and clearing a stage is deleting its folder. The rule against unlisted tables remains, confined to one folder, and still catches a stage that did not finish.
  - *Parallel testing stages write into separate folders*, which suits a pipeline.
  - *Cost.* Results are in three folders, not one `results/` folder.
  - *One name chosen with care:* the QC verdicts are `discovery_pair_qc.parquet` and `positive_control_pair_qc.parquet`, not `discovery_pairs.parquet`, so that a file copied out of its folder is not mistaken for the pairs-to-test table of `parameters/`.
- **Stage names.** `parameters`, `assignment`, `qc`, `calibration`, `power`, `discovery`: neutral words, not R function names.
- **`provenance.language` is free text** *(Adopted)*. The two libraries here write `R` and `Python`. A list of allowed values would shut out a third implementation. A saved fit is reused only between languages that share the fitting method (T6), so a fit from an unknown language is simply not reused.
- **IGVF.** Settings are spread over several places: a resolved run configuration (`full_run_config.json`), the resolved control regime (`inference_mode.json`), a list of covariate options in the file's `uns`, and per-step JSON files in the work folders. They describe one pipeline run. This format's stage records describe one analysis stage and are what a reader in either language needs to continue it.

### A3. What a stage record contains — *Adopted*

- **Settings** are what the user chose, with defaults filled in: everything needed to rerun the stage. They use plain words, not R's internal codes:

  | Setting | Values | Replaces |
  |---|---|---|
  | `side` | `left`, `both`, `right` | `side_code` −1, 0, 1 |
  | `control_group` | `complement`, `nt_cells` | `control_group_complement` true/false |
  | `resampling_mechanism` | `crt`, `permutations` | `run_permutations` true/false |
  | `resampling_approximation` | the user-facing name | internal names such as `crt_spa_fast` |
  | `formula` | text | an R formula with its environment |

- **Names inside settings.** A formula names covariates by the modality's name (`rna_n_umis`, T7). What stays in role words: table columns (`response_id`, `grna_id`, `grna_target`) and setting names (`response_fit_method`, `response_n_umis_range`). A pair is a target and a response whatever was measured.
- **`moi`** *(Decided)*. Whether a screen is treated as low or high multiplicity of infection selects an analysis's defaults: the control group, and what is done with cells that carry several gRNAs. The same data can reasonably be analysed either way, and real screens are not all one or the other. So it is a setting of the parameters record, and the data file has no such entry.
- **Summaries** are small derived counts kept for display (cells removed by each rule). A reader must not need them.
- **Not included.** R-only values: the design matrix's `assign` and `contrasts` attributes, table classes, and the old pipeline flags (`nuclear`, `nf_pipeline`, `integer_id`, `elements_to_analyze`).
- **JSON rules.** UTF-8. Missing is `null`. No NaN or infinity. Each field's type is fixed by `SPEC.md` and its schemas; a reader does not infer integer-versus-double from how a number is printed. Doubles are written with 17 significant digits, so that the number read is the number written. *(measured)* With 15 digits, the stored significance threshold (T5) came back a hair smaller and misclassified the boundary pair.
- **The resampling budget** (`B1`, `B2`, `B3`) is in the parameters record. QC and the calibration check may reset `B3` *(verified, `R/s4_analysis_functs_1.R`)*, so their records carry `B3` too, and the latest stage's value is the one in force.
- **IGVF.** The control group is recorded with the same two words, `complement` and `nt_cells` (`bin/inference_sceptre.R`). The formats agree.

### A4. Which stages have run, and rerunning — *Adopted*

- **Decision.** A stage has run if and only if its record exists. To rerun a stage, a writer first deletes the records and files of that stage and of every later stage, then writes.
- **Order.** parameters → assignment → qc → {calibration, power, discovery}. The last three are independent.
- **Why.** The two codebases differ today. The fork clears later stages *(verified, `analysis/object.py`)*. R only marks them as not run and leaves their values in the object *(verified, `R/s4_helpers.R`)*. On disk, stale files beside a missing record would be read by a careless reader; deleting them makes the folder say exactly what is valid.

### A5. Integrity: row counts in the stage record, no checksum file — *Adopted*

- **Decision.** No `checksums.tsv`. Each stage record lists its files and their row counts.
- **Why.** What a checksum file catches in practice is a truncated or half-written folder: in the fork's earlier layout, sizes were checked on every read and checksums only on request *(verified, fork `io/analysis_files.py`)*. A record written last, with row counts that a reader checks against each table's own footer, catches the same failures without hashing large files on every write.
- **Left out on purpose.** Detection of silent corruption with the size unchanged. A stage record may carry an optional SHA-256 per file for anyone who wants it (principle 6).

### A6. Writing safely — *Adopted*

A writer writes each file under a temporary name and renames it into place.
Parallel jobs each write their own part file into a table's folder. Public
pysceptre's shard scripts already work this way *(verified,
`scripts/run_trans_background.py`)*. In the Nextflow pipeline each stage writes
a new object file in its own work folder and hands it on, so no job appears to
modify a file another job reads *(inferred from `main.nf` and the `bin/` scripts)*.

---

## 4. The tables

### T0. Rules for every table — *Decided*

1. Parquet. Types as in section 2.
2. Rows sorted by the table's key, in byte order.
3. Row groups of at most one million rows. *(measured)* R's Parquet writer defaults to 122,880 and takes the setting.
4. A table is one `.parquet` file or a folder of part files with the same columns. Each part is sorted; parts from parallel jobs need not be sorted across each other.
5. Cells, gRNAs and responses are named by identifier.
6. A reader takes the columns it knows and ignores the rest.
7. An identifier column may be plain text or dictionary-encoded.

Two details that are **not** rules, because the libraries differ: R's writer
marks a column without missing values as required and Python's marks every
column optional *(measured)*; and any of snappy, zstd or gzip compression may
be used. Readers handle all of these already.

Deferred to reader code, not the format: chunked reads, filtered reads.

### T1. `assignment/assignments.parquet` — *Decided*

- **Columns.** `grna_id`, `cell_id`. One row per assigned pair. Key: both.
- **Alternatives.** A 0/1 matrix in the data file (public pysceptre); a matrix in its own `.h5ad`; lists inside the settings record.
- **Why, independent of today's code.** Assignments are tiny next to the counts. Target unions and control groups need many gRNAs at once. Parallel jobs can each write a part. Another tool can hand in a plain two-column table.
- **Supporting evidence from the code.** All four codebases load every assignment: the pipeline merges per-job lists *(verified)*, the fork rebuilds everything from the per-gRNA lists *(verified, `pipeline/r_state.py`)*, and R holds them in the object. Loading everything is sound at realistic sizes: *(measured)* a 10-million-row table reads in 7.8 s in Python and 10.1 s in R, 76 MB on disk (its measurement, not repeated here).
- **Derived on read.** Per-target unions, per-control lists, and gRNAs per cell.
- **Room to grow.** A score column (a posterior probability, a p-value) can be added under rule 6.
- **Deferred: several gRNAs on one vector.** The data file holds per-gRNA counts and is written once (S8, S12), so sceptre's sum-within-vector has to become a step of an analysis. That implies a setting in the assignment stage's record, and an assignment table that names the unit it assigned, gRNA or vector, as results tables already switch their key by strategy. Nothing in the data file depends on the answer. It is listed in [section 7](#suggested-changes-to-what-the-packages-produce). The prototype's readers do not use `vector_id` when they build a package object.
- **IGVF.** A 0/1 cells × guides matrix in the data file, `guide.layers['guide_assignment']`, with an optional second layer of posterior probabilities. pertpy does the same (`layers["assigned_guides"]`). The matrix and the table hold the same information. This format uses the table because the assignment belongs to an analysis, not to the data (S12 of the data file format), and because parallel jobs can each write a part. Since both precedents use the matrix, an export to that shape is worth providing; it is package code, not format.

### T1b. `assignment/flagged_cells.parquet` (optional) — *Adopted*

- **Column.** `cell_id`: cells the assignment method itself judges not to carry exactly one gRNA, beyond what the assignment table shows.
- **Why it exists.** In low MOI, QC removes cells with no gRNA or several. For most methods that follows from the assignment table (a cell with zero rows or more than one). sceptre's "maximum" method is the exception: it assigns every cell its top gRNA and separately flags cells whose top gRNA has too small a share of the cell's gRNA UMIs *(verified, `R/assign_grna_functions.R`)*. Those cells have one assignment row and are still multiplets.
- **Alternatives.** Have each reader re-apply the method's thresholds to the gRNA counts: that puts method logic in every reader, and a new method would need new reader code. Or drop the flagged cells from the assignment table: that loses which gRNA they were assigned.
- **Derived on read.** The full set of cells without exactly one gRNA is the cells with zero or several assignment rows, plus the cells in this table.

### T2. `qc/cells.parquet` — *Adopted*

- **Columns.** `cell_id`, `in_use`. One row per cell. Key: `cell_id`.
- **No `n_grnas` column.** It is the number of assignment rows for the cell (principle 7).
- **Flag.** sceptre records only how many cells each QC rule removed, not which rule removed which cell *(verified, `R/cellwise_qc_functs.R`)*. A `removed_by` column would be more useful and can be added later under rule 6.
- **IGVF.** Cells that fail QC are removed from the file; only the retained cells are stored, with aggregate before-and-after counts in `uns`. A second analysis with different QC then needs a second copy of the data. Here every cell stays in the data file and each analysis says which it used.

### T3. The pairs to test, and their QC — *Adopted*

- **Decision.** The pairs to test are inputs, written once by the parameters stage: `parameters/discovery_pairs.parquet` and `parameters/positive_control_pairs.parquet`, columns `grna_target`, `response_id`. Pairwise QC writes its verdict separately: `qc/discovery_pair_qc.parquet` and `qc/positive_control_pair_qc.parquet`, the same columns plus `pass_qc`.
- **Under the singleton and bonferroni strategies** each gRNA is tested on its own *(verified, `R/preprocessing_functions.R`)*, so the QC verdict is per gRNA and the pair-QC tables have a `grna_id` column as well. The pairs themselves stay as the user gave them, per target; a reader expands them through the gRNA table.
- **The `grna_id` column is a rule, not an option** *(Adopted)*. Whether a pair-QC or results table has the column follows from the strategy and the stage: never under `union`; always under `singleton`; under `bonferroni` in pair QC and calibration, and not in power or discovery, whose results sceptre reports per target. A reader can therefore know a table's key without opening it. The bonferroni irregularity is sceptre's, and the rule changes when sceptre's output does ([section 7](#suggested-changes-to-what-the-packages-produce)).
- **Alternative.** One table with a `pair_type` column and a `pass_qc` column filled in later (the first draft). That makes the QC stage rewrite a file from the parameters stage, against principle 3.
- **Cost.** The key columns are stored twice. *(measured)* A sorted key-only table takes 2 bytes per row, so 0.7 GB at 339 million pairs.
- **Row order.** Not kept. Duplicate pairs are an error *(verified, `R/check_functions.R`)* and results are re-sorted, so order carries no meaning.
- **Not stored.** `grna_group` (S8).
- **IGVF.** The pairs are a table in the data file's `uns` (`pairs_to_test`), one row per (guide, gene). sceptre's unit is the target, so IGVF's sceptre step converts guides to targets before testing (`bin/inference_sceptre.R`). This format stores the pairs as sceptre takes them, per target, in the analysis.

### T4. `calibration/negative_control_groups.parquet` — *Adopted*

- **Columns.** `grna_group` (the group's name), `grna_id`. One row per member.
- **Why.** Today a control group's members exist only inside its name, such as `non-targeting_grna3&non-targeting_grna4`, and both codebases split the name on `&` to recover them *(verified, `R/association_testing_functions_v2.R`; fork `pipeline/calibration.py`)*. That is why `&` is forbidden in gRNA identifiers. Encoding a set in a name is a fragile pattern (principle 8), and an explicit table costs a few rows.
- **Kept.** The names themselves, so results look as they do today. A group's name appears in the calibration results in the column that names the tested unit: `grna_target` under the union strategy, `grna_id` otherwise (where each "group" is one non-targeting gRNA).
- **Not stored.** A separate table of negative-control pairs. `calibration/results.parquet` already has one row per such pair.
- **Honest note.** The prototype's writers fill this table by splitting the names, because the packages hold nothing else today. Readers never split.
- **IGVF.** No counterpart. Its calibration is a separate pass over non-targeting "pseudo-elements" with its own output tables (`pipeline_outputs/calibration/`).

### T5. Results: `<stage>/results.parquet` — *Adopted*

- **Decision.** Three tables, one per testing stage. Key: `grna_target`, `response_id` (and `grna_id` under the singleton strategy). Columns as sceptre reports them, less what follows from others: `pass_qc`, `p_value`, `fold_change`, `se_fold_change`. Diagnostic columns may follow (rule 6).
- **`log_2_fold_change` is always stored** *(Decided)*.
  - It follows from `fold_change`, and principle 7 would leave it out. It is stored because rebuilding it does not give the same number everywhere: the last bit of a logarithm differs between platforms. On the example, 136 of 15,768 values rebuilt on Linux differ from those computed on macOS, by at most 3 parts in 10^16 *(measured)*; and `log2(x)` differs from `log(x) / log(2)`, sceptre's expression *(verified, `R/aux_functions.R`)*, in the last bit of about 30% of rows *(measured)*. Stored, the value a user sees is the one the writer computed, in both languages and on every platform.
  - Under the bonferroni strategy, sceptre's discovery and power results have no `fold_change` and no `se_fold_change`, only `log_2_fold_change` *(verified, `R/s4_analysis_functs_2.R`)*. So it is the one effect size every results table has, and the required one.
  - Where a table has both, the validators check that they agree to a relative difference of 10^-9 (T5.2). That is the guard against two copies of one fact disagreeing, which is what principle 7 is for.
  - **Rejected.** Storing it only where there is no `fold_change`, and having readers rebuild it: that needs a named expression, a tolerance and an equality test in every writer, to save 8 bytes a row. Storing only `log_2_fold_change`: the standard error sceptre reports is on the fold-change scale.
- **Sorted by key, not by p-value.** Sorting by key lets parts be written independently and lets a reader find one target's rows. Packages sort by p-value when they hand the table to a user. This also removes a dependence on job order among tied rows *(inferred from `bin/process_association_analysis_results.R`)*.
- **Size.** *(measured, synthetic, before `log_2_fold_change` was always stored, which adds 8 bytes)* 33 bytes per row, so about 11 GB at 339 million pairs. Reading two columns of 20 million rows takes 0.05 s in Python and 0.8 s in R; as 20 part files, 1.7 s in R.
- **Resampled statistics** (the largest output setting) go in `<stage>/resamples.parquet`: the key columns and one list column, `z_null`. Today packages hold them as thousands of columns `z_null_1`, `z_null_2`, ..., padded with missing values *(verified, both codebases)*. *(measured)* Asked to keep them for the example's 1,406 pairs, the R build used here failed outright; the example keeps them for 40 pairs.
- **The `significant` column is not stored.** It depends on the whole table, so parallel jobs could not write it. The stage record stores `significance_threshold`, the largest p-value called significant (or `null` if none was), and a reader computes `significant` as `p_value <= threshold`. This is exact for every correction sceptre offers, because each calls significant exactly the pairs up to some p-value. *(measured)* It reproduces R's column for every pair of every example.
- **Very small p-values.** A p-value below about 1e-308 becomes zero in a double. A `log_p_value` column can be added when a method needs it (rule 6); none does today.
- **IGVF.** Two tables per scope, per guide and per element, which correspond to this format's singleton and union strategies. They differ from this format in four ways:
  - *Columns carry the method's name* (`sceptre_p_value`, `sceptre_log2_fc`, `sceptre_fc_se`), because one table holds two methods' results. Here an analysis is one method's, so the names are plain.
  - *Annotations are copied into every row* (guide sequence, coordinates, gene symbol, number of perturbed cells). Here a row has keys and results, and annotations are joined from the data file.
  - *Adjusted p-values are stored* (`sceptre_q_value`) along with `sceptre_negLog10p`. Both follow from the p-values. Here neither is stored; the threshold plays the role of the adjusted column.
  - *Each table is stored three times*: in the file's `uns`, as Parquet and as TSV. Here once, as Parquet.

  IGVF's working table format is Parquet with zstd, the same as here. Its per-row count of perturbed cells is a useful diagnostic that these tables do not carry; it is listed in [section 7](#7-fixed-now-free-to-grow).

### T6. The fit cache: `fits/response/` — *Decided (placement; where the identity lives), Adopted (details)*

- **What it is.** A cache of fitted response models. A reader may ignore it, and deleting it changes no result beyond numerical tolerance.
- **Columns.** `response_id`, `coefficients` (a list of doubles, in the order of the design's terms), `theta`. Key: `response_id`.
- **Alternatives for the columns.**
  - One row per coefficient, with the dispersion as the row whose term is missing: a trick.
  - One column per term: a covariate named `theta` would collide, and the columns change with the formula.
- **Outside the stages** *(Decided)*. A response fit depends on the data, the formula, the fitting method, the control group and which cells passed QC. All of that is fixed once the parameters and QC stages have run. Which testing stage happens to compute a fit first is an accident. Keeping fits in the folder of the stage that first made them would:
  - delete still-valid fits when that stage was rerun;
  - stop stages running in parallel from sharing them;
  - name them after a stage that does not own them.

- **How the cache works.**
  - It is a folder of part files. Each job adds a part of its own and never modifies another's, which keeps "every file is written once" and uses the rule that a table may be a folder of parts.
  - Two parts may hold the same response, when two parallel jobs fitted it. A reader takes either. So the rule that a key is unique does not hold across the parts of this table, and the specification says so.
  - Fits are listed in no stage's record. What makes them valid is their **identity**: the dataset `id`, the formula, a hash of the identifiers of the cells in use, the fitting method and the control group. Each part carries the identity of its own fits, inside the file (next point).
  - The cache is deleted when the QC stage or an earlier one is rerun, and kept when a testing stage is rerun. *(measured)* Both libraries do this, and a conformance case covers it. Since the identity moved into the parts, deleting is housekeeping (SHOULD): correctness does not depend on it. It is kept because a plain read of the whole folder would otherwise mix old parts with current ones.
- **The identity is inside each part** *(Decided)*. It is JSON text under the key `sceptre_fit_identity` in the Parquet file's own key-value metadata. There is no `identity.json`, neither in the folder nor beside it; the folder holds part files and nothing else.
  - **Why.** The identity exists to stop a fit being reused when it was computed for other cells, another formula or another method. If that check passes wrongly, the result is wrong p-values with no error. With a separate file, the check reads one file and trusts that it describes the others, and they can come apart: QC rerun with the clearing interrupted; a job adding parts without checking the file; two parallel jobs both finding a stale cache and both clearing and rewriting it. Each is unlikely, and the failure is silent. With the identity inside each part there is nothing to keep in step, a stale part is harmless, and no file is shared between writers, which fits "every file is written once". No writer has to write a shared record or verify one.
  - **Why the other tables do not need this.** A stage finishes. Its record, written last, lists its files with their row counts and so marks the stage as complete, and rerunning the stage deletes the folder first. A record written last cannot work for something open-ended that jobs keep adding to. That structural difference is the whole reason the cache is treated differently.
  - **Alternatives.** A record file inside the folder (`fits/response/identity.json`), or beside it (`fits/response.json`). Both leave one file describing many, which is the weakness above; the first also breaks folder readers, as measured next.
  - **Evidence** *(measured; the first and third lines are repeated in the suite on every run)*:
    - A JSON file inside a folder of parts: `pyarrow.dataset`, `pyarrow.parquet.read_table`, `pandas.read_parquet`, R's `arrow::open_dataset`, `polars.scan_parquet` on the folder and duckdb's `read_parquet('folder/*')` all fail. Only an explicit `*.parquet` pattern works.
    - The same file named `_identity.json`: the Arrow family ignores it and works; polars and duckdb still fail. So an underscore is not enough.
    - Parts only: everything works. All six readers open the cache folder of every example directly.
    - Both directions across languages, with the libraries already in use: pyarrow writes the key and nanoparquet reads it; nanoparquet writes it and pyarrow reads it. Arrow adds a key of its own, `ARROW:schema`, which is ignored. duckdb and R's arrow show the key too.
  - **Costs.**
    - It cannot be read with a text editor, so the validators show it, and `fit_parts()` in both libraries returns it.
    - It is a mechanism the format uses nowhere else.
    - A tool that rewrites a part may drop the metadata. The part is then ignored, which is the safe direction. *(measured)* pyarrow's `read_table` exposes only the metadata that Arrow's own schema record holds, so a part written by R, read and written back with pyarrow, loses the key unless the code copies it from the file's metadata. A reader must take the key from the file's metadata, not from the table.
  - **Parts that are not the analysis's are not reported.** A part with another identity is legal and unused, so the validators say nothing about it. They do report a part with no identity, or one that does not conform, and anything in the folder that is not a part file.
  - **One limit of reading the folder alone.** Four of the five identity fields can be read from the analysis folder. The fifth, the hash of the cells in use, needs the data file's order of cells. So a reader given only the folder returns the fits made for the analysis's parameters when all such parts agree on the cells, and nothing when they do not (parts left by an earlier QC beside current ones). A reader that has the data file gives the full identity and gets exactly its parts. Both libraries work this way.
  - *(not verified)* GeoParquet is believed to keep its metadata the same way, as JSON under one key.
- **When a cached fit is reused.** Only from a part whose identity equals, field by field, what the reader would use itself, and only if the fitting method is one the reader's language shares. Identities are compared as parsed values, with white space in the formula ignored. The identity does not name the language that made a fit: a method name means the same estimator in every language that lists it, and parts may come from several languages.
- **It is the part of the layout a new reader is most likely to misunderstand**, because it looks like a seventh stage and is not: no record, it fills gradually, and it is cleared at a different moment. The specification says so at the head of its section.
- **Room for gRNA fits** *(Decided)*. These will be wanted eventually, and this version must not step on them. The folder is named for its kind, `fits/response/`, and `fits/grna/` is reserved. The two would not be one table with a kind column: they have different keys (a response, against a target or a gRNA depending on the strategy), different columns (no dispersion for a gRNA fit) and different dependencies. A gRNA fit would also depend on the assignments, so its cache would be cleared when assignment is rerun *(inferred from how conditional resampling works, not read in the code; to verify before the folder is specified)*.
- **`response` is unambiguous** because an analysis names one response modality (A1). The folder is not named for the modality (`fits/rna/`), which would put a user-chosen name beside the reserved `grna`.
- **Fitting methods are shared, named definitions** *(Decided in principle)*. The methods are meant to be consistent across R and Python.
  - A method name means the same estimator in both languages, and the names written to files come from one shared list.
  - A method that exists in one language only (glmGamPoi, in R) is marked as such in that list, and its fits are never reused by the other language.
  - How closely two implementations of one method must agree belongs in the hub's parity contract. The format's part is to make it checkable, by recording the method with every saved fit.
- **This closes a gap in both codebases.** Neither compares a saved fit's method with the one about to be used at the moment of reuse. R clears its cache when the method is changed through `set_analysis_parameters()` *(verified, `R/s4_analysis_functs_1.R`)*, but at test time it reuses any saved entry for the response *(verified, `R/medium_level_functs_v2.R`)*, so fits loaded from files are trusted. The fork's reuse depends on a global setting, not on the method the analysis recorded *(verified, `pipeline/discovery.py`)*.
- **How far the methods differ today** *(measured on the example: 100 responses, 1,406 pairs)*:

  | Compared | Dispersion (theta) | p-values | Discoveries |
  |---|---|---|---|
  | R `sceptre` vs R `pysceptre` | agree to 1e-9 | agree to 2.5e-8 | same 11 |
  | R `batched` vs Python's default (`batched`) | identical | agree to 8e-12 | same 11 |
  | R `batched` vs R `sceptre` | 98 responses agree to within 1e-6 (median difference 4e-14). **Two differ completely:** `batched` gives 1000 (its upper bound) where `sceptre` gives 0.07 and 0.03 | 22 pairs differ, all in those two responses; 9 by more than 10%, the largest by a factor of 2 | same 11 |

  So there are two families: `sceptre` and `pysceptre` agree, and `batched` (the current default in both languages) agrees with itself across languages but treats some responses differently. glmGamPoi was not installed and was not measured. Which estimator is the default is a scientific decision, outside the format ([section 7](#suggested-changes-to-what-the-packages-produce)).
- **Not adopted.** Public pysceptre's cache also keys on library versions and CPU features, for bit-exact reuse *(verified, `pipeline/fingerprint.py`)*. That is stricter than a cross-language cache needs.

### T7. The design matrix — *Decided*

- **Decision.** Not stored. The parameters record holds the formula text and the names of the matrix's columns; each language rebuilds the matrix. Only additive formulas of covariate names, `log(x)` and `log(x + 1)` are allowed, and a writer refuses others. No stored-matrix fallback.
- **The covariate table** *(Decided)*. The model's covariate table is the data file's shared cell table together with the cell tables of `grna` and of the response modality, a modality's columns prefixed with the modality's name and an underscore: `rna_n_umis`, `grna_n_umis`. The formula text and the recorded term names use these names (`log(rna_n_umis)`). The data file keeps a shared column from taking such a name (S5).
  - *Alternative.* A role prefix (`response_n_umis`) whatever was measured. A role prefix cannot name the RNA depth in an analysis of protein.
  - *The cost.* Today's covariate names are `response_n_umis` and so on, and users write them in formulas. In a file they become `rna_n_umis`. Most users never type them, because the default formula is built automatically, and a package can keep accepting the old names; the prototype translates in both directions.
  - *Dots in names.* The data file allows a dot in a column name (S5). *(measured)* R takes a dotted name in a formula as it is. pysceptre's formula parser rejects it, so the prototype's Python reader gives the column a dot-free alias for the call; with that, R and Python build bitwise-equal matrices for `~log(rna_n_umis) + S.Score + batch`. The parser change is listed in [section 7](#suggested-changes-to-what-the-packages-produce).
  - *Locks in.* The prefix rule and the underscore. MuData's own `rna:n_umis` colon style is not used, because a colon is not valid in a formula.
- **The baseline of a factor.** The data file stores a categorical's levels in order (S7). The first stored category is the baseline of the model. This is R's default contrasts; nothing is invented.
- **Evidence.** *(measured)* On the example, R's `model.matrix()` and pysceptre's `model_matrix()` give bitwise-equal matrices on one machine.
- **What a reader checks.** The rebuilt matrix's term names. The names (`batchbatch_2`) encode the factor levels, the baseline and the term order, which are the ways two languages can really disagree. That the data is the same data is the job of the link to the data file (A1).
- **IGVF.** The design is recorded as a list of covariate options in the file's `uns`, and each method builds its own design from them. There is no check that two methods, or two languages, built the same matrix.
- **No hash of the matrix.** A hash of exact values fails on a correct matrix, because a logarithm's last bit differs between platforms: the natural logarithms of the whole numbers 1 to 200,000 differ in 16 values between macOS and Linux *(measured, numpy on both)*. With a logarithm in the default formula, an analysis written on one would be refused on the other. A check that survives this, such as a hash of the covariates before any arithmetic, can be added later ([section 7](#7-fixed-now-free-to-grow)).
- **Restriction to state plainly.** R sceptre accepts any formula today. Under this rule an analysis with an interaction or a spline cannot be saved.

---

## 5. Code changes this implies

The format is chosen on its merits; these are the places where package code
would follow, for the data file as well as for the analysis. The prototype's
writers already do the file side of each.

| Where | Change | Why |
|---|---|---|
| R writer | Write counts as 32-bit integers and both index arrays in the index type; write file metadata as a native record | S4, S10 |
| Python writer | Turn off `nullable-string-array`; write 32-bit counts (scipy already gives the index type) | S7, S4 |
| R sceptre | Check at import that counts are integers | S4 |
| Both | Convert text and true/false covariates to categoricals at import | S7 |
| Both | Keep control-group membership as data instead of splitting names on `&` | T4 |
| Both | Check a saved fit's recorded method before reusing it; keep one shared list of fitting-method names | T6 |
| Both | Sort results by p-value when returning them, not when storing them | T5 |
| Both | Hold resampled statistics as one list per pair, not thousands of columns | T5 |
| R sceptre | Clear later stages' values when a stage is rerun, not only their flags | A4 |
| R sceptre | Record which QC rule removed each cell | T2 |
| pysceptre | Read and write these files in the package, not only in scripts | S1 |
| Pipeline | Replace the `.rds` object passed between stages with stage files | A2, A6 |
| Both | Translate between the packages' `response_` covariate names and the file's modality names, as the prototype does; later, accept the modality names directly | A1, T7 |
| Both | Translate the gRNA table's `chr`, `start`, `end` to and from `target_chr`, `target_start`, `target_end` | S8 |
| Both | Store `log_2_fold_change` as computed, and take the stored value on read | T5 |

---

## 6. Compared with IGVF and pertpy

How the data files compare is in [the data file format's
rationale](../perturb-seq-data-format/RATIONALE.md#5-compared-with-igvf-and-pertpy), which also says where the
facts about IGVF and pertpy come from.

### IGVF

**Where the two agree**

| Area | Both |
|---|---|
| Identifiers | Table rows are named by identifier, not position |
| Control group | `complement` or `nt_cells` |
| Two levels of testing | Per element and per guide there; union and singleton here |
| Result tables | Parquet with zstd as the working format |
| Continuing an analysis | The parent data is never overwritten, and a hash ties the continuation to it |

**Where they differ**

| Area | IGVF | This format | Which is better, and why | See |
|---|---|---|---|---|
| Data and analysis | One file holds both | Separate | This format: one dataset, many analyses. IGVF's own code shows the cost of patching results into the data file | S12 |
| Cells | Only cells that passed QC | All cells; each analysis says which it used | This format: QC is an analysis choice | T2 |
| Assignments | 0/1 matrix layer in the data file | Table in the analysis | This format, for the same reason. An export to the matrix is worth having | T1 |
| MOI | `high` or `low`, in the data file | `high` or `low`, in the analysis's parameters | This format: it selects an analysis's defaults, and one dataset can be analysed either way | S10, A3 |
| Results | Three copies; method-prefixed columns; annotations in every row; adjusted p-values stored | One copy; keys and results only | This format: nothing stored twice | T5 |
| Settings | Several run-level files | One record per stage | This format, for continuing an analysis | A2 |
| Link to the data | SHA-256 of the whole file | Dataset id and identifier hashes | Mixed: IGVF's catches in-place edits; this format's survives recompression and is fast | A1 |
| Names | `guide_id`, `gene_id`, `intended_target_name` | `grna_id`, `response_id`, `grna_target` | Neither is a standard. An importer maps one to the other | A3 |

**What an IGVF importer has to do:** turn the assignment layer into a table,
and move the pairs and results into an analysis folder.

### pertpy

pertpy defines no file format and works on in-memory AnnData objects.

| Area | pertpy | This format |
|---|---|---|
| What a perturbation is | A text label on each cell (`perturbation`, control label `control` or `NT`) | Rows of (gRNA, cell) in a table |
| Several gRNAs in a cell | Names joined with `+`, capped at five, then `multiple` | One row each |
| Assignments | 0/1 matrix layer, `assigned_guides` | Table in the analysis |
| Results | In-memory tables: `variable`, `p_value`, `adj_p_value`, `log_fc`, `statistic` | Parquet tables |

The per-cell label works at low MOI and breaks down at high MOI; a label can
be derived from this format's table when a tool wants one. pertpy's `log_fc`
is base 2 but does not say so.

---

## 7. Fixed now, free to grow

### What would be expensive to change after the freeze

| Kind | What |
|---|---|
| Structure | One record per stage; tables keyed by identifier; a folder per stage |
| Names | The table and column names in sections 3 and 4; the stage names; the covariate prefix rule (T7) |
| Removals | A column or file that readers rely on cannot be taken away |

Everything else can be added with a minor version, because readers ignore what
they do not know (principle 6), or loosened, because old files stay valid. The
table of planned additions below says which additions are not safe to ignore.

### To settle before the freeze

- **Review** (`overview.html` says where to look).
- **An analysis at scale.** The tables of an analysis have been written and read at real size with made-up content (`scale-test/README.md`): results tables are no problem at 100 million pairs; resampled statistics take 4 GB per million pairs, and a row group of a million rows (the limit of T0.3) makes reading one pair cost 7 seconds and 11 GB where groups of 10,000 rows cost 0.1 seconds; reading the fit cache costs about 9 milliseconds a part file, so 10,000 parts take a minute and a half. What to change in T0.3, T5.5 and T6.4 because of this is open. No real analysis larger than the 1,000-cell examples has been written or read in this format.

### Planned additions to the format

None of these is in format 0.1. The last column says what kind of version
change each would be once the format is frozen: **minor** if an older
reader can safely ignore it, **major** if it could misread the file.

| Addition | Where | Why wait | Version |
|---|---|---|---|
| A hash of the data file's contents: the cell and feature tables, and perhaps the counts | `analysis.json`, beside the identifier hashes | The link (A1) catches a data file whose cells or features changed, not one whose values were edited in place. Hashing the tables is cheap; hashing the counts is a full read of the file | minor |
| A check that a rebuilt design matrix is the writer's | Parameters record | A hash of exact values is not portable (T7). A hash of the covariates before any arithmetic would be | minor |
| A small stamp in every stage table's Parquet metadata: format version, dataset `id`, stage, table name | Key-value metadata of each table file, as the fit cache's parts carry their identity (T6) | A convenience, not needed for correctness: a results file copied out of its folder would identify itself, and a validator could confirm each file belongs. It is a second copy of facts the stage record holds, so the two would have to be kept in agreement. | minor: an old reader ignores the metadata |
| Why each cell was removed (`removed_by`) | `qc/cells.parquet` | sceptre does not record it (T2) | minor |
| An assignment score (posterior probability or p-value) | `assignment/assignments.parquet` | Depends on the assignment method | minor |
| Per-gRNA values an assignment method computes, such as fitted thresholds | A new optional table of the assignment stage, keyed by `grna_id` | They are computed by an analysis, so they are not data and not settings. Nothing needs them now: sceptre's thresholding uses one threshold, stored as a setting | minor |
| `log_p_value`, for p-values below a double's range | Results tables | No method needs it today | minor |
| Per-pair counts: cells carrying the perturbation, and cells with a nonzero count in each group | Results or pair-QC tables | The development builds do not put them in the result tables | minor |
| A SHA-256 per file, and one of the data file | Stage records, `analysis.json` | For those who want IGVF's stronger integrity check | minor |
| Cached fits of the gRNA models | `fits/grna/`, a reserved folder | Not needed yet; its key, columns and clearing rule differ from the response cache's (T6) | minor |
| Two measured modalities tested in one analysis | Analysis folder | No dataset needs it yet | major: it changes what `response_modality` means |

### Suggested changes to what the packages produce

These are changes to sceptre and pysceptre, not to the format. None is made
here: the packages' functionality is frozen while their foundations are rebuilt. Each would
remove a workaround in the format or its readers.
The suggestions about the data file are in [its
rationale](../perturb-seq-data-format/RATIONALE.md#suggested-changes-to-what-the-packages-produce).

| Change | What it would remove | Evidence |
|---|---|---|
| Bonferroni results carry `fold_change` and `se_fold_change` of the reported gRNA, like the other strategies | Results tables with `log_2_fold_change` and no `fold_change` (T5) | `R/s4_analysis_functs_2.R` keeps only the log2 value |
| One expression for `log_2_fold_change` everywhere | A last-bit difference waiting to happen: pysceptre has `np.log2(fc)` in `pipeline/analysis_result.py` and `log(fc)/log(2)` in `pipeline/result_fields.py` | *(verified)* |
| Keep control-group membership as data | Writers splitting group names on `&`; the ban on `&` in gRNA identifiers (T4) | Both codebases split names |
| Record which QC rule removed each cell | The summaries-only QC record (T2) | `R/cellwise_qc_functs.R` |
| Hold resampled statistics as one list per pair | Thousands of padded columns; the failure when keeping them for 1,406 pairs (T5) | *(measured)* |
| Sort results by p-value when showing them, not when storing them | Order depending on job order among ties (T5) | *(inferred)* |
| Check a saved fit's recorded method before reusing it; keep one shared list of method names | Fits trusted without a check (T6) | Both codebases |
| Decide the shared default fitting method | `batched` gives a dispersion of 1000 on 2 of 100 responses where the other two methods give 0.03 to 0.07 (T6) | *(measured)* |
| Clear later stages' values when a stage is rerun | Stale values behind a "not run" flag in R (A4) | `R/s4_helpers.R` |
| Resolve a dataset's `id` to a location from the user's configuration (a search path or a small registry) | Passing the data file's path by hand every time (A1) | Not built |
| Accept modality-name covariates (`rna_n_umis`), keeping `response_n_umis` as an alias; add a response-modality argument and a modality name at import | The translation layer in the readers and writers (A1, T7) | — |
| Let pysceptre's formula parser take covariate names with dots | An alias where pysceptre reads the file; a pysceptre analysis with such a covariate fails today (T7) | *(measured)* `model_matrix` raises on `S.Score` |
| Combine a vector's gRNAs in the assignment stage, not at import: a setting in the assignment record, and an assignment table that names the unit assigned | Summing at import, which a written-once data file of per-gRNA counts cannot hold (T1) | sceptre sums within a vector in `import_data()` |

---

## 8. How the format is defined and checked

The method is the data file format's: numbered rules of three kinds,
definitions a machine can read, conformance files and validators. [Its
rationale](../perturb-seq-data-format/RATIONALE.md#7-how-the-format-is-defined-and-checked) explains it, with
what counts as evidence, what is taken from Apache Arrow's practice and where
the code lives. This section gives what is particular to analyses.

### The pieces

| Piece | Where | What it is for |
|---|---|---|
| The rules | `SPEC.md` | Numbered statements, each a file rule, a behaviour rule or guidance, with a scope, a glossary and the library versions it is defined against. A validator's report cites them (`T5.2`) |
| The introduction | `overview.html` | One page: what is in an analysis folder, the choices that are expensive to change with their reasons in short, the questions that need a view, and what has been tested |
| The reasons | this document | Why each rule is what it is |
| Table definitions | `spec/tables.json` | Every table's stage, key, columns and types. Both prototype libraries and both validators read this one file, so the two languages cannot drift apart |
| Record schemas | `spec/analysis.schema.json`, `spec/stage.schema.json`, `spec/fits.schema.json` | JSON Schemas (draft 7) for `analysis.json`, the six stage records, and the identity each part of the fit cache carries |
| Conformance files | `conformance/`, built by `verify/make_conformance.py` | 16 conforming analysis folders, written by R and by Python, and 34 that each break a rule, with the rules a validator must report. Each names its data file, which is among the data file format's files |
| Validator | `validate_analysis()` in `format/sceptre_analysis.R` and `format/sceptre_analysis.py` | Checks an analysis folder, with its data file when one is given, without building a package object, and reports every problem, not only the first |

**What the specification leaves out.** It fixes the form of results and not
their meaning: what a p-value tests is defined by sceptre's method
documentation. Left for later: each result column's meaning, tighter wording
of individual rules, schemas for method-specific settings (an assignment
method's `hyperparameters` is any JSON object), and an implementation by
someone who has only the specification.

### What counts as evidence

Beyond the points in the data file format's rationale:

1. **Whether sceptre and pysceptre compute the same thing from the same files
   is not the format's job.** That is parity between the packages. The
   format's part is to make the comparison possible: it records the formula,
   the methods, the seed, the settings and the writer; it stores the design's
   term names; and it stops fits crossing where the method differs.
2. **Sufficiency is the format's job, and needs a package.** Whether the files
   hold everything an analysis needs can only be tested by rebuilding a
   package's object from them.
3. **What a standard reader returns still has to be interpreted.** The link
   between an analysis and its data and the covariate table with modality
   prefixes are each package's own code, as the data file format's rationale
   says of the data file's interpretation.

### What is checked

**The verification suite** (`verify/run.sh`) loads neither sceptre nor
pysceptre and runs on the committed files with standard libraries only.
*(measured, on macOS, Linux and Windows)* Each group of checks runs in R and in
Python, all passing:

- the validators on the conformance set. Both validators accept all 16
  conforming analyses and report exactly the expected rules for all 34 broken
  ones, the same in both languages. 22 rules are exercised by at least one
  broken analysis. One more is checked with no file: T0.3 (row groups of at
  most a million rows), because the file would be large;
- two writers on seeded random content, a data file and an analysis of it,
  compared array by array and with the source;
- readers written elsewhere, on the tables: pandas, polars, duckdb, arrow;
- an encoding matrix for the tables: each kind of column written by each
  language's standard Parquet library and read by the other's;
- what a writer does beyond producing a conforming analysis: a stage written
  again leaves nothing of its first attempt;
- what a reader derives from the files, against values sceptre computed,
  stored beside the examples.

**The overview** (`overview.qmd`, rendered to `overview.html`) is the
introduction, not a check: it shows one analysis as an annotated tree, with
descriptions read from `spec/`, and reports what the verification suite found
when it was rendered. It renders R and Python in one process *(measured,
macOS only)*. R's arrow and Python's pyarrow cannot both be loaded in that
process, so its R side reads Parquet with nanoparquet.

**What these cannot show.**

- *Sufficiency.* Appendix B records the one time it was checked, with both
  packages. Nothing here repeats that check, so a change to the format is
  covered only when the packages test it.
- *Anything about results.*

**The example files are fixed.** Nothing here regenerates them, because doing
so needs sceptre. A format change converts them in place, with code kept
beside the change, which is also a direct test that old files can be carried
forward.

### How adoption would go

None of this has started, and it is not a format question.

1. **Finish the format side first**: the verification suite, the harness, and
   a review of structure and names. Package code written before that review
   risks being rewritten.
2. **Settle what the package branches start from.** This is the real blocker.
3. **In each package, two layers**: the plain layer, lifted from the prototype,
   and a small new adapter to the package's own object.
4. **Small changes in this order**: R writer, R reader, Python reader, Python
   writer. The writer comes first in R because R is the reference
   implementation and a writer has independent evidence at once (the validators
   and the standard Python libraries), while a reader has nothing real to read
   until a writer exists.

---

## Appendix A. Measurements

One run each, on one Mac, with warm file caches. Synthetic data unless noted.
The library versions are those of [the data file format's
measurements](../perturb-seq-data-format/RATIONALE.md#appendix-a-measurements).

**Tables:**

| Test | Result |
|---|---|
| Results table, 20 million rows, 8 columns | 33 bytes per row with zstd, 35 with snappy |
| Read 2 columns of it | 0.05 s Python; 0.8 s R; 1.7 s R from 20 part files; identical to the single file |
| Read all 8 columns in R | 1.9 s, 917 MB in memory |
| Sorted key-only pairs table, 20 million rows | 2.0 bytes per row |
| List-of-doubles column | Python → R and R → Python both correct |
| 64-bit integer column read in R | arrives as a double (exact up to 2⁵³) |
| R writer's default row group | 122,880 rows; 1 million when set |
| SHA-256 of 10 million cell identifiers | 0.3 s Python; 2.4 s R |
| Sort order | R radix sort = Python byte sort; R default sort differs |

**Records, formulas and results:**

| Test | Result |
|---|---|
| Design matrix, R against pysceptre, on the example | bitwise equal |
| Fitting methods on the example (real data, 500 cells) | see the table in T6 |
| A double written to JSON with 15 significant digits | the significance threshold came back smaller; 1 pair misclassified. With 17 digits: exact |
| R sceptre keeping resampled statistics as columns | fails for 1,406 pairs; works for 40 (4,999 columns) |
| `log_2_fold_change` recomputed as `log(x) / log(2)`, 17,219 result rows of the examples | identical to the stored value on every row, in R and in Python |
| The same recomputed with `log2(x)` | a different last bit on about 30% of rows |
| `fold_change` recomputed from the log2 value | a different last bit on a few rows per table |

**At real size.** In the scale run on a real screen of 207,324 cells, which
[the data file format's rationale](../perturb-seq-data-format/RATIONALE.md#scale-test-on-real-data) reports:

| Check | Result |
|---|---|
| An assignment table of 6.6 million rows, written by Python and rewritten by R | Equal |
| Read the assignment table | 5 s in Python; 1 s in R |
| Write the assignment table | Python: 111 s, sorting rows one at a time; not rerun with the numpy sort the writer uses now. R: 1 s |

The tables of an analysis were also measured at real size with made-up
content; `scale-test/README.md` has those numbers. No real analysis at that
size has been written or read.

## Appendix B. What the prototype implements

Everything in sections 2 to 4, and the validator of section 8. The code is
`format/sceptre_analysis.R` and `format/sceptre_analysis.py`, which load the
data file format's code; `example/` holds analyses in the format.

**Sufficiency, checked once with both packages.** Each example analysis was
rebuilt as a sceptre object and as a pysceptre `Analysis` from the files. The
code that did so needed development builds of both packages and is not in this
folder (commit `7e10af4` holds it). Nothing here reruns these checks:

| What | Result |
|---|---|
| R writes, R reads back as a `sceptre_object` | Every field identical or equal for every example analysis, on high MOI and low MOI |
| R writes, Python reads into pysceptre's `Analysis` | Equal, part by part, to what pysceptre's own reader gives for the same analyses in its own layout |
| Python rewrites everything, R reads the Python copy | Objects as above; all stage records equal; all tables identical; same number types in the data file |
| Python continues an analysis R stopped after QC (parity between the packages, not evidence about the format) | p-values agree with R's to a relative 8e-12; same discoveries |
| R resumes from Python-written files | Results identical to an analysis that was never saved |
| Per-target unions, per-control lists, gRNAs per cell derived from assignments | Equal to what the packages hold |
| Design matrix rebuilt from the formula | Same term names and the same hash in R and Python |
| Resampled statistics as a list column | 40 pairs, 499 to 4,999 statistics each, restored as the packages' columns |
| The three gRNA strategies, and low MOI with the "maximum" method, flagged cells and NT cells as controls | All round-trip |
| A table as a folder of part files, with an unknown column and dictionary-encoded identifiers | Read to the same object, from parts written by either language |
| The same data under the modality name `atac`, with a `hashing` modality beside it | Read to the same object by both; without an analysis, a reader asks which modality is the response |
| A subset data file, and a re-imported one | Refused |
| A table with fewer rows than its stage wrote | Refused |
| A covariate with a dot in its name | R and Python build the same design matrix (Python through an alias) |

**Limits that remain.** These are limits of the prototype or of what was
tested, not open choices:

- Nothing here writes a stage while an analysis runs: the writers are given
  finished tables. Writing stage by stage belongs in the packages.
- A value may come back equal but not identical: a whole-number setting such
  as a threshold of `5` is an integer in JSON.
- Refused with a message: a formula outside the additive subset.
- The packages call the response modality `response`; a package's reader and
  writer translate. The check of a response modality under another name
  (`atac`) reused the RNA counts.
- The tables of an analysis were measured at real size with made-up content,
  and the suite passes on Linux as on macOS. Not tested: a real analysis at
  scale, and methods other than those in the examples.

## Appendix C. The examples

| Dataset | Analyses |
|---|---|
| `example/highmoi` | Of sceptre's built-in high-MOI data (500 cells, 100 responses, 60 gRNAs): `left_sided`, `two_sided`, `singleton`, `bonferroni`, `stopped_after_qc` (through QC only), and `with_resamples` (discovery only, 40 pairs, every resampled statistic) |
| `example/lowmoi` | Of sceptre's built-in low-MOI data (1,000 cells, 100 responses, 50 gRNAs): `default`: the "maximum" assignment method and the NT cells as controls |

The data files are the examples of the data file format, under
`../perturb-seq-data-format/example/`. These files are fixed: nothing here
regenerates them, and `example/*/expected/` holds values sceptre computed for
each analysis.
