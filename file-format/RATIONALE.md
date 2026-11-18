# The sceptre file format: rationale

[`SPEC.md`](SPEC.md) states the rules of the format, numbered, with no reasons.
This document gives the reasons: for each choice, the alternatives, the
evidence, what it locks in, and how it compares with the IGVF pipeline's
format. Sections here and there share their letters and numbers (S3, A1, T5,
...), so rule `S4.1` of the specification is explained in section S4 below.
[Section 11](#11-how-the-format-is-defined) says how the format is formally
defined and checked.

**Status:** draft (format 0.1), 18 November 2026. Every choice below is settled
for this draft and implemented in the prototype beside this file. None has yet
been reviewed. [`REVIEW.md`](REVIEW.md) is a one-page guide for reviewers.

**Major version 0 is a draft period.** Until sceptre's first release that reads
and writes these files, the format may change and no migration code is
promised. At that release it is frozen as format 1.0 (S10 explains the
two-part version). [Section 10](#10-fixed-now-free-to-grow)
says which parts would be expensive to change after that and which can grow.

Each choice is labelled:

- **Decided**: Eugene has ruled on it.
- **Adopted**: my recommendation, taken up so that the draft is complete.
  Eugene can overturn any of them; [section 8](#8-questions-settled-in-this-draft)
  lists the ones that were open questions and how each was settled.

**How to read this.** Sections 1 and 2 give the whole format in two pages.
Sections 3 to 6 defend each choice in a fixed form: the decision, the
alternatives, why, and what it locks into files. Where this document and
`SPEC.md` disagree, `SPEC.md` is the rule and the disagreement is a mistake to
fix. Evidence is tagged:

- *(verified)*: read in the code, with the file named;
- *(measured)*: a test run for this document, listed in [Appendix A](#appendix-a-measurements);
- *(inferred)*: my reasoning, not checked.

**Precedents.** Each decision that has a counterpart in the IGVF pipeline's
output ends with an **IGVF** note saying what IGVF does and why this format
does or does not follow it. [Section 9](#9-compared-with-igvf-and-pertpy) puts
those notes in one table and adds pertpy. The IGVF facts come from
`docs/mudata_schema.md`, `docs/output.md` and the writer scripts in `bin/` on the
`dev` branch of `pinellolab/CRISPR_Pipeline` at commit `ead4db3` (28 October
2026) *(verified by reading them)*. I did not inspect an IGVF output file, and
I did not read the Google Sheet and Doc that repository links to.

Most code facts come from four read-throughs of the codebases (R sceptre on
the `sceptre3` branch and on `main`, the pysceptre fork, public pysceptre, and
the Nextflow pipeline). I re-read the ones the main decisions rest on myself:
the splitting of control-group names in R, the cell-keyed random numbers, R's
stage-status update, the duplicate-gRNA warning, public pysceptre's count type
and compression notes, the fork's use of saved fits, and the pipeline's result
sort. The rest are as the read-throughs reported them.

---

## 1. The format at a glance

```
my_screen.h5mu                        THE DATA. Written once at import. Never modified.
                                      Modalities: grna, plus one per thing measured (rna, ...).
                                      Carries its own format version and a dataset id.

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
data file is read on its own and the two parts can change independently.

## 2. The eight principles

Every choice below follows from these. Where a choice needed a trade-off, the
earlier principle wins.

| # | Principle | Consequence |
|---|---|---|
| 1 | **The same files from R and Python.** | Only encodings that both languages read *and* write correctly today are allowed. That subset was measured, not assumed (Appendix A). |
| 2 | **Data and analysis are separate.** One dataset, many analyses. | Nothing that depends on an analysis choice is in the data file. |
| 3 | **Every file is written once.** No stage edits a file an earlier stage wrote. | Stages add files. Rerunning a stage deletes its files and those of later stages. Parallel jobs write part files. |
| 4 | **Identifiers, never positions.** | Tables survive subsetting, reordering and being written in pieces. There is no 0-versus-1 question. |
| 5 | **Existing conventions before invented ones.** | AnnData, MuData and Parquet encodings are used as specified. sceptre adds names and rules, not encodings. |
| 6 | **Readers ignore what they do not know.** | Columns, keys and files that do not change the meaning of what is already there can be added with a minor version (S10). |
| 7 | **What follows from other stored values is not stored.** | Unless recomputing it needs a full pass over the counts and nearly every analysis needs it. |
| 8 | **Files are harder to change than code.** | How today's code reads data is evidence when the pattern is sound. Where it is not, the format is chosen on its merits and the code change is listed in [section 7](#7-code-changes-this-implies). |

---

## 3. Number and text types

One rule per kind of value, the same in both languages.

| Value | Type in the file | Why |
|---|---|---|
| Counts (`X/data`) | 32-bit signed integer | Counts are integers. One count never approaches 2³¹. R has no unsigned type. *(measured)* 33% smaller than doubles uncompressed. |
| The two index arrays of the sparse matrix (`X/indices`, `X/indptr`) | one signed integer type for both: 32-bit, or 64-bit past the limit | 32-bit while the stored values, the cells and the features each number at most 2,147,483,647; 64-bit otherwise. It is the type scipy chooses, so the standard Python tools read the counts with nothing done to them. See S4. |
| Integers in tables | 32-bit signed integer | R's integer is 32-bit. *(measured)* R's Parquet reader returns 64-bit integers as doubles, and R's MuData reader returns a 64-bit column as all-missing. No value in the format needs more than 32 bits. |
| Real numbers | 64-bit float (double) | p-values need the range. Never 32-bit. |
| True/false in Parquet tables | Parquet boolean | Round-trips both ways, with missing values. *(measured)* |
| True/false in the data file | categorical with levels `FALSE`, `TRUE` | *(measured)* R writes a logical as an integer and reads a Python boolean as text. See S7. |
| Text | UTF-8 | In the data file, AnnData's plain `string-array` only (S7). In Parquet, plain or dictionary-encoded. |
| A list of numbers per row | Parquet list of doubles | *(measured)* Round-trips both ways between the two Parquet libraries. Used for fit coefficients and resampled statistics. |
| Missing values | float: NaN. Categorical: code −1. Parquet: null. | No integer or plain-text column in the data file may have missing values (S7). |
| Sort order | byte order of the UTF-8 text | *(measured)* R's radix sort and Python's byte sort agree; R's default sort follows the locale and does not. |

**Readers are lenient, writers are strict.** A reader accepts any integer or
float type for a number. A writer writes exactly the types above. That way a
file from another tool still opens, and files from our two writers are the same.

**IGVF.** Counts are stored in the smallest unsigned integer type that holds
them, usually 16-bit (`bin/count_matrix_utils.py`); see S4. Effect sizes in
result tables may be rounded to 32-bit floats, with p-values left as doubles
(`bin/compact_result_dtypes.py`). This format keeps every real number a double:
one rule, and no rounding between what a package computed and what it reads back.

---

## 4. The data file

### S1. One MuData `.h5mu` file — *Adopted*

- **Decision.** The dataset is a single HDF5 file in the MuData layout.
- **Alternatives.** A Zarr store; two `.h5ad` files; the `ondisc` `.odm` format; a folder of Parquet and matrix files.
- **Why.** MuData is the scverse standard for one set of cells measured in several ways. Both languages open it with standard libraries *(measured)*. One file is easy to move and to publish. HDF5 allows one writer at a time, which costs nothing for a file written once.
- **Locks in.** HDF5. Zarr suits cloud object storage better; AnnData defines the same encodings for Zarr, so a Zarr variant could be added later without changing anything else in this document. Public pysceptre's roadmap asks for Zarr *(verified, `docs/ROADMAP.md`)*; I did not test R's Zarr support.
- **IGVF.** The same container: one MuData `.h5mu` file, `inference_mudata.h5mu`. The two formats agree here.

### S2. Modalities: `grna`, plus one per thing measured — *Decided (18 November)*

- **Decision.**
  - The gRNA modality is named `grna`. This name is fixed.
  - Every other modality is named for what was measured: `rna` for RNA (recommended), and, for example, `atac` or `protein`. These names are free. A name starts with a letter and uses letters, digits and underscores, so that it can prefix a covariate name (S5).
  - A file has `grna` and at least one other modality. All modalities hold the same cells in the same order.
  - Which modality an analysis tests is the analysis's choice, recorded in `analysis.json` (A1). A reader finds the response modality there, never by a hard-coded name.
  - A modality that an analysis does not name is ignored. A file may carry a hashing modality, for instance.
  - Every modality follows the rules for counts, cell tables and feature tables, whether or not an analysis uses it *(Adopted, 18 November)*. A file is then conforming or not on its own, whatever analysis looks at it, and any modality can become a response later.
- **Alternatives.**
  - Name the modality for its role: `response` (the earlier drafts, and Eugene's ruling of 17 November, which he reopened on 18 November).
  - IGVF's names, `gene` and `guide`.
  - A fixed vocabulary of assay names.
- **Why.**
  - Which modality is the response is an analysis choice. In an ECCITE-seq screen one analysis tests RNA and another tests protein, on the same data file. Naming a modality `response` would put an analysis choice in the data file, against principle 2.
  - `response` cannot name two modalities, and it does not say what a chromatin-only file contains.
  - No vocabulary is settled beyond `rna`: IGVF uses `gene`, `guide` and `hashing`; pertpy's one MuData example uses `rna`, `adt`, `hto` and `gdo` *(verified, pertpy tutorials at commit `4924c25`)*. So the format fixes only `grna` and leaves the rest free.
  - `grna` is kept because sceptre's vocabulary uses it throughout (`grna_id`, `grna_target`).
- **What stays in role words.** Table columns (`response_id`, `grna_id`, `grna_target`) and setting names (`response_fit_method`, `response_n_umis_range`). A pair is a target and a response whatever was measured.
- **In the packages.** Both packages call the response modality `response` in memory today, and users write `response_n_umis` in formulas. The prototype's readers and writers translate between that and the file's names, so no package behaviour changes *(measured: every round trip in Appendix B still passes)*.
- **Locks in.** The name `grna`, and the rule that names are otherwise free.
- **Not in this format.** Testing two modalities in one analysis. Under this naming it is two analyses of one file.
- **IGVF.** Modalities `gene` and `guide`, with `hashing` when present. This format takes IGVF's idea of a file that may carry more than two modalities, and not its names. An importer from IGVF files renames the two modalities.

### S3. Counts: cells × features, compressed by feature, canonical — *Decided (orientation)*

- **Decision.** `X` holds raw counts as an AnnData `csc_matrix` of shape cells × features. Within a feature, cell positions are strictly increasing, with no repeats and no stored zeros.
- **Alternatives.** Compressed by cell (`csr_matrix`, the common default); dense.
- **Why.** One feature is one contiguous read, which is what fitting one gene needs. Public pysceptre's disk-backed reader requires this layout *(verified, `scripts/sceptre_io.py`)*, the fork reads it directly and has to re-orient the other one *(verified, `execution/scverse_inputs.py`)*, and R's backed reader relies on it *(verified, `R/backed_counts.R`)*. The pattern is sound: every stage works gene by gene.
- **Locks in.** The orientation. A file from another tool is usually compressed by cell and is converted once at import.
- **IGVF.** Compressed by cell (`csr_matrix`), which is what the upstream mapping tools and scanpy produce *(verified, `bin/count_matrix_utils.py`)*. To fit a gene's model, the pipeline splits the whole file into per-gene chunk files first (`bin/chunk_mudata_sceptre.py`). Compression by feature makes that split unnecessary, so this format does not follow IGVF here.

### S4. Count type: 32-bit integers; one index type for both index arrays — *Decided by Eugene (index type, 18 November); otherwise Adopted*

- **Decision.** As in section 3.
- **Alternatives.**
  - Doubles: what R's MuData writes by default *(measured)*. 588 MB against 395 MB uncompressed on a synthetic 200,000 × 4,000 matrix *(measured)*, and it admits non-integer "counts".
  - The smallest unsigned type that fits: what public pysceptre does *(verified)*. 298 MB on the same matrix *(measured)*, but the type then differs from file to file and R has no unsigned integers.
- **Why.** A fixed type makes the two writers agree. With gzip the choice barely matters (98 MB against 103 MB, *measured*), so the smaller unsigned types buy little. R can write these types with six lines after the MuData writer *(measured)*.
- **Locks in.** Counts are integers. R sceptre does not check this today *(verified, `R/check_functions.R`)*; the fork does *(verified, `storage/count_identity.py`)*.
- **One level: a file conforms or it does not** *(Decided by Eugene, 18 November)*. A file whose counts or index arrays are in another integer or floating-point type, with whole values, does not conform, and a reader must still read it (S4.4). That is "readers lenient, writers strict" as one behaviour rule. *(measured)* Both readers read the five such conformance files to the same counts as the conforming file.
  - **Tried and dropped: a middle level, "readable".** Until 18 November the validators reported three levels: conforming, readable (not conforming, but a reader must accept it) and not readable. It was added for counts in another number type, as IGVF writes them. It was dropped because no real file from another tool would land in it: an IGVF file also has other modality names, the other matrix layout and assignments as a layer, so it breaks many rules and needs an importer regardless. The level existed for a file that does not occur, and the index-type decision was adding more cases to it. The validators now return rule, place and message, and a file conforms when they return nothing.
- **The index type** *(Decided by Eugene, 18 November)*. Both index arrays of a modality are stored in one signed integer type: 32-bit while the stored values, the cells and the features each number at most 2,147,483,647, and 64-bit otherwise. This is scipy's own rule *(verified in scipy 1.18.1: `get_index_dtype`; a matrix with a dimension of 2³¹ − 1 gets 32-bit arrays and one of 2³¹ gets 64-bit)*.
  - **What it replaced.** Until 18 November the rule was 32-bit cell positions with 64-bit column pointers, always. That rule was invented here; scipy's is the existing convention. *(measured)* Mixed widths made one common scanpy function, `calculate_qc_metrics`, fail on a plain anndata read, because scipy rejects mixed widths in some operations, and they needed a step in our own reader. The case for 64-bit pointers was arithmetic only: the scale test's 699 million nonzeros are under the limit.
  - **Cost.** Past about 2.1 billion nonzeros both arrays are 64-bit, so a cell position costs 8 bytes, not 4: roughly 12 bytes per nonzero against 8, half as much again before compression, on the datasets where size matters most. And the type now depends on the data. It is still determined by the data alone, so R and Python write the same types; both validators compute the type the size requires.
  - **Rejected.** Always 64-bit for both: one fixed type, but every file half as much again larger before compression.
  - **Wider than needed, or two widths, does not conform**, like counts in another type, and a reader still reads it (S4.4). The Python reader keeps its step that brings two widths to one, as leniency for such files only; it does nothing to a conforming file.
  - **What was checked** *(measured, 18 November, `verify/wide_index.*` and the conformance set)*. scanpy's `calculate_qc_metrics` now runs on every example and conformance data file read with anndata alone, for files written by R and by Python. No test file can reach the limit, so the 64-bit branch is tested with the limit lowered, in both writers and both validators: each writer then stores both arrays as 64-bit, the two files hold the same arrays, and each validator passes them under the lowered limit and reports S4.2 and S4.3 under the real one. **64-bit cell positions in R:** the high-MOI example with both arrays 64-bit reads to the same counts through the prototype's reader in memory, through MuData's backed reader (one feature at a time and whole), and through HDF5Array, with no warning beyond those on the 32-bit file.
  - **What was not checked.** Values above 2,147,483,647 in either array. That needs a file of more than 2.1 billion nonzeros, or cells, and none was built. R holds such values as doubles when read with rhdf5, exactly up to 2⁵³; whether MuData's backed reader and HDF5Array do the same at that size is unmeasured. The test above shows the wider type is handled, not the larger values.
- **Flag.** R's in-memory sparse matrix cannot hold more than 2³¹ − 1 nonzeros. Past that size R must read disk-backed. That is a limit of R, not of the file.
- **IGVF.** The smallest unsigned type that fits, chosen from the data's maximum, as public pysceptre does. This is the closest call in this document: IGVF's rule saves about a quarter of the uncompressed size (298 MB against 395 MB above) and nothing meaningful with gzip. This format keeps one fixed type so that the type never depends on the data. Readers here accept IGVF's types as they are. IGVF's code independently works around the same scipy problem with mixed index types (`normalize_sparse_index_dtypes`).

### S5. Cell covariates: shared table for cell facts, modality tables for measurements — *Decided*

- **Decision.** Facts about the cell (batch, donor, anything the user supplies) are in MuData's shared cell table. Measurements of one modality (`n_umis`, `n_nonzero`, `p_mito`) are in that modality's cell table. The model's covariate table is all of them, a modality's columns prefixed with the modality's name and an underscore: `rna_n_umis`, `grna_n_umis`. Formula text stored in an analysis uses these names (T7).
- **Alternatives.** Everything on the response modality (both earlier drafts); everything in the shared table; a role prefix (`response_n_umis`) whatever was measured.
- **Why.** A batch belongs to the cell, not to its RNA. With two measured modalities, "everything on the response" has no answer for where batch goes, and a role prefix cannot name the RNA depth in an analysis of protein. *(measured)* The shared table round-trips in both directions between R and Python.
- **The cost.** Today's covariate names are `response_n_umis` and so on, and users write them in formulas. In a file they become `rna_n_umis`. Most users never type them, because the default formula is built automatically, and a package can keep accepting the old names; the prototype translates in both directions.
- **Covariate names** *(Adopted, 18 November)*. A name starts with a letter and uses letters, digits, underscores and dots. Dots are allowed because they are common: Seurat's defaults, which pertpy's example dataset carries, include `S.Score` and `percent.mito`. Without this rule such a column could sit in a file and never be used in a formula. *(measured)* R takes a dotted name in a formula as it is. pysceptre's formula parser rejects it, so the prototype's Python reader gives the column a dot-free alias for the call; with that, R and Python build bitwise-equal matrices for `~log(rna_n_umis) + S.Score + batch`. The parser change is listed in section 10.
- **Locks in.** The prefix rule and the underscore. MuData's own `rna:n_umis` colon style is not used, because a colon is not valid in a formula.
- **IGVF.** Per-cell measurements are on the modality's cell table, as here, under scanpy-derived names (`total_gene_umis`, `num_expressed_genes`, `total_guide_umis`, `percent_mito`); some are also copied to the shared table. IGVF's `percent_mito` is a percentage and sceptre's `p_mito` is a proportion, so an importer must convert it.

### S6. Which computed columns are stored — *Adopted*

- **Decision.** Every modality's cell table has `n_umis` and `n_nonzero`. `p_mito` is optional, for a modality that has mitochondrial genes. Not stored: each cell's top gRNA (`max_grna`, `max_grna_frac_umis`).
- **Why required, not recommended** *(Adopted, 18 November)*. A writer can always compute the two columns, the default formula needs them, and "recommended" would make every reader handle their absence. The prototype's writers compute them when a modality lacks them. `p_mito` cannot be required: it exists only where mitochondrial genes do.
- **Why these and no others.** They are used by the default formula and by cell QC in nearly every analysis, and recomputing them is a full pass over the counts (principle 7). `p_mito` also records an import decision: which genes counted as mitochondrial. The top gRNA is read by one assignment method only *(verified, `R/assign_grna_functions.R`)* and follows from the gRNA counts, which are small. The prototype rebuilds it and matches R's values exactly *(measured)*.
- **Alternatives.** scanpy's names (`total_counts`, `n_genes_by_counts`). sceptre's names are kept because users write them in formulas *(verified, default formula builder in `R/aux_functions.R`)*.
- **Not checked.** A validator checks that the columns exist and are numbers, not that they equal what the counts give.
- **IGVF.** Stores about twenty per-cell QC columns from scanpy, including logarithms and percentages of other stored columns (`log1p_total_counts`, `pct_counts_in_top_50_genes`). Principle 7 excludes those here. A file may still carry them: readers ignore columns they do not know.

### S7. Column types in the cell and feature tables — *Adopted*

- **Decision.** A column is one of four kinds: double; 32-bit integer with no missing values; categorical; plain text with no missing values. In addition:
  - a covariate that is text or true/false becomes a categorical once, at import;
  - an integer column with missing values is stored as doubles;
  - strings use AnnData's `string-array` encoding, never `nullable-string-array`.
- **Why.** These are the encodings that survive both directions today *(all measured, Appendix A)*:

  | Written by | Column | Read by the other language as |
  |---|---|---|
  | Python (default settings) | any file at all | **R cannot open the file**: current pandas writes identifiers as `nullable-string-array` |
  | Python | boolean | text (`"TRUE"`, `"FALSE"`) |
  | Python | 64-bit integer | all missing, silently |
  | R | logical | 8-bit integer |
  | R | integer with missing values | **Python cannot open the file** |
  | R | text with missing values | empty strings |
  | either | double, 32-bit integer, categorical, plain text | correct |

- **Alternatives.** Allow every AnnData type and fix the libraries. That is the right long-term direction and the fixes belong upstream (see "Outside the format" in section 8). Until they land, a wider set would break principle 1.
- **Factors.** A categorical's stored category order is the level order, and the first category is the baseline. This is AnnData's convention and R's default contrasts; nothing is invented. Because every text covariate is a categorical, R and Python cannot disagree about level order (R sorts by locale, Python by code point).
- **Locks in.** No boolean or 64-bit columns in the data file. Lifting this later is allowed by principle 6.
- **IGVF.** The file is written by Python only and uses booleans (`targeting`, `mt`) and whatever integer width pandas chose. R reads only the count matrices, the MOI and the pairs from it (`bin/inference_sceptre.R`), so the encodings R misreads never matter there. This format's rule is stricter because both languages must read and write every column. IGVF's writers also replace missing strings with a token and store low-cardinality strings as categoricals (`make_h5mu_safe_dataframe`), for the reason S7 gives.

### S8. The gRNA table — *Adopted*

- **Decision.** One row per gRNA in the gRNA modality's feature table: `grna_target` (required), and optionally `target_chr`, `target_start`, `target_end`, the position of the *target*. A non-targeting gRNA has `grna_target` equal to `non-targeting`.
- **Why.** These are facts about the gRNA, known at import. `non-targeting` is sceptre's long-standing convention *(verified, `R/assign_grna_functions.R`)*; a second marker could contradict it.
- **Why `target_chr` and not `chr`.** sceptre's gRNA table calls these columns `chr`, `start` and `end` and uses them as the target's position when it builds cis pairs. The bare names do not say whether they locate the gRNA or its target, and IGVF's files carry both. The file names say which; the packages keep their own names and the prototype translates.
- **Not stored.** `grna_group`. It depends on the analysis's integration strategy *(verified, `R/preprocessing_functions.R`)*, so it is derived (principle 2).
- **A gRNA under two targets is refused.** It cannot be one row, and a writer stops with a message. R accepts such tables today, with a warning on the development branch *(verified, `R/accessors.R`)*, and public pysceptre counts them and keeps the full table separately *(verified, `scripts/sceptre_io.py`)*. What pooling means for such a gRNA is a scientific question first; a separate gRNA-to-target table can be added later without changing anything else (question 3).
- **Several gRNAs on one vector: `vector_id`** *(Decided by Eugene, 18 November)*.
  - *The ruling.* A vector has one target. Several gRNAs delivered together against the same element are in scope. A vector whose gRNAs have different targets is a combinatorial screen and is out of scope: sceptre is not equipped to analyse those, and the format should be future-proof without accommodating every possibility. This is the rule sceptre already enforces, that `vector_id` is nested within `grna_target`.
  - *The column.* An optional text column `vector_id` in the gRNA table, beside `grna_target`. Which gRNA sits on which vector is fixed at cloning and known at import, so it belongs in the data file by the format's dividing rule. sceptre keeps it today as an optional column of its gRNA table under the same name, so this continues a convention.
  - *Counts stay per gRNA.* A writer never sums over a vector's gRNAs. sceptre sums them at import today, which destroys the one piece of information needed to tell a cell carrying two gRNAs of one vector, as intended, from a cell carrying gRNAs of two vectors, a recombination product. Dual-guide libraries report such recombinants at roughly a quarter of reads *(as reported by the side session from Replogle et al. 2022; not re-read by me)*. Summing within a vector also raises sensitivity when per-gRNA counts are low, and a file without `vector_id` cannot reproduce it.
  - *Why a reader must refuse two targets on a vector (S8.7).* It keeps the door open cheaply. A reader must either read a file correctly or detect that it cannot (S10). With the refusal stated now, allowing combinatorial screens later is a minor version change, because today's readers reject those files cleanly. Without it, the same change would be major.
  - *Why a vector may not be named like another vector's gRNA (S8.8).* sceptre's collapse renames `vector_id` to `grna_id`, so a stray match would silently merge two different things. A vector with one gRNA may still be named after it.
  - *What this covers without anything more.* Under one target per vector, barcode-based designs need no new structure: a barcode standing for several gRNAs with one target is a feature with counts and a target, already a row of the gRNA table. *(As reported by the side session from the papers, not re-read by me:)* Adamson et al. 2016 and Norman et al. 2019 identify a multi-guide vector by one barcode; Replogle et al. 2022 count each guide separately by direct capture; CaRPool-seq, CROPseq-multi, Cas12a arrays and VIPerturb-seq each use a barcode per array or per guide. Both schemes fit.
  - *Not added, on purpose:* a vector-to-gRNA table (needed only when a vector carries several targets), targets as lists, a position within the vector, a field for the detection scheme, and a field for the target molecule.
  - *Deferred: the analysis side.* Because the data file holds per-gRNA counts and is written once, sceptre's sum-within-vector has to become a step of an analysis. That implies a setting in the assignment stage's record, and an assignment table that names the unit it assigned, gRNA or vector, as results tables already switch their key by strategy. Nothing in the data file depends on the answer. It is listed in section 10. The prototype's readers do not use `vector_id` when they build a package object.
- **No entry for the target molecule** *(Decided by Eugene, 18 November)*. Target coordinates are genomic positions; a screen that knocks down RNA names its transcript in `grna_target` and may leave coordinates out. Nothing in an analysis depends on the molecule: counts, assignments, targets and tests are the same, and the coordinates are used only when the user builds cis pairs. The optional `perturbation_type` (S10) records which kind of screen it was, and its knockdown value is the RNA-targeting case.
- **The name `grna` is a slight wart.** The modality holds whatever identifies the perturbation, which for some designs is a barcode and could in principle be something other than a gRNA. It is kept because sceptre's vocabulary uses it throughout and another name buys nothing within this scope.
- **Positions are 1-based with both ends included** *(Adopted, 18 November, on the recommendation of a side session; Eugene can overturn it)*. `target_start` is the first base of the target and `target_end` the last. One fixed convention for the whole format, with no entry saying which is in use.
  - *sceptre's own data are already in this convention* *(verified by the side session on the `sceptre3` build; I re-checked the first rows in this folder's example)*. In sceptre's high-MOI example every target has `end − start = 49`, and a target's second gRNA starts at the first one's end plus 1 (23567064 to 23567113, then 23567114 to 23567163). Adjacent 50-base tiles are consistent only with closed intervals. The gene positions sceptre pairs them with come from a GTF file, which is 1-based with both ends included.
  - *Nothing has to be converted.* Files written from existing sceptre objects are already right, and target and gene positions are in one system. It is also the convention of Ensembl and of Bioconductor's `GRanges`, which the planned genome-build check would compare.
  - *Why it was never pinned down.* sceptre builds cis pairs from the midpoint of a target and a distance threshold of hundreds of kilobases, so being one base off changes nothing. It matters for exchanging files, not for sceptre's results.
  - *The cost.* Python tools and BED files usually count from 0 and exclude the end. A fixed, stated rule and one sentence ("add 1 to a BED start at import") handle it.
  - *Checked.* Both validators report a start below 1, an end before its start, and a gRNA with some of the three position columns and not all.
  - *Chromosome names* are text; the format does not choose between `chr1` and `1`. sceptre's own data use the `chr` prefix.
  - *IGVF* does not fix a convention: its schema has a per-guide column `coordinate_system`, with zero-based half-open and one-based inclusive given as examples. That passes the question to every reader, which this format avoids. crisprVerse uses `GRanges`, but its guide coordinate is a single position at the PAM, a different thing from the targeted element *(as reported by the side session)*.
- **IGVF.** Much richer, and better than this format today. Each gRNA has a `targeting` flag, a `type` (targeting, non-targeting, safe-targeting, negative control), its `spacer` sequence, its own position (`guide_chr`, `guide_start`, `guide_end`, `strand`, `pam`), its target's name and position (`intended_target_name`, `intended_target_chr`, `intended_target_start`, `intended_target_end`) and a collision-safe target key; source columns are passed through. This format takes the distinction between gRNA position and target position now, through the column names. The rest are optional columns that can be added without a format change and are listed in [section 10](#10-fixed-now-free-to-grow). A file may already carry them.

### S9. Identifiers — *Adopted*

- **Decision.** Every cell and every feature has a unique text identifier. Feature identifiers are unique across all modalities of a file, so no table needs a modality column.
- **Why.** Principle 4. Cell identifiers are also load-bearing for results: R's current random-number protocol keys each cell's draws on its identifier *(verified, `R/rng_v2.R`)*.
- **Screens without barcodes.** The writer assigns the positions `"1"`, `"2"`, ... as identifiers and stores them, as R does in memory today *(verified)*. Because identifiers affect results, they are fixed in the file, so both languages see the same ones.
- **No line feeds.** An identifier must not contain a line feed, because the hash that links an analysis to its data joins identifiers with line feeds (A1).
- **IGVF.** The same: cells, genes and guides are named by identifier (the cell barcode with a measurement-set suffix, the Ensembl gene identifier without its version, the guide name), and tables are keyed by them.

### S10. File-level metadata: a native record, not JSON text — *Adopted*

- **Decision.** `uns/sceptre` is an AnnData dictionary with scalar text entries: `format` (the version, such as `0.1`), `moi` (`high` or `low`), and `id` (a random identifier assigned when the file is written, used in A1).
- **A two-part version** *(Decided by Eugene, 18 November)*. The version is `MAJOR.MINOR`. A minor increase adds something an older reader can safely ignore; a major increase is a change after which an older reader could misread a file. A reader reads any file of the major version it knows and refuses another.
  - *Why.* With one number, every addition either goes unrecorded or forces every reader to refuse the file. Two parts let a file say what it uses while old readers keep working. The stamp is in every file, so this had to be settled before the freeze.
  - *Precedent.* Apache Arrow versions its format this way and states the rule as: an older reader must either read newer data or detect that it cannot *(verified by the side session that read Arrow's `Versioning.rst` and `Changing.rst`; I did not re-read them)*.
  - *What it corrects.* Principle 6, "readers ignore what they do not know", is safe only for additions that do not change the meaning of what was already there. Section 10 marks each planned addition as minor or major.
  - *Stored as text*, the same in HDF5 and in JSON, and compared as two integers, never as text.
- **Alternatives.** One JSON string (the first draft of this format); HDF5 attributes.
- **Why.** Principle 5. *(measured)* Both languages write and read a nested record with true scalars; R needs one attribute set on each value to write a scalar. A JSON string needed a low-level HDF5 write in R and is opaque to HDF5 tools.
- **No list of text columns.** The first draft kept one, to mark text stored as a categorical. Under S7 every covariate that was text *is* a categorical, and the few text columns of the feature tables are named in `SPEC.md`.
- **The version here is the data file's own.** An analysis has a separate one in `analysis.json`.
- **What kind of screen it was: `perturbation_type`, `organism`, `genome_build`** *(Decided by Eugene, 18 November)*. Three optional text entries. No reader's behaviour depends on them, and a file without them conforms.
  - *Why they are in now, though nothing reads them.* "Nothing reads it, so it can be added later" is the right test for a value that follows from other stored values (principle 7). It is the wrong test for provenance. What kind of CRISPR perturbation a screen used is known at import, cannot be worked out from anything in the file, and is lost for good if not written down. It is also needed to read a results table at all: a log2 fold change of −1 is the expected result under CRISPRi and a surprising one under CRISPRa. The dividing rule already settles where it goes: the data file holds what is known at import.
  - *The vocabulary* for `perturbation_type` is crisprVerse's five kinds: `CRISPRko`, `CRISPRi`, `CRISPRa`, `CRISPRbe`, `CRISPRkd` *(verified by the side session in the crisprDesign README; that `CRISPRkd` means knockdown of RNA with an RNA-targeting enzyme is from memory and should be confirmed)*. CELLxGENE's schema, the EMBL-EBI Perturbation Catalogue and scPerturb record the same concept *(as reported by the side session)*. A new kind is a minor addition.
  - *`genome_build`* is the more pressing of the other two: a coordinate without a build is only half-interpretable. It does not settle whether coordinates are 0-based or 1-based, which is still open (S8).
  - *It replaces a separate entry for the target molecule.* The knockdown kind is the RNA-targeting case, so one entry covers both distinctions.
- **IGVF.** `moi` (`high` or `low`, the same two words) and `capture_method` are in the guide modality's `uns`, each as a one-element string array. There is no format version and no dataset identifier; a file is identified by a SHA-256 of all its bytes in a separate manifest (A1). This format adds the version and the identifier because files outlive the code that wrote them.

### S11. Compression — *Adopted: writer's choice*

- **Decision.** The count arrays may be uncompressed or gzip-compressed. No other filter.
- **Why.** HDF5 decompresses transparently, so this is not locked in (principle 8). *(measured)* gzip made a synthetic matrix 5.7 times smaller and 15 times slower to read (588 MB and 0.07 s, against 103 MB and 1.2 s). Other filters need plugins that one language or the other lacks. *(inferred; not tested)*
- **Default for our writers.** Uncompressed, for files that are analysed; gzip for archiving. Public pysceptre chose the same for the same reason *(verified)*.

### S12. What is not in the data file — *Decided*

Nothing an analysis produces: no assignments, no kept-or-removed flag, no
gRNA groups, no results. Public pysceptre keeps a kept-or-removed flag in the
data file and then needs a script that rewrites the whole file to work
disk-backed on the kept cells *(verified, `scripts/write_in_use_h5mu.py`)*.
That is the cost this rule avoids.

**IGVF.** This is the largest difference between the two formats. IGVF's one
file holds the data, the assignments (a layer), the pairs to test and the
results (tables in `uns`), and only the cells that passed QC. That gives one
self-contained file to upload to a portal, which is IGVF's purpose. For one
dataset with several analyses it costs what its own code documents
*(verified)*:

- a dedicated script patches results into the file without rewriting the
  matrices, and records a silent-corruption bug it had to work around
  (`bin/mudata_uns_io.py`);
- sceptre's own write of the file is switched off by default, because it was
  "tens of gigabytes of writes per run" (`bin/inference_sceptre.R`);
- parallel work means splitting the file into chunk files;
- the newer "checkpoint" mode keeps an immutable parent file with a hash and
  writes each continuation to a new folder (`docs/preinference_checkpoint.md`),
  which is a step toward the separation this format starts from.

A single-file export for a portal can be produced from this format; the
reverse, several analyses sharing one data file, cannot be recovered from a
merged file.

---

## 5. The analysis folder

### A1. Linking an analysis to its dataset and its response modality — *Adopted*

- **Decision.** `analysis.json` records:
  - `format`, the analysis's format version (S10);
  - `response_modality`, the name of the modality this analysis tests (S2). It is always written, even when the file has only one candidate;
  - `screen`: the dataset's `id` (S10), its numbers of cells, responses and gRNAs, a SHA-256 of each ordered identifier list, and a relative path as a hint.

  A reader refuses a dataset whose `id` or hashes differ, or that has no modality of that name.
- **Why the response modality is here and not in the parameters record.** The link's response count and hash are those of the response modality, so the two belong together, and nothing in an analysis survives a change of response modality.
- **Alternatives for the link.**
  - A path only: breaks when either side is moved.
  - Size and checksum of the data file (the fork's newest branch, *verified, `io/analysis_files.py` on `sceptre3-dsfix`*): breaks if the file is recompressed, and a checksum of a 100 GB file is slow.
  - Requiring the analysis to live inside the dataset's folder: impossible for a shared, read-only dataset.
- **Why two checks.** The `id` survives moving and copying and changes when the dataset is re-imported. The identifier hashes catch a file that was subset by a generic tool that kept the `id`. *(measured)* Hashing 10 million cell identifiers takes 0.3 s in Python and 2.4 s in R.
- **Not caught.** Counts edited in place with the identifiers and `id` unchanged.
- **Defaults are for packages.** A package may choose the response modality for the user when the file has one candidate. The file is always explicit.
- **Convention.** Analyses live in `analyses/<name>/` beside the data file unless there is a reason not to.
- **Analyses far from their data** *(Decided by Eugene, 18 November; details Adopted)*. A data file may sit read-only on a server while its analyses live on a laptop or elsewhere on the server. The link already allows this: it is the `id` and the hashes, a reader is told where the data file is and never follows `screen.path`, the data file knows nothing about its analyses, and an analysis folder can be read for its results with no data file at all. Four things were tightened:
  - *The hint is a relative path or a URL, never an absolute local path.* An absolute path names a user and a server's layout, and means nothing on another machine. The schema enforces the form and the validators report it (A1.4).
  - *A writer records a path only when the data file is near*: in the analysis folder's parent or grandparent, or below one of them. Farther away, the relative path is a chain of `..` that names directories and helps nobody, so the writer records nothing.
  - *A dataset may have a `name`*, a label for people, stored in the data file and copied into the link. An error can then say which dataset an analysis was run on instead of quoting a random `id`. The name is not part of the identity check.
  - *A dataset may have a `source`*: a URL, DOI or accession, such as a GEO accession for a published screen.
- **No dataset version** *(Decided by Eugene, 18 November)*: "too much bookkeeping". The `id` changes when a dataset is imported again, which is the protection an analysis needs.
- **Finding the data is package behaviour, not format.** A package can resolve a dataset `id` to a location through a search path or a small table in the user's own configuration. Nothing specific to one machine goes in `analysis.json`. If a package ever follows the hint by itself, it should do so only for locations the user has configured: a folder from an untrusted source could otherwise point a reader at any file or URL. None of this is built.
- **Precedents** *(read by the side session of 18 November; I did not re-read them)*. The common pattern is to identify data by a name or a hash, keep its location as a separate, replaceable mapping, allow only relative paths or URLs, and carry a label for people.
  - *Frictionless Data Package*: a file reference is a relative POSIX path or an HTTP(S) URL; absolute paths are disallowed, as a security risk and as rarely useful.
  - *BIDS*: references name a dataset, not a location; a separate table maps each name to a local or remote address; source datasets are recorded by URL, DOI and version.
  - *DVC*: data is identified by a content hash with a relative path; the physical location is a named remote configured outside the pointer file.
- **IGVF.** Its checkpoint manifest links a continuation to its parent by a SHA-256 of the whole parent file, plus cell and feature counts, the pipeline commit and a hash of the gene annotation file (`bin/preinference_checkpoint.py`). That catches an in-place edit of the counts, which this format's link does not. It costs a full read of the file and breaks if the file is recompressed. An optional whole-file checksum could be added here for those who want that guarantee ([section 10](#10-fixed-now-free-to-grow)).

### A2. One record per stage, written last — *Adopted (changed from my earlier advice)*

- **Decision.** Each completed stage has a folder named for it, holding its tables and `record.json`, which is written after the tables. The record holds the settings that stage used, small summaries, provenance (package, version, language, time), and the list of files the stage wrote with their row counts.
- **Alternative.** One `analysis.json` holding everything, rewritten by each stage. I recommended this earlier.
- **Why I changed.**
  - It keeps principle 3: no file is ever rewritten. The single record is the one file every stage would modify.
  - sceptre's settings arrive stage by stage (`assign_grnas(method = ...)`, `run_qc(...)`) *(verified, both codebases)*, so they have a natural home.
  - The record is a commit marker. A stage whose record is missing did not finish, whatever files exist.
  - The calibration check, power check and discovery analysis do not depend on one another's results *(verified, same rank in both codebases)*. A pipeline could run them at once, and they would race to rewrite a single record.
- **Cost.** Up to six small files instead of one. Reading all settings means reading several.
- **A folder per stage** *(Decided by Eugene, 18 November)*. Earlier drafts kept the records in `stages/` and the tables elsewhere: some grouped by kind (`pairs/`, `pair_qc/`, `results/`, `fits/`), some loose at the top, with stage ownership recorded only inside the records. Folder and file names are expensive to change later, so this was settled now.
  - *It mixed two organising ideas*, and used two vocabularies for one three-way split: pair sets were named by pair type and results by stage, so `pairs/positive_control` went with `results/power`. Some names said what they were (`pair_qc`) and some did not (`cells.parquet`).
  - *Rules existed only to police the layout.* A record could list another stage's table, and clearing a stage meant finding its files through its record. With a folder per stage, a record names files in its own folder and nothing else (the schema enforces it), and clearing a stage is deleting its folder. Built and checked: the rule against listing another stage's table became a property of the schema; the rule against unlisted tables remains, now confined to one folder, and still catches a stage that did not finish.
  - *Parallel testing stages write into separate folders*, which suits a pipeline.
  - *Cost.* Results are in three folders, not one `results/` folder.
  - *One name chosen with care:* the QC verdicts are `discovery_pair_qc.parquet` and `positive_control_pair_qc.parquet`, not `discovery_pairs.parquet`, so that a file copied out of its folder is not mistaken for the pairs-to-test table of `parameters/`.
- **Stage names.** `parameters`, `assignment`, `qc`, `calibration`, `power`, `discovery`: neutral words, not R function names.
- **`provenance.language` is free text** *(Adopted, 18 November)*. The two libraries here write `R` and `Python`. A list of allowed values would shut out a third implementation. A saved fit is reused only between languages that share the fitting method (T6), so a fit from an unknown language is simply not reused.
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

- **Names inside settings.** A formula names covariates as the file does (`rna_n_umis`, S5). Setting names themselves keep role words (`response_n_umis_range`, `response_fit_method`).
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
- **Why.** What the October checksum file caught in practice was a truncated or half-written folder: sizes were checked on every read, checksums only on request *(verified, fork `io/analysis_files.py`)*. A record written last, with row counts that a reader checks against each table's own footer, catches the same failures without hashing large files on every write.
- **Left out on purpose.** Detection of silent corruption with the size unchanged. A stage record may carry an optional SHA-256 per file for anyone who wants it (principle 6).

### A6. Writing safely — *Adopted*

A writer writes each file under a temporary name and renames it into place.
Parallel jobs each write their own part file into a table's folder. Public
pysceptre's shard scripts already work this way *(verified,
`scripts/run_trans_background.py`)*. In the Nextflow pipeline each stage writes
a new object file in its own work folder and hands it on, so no job appears to
modify a file another job reads *(inferred from `main.nf` and the `bin/` scripts)*.

---

## 6. The tables

### T0. Rules for every table — *Decided*

1. Parquet. Types as in section 3.
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
- **Alternatives.** A 0/1 matrix in the data file (public pysceptre); a matrix in its own `.h5ad`; lists inside the settings record (October).
- **Why, independent of today's code.** Assignments are tiny next to the counts. Target unions and control groups need many gRNAs at once. Parallel jobs can each write a part. Another tool can hand in a plain two-column table.
- **Supporting evidence from the code.** All four codebases load every assignment: the pipeline merges per-job lists *(verified)*, the fork rebuilds everything from the per-gRNA lists *(verified, `pipeline/r_state.py`)*, and R holds them in the object. Loading everything is sound at realistic sizes: the side session measured a 10-million-row table at 7.8 s in Python and 10.1 s in R, 76 MB on disk (its measurement, not repeated here).
- **Derived on read.** Per-target unions, per-control lists, and gRNAs per cell.
- **Room to grow.** A score column (a posterior probability, a p-value) can be added under rule 6.
- **IGVF.** A 0/1 cells × guides matrix in the data file, `guide.layers['guide_assignment']`, with an optional second layer of posterior probabilities. pertpy does the same (`layers["assigned_guides"]`). The matrix and the table hold the same information. This format uses the table because the assignment belongs to an analysis, not to the data (S12), and because parallel jobs can each write a part. Since both precedents use the matrix, an export to that shape is worth providing; it is package code, not format.

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
- **The `grna_id` column is a rule, not an option** *(Adopted, 18 November)*. Whether a pair-QC or results table has the column follows from the strategy and the stage: never under `union`; always under `singleton`; under `bonferroni` in pair QC and calibration, and not in power or discovery, whose results sceptre reports per target. A reader can therefore know a table's key without opening it. The bonferroni irregularity is sceptre's, and the rule changes when sceptre's output does (section 10).
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
- **`log_2_fold_change` is stored only where there is no `fold_change`** *(Decided, 18 November)*.
  - It follows from `fold_change` (principle 7). A reader computes it as `log(fold_change) / log(2)`, which is sceptre's own expression *(verified, `R/aux_functions.R`)*. The expression matters: `log2()` gives a different last bit on about 30% of rows *(measured)*.
  - *(measured)* On all 17,219 result rows of the examples, the derived value is identical to the one sceptre computed, in R and in Python.
  - Under the bonferroni strategy, sceptre's discovery and power results have no `fold_change` and no `se_fold_change`, only `log_2_fold_change` *(verified, `R/s4_analysis_functs_2.R`)*. There it is the stored column. The reverse derivation is not exact, so it cannot be dropped.
  - A reader uses a stored `log_2_fold_change` when it is present and derives it otherwise. A writer leaves it out only when it equals the derivation.
  - This changes nothing a user sees: the packages' result tables keep the column.
- **Sorted by key, not by p-value.** Sorting by key lets parts be written independently and lets a reader find one target's rows. Packages sort by p-value when they hand the table to a user. This also removes a dependence on job order among tied rows *(inferred from `bin/process_association_analysis_results.R`)*.
- **Size.** *(measured, synthetic)* 33 bytes per row, so about 11 GB at 339 million pairs. Reading two columns of 20 million rows takes 0.05 s in Python and 0.8 s in R; as 20 part files, 1.7 s in R.
- **Resampled statistics** (the largest output setting) go in `<stage>/resamples.parquet`: the key columns and one list column, `z_null`. Today packages hold them as thousands of columns `z_null_1`, `z_null_2`, ..., padded with missing values *(verified, both codebases)*. *(measured)* Asked to keep them for the example's 1,406 pairs, the R build used here failed outright; the example keeps them for 40 pairs.
- **The `significant` column is not stored.** It depends on the whole table, so parallel jobs could not write it. The stage record stores `significance_threshold`, the largest p-value called significant (or `null` if none was), and a reader computes `significant` as `p_value <= threshold`. This is exact for every correction sceptre offers, because each calls significant exactly the pairs up to some p-value. *(measured)* It reproduces R's column for every pair of every example.
- **Very small p-values.** A p-value below about 1e-308 becomes zero in a double. A `log_p_value` column can be added when a method needs it (rule 6); none does today.
- **IGVF.** Two tables per scope, per guide and per element, which correspond to this format's singleton and union strategies. They differ from this format in four ways:
  - *Columns carry the method's name* (`sceptre_p_value`, `sceptre_log2_fc`, `sceptre_fc_se`), because one table holds two methods' results. Here an analysis is one method's, so the names are plain.
  - *Annotations are copied into every row* (guide sequence, coordinates, gene symbol, number of perturbed cells). Here a row has keys and results, and annotations are joined from the data file.
  - *Adjusted p-values are stored* (`sceptre_q_value`) along with `sceptre_negLog10p`. Both follow from the p-values. Here neither is stored; the threshold plays the role of the adjusted column.
  - *Each table is stored three times*: in the file's `uns`, as Parquet and as TSV. Here once, as Parquet.

  IGVF's working table format is Parquet with zstd, the same as here. Its per-row count of perturbed cells is a useful diagnostic that these tables do not carry; it is listed in [section 10](#10-fixed-now-free-to-grow).

### T6. The fit cache: `fits/response/` — *Decided (placement; where the identity lives), Adopted (details)*

- **What it is.** A cache of fitted response models. A reader may ignore it, and deleting it changes no result beyond numerical tolerance.
- **Columns.** `response_id`, `coefficients` (a list of doubles, in the order of the design's terms), `theta`. Key: `response_id`.
- **Alternatives for the columns.**
  - One row per coefficient, with the dispersion as the row whose term is missing (the October layout): a trick.
  - One column per term: a covariate named `theta` would collide, and the columns change with the formula.
- **Outside the stages** *(Decided by Eugene, 18 November)*. A response fit depends on the data, the formula, the fitting method, the control group and which cells passed QC. All of that is fixed once the parameters and QC stages have run. Which testing stage happens to compute a fit first is an accident. Earlier drafts kept fits in the folder of the stage that first made them, which would:
  - delete still-valid fits when that stage was rerun;
  - stop stages running in parallel from sharing them;
  - name them after a stage that does not own them.

  (The earlier prototype wrote every fit into the first testing stage's file, an artefact of writing a finished object in one go.)
- **How the cache works.**
  - It is a folder of part files. Each job adds a part of its own and never modifies another's, which keeps "every file is written once" and uses the rule that a table may be a folder of parts.
  - Two parts may hold the same response, when two parallel jobs fitted it. A reader takes either. So the rule that a key is unique does not hold across the parts of this table, and the specification says so.
  - Fits are listed in no stage's record. What makes them valid is their **identity**: the dataset `id`, the formula, a hash of the identifiers of the cells in use, the fitting method and the control group. Each part carries the identity of its own fits, inside the file (next point).
  - The cache is deleted when the QC stage or an earlier one is rerun, and kept when a testing stage is rerun. *(measured)* Both libraries do this, and a conformance case covers it. Since the identity moved into the parts, deleting is housekeeping (SHOULD): correctness does not depend on it. It is kept because a plain read of the whole folder would otherwise mix old parts with current ones.
- **The identity is inside each part** *(Decided by Eugene, 18 November)*. It is JSON text under the key `sceptre_fit_identity` in the Parquet file's own key-value metadata. There is no `identity.json`, neither in the folder nor beside it; the folder holds part files and nothing else.
  - **Why.** The identity exists to stop a fit being reused when it was computed for other cells, another formula or another method. If that check passes wrongly, the result is wrong p-values with no error. With a separate file, the check reads one file and trusts that it describes the others, and they can come apart: QC rerun with the clearing interrupted; a job adding parts without checking the file; two parallel jobs both finding a stale cache and both clearing and rewriting it. Each is unlikely, and the failure is silent. With the identity inside each part there is nothing to keep in step, a stale part is harmless, and no file is shared between writers, which fits "every file is written once". It also removes the earlier rule that the first writer writes the record and later writers verify it.
  - **Why the other tables do not need this.** A stage finishes. Its record, written last, lists its files with their row counts and so marks the stage as complete, and rerunning the stage deletes the folder first. A record written last cannot work for something open-ended that jobs keep adding to. That structural difference is the whole reason the cache is treated differently.
  - **What it replaced.** `fits/response/identity.json` (drafts up to 18 November), and the recommendation, mine and a side session's, to move that file beside the folder as `fits/response.json`.
  - **Evidence** *(measured by a side session on 18 November; the first and third lines are repeated in the suite on every run)*:
    - A JSON file inside a folder of parts: `pyarrow.dataset`, `pyarrow.parquet.read_table`, `pandas.read_parquet`, R's `arrow::open_dataset`, `polars.scan_parquet` on the folder and duckdb's `read_parquet('folder/*')` all fail. Only an explicit `*.parquet` pattern works.
    - The same file named `_identity.json`: the Arrow family ignores it and works; polars and duckdb still fail. So an underscore is not enough.
    - Parts only: everything works. *(measured here)* All six readers open the cache folder of every example directly.
    - Both directions across languages, with the libraries already in use: pyarrow writes the key and nanoparquet reads it; nanoparquet writes it and pyarrow reads it. Arrow adds a key of its own, `ARROW:schema`, which is ignored. duckdb and R's arrow show the key too.
  - **Costs.**
    - It cannot be read with a text editor, so the validators show it, and `fit_parts()` in both libraries returns it.
    - It is a mechanism the format uses nowhere else.
    - A tool that rewrites a part may drop the metadata. The part is then ignored, which is the safe direction. *(measured here)* This happened while building the conformance set: pyarrow's `read_table` exposes only the metadata that Arrow's own schema record holds, so a part written by R, read and written back with pyarrow, lost the key until the builder copied it from the file's metadata. A reader must take the key from the file's metadata, not from the table.
  - **Parts that are not the analysis's are not reported.** A part with another identity is legal and unused, so the validators say nothing about it. They do report a part with no identity, or one that does not conform, and anything in the folder that is not a part file.
  - **One limit of reading the folder alone.** Four of the five identity fields can be read from the analysis folder. The fifth, the hash of the cells in use, needs the data file's order of cells. So a reader given only the folder returns the fits made for the analysis's parameters when all such parts agree on the cells, and nothing when they do not (parts left by an earlier QC beside current ones). A reader that has the data file gives the full identity and gets exactly its parts. Both libraries work this way.
  - *(not verified)* GeoParquet is believed to keep its metadata the same way, as JSON under one key. A side session recalled this and did not confirm it, and neither did I.
- **When a cached fit is reused.** Only from a part whose identity equals, field by field, what the reader would use itself, and only if the fitting method is one the reader's language shares. Identities are compared as parsed values, with white space in the formula ignored. The identity does not name the language that made a fit: a method name means the same estimator in every language that lists it, and parts may come from several languages.
- **It is the part of the layout a new reader is most likely to misunderstand**, because it looks like a seventh stage and is not: no record, it fills gradually, and it is cleared at a different moment. The specification says so at the head of its section.
- **Room for gRNA fits** *(Eugene, 18 November: these will be wanted eventually, and this version must not step on them)*. The folder is named for its kind, `fits/response/`, and `fits/grna/` is reserved. The two would not be one table with a kind column: they have different keys (a response, against a target or a gRNA depending on the strategy), different columns (no dispersion for a gRNA fit) and different dependencies. A gRNA fit would also depend on the assignments, so its cache would be cleared when assignment is rerun *(inferred by the side session from how conditional resampling works, not read in the code; to verify before the folder is specified)*.
- **`response` is unambiguous** because an analysis names one response modality (A1). The folder is not named for the modality (`fits/rna/`), which would put a user-chosen name beside the reserved `grna`.
- **Fitting methods are shared, named definitions** *(Decided in principle by Eugene, 17 November: "We want these methods to be consistent across R and Python.")*
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

  So there are two families: `sceptre` and `pysceptre` agree, and `batched` (the current default in both languages) agrees with itself across languages but treats some responses differently. glmGamPoi was not installed and was not measured. Which estimator is the default is a scientific decision, outside the format (section 8).
- **Not adopted.** Public pysceptre's cache also keys on library versions and CPU features, for bit-exact reuse *(verified, `pipeline/fingerprint.py`)*. That is stricter than a cross-language cache needs.

### T7. The design matrix — *Decided*

- **Decision (Eugene, 17 November).** Not stored. The parameters record holds the formula text and a hash of the matrix; each language rebuilds the matrix. Only additive formulas of covariate names, `log(x)` and `log(x + 1)` are allowed, and a writer refuses others. No stored-matrix fallback.
- **Names.** The formula text and the recorded term names use the file's covariate names (`log(rna_n_umis)`, S5). The hash is of the matrix's values, so it does not depend on names.
- **Evidence.** *(measured)* On the example, R's `model.matrix()` and pysceptre's `model_matrix()` give bitwise-equal matrices and the same hash.
- **What a reader checks.** Both the rebuilt matrix's term names and its hash. The names (`batchbatch_2`) encode the factor levels, the baseline and the term order, which are the ways two languages can really disagree; the hash catches everything else.
- **IGVF.** The design is recorded as a list of covariate options in the file's `uns`, and each method builds its own design from them. There is no check that two methods, or two languages, built the same matrix.
- **Known limit.** A bit-exact hash can fail on a correct matrix if a logarithm differs in its last bit between platforms or math libraries. On this machine R and Python agree to the last bit; I could not test another (question 4).
- **Restriction to state plainly.** R sceptre accepts any formula today. Under this rule an analysis with an interaction or a spline cannot be saved.

---

## 7. Code changes this implies

The format is chosen on its merits; these are the places where package code
would follow. The prototype's readers and writers already do the file side of
each.

| Where | Change | Why |
|---|---|---|
| R writer | Write counts as 32-bit integers and both index arrays in the index type; write file metadata as a native record | S4, S10 |
| Python writer | Turn off `nullable-string-array`; write 32-bit counts (scipy already gives the index type) | S7, S4 |
| Python reader | Bring a sparse matrix's two index arrays to one integer type after reading | S4 |
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
| Both | Translate between the packages' `response_` covariate names and the file's modality names, as the prototype does; later, accept the modality names directly | S2, S5 |
| Both | Translate the gRNA table's `chr`, `start`, `end` to and from `target_chr`, `target_start`, `target_end` | S8 |
| Both | Derive `log_2_fold_change` on read with `log(x) / log(2)` | T5 |

## 8. Questions settled in this draft

These were open when the review began. Each is now settled as shown and
implemented. Any of them can be reopened; the last column says what reopening
would cost.

| # | Question | Settled as | To change it later |
|---|---|---|---|
| 1 | Per-stage records or one record? (A2) | Per-stage records, written last | A format change: it moves where settings live |
| 2 | Count type (S4) | 32-bit integer counts; both index arrays in one type, 32-bit until the size needs 64 (row 29) | Readers already accept any numeric type; only writers change |
| 3 | A gRNA under two targets (S8) | Refused, with a message | Add a gRNA-to-target table; nothing else changes |
| 4 | How to check the rebuilt design matrix (T7) | Term names and the bit-exact hash Eugene decided on | Relax the hash to a warning if a second platform shows last-bit differences |
| 5 | The `significant` column (T5) | Not stored; one threshold in the stage record | Storing the column too is an added column |
| 6 | p-values too small for a double (T5) | Nothing now | An added `log_p_value` column |
| 7 | Stage names (A2) | `parameters`, `assignment`, `qc`, `calibration`, `power`, `discovery` | A format change: file names |
| 8 | Default compression of counts (S11) | Uncompressed; gzip allowed | Nothing: readers do not care |
| 9 | Several measured modalities (S2) | A file may hold several; an analysis names one as its response. Testing two in one analysis is not in this format | A format addition when a dataset needs it |
| 10 | Identifiers for cells without barcodes (S9) | Positions as text, written into the file | Nothing for existing files |
| 11 | Cells flagged by the assignment method (T1b) | An optional table of cell identifiers | Add a reason column |
| 12 | Zarr (S1) | HDF5 only | A Zarr variant with the same encodings |
| 13 | Modality names (S2) | `grna` fixed; others named for what was measured, `rna` recommended; the analysis names its response modality (Eugene, 18 November) | A format change: names |
| 14 | The gRNA table's coordinate columns (S8) | `target_chr`, `target_start`, `target_end` | A format change: names |
| 15 | `fold_change` and `log_2_fold_change` both stored (T5) | The log2 column only where there is no `fold_change` | Storing it everywhere again is an added column |
| 16 | Modalities an analysis does not use (S2) | Allowed and ignored; they follow the same rules as the others | Nothing |
| 17 | The form of the version (S10) | Two parts, `MAJOR.MINOR` (Eugene, 18 November) | A format change: it is the stamp itself |
| 18 | `n_umis` and `n_nonzero` (S6) | Required for every modality | Loosening it is safe for old files |
| 19 | Covariate names (S5) | Letters, digits, underscores, dots | Widening it is safe for old files |
| 20 | When tables have `grna_id` (T3, T5) | Fixed by strategy and stage | Follows sceptre's output |
| 21 | Counts in another number type (S4) | Not conforming; a reader still reads them (S4.4). See row 31 | Nothing |
| 22 | The hint to the data file (A1) | A relative path or a URL, never absolute; recorded only when the file is near; an optional dataset name and source; no dataset version (Eugene, 18 November) | Loosening the form is safe for old files |
| 23 | Several gRNAs on one vector (S8) | An optional `vector_id`; one target per vector, and readers refuse otherwise; counts per gRNA (Eugene, 18 November) | Allowing several targets later is a minor change |
| 24 | The target molecule (S8) | No entry of its own; `perturbation_type` covers it (Eugene, 18 November) | A minor addition |
| 28 | The coordinate convention (S8) | 1-based, both ends included; all three position columns or none | A format change: it changes what stored numbers mean |
| 29 | The index type of the counts (S4) | Both index arrays in one width: 32-bit until a dataset needs 64, scipy's rule (Eugene, 18 November). Replaces 32-bit positions with 64-bit pointers | A format change: types |
| 30 | Where the fit cache's identity lives (T6) | Inside each part file, in its Parquet key-value metadata; no `identity.json`; clearing the cache becomes housekeeping (Eugene, 18 November) | A format change: where a record lives |
| 31 | Conforming and "readable" as two levels (S4) | Dropped: a file conforms or it does not, and one behaviour rule says what a reader still accepts (Eugene, 18 November) | The validators' output: no `readable` column, no `conformance_level()` |
| 26 | The layout of the analysis folder (A2) | One folder per stage, holding its record and its tables (Eugene, 18 November) | A format change: paths |
| 27 | Where fits live (T6) | A cache outside the stages, `fits/response/`, filled by part files, with its own clearing rule; `fits/grna/` reserved (Eugene, 18 November) | A format change: paths |
| 25 | What kind of screen it was (S10) | Optional `perturbation_type`, `organism`, `genome_build`, recorded now because they cannot be recovered later (Eugene, 18 November) | A new kind of perturbation is a minor addition |

### Outside the format

Three matters came up that the format does not decide and should not
([section 10](#10-fixed-now-free-to-grow) lists them with the other suggested
changes to the packages):

- **Which response-fitting method is the shared default.** The format records
  the method with every saved fit; which one to use is a scientific decision
  for the method's authors.
  The measurement in T6 shows the size of the difference.
- **Upstream library fixes.** S7's restrictions exist because of measured
  problems: R's MuData writes a logical as an integer, reads a boolean as text,
  writes a missing string as empty, and writes missing integers in a form
  Python cannot read; and current pandas writes strings in a form R cannot
  read. Reporting them upstream would let a later format lift the restrictions.
  Nothing has been reported yet.
- **Review.** The format has not been reviewed.

---

## 9. Compared with IGVF and pertpy

### IGVF

The IGVF pipeline's output is the closest existing format: one `.h5mu` file per
run, `inference_mudata.h5mu`, with result tables beside it. Sources are given at
the top of this document.

**Where the two agree**

| Area | Both |
|---|---|
| Container | MuData `.h5mu`; cells × features; raw integer counts; sparse |
| Identifiers | Cells, features and table rows are named by identifier, not position |
| MOI | Stored in the file as `high` or `low` |
| Control group | `complement` or `nt_cells` |
| Two levels of testing | Per element and per guide there; union and singleton here |
| Result tables | Parquet with zstd as the working format |
| Continuing an analysis | The parent data is never overwritten, and a hash ties the continuation to it |
| Library problems met | Mixed index types in scipy; missing strings in HDF5 |

**Where they differ**

| Area | IGVF | This format | Which is better, and why | See |
|---|---|---|---|---|
| Data and analysis | One file holds both | Separate | This format: one dataset, many analyses. IGVF's own code shows the cost of patching results into the data file | S12 |
| Cells | Only cells that passed QC | All cells; each analysis says which it used | This format: QC is an analysis choice | T2 |
| Assignments | 0/1 matrix layer in the data file | Table in the analysis | This format, for the same reason. An export to the matrix is worth having | T1 |
| R and Python | Python writes; R reads the matrices, MOI and pairs | Both read and write everything | This format: it is the first requirement | S7 |
| Format number, dataset id | Neither | Both | This format | S10 |
| Matrix layout | Compressed by cell | Compressed by feature | This format: one gene is one read, with no chunk files | S3 |
| Results | Three copies; method-prefixed columns; annotations in every row; adjusted p-values stored | One copy; keys and results only | This format: nothing stored twice | T5 |
| Settings | Several run-level files | One record per stage | This format, for continuing an analysis | A2 |
| Link to the data | SHA-256 of the whole file | Dataset id and identifier hashes | Mixed: IGVF's catches in-place edits; this format's survives recompression and is fast | A1 |
| Count type | Smallest unsigned type that fits | Fixed 32-bit | A close call; kept for a type that never depends on the data | S4 |
| gRNA annotation | Type, spacer, gRNA and target positions, target key | Target and target position | **IGVF.** Taken in part now (S8); the rest planned | S8, 10 |
| Gene annotation | Symbol, position, mitochondrial and ribosomal flags | Unspecified | **IGVF.** Planned | 10 |
| Reference genome | Build per gRNA; annotation file hash in the manifest | Not recorded | **IGVF.** Planned | 10 |
| Other modalities | Optional `hashing` | Any, ignored unless named | Taken from IGVF | S2 |
| Names | `gene`, `guide`, `guide_id`, `gene_id`, `intended_target_name` | `rna`, `grna`, `grna_id`, `response_id`, `grna_target` | Neither is a standard. An importer maps one to the other | S2 |

**What this format took from IGVF:** a file may carry modalities an analysis
does not use; the gRNA table says whose position its coordinates are; and the
list of annotation to add in [section 10](#10-fixed-now-free-to-grow).

**What an IGVF importer has to do:** rename the modalities and columns,
re-orient the matrices, turn the assignment layer into a table, convert
`percent_mito` to a proportion, and move the pairs and results into an analysis
folder. The cells IGVF removed in QC cannot be recovered from its file.

### pertpy

pertpy defines no file format. It works on in-memory AnnData objects, and the
user passes the name of every column *(verified, source and tutorials at commit
`4924c25`, 18 November 2026; no pertpy dataset was downloaded)*. Its habits:

| Area | pertpy | This format |
|---|---|---|
| Container | 54 of its 55 bundled datasets are single `.h5ad` files; one is MuData with `rna`, `adt`, `hto`, `gdo` | MuData |
| What a perturbation is | A text label on each cell (`perturbation`, control label `control` or `NT`) | Rows of (gRNA, cell) in a table |
| Several gRNAs in a cell | Names joined with `+`, capped at five, then `multiple` | One row each |
| gRNA-to-target mapping | A per-cell column | A property of the gRNA |
| Assignments | 0/1 matrix layer, `assigned_guides` | Table in the analysis |
| Results | In-memory tables: `variable`, `p_value`, `adj_p_value`, `log_fc`, `statistic` | Parquet tables |

The per-cell label works at low MOI and breaks down at high MOI; a label can
be derived from this format's table when a tool wants one. pertpy's `log_fc`
is base 2 but does not say so. The one thing pertpy settles for this format is
that `rna` is the shared name for the RNA modality; its `adt` and `gdo` are
Seurat's assay names and appear in one dataset.

---

## 10. Fixed now, free to grow

### What would be expensive to change after the freeze

| Kind | What |
|---|---|
| Structure | Data separate from analysis; one record per stage; tables keyed by identifier; the matrix layout; the number types |
| Names | `grna`; the table and column names in sections 4 to 6; the stage names; the covariate prefix rule |
| Removals | A column or file that readers rely on cannot be taken away |

Everything else can be added with a minor version, because readers ignore what
they do not know (principle 6), or loosened, because old files stay valid. The
table of planned additions below says which additions are not safe to ignore.

### To settle before the freeze

- **Review** ([`REVIEW.md`](REVIEW.md) says where to look).
- **An analysis at scale.** The data file has been through one real large dataset, in both languages and disk-backed (Appendix A). No analysis larger than the 1,000-cell examples has been written or read in this format.

### Planned additions to the format

None of these is in format 0.1. The last column says what kind of version
change each would be once the format is frozen (S10): **minor** if an older
reader can safely ignore it, **major** if it could misread the file.

| Addition | Where | Why wait | Version |
|---|---|---|---|
| gRNA type: targeting, non-targeting, safe-targeting, other negative control | gRNA table | It overlaps with marking controls by `grna_target` equal to `non-targeting`, which both codebases depend on, and a safe-targeting gRNA has a real target and is a control, which that marker cannot say. This needs a decision, not only a column | minor if the column is information and the marker stays the rule; major if the column becomes the rule |
| Spacer sequence and PAM | gRNA table | Known at import and not derivable, like the entries of S10, but they serve guide design and reanalysis, not reading results | minor |
| The gRNA's own position and strand | gRNA table | sceptre does not read it | minor |
| Gene symbol, gene position, mitochondrial flag | Feature table of the measured modality | The flag would record which genes `p_mito` counted; it needs the importers to keep it | minor |
| A hash of the gene annotation used | File metadata | The genome build is recorded now (S10); the annotation is not | minor |
| How the gRNAs were captured (IGVF's `capture_method`) | File metadata | sceptre does not use it | minor |
| A declared type for each modality, beside its name | File metadata | No dataset needs it yet | minor |
| A gRNA-to-target table, for a gRNA with two targets | Data file | A scientific question first (S8) | major: a file using it has gRNAs an old reader would assign to one target |
| A small stamp in every stage table's Parquet metadata: format version, dataset `id`, stage, table name | Key-value metadata of each table file, as the fit cache's parts carry their identity (T6) | A convenience, not needed for correctness: a results file copied out of its folder would identify itself, and a validator could confirm each file belongs. It is a second copy of facts the stage record holds, so the two would have to be kept in agreement. Agreed with Eugene on 18 November that it waits | minor: an old reader ignores the metadata |
| Why each cell was removed (`removed_by`) | `qc/cells.parquet` | sceptre does not record it (T2) | minor |
| An assignment score (posterior probability or p-value) | `assignment/assignments.parquet` | Depends on the assignment method | minor |
| Per-gRNA values an assignment method computes, such as fitted thresholds | A new optional table of the assignment stage, keyed by `grna_id` | They are computed by an analysis, so they are not data and not settings. Nothing needs them now: sceptre's thresholding uses one threshold, stored as a setting | minor |
| `log_p_value`, for p-values below a double's range | Results tables | No method needs it today | minor |
| Per-pair counts: cells carrying the perturbation, and cells with a nonzero count in each group | Results or pair-QC tables | The development builds do not put them in the result tables | minor |
| A SHA-256 per file, and one of the data file | Stage records, `analysis.json` | For those who want IGVF's stronger integrity check | minor |
| Cached fits of the gRNA models | `fits/grna/`, a reserved folder | Not needed yet; its key, columns and clearing rule differ from the response cache's (T6) | minor |
| A screen that combines targets on one vector (combinatorial perturbations) | Data file, S8.7 | Out of scope: sceptre cannot analyse them | minor, because today's readers refuse such files |
| Which package and version wrote the data file | File metadata | Stage records already say this for each stage; the data file does not. The format and the packages are versioned separately, so this is information, not a rule | minor |
| True booleans and nullable integers in the data file | S7 | Needs the R MuData and pandas fixes upstream | major: an old R reader misreads them |
| A Zarr variant | S1 | Not tested | a second container, not a version of this one |
| Two measured modalities tested in one analysis | Analysis folder | No dataset needs it yet | major: it changes what `response_modality` means |

### Suggested changes to what the packages produce

These are changes to sceptre and pysceptre, not to the format. None is made
now: functionality is frozen while the foundations are rebuilt. Each would
remove a workaround in the format or its readers.

| Change | What it would remove | Evidence |
|---|---|---|
| Bonferroni results carry `fold_change` and `se_fold_change` of the reported gRNA, like the other strategies | The one case where `log_2_fold_change` must be stored (T5) | `R/s4_analysis_functs_2.R` keeps only the log2 value |
| One expression for `log_2_fold_change` everywhere | A last-bit difference waiting to happen: pysceptre has `np.log2(fc)` in `pipeline/analysis_result.py` and `log(fc)/log(2)` in `pipeline/result_fields.py` | *(verified)* |
| Keep control-group membership as data | Writers splitting group names on `&`; the ban on `&` in gRNA identifiers (T4) | Both codebases split names |
| Record which QC rule removed each cell | The summaries-only QC record (T2) | `R/cellwise_qc_functs.R` |
| Hold resampled statistics as one list per pair | Thousands of padded columns; the failure when keeping them for 1,406 pairs (T5) | *(measured)* |
| Sort results by p-value when showing them, not when storing them | Order depending on job order among ties (T5) | *(inferred)* |
| Check a saved fit's recorded method before reusing it; keep one shared list of method names | Fits trusted without a check (T6) | Both codebases |
| Decide the shared default fitting method | `batched` gives a dispersion of 1000 on 2 of 100 responses where the other two methods give 0.03 to 0.07 (T6) | *(measured)* |
| Clear later stages' values when a stage is rerun | Stale values behind a "not run" flag in R (A4) | `R/s4_helpers.R` |
| Check at import that counts are integers | Non-integer "counts" reaching a file (S4) | `R/check_functions.R` |
| Accept modality-name covariates (`rna_n_umis`), keeping `response_n_umis` as an alias; add a response-modality argument and a modality name at import | The translation layer in the readers and writers (S2, S5) | — |
| Rename `chr`, `start`, `end` in the gRNA table, or accept both | The same (S8) | — |
| Resolve a dataset's `id` to a location from the user's configuration (a search path or a small registry) | Passing the data file's path by hand every time (A1) | Not built |
| Let pysceptre's formula parser take covariate names with dots | The alias in the prototype's Python reader; a pysceptre analysis with such a covariate fails today (S5) | *(measured)* `model_matrix` raises on `S.Score` |
| Decide what a gRNA under two targets means | The refusal (S8) | R warns, pysceptre counts |
| Combine a vector's gRNAs in the assignment stage, not at import: a setting in the assignment record, and an assignment table that names the unit assigned | Summing at import, which a written-once data file of per-gRNA counts cannot hold (S8) | sceptre sums within a vector in `import_data()` |
| Report the R MuData and pandas problems upstream | The restrictions of S7 | Appendix A |

---

## 11. How the format is defined

### The pieces

| Piece | Where | What it is for |
|---|---|---|
| The rules | `SPEC.md` | Numbered statements, each a file rule, a behaviour rule or guidance, with a scope, a glossary and the library versions it is defined against. A validator's report cites them (`S4.1`) |
| The reviewer's page | `REVIEW.md` | One page: the format in ten lines, the choices that are expensive to change, the questions that need a view, and what can be skipped |
| The reasons | this document | Why each rule is what it is |
| Data-file definition | `spec/screen.json` | The structure of a data file as data: required groups, number types, encodings, allowed column kinds, required columns and metadata. Both validators read it |
| Table definitions | `spec/tables.json` | Every table's stage, key, columns and types. Both prototype libraries and both validators read this one file, so the two languages cannot drift apart |
| Record schemas | `spec/analysis.schema.json`, `spec/stage.schema.json`, `spec/fits.schema.json` | JSON Schemas for `analysis.json`, the six stage records, and the identity each part of the fit cache carries |
| Data-file rules that need computing | `SPEC.md`, section 2 | The order of cells within a feature, uniqueness, the same cells in every modality: code in each language |
| Conformance files | `conformance/`, built by `make_conformance.py` | 23 conforming files, written by R and by Python, and 73 that each break a rule, with the rules a validator must report. Five of the 73 store counts or index arrays in another number type, which a reader must still read (S4.4) |
| Validators | `validate_screen()`, `validate_analysis()` in both prototype libraries | Check a data file or an analysis folder without building a package object, and report every problem, not only the first |

**Why define it this way.**

- Writing the schemas exposes ambiguities in the prose while they are still
  free to fix. The list below is what it exposed.
- Reviewers get something exact to disagree with.
- People can work with the format through AI assistants, which need rules a
  machine can read and a validator to check their own output against.

### Three kinds of statement

A first version of the specification mixed three things: facts about files,
duties of software, and advice. Some of its MUSTs could not be tested by any
program (that "facts about the cell" are in the shared table, for instance).
Every statement is now one of:

- a **file rule**, true of every conforming file and tested by the validators;
- a **behaviour rule**, about what a reader or writer does, tested by the
  round-trip checks and the conformance run;
- **guidance**, which is advice or a definition and binds nobody.

`SPEC.md` gives the counts and names the behaviour rules that nothing tests
yet. The sorting is meant to be final. The wording of individual rules is not,
and the specification says so.

**What the specification deliberately leaves out.** It fixes the form of
results and not their meaning: what a p-value tests is defined by sceptre's
method documentation. Its status line says plainly how little it has been
tested. Left for later: defining each result column's meaning, tightening
each rule's wording, schemas for method-specific settings, and an
implementation by someone who has only the specification.

### What was checked

*(measured)* Both validators accept all 23 conforming files and report exactly the expected
rules for all 73 broken ones. The two languages report the same rules for
every case. Both readers read each of the five files of S4.4 to the same
counts as the conforming file it was made from. 57 rules are exercised by at least one broken
file. The validators check two more rules that have no conformance file: S11.1
(no filter other than gzip), because R cannot read a file compressed otherwise
at all, and T0.3 (row groups of at most a million rows), because the file would
be large.

The JSON Schemas are checked with standard libraries: `jsonschema` in Python
and `jsonvalidate` (the `ajv` engine) in R. The schemas use draft 7, which both
support.

### What the formal definition exposed

Writing the rules down precisely raised twelve questions the prose had left
open. Each is now settled as below. Those marked *Adopted* follow the
recommendations of the side session of 18 November, which Eugene read without
objecting; he can overturn any of them.

| # | Question | Settled as | Rule |
|---|---|---|---|
| 1 | Is a data file with another count type (16-bit unsigned, as IGVF writes) a conforming file? | No, and a reader must still read it. A separate "readable" level was tried and dropped (section 8, row 31) | S4.1, S4.4 |
| 2 | Are `n_umis`, `n_nonzero` and `p_mito` required? | `n_umis` and `n_nonzero` are required for every modality; `p_mito` is optional (*Adopted*) | S6 |
| 3 | Must a required table exist when it has no rows? | Yes | A5.4 |
| 4 | When do the pair-QC and results tables have a `grna_id` column? | A rule, fixed by strategy and stage, and checked (*Adopted*) | T3.3, T5.6 |
| 5 | A file in the folder that no record lists: allowed or stale? | Allowed if the format does not define it; an error if it is one of the format's tables | A4.1, A4.3 |
| 6 | Is `seed` required? | Yes, an integer | A2.1 |
| 7 | Which characters may a covariate name have? | Letters, digits, underscores and dots, starting with a letter; formulas accept the same (*Adopted*) | S5.4, T7.1 |
| 8 | Do the rules apply to a modality no analysis uses? | Yes, all of them, as a plain rule (*Adopted*) | S2.7 |
| 9 | Method-specific settings (`hyperparameters` of an assignment method) | Any JSON object; no schema for now. This is why a threshold of `5` comes back as an integer | A2.1 |
| 10 | The form of `written_at` | UTC, `YYYY-MM-DDTHH:MM:SSZ` | A2.4 |
| 11 | Exactly what is hashed | Identifiers: UTF-8, joined by line feeds, which identifiers must not contain (*Adopted*). Design matrix: little-endian doubles, column by column | A1.3, S9.4, T7.3 |
| 12 | Signed or unsigned 32-bit integers | Signed; both validators can tell | S4, S7.1 |

### What counts as evidence about the format

Five points, agreed with Eugene on 18 November, that shape how the format is
checked:

1. **There is no canonical writer.** The standard writers (R's MuData, Python's
   mudata and anndata, nanoparquet, pyarrow) do not produce conforming files as
   they come; the writers that do are the prototype's. So the format is defined
   by `SPEC.md` and `spec/`, and the validators are the arbiter. Any code whose
   files pass the validators and read back correctly in the standard libraries
   is a conforming writer.
   *Custom writers, generic readers plus a thin layer* (accepted by Eugene,
   18 November). Each conforming writer is a wrapper around the standard
   writer, not a replacement. R needs the most: counts written by hand, and
   column conversions. Python needs a setting, types set before writing, and
   the same conversions. Readers can be generic: the standard libraries open a
   conforming file and return the right values. What they do not do is
   interpret it: the version check, the link between an analysis and its data,
   the refusals the specification requires, the covariate table with modality
   prefixes, true/false arriving as a categorical, and derived values. (Until
   18 November a Python reader also had to bring the two index arrays to one
   width; S4 now stores them so.) The reason is that the format allows only
   what both languages read correctly today and makes the writers stay inside
   that set. The burden falls on the two writers we control and not on
   readers, which include tools we do not. Every R-side restriction in S7
   traces to MuData's reader or writer, not to SingleCellExperiment: it cannot
   open pandas' nullable strings, reads a Python boolean as text and a 64-bit
   integer column as all missing, and writes a logical as an integer, missing
   integers in a form Python cannot open, and missing text as empty. If MuData
   is fixed or replaced, several restrictions could be lifted. None has been
   reported upstream.
2. **A round trip through one's own reader and writer is weak evidence.** It
   tests that the pair is self-consistent and cannot catch an error they share.
   It is a fair requirement of a package that adopts the format. Stronger
   evidence is independent: a file read by a reader we did not write; data read
   back and compared with its source, not with a rewrite; R's file compared
   with Python's.
3. **Whether sceptre and pysceptre compute the same thing from the same files
   is not the format's job.** That is parity between the packages. The
   format's part is to make the comparison possible: it records the formula,
   the methods, the seed, the settings and the writer; it stores the design's
   term names and hash; and it stops fits crossing where the method differs.
4. **So the cross-language continuation** (an analysis stopped in R and
   finished in pysceptre, agreeing to 8e-12) **is a parity result**, not
   evidence about the format.
5. **One package-dependent question is the format's concern: sufficiency.** Do
   the files hold everything an analysis needs? It is a property of the
   format's design that can only be tested by rebuilding a package's object.

**Two suites follow from this.**

- *The package-free suite* (`verify/`) loads neither sceptre nor pysceptre and
  runs on committed files with standard libraries only. *(measured, one Mac)*
  About 790 checks in five groups, each in R and in Python, all passing:
  - an encoding matrix: each kind of value written by each language's standard
    libraries and read by the other's;
  - two writers on seeded random content, compared array by array and with the
    source;
  - readers we did not write: anndata, scanpy, polars, duckdb, HDF5Array, arrow;
  - the validators and readers on the conformance set;
  - what a reader derives from the files, against values sceptre held, captured
    once and stored beside the examples.
- *The package suite* rebuilt each package's object from the files and
  compared it with the original, which is the sufficiency test. It found most
  of the gaps in earlier drafts. It also held the parity checks.

**The package suite was removed for now** *(Decided by Eugene, 18 November)*.
Nobody else could run it: it needed the unreviewed `sceptre3` builds of both
packages on one machine, and its code reached into internals that the rebuild
of the packages will change. It is in git history up to commit `7e10af4` and
returns inside the packages when they adopt the format. The costs are real:

- **No sufficiency check.** Rebuilding a real sceptre object is the only test
  that the files hold everything an analysis needs. The format passed it at
  every step up to commit `7e10af4`. Changes after that are not covered until the
  packages test them.
- **Frozen examples.** The example generator ran sceptre, so the committed
  example files can no longer be regenerated. A format change from now on
  converts them in place with package-free code, written with the change and
  kept. That is also a direct test that old files can be carried forward.
- **Frozen expectations.** `example/*/expected/` was produced by sceptre once,
  by the removed code, and stays.
- **No parity checks here.** The continuation across languages and the
  agreement of the design matrices went with it. They compare the packages,
  which point 3 above already puts outside the format.

**What the package-free suite found.** Writing it changed one thing and raised
two questions, both since settled:

- *Fixed.* The two writers stored a feature table's text differently: R as
  text, Python as a categorical, because mudata converts text with repeated
  values as it writes. Both conformed and both read back the same, but the
  files were not the same. The Python writer now stores it as R does.
- *Settled: one common scanpy function failed on the counts as anndata read
  them.* The format stored column pointers as 64-bit integers and cell
  positions as 32-bit, and scipy rejects mixed index widths in some operations.
  *(measured by a side session on 18 November: scanpy 1.12.5, scipy 1.18.1,
  anndata 0.13.4)* Most of scanpy worked on a raw read;
  `scanpy.pp.calculate_qc_metrics` failed with "Output dtype not compatible
  with inputs", as it does for any scipy matrix with those widths. Eugene
  decided on 18 November to store both arrays in one width by scipy's rule; see
  S4 for the reason, the cost and what was checked. The suite now runs that
  function on every file as a check.
- *Settled: the fit cache's record sat inside its folder of parts.* A reader
  pointed at `fits/response/` as a Parquet dataset met `identity.json` and
  failed. Eugene decided on 18 November that the identity goes inside each
  part file and the folder holds parts only; see T6. The suite now opens the
  folder itself with six readers as a check.

**What the package-free suite cannot show.** Sufficiency, and anything about
results. The validators and both writers also still come from one author; the
independent readers reduce that weakness and only an outside implementer
removes it. The suite has run on one operating system: a workflow file for
Linux and macOS is in the hub and has not been run.

**The prototype's reader and writer are temporary** *(Decided by Eugene, 18
November)*. See "Where this lives" below.

### What a code review found

A review of the specification and the code that checks it (18 November) found
nine faults. None was in the layout. All are fixed, and each fix has a check:

| Found | Fixed by |
|---|---|
| Cached fits were matched on the formula as exact text, so two writers that space a formula differently could not share fits | Formulas are compared with white space removed, in both libraries and both validators (T6.5) |
| A stage written again after a crash could leave a table of the first attempt that its new record did not list | A writer deletes a stage folder that has no record before writing the stage (A2.3) |
| A file with a modality no analysis uses could not be written again after Python read it, and R dropped that modality on a round trip | Both readers carry every modality, and a copy keeps them (S2.6) |
| The readers followed a record's `stage` and file paths unchecked | Both readers refuse a record of another stage or a path outside the stage's folder (A5.3) |
| The readers refused a vector named like another vector's gRNA; the specification asked this only of the other vector rule | The specification now says so (S8.9) |
| `p_mito` stored as integers was rejected though the specification says a number | Integers are allowed |
| The Python validator did not decode a table's index attribute | It does |
| Both validators read every table twice | Once |
| The conformance file for an extra modality was made by copying raw HDF5 groups, which R could not read | It is written as a proper MuData file; the review's own fix exposed this |

### Taken from Apache Arrow's practice

Four things, on the side session's recommendation:

- **A two-part version** (S10). Decided by Eugene.
- **The forward-compatibility rule**: a reader either reads a file correctly or
  detects that it cannot. It is at the top of `SPEC.md`.
- **The change rule**: nothing enters the specification until both libraries
  implement it and a conformance file exercises it.
- **A table of what each library supports** (`SPEC.md`, section 7).

Not taken: a voting process, a separate repository of test files, a folder of
proposals. An archive of files written by each released version is worth
starting at the freeze, so that later readers can be tested against them.

### Where this lives

*(Decided by Eugene, 18 November.)*

- **Now, while the format is being decided:** the readers and writers live
  outside the packages, in `format/`. That is the right place while the format
  still changes: the layout, the names, the version and the coordinate
  convention each changed in one edit, with no package release or review.
- **They are a demonstration, not a reference.** Prototype code by one author.
  Nobody should depend on it as a package's reader.
- **When this phase ends:** the plain layer (files to standard containers)
  moves into sceptre and pysceptre and becomes their real reader and writer.
  The layer that reached into package internals (the package suite, removed
  for now) is replaced
  by a small adapter in each package. The copy here is retired once both
  packages pass the conformance suite.
- **sceptre and pysceptre are then the two reference implementations.** That
  fits the change rule of the specification: nothing enters it until both
  implement it and a conformance file exercises it.
- **What stays with the specification for good:** the specification, `spec/`,
  the conformance files with their expected results, the validators, and a
  harness that runs any reader or writer against them (`verify/README.md`).

**Precedents** *(read by a side session on 18 November from each project's own
repository; I did not re-read them)*. In all seven, the reference
implementations are the real libraries, and what sits beside the specification
is neutral: test files, expected values, a harness, sometimes a validator.
None keeps a separate, non-production reference reader and writer once real
libraries exist. Seven is a small sample.

| Project | What it does |
|---|---|
| Apache Parquet | A format change needs two implementations that interoperate, one of them the Java library; a prototype is allowed while a change is discussed. A separate repository holds test files |
| Apache Arrow | At least two reference implementations, each inside a complete library, with integration tests. Test data is described in a neutral JSON form, and every producer is tested against every consumer |
| CommonMark | Two named reference implementations, which are also the libraries people use; the examples in the specification are the conformance tests, with a runner that takes any program |
| JSON Schema, TOML | No reference implementation; a test suite that does not depend on a language |
| BIDS | No reference writer; one official validator |
| SOMA | A specification and its first implementation, with Python and R interfaces |
| GeoParquet | A specification, a folder of examples, and a list of implementations |

Two things here go beyond those sources and are a judgement, not a precedent:
the harness takes any reader or writer as an input, so that a package can be
tested against this suite without the suite importing it; and there is one
validator per language, kept in agreement by the conformance set, with no
third to be added.

### How adoption would go

For the record. None of this has started, and it is not a format question.

1. **Finish the format side first**: the package-free suite, the harness, and
   a review of structure and names. Package code written before that review
   risks being rewritten.
2. **Settle what the package branches start from.** The package layer of the
   prototype is written against unreviewed `sceptre3` builds, and the plan for
   rebuilding the packages is not decided. This is the real blocker.
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
R 4.6.0 with MuData 1.16.0, rhdf5 2.56.0, nanoparquet 0.5.2. Python 3.12 with
mudata 0.4.1, anndata 0.13.4, pandas 3.1.0, pyarrow 26.0.0, h5py 3.16.0.

**Counts** (200,000 cells × 4,000 features, 48 million nonzeros, 6% dense):

| Stored as | File | Read all, 200 features at a time | 200 single features |
|---|---|---|---|
| double, uncompressed | 588 MB | 0.07 s | 0.02 s |
| double, gzip | 103 MB | 1.17 s | 0.17 s |
| 32-bit integer, uncompressed | 395 MB | 0.04 s | 0.02 s |
| 32-bit integer, gzip | 98 MB | 1.00 s | 0.16 s |
| 16-bit unsigned, uncompressed | 298 MB | 0.03 s | 0.01 s |

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

**Data-file encodings:** the table in S7, plus:

| Test | Result |
|---|---|
| Python writes 32-bit counts and 64-bit column pointers (the rule before 18 November) | R reads it in memory and disk-backed |
| R writes the same types (six extra lines) | Python reads it with those types |
| Nested metadata record with scalars | both directions correct; R needs `encoding-scalar` set to write a scalar |
| Shared cell table | both directions correct |
| gzip-compressed counts | R reads them |
| Python nullable boolean and nullable integer columns | R reads both correctly; R cannot write either in a form Python reads |
| Design matrix, R against pysceptre, on the example | bitwise equal |
| Fitting methods on the example (real data, 500 cells) | see the table in T6 |
| A sparse matrix with 64-bit pointers and 32-bit indices in scipy | conversion to coordinate form fails until both are one type |
| A double written to JSON with 15 significant digits | the significance threshold came back smaller; 1 pair misclassified. With 17 digits: exact |
| R sceptre keeping resampled statistics as columns | fails for 1,406 pairs; works for 40 (4,999 columns) |
| `log_2_fold_change` recomputed as `log(x) / log(2)`, 17,219 result rows of the examples | identical to the stored value on every row, in R and in Python |
| The same recomputed with `log2(x)` | a different last bit on about 30% of rows |
| `fold_change` recomputed from the log2 value | a different last bit on a few rows per table |

Not measured: Zarr; other HDF5 filters; a second platform.

### Scale test on real data

One run, 18 November 2026, on the Gasperini et al. (2019) at-scale screen:
207,324 cells, 13,135 genes, 13,077 gRNAs, 699 million nonzero RNA counts. One
cluster node, 16 cores, 90 GB of memory. R 4.5.2 with MuData 1.14.0,
rhdf5 2.54.1 and HDF5Array 1.38.0; Python with the versions above. The scripts
and their recorded results are in `scale-test/`.

**What was checked, and held:**

| Check | Result |
|---|---|
| The file R wrote against the file Python wrote | Every count array identical; every cell and feature table identical; same number types |
| Both validators on both files | Conforming |
| Read back in memory, R and Python | Same shape, nonzeros and count sums as the source |
| 200 random genes read from disk one at a time, R (rhdf5) and Python (h5py) | Equal to the in-memory values |
| The reader the `sceptre3` branch uses: MuData's backed reader, one response at a time | Sums identical to direct HDF5 reads; the 64-bit column pointers are read correctly |
| An assignment table of 6.6 million rows, written by Python and rewritten by R | Equal |

**Sizes and times.** The file is 5.7 GB from either language, against 8.6 GB
for the same data with counts as doubles. Times are from a cold start;
memory is the peak of the whole process, which for later steps includes data
already held.

| Step | Python | R |
|---|---|---|
| Write the data file | 12 s, 15 GB | 35 s, 25 GB |
| Validate the data file | 116 s, 29 GB | 79 s, 21 GB |
| Read the data file in memory | 37 s, 15 GB | 106 s, 25 GB |
| Read 200 genes from disk, one at a time | 3 s | 21 s with rhdf5; 3 s with HDF5Array |
| Open disk-backed with MuData and read 200 responses | — | 5 s to open, 11 s to read, under 1 GB |
| Read the assignment table | 5 s | 1 s |
| Write the assignment table | 111 s (before the fix below) | 1 s |

**What the run found.**

- *Dead space.* In the first run the Python-written file was 14.3 GB. Both
  prototype writers wrote the counts as doubles and then rewrote the arrays in
  their final types, and HDF5 does not give the first copy's space back. The
  500-cell examples hid this. Both writers now write each array once, in its
  final type: the file fell to 5.7 GB, Python's write from 106 s to 12 s and
  R's from 399 s to 35 s. This is a lesson about writers, not a change to the
  format.
- *A slow sort.* The prototype's Python table writer sorted rows one at a time
  in Python: 111 s for 6.6 million rows, against 1 s in R. It now sorts with
  numpy (2 s for 2 million rows on a laptop); this was not rerun at scale.
- *A mistake in the test itself.* The first run reported that R's disk-backed
  sums disagreed with the in-memory ones. HDF5Array presents the stored
  cells × features matrix as features × cells, and the test took columns. With
  the orientation handled, the sums agree.
- *The validators are not built for size.* They read whole arrays into memory:
  about 29 GB in Python and 21 GB in R for this file. That is acceptable for a
  check run once per dataset and would need chunking for larger ones.

**What this does not cover.** An analysis at this scale: only the data file and
one table were tested, not the stage records, results or fits. A second
dataset. Datasets large enough to pass 2³¹ nonzeros: this one has 699 million,
well under that limit, so the file never needed its 64-bit column pointers.
The run shows that readers, including the backed reader, handle the wider type
correctly. That the wider type is needed for larger data remains arithmetic,
not a measurement, and R's in-memory matrix would not fit such data.

**After the index-type decision (18 November).** These files were written under
the earlier rule, with 64-bit column pointers beside 32-bit cell positions.
Under the present rule both arrays of this dataset are 32-bit, and the file is
the same size to within the pointers (8 bytes per feature). The test was not
run again.

## Appendix B. What the prototype implements

Everything in sections 3 to 6, and the validators of section 11. The code is
`format/sceptre_format.R` and `format/sceptre_format.py`; `example/` holds
files in the format.

The checks that remain are the package-free suite, `verify/`, described in
section 11. The table below is a record of what the package suite showed when
it last ran, at commit `7e10af4`, before it was removed. Nothing reruns it.

**Checked in both languages by the package suite, as of commit `7e10af4`** (notebook section in brackets):

| What | Result |
|---|---|
| R writes, R reads back as a `sceptre_object` [6] | Every field identical or equal for all 7 analyses, on high MOI and low MOI |
| R writes, Python reads into pysceptre's `Analysis` [7] | Equal, part by part, to what pysceptre's own reader gives for the same analyses in the October layout |
| Python rewrites everything, R reads the Python copy [6, 7] | Objects as above; all stage records equal; all tables identical; same number types in the data file |
| Python continues an analysis R stopped after QC [8] (a parity result between the packages, not evidence about the format) | p-values agree with R's to a relative 8e-12; same discoveries |
| R resumes from Python-written files [8] | Results identical to an analysis that was never saved |
| Negative-control groups read from the table [4] | Every member is a non-targeting gRNA; every tested group has members |
| Per-target unions, per-control lists, gRNAs per cell derived from assignments [4] | Equal to what the packages hold |
| `significant` derived from the stored threshold [4] | Equal to R's column for every pair |
| `log_2_fold_change` derived from `fold_change` [4] | Identical to R's column for every pair; stored only in the bonferroni tables |
| Design matrix rebuilt from the formula [4] | Same term names and the same hash in R and Python |
| Resampled statistics as a list column [4] | 40 pairs, 499 to 4,999 statistics each, restored as the packages' columns |
| The three gRNA strategies: union, singleton, bonferroni [5–7] | All round-trip |
| Low MOI: the "maximum" method, flagged cells, NT cells as controls [5–7] | Round-trips |
| A table as a folder of part files, with an unknown column and dictionary-encoded identifiers [9] | Read to the same object, from parts written by either language |
| The same data under the modality name `atac`, with a `hashing` modality beside it, written by each language [9] | Read to the same object by both; without an analysis, a reader asks which modality is the response |
| A subset data file, and a re-imported one [9] | Refused |
| A table with fewer rows than its stage wrote [9] | Refused |
| Clearing a stage [9] | Its record and files and all later ones are removed |
| Saved fits whose recorded method differs [9] | Not reused, with a message |
| Every allowed column type, written by each language and read by the other [9] | Same result both ways |
| The validators on the conformance files [10] | Both languages report the same rules for every file |
| A covariate with a dot in its name [9] | R and Python build the same design matrix (Python through an alias) |

**Limits that remain.** These are limits of the prototype or of what was
tested, not open choices:

- The writers save a finished object in one go. They produce the same files a
  stage-by-stage writer would, but nothing here writes a stage while it runs.
  That belongs in the packages.
- The packages are unchanged: they still call the response modality `response`
  and the gRNA table's coordinates `chr`, `start`, `end`. The prototype
  translates on the way in and out.
- The readers build each package's object through its internal accessors, on
  the `sceptre3` builds of R sceptre and of the pysceptre fork. In R, one
  per-response count that sceptre computes during pairwise QC is recomputed on
  read with sceptre's own routine.
- Three object fields come back equal but not identical in R: the gRNA table's
  `start` column (`target_start` in the file; an integer with missing values, stored as doubles), an
  assignment threshold (a double `5` that JSON returns as an integer), and one
  table's class.
- Refused with a message: a gRNA under two targets; a formula outside the
  additive subset.
- Tested once at scale, for the data file only (Appendix A): 207,324 cells.
  Not tested: an analysis at scale, a second platform, Zarr, the R readers
  anndataR and zellkonverter, methods other than
  those in the examples, and a file whose second modality
  is real data (the `atac` and `hashing` test reuses the RNA counts and a
  synthetic matrix).

## Appendix C. The examples

| Dataset | What it is | Analyses |
|---|---|---|
| `example/highmoi` | sceptre's built-in high-MOI data: 500 cells, 100 responses, 60 gRNAs | `left_sided`, `two_sided`, `singleton`, `bonferroni`, `stopped_after_qc` (through QC only), `with_resamples` (discovery only, 40 pairs, every resampled statistic), and `continued_in_python` (R through QC, then pysceptre) |
| `example/lowmoi` | sceptre's built-in low-MOI data: 1,000 cells, 100 responses, 50 gRNAs | `default`: the "maximum" assignment method and the NT cells as controls |
| `example/types` | 6 cells, synthetic | none: a data file with a factor whose levels are not alphabetical, a text covariate, a true/false covariate, a covariate with missing values, and a gRNA table with missing text and a missing integer |

These files are frozen: the code that generated them needed sceptre and was
removed after commit `7e10af4`.
