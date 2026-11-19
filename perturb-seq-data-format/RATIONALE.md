# The Perturb-seq data file format: rationale

[`SPEC.md`](SPEC.md) states the rules of the Perturb-seq data file format,
numbered, with no reasons. This document gives the reasons: for each choice,
the alternatives, the evidence, what it locks in, and how it compares with the
IGVF pipeline's format and with pertpy. Sections here and there share their
letters and numbers (S3, S8, ...), so rule `S4.1` of the specification is
explained in section S4 below.

The format is meant for any software that handles Perturb-seq data. Much of
the evidence below comes from the code of sceptre and pysceptre, because
those are the tools whose needs were read closely; it is evidence about what
an analysis needs from a data file, not a requirement to use them. Analyses
of a data file are stored apart from it; the [sceptre analysis
format](../sceptre-analysis-format/) is one way, and its rules are numbered A
and T.

**Status:** draft (format 0.1), November 2026. Every choice below is
implemented in the prototype beside this file. None has been reviewed;
`overview.html` is the introduction, and says what review is asked for. Until
the format is frozen as 1.0 it may change and no migration code is promised
(S10). [Section 6](#6-fixed-now-free-to-grow) says which parts would be
expensive to change after that and which can grow.

Each choice is labelled **Decided** (the maintainer has ruled on it) or
**Adopted** (the author's recommendation, taken up so that the draft is
complete, and more open to challenge).

**How to read this.** Sections 1 and 2 give the whole format in a page.
Sections 3 and 4 defend each choice in a fixed form: the decision, the
alternatives, why, and what it locks into files. Where this document and
`SPEC.md` disagree, `SPEC.md` is the rule and the disagreement is a mistake to
fix. Evidence is tagged:

- *(verified)*: read in the code, with the file named;
- *(measured)*: a test run for this document, listed in [Appendix A](#appendix-a-measurements);
- *(inferred)*: reasoning, not checked;
- *(not independently checked)*: taken from a source the author did not read.

**Precedents.** Each decision that has a counterpart in the IGVF pipeline's
output ends with an **IGVF** note saying what IGVF does and why this format
does or does not follow it. [Section 5](#5-compared-with-igvf-and-pertpy) puts
those notes in one table and adds pertpy. The IGVF facts come from
`docs/mudata_schema.md`, `docs/output.md` and the writer scripts in `bin/` on the
`dev` branch of `pinellolab/CRISPR_Pipeline` at commit `ead4db3` *(verified by
reading them)*. No IGVF output file was inspected, nor the Google Sheet and Doc
that repository links to.

The code facts come from read-throughs of the codebases: R sceptre on the
`sceptre3` branch and on `main`, the pysceptre fork, public pysceptre, and the
Nextflow pipeline.

---

## 1. The format at a glance

```
my_screen.h5mu                        THE DATA. Written once at import. Never modified.
                                      Modalities: grna, plus one per thing measured (rna, ...).
                                      Carries its own format version and a dataset id.
```

The data file carries its own format version, because it is read on its own.
Nothing an analysis produces is in it. An analysis is stored apart, names the
data file by its identifier, and has a format version of its own; the two can
change independently.

## 2. The eight principles

Every choice below follows from these. Where a choice needed a trade-off, the
earlier principle wins.

| # | Principle | Why | What follows |
|---|---|---|---|
| 1 | **The same file from R and Python.** | Screens are analysed in both languages, and what a file means should not depend on which language wrote or read it. | Only encodings that both languages read *and* write correctly today are allowed; that subset was measured, not assumed (the measurements are in Appendix A of this document). Both languages read the file with standard libraries, so the format has no reader of its own. |
| 2 | **Data only: nothing in the file depends on an analysis.** | One dataset is analysed many times, with different settings. The data is large and fixed at import; an analysis is small and gets redone. Anything in the file that depended on an analysis choice would fit one analysis only, or force the file to be rewritten. | The file holds no assignments, no quality-control flags and no grouping of gRNAs (S8.5, S12). Counts are stored per gRNA, and are never summed over the gRNAs that one vector delivers together (S8.9): combining them is a step of an analysis. There is no label for low or high multiplicity of infection, which is a choice about how to analyse a screen (S10). No modality is named for a role such as "response": which one an analysis tests is that analysis's choice (S2). Covariates that a user supplies for an analysis, such as a score, a cluster or a cell-type label that some method computed, are part of that analysis and are saved with it, not in the data file (S12.3). The file is written once and software that analyses it does not modify it (S12.2). |
| 3 | **It belongs to no tool.** | A file that many tools read and write cannot carry one tool's name or assumptions, and data outlives the tool that imported it. | The format is a set of names and rules on a standard MuData file. Its own entries sit at the top level of `uns`, under no tool's name (S10). Entries, modalities and columns that other tools add are left alone (S10.2, S2). The validator needs no analysis tool. |
| 4 | **Every cell and feature has an identifier.** | Other files refer to the cells and features of a data file: analyses, annotations, results. A row number means something only for one file in one order, and R counts from 1 where Python counts from 0. | Cell identifiers are unique, and feature identifiers are unique across the modalities of a file (S9.1, S9.2). Cells without barcodes are given identifiers, which are stored (S9.3). The gRNA table gives each gRNA's target, and its vector where there is one, by name. |
| 5 | **Existing conventions before invented ones.** | An existing format comes with readers, tools and users. An invented encoding needs code in two languages and has no other reader. | AnnData and MuData encodings are used as specified. The format adds names and rules, not encodings. Target positions follow the convention of GTF files (S8.10). |
| 6 | **Readers ignore what they do not know.** | The format will grow, and files outlive the code that wrote them. Older software has to keep working on newer files where it can, and fail cleanly where it cannot. | Columns, modalities and entries that do not change the meaning of what is already there can be added with a minor version (S10). Software refuses a file of a major version it does not know (S10.3). |
| 7 | **What follows from other stored values is not stored.** | Two copies of one fact can disagree, and the file then does not say which is right. | There are no logarithms or percentages of stored columns, and no per-cell top gRNA (S6). The exception is a value that nearly every analysis needs and that takes a full pass over the counts to recompute: `n_umis`, `n_nonzero` and `p_mito` may be stored, and a writer should store the first two (S6). |
| 8 | **Files are harder to change than code.** | Code is corrected in the next release. A file already written stays as it is for years. | Each choice is made on its merits and before the freeze. How existing code reads data is evidence when the pattern is sound, and never the only reason. What does not lock anything in is left open: compression is the writer's choice, because HDF5 decompresses transparently (S11). The section [Fixed now, free to grow](#6-fixed-now-free-to-grow) says which parts would be expensive to change. |

The [sceptre analysis format](../sceptre-analysis-format/RATIONALE.md#the-principles)
has principles of its own. Several are these, applied to an analysis folder.

---

## 3. Number and text types

One rule per kind of value, the same in both languages.

| Value | Type in the file | Why |
|---|---|---|
| Counts (`X/data`) | 32-bit signed integer | Counts are integers. One count never approaches 2³¹. R has no unsigned type. *(measured)* 33% smaller than doubles uncompressed. |
| The two index arrays of the sparse matrix (`X/indices`, `X/indptr`) | one signed integer type for both: 32-bit, or 64-bit past the limit | 32-bit while the stored values, the cells and the features each number at most 2,147,483,647; 64-bit otherwise. It is the type scipy chooses, so the standard Python tools read the counts with nothing done to them. See S4. |
| Other integers | 32-bit signed integer | R's integer is 32-bit. *(measured)* R's MuData reader returns a 64-bit column as all-missing. No value in the format needs more than 32 bits. |
| Real numbers | 64-bit float (double) | One type for every real number. Never 32-bit. |
| True/false in the data file | categorical with levels `FALSE`, `TRUE` | *(measured)* R writes a logical as an integer and reads a Python boolean as text. See S7. |
| Text | UTF-8 | AnnData's plain `string-array` only (S7). |
| Missing values | float: NaN. Categorical: code −1. | No integer or plain-text column in the data file may have missing values (S7). |

**Readers are lenient, writers are strict.** A reader accepts any integer or
float type for a number. A writer writes exactly the types above. That way a
file from another tool still opens, and files from our two writers are the same.

**IGVF.** Counts are stored in the smallest unsigned integer type that holds
them, usually 16-bit (`bin/count_matrix_utils.py`); see S4.

---

## 4. The data file

### S1. One MuData `.h5mu` file — *Adopted*

- **Decision.** The dataset is a single HDF5 file in the MuData layout.
- **Alternatives.** A Zarr store; two `.h5ad` files; the `ondisc` `.odm` format; a folder of Parquet and matrix files.
- **Why.** MuData is the scverse standard for one set of cells measured in several ways. Both languages open it with standard libraries *(measured)*. One file is easy to move and to publish. HDF5 allows one writer at a time, which costs nothing for a file written once.
- **Locks in.** HDF5. Zarr suits cloud object storage better; AnnData defines the same encodings for Zarr, so a Zarr variant could be added later without changing anything else in this document. Public pysceptre's roadmap asks for Zarr *(verified, `docs/ROADMAP.md`)*; I did not test R's Zarr support.
- **IGVF.** The same container: one MuData `.h5mu` file, `inference_mudata.h5mu`. The two formats agree here.

### S2. Modalities: `grna`, plus one per thing measured — *Decided*

- **Decision.**
  - The gRNA modality is named `grna`. This name is fixed.
  - Every other modality is named for what was measured: `rna` for RNA (recommended), and, for example, `atac` or `protein`. These names are free. A name starts with a letter and uses letters, digits and underscores, so that it can prefix a column name where software sets the cell tables side by side (S5).
  - A file has `grna` and at least one other modality. All modalities hold the same cells in the same order.
  - The file marks no modality as the one to analyse. Which modality an analysis tests is that analysis's choice, recorded with the analysis and not in the file.
  - A file may carry more modalities than any one use needs: a hashing modality, for instance. A tool that reads a file and writes it again keeps them all.
  - Every modality follows the rules for counts, cell tables and feature tables, whichever of them a piece of software uses *(Adopted)*. A file is then conforming or not on its own, and any modality can be analysed later.
- **Alternatives.**
  - Name the modality for its role: `response`.
  - IGVF's names, `gene` and `guide`.
  - A fixed vocabulary of assay names.
- **Why.**
  - Which modality is the response is an analysis choice. In an ECCITE-seq screen one analysis tests RNA and another tests protein, on the same data file. Naming a modality `response` would put an analysis choice in the data file, against principle 2.
  - `response` cannot name two modalities, and it does not say what a chromatin-only file contains.
  - No vocabulary is settled beyond `rna`: IGVF uses `gene`, `guide` and `hashing`; pertpy's one MuData example uses `rna`, `adt`, `hto` and `gdo` *(verified, pertpy tutorials at commit `4924c25`)*. So the format fixes only `grna` and leaves the rest free.
  - `grna` is kept because sceptre's vocabulary uses it throughout (`grna_id`, `grna_target`).
- **Locks in.** The name `grna`, and the rule that names are otherwise free.
- **IGVF.** Modalities `gene` and `guide`, with `hashing` when present. This format takes IGVF's idea of a file that may carry more than two modalities, and not its names. An importer from IGVF files renames the two modalities.

### S3. Counts: cells × features, compressed by feature, canonical — *Decided (orientation)*

- **Decision.** `X` holds raw counts as an AnnData `csc_matrix` of shape cells × features. Within a feature, cell positions are strictly increasing, with no repeats and no stored zeros.
- **Alternatives.** Compressed by cell (`csr_matrix`, the common default); dense.
- **Why.** One feature is one contiguous read, which is what fitting one gene needs, and testing for association goes gene by gene. Three codebases show it. Public pysceptre's disk-backed reader requires this layout *(verified, `scripts/sceptre_io.py`)*, the fork reads it directly and has to re-orient the other one *(verified, `execution/scverse_inputs.py`)*, and R's backed reader relies on it *(verified, `R/backed_counts.R`)*. IGVF's pipeline, which stores the other layout, splits its file into per-gene chunk files before fitting (below).
- **Locks in.** The orientation. A file from another tool is usually compressed by cell and is converted once at import.
- **IGVF.** Compressed by cell (`csr_matrix`), which is what the upstream mapping tools and scanpy produce *(verified, `bin/count_matrix_utils.py`)*. To fit a gene's model, the pipeline splits the whole file into per-gene chunk files first (`bin/chunk_mudata_sceptre.py`). Compression by feature makes that split unnecessary, so this format does not follow IGVF here.

### S4. Count type: 32-bit integers; one index type for both index arrays — *Decided (index type); otherwise Adopted*

- **Decision.** As in section 3.
- **Alternatives.**
  - Doubles: what R's MuData writes by default *(measured)*. 588 MB against 395 MB uncompressed on a synthetic 200,000 × 4,000 matrix *(measured)*, and it admits non-integer "counts".
  - The smallest unsigned type that fits: what public pysceptre does *(verified)*. 298 MB on the same matrix *(measured)*, but the type then differs from file to file and R has no unsigned integers.
- **Why.** A fixed type makes the two writers agree. With gzip the choice barely matters (98 MB against 103 MB, *measured*), so the smaller unsigned types buy little. R can write these types with six lines after the MuData writer *(measured)*.
- **Locks in.** Counts are integers. R sceptre does not check this today *(verified, `R/check_functions.R`)*; the fork does *(verified, `storage/count_identity.py`)*.
- **One level: a file conforms or it does not** *(Decided)*. A file whose counts or index arrays are in another integer or floating-point type, with whole values, does not conform, and a reader must still read it (S4.4). That is "readers lenient, writers strict" as one behaviour rule. *(measured)* Both readers read the five such conformance files to the same counts as the conforming file.
  - **Alternative: a middle level, "readable"**, for a file that does not conform and that a reader must accept, reported apart from files that are damaged. Rejected because no real file from another tool would land in it: an IGVF file also has other modality names, the other matrix layout and assignments as a layer, so it breaks many rules and needs an importer regardless. The validators return rule, place and message, and a file conforms when they return nothing.
- **The index type** *(Decided)*. Both index arrays of a modality are stored in one signed integer type: 32-bit while the stored values, the cells and the features each number at most 2,147,483,647, and 64-bit otherwise. This is scipy's own rule *(verified in scipy 1.18.1: `get_index_dtype`; a matrix with a dimension of 2³¹ − 1 gets 32-bit arrays and one of 2³¹ gets 64-bit)*.
  - **Alternative: 32-bit cell positions with 64-bit column pointers, always.** One fixed pair of types, but invented here, where scipy's rule is the existing convention. *(measured)* Mixed widths make one common scanpy function, `calculate_qc_metrics`, fail on a plain anndata read, because scipy rejects mixed widths in some operations, and a Python reader then needs a step of its own. The case for 64-bit pointers at every size is arithmetic only: the scale test's 699 million nonzeros are under the limit.
  - **Cost.** Past about 2.1 billion nonzeros both arrays are 64-bit, so a cell position costs 8 bytes, not 4: roughly 12 bytes per nonzero against 8, half as much again before compression, on the datasets where size matters most. And the type depends on the data. It is still determined by the data alone, so R and Python write the same types; both validators compute the type the size requires.
  - **Rejected.** Always 64-bit for both: one fixed type, but every file half as much again larger before compression.
  - **Wider than needed, or two widths, does not conform**, like counts in another type, and a reader still reads it (S4.4). The Python reader brings two widths to one for such files; the step does nothing to a conforming file.
  - **What was checked** *(measured, `verify/wide_index.*` and the conformance set)*. scanpy's `calculate_qc_metrics` runs on every example and conformance data file read with anndata alone, for files written by R and by Python. No test file can reach the limit, so the 64-bit branch is tested with the limit lowered, in both writers and both validators: each writer then stores both arrays as 64-bit, the two files hold the same arrays, and each validator passes them under the lowered limit and reports S4.2 and S4.3 under the real one. **64-bit cell positions in R:** the high-MOI example with both arrays 64-bit reads to the same counts through MuData's reader in memory, through its backed reader (one feature at a time and whole), and through HDF5Array, with no warning beyond those on the 32-bit file.
  - **What was not checked.** Values above 2,147,483,647 in either array. That needs a file of more than 2.1 billion nonzeros, or cells, and none was built. R holds such values as doubles when read with rhdf5, exactly up to 2⁵³; whether MuData's backed reader and HDF5Array do the same at that size is unmeasured. The test above shows the wider type is handled, not the larger values.
- **Flag.** R's in-memory sparse matrix cannot hold more than 2³¹ − 1 nonzeros. Past that size R must read disk-backed. That is a limit of R, not of the file.
- **IGVF.** The smallest unsigned type that fits, chosen from the data's maximum, as public pysceptre does. This is the closest call in this document: IGVF's rule saves about a quarter of the uncompressed size (298 MB against 395 MB above) and nothing meaningful with gzip. This format keeps one fixed type so that the type never depends on the data. Readers here accept IGVF's types as they are. IGVF's code independently works around the same scipy problem with mixed index types (`normalize_sparse_index_dtypes`).

### S5. Cell covariates: shared table for cell facts, modality tables for measurements — *Decided*

- **Decision.** Facts about the cell (batch, donor, anything the user supplies) are in MuData's shared cell table. Measurements of one modality (`n_umis`, `n_nonzero`, `p_mito`) are in that modality's cell table, under their plain names. A shared column is not named with a modality's name and an underscore, so that software which sets the tables side by side, each modality's columns under its name (`rna_n_umis`, `grna_n_umis`), gets distinct names. The sceptre analysis format builds a model's covariate table that way (T7 there).
- **Alternatives.** Everything on one modality's cell table; everything in the shared table.
- **Why.** A batch belongs to the cell, not to its RNA. With two measured modalities, "everything on one modality" has no answer for where batch goes. *(measured)* The shared table round-trips in both directions between R and Python.
- **Column names** *(Adopted)*. A name starts with a letter and uses letters, digits, underscores and dots. Dots are allowed because they are common: Seurat's defaults, which pertpy's example dataset carries, include `S.Score` and `percent.mito`. Without this rule such a column could sit in a file and never be used in a formula. *(measured)* R takes a dotted name in a formula as it is.
- **Locks in.** That a shared column does not start with a modality's name and an underscore.
- **IGVF.** Per-cell measurements are on the modality's cell table, as here, under scanpy-derived names (`total_gene_umis`, `num_expressed_genes`, `total_guide_umis`, `percent_mito`); some are also copied to the shared table. IGVF's `percent_mito` is a percentage and this format's `p_mito` is a proportion, so an importer must convert it.

### S6. Which computed columns are stored — *Adopted*

- **Decision.** A modality's cell table may have `n_umis` and `n_nonzero`, and a writer should store them. `p_mito` is optional too, for a modality that has mitochondrial genes. Not stored: each cell's top gRNA (`max_grna`, `max_grna_frac_umis`).
- **Why optional** *(Decided)*. Not every tool needs them, and a file that comes straight from another pipeline should conform without them. A tool that needs them and finds them absent computes them in one pass over the counts; the sceptre analysis format lets an analysis carry them itself. They are defined here so that every file that has them means the same thing by them.
- **Why these and no others.** They are used by sceptre's default formula and by cell QC in nearly every analysis, and recomputing them is a full pass over the counts (principle 7). `p_mito` also records an import decision: which genes counted as mitochondrial. The top gRNA is read by one assignment method only *(verified, `R/assign_grna_functions.R`)* and follows from the gRNA counts, which are small. The prototype rebuilds it and matches R's values exactly *(measured)*.
- **Alternatives.** scanpy's names (`total_counts`, `n_genes_by_counts`). sceptre's names are kept because users write them in formulas *(verified, default formula builder in `R/aux_functions.R`)*.
- **Not checked.** A validator checks that the columns, where present, are numbers, not that they equal what the counts give.
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

- **Alternatives.** Allow every AnnData type and fix the libraries. That is the right long-term direction and the fixes belong upstream (see the suggested changes in [section 6](#6-fixed-now-free-to-grow)). Until they land, a wider set would break principle 1.
- **Factors.** A categorical's stored category order is the level order. This is AnnData's convention; nothing is invented. Because every text covariate is a categorical, R and Python cannot disagree about level order (R sorts by locale, Python by code point).
- **Locks in.** No boolean or 64-bit columns in the data file. Lifting this later is allowed by principle 6.
- **IGVF.** The file is written by Python only and uses booleans (`targeting`, `mt`) and whatever integer width pandas chose. R reads only the count matrices, the MOI and the pairs from it (`bin/inference_sceptre.R`), so the encodings R misreads never matter there. This format's rule is stricter because both languages must read and write every column. IGVF's writers also replace missing strings with a token and store low-cardinality strings as categoricals (`make_h5mu_safe_dataframe`), for the reason S7 gives.

### S8. The gRNA table — *Adopted*

- **Decision.** One row per gRNA in the gRNA modality's feature table: `grna_target` (required), and optionally `target_chr`, `target_start`, `target_end`, the position of the *target*. A non-targeting gRNA has `grna_target` equal to `non-targeting`.
- **Why.** These are facts about the gRNA, known at import. `non-targeting` is sceptre's long-standing convention *(verified, `R/assign_grna_functions.R`)*; a second marker could contradict it.
- **Why `target_chr` and not `chr`.** sceptre's gRNA table calls these columns `chr`, `start` and `end` and uses them as the target's position when it builds cis pairs. The bare names do not say whether they locate the gRNA or its target, and IGVF's files carry both. The file names say which; the packages keep their own names and the prototype translates.
- **Not stored.** `grna_group`. It depends on a choice an analysis makes, how gRNAs are combined *(verified, `R/preprocessing_functions.R`)*, so it is not data (principle 2).
- **A gRNA under two targets is refused.** It cannot be one row, and a writer stops with a message. R accepts such tables today, with a warning on the development branch *(verified, `R/accessors.R`)*, and public pysceptre counts them and keeps the full table separately *(verified, `scripts/sceptre_io.py`)*. What pooling means for such a gRNA is a scientific question first; a separate gRNA-to-target table can be added later without changing anything else ([section 6](#6-fixed-now-free-to-grow)).
- **The kind of gRNA: `grna_type`** *(Decided)*.
  - *The column.* An optional text column of the gRNA table with four values: `targeting`, `non-targeting`, `safe-targeting`, `negative control`. These are the values of the `type` column in the IGVF pipeline's guide table *(verified, `docs/mudata_schema.md` on its `dev` branch)*.
  - *Controls keep one marker.* A gRNA of any of the three control kinds has `grna_target` equal to `non-targeting` and no target position. That is how the IGVF pipeline treats them: it takes a guide as a control if its type is one of the three or its target name is `non-targeting`, and its guide documentation says controls should carry no invented positions *(verified, `bin/intended_target_key_utils.py` and `docs/guide-metadata.md`)*. So whether a gRNA is a control is still decided by `grna_target` alone (S8.3), and the column adds one thing: which kind of control it is. A safe-targeting gRNA does cut somewhere; the format records that it is a control and not where it cuts.
  - *The two must agree* (S8.16). A file in which a control kind has a real target, or a `targeting` gRNA has the target `non-targeting`, says two different things about one gRNA.
  - *Why optional.* Most existing data says only "non-targeting", including libraries whose controls are in fact safe-targeting, so a required column would make an importer state a guess as fact. For a targeting gRNA the value follows from the target, which principle 7 says not to store.
  - *Not added:* columns for the gRNA's own position, spacer, PAM and strand. They are planned ([section 6](#6-fixed-now-free-to-grow)); a file that carries them keeps IGVF's names.
- **Several gRNAs on one vector: `vector_id`** *(Decided)*.
  - *The rule.* A vector has one target. Several gRNAs delivered together against the same element are in scope. A vector whose gRNAs have different targets is a combinatorial screen and is out of scope: sceptre is not equipped to analyse those, and the format should be future-proof without accommodating every possibility. This is the rule sceptre already enforces, that `vector_id` is nested within `grna_target`.
  - *The column.* An optional text column `vector_id` in the gRNA table, beside `grna_target`. Which gRNA sits on which vector is fixed at cloning and known at import, so it belongs in the data file by the format's dividing rule. sceptre keeps it today as an optional column of its gRNA table under the same name, so this continues a convention.
  - *Counts stay per gRNA.* A writer never sums over a vector's gRNAs. sceptre sums them at import today, which destroys the one piece of information needed to tell a cell carrying two gRNAs of one vector, as intended, from a cell carrying gRNAs of two vectors, a recombination product. Dual-guide libraries report such recombinants at roughly a quarter of reads *(Replogle et al. 2022; not independently checked)*. Summing within a vector also raises sensitivity when per-gRNA counts are low, and a file without `vector_id` cannot reproduce it.
  - *Why a reader must refuse two targets on a vector (S8.7).* It keeps the door open cheaply. A reader must either read a file correctly or detect that it cannot (S10). With the refusal stated now, allowing combinatorial screens later is a minor version change, because today's readers reject those files cleanly. Without it, the same change would be major.
  - *Why a vector may not be named like another vector's gRNA (S8.8).* sceptre's collapse renames `vector_id` to `grna_id`, so a stray match would silently merge two different things. A vector with one gRNA may still be named after it.
  - *What this covers without anything more.* Under one target per vector, barcode-based designs need no new structure: a barcode standing for several gRNAs with one target is a feature with counts and a target, already a row of the gRNA table. *(from the papers as reported to me; not independently checked)* Adamson et al. 2016 and Norman et al. 2019 identify a multi-guide vector by one barcode; Replogle et al. 2022 count each guide separately by direct capture; CaRPool-seq, CROPseq-multi, Cas12a arrays and VIPerturb-seq each use a barcode per array or per guide. Both schemes fit.
  - *Not added, on purpose:* a vector-to-gRNA table (needed only when a vector carries several targets), targets as lists, a position within the vector, a field for the detection scheme, and a field for the target molecule.
  - *Left to analyses.* Because the data file holds per-gRNA counts and is written once, summing within a vector is a step of an analysis. Nothing in the data file depends on how an analysis does it. What it implies for sceptre's analyses is in [the sceptre analysis format's rationale](../sceptre-analysis-format/RATIONALE.md#t1-assignmentassignmentsparquet--decided).
- **No entry for the target molecule** *(Decided)*. Target coordinates are genomic positions; a screen that knocks down RNA names its transcript in `grna_target` and may leave coordinates out. Nothing in an analysis depends on the molecule: counts, assignments, targets and tests are the same, and the coordinates are used only when the user builds cis pairs. The optional `perturbation_type` (S10) records which kind of screen it was, and its knockdown value is the RNA-targeting case.
- **The name `grna` is a slight wart.** The modality holds whatever identifies the perturbation, which for some designs is a barcode and could in principle be something other than a gRNA. It is kept because sceptre's vocabulary uses it throughout and another name buys nothing within this scope.
- **Positions are 1-based with both ends included** *(Adopted)*. `target_start` is the first base of the target and `target_end` the last. One fixed convention for the whole format, with no entry saying which is in use.
  - *sceptre's own data are already in this convention* *(verified in this folder's example)*. In sceptre's high-MOI example every target has `end − start = 49`, and a target's second gRNA starts at the first one's end plus 1 (23567064 to 23567113, then 23567114 to 23567163). Adjacent 50-base tiles are consistent only with closed intervals. The gene positions sceptre pairs them with come from a GTF file, which is 1-based with both ends included.
  - *Nothing has to be converted.* Files written from existing sceptre objects are already right, and target and gene positions are in one system. It is also the convention of Ensembl and of Bioconductor's `GRanges`, which the planned genome-build check would compare.
  - *Why it was never pinned down.* sceptre builds cis pairs from the midpoint of a target and a distance threshold of hundreds of kilobases, so being one base off changes nothing. It matters for exchanging files, not for sceptre's results.
  - *The cost.* Python tools and BED files usually count from 0 and exclude the end. A fixed, stated rule and one sentence ("add 1 to a BED start at import") handle it.
  - *Checked.* Both validators report a start below 1, an end before its start, and a gRNA with some of the three position columns and not all.
  - *Chromosome names* are text; the format does not choose between `chr1` and `1`. sceptre's own data use the `chr` prefix.
  - *IGVF* does not fix a convention: its schema has a per-guide column `coordinate_system`, with zero-based half-open and one-based inclusive given as examples. That passes the question to every reader, which this format avoids. crisprVerse uses `GRanges`, but its guide coordinate is a single position at the PAM, a different thing from the targeted element *(not independently checked)*.
- **IGVF.** Much richer, and better than this format today. Each gRNA has a `targeting` flag, a `type` (targeting, non-targeting, safe-targeting, negative control), its `spacer` sequence, its own position (`guide_chr`, `guide_start`, `guide_end`, `strand`, `pam`), its target's name and position (`intended_target_name`, `intended_target_chr`, `intended_target_start`, `intended_target_end`) and a collision-safe target key; source columns are passed through. This format takes the distinction between gRNA position and target position now, through the column names. The rest are optional columns that can be added without a format change and are listed in [section 6](#6-fixed-now-free-to-grow). A file may already carry them.

### S9. Identifiers — *Adopted*

- **Decision.** Every cell and every feature has a unique text identifier. Feature identifiers are unique across all modalities of a file, so no table needs a modality column.
- **Why.** Principle 4. Cell identifiers are also load-bearing for results: R's current random-number protocol keys each cell's draws on its identifier *(verified, `R/rng_v2.R`)*.
- **Screens without barcodes.** The writer assigns the positions `"1"`, `"2"`, ... as identifiers and stores them, as R does in memory today *(verified)*. Because identifiers affect results, they are fixed in the file, so both languages see the same ones.
- **No line feeds.** An identifier must not contain a line feed, so that a list of identifiers joined by line feeds is unambiguous. The sceptre analysis format hashes such a list to link an analysis to its data (A1 there).
- **IGVF.** The same: cells, genes and guides are named by identifier (the cell barcode with a measurement-set suffix, the Ensembl gene identifier without its version, the guide name), and tables are keyed by them.

### S10. File-level metadata: a native record, not JSON text — *Adopted*

- **Decision.** Scalar text entries at the top level of `uns`: `perturb_seq_format` (the version, such as `0.1`) and `dataset_id` (a random identifier assigned when the file is written, by which an analysis names its data file), and the optional entries below.
- **At the top level, under no tool's name** *(Decided)*. These entries are facts about the dataset, and other tools should be able to read and write them. Under a key named for one tool, another tool would have to write into that tool's space or invent its own key for the same facts. Dataset facts at the top level of `uns` is what CELLxGENE's schema does *(not independently checked)*. The two entries that are the format's own bookkeeping carry distinctive names, because every tool shares this level: the version entry names the format, and `dataset_id` and `title` are not `id` and `name`.
- **No `moi`** *(Decided)*. Whether a screen is treated as low or high multiplicity of infection is a choice about how to analyse it. The same data can reasonably be analysed either way, and real screens are not all one or the other. In the sceptre analysis format it is a setting of the analysis's parameters record (A3.6 there), not an entry of the data file.
- **A two-part version** *(Decided)*. The version is `MAJOR.MINOR`. A minor increase adds something an older reader can safely ignore; a major increase is a change after which an older reader could misread a file. A reader reads any file of the major version it knows and refuses another.
  - *Why.* With one number, every addition either goes unrecorded or forces every reader to refuse the file. Two parts let a file say what it uses while old readers keep working. The stamp is in every file, so this had to be settled before the freeze.
  - *Precedent.* Apache Arrow versions its format this way and states the rule as: an older reader must either read newer data or detect that it cannot *(Arrow's `Versioning.rst` and `Changing.rst`; not independently checked)*.
  - *What it corrects.* Principle 6, "readers ignore what they do not know", is safe only for additions that do not change the meaning of what was already there. [Section 6](#6-fixed-now-free-to-grow) marks each planned addition as minor or major.
  - *Stored as text*, the same in HDF5 and in JSON, and compared as two integers, never as text.
- **Alternatives.** One JSON string (the first draft of this format); HDF5 attributes.
- **Why.** Principle 5. *(measured)* Both languages write and read a nested record with true scalars; R needs one attribute set on each value to write a scalar. A JSON string needed a low-level HDF5 write in R and is opaque to HDF5 tools.
- **No list of text columns.** The first draft kept one, to mark text stored as a categorical. Under S7 every covariate that was text *is* a categorical, and the few text columns of the feature tables are named in `SPEC.md`.
- **The version here is the data file's own.** An analysis has a separate one in `analysis.json`.
- **What kind of screen it was: `perturbation_type`, `organism`, `genome_build`** *(Decided)*. Three text entries. No reader's behaviour depends on them. The first two are optional; `genome_build` is required in a file whose gRNA table has positions (S8.14), and optional otherwise.
  - *Why they are in now, though nothing reads them.* "Nothing reads it, so it can be added later" is the right test for a value that follows from other stored values (principle 7). It is the wrong test for provenance. What kind of CRISPR perturbation a screen used is known at import, cannot be worked out from anything in the file, and is lost for good if not written down. It is also needed to read a results table at all: a log2 fold change of −1 is the expected result under CRISPRi and a surprising one under CRISPRa. The dividing rule already settles where it goes: the data file holds what is known at import.
  - *The vocabulary* for `perturbation_type` is crisprVerse's five kinds: `CRISPRko`, `CRISPRi`, `CRISPRa`, `CRISPRbe`, `CRISPRkd` *(crisprDesign README, not independently checked; that `CRISPRkd` means knockdown of RNA with an RNA-targeting enzyme should be confirmed)*. CELLxGENE's schema, the EMBL-EBI Perturbation Catalogue and scPerturb record the same concept *(not independently checked)*. A new kind is a minor addition.
  - *`genome_build` is required with positions.* Positions are used to pair targets with nearby genes. If the two sets of positions are in different builds the pairs are wrong and nothing fails, and a file without the build gives no way to notice. The value is free text, such as `GRCh38`: the format keeps no list of builds, so software comparing two builds has to know that `GRCh38` and `hg38` name the same one. The cost is that positions of unknown build cannot be stored.
  - *It replaces a separate entry for the target molecule.* The knockdown kind is the RNA-targeting case, so one entry covers both distinctions.
- **IGVF.** `moi` (`high` or `low`) and `capture_method` are in the guide modality's `uns`, each as a one-element string array; this format has no `moi` entry, because it is a choice of an analysis. There is no format version and no dataset identifier; a file is identified by a SHA-256 of all its bytes in a separate manifest. This format adds the version and the identifier because files outlive the code that wrote them.

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

## 5. Compared with IGVF and pertpy

### IGVF

The IGVF pipeline's output is the closest existing format: one `.h5mu` file per
run, `inference_mudata.h5mu`, with result tables beside it. Sources are given at
the top of this document.

**Where the two agree**

| Area | Both |
|---|---|
| Container | MuData `.h5mu`; cells × features; raw integer counts; sparse |
| Identifiers | Cells and features are named by identifier, not position |
| Library problems met | Mixed index types in scipy; missing strings in HDF5 |

**Where they differ**

| Area | IGVF | This format | Which is better, and why | See |
|---|---|---|---|---|
| Data and analysis | One file holds both | Separate | This format: one dataset, many analyses. IGVF's own code shows the cost of patching results into the data file | S12 |
| Cells | Only cells that passed QC | All cells | This format: which cells are kept is an analysis choice | S12 |
| Assignments | 0/1 matrix layer in the data file | Not in the data file | This format, for the same reason | S12 |
| R and Python | Python writes; R reads the matrices, MOI and pairs | Both read and write everything | This format: it is the first requirement | S7 |
| Format number, dataset id | Neither | Both | This format | S10 |
| MOI | `high` or `low`, in the data file | Not in the data file | This format: it is a choice about how to analyse a screen, and one dataset can be analysed either way | S10 |
| Matrix layout | Compressed by cell | Compressed by feature | This format: one gene is one read, with no chunk files | S3 |
| Count type | Smallest unsigned type that fits | Fixed 32-bit | A close call; kept for a type that never depends on the data | S4 |
| gRNA annotation | Type, spacer, gRNA and target positions, target key | Target, target position, and optionally the type, with IGVF's four values | **IGVF** for the spacer and the gRNA's own position, which are planned | S8, section 6 |
| Gene annotation | Symbol, position, mitochondrial and ribosomal flags | Unspecified | **IGVF.** Planned | section 6 |
| Reference genome | Build per gRNA; annotation file hash in the manifest | One build for the file, required when the gRNA table has positions; no annotation hash | **IGVF** for the annotation hash, which is planned | S8, S10 |
| Other modalities | Optional `hashing` | Any number, under free names | Taken from IGVF | S2 |
| Names | `gene`, `guide`, `intended_target_name` | `rna`, `grna`, `grna_target` | Neither is a standard. An importer maps one to the other | S2 |

In the rows on data and analysis, cells, assignments and MOI, this format
keeps out of the data file what IGVF keeps in it; S12 gives the reason. Where
sceptre's analyses keep them is compared in [the sceptre analysis format's
rationale](../sceptre-analysis-format/RATIONALE.md#igvf).

**What this format took from IGVF:** a file may carry more than two
modalities; the gRNA table says whose position its coordinates are; and the
list of annotation to add in [section 6](#6-fixed-now-free-to-grow).

**What an IGVF importer has to do:** rename the modalities and columns,
re-orient the matrices, convert `percent_mito` to a proportion, and take the
assignment layer, the pairs and the results out of the file, into an analysis.
The cells IGVF removed in QC cannot be recovered from its file.

### pertpy

pertpy defines no file format. It works on in-memory AnnData objects, and the
user passes the name of every column *(verified, source and tutorials at commit
`4924c25`; no pertpy dataset was downloaded)*. Its habits:

| Area | pertpy | This format |
|---|---|---|
| Container | 54 of its 55 bundled datasets are single `.h5ad` files; one is MuData with `rna`, `adt`, `hto`, `gdo` | MuData |
| gRNA-to-target mapping | A per-cell column | A property of the gRNA |

The one thing pertpy settles for this format is
that `rna` is the shared name for the RNA modality; its `adt` and `gdo` are
Seurat's assay names and appear in one dataset. How pertpy holds assignments
and results is compared in [the sceptre analysis format's
rationale](../sceptre-analysis-format/RATIONALE.md#pertpy).

---

## 6. Fixed now, free to grow

### What would be expensive to change after the freeze

| Kind | What |
|---|---|
| Structure | Data separate from analysis; the matrix layout; the number types |
| Names | `grna`; the column and entry names in section 4 |
| Removals | A column or entry that readers rely on cannot be taken away |

Everything else can be added with a minor version, because readers ignore what
they do not know (principle 6), or loosened, because old files stay valid. The
table of planned additions below says which additions are not safe to ignore.

### To settle before the freeze

- **Review** (`overview.html` says where to look).

### Planned additions to the format

None of these is in format 0.1. The last column says what kind of version
change each would be once the format is frozen (S10): **minor** if an older
reader can safely ignore it, **major** if it could misread the file.

| Addition | Where | Why wait | Version |
|---|---|---|---|
| Spacer sequence and PAM | gRNA table | Known at import and not derivable, like the entries of S10, but they serve guide design and reanalysis, not reading results | minor |
| The gRNA's own position and strand | gRNA table | sceptre does not read it | minor |
| Gene symbol, gene position, mitochondrial flag | Feature table of the measured modality | The flag would record which genes `p_mito` counted; it needs the importers to keep it | minor |
| A hash of the gene annotation used | File metadata | The genome build is recorded now (S10); the annotation is not | minor |
| How the gRNAs were captured (IGVF's `capture_method`) | File metadata | sceptre does not use it | minor |
| A declared type for each modality, beside its name | File metadata | No dataset needs it yet | minor |
| A gRNA-to-target table, for a gRNA with two targets | Data file | A scientific question first (S8) | major: a file using it has gRNAs an old reader would assign to one target |
| A screen that combines targets on one vector (combinatorial perturbations) | Data file, S8.7 | Out of scope: sceptre cannot analyse them | minor, because today's readers refuse such files |
| Which software and version wrote the data file | File metadata | The data file does not say this. The format and the software that writes it are versioned separately, so this is information, not a rule | minor |
| True booleans and nullable integers in the data file | S7 | Needs the R MuData and pandas fixes upstream | major: an old R reader misreads them |
| A Zarr variant | S1 | Not tested | a second container, not a version of this one |

### Suggested changes to what the packages produce

These are changes to sceptre and pysceptre, not to the format. None is made
here. Each would remove a workaround in the format or in the code that reads
and writes it. The suggestions about analyses are in [the sceptre analysis
format's rationale](../sceptre-analysis-format/RATIONALE.md#suggested-changes-to-what-the-packages-produce).

| Change | What it would remove | Evidence |
|---|---|---|
| Check at import that counts are integers | Non-integer "counts" reaching a file (S4) | `R/check_functions.R` |
| Rename `chr`, `start`, `end` in the gRNA table, or accept both | The same (S8) | — |
| Decide what a gRNA under two targets means | The refusal (S8) | R warns, pysceptre counts |
| Report the R MuData and pandas problems upstream | The restrictions of S7 | Appendix A |

---

## 7. How the format is defined and checked

### The pieces

| Piece | Where | What it is for |
|---|---|---|
| The rules | `SPEC.md` | Numbered statements, each a file rule, a behaviour rule or guidance, with a scope, a glossary and the library versions it is defined against. A validator's report cites them (`S4.1`) |
| The introduction | `overview.html` | One page: the principles, what is in a file, the choices that are expensive to change with their reasons in short, the questions that need a view, and what has been tested |
| The reasons | this document | Why each rule is what it is |
| Data-file definition | `spec/screen.json` | The structure of a data file as data: required groups, number types, encodings, allowed column kinds, required columns and metadata. Both validators read it |
| Conformance files | `conformance/`, built by `verify/make_conformance.py` | 10 conforming files, written by R and by Python, and 43 that each break a rule, with the rules a validator must report. Five of the 41 store counts or index arrays in another number type, which a reader must still read (S4.4) |
| Validator | `validate_screen()` in `format/perturb_seq_data.R` and `format/perturb_seq_data.py` | Checks a data file without building a package object, and reports every problem, not only the first |

**Why define it this way.** Writing the schemas exposes ambiguities in the
prose while they are still free to fix. Reviewers get something exact to
disagree with. And people can work with the format through AI assistants,
which need rules a machine can read and a validator to check their own output
against.

**Three kinds of statement.** Facts about files, duties of software and advice
are kept apart, because a MUST that no program can test is not a rule. Every
statement is a **file rule** (true of every conforming file, tested by the
validators), a **behaviour rule** (what a reader or writer does, tested by the
conformance run or the verification suite) or **guidance**
(advice or a definition, binding nobody). `SPEC.md` gives the counts and names
the rules that nothing tests.

**What the specification leaves out.** Left for later: tighter wording of
individual rules, and an implementation by someone who has only the
specification.

### What counts as evidence

1. **There is no canonical writer.** The standard writers (R's MuData,
   Python's mudata and anndata) do not produce conforming
   files as they come. So the format is defined by `SPEC.md` and `spec/`, and
   the validators are the arbiter. Any code whose files pass the validators and
   read back correctly in the standard libraries is a conforming writer.
2. **Writers are custom; readers can be generic, plus a thin layer.** Each
   conforming writer is a wrapper around the standard writer. R needs the
   most: counts written by hand, and column conversions. Python needs a
   setting, types set before writing, and the same conversions. The standard
   libraries open a conforming file and return the right values. What they do
   not do is interpret it: the version check, the refusals the specification
   requires, true/false arriving as a categorical, and derived values.
   That interpretation is each package's own, so the format has no reader.
   The reason is that the format allows only what both languages read
   correctly today and makes the writers stay inside that set. The burden
   falls on the writers, which are few, and not on readers, which include
   tools outside this project.
3. **The R-side restrictions come from one library.** Every restriction in S7
   traces to MuData's reader or writer, not to SingleCellExperiment: it cannot
   open pandas' nullable strings, reads a Python boolean as text and a 64-bit
   integer column as all missing, and writes a logical as an integer, missing
   integers in a form Python cannot open, and missing text as empty. If MuData
   is fixed or replaced, several restrictions could be lifted. None has been
   reported upstream.
4. **A round trip through one's own reader and writer is weak evidence.** It
   shows that the pair is self-consistent and cannot catch an error they share.
   Stronger evidence is independent: a file read by a reader written elsewhere;
   data read back and compared with its source, not with a rewrite; R's file
   compared with Python's.

### What is checked

**The verification suite** (`verify/run.sh`) needs no package that analyses
the data and runs on the committed files with standard libraries only.
*(measured, on macOS, Linux and Windows)* Each group of checks runs in R and in
Python, all passing:

- the validators and readers on the conformance set. Both validators accept
  all 10 conforming files and report exactly the expected rules for all 43
  broken ones, the same in both languages. 36 rules are exercised by at least
  one broken file. One more is checked with no file: S11.1 (no filter other
  than gzip), because R cannot read a file compressed otherwise at all;
- an encoding matrix: each kind of value written by each language's standard
  libraries and read by the other's;
- readers written elsewhere: anndata, scanpy, HDF5Array;
- the 64-bit branch of the index type, with the limit lowered (S4);
- what a writer does beyond producing a conforming file: a copy keeps the
  dataset's identifier and every modality, and past the index limit
  (lowered) both index arrays are 64-bit.

The two temporary writers are compared with each other, on seeded random
content, in the sceptre analysis format's folder, because that check writes
an analysis as well.

**The overview** (`overview.qmd`, rendered to `overview.html`) is the
introduction, not a check: it shows the example data file as an annotated
tree, with descriptions read from `spec/`, and reports what the verification
suite found when it was rendered. It renders R and Python in one process
*(measured, macOS only)*.

**What these cannot show.**

- *Independence of authorship.* The validators and both writers come from one
  author. The readers written elsewhere reduce that weakness and only an
  outside implementer removes it.

**The example files are fixed.** Nothing here regenerates them. A format
change converts them in place, with code kept beside the change, which is
also a direct test that old files can be carried forward.

### Taken from Apache Arrow's practice

- **A two-part version** (S10).
- **The forward-compatibility rule**: a reader either reads a file correctly or
  detects that it cannot. It is at the top of `SPEC.md`.
- **The change rule**: nothing enters the specification until both libraries
  implement it and a conformance file exercises it.
- **A table of what each library supports** (`SPEC.md`, section 5).

Not taken: a voting process, a separate repository of test files. An archive
of files written by each released version is worth starting at the freeze, so
that later readers can be tested against them.

### Where the code lives

- **The prototype in `format/` is a demonstration, not a reference.** It lives
  outside the packages because the format can then change in one edit, with no
  package release.
- **Its writers move into sceptre and pysceptre**, each taking its package's
  own object. The data file has no reader here: it is a standard `.h5mu`, and
  each package builds its object from what the standard reader returns. The copy here is retired once both packages pass the conformance
  suite.
- **sceptre and pysceptre are then the two reference implementations**, which
  fits the change rule above.
- **What stays with the specification:** `spec/`, the conformance files with
  their expected results, the validators, and a harness that runs any reader
  or writer against them (`verify/README.md`).

**Precedents** *(from each project's own repository; not independently
checked)*. In all seven, the reference implementations are the real libraries,
and what sits beside the specification is neutral: test files, expected
values, a harness, sometimes a validator. None keeps a separate,
non-production reference reader and writer once real libraries exist. Seven is
a small sample.

| Project | What it does |
|---|---|
| Apache Parquet | A format change needs two implementations that interoperate, one of them the Java library; a prototype is allowed while a change is discussed. A separate repository holds test files |
| Apache Arrow | At least two reference implementations, each inside a complete library, with integration tests. Test data is described in a neutral JSON form, and every producer is tested against every consumer |
| CommonMark | Two named reference implementations, which are also the libraries people use; the examples in the specification are the conformance tests, with a runner that takes any program |
| JSON Schema, TOML | No reference implementation; a test suite that does not depend on a language |
| BIDS | No reference writer; one official validator |
| SOMA | A specification and its first implementation, with Python and R interfaces |
| GeoParquet | A specification, a folder of examples, and a list of implementations |

Two things here go beyond those sources and are a judgement: the harness takes
any reader or writer as an input, so that a package can be tested without the
suite importing it; and there is one validator per language, kept in agreement
by the conformance set.

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

**Data-file encodings:** the table in S7, plus:

| Test | Result |
|---|---|
| Python writes 32-bit counts and 64-bit column pointers | R reads it in memory and disk-backed |
| R writes the same types (six extra lines) | Python reads it with those types |
| Nested metadata record with scalars | both directions correct; R needs `encoding-scalar` set to write a scalar |
| Shared cell table | both directions correct |
| gzip-compressed counts | R reads them |
| Python nullable boolean and nullable integer columns | R reads both correctly; R cannot write either in a form Python reads |
| A sparse matrix with 64-bit pointers and 32-bit indices in scipy | conversion to coordinate form fails until both are one type |

Not measured: Zarr; other HDF5 filters; a second platform.

### Scale test on real data

One run on the Gasperini et al. (2019) at-scale screen:
207,324 cells, 13,135 genes, 13,077 gRNAs, 699 million nonzero RNA counts. One
cluster node, 16 cores, 90 GB of memory. R 4.5.2 with MuData 1.14.0,
rhdf5 2.54.1 and HDF5Array 1.38.0; Python with the versions above. The
recorded results are in `scale-test/`. The files of this run stored 64-bit
column pointers beside 32-bit cell positions; under S4 both arrays of this
dataset are 32-bit, and the run was not repeated.

**What was checked, and held:**

| Check | Result |
|---|---|
| The file R wrote against the file Python wrote | Every count array identical; every cell and feature table identical; same number types |
| Both validators on both files | Conforming |
| Read back in memory, R and Python | Same shape, nonzeros and count sums as the source |
| 200 random genes read from disk one at a time, R (rhdf5) and Python (h5py) | Equal to the in-memory values |
| MuData's backed reader, one feature at a time | Sums identical to direct HDF5 reads; the 64-bit column pointers are read correctly |

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
| Open disk-backed with MuData and read 200 features | — | 5 s to open, 11 s to read, under 1 GB |

An assignment table of real size was written and read in the same run; [the
sceptre analysis format's rationale](../sceptre-analysis-format/RATIONALE.md#appendix-a-measurements) reports it.

**What the run showed about writers and validators.**

- *Dead space.* A writer that stores the counts and then rewrites the arrays in
  their final types leaves the first copy as dead space, because HDF5 does not
  give it back: 14.3 GB against 5.7 GB for this dataset, and 106 s against
  12 s to write in Python (399 s against 35 s in R). A writer stores each
  array once, in its final type.
- *Orientation in R.* HDF5Array presents the stored cells × features matrix as
  features × cells.
- *The validators are not built for size.* They read whole arrays into memory:
  about 29 GB in Python and 21 GB in R for this file. That is acceptable for a
  check run once per dataset and would need chunking for larger ones.

**What this does not cover.** A second dataset. A dataset past 2³¹ − 1
nonzeros, where both index arrays are 64-bit: this one has 699 million. R's
in-memory matrix would not fit such data.

## Appendix B. What the prototype implements

Everything in sections 3 and 4, and the validator of section 7. The code is
`format/perturb_seq_data.R` and `format/perturb_seq_data.py`; `example/` holds
files in the format. Whether the files hold everything an analysis needs was
checked once with both packages, together with the analysis format; [its
rationale](../sceptre-analysis-format/RATIONALE.md#appendix-b-what-the-prototype-implements) records that.

**Limits that remain.** These are limits of the prototype or of what was
tested, not open choices:

- sceptre and pysceptre call the gRNA table's coordinates `chr`, `start`,
  `end`; a package's reader and writer translate.
- A value may come back equal but not identical: an integer column with
  missing values is stored as doubles.
- Refused with a message: a gRNA under two targets.
- Tested once at scale (Appendix A): 207,324 cells. The suite passes on Linux
  as on macOS. Not tested: Zarr, the R readers anndataR and zellkonverter, and
  a file whose second measured modality is real data (the conformance file
  with a `hashing` modality holds a synthetic matrix).

## Appendix C. The examples

| Dataset | What it is |
|---|---|
| `example/highmoi` | sceptre's built-in high-MOI data: 500 cells, 100 genes, 60 gRNAs |
| `example/lowmoi` | sceptre's built-in low-MOI data: 1,000 cells, 100 genes, 50 gRNAs |
| `example/types` | 6 cells, synthetic: a data file with a factor whose levels are not alphabetical, a text covariate, a true/false covariate, a covariate with missing values, and a gRNA table with missing text and a missing integer |

These files are fixed: nothing here regenerates them. Analyses of the first
two are the examples of the sceptre analysis format.
