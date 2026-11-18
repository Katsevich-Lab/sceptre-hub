# Reviewing the sceptre file format: one page

The specification and its rationale are about 1,400 lines. This page says where
review time is best spent. It should take ten minutes to read.

## The format in ten lines

1. A dataset is **one MuData `.h5mu` file**. It is written once and never modified.
2. It holds raw counts for the gRNAs (modality `grna`) and for whatever else was measured (`rna`, ...), over the same cells.
3. Counts are stored by feature, so one gene is one contiguous read.
4. An **analysis is a folder**, separate from the data. One dataset can have any number of analyses.
5. `analysis.json` says which dataset the analysis belongs to and which modality it tests.
6. Each completed stage has **its own folder**: Parquet tables, then a small JSON record of its settings and files. The record is written last, so a stage without a record did not finish.
7. Every table row names cells, gRNAs and responses by identifier, never by position.
8. Nothing is stored that follows from something else stored.
9. R and Python read and write exactly the same files.
10. Every file carries a version, `MAJOR.MINOR`. Version 0 is a draft.

`SPEC.md`, section 1, has the layout as a diagram.

## Where to look hardest: what is expensive to change later

After the format is frozen, these cannot change without breaking existing
files. Disagreement is cheap now and costly then.

| Choice | Where it is argued |
|---|---|
| Data and analysis are separate files | `RATIONALE.md`, S12 |
| One record per stage, written last, instead of one settings file | A2 |
| The layout of the analysis folder: one folder per stage, holding its record and its tables | A2 |
| Fits as a cache outside the stages (`fits/response/`), shared by the testing stages and cleared only when QC is rerun | T6 |
| Tables keyed by identifier, sorted by key, not by p-value | T0, T5 |
| Counts stored by feature (column-compressed) | S3 |
| Modality names: `grna` fixed, others named for what was measured; the analysis names its response modality | S2, A1 |
| Covariate names carry the modality's name: `rna_n_umis`, not `response_n_umis` | S5 |
| Column and table names throughout (`grna_target`, `target_chr`, `pass_qc`, ...) | sections 4 to 6 |
| The stage names: `parameters`, `assignment`, `qc`, `calibration`, `power`, `discovery` | A2 |
| The design matrix is not stored: formula text plus a hash, additive formulas only | T7 |
| Assignments as a table of (gRNA, cell) rows, not a matrix | T1 |
| The two-part version and what counts as a breaking change | S10, and the top of `SPEC.md` |

## Questions that need your view

1. **Is the scope right?** The specification fixes the form of results and leaves their meaning (what a p-value tests) to sceptre's method documentation.
2. **Is anything missing that a real screen needs?** In particular for gRNA annotation: the format has a target, its position and an optional vector, and no gRNA type or spacer sequence yet (`RATIONALE.md`, section 10).
3. **Coordinates.** Target positions are 1-based with both ends included, as in GTF files and in sceptre's own data, and the genome build can be recorded. Is that the right convention, and should the build be required when positions are present?
4. **Scope.** A gRNA has one target, and so does a vector carrying several gRNAs; combinatorial screens are refused. DNA-targeting is assumed. Are those limits acceptable for now?
5. **The fitting-method default.** Saved fits are shared between R and Python only when the method matches, and the current default differs from the alternatives on a few responses (T6). This is a scientific decision the format does not make.
6. **Is the restriction to additive formulas acceptable?** An analysis with an interaction or a spline cannot be saved.
7. **Twelve smaller questions** that writing the rules down exposed are listed in `RATIONALE.md`, section 11, with how each was settled.

## What you can skip

These can be added or changed later without breaking any file, so they do not
need review now:

- **Optional columns and tables not yet in the format**: gRNA type, spacer, gene annotation, per-cell QC reasons, assignment scores. Section 10 of the rationale lists them.
- **The wording of individual rules.** It will change. The sorting of rules into file rules, behaviour rules and guidance is meant to be final; the sentences are not.
- **Compression settings and Parquet details.** Readers do not depend on them.
- **The prototype code** in `format/`. It exists to show the format works in both languages and will be replaced by code in the packages.
- **Comparisons with IGVF and pertpy** (section 9), unless you know those formats and see an error.

## What has and has not been tested

- **Tested without either package** (`verify/run.sh`): R and Python write identical files from the same content; each reads the other's and gets the source back; tools outside this code (anndata, polars, duckdb, HDF5Array, arrow) read them; two validators agree on a set of valid and deliberately broken files.
- **A notebook** (`notebook.qmd`, rendered with `./render.sh`) walks through the format in R and Python with standard libraries only. It was rendered on macOS only.
- **Tested with the packages, up to commit `7e10af4` only:** every example analysis was rebuilt as a sceptre and a pysceptre object and matched the original, which showed the files held everything an analysis needs. That suite was then removed, because it needed unreviewed builds of both packages on one machine. Changes to the format after that commit have no such check until the packages adopt the format.
- **No open finding** from the independent readers. Two were settled on 18 November: one scanpy function failed on the counts as read, fixed by storing both index arrays in one width (S4); and a dataset reader pointed at the fit cache's folder tripped on a JSON record inside it, fixed by putting the identity of the fits inside each part file (T6).
- **Tested once at scale:** the data file, on a real screen of 207,324 cells, written by both languages to identical files and read back in memory and from disk (`RATIONALE.md`, Appendix A).
- **Not tested:** an analysis larger than 1,000 cells; a second operating system; real data in a second measured modality.
- **Known limit:** everything was written and implemented by one author. An independent implementation from the specification alone has not been attempted.
