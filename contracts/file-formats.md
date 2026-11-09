# File format registry

Every file that SCEPTRE writes for later reading (`screen.h5mu` and, in
sceptre3, the analysis directory around it) carries a stamp in its metadata:
a JSON object with an integer `format`. A reader compares that number with the
highest format it knows and refuses a newer one with a message to upgrade.
This file is the only registry of what each number means. A number is never
reused, and a new layout gets a new number here before any code writes it.

| Format | Layout | First written by | Read by |
| --- | --- | --- | --- |
| 1 | `screen.h5mu` alone: two MuData modalities, `response` and `grna`, counts as sparse `X` (cells × features), cell covariates in `/mod/response/obs`, feature names as `var` names, metadata `sceptre_format` and the optional entries `covariate_text_columns` and `cell_names_absent`. Everything else (gRNA targets, analysis parameters, assignments, results) stays in `sceptre_object.rds` beside it. Logical covariates are stored as int8 (a known limitation; Python reads them as integers). | sceptre 0.99.2 (Bioconductor 3.24) | sceptre ≥ 0.99.2; sceptre3 and pysceptre must keep reading it |
| 2 | Reserved for sceptre3's analysis directory (`screen.h5mu` plus `state.h5ad` and Parquet tables), once it is settled. Not yet written by any release. | — | — |

## History

The sceptre3 candidate used the internal numbers 1 to 3 during October 2026
for three successive layouts of its analysis directory (no stamp; state as JSON
per key; gRNA assignments as Parquet edge tables). None was released, and the
candidate's reader accepts only the last of them. Those numbers are retired:
the 3.24 release of the R package is the first public format and is numbered
1 (ruling R13, 2026-11-09). The candidate's next layout takes 2 and must read
1. Directories written by the candidate with the old numbers are regenerated.

## Rules

- A format number changes only when a reader of the previous number would
  misread the file. Additive metadata entries do not change it.
- A release that introduces a number records it here in the same change, with
  the layout in one paragraph and a pointer to the test that reads a frozen
  file of that format.
- Readers accept every number up to their own and refuse higher ones with a
  message naming the release that can read the file.
