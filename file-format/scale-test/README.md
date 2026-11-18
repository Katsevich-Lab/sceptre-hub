# Scale test

One run of the data-file half of the format on a real large screen: the
Gasperini et al. (2019) at-scale data, 207,324 cells, 13,135 genes, 13,077
gRNAs and 699 million nonzero counts. It was run on 18 November 2026 on one
node of a cluster with 16 cores and 90 GB of memory. `RATIONALE.md`, Appendix
A, reports the results.

| File | What it does |
|---|---|
| `scale_py.py convert` | Reads the data from the October layout, writes it as a data file of the format, validates it, reads it back in memory and gene by gene, and writes an assignment table of real size |
| `scale_r.R` | Validates the Python-written file, reads it in memory and from disk, writes it again from R, and reads and rewrites the assignment table |
| `scale_py.py compare` | Checks that the R-written file equals the Python-written one |
| `scale_r_backed.R` | Reads the file the way the `sceptre3` branch does: MuData's backed reader, one response at a time |
| `results/` | What the three scripts recorded |

Run from the folder that holds `format/` and `spec/`, with
`SCEPTRE_SCALE_SOURCE` set to a folder holding `screen.h5mu` and
`grna_targets.parquet` in the October layout.

Not covered: an analysis at this scale (only the data file and one table), a
second dataset, a second operating system.
