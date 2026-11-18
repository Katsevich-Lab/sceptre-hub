# Scale test

Recorded results of one run of the data-file half of the format on a real
large screen: the Gasperini et al. (2019) at-scale data, 207,324 cells, 13,135
genes, 13,077 gRNAs and 699 million nonzero counts, on one cluster node with
16 cores and 90 GB of memory. `RATIONALE.md`, Appendix A, reports them.

| File | What it records |
|---|---|
| `results/results-python.json` | Python writes the data as a data file of the format, validates it, reads it back in memory and gene by gene, and writes an assignment table of real size |
| `results/results-r.json` | R validates the Python-written file, reads it in memory and from disk, writes it again, and reads and rewrites the assignment table |
| `results/results-r-backed.json` | R reads the file through MuData's backed reader, one response at a time |

The results were produced by commit `e490b1c`, whose `scale-test/` holds the
scripts. They predate the index type of S4: those files stored 64-bit column
pointers beside 32-bit cell positions. The scripts were not kept current and
are not in this folder.

Not covered: an analysis at this scale (only the data file and one table), a
second dataset, a second operating system.
