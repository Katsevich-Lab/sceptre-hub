# Scale test

A recorded run at real size, on one cluster node (Linux). It is too large for
the verification suite and is run by hand.

One run on a real large screen: the Gasperini et al. (2019) at-scale data,
207,324 cells, 13,135 genes, 13,077 gRNAs and 699 million nonzero counts, on
16 cores and 90 GB of memory. `RATIONALE.md`, Appendix A, reports it.

| File | What it records |
|---|---|
| `results/results-python.json` | Python writes the data as a data file of the format, validates it, reads it back in memory and gene by gene. It also records an assignment table of real size, which the sceptre analysis format's rationale reports |
| `results/results-r.json` | R validates the Python-written file, reads it in memory and from disk, and writes it again. It also records reading and rewriting the assignment table |
| `results/results-r-backed.json` | R reads the file through MuData's backed reader, one feature at a time |

The results were produced by commit `e490b1c`, whose `scale-test/` holds the
scripts. They predate the index type of S4: those files stored 64-bit column
pointers beside 32-bit cell positions. The scripts were not kept current and
are not in this folder.

## Not covered

A second dataset; a dataset past 2³¹ − 1 nonzero counts; a second run.
