# Scale test

Two recorded runs at real size, each on one cluster node (Linux). They are
too large for the verification suite and are run by hand.

## A real large screen

One run on the Gasperini et al. (2019) at-scale data,
207,324 cells, 13,135 genes, 13,077 gRNAs and 699 million nonzero counts, on
16 cores and 90 GB of memory. `RATIONALE.md`, Appendix A, reports it.

| File | What it records |
|---|---|
| `results/results-python.json` | Python writes the data as a data file of the format, validates it, reads it back in memory and gene by gene. It also records an assignment table of real size, which the sceptre analysis format's rationale reports |
| `results/results-r.json` | R validates the Python-written file, reads it in memory and from disk, and writes it again. It also records reading and rewriting the assignment table |
| `results/results-r-backed.json` | R reads the file through MuData's backed reader, one feature at a time |

The results were produced by commit `e490b1c`, whose `scale-test/` holds the
scripts. They predate the index type of D4: those files stored 64-bit column
pointers beside 32-bit cell positions. The scripts were not kept current and
are not in this folder.

## Past the 32-bit index limit

A second run, on a synthetic file made to cross the limit of D4: 1,100,000
cells, 8,000 genes and 2,199,990,946 nonzero counts in `rna`, with a small
`grna` modality of 400 gRNAs. The file is 26.6 GB. It ran on one cluster node
(Linux) with 32 cores and 496 GB of memory. `past-limit/` holds the scripts
and `past-limit/results/` what each step recorded.

| Step | Outcome | Time | Peak memory |
|---|---|---|---|
| Python writes the file: `rna` with both index arrays 64-bit, `grna` with both 32-bit | as D4 requires | 58 s | 53 GB |
| Python validates it | conforms | 48 s | 71 GB |
| Python reads it in memory with mudata, and 20 genes from disk with h5py | equal to what was generated | 28 s; 4 s | 27 GB |
| anndata reads it, in memory and from disk; scanpy computes its count summaries | equal to what was generated | 40 s; 4 s; 18 s | 61 GB |
| R validates the Python-written file | conforms | 242 s | 138 GB |
| R reads 20 genes from disk with HDF5Array, and 3 with MuData's backed reader | equal to what was generated | about 0.1 s a gene | 36 GB |
| R reads it in memory with MuData | stops with an error | 60 s | 27 GB |
| R writes it | not possible | | |

- **The genes checked** include the one whose counts cross position 2³¹ − 1 in
  the stored arrays, and the last ones.
- **R in memory.** R's sparse matrix holds at most 2³¹ − 1 nonzero counts.
  MuData's in-memory reader stops with an error about the column pointers,
  which does not say that the matrix is too large. It returns no wrong counts.
- **R writing.** The prototype's writer takes an in-memory matrix, so it has
  no way to be given this data.
- **Memory.** Both validators read whole arrays at once, which is where their
  peaks come from.

## Not covered

A second real dataset; a file past the limit written by R; a second run of either.
