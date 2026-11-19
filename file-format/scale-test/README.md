# Scale tests

Recorded runs at real size, each on one cluster node (Linux). They are too
large for the verification suite and are run by hand.

## The data file

One run on a real large screen: the Gasperini et al. (2019) at-scale data,
207,324 cells, 13,135 genes, 13,077 gRNAs and 699 million nonzero counts, on
16 cores and 90 GB of memory. `RATIONALE.md`, Appendix A, reports it.

| File | What it records |
|---|---|
| `results/results-python.json` | Python writes the data as a data file of the format, validates it, reads it back in memory and gene by gene, and writes an assignment table of real size |
| `results/results-r.json` | R validates the Python-written file, reads it in memory and from disk, writes it again, and reads and rewrites the assignment table |
| `results/results-r-backed.json` | R reads the file through MuData's backed reader, one response at a time |

The results were produced by commit `e490b1c`, whose `scale-test/` holds the
scripts. They predate the index type of S4: those files stored 64-bit column
pointers beside 32-bit cell positions. The scripts were not kept current and
are not in this folder.

## The tables of an analysis

Made-up tables with the columns, types and sizes of a real analysis: results,
resampled statistics and the fit cache. No statistics are computed and no
package is used, so this measures the layout of the files and nothing else.
One run of each measurement, on 24 cores and 132 GB of memory, with files on a
network file system; pyarrow 26.0.0 in Python, nanoparquet and arrow 26.0.1 in R.

```sh
python analysis/analysis_scale.py run   <work folder> tables.json     # every table, Python and nanoparquet
python analysis/analysis_scale.py run_r <work folder> r-arrow.json    # R's arrow: whole reads and chosen rows
```

`analysis/results/` holds what the two commands wrote. Times are seconds and
memory is the peak of the process, in GB. Before each read the operating
system was asked to drop the file from its cache; on a network file system it
may not have.

### Results tables: no problem at 100 million pairs

One file, zstd, row groups of a million rows, sorted by target then response.

| | 20 million pairs | 100 million pairs |
|---|---|---|
| Size on disk | 0.5 GB | 2.5 GB |
| Write (Python) | 3 s | 15 s |
| Read all: Python / R nanoparquet / R arrow | 0.4 s / 2.5 s / not run | 17 s, 9 GB / 15 s, 8 GB / 11 s, 12 GB |
| Read one target or one pair: Python / R arrow | 0.1 s | 0.1 s / 2 to 3 s |
| Read one response: Python / R arrow | 0.5 s | 6 s, 3 GB / 5 s, 4 GB |

### Resampled statistics: large, and the row group limit is too loose

499 resampled statistics per pair, as random normal values.

| Layout (1 million pairs) | Size on disk |
|---|---|
| The format: a list of doubles, zstd | 3.8 GB |
| The same with snappy | 4.0 GB |
| A list of 32-bit values, zstd | 1.9 GB |
| One column per resample, zstd | 4.7 GB |

Random values do not compress, so the size is that of the numbers: 4 GB per
million pairs, and 38 GB for 10 million pairs.

The size of a row group decides what a selective read costs, because a reader
decodes whole row groups:

| Read (zstd, list of doubles) | Row groups of 1,000,000 rows (the limit of T0.3) | Row groups of 10,000 rows |
|---|---|---|
| One pair, 1 million pairs, Python | 7 s, 11 GB | 0.1 s, 0.3 GB |
| One pair, 1 million pairs, R arrow | 9 s, 11 GB | 1.4 s, 0.4 GB |
| One pair, 10 million pairs, Python | 7 s, 11 GB | 0.1 s, 0.3 GB |
| One response, 10 million pairs, Python | 191 s, 57 GB | 14 s, 43 GB |
| Everything, 1 million pairs: Python / R nanoparquet / R arrow | 10 s / 18 s / 16 s, each 11 GB | 1.3 s / 22 s / 20 s, each 9 to 11 GB |
| Everything, 10 million pairs, Python | 145 s, 81 GB | 35 s, 83 GB |
| Memory to write, 10 million pairs | 14 GB | 0.3 GB |

Reading one response is slow in both because the table is sorted by target
first: a response's rows are spread over every row group.

### The fit cache: time grows with the number of part files

20,000 fits of 20 coefficients each, 4 MB in one file, split over more and
more part files. Reading costs about 9 milliseconds a part, the same whether
or not the files were read just before.

| Parts | Size | Write | Identities only | Read all: Python, part by part | Python, as one dataset | R nanoparquet |
|---|---|---|---|---|---|---|
| 1 | 4 MB | 0.04 s | 0.01 s | 0.03 s | 0.02 s | 0.1 s |
| 100 | 4 MB | 0.9 s | 0.3 s | 0.8 s | 0.2 s | 0.9 s |
| 1,000 | 6 MB | 10 s | 3 s | 9 s | 3 s | 10 s |
| 10,000 | 22 MB | 99 s | 36 s | 90 s | 29 s | 96 s |
| 20,000 | 41 MB | 206 s | 78 s | 176 s | 56 s | 188 s |

## Not covered

A real analysis at this size, written by a package; real resampled statistics,
which may compress differently from random values; a second run of any
measurement; a local disk; a second dataset for the data file.
