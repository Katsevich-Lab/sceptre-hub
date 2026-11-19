# The verification suite

Every check of the sceptre analysis format that needs neither sceptre nor
pysceptre; the checks that need sceptre are in `../sceptre-suite`. It uses the data file format in `../../perturb-seq-data-format`:
its code, its example files and its conformance files. `verify/run.sh` runs
all of it and prints a summary; `verify/run.sh format` runs the part that uses
no writer (the validators and standard readers on the committed files, and the
standard Parquet libraries on each kind of column), and
`verify/run.sh writers` the checks of the two writers in `format/`.
`../overview.html` says what each check shows, and `../README.md` which
libraries it needs.

## Testing another implementation

The suite is built so that a reader or writer can be tested without this
folder importing it. Nothing here needs to know what language it is written in.

**A writer.** `verify/random_inputs.py` writes seeded random content in a
neutral form, `verify/out/random/seed-N/input.json`: identifiers, counts as
(cell, feature, count) triplets, covariates, the gRNA table, and the tables of
an analysis. An implementation reads that file and writes, beside it,
`screen-<name>.h5mu` and `analysis-<name>/`. Then

```sh
python verify/two_writers.py <name> python
```

compares its files with the Python library's, array by array and column by
column, runs the validator on them, and reads them back to compare with the
source. `verify/two_writers.R` is the R library doing exactly this, and shows
what has to be written.

**A reader.** `conformance/cases.json` lists every conformance file with what
must be concluded about it: the rules broken, whether it conforms, the data
file it is checked against, and whether a reader must refuse it.
`verify/check_conformance.py` and `verify/check_conformance.R` run the validators
and the standard readers over that list; another reader is tested by running it over the same
list and comparing with the same expectations. `example/*/expected/` holds
values a reader should be able to derive from the example analyses.

**What "passing" means.** There is no canonical writer. A writer conforms when
its files pass the validators and read back correctly in the standard
libraries; a reader conforms when it accepts every conforming folder, gets the
source back, and refuses what must be refused.
