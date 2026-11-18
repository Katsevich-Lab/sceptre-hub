# The verification suite

Every check of the format. It needs neither sceptre nor pysceptre.
`verify/run.sh` runs all of it and prints a summary; `../README.md` says what
each check shows and which libraries it needs.

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
must be concluded about it: the rules broken, whether it conforms, whether a
reader must still read it (S4.4), and whether a reader must refuse it.
`verify/check_conformance.py` and `verify/check_conformance.R` run the validators
and the standard readers over that list; another reader is tested by running it over the same
list and comparing with the same expectations. `example/*/expected/` holds
values a reader should be able to derive from the example analyses.

**What "passing" means.** There is no canonical writer. A writer conforms when
its files pass the validators and read back correctly in the standard
libraries; a reader conforms when it accepts every conforming file and every
file S4.4 names, gets the source back, and refuses what must be refused.

## To do

Two more independent readers in R: **anndataR** and **zellkonverter**. Notes for whoever adds them, neither verified against
current versions: both are believed to read single `.h5ad` files and not
`.h5mu`, so the test would be to export one modality to `.h5ad` and open it;
and zellkonverter reads through Python's anndata, so it largely repeats a check
that exists. They matter once R users open these files outside
sceptre.
