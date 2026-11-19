# The verification suite

Every check of the Perturb-seq data file format. It needs no analysis tool.
`verify/run.sh` runs all of it and prints a summary. [The overview](../overview.qmd) says
what each check shows, and `../README.md` which libraries it needs.

## Testing another implementation

**A reader or a validator.** `conformance/cases.json` lists every conformance
file with what must be concluded about it: the rules broken, whether it
conforms, whether a reader must still read it (D4.4), and whether a reader must
refuse it. `verify/check_conformance.py` and `verify/check_conformance.R` run
the validators and the standard readers over that list; another implementation
is tested by running it over the same list and comparing with the same
expectations.

**A writer.** There is no canonical writer. A writer conforms when its files
pass the validators and read back correctly in the standard libraries.
`verify/writers.py` and `verify/writers.R` check what the two writers in
`format/` do beyond that: a copy keeps the dataset's identifier and every
modality, and past the index limit both index arrays are 64-bit. The two
writers are compared with each other on the same random content in the
sceptre analysis format's suite, which writes an analysis as well.

## To do

Two more independent readers in R: **anndataR** and **zellkonverter**. Notes for whoever adds them, neither verified against
current versions: both are believed to read single `.h5ad` files and not
`.h5mu`, so the test would be to export one modality to `.h5ad` and open it;
and zellkonverter reads through Python's anndata, so it largely repeats a check
that exists.
