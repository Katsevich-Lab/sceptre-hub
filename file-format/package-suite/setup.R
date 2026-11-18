# Libraries for the package suite: the format code, plus the sceptre3 build of
# sceptre from the tour library (see README.md) and the code that builds a
# sceptre_object from the files. Run from the folder that holds format/.
source("format/setup.R")
.libPaths(c(Sys.getenv("SCEPTRE_TOUR_RLIB", path.expand("~/code/research/sceptre-tour-rlib")), .libPaths()))
suppressPackageStartupMessages(library(sceptre))
source(file.path(proposal_dir, "package-suite", "sceptre_objects.R"))
