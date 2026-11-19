# Libraries and code for the sceptre suite: the format code, the reader of a
# data file as a SingleCellExperiment, sceptre, and the code that takes a
# sceptre_object to and from the files. Run from the folder that holds format/.
source("format/setup.R")
source(file.path(Sys.getenv("PERTURB_SEQ_DATA_FORMAT_DIR"), "verify", "as_sce.R"))
suppressPackageStartupMessages(library(sceptre))
SCEPTRE_VERSION <- "0.99.1"
if (as.character(utils::packageVersion("sceptre")) != SCEPTRE_VERSION) {
    stop("The sceptre suite needs sceptre ", SCEPTRE_VERSION, "; found ", utils::packageVersion("sceptre"), ".")
}
source("sceptre-suite/sceptre_objects.R")
source("sceptre-suite/objects.R")
