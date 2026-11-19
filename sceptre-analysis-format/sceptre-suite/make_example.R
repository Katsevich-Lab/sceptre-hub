# Writes example/ again from sceptre: the example analyses and, in
# example/<dataset>/expected/, what sceptre itself held for each. Run by hand
# when the format or the examples change, from the folder that holds format/:
#     Rscript sceptre-suite/make_example.R
# The conformance set is built from example/, so build it again afterwards.
source("sceptre-suite/setup.R")
unlink("example", recursive = TRUE)
write_examples(suppressMessages(example_objects()), "example")
print(file_listing("example"))
