# Source the contents of the R/ subirectory when this repo is cloned as a
# submodule _at the root_ of another repo
sapply(list.files("./rep2si-common/R", full.names = TRUE), source)
