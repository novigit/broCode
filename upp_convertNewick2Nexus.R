#!/usr/local/bin/Rscript
##!/pdc/vol/R/3.2.0/bin/Rscript

## Author: Jimmy Saw
## Date modified: 2015-06-08
## Description: This script converts newick to nexus format

library(ape)
args <- commandArgs(TRUE)
tree <- read.tree(args[1])
outfile <- paste0(make.names(args[1]), ".nex")
write.nexus(tree, file=outfile)