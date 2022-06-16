#!/usr/bin/env Rscript
library(data.table)
args <- commandArgs(trailingOnly = TRUE)
dt <- fread(input = args[1], sep = " ", header = FALSE, showProgress = FALSE)
s <- summary(dt$V3)
cat(s, "\n")
