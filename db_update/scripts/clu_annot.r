#!/usr/bin/env Rscript

# Check if basic packages are installed -----------------------------------

is.installed <- function(pkg) {
  is.element(pkg, installed.packages()[, 1])
}

if (!is.installed("maditr")) {
  cat("We will try to install the package... (this will be only be done once)\n")
  Sys.sleep(5)
  if (!is.installed("maditr")) {
    suppressMessages(install.packages("maditr", repos = "http://cran.us.r-project.org"))
  }
}

library(tidyverse)
library(data.table)
library(maditr)
library(optparse)

# Script command line options ---------------------------------------------

option_list <- list(
  make_option(c("--pfam_annot"),
    type = "character", default = NULL,
    help = "Pfam multi-domain annotations path", metavar = "character"
  ),
  make_option(c("--clusters"),
    type = "character", default = NULL,
    help = "Clustering results path", metavar = "character"
  ),
  make_option(c("--partial"),
    type = "character", default = 1,
    help = "ORF partial info file", metavar = "character"
  ),
  make_option(c("--output_annot"),
    type = "character", default = NULL,
    help = "Output annotated clusters", metavar = "character"
  ),
  make_option(c("--output_noannot"),
    type = "character", default = NULL,
    help = "Output not annotated clusters", metavar = "character"
  )
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$pfam_annot) | is.null(opt$clusters) |
  is.null(opt$partial) | is.null(opt$output_annot) |
  is.null(opt$output_noannot)) {
  print_help(opt_parser)
  stop("You need to provide the path to the annotation and clustering results, and to the partial-info file, and to the output files\n", call. = FALSE)
}

setDTthreads(opt$threads)
options(datatable.verbose = FALSE)
options(readr.num_columns = 0)

# Pfam annotation results
annot <- fread(opt$pfam_annot, stringsAsFactors = F, header = F) %>% setNames(c("orf", "name", "acc", "clan"))
# Clustering results
clu <- fread(opt$clusters, stringsAsFactors = F, header = F) %>% setNames(c("cl_name","rep", "orf"))

clu_annot <- clu %>% dt_left_join(annot, by = "orf")

# ORF completeness info:
# from the gene prediction (Prodigal) with obtain an indicator of if a gene runs off the edge of a sequence or into a gap.
# "0" indicates the gene has a true boundary (a start or a stop),
# "1" indicates the gene is "unfinished" at that edge (i.e. a partial gene).
# For example: "01" means a gene is partial at the right boundary, "11" indicates both edges are incomplete,
# and "00" indicates a complete gene with a start and stop codon.
partial <- fread(opt$partial, header = F, colClasses = c("character", "character"), sep="\t") %>% setNames(c("orf", "partial"))

# Annotated clusters
clu_annot %>%
  group_by(cl_name,rep) %>%
  mutate(annot = ifelse(any(!is.na(name)), "annot", "noannot")) %>%
  dt_filter(annot == "annot") %>%
  dt_select(-annot) %>%
  dt_left_join(partial, by = "orf") %>%
  select(cl_name,rep,orf,name,acc,clan,partial) %>%
  write_tsv(path = opt$output_annot, col_names = F)

# Not annotated clusters
clu_annot %>%
  group_by(cl_name,rep) %>%
  mutate(annot = ifelse(any(!is.na(name)), "annot", "noannot")) %>%
  dt_filter(annot == "noannot") %>%
  dt_select(-annot) %>%
  dt_left_join(partial, by = "orf") %>%
  select(cl_name,rep,orf,name,acc,clan,partial) %>%
  write_tsv(path = opt$output_noannot, col_names = F)
