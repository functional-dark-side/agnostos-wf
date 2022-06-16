#!/usr/bin/env Rscript

# Check if basic packages are installed -----------------------------------

is.installed <- function(pkg){
  is.element(pkg, installed.packages()[,1])
}

if (!is.installed("valr")){
  cat("We will try to install the package... (this will be only be done once)\n")
  Sys.sleep(5)
  if (!is.installed("valr")){
    suppressMessages(install.packages("BiocManager", repos = "http://cran.us.r-project.org"))
    suppressMessages(BiocManager::install("rtracklayer"))
    suppressMessages(install.packages("valr", repos = "http://cran.us.r-project.org"))
  }
}

library(tidyverse)
library(data.table)
library(valr)
library(parallel)
library(optparse)

# Script command line options ---------------------------------------------

option_list <- list(
  make_option(c("-o", "--orfs"), type="character", default=NULL,
              help="Predicted ORFs path", metavar="character"),
  make_option(c("-s", "--shadows"), type="character", default=NULL,
              help="Shadow ORFs results path", metavar="character"),
  make_option(c("-t", "--threads"), type="numeric", default=NULL,
              help="Number of threads", metavar="numeric")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$orfs) | is.null(opt$shadows)){
  print_help(opt_parser)
  stop("You need to provide the path to the ORFs and to the output folder\n", call.=FALSE)
}

setDTthreads(opt$threads)
options(datatable.verbose = FALSE)

# Read table with orfs names
orfs <- fread(opt$orfs,stringsAsFactors = F, header = F, sep="\t") %>% select(V1)

orfs1 <- orfs %>% mutate(info=str_match(V1,"\\+_(.*?)_orf")[,1]) %>%
  separate(col=info,into=c("strand","start","end","orf"), sep="_") %>%
  select(-orf) %>% drop_na() %>%
  mutate(contig=gsub("_\\+.*","",V1))
orfs2 <- orfs %>% mutate(info=str_match(V1,"\\-_(.*?)_orf")[,1]) %>%
  separate(col=info,into=c("strand","start","end","orf"), sep="_") %>%
  select(-orf) %>% drop_na() %>%
  mutate(contig=gsub("_\\-.*","",V1))

orfs <- rbind(orfs1,orfs2) %>% setNames(c("orf","strand","start","end","contig")) %>% #select(contig) %>% distinct()
  group_by(contig) %>% add_count(contig) %>% filter(n>1) %>% select(-n)

orfs <- orfs %>% ungroup() %>% data.table()
orfs$start <- as.integer(orfs$start)
orfs$end <- as.integer(orfs$end)
contigs <- orfs[, .(count=.N), by = contig]
contigs <- contigs[count > 1]

overlaps_dt_f <- function(X,tbl){
  tbl <- tbl[contig == X]

  setkey(tbl, contig, start, end)
  over <- foverlaps(tbl, tbl, nomatch = 0)
  over <- over[orf != i.orf]
  over <- over[!duplicated(data.table(pmin(orf,i.orf),pmax(orf,i.orf)))]
  over[, overlap := abs(over[, ifelse(start > i.start, start, i.start)] - over[, ifelse(end < i.end, end, i.end)])]
  over <- over[overlap >= 50]

  same_str <- over %>% filter(strand == i.strand) %>% filter(overlap>=60)
  opp_str1 <- over %>% filter(strand != i.strand) %>% filter(overlap>=50) %>%
    filter(end <= i.end & i.start >= start || i.end <= end & start >= i.start) # both 3' included
  opp_str2 <- over %>% filter(strand != i.strand) %>% filter(overlap>=120) %>%
    filter(start >= i.start || i.end <= end) # at least one 5' included
  overl_orfs <- rbind(same_str, opp_str1, opp_str2)
  if (dim(overl_orfs)[1] == 0) {
    return(NULL)
  }else{
    overl_orfs
  }
}

shadows <- mclapply(contigs$contig,overlaps_dt_f,tbl=orfs, mc.cores = opt$threads - 2)
shadows <- Filter(Negate(is.null),shadows)
shadows <- plyr::ldply(shadows,data.frame)
write.table(shadows,opt$shadows,col.names = F, row.names = F, sep = "\t", quote = F, append = T)
