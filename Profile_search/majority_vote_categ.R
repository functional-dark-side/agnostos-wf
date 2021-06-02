#!/usr/bin/env Rscript

is.installed <- function(pkg){
  is.element(pkg, installed.packages()[,1])
}

if (!is.installed("tidyverse") || !is.installed("data.table") || !is.installed("maditr")){
  cat("We will try to install the packages... (this will be only be done once)\n")
  Sys.sleep(5)
  if (!is.installed("tidyverse")){
    suppressMessages(install.packages("tidyverse", repos = "https://cloud.r-project.org/"))
  }
  if (!is.installed("data.table")){
    suppressMessages(install.packages("data.table", repos = "https://cloud.r-project.org/"))
  }
  if (!is.installed("maditr")){
    suppressMessages(install.packages("maditr", repos = "https://cloud.r-project.org/"))
  }
}

library(tidyverse)
library(data.table)
library(parallel)
library(maditr)

# Functions
majority_vote <- function (x, seed = 12345) {
  set.seed(seed)
  whichMax <- function(x) {
    m <- seq_along(x)[x == max(x, na.rm = TRUE)]
    if (length(m) > 1)
      sample(m, size = 1)
    else m
  }
  x <- as.vector(x)
  tab <- table(x)
  m <- whichMax(tab)
  out <- list(table = tab, ind = m, majority = names(tab)[m])
  return(out)
}

apply_majority <- function(X){
  DT.1 <- X[,majority:=majority_vote(category)$majority, by="gene"]
  df <- DT.1 %>% as_tibble() %>% distinct()
}

get_majority <- function(X){

  list_genes <- X %>%                                        # Split into groups by gene-caller-id
    split(.$gene)

  maj_l <- mclapply(list_genes,apply_majority, mc.cores=8) # run majority_vote function
  maj_df <- plyr::ldply(maj_l, data.frame) %>%  # bind list rowwise and get distint votes for gene category
    select(gene,majority) %>%
    distinct() %>%
    as_tibble()
}

# Arguments from terminal: tsv-table with seaerch results
args = commandArgs(trailingOnly=TRUE)

res=basename(args[1])
res=gsub(".tsv.gz","",res)
dir=dirname(args[1])

# Find consensus category
# 1. read e-value filtered results
efilters <- fread(args[1], stringsAsFactors = F, header = F) %>% dt_select(V1,V2,V11,V13,V14,V15) %>%
  setNames(c("cl_name","gene","evalue","qcov","tcov","category"))
# get major categories
votes <- get_majority(efilters)
# Summary of consensus categories
df_summ <- efilters %>% left_join(votes) %>% group_by(gene) %>%
  arrange(evalue, desc(qcov),desc(tcov)) %>% mutate(is_best=ifelse(row_number()==1,TRUE,FALSE)) %>%
  select(gene,majority,cl_name,category,is_best,evalue)
#write.table(df_summ,paste(getwd(),"/",res,"_summary.tsv",sep=""), col.names = T, row.names = F, quote = F, sep = "\t")
df_best <- df_summ %>% filter(is_best==TRUE) %>% select(-is_best,-category) %>% rename(category=majority) %>%
  mutate(gene_callers_id=gsub(".*_","",gene) %>% select(gene,cl_name,category,evalue)
write.table(df_best,paste(dir,"/",res,"_best-hits.tsv",sep=""), col.names = T, row.names = F, quote = F, sep = "\t")

# Join with info about genomes/MAGs
# Must be a table with the corrispondence of the genes to the contigs (and in case also to the genomes/MAGs)
#The format should be gene - genome (or sample_ID) etc..
#info <- fread(args[2], stringsAsFactors = F, header = T) %>%
#setNames(c("gene","sample"))

#info <- info %>% group_by(sample) %>% add_count() %>%
#  rename(total_ngenes=n) %>% ungroup()

#res_info_class <- df_best %>% left_join(info,by="gene") %>%
#  mutate(class=case_when(grepl('U',category) ~ "Unknown",
#                         TRUE ~ "Known")) %>%
#  group_by(sample, total_ngenes, class) %>% count()
#write.table(res_info_class,paste0(dir,"/",res,"_summary-classes.tsv"), col.names = T, row.names = F, quote = F, sep = "\t")

#res_info_categ <- df_best %>% left_join(info,by="gene") %>%
#  group_by(sample, total_ngenes, category) %>% count()
#write.table(res_info_categ,paste0(dir,"/",res,"_summary-categ.tsv"), col.names = T, row.names = F, quote = F, sep = "\t")
