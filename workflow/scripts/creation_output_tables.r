#!/usr/bin/env Rscript

# Check if basic packages are installed -----------------------------------

library(tidyverse)
library(data.table)
library(maditr)
library(optparse)

# Script command line options ---------------------------------------------

option_list <- list(
  make_option(c("--clu_or"), type="character", default=NULL,
              help="clusterDB origin size table path", metavar="character"),
  make_option(c("--cat"), type="character", default=NULL,
              help="cluster categories path", metavar="character"),
  make_option(c("--clu_info"), type="character", default=NULL,
              help="clusterDB info table", metavar="character"),
  make_option(c("--comm"), type="character", default=NULL,
              help="cluster communities", metavar="character"),
  make_option(c("--hq_clu"), type="character", default=NULL,
              help="HQ clusters", metavar="character"),
  make_option(c("--k_annot"), type="character", default=NULL,
              help="K annotations", metavar="character"),
  make_option(c("--is_singl"), type="character", default=NULL,
              help="use singleton or not", metavar="character"),
  make_option(c("--s_categ"), type="character", default=NULL,
              help="singelton categories", metavar="character"),
  make_option(c("--res"), type="character", default=NULL,
              help="output results", metavar="character"),
  make_option(c("--threads"), type="numeric", default=1,
              help="number of threads", metavar="numeric")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

dir <- dirname(opt$cat)

if (is.null(opt$clu_or) |
    is.null(opt$cat) | is.null(opt$clu_info) |
    is.null(opt$comm) | is.null(opt$hq_clu) |
    is.null(opt$k_annot) | is.null(opt$threads) |
    is.null(opt$is_singl) | is.null(opt$res)){
  print_help(opt_parser)
  stop("You need to provide the path to the previous validation step results and output files paths\n", call.=FALSE)
}

setDTthreads(opt$threads)
options(datatable.verbose = FALSE)


## GC DB information ################################################################
DB_cl <- fread(opt$clu_or) %>%
  setNames(c("cl_name","db","size")) %>% mutate(cl_name=as.character(cl_name))

# Load cluster DB results
summary_DB <- fread(opt$cat, stringsAsFactors = F, header = F, nThread = 32) %>%
  setNames(c("cl_name","category")) %>% mutate(cl_name=as.character(cl_name))
# Table with gene - cluster - community - category
Communities <- fread(opt$comm,stringsAsFactors = F, header = T) %>% select(cl_name,com) %>%
  mutate(cl_name=as.character(cl_name))
DB_clu_info <- fread(opt$clu_info,stringsAsFactors = F, header = F, nThread = 32) %>%
  setNames(c("cl_name","rep","gene","length","size")) %>% select(-rep,-length,-size) %>%
  mutate(cl_name=as.character(cl_name)) %>% mutate(gene=as.character(gene)) %>%
  dt_inner_join(DB_cl) %>%
  dt_left_join(summary_DB) %>%
  dt_left_join(Communities)

DB_clu_info <- DB_clu_info %>%
  mutate(category=ifelse(size==1,"SINGL",
                         ifelse(size>1 & is.na(category),"DISC",category)))

if(opt$is_singl=="true"){
  s_cat <- fread(opt$s_cat, header=F, sep="\t") %>%
  setNames(c("cl_name","category_s","gene")) %>% select(-cl_name) %>% mutate(gene=as.character(gene))

  DB_clu_info <- DB_clu_info %>% mutate(gene=as.character(gene)) %>%
    dt_left_join(s_cat, by="gene") %>%
    mutate(is.singleton=ifelse(category=="SINGL",TRUE,FALSE),
           category=ifelse(is.na(category_s),category,category_s)) %>% select(-category_s) %>%
           rename(n_genes=size,community=com, origin_db=db) %>%
           select(cl_name, origin_db, n_genes, category, community, is.singleton, gene)
}

write.table(DB_clu_info, paste0(dir,"/DB_cluster_information.tsv"), sep="\t", row.names=F, quote=F)

# Expanded summary info table
## gene - cl_name - category - db - is.HQ - community - pfam
# HQ clusters
HQ_cl <- fread(opt$hq_clu) %>%
  mutate(cl_name=as.character(cl_name), is.HQ=TRUE) %>% select(cl_name,is.HQ)
DB_info_exp <- DB_clu_info %>%
  dt_left_join(HQ_cl)

# Join with Pfam annotations
K_annot <- fread(opt$k_annot, stringsAsFactors = F, header = F) %>%
        select(V1,V2) %>% setNames(c("gene","pfam")) %>% mutate(gene=as.character(gene))

DB_info_exp <- DB_info_exp  %>%
      dt_left_join(K_annot)

DB_info_exp <- DB_info_exp %>%
    distinct() %>% mutate(is.HQ=ifelse(is.na(is.HQ),FALSE,is.HQ)) %>%
    mutate(gene_callers_id=gsub(".*_","",gene))

if(opt$is_singl=="true"){
    DB_info_exp <- DB_info_exp %>%
    select(gene_callers_id,cl_name, community, n_genes,category,is.singleton,is.HQ,pfam)
}else{
   DB_info_exp <- DB_info_exp %>%
   select(gene_callers_id,cl_name,community,n_genes,category,is.HQ,pfam)
}

write.table(DB_info_exp, opt$res, col.names = T, row.names = F, sep="\t",quote = F)
