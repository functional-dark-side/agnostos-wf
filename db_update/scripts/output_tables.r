#!/usr/bin/env Rscript

# Check if basic packages are installed -----------------------------------

library(tidyverse)
library(data.table)
library(maditr)
library(optparse)

# Script command line options ---------------------------------------------

option_list <- list(
  make_option(c("-o", "--clu_or"), type="character", default=NULL,
              help="clusterDB origin size table path", metavar="character"),
  make_option(c("-g", "--contig"), type="character", default=NULL,
              help="contig genes table path", metavar="character"),
  make_option(c("-c", "--cat"), type="character", default=NULL,
              help="cluster categories path", metavar="character"),
  make_option(c("-i", "--info"), type="character", default=1,
              help="clusterDB info table", metavar="character"),
  make_option(c("-n", "--name"), type="character", default=1,
              help="new dataset name", metavar="character"),
  make_option(c("-m", "--comm"), type="character", default=1,
              help="cluster communities", metavar="character"),
  make_option(c("-h", "--hq_cl"), type="character", default=1,
              help="HQ clusters", metavar="character"),
  make_option(c("-k", "--k_annot"), type="character", default=1,
              help="K annotations", metavar="character"),
  make_option(c("-w", "--kwp_annot"), type="character", default=1,
              help="KWP annotations", metavar="character"),
  make_option(c("-u", "--gu_annot"), type="character", default=1,
              help="GU annotations", metavar="character")
  make_option(c("-t", "--threads"), type="character", default=1,
              help="number of threads", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

dir <- dirname(opt$contig)

if (is.null(opt$clu_or) | is.null(opt$contig) |
    is.null(opt$cat) | is.null(opt$info) |
    is.null(opt$name) | is.null(opt$comm) |
    is.null(opt$hq) | is.null(opt$annot) |
    is.null(opt$kwp) | is.null(opt$gu) | is.null(opt$threads)){
  print_help(opt_parser)
  stop("You need to provide the path to the previous validation step results and output files paths\n", call.=FALSE)
}

setDTthreads(opt$threads)
options(datatable.verbose = FALSE)


## AgnostosDB update (integration) ################################################################

DB_cl <- fread(opt$clu_or) %>%
  setNames(c("cl_name","db","size")) %>% filter(grepl(opt$name,db)) %>% mutate(cl_name=as.character(cl_name))

DB_cl %>% filter(db==opt$name) %>% nrow()
DB_cl %>% filter(db==opt$name) %>% .$size %>% sum()

# Load gene - contig table
gene_info <- fread(opt$contig,header=F,sep="\t") %>% setNames(c("contig","gene")) %>%
  group_by(contig) %>% add_count() %>%
  rename(gene_x_contig=n) %>% ungroup()

# Load DB integration results
summary_DB <- fread(opt$cat, stringsAsFactors = F, header = F, nThread = 32) %>%
  setNames(c("cl_name","category")) %>% mutate(cl_name=as.character(cl_name))
DB_clu_info <- fread(opt$info,stringsAsFactors = F, header = F, nThread = 32) %>%
  setNames(c("cl_name","rep","gene","length","size")) %>% select(-rep,-length,-size) %>%
  mutate(cl_name=as.character(cl_name)) %>%
  dt_inner_join(DB_cl) %>%
  dt_left_join(summary_DB)

DB_info <- gene_info %>% dt_left_join(DB_clu_info) %>%
  mutate(category=ifelse(size==1,"SINGL",
                         ifelse(size>1 & is.na(category),"DISC",category)))

# Table with gene - cluster - community - category
Communities <- fread(opt$comm,stringsAsFactors = F, header = T)

DB_GC_GCCs <- DB_info %>% select(gene,cl_name,category) %>%
  left_join(Communities %>% mutate(cl_name=as.character(cl_name)) %>% select(-category))
write.table(DB_GC_GCCs,paste0(dir,"/DB_genes_clusters_communities.tsv"), col.names = T, row.names = F, quote = F, sep = "\t")

# Minimal summary info table
DB_info_red <- DB_info %>% select(contig,gene,cl_name,category,db)
write.table(DB_info_red,paste0(dir,"/DB_genes_summary_info_red.tsv"), col.names = T, row.names = F, quote = F, sep = "\t")

# Expanded summary info table
## contig - gene - cl_name - category - db - is.HQ - is.LS - project - Niche-breadth

# HQ clusters
HQ_cl <- fread(opt$hq_clu) %>%
  mutate(cl_name=as.character(cl_name), is.HQ=TRUE) %>% select(cl_name,is.HQ)
DB_info_exp <- DB_info %>%
  dt_left_join(HQ_cl)

# join with  lineage-specificity)
ls_clu <- fread(paste0(dir,"/agnostosDB_dbf02445-20200519_phylogenetic/mg_gtdb_lineage_specific_clusters.tsv.gz"),nThread = opt$threads)

DB_info_exp <- DB_info_exp %>% dt_left_join(ls_clu %>% mutate(cl_name=as.character(cl_name)) %>%
                                              mutate(is.LS=TRUE) %>% select(cl_name,is.LS,lowest_rank,lowest_level))

#join with environmental/metagenomic projecs
mg <- fread(paste0(dir,"/agnostosDB_dbf02445-20200519_environmental/mg_cluster_sample_norfs_coverage.tsv.gz"),nThread = opt$threads ) #environmental

mg_project <- mg %>% mutate(project=case_when(grepl('TARA',sample) ~ "TARA",
                                              grepl('SRS',sample) ~ "HMP",
                                              grepl('MP',sample) ~ "Malaspina",
                                              grepl('OSD',sample) ~ "OSD",
                                              TRUE ~ "GOS")) %>%
  group_by(cl_name,project) %>% summarise(norfs=sum(norfs),coverage=sum(coverage))

DB_info_exp <- DB_info_exp %>% dt_left_join(mg_project %>% ungroup() %>% mutate(cl_name=as.character(cl_name)))

# Join with Nieche breadth info
load(paste0(dir,"/agnostosDB_dbf02445-20200519_environmental/niche_breadth/gCl_nb_all_mv.Rda"))
DB_info_exp <- DB_info_exp %>%
  dt_left_join(gCl_nb_all_mv %>%
                 rename(niche_breadth_sign=sign_mv,cl_name=gCl_name) %>%
                 mutate(cl_name=as.character(cl_name)) %>% select(-categ))
# Join with Pfam annotations
K_annot <- fread(opt$k_annot, stringsAsFactors = F, header = F) %>%
  select(V1,V4) %>% setNames(c("cl_name","pfam"))

DB_info_exp <- DB_info_exp  %>%
  dt_left_join(K_annot %>% mutate(cl_name=as.character(cl_name)) %>% dt_filter(pfam!="mono", pfam!="multi"))


DB_info_exp <- DB_info_exp %>%
  distinct() %>% mutate(is.HQ=ifelse(is.na(is.HQ),FALSE,is.HQ), is.LS=ifelse(is.na(is.LS),FALSE,is.LS)) %>%
  rename(cl_size=size,mg_project=project)

write.table(DB_info_exp,paste0(dir,"/DB_genes_summary_info_exp.tsv"), col.names = T, row.names = F, sep="\t",quote = F)

# Create table with all annotations
kwp_annot <- fread(opt$kwp_annot, stringsAsFactors = F, header = F) %>%
  select(V1,V4) %>% setNames(c("cl_name","annot"))
kwp_annot <- kwp_annot %>% dt_filter(!grepl("Unchar|Hypo|hypo|unchar",annot))

gu_annot <- fread(opt$gu_annot , stringsAsFactors = F, header = F) %>%
  select(V1,V4) %>% setNames(c("cl_name","annot"))

other_annot <- bind_rows(kwp_annot,gu_annot) %>% distinct()
DB_annot <- DB_info_exp %>% select(cl_name,category,pfam) %>% distinct() %>%
  dt_left_join(other_annot %>% mutate(cl_name=as.character(cl_name)) %>% rename(other_annot=annot)) %>%
  distinct()
write.table(DB_annot,paste0(dir,"/DB_cluster_annotations.tsv"), col.names = T, row.names = F, sep="\t",quote = F)

# Two tables with presence in GTDB genomes and MG samples
gtdb <- fread(paste0(dir,"agnostosDB_dbf02445-20200519_phylogenetic/gtdb_cluster_genome_norfs.tsv.gz"), nThread = opt$threads) #phylogenetic

DB_genomes <- DB_info %>% select(cl_name,category) %>% distinct() %>%
  inner_join(gtdb %>% mutate(cl_name=as.character(cl_name)), by="cl_name")
write.table(DB_genomes,paste0(dir,"/DB_clusters_in_gtdb_genomes.tsv"), col.names = T, row.names = F, sep="\t",quote = F)

DB_metagenomes <- DB_info %>% select(cl_name,category) %>% distinct() %>%
  inner_join(mg %>% ungroup() %>% mutate(cl_name=as.character(cl_name))) %>% distinct() %>%
  mutate(project=case_when(grepl('TARA',sample) ~ "TARA",
                           grepl('SRS',sample) ~ "HMP",
                           grepl('MP',sample) ~ "Malaspina",
                           grepl('OSD',sample) ~ "OSD",
                           TRUE ~ "GOS"))
write.table(DB_metagenomes,paste0(dir,"/DB_clusters_in_metagenomes.tsv"), col.names = T, row.names = F, sep="\t",quote = F)
