#!/usr/bin/env Rscript

# Check if basic packages are installed -----------------------------------

is.installed <- function(pkg) {
  is.element(pkg, installed.packages()[, 1])
}

if (!is.installed("entropy") || !is.installed("ggridges")) {
  cat("We will try to install the packages... (this will be only be done once)\n")
  Sys.sleep(5)
  if (!is.installed("entropy")) {
    suppressMessages(install.packages("RSQlite", repos = "http://cran.us.r-project.org"))
  }
  if (!is.installed("ggridges")) {
    suppressMessages(install.packages("dbplyr", repos = "http://cran.us.r-project.org"))
  }
}

library(data.table)
library(tidyverse)
library(maditr)
library(entropy)
library(ggridges)
library(optparse)
library(parallel)

# Script command line options ---------------------------------------------

option_list <- list(
  make_option(c("-r", "--ref_clu"),
    type = "character", default = NULL,
    help = "Refined set of clusters", metavar = "character"
  ),
  make_option(c("-a", "--clu_categ"),
    type = "character", default = NULL,
    help = "Cluster ids-categories", metavar = "character"
  ),
  make_option(c("-t", "--mmseqs_tax"),
    type = "character", default = NULL,
    help = "Cluster taxonomy results", metavar = "character"
  ),
  make_option(c("-k", "--kaiju_tax"),
    type = "character", default = NULL,
    help = "Cluster taxonomy results", metavar = "character"
  ),
  make_option(c("-d", "--clu_dark"),
    type = "character", default = NULL,
    help = "Cluster darkness results", metavar = "character"
  ),
  make_option(c("-i", "--dpd_info"),
    type = "character", default = NULL,
    help = "DPD info table", metavar = "character"
  ),
  make_option(c("-p", "--compl"),
    type = "character", default = NULL,
    help = "output cluster completeness", metavar = "character"
  ),
  make_option(c("-q", "--hq_clu"),
    type = "character", default = NULL,
    help = "output HQ clusters", metavar = "character"
  ),
  make_option(c("-s", "--summ_stats"),
    type = "character", default = NULL,
    help = "output summary stats", metavar = "character"
  ),
  make_option(c("-o", "--output"),
    type = "character", default = NULL,
    help = "Stats output folder", metavar = "character"
  )
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$ref_clu) | is.null(opt$clu_categ) |
  is.null(opt$mmseqs_tax) | is.null(opt$kaiju_tax) | is.null(opt$clu_dark) |
  is.null(opt$dpd_info) | is.null(opt$compl) |
  is.null(opt$hq_clu) | is.null(opt$summ_stats) |
  is.null(opt$output)) {
  print_help(opt_parser)
  stop("You need to provide the path to the refiend clusters, taxonomy and darkness results and output paths\n", call. = FALSE)
}

options(datatable.verbose = FALSE)

# Refined cluster with size orf length and partial info
ref_clu <- fread(opt$ref_clu,
  stringsAsFactors = F, header = F,
  colClasses = c(V1 = "character", V6 = "character")
) %>%
  setNames(c("cl_name", "rep", "orf", "size", "length", "partial"))

# Join with category info
category <- fread(opt$clu_categ,
  stringsAsFactors = F, header = F,
  colClasses = c(V1 = "character")
) %>%
  setNames(c("cl_name", "category"))

ref_clu <- ref_clu %>% dt_inner_join(category)

# Start caluclating stats per cluster catgory (or first xcluster and then xcategory????)
# cluster ORF length stats
clu_stats <- ref_clu %>%
  group_by(cl_name) %>%
  mutate(size=max(size)) %>%
  mutate(
    min_len = min(length), mean_len = mean(length),
    median_len = median(length), max_len = max(length), sd_len = sd(length)
  )

# Cluster completeness stats (then summarize by cluster category)
clu_compl <- clu_stats %>%
  group_by(cl_name, partial, size) %>%
  count() %>%
  mutate(p = n / size) %>%
  ungroup() %>%
  dt_select(cl_name, partial, p) %>%
  spread(partial, p, fill = 0)

if("01" %in% colnames(clu_compl)){
  clu_compl <- clu_compl %>%
    rename(p00 = `00`, p10 = `10`, p01 = `01`, p11 = `11`)
} else {
  clu_compl <- clu_compl %>%
    mutate(p00 =`00`, p01=0, p10=0, p11= `11`) %>% select(-`00`,-`11`)
}

write_tsv(clu_compl, path=opt$compl, col_names = T)

clu_stats <- clu_stats %>% dt_left_join(clu_compl) %>%
  mutate(rep_compl=ifelse(rep==orf & partial=="00", TRUE,FALSE)) %>%
  select(-orf,-length,-partial) %>% distinct()

# HQ cluster set
# Completeness: Percentage of complete ORFs ("00") (tranformed with sinus function to be between 0-0.5) +
# 0.5 if repres is complete or 0 if it's not

pcompl <- clu_stats %>%
  select(cl_name, p00, rep_compl) %>%
  mutate(rep_compl_v = ifelse(rep_compl == T, 0.5, 0), rep.is.compl = ifelse(rep_compl == T, 1, 0))

# BrokenStick to p00 distribution (where representative is complete)
brStick <- function(X) {
  x <- X[[2]]
  m <- 0
  out <- matrix(NA, ncol = 2, nrow = length(x))
  colnames(out) <- c("Observed", "BSM")

  # colnames(out) <- c("% of Variability", "B-Stick Threshold")
  for (i in 1:length(x)) {
    for (k in i:length(x)) {
      m <- m + ((1 / length(x)) * (1 / k))
    }
    out[i, ] <- c((x[i] / sum(x)), m)
    m <- 0
  }
  out <- as_tibble(out) %>% mutate(thresh = X[[1]])
  out_w <- out %>%
    gather(class, value, -thresh) %>%
    mutate(
      thresh = as.character(thresh),
      class = fct_rev(class)
    )
  plot <- ggplot(out_w, aes(thresh, value, fill = class)) +
    geom_col(position = "dodge", color = "black", alpha = 0.7) +
    geom_line(aes(group = class, color = class), position = position_dodge(width = 0.9)) +
    geom_point(position = position_dodge(width = 0.9), colour = "black", shape = 21) +
    theme_light() +
    theme(
      legend.position = "top",
      legend.title = element_blank()
    ) +
    scale_y_continuous(labels = scales::percent) +
    xlab("Filtering threshold") +
    ylab("Variability")

  h_bsm <- out %>%
    filter(Observed > BSM) %>%
    .$thresh

  return(list(bsm_table = out, plot = plot, thresh_by_bsm = h_bsm))
}
br_compl_cl <- plyr::ldply(seq(0, 1, 0.05), function(x) {
  data.frame(threshold = x, clusters = dim(pcompl %>% filter(p00 >= x))[1])
})
br_compl <- brStick(br_compl_cl)
br_compl$plot
lag <- br_compl$thresh_by_bsm %>%
  enframe() %>%
  mutate(lag = round(value - dplyr::lag(value), 2)) %>%
  filter(lag > .01) %>%
  top_n(1) %>%
  .$name %>%
  head(1)
if (length(lag) != 0) {
  rej_threshold <- br_compl$thresh_by_bsm[lag - 1]
} else {
  rej_threshold <- br_compl$thresh_by_bsm[length(br_compl$thresh_by_bsm)]
}

# We define as HQ clusters the clusters with a percentage of complete ORFs higher than the identified threshold
# and a complete representative
HQ_clusters <- clu_stats %>%
  filter(p00 >= rej_threshold, rep_compl == T) %>%
  select(cl_name, category)
write_tsv(HQ_clusters, path = opt$hq_clu, col_names = T)

# Category stats
# category ORF length stats
cat_stats <- ref_clu %>%
  group_by(category) %>%
  mutate(min_len = min(length), mean_len = mean(length), median_len = median(length), max_len = max(length), sd_len = sd(length))
# category size stats
cat_stats <- cat_stats %>%
  group_by(category) %>%
  mutate(min_size = min(size), mean_size = mean(size), median_size = median(size), max_size = max(size), sd_size = sd(size))
# category completeness stats
cat_compl <- cat_stats %>%
  group_by(category, partial) %>%
  add_count() %>%
  mutate(p = n / sum(size)) %>%
  ungroup() %>%
  dt_select(category, partial, p) %>%
  distinct() %>%
  spread(partial, p, fill = 0)

  if("01" %in% colnames(cat_compl)){
    cat_compl <- cat_compl %>%
      rename(p00 = `00`, p10 = `10`, p01 = `01`, p11 = `11`)
  } else {
    cat_compl <- cat_compl %>%
      mutate(p00 =`00`, p01=0, p10=0, p11= `11`) %>% select(-`00`,-`11`)
  }

# Cluster taxonomy stats
cl_mmseqs <- fread(opt$mmseqs_tax,
  stringsAsFactors = F, header = F, fill = T, sep="\t") %>%
  setNames(c("orf", "cl_name", "category", "taxid", "rank", "classific", "lineage"))
# Filter out the unclassified ones and split the taxonomic levels in the lineage
cl_mmseqs <- cl_mmseqs %>%
  filter(classific != "unclassified") %>%
  mutate(lineage = gsub("-_cellular__organisms;", "", lineage)) %>%
  mutate(
    domain = str_match(lineage, "^d_(.*?);.*")[, 2],
    phylum = str_match(lineage, ";p_(.*?);.*")[, 2],
    class = str_match(lineage, ";c_(.*?);.*")[, 2],
    order = str_match(lineage, ";o_(.*?);.*")[, 2],
    family = str_match(lineage, ";f_(.*?);.*")[, 2],
    genus = str_match(lineage, ";g_(.*?);.*")[, 2],
    species = str_match(lineage, ";s_(.*?)$")[, 2]
  )  %>%
  mutate(domain=case_when(classific=="Bacteria" ~ "Bacteria",
                          classific=="Archaea" ~ "Archaea",
                          classific=="Eukaryota" ~ "Eukaryota",
                          classific=="Viruses" ~ "Viruses",
                          TRUE ~ domain)) %>%
  unite(mmseqs_tax, domain:species,sep=";") %>% select(orf,cl_name,category,mmseqs_tax)

cl_kaiju <- fread(opt$kaiju_tax,
  header = F, fill=T, sep="\t") %>%
  setNames(c("orf","domain","phylum","class","order","family","genus","species")) %>%
  unite(kaiju_tax, domain:species,sep=";")

cl_kaiju <- cl_kaiju %>% dt_left_join(ref_clu %>% select(cl_name,category,orf))

cluster_taxonomy <- ref_clu %>% select(cl_name,category,orf) %>%
 dt_left_join(cl_mmseqs %>% select(-category,-cl_name)) %>%
 dt_left_join(cl_kaiju %>% select(-category,-cl_name))
write_tsv(cluster_taxonomy,path = paste0(opt$output, "/cluster_category_taxonomies.tsv"), col_names = TRUE)

# Majority vote functions
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
  DT.1 <- X[,majority:=majority_vote(mmseqs_tax)$majority, by="cl_name"]
  df <- DT.1 %>% as_tibble() %>% distinct()
}

get_majority <- function(X){

  list_genes <- X %>%                                        # Split into groups by gene-caller-id
    split(.$cl_name)

  maj_l <- mclapply(list_genes,apply_majority, mc.cores=8) # run majority_vote function
  maj_df <- plyr::ldply(maj_l, data.frame) %>%  # bind list rowwise and get distint votes for gene category
    select(cl_name,majority) %>%
    distinct() %>%
    as_tibble()
}

if(dim(cl_mmseqs)[1] > 0) {
  mmseqs_maj <- get_majority(cl_mmseqs %>% select(cl_name,mmseqs_tax) %>% distinct())
}
apply_majority <- function(X){
  DT.1 <- X[,majority:=majority_vote(kaiju_tax)$majority, by="cl_name"]
  df <- DT.1 %>% as_tibble() %>% distinct()
}
kaiju_maj <- get_majority(cl_kaiju %>% select(cl_name,kaiju_tax) %>% distinct())

if(dim(cl_mmseqs)[1] > 0){
  clu_majority_tax <- bind_rows(mmseqs_maj %>% mutate(cl_name=as.character(cl_name)),kaiju_maj)
}else{
  clu_majority_tax <- kaiju_maj
}

# Cluster darkness stats
dpd_res <- fread(opt$clu_dark,
  stringsAsFactors = F, header = F
) %>%
  setNames(c("orf", "dpd_acc", "cl_name", "category"))

dpd_info <- fread(opt$dpd_info, stringsAsFactors = F, header = T, sep = "\t", fill = T) %>%
  mutate(Darkness = Darkness / 100, Disorder = Disorder / 100, Compositional_Bias = Compositional_Bias / 100, Transmembrane = Transmembrane / 100)

cl_dark <- dpd_res %>%
  left_join(dpd_info, by = c("dpd_acc" = "Primary_Accession")) %>%
  group_by(cl_name, category) %>%
  summarise(
    mean_darkness = mean(Darkness), median_darkness = median(Darkness),
    mean_disorder = mean(Disorder), median_disorder = median(Disorder)
  )
write_tsv(cl_dark, path = paste0(opt$output, "/cluster_dpd_perc.tsv"), col_names = TRUE)
cat_dark <- dpd_res %>%
  left_join(dpd_info, by = c("dpd_acc" = "Primary_Accession")) %>%
  group_by(category) %>%
  replace_na(list(Darkness=0, Disorder=0,Compositional_Bias=0,Transmembrane=0)) %>%
  summarise(
    mean_darkness = mean(Darkness), median_darkness = median(Darkness),
    mean_disorder = mean(Disorder), median_disorder = median(Disorder)
  )
write_tsv(cat_dark, path = paste0(opt$output, "/cluster_category_dpd_perc.tsv"), col_names = TRUE)

## Cluster general stats
clu_stats <- clu_stats %>% select(-rep_compl) %>% mutate(cl_name = as.character(cl_name)) %>%
  dt_left_join(clu_majority_tax %>% ungroup() %>% mutate(cl_name = as.character(cl_name))) %>%
  dt_left_join(cl_dark %>% ungroup() %>% mutate(cl_name = as.character(cl_name))) %>%
  dt_left_join(HQ_clusters %>% mutate(is.HQ = TRUE) %>% mutate(cl_name = as.character(cl_name))) %>% distinct()
write_tsv(clu_stats, path = opt$summ_stats, col_names = TRUE)

## Category general stats
cat_stats <- cat_stats %>%
  dt_left_join(cat_compl) %>%
  dt_select(-cl_name, -orf, -rep, -size, -length, -partial) %>%
  distinct() %>%
  dt_left_join(cat_dark)
write_tsv(cat_stats, path = paste0(opt$output, "/only_category_summary_stats.tsv"), col_names = TRUE)
