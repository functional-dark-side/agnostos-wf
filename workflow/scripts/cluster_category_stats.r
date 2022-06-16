#!/usr/bin/env Rscript

# Check if basic packages are installed -----------------------------------

is.installed <- function(pkg) {
  is.element(pkg, installed.packages()[, 1])
}

if (!is.installed("entropy") || !is.installed("ggridges")) {
  cat("We will try to install the packages... (this will be only be done once)\n")
  Sys.sleep(5)
  if (!is.installed("entropy")) {
    suppressMessages(install.packages("entropy", repos = "http://cran.us.r-project.org"))
  }
  if (!is.installed("ggridges")) {
    suppressMessages(install.packages("ggridges", repos = "http://cran.us.r-project.org"))
  }
}

library(data.table)
library(tidyverse)
library(entropy)
library(RSQLite)
library(dbplyr)
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
  make_option(c("-f", "--mutantDB"),
    type = "character", default = NULL,
    help = "feba.db for mutants", metavar = "character"
  ),
  make_option(c("-u", "--mutants"),
    type = "character", default = NULL,
    help = "mutants results", metavar = "character"
  ),
  make_option(c("-e", "--eggnog"),
    type = "character", default = NULL,
    help = "eggNOG annotations", metavar = "character"
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
  is.null(opt$kaiju_tax) | is.null(opt$clu_dark) |
  is.null(opt$dpd_info) | is.null(opt$compl) |
  is.null(opt$hq_clu) | is.null(opt$mutants) |
  is.null(opt$eggnog) | is.null(opt$mutantDB) |
  is.null(opt$summ_stats) |
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

ref_clu <- ref_clu %>% inner_join(category)

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
  select(cl_name, partial, p) %>%
  spread(partial, p, fill = 0)

if ("01" %in% colnames(clu_compl) & "11" %in% colnames(clu_compl)){
  clu_compl <- clu_compl %>%
    rename(p00 = `00`, p10 = `10`, p01 = `01`, p11 = `11`)
} else if ("11" %in% colnames(clu_compl)){
  clu_compl <- clu_compl %>%
    mutate(p00 =`00`, p01=0, p10=0, p11= `11`) %>% select(-`00`,-`11`)
} else {
  clu_compl <- clu_compl %>%
    mutate(p00 =`00`, p01=0, p10=0, p11=0) %>% select(-`00`)
}

write_tsv(clu_compl, path=opt$compl, col_names = T)

clu_stats <- clu_stats %>% left_join(clu_compl) %>%
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
  select(category, partial, p) %>%
  distinct() %>%
  spread(partial, p, fill = 0)

  if ("01" %in% colnames(cat_compl) & "11" %in% colnames(cat_compl)){
    cat_compl <- cat_compl %>%
      rename(p00 = `00`, p10 = `10`, p01 = `01`, p11 = `11`)
  } else if ("11" %in% colnames(cat_compl)){
    cat_compl <- cat_compl %>%
      mutate(p00 =`00`, p01=0, p10=0, p11= `11`) %>% select(-`00`,-`11`)
  } else {
    cat_compl <- cat_compl %>%
      mutate(p00 =`00`, p01=0, p10=0, p11=0) %>% select(-`00`)
  }

# Cluster taxonomy stats
cl_kaiju <- fread(opt$kaiju_tax,
  header = F, fill=T, sep="\t") %>%
  setNames(c("orf","domain","phylum","class","order","family","genus","species"))

cl_kaiju <- cl_kaiju %>% left_join(ref_clu %>% select(cl_name,category,orf))

cl_kaiju_tbl <- cl_kaiju %>% unite(kaiju_tax, domain:species,sep=";")

write_tsv(cl_kaiju_tbl,
  path = paste0(opt$output, "/cluster_kaiju_taxonomy.tsv"), col_names = TRUE)

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

get_majority <- function(X){

  list_genes <- X %>%                                        # Split into groups by gene-caller-id
    split(.$cl_name)

  maj_l <- mclapply(list_genes,apply_majority, mc.cores=8) # run majority_vote function
  maj_df <- plyr::ldply(maj_l, data.frame) %>%  # bind list rowwise and get distint votes for gene category
    select(cl_name,majority) %>%
    distinct() %>%
    as_tibble()
}

apply_majority <- function(X){
  DT.1 <- X[,majority:=majority_vote(kaiju_tax)$majority, by="cl_name"]
  df <- DT.1 %>% as_tibble() %>% distinct()
}
kaiju_maj <- get_majority(cl_kaiju_tbl %>% select(cl_name,kaiju_tax) %>% distinct())

# Cluster taxonomic entropy
# Calculate cluster taxonomic entropy
d <- cl_kaiju %>% group_by(cl_name, category,domain) %>%
  summarise (n=n()) %>%
  mutate(df = n / sum(n)) %>%
  group_by(cl_name, category) %>%
  mutate(de = entropy.empirical(df, unit="log2")) %>%
  select(1,2,6) %>% distinct()
p <- cl_kaiju %>% group_by(cl_name, category,phylum) %>%
  summarise (n=n()) %>%
  mutate(pf = n / sum(n)) %>%
  group_by(cl_name, category) %>%
  mutate(pe = entropy.empirical(pf, unit="log2")) %>%
  select(1,2,6) %>% distinct()
o <- cl_kaiju %>% group_by(cl_name, category,order) %>%
  summarise (n=n()) %>%
  mutate(of = n / sum(n)) %>%
  group_by(cl_name, category) %>%
  mutate(oe = entropy.empirical(of, unit="log2")) %>%
  select(1,2,6) %>% distinct()
c <- cl_kaiju %>% group_by(cl_name, category,class) %>%
  summarise (n=n()) %>%
  mutate(cf = n / sum(n)) %>%
  group_by(cl_name, category) %>%
  mutate(ce = entropy.empirical(cf, unit="log2")) %>%
  select(1,2,6) %>% distinct()
f <- cl_kaiju %>% group_by(cl_name, category,family) %>%
  summarise (n=n()) %>%
  mutate(ff = n / sum(n)) %>%
  group_by(cl_name, category) %>%
  mutate(fe = entropy.empirical(ff, unit="log2")) %>%
  select(1,2,6) %>% distinct()
g <- cl_kaiju %>% group_by(cl_name, category,genus) %>%
  summarise (n=n()) %>%
  mutate(gf = n / sum(n)) %>%
  group_by(cl_name, category) %>%
  mutate(ge = entropy.empirical(gf, unit="log2")) %>%
  select(1,2,6) %>% distinct()
s <- cl_kaiju %>% group_by(cl_name, category,species) %>%
  summarise (n=n()) %>%
  mutate(sf = n / sum(n)) %>%
  group_by(cl_name, category) %>%
  mutate(se = entropy.empirical(sf, unit="log2")) %>%
  select(1,2,6) %>% distinct()

homo_tax <- cl_kaiju %>% select(-orf,-category) %>%
  group_by(cl_name) %>%
  summarise(Num_d=length(unique(domain[!is.na(domain)])),
            Num_p=length(unique(phylum[!is.na(phylum)])),
            Num_o=length(unique(`order`[!is.na(`order`)])),
            Num_c=length(unique(class[!is.na(class)])),
            Num_f=length(unique(family[!is.na(family)])),
            Num_g=length(unique(genus[!is.na(genus)])),
            Num_s=length(unique(species[!is.na(species)])))

tax_entropy <- d %>% left_join(p) %>% left_join(c) %>% left_join(o) %>%
  left_join(f) %>% left_join(g) %>% left_join(s) %>% left_join(homo_tax)
rm(d,p,c,o,f,g,s)
write.table(tax_entropy, paste0(opt$output, "/cluster_taxonomic_entropy.tsv"),
 row.names = F, quote =  F, sep = "\t")

# EggNOG annotations
nog_res <- fread(opt$eggnog,
  stringsAsFactors = F, header = F
) %>%
  setNames(c("orf", "cog_categ","desc","cl_name", "category"))

nog_res <- nog_res %>%
  separate_rows(cog_categ,sep=",") %>%
  filter(cog_categ!="") %>%
  distinct() %>%
  select(cl_name,category,orf,cog_categ,desc) %>%
  rename(description=desc)

write.table(nog_res, paste0(opt$output, "/COG_categories_per_gene.tsv"),
   row.names = F, quote =  F, sep = "\t")

# Eggnog annotations entropy
eggnog_entropy <- nog_res %>%
  group_by(cl_name,category,cog_categ) %>%
  summarise(n=n()) %>%
  mutate(cat_f = n / sum(n)) %>%
  group_by(cl_name, category) %>%
  mutate(cat_entropy = entropy.empirical(cat_f, unit="log2")) %>%
  select(cl_name,category,cat_entropy) %>% distinct()

  write.table(nog_res, paste0(opt$output, "/cluster_COG_category_entropy.tsv"),
     row.names = F, quote =  F, sep = "\t")

# Cluster mutant information
mutants_db <- dbConnect(drv=SQLite(), dbname=opt$mutantDB)
experim <-  tbl(mutants_db, "Experiment") %>%
  collect(n = Inf) %>%
  select(orgId, expName, expDesc, expGroup) %>% distinct()
gene_fit <- tbl(mutants_db, "GeneFitness") %>% collect(n = Inf) %>%
    inner_join(experim) %>%
    mutate(mutantId=paste(orgId,locusId,sep=":"))
mutants <- fread(opt$mutants,
   stringsAsFactors = F, header = F, sep="\t") %>%
   setNames(c("orf", "mutantId", "cl_name", "category"))

mutants <- mutants %>% left_join(gene_fit)

write.table(mutants, paste0(opt$output, "/mutant_phenotypes_per_gene.tsv"),
   row.names = F, quote =  F, sep = "\t")

mutants <- mutants %>%
  mutate(class = ifelse(category == "K" | category == "KWP", "Known", "Unknown")) %>%
  select(orgId, locusId, cl_name, category, class, expDesc, expGroup, fit) %>%
  group_by(orgId, locusId, cl_name, category, class, expDesc, expGroup) %>%
  summarise(fit = mean(fit))

write.table(mutants, paste0(opt$output, "/cluster_mutant_phenotypes.tsv"),
     row.names = F, quote =  F, sep = "\t")

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
  left_join(kaiju_maj %>% ungroup() %>% mutate(cl_name = as.character(cl_name))) %>%
  left_join(cl_dark %>% ungroup() %>% mutate(cl_name = as.character(cl_name))) %>%
  left_join(HQ_clusters %>% mutate(is.HQ = TRUE) %>% mutate(cl_name = as.character(cl_name))) %>% distinct()
write_tsv(clu_stats, path = opt$summ_stats, col_names = TRUE)
