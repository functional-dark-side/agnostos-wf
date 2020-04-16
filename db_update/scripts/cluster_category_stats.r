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

library(tidyverse)
library(data.table)
library(maditr)
library(entropy)
library(ggridges)
library(optparse)

# Script command line options ---------------------------------------------

option_list <- list(
  make_option(c("-r", "--ref_clu"),
    type = "character", default = NULL,
    help = "Refined set of clusters", metavar = "character"
  ),
  make_option(c("-c", "--clu_categ"),
    type = "character", default = NULL,
    help = "Cluster ids-categories", metavar = "character"
  ),
  make_option(c("-t", "--clu_tax"),
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
  is.null(opt$clu_tax) | is.null(opt$clu_dark) |
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
cl_tax <- fread(opt$clu_tax,
  stringsAsFactors = F, header = F, fill = T, sep="\t") %>%
  setNames(c("orf", "cl_name", "category", "taxid", "rank", "classific", "lineage"))
# Filter out the unclassified ones and split the taxonomic levels in the lineage
cl_tax <- cl_tax %>%
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
                          TRUE ~ domain))

# Taxonomic homogeneity
homo_tax <- cl_tax %>%
  select(-orf, -taxid, -classific, -rank, -lineage) %>%
  group_by(cl_name, category) %>%
  summarise(
    homo_d = length(unique(domain[!is.na(domain)])),
    homo_p = length(unique(phylum[!is.na(phylum)])),
    homo_o = length(unique(`order`[!is.na(`order`)])),
    homo_c = length(unique(class[!is.na(class)])),
    homo_f = length(unique(family[!is.na(family)])),
    homo_g = length(unique(genus[!is.na(genus)])),
    homo_s = length(unique(species[!is.na(species)]))
  )

### Plot:
plot_tax <- gather(homo_tax, ranks, n, homo_p:homo_s) %>%
  mutate(ranks = case_when(
    ranks == "homo_p" ~ "Phylum",
    ranks == "homo_c" ~ "Class",
    ranks == "homo_o" ~ "Order",
    ranks == "homo_f" ~ "Family",
    ranks == "homo_g" ~ "Genus",
    TRUE ~ "Species"
  )) %>%
  mutate(n = as.numeric(n), ranks = as.factor(ranks)) %>%
  as_tibble()
plot_tax$category <- factor(as.factor(plot_tax$category), levels = c("K", "KWP", "GU", "EU"))
plot_tax$ranks <- factor(as.factor(plot_tax$ranks), levels = c("Phylum", "Class", "Order", "Family", "Genus", "Species"))
# joy_div plot
taxp <- ggplot(plot_tax %>% filter(n <= 50 & n > 0), aes(y = category, x = n, fill = category)) +
  geom_density_ridges(scale = 4, bandwidth = .3, alpha = .8, size = .3) +
  scale_x_continuous(breaks = c(1, 5, 10, 15, 20)) +
  coord_cartesian(xlim = c(0, 20)) +
  facet_wrap(. ~ ranks) +
  scale_fill_manual(values = c("#233B43", "#556c74", "#65ADC2", "#E84646")) + ## 4B636A K
  theme_light() +
  xlab("Number of different taxonomies inside each cluster") +
  ylab("") +
  theme(
    axis.text.x = element_text(size = 6),
    axis.text.y = element_text(size = 6),
    axis.title.x = element_text(size = 7),
    legend.text = element_text(size = 6),
    legend.title = element_blank(),
    legend.key.size = unit(.39, "cm"),
    strip.text = element_text(size = 7, colour = "white"),
    strip.background = element_rect(fill = "#273133")
  )
save(taxp, file = paste0(opt$output, "/cluster_category_taxonomy_plot.rda"))

# Calculate cluster taxonomic entropy
d <- cl_tax %>%
  group_by(cl_name, category, domain) %>%
  summarise(n = n()) %>%
  mutate(df = n / sum(n)) %>%
  group_by(cl_name, category) %>%
  mutate(de = entropy.empirical(df, unit = "log2")) %>%
  select(1, 2, 6) %>%
  distinct()
p <- cl_tax %>%
  group_by(cl_name, category, phylum) %>%
  summarise(n = n()) %>%
  mutate(pf = n / sum(n)) %>%
  group_by(cl_name, category) %>%
  mutate(pe = entropy.empirical(pf, unit = "log2")) %>%
  select(1, 2, 6) %>%
  distinct()
o <- cl_tax %>%
  group_by(cl_name, category, order) %>%
  summarise(n = n()) %>%
  mutate(of = n / sum(n)) %>%
  group_by(cl_name, category) %>%
  mutate(oe = entropy.empirical(of, unit = "log2")) %>%
  select(1, 2, 6) %>%
  distinct()
c <- cl_tax %>%
  group_by(cl_name, category, class) %>%
  summarise(n = n()) %>%
  mutate(cf = n / sum(n)) %>%
  group_by(cl_name, category) %>%
  mutate(ce = entropy.empirical(cf, unit = "log2")) %>%
  select(1, 2, 6) %>%
  distinct()
f <- cl_tax %>%
  group_by(cl_name, category, family) %>%
  summarise(n = n()) %>%
  mutate(ff = n / sum(n)) %>%
  group_by(cl_name, category) %>%
  mutate(fe = entropy.empirical(ff, unit = "log2")) %>%
  select(1, 2, 6) %>%
  distinct()
g <- cl_tax %>%
  group_by(cl_name, category, genus) %>%
  summarise(n = n()) %>%
  mutate(gf = n / sum(n)) %>%
  group_by(cl_name, category) %>%
  mutate(ge = entropy.empirical(gf, unit = "log2")) %>%
  select(1, 2, 6) %>%
  distinct()
s <- cl_tax %>%
  group_by(cl_name, category, species) %>%
  summarise(n = n()) %>%
  mutate(sf = n / sum(n)) %>%
  group_by(cl_name, category) %>%
  mutate(se = entropy.empirical(sf, unit = "log2")) %>%
  select(1, 2, 6) %>%
  distinct()
tax_entropy <- d %>%
  left_join(p) %>%
  left_join(c) %>%
  left_join(o) %>%
  left_join(f) %>%
  left_join(g) %>%
  left_join(s)
rm(d,p, c, o, f, g, s)
write_tsv(tax_entropy, path = paste0(opt$output, "/cluster_category_taxonomy_entropy.tsv"), col_names = TRUE)

## Retrieve the prevalent taxa for each cluster
prev_tax <- cl_tax %>%
  group_by(cl_name, category, domain, phylum, class, order, family, genus, species) %>%
  count() %>%
  ungroup() %>%
  group_by(cl_name, category) %>%
  arrange(desc(n)) %>%
  slice(1) %>%
  select(-n)
write_tsv(prev_tax, path = paste0(opt$output, "/cluster_category_prevalent_tax.tsv"), col_names = TRUE)

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
clu_stats <- clu_stats %>% select(-rep_compl) %>%
  dt_left_join(tax_entropy %>% ungroup() %>% mutate(cl_name = as.character(cl_name))) %>%
  dt_left_join(prev_tax %>% ungroup() %>% mutate(cl_name = as.character(cl_name))) %>%
  dt_left_join(cl_dark %>% ungroup() %>% mutate(cl_name = as.character(cl_name))) %>%
  dt_left_join(HQ_clusters %>% mutate(is.HQ = TRUE)) %>% distinct()
write_tsv(clu_stats, path = opt$summ_stats, col_names = TRUE)

## Category general stats
cat_stats <- cat_stats %>%
  dt_left_join(cat_compl) %>%
  dt_select(-cl_name, -orf, -rep, -size, -length, -partial) %>%
  distinct() %>%
  dt_left_join(tax_entropy %>% ungroup() %>% select(-cl_name) %>% group_by(category) %>% summarise_all(sum)) %>%
  dt_left_join(cat_dark)
write_tsv(clu_stats, path = paste0(opt$output, "/only_category_summary_stats.tsv"), col_names = TRUE)
