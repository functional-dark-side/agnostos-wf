#!/usr/bin/env Rscript

# Check if basic packages are installed -----------------------------------

is.installed <- function(pkg){
  is.element(pkg, installed.packages()[,1])
}

if (!is.installed("RSQLite") || !is.installed("dbplyr") || !is.installed("cowplot")){
  cat("We will try to install the packages... (this will be only be done once)\n")
  Sys.sleep(5)
  if (!is.installed("RSQLite")){
    suppressMessages(install.packages("RSQlite", repos = "https://cloud.r-project.org/"))
  }
  if (!is.installed("dbplyr")){
    suppressMessages(install.packages("dbplyr", repos = "https://cloud.r-project.org/"))
  }
  if (!is.installed("cowplot")){
    suppressMessages(install.packages("cowplot", repos = "https://cloud.r-project.org/"))
  }
}

library(tidyverse)
library(data.table)
library(RSQLite)
library(dbplyr)
library(ggrepel)
library(ggsci)
library(cowplot)
library(optparse)

# Script command line options ---------------------------------------------

option_list <- list(
  make_option(c("-d", "--valdb"), type="character", default=NULL,
              help="Validation sqlite results path", metavar="character"),
  make_option(c("-f", "--fval_res"), type="character", default=NULL,
              help="Functional validation results path", metavar="character"),
  make_option(c("-c", "--cval_res"), type="character", default=NULL,
              help="Compositional validation results path", metavar="character"),
  make_option(c("-v", "--val_annot"), type="character", default=1,
              help="Validation annotations", metavar="character"),
  make_option(c("-r", "--val_res"), type="character", default=1,
              help="Validation summary results", metavar="character"),
  make_option(c("-s", "--val_stats"), type="character", default=1,
              help="Validation results stats", metavar="character"),
  make_option(c("-g", "--good"), type="character", default=1,
              help="Validation good clusters", metavar="character"),
  make_option(c("-x", "--fplots"), type="character", default=1,
              help="Funct. validation plots", metavar="character"),
  make_option(c("-y", "--cplots"), type="character", default=1,
              help="Comp. validation plots", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$valdb) | is.null(opt$fval_res) |
    is.null(opt$cval_res) | is.null(opt$val_annot)|
    is.null(opt$val_res) | is.null(opt$val_stats) |
    is.null(opt$good) | is.null(opt$fplots) |
    is.null(opt$cplots) ){
  print_help(opt_parser)
  stop("You need to provide the path to the previous validation step results and output files paths\n", call.=FALSE)
}

setDTthreads(opt$threads)
options(datatable.verbose = FALSE)

db_file = opt$valdb

con <- dbConnect(SQLite(), db_file)

dbGetQuery(con, 'PRAGMA foreign_keys = ON')
dbGetQuery(con, 'PRAGMA auto_vacuum = ON')

# Functional validation
dbGetQuery(con, 'create table if not exists funct_val
           (old_repres text,
           jacc_median_raw numeric,
           jacc_median_sc numeric,
           annot_type text,
           prop_type numeric,
           prop_partial numeric,
           annot_categ text
         )'
)
dbGetQuery(con,"CREATE INDEX if not exists index_repres ON funct_val (old_repres)")

# Load the table with both validation results
funct_val_df <- read_tsv(opt$fval_res, col_names = T,cols(
  rep = col_character(),
  jacc_median_raw = col_double(),
  jacc_median_sc = col_double(),
  type = col_character(),
  prop_type = col_double(),
  prop_partial = col_double(),
  annot_categ = col_character()
)) %>%
 setNames(c("old_repres","jacc_median_raw","jacc_median_sc", "annot_type","prop_type","prop_partial","annot_categ"))

dbWriteTable(con, "funct_val", funct_val_df, append = TRUE)

# Compostional validation
dbGetQuery(con, 'create table if not exists comp_val
           (cl_name text primary key,
           new_repres text,
           n_orfs integer,
           n_vertices integer,
           n_edges integer,
           density numeric,
           cut_w numeric,
           connected logic,
           n_compon integer,
           tr_min_id numeric,
           tr_mean_id numeric,
           tr_median_id numeric,
           tr_max_id numeric,
           raw_min_id numeric,
           raw_mean_id numeric,
           raw_median_id numeric,
           raw_max_id numeric,
           min_len integer,
           mean_len numeric,
           median_len numeric,
           max_len integer,
           rejected integer,
           core integer,
           prop_rejected numeric
           )'
)
dbGetQuery(con,"CREATE INDEX if not exists index_cl_name ON comp_val (cl_name)")

comp_val_df <- read_tsv(opt$cval_res, col_names = F, cols(
  X1 = col_character(),
  X2 = col_character(),
  X3 = col_double(),
  X4 = col_double(),
  X5 = col_double(),
  X6 = col_character(),
  X7 = col_logical(),
  X8 = col_logical(),
  X9 = col_double(),
  X10 = col_double(),
  X11 = col_double(),
  X12 = col_double(),
  X13 = col_double(),
  X14 = col_double(),
  X15 = col_double(),
  X16 = col_double(),
  X17 = col_double(),
  X18 = col_double(),
  X19 = col_double(),
  X20 = col_double(),
  X21 = col_double(),
  X22 = col_double(),
  X23 = col_double(),
  X24 = col_double()
)) %>%
 setNames(c("cl_name","new_repres", "n_orfs","n_vertices","n_edges", "density","cut_w","connected","n_compon","tr_min_id",
            "tr_mean_id","tr_median_id","tr_max_id","raw_min_id","raw_mean_id","raw_median_id","raw_max_id","min_len","mean_len",
            "median_len","max_len","rejected","core","prop_rejected"))

dbWriteTable(con, "comp_val", comp_val_df, append = TRUE)

# Cluster old new representatives
annot_noannot_df <- read_tsv(opt$val_annot, col_names = F,
                             cols(X1 = col_character(),
                                  X2 = col_character(),
                                  X3 = col_character(),
                                  X4 = col_character(),
                                  X5 = col_character()
                             )) %>%
setNames(c("cl_name","new_repres","old_repres","new_repres_annot","funct_annot"))

# Join compositional and functional results
cl_val_df <-  comp_val_df %>%
  inner_join(annot_noannot_df,by=c("cl_name","new_repres")) %>%
  left_join(funct_val_df, by=c("old_repres")) %>%
  dplyr::select(cl_name,new_repres,new_repres_annot,funct_annot,old_repres,n_orfs,n_vertices,n_edges,density,cut_w,connected,
    n_compon,tr_min_id,tr_mean_id,tr_median_id,tr_max_id,raw_min_id,raw_mean_id,raw_median_id,raw_max_id,min_len,mean_len,median_len,max_len,
    rejected,core, prop_rejected,jacc_median_raw,jacc_median_sc,annot_type, prop_type,prop_partial,annot_categ)

#Write results to the database
dbGetQuery(con, 'create table if not exists cluster_val_res
           (cl_name text primary key,
           new_repres text,
           new_repres_annot text,
           funct_annot text,
           old_repres text,
           n_orfs integer,
           n_vertices integer,
           n_edges integer,
           density numeric,
           cut_w numeric,
           connected logic,
           n_compon integer,
           tr_min_id numeric,
           tr_mean_id numeric,
           tr_median_id numeric,
           tr_max_id numeric,
           raw_min_id numeric,
           raw_mean_id numeric,
           raw_median_id numeric,
           raw_max_id numeric,
           min_len integer,
           mean_len numeric,
           median_len numeric,
           max_len integer,
           rejected integer,
           core integer,
           prop_rejected numeric,
           jacc_median_raw numeric,
           jacc_median_sc numeric,
           annot_type text,
           prop_type numeric,
           prop_partial numeric,
           annot_categ text
           )'
)
dbGetQuery(con,"CREATE INDEX if not exists index_cl_name ON cluster_val_res (cl_name)")
dbWriteTable(con, "cluster_val_res", cl_val_df, append = TRUE)


#Write a table
write.table(cl_val_df,file=opt$val_res,col.names=T,row.names=F,quote=F,sep="\t")

#Good and bad clusters stats
p_rej_cl <- plyr::ldply(seq(0,1, 0.01), function(x) {data.frame(threshold = x, clusters = dim(cl_val_df %>% filter(rejected>0) %>% filter(prop_rejected >= x))[1])})
brStick <- function (X) {
  x <- X[[2]]
  m <- 0
  out <- matrix(NA, ncol = 2, nrow = length(x))
  colnames(out) <- c("Observed", "BSM")

  #colnames(out) <- c("% of Variability", "B-Stick Threshold")
  for (i in 1:length(x)) {
    for (k in i:length(x)) {
      m <- m + ((1 / length(x)) * (1 / k))
    }
    out[i, ] <- c((x[i] / sum(x)), m)
    m <- 0
  }
  out <- as_tibble(out) %>% mutate(thresh = X[[1]])
  out_w <- out %>% gather(class, value, -thresh) %>%
    mutate(thresh = as.character(thresh),
           class = fct_rev(class))
  plot <- ggplot(out_w, aes(thresh, value, fill = class)) +
    geom_col(position = "dodge", color = "black", alpha = 0.7) +
    geom_line(aes(group = class, color = class), position=position_dodge(width=0.9)) +
    geom_point(position=position_dodge(width=0.9), colour="black",  shape = 21) +
    theme_light() +
    theme(legend.position = "top",
          legend.title = element_blank()) +
    scale_y_continuous(labels = scales::percent) +
    xlab("Filtering threshold") +
    ylab("Variability")

  h_bsm <- out %>% filter(Observed > BSM) %>% .$thresh

  return(list(bsm_table = out, plot = plot, thresh_by_bsm = h_bsm))
}
lag <- brStick(p_rej_cl)$thresh_by_bsm %>% enframe() %>% mutate(lag = round(value - lag(value), 2)) %>%
  filter(lag > .01) %>% top_n(1) %>% .$name
if (length(lag)!=0){
  rej_threshold <- brStick(p_rej_cl)$thresh_by_bsm[lag - 1]
} else {
  rej_threshold <- brStick(p_rej_cl)$thresh_by_bsm[length(brStick(p_rej_cl)$thresh_by_bsm)]
}
val_stats <- data.frame(total_clusters = dim(cl_val_df)[1],
                        total_orfs = sum(cl_val_df$n_orfs),
                        good_cl_n = dim(cl_val_df %>%
                                          filter(funct_annot=="noannot" & prop_rejected < rej_threshold | funct_annot!="noannot" & prop_rejected<rej_threshold & jacc_median_raw==1))[1],
                        good_cl_orfs = sum(cl_val_df %>%
                                          filter(funct_annot=="noannot" & prop_rejected < rej_threshold | funct_annot!="noannot" & prop_rejected<rej_threshold & jacc_median_raw==1) %>%
                                          select(n_orfs)),
                        bad_cl_n = dim(cl_val_df %>%
                                         filter(prop_rejected >= rej_threshold | funct_annot!="noannot" & jacc_median_raw<1 | funct_annot!="noannot" & is.na(jacc_median_raw)))[1],
                        bad_cl_orfs = sum(cl_val_df %>%
                                         filter(prop_rejected >= rej_threshold | funct_annot!="noannot" & jacc_median_raw<1 | funct_annot!="noannot" & is.na(jacc_median_raw)) %>%
                                           select(n_orfs)),
                        cl_with_rej = dim(cl_val_df %>% filter(rejected>0))[1],
                        orfs_cl_with_rej = sum(cl_val_df %>% filter(rejected>0) %>% dplyr::select(n_orfs)),
                        cl_without_rejected = dim(cl_val_df %>% filter(rejected==0))[1],
                        orfs_cl_without_rej = sum(cl_val_df %>% filter(rejected==0) %>% dplyr::select(n_orfs)),
                        rejected_orfs = sum(cl_val_df$rejected),
                        comp_bad = dim(cl_val_df %>% filter(prop_rejected >= rej_threshold))[1],
                        comp_bad_orfs = sum(cl_val_df %>% filter(prop_rejected >= rej_threshold) %>% dplyr::select(n_orfs)),
                        bad_cl_rej_orfs = sum(cl_val_df %>% filter(prop_rejected >= rej_threshold) %>% dplyr::select(rejected)),
                        comp_good = dim(cl_val_df %>% filter(prop_rejected < rej_threshold))[1],
                        comp_good_orfs = sum(cl_val_df %>% filter(prop_rejected < rej_threshold) %>% dplyr::select(n_orfs)),
                        good_cl_rej_orfs = sum(cl_val_df %>% filter(prop_rejected < rej_threshold) %>% dplyr::select(rejected)),
                        func_good = dim(cl_val_df %>% filter(funct_annot!="noannot" & jacc_median_raw==1))[1],
                        func_bad = dim(cl_val_df %>%
                           filter(funct_annot!="noannot" & jacc_median_raw<1 | funct_annot!="noannot" & is.na(jacc_median_raw)))[1],
                        HA = dim(cl_val_df %>% filter(annot_categ=="HA"))[1],
                        MoDA = dim(cl_val_df %>% filter(annot_categ=="MoDA"))[1],
                        MuDA = dim(cl_val_df %>% filter(annot_categ=="MuDA"))[1],
                        stringsAsFactors = F)
write.table(val_stats, opt$val_stats, col.names = T, row.names = F, quote = F, sep = "\t")

# Good clusters name and representatives
good_cl <- cl_val_df %>%
  filter(funct_annot=="noannot" & prop_rejected < rej_threshold | funct_annot!="noannot" & prop_rejected < rej_threshold & jacc_median_raw==1) %>%
  dplyr::select(cl_name, new_repres, new_repres_annot,funct_annot, old_repres)
write.table(good_cl, opt$good, col.names = T, row.names = F, quote = F, sep = "\t")

# PLOTS
# Functional validation results
cl_val_df <- cl_val_df %>% mutate(funct_annot=ifelse(funct_annot=="new_noannot","noannot","annot"))
cl_val_func <- cl_val_df %>% filter(funct_annot != "noannot")

#plot with trasparent background (change color parameters for white background)
#tiff("clstr_jacc_shingl_raw.tiff", width=2500,height=2000,units = 'px',res = 500, compression = 'lzw', bg = "transparent")
f_raw <- ggplot(cl_val_func, aes(jacc_median_raw)) +
  theme(axis.title.x  = element_text(size = 22),
        axis.text.x  = element_text(size = 20)) +
  theme_bw() + xlab("Similarity") + ylab("Density") +
  geom_density(fill="#18BE8C", colour="#068666", alpha=.8, adjust=0.4) +
  theme(axis.title = element_text(size = 16, colour="white"),
        axis.text = element_text(size=14, colour="white"),
        panel.border = element_rect(colour = "grey"),
        panel.background = element_rect(fill = "transparent",colour = NA), # or theme_blank()
        plot.background = element_rect(fill = "transparent",colour = NA))
#print(all)
#dev.off()
f_sc <- ggplot(cl_val_func, aes(jacc_median_sc)) +
  theme(axis.title.x  = element_text(size = 22),
        axis.text.x  = element_text(size = 20)) +
  theme_bw() + xlab("Similarity") + ylab("Density") +
  geom_density(fill="#18BE8C", colour="#068666", alpha=.8, adjust=0.4) +
  theme(axis.title = element_text(size = 16, colour="black"),
        axis.text = element_text(size=14, colour="black"))

save(f_raw,f_sc,file=opt$fplots)

# Compositional validation results
#With numbers in the labels
rej_desc <- data.frame(class = c(paste0("With rejected (", dim(cl_val_df %>% filter(rejected>0))[1], ")"),
                                 paste0("Without rejected (",dim(cl_val_df %>% filter(rejected==0))[1], ")")),
                                 num = c(dim(cl_val_df %>% filter(rejected>0))[1], dim(cl_val_df %>% filter(rejected==0))[1]))

p_desc <- ggplot(rej_desc, aes(class, num)) +
  geom_bar(stat = "identity") +
  scale_y_continuous(labels = scales::comma) +
  ylab("Number of clusters") +
  xlab("") +
  ggtitle("Clusters with rejected ORFs") +
  theme_light() +
  theme(plot.title = element_text(size=13),
        axis.text = element_text(size=11),
        axis.title = element_text(size=12))

cl_val_df$funct_annot <- factor(cl_val_df$funct_annot, levels = c("annot", "noannot"))
cl_val_df$annot_categ <- factor(cl_val_df$annot_categ, levels = c("HA", "MoDA", "MuDA", NA))

class_names <- c(
  `annot` = "Annotated",
  `noannot` = "Not annotated"
)

# Prop. rejected ORFs vs number of clusters
# Number of remaining clusters after applying different thresholds
# based on the number of bad aligned ORFs per cluster.
#With BSM
bsm <- brStick(p_rej_cl)$plot[[1]] %>% filter(class=="BSM") %>%
  mutate(clusters=value*sum(p_rej_cl$clusters)) %>%
  select(-value) %>% rename(threshold=thresh)
p_rej_cl_bsm <- rbind(p_rej_cl %>% mutate(class="Observed",threshold=as.character(threshold)), bsm) %>%
  as.tibble() %>% mutate(class=as.factor(class))

p_rej_clp <- ggplot(p_rej_cl_bsm %>% as.tibble(), aes(threshold, clusters, fill=class)) +
  geom_col(position ="dodge", color = "black", alpha = 0.7) +
  geom_line(aes(group = class, color = class), position=position_dodge(width=0.9)) +
  geom_point(position=position_dodge(width=0.9), colour="black",  shape = 21) +
  ggrepel::geom_label_repel(data=(p_rej_cl_bsm %>% filter(threshold==rej_threshold) %>% mutate(label=ifelse(class=="Observed",paste0("T=",threshold,"; N=",clusters),NA))),
                            aes(label = label), size = 4, position = position_dodge(width=0.9),
                            box.padding = unit(0.4, "lines"),point.padding = unit(0.45, "lines"),show.legend = F,fill="white") +
  theme_light() +
  xlab("Proportion of non-homolog ORFs per cluster") +
  ylab("Number of clusters") +
  scale_color_manual(values=c("#94bacc","#d95472")) +
  scale_fill_manual(values=c("#94bacc","#d95472")) +
  scale_y_continuous(labels = scales::comma) +
  scale_x_discrete(breaks = unique(p_rej_cl_bsm$threshold[seq(1, length(unique(p_rej_cl_bsm$threshold)), by = 10)])) +
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=13),
        legend.title = element_blank())

# Prop. rejected ORFs vs average ORFs similarity (divided by annotation type)
# Relationship between the proportion of rejected ORFs identified by LEON-BIS
# and the average ORF similarity in each cluster (In red rejected clusters).
p_rej_simil <- ggplot(cl_val_df, aes(raw_mean_id/100, prop_rejected)) +
  geom_jitter(alpha = 0.5) +
  geom_jitter(data = cl_val_df %>% filter(prop_rejected >= rej_threshold), aes(raw_mean_id/100, prop_rejected, size = n_orfs), color = "#C84359", alpha = 0.5) +
  geom_rug(data = cl_val_df, aes(color=funct_annot), alpha = 1/2, position = "jitter") +
  theme_light() +
  ylab("Proportion of rejected ORFs per cluster") +
  xlab("Average ORF similarity per cluster") +
  scale_x_continuous(labels = scales::percent) +
  scale_y_continuous(labels = scales::percent) +
  facet_wrap(~annot_categ, scales = "free", nrow = 1) +
  ggsci::scale_color_jco(name = "Cluster type", labels = c("Annot.","Not annot.")) +
  scale_size_continuous(name = "Number of ORFs") +
  theme(axis.title = element_text(size = 13),
      axis.text = element_text(size = 11),
      strip.text = element_text(size=13),
      legend.text = element_text(size=11),
      legend.key.size = unit(1,"cm"),
      legend.title = element_text(size=13))

# Size distribution of the kept and rejected clusters.
p_size_rej <- ggplot(cl_val_df %>% filter(rejected>0,prop_rejected >= 0.1), aes(n_orfs, fill = "Rejected")) +
  geom_histogram(data = cl_val_df %>% filter(prop_rejected < 0.1), aes(n_orfs, fill = "Kept"), color = "white", size = 0.1) +
  geom_histogram(color = "white", size = 0.1) +
  theme_light() +
  xlab("Cluster size (log10)") +
  ylab("Number of clusters") +
  scale_x_log10() +
  #scale_y_log10() +
  scale_y_continuous(labels = scales::comma) +
  scale_fill_manual(values = c("#4A4A4A", "#F0D999"), name = "") +
  theme(legend.position = c(0.88, 0.85),
        axis.title = element_text(size=13),
        axis.text = element_text(size=11))

# Plots together
p_panel <- ggdraw() +
            draw_plot(p_desc, x = 0, y = .5, width = .30, height = .5) +
            draw_plot(p_rej_clp, x = .30, y = .5, width = .40, height = .5) +
            draw_plot(p_size_rej, x = .70, y = .5, width = .30, height = .5) +
            draw_plot(p_rej_simil, x = 0, y = 0, width = 1, height = .5) +
            draw_plot_label(label = c("a", "b", "c", "d"), size = 15,
                  x = c(0, 0.30, 0.70, 0), y = c(1, 1, 1, 0.53))

save(cl_val_df,p_rej_cl_bsm,p_desc,p_size_rej, p_rej_clp, p_rej_simil, p_panel, file=opt$cplots)
