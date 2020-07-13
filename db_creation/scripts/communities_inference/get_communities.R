#!/usr/bin/env Rscript

# Check if basic packages are installed -----------------------------------

is.installed <- function(pkg){
  is.element(pkg, installed.packages()[,1])
}

if (!is.installed("crayon") || !is.installed("unixtools") || !is.installed("config")){
  cat("We will try to install the packages crayon, optparse and config... (this will be only be done once)\n")
  Sys.sleep(5)
  if (!is.installed("crayon")){
    suppressMessages(install.packages("crayon", repos = "http://cran.us.r-project.org"))
  }
  if (!is.installed("config")){
    suppressMessages(install.packages("config", repos = "http://cran.us.r-project.org"))
  }
  if (!is.installed("unixtools")){
    suppressMessages(install.packages("unixtools",,"http://rforge.net/",type="source"))
  }
}

suppressMessages(library(crayon))
suppressMessages(library(optparse))


# Check if packaged are installed -----------------------------------------

cat("\nChecking if all packages are installed...\n\n")


needed <- c("igraph", "tidyverse", "pbmcapply", "maditr", "data.table",
            "tidygraph", "bigreadr", "unixtools", "stringi", "foreach",
            "doSNOW", "itertools", "parallel","reshape2")

missing_package <- FALSE
# For loop to run through each of the packages
for (p in 1:length(needed)){
  if(is.installed(needed[p])){
    cat(sprintf("%-10s: %s", needed[p], green("Installed\n")))
  }else{
    cat(sprintf("%-10s: %s", needed[p], red("Not installed\n")))
    missing_package <- TRUE
  }
}

quit_not_installed <- function(){
  cat("\nMissing packages, please install them.\n")
  quit(save = "no", status = 1)
}

if (missing_package) {
  quit_not_installed()
}else{
  cat("\nAll packages installed.\n")
}

Sys.sleep(2)
system("clear")



# Load libraries ----------------------------------------------------------
cat("Loading libraries...")
silent <- suppressMessages(lapply(needed, function(X) {require(X, character.only = TRUE)}))
rm(silent)
cat(" done\n\n")

# Some functions ----------------------------------------------------------

msg <- function(X){
  cat(crayon::white(paste0("[",format(Sys.time(), "%T"), "]")), X)
}

msg_sub <- function(X){
  cat(crayon::white(paste0("  [",format(Sys.time(), "%T"), "]")), X)
}

# Script command line options ---------------------------------------------

option_list <- list(
  make_option(c("-c", "--config"), type="character", default=NULL,
              help="YAML config file", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list);
opt <- parse_args(opt_parser);


if (is.null(opt$config)){
  print_help(opt_parser)
  stop("At least one arguments must be supplied (YAML config file).n", call.=FALSE)
}

msg("Getting config parameters...")
cfg <- config::get(file = opt$config, use_parent = FALSE)
cat(" done\n")

# Create a new environment
lo_env <- new.env()

# Define data.table parameters
options(datatable.optimize=Inf)
options(datatable.verbose = cfg$dt_verbose)
setDTthreads(cfg$dt_cores)
# getDTthreads()



# Read data ---------------------------------------------------------------

#load("/bioinf/projects/megx/UNKNOWNS/chiara/architect/communities/knowns_DA.rda",
#    verbose = TRUE)


# Number of k clusters
# 1,122,105
msg("Reading cluster categories...")
cl_cat <- fread(file = cfg$cl_cat, header = FALSE, showProgress = FALSE, nThread = cfg$dt_cores)
setnames(cl_cat, names(cl_cat), c("cl_name", "category"))
cat(" done\n")

msg("Reading cluster completeness...")
cl_compl <- fread(file = cfg$cl_comp, header = TRUE, showProgress = FALSE, nThread = cfg$dt_cores)
cat(" done\n")

msg("Reading PFAM31 clan information...")
p_clan <- fread(file = cfg$p_clan, header = FALSE, showProgress = FALSE, nThread = cfg$dt_cores) %>%
  dt_select(V2, V4) %>%
  setnames(c("V2", "V4"), c("clan", "pfam"))
cat(" done\n")

msg("Reading simplified domain architecture information...")
p_doms <- fread(input = cfg$p_doms, header = TRUE, showProgress = FALSE, nThread = cfg$dt_cores) %>%
  mutate(cl_name = as.character(cl_name))
cat(" done\n\n")

k_comp <- p_doms %>% dt_filter(cl_name %in% (cl_cat %>% dt_filter(category == "K") %>% .$cl_name))


tmp <- file.path(cfg$wd, "tmp")
msg(paste0("Creating TMP folder at ", tmp, "..."))
if (dir.exists(tmp)) unlink(tmp, recursive = TRUE, force = TRUE)
dir.create(tmp, recursive = TRUE)
# Set a new tempdir
set.tempdir(tmp)
cat(" done\n")

time_string <- paste(format(Sys.Date()), format(Sys.time(), "%H%M%S"), sep = "-")
results <- file.path(cfg$wd, time_string)
msg(paste0("Creating RESULTS folder at ", results, "..."))
if (dir.exists(results)) unlink(results, recursive = TRUE, force = TRUE)
dir.create(results, recursive = TRUE)
cat(" done\n\n")

msg("Reading Known HHBLITS all-vs-all results...")
lo_env$k_hhblits_all <- big_fread1(cfg$hhblits_results,
                                   print_timings = FALSE,
                                   every_nlines = 2e9,
                                   data.table = TRUE,
                                   .combine = list,
                                   nThread = cfg$dt_cores,
                                   header = FALSE
)
k_hhblits_all_nrow <- map(unlist(lo_env$k_hhblits_all, recursive = FALSE), nrow) %>% unlist() %>% sum()
cat(paste(" Read", green(scales::comma(k_hhblits_all_nrow)), "entries\n"))

# Let's filter for prob >= 50, and cov > 0.6
msg("Filtering Known HHBLITS all-vs-all results with P>=50 and cov > 0.6...")
lo_env$k_hhblits <-  data.table::rbindlist(map(unlist(lo_env$k_hhblits_all, recursive = FALSE), function(X){
  X <- X %>% dt_filter(V1 != V2, V3 >= 50, V13 > 0.6, V14 > 0.6)
  setnames(X, names(X),
           c("cl_name1", "cl_name2", "probability", "e-value", "Score",
             "Cols", "q_start", "q_stop", "t_start", "t_stop", "q_len",
             "t_len", "q_cov", "t_cov"))
  X %>% maditr::dt_mutate(cl_name1 = as.character(cl_name1),
                          cl_name2 = as.character(cl_name2),
                          score_col = Score/Cols)
}))
# How many clusters do we have (removed self-hits)
# 1,122,105
# 1,100,297 p50;c0.6
unique_hhblits_k_cl <- c(lo_env$k_hhblits$cl_name1, lo_env$k_hhblits$cl_name2) %>% unique()
cat(paste0(" Found ", scales::comma(length(unique_hhblits_k_cl)), " clusters\n"))

# We are missing 21,808
# Which ones
# Find missing clusters ---------------------------------------------------
msg("Identifying missing clusters...")
k_missing_ids <- setdiff(cl_cat %>% dt_filter(category == "K") %>% .$cl_name, unique_hhblits_k_cl)
lo_env$k_hhblits_missing <-  data.table::rbindlist(map(unlist(lo_env$k_hhblits_all, recursive = FALSE), function(X){
  X <- X %>% dt_filter(V1 %in% k_missing_ids | V2 %in% k_missing_ids) %>%
    dt_filter(V1 != V2)
  setnames(X, names(X),
           c("cl_name1", "cl_name2", "probability", "e-value", "Score",
             "Cols", "q_start", "q_stop", "t_start", "t_stop", "q_len",
             "t_len", "q_cov", "t_cov"))
  X %>% maditr::dt_mutate(cl_name1 = as.character(cl_name1),
                          cl_name2 = as.character(cl_name2),
                          score_col = Score/Cols)
}))
cat(paste0("Found ", scales::comma(length(k_missing_ids))), "clusters missing\n")


# Prepare data for network analysis ---------------------------------------
msg("Inferring graph from HHBLITS results...")
lo_env$k_hhb_bh_score <- lo_env$k_hhblits %>% dt_select(cl_name1, cl_name2, score_col) %>% as_tibble()

k_hh_g <- lo_env$k_hhblits %>%
  dt_select(cl_name1, cl_name2, score_col) %>%
  as_tibble() %>%
  rename(weight = score_col) %>%
  igraph::graph_from_data_frame(directed = FALSE) %>%
  igraph::simplify(remove.multiple = TRUE, remove.loops = TRUE, edge.attr.comb = list("max")) %>%
  as_tbl_graph()

k_hh_g <- k_hh_g %>%
  activate(nodes) %>%
  #select(-archit) %>%
  inner_join(k_comp %>% as_tibble() %>% rename(name = cl_name) %>% mutate(name = as.character(name)) %>% select(name, archit, original), by = "name") %>%
  filter(!node_is_isolated())

cat(paste0(" Graph contains ", scales::comma(vcount(k_hh_g)), " vertices and ", scales::comma(ecount(k_hh_g)), " edges\n"))

# Identify communities using MCL ------------------------------------------
source(cfg$graph_lib)

if (cfg$k_g_mcl_list != ""){
  msg("Loading already computed MCL clusters...")
  load(cfg$k_g_mcl_list)
  inflation_list <- seq(cfg$mcl_inflation_min, cfg$mcl_inflation_max, cfg$mcl_inflation_step)
  cat(" done\n")
}else{
  msg("Preparing graph for MCL clustering...")
  gx <- k_hh_g
  lo_env$w <- E(gx)$weight
  lo_env$w  <- lo_env$w - min(lo_env$w ) + 0.001
  E(gx)$weight <- lo_env$w
  cat(" done\n")

  msg(paste0("Running MCL clustering with inflation values from ", cfg$mcl_inflation_min, " to ", cfg$mcl_inflation_max, " with ", cfg$mcl_inflation_step, " steps...\n"))
  inflation_list <- seq(cfg$mcl_inflation_min, cfg$mcl_inflation_max, cfg$mcl_inflation_step)
  mcl_bin <- cfg$mcl_bin
  k_g_mcl_list <- parallel::mclapply(inflation_list,
                                        optimal_mcl,
                                        G = k_hh_g,
                                        mcl_cores = cfg$mcl_cores,
                                        #ignore.interactive = TRUE,
                                        Gx = gx,
                                        #max.vector.size = 1e+07,
                                        mc.cores = cfg$mcl_jobs,
                                        mc.cleanup = TRUE,
                                        mc.silent = TRUE,
                                        tmp_dir = tmp)

  #k_g_mcl_list <- lapply(inflation_list, optimal_mcl, G = k_hh_g, Gx = gx, mcl_cores = cfg$mcl_cores, tmp_dir = tmp)
  names(k_g_mcl_list) <- inflation_list

  k_g_mcl_list_file <- file.path(results, paste0("k_g_mcl_list_", time_string, ".Rda"))
  msg(paste0("Saving computed MCL communities in ", k_g_mcl_list_file, "..."))
  save(k_g_mcl_list, file = k_g_mcl_list_file, compress = FALSE)
  cat(" done\n")
}

# Contract identified communities -----------------------------------------
if (cfg$k_gc != ""){
  msg("Loading already contracted MCL communities...")
  load(cfg$k_gc)
  cat(" done\n")
}else{
  msg("Contracting MCL communities...\n")
  k_gc <- parallel::mclapply(k_g_mcl_list,
                                contract_graphs,
                                G = k_hh_g,
                                #max.vector.size = 3e+09,
                                mc.cores = cfg$max_gc_jobs,
                                #ignore.interactive = TRUE,
                                mc.cleanup = TRUE,
                                mc.silent = TRUE
  )
  names(k_gc) <- inflation_list
  msg("Contracting MCL communities... done\n")

  k_gc_file <- file.path(results, paste0("k_gc_", time_string, ".Rda"))
  msg(paste0("Saving contracted MCL communities results in ", k_gc_file, "..."))
  save(k_gc, file = k_gc_file, compress = FALSE)
  cat(" done\n")
}

msg("Computing modularity for each MCL community result...")
k_hh_g_modularity <- map_df(k_g_mcl_list, function(X){
  vnames <- V(k_hh_g)$name
  tibble(modularity = modularity(k_hh_g, X[match(vnames,X$vertex),]$com))
}, .id = "inflation")
cat(" done\n")

# Number of different DAs in each MCL cluster
if (cfg$k_gc_da_sg != ""){
  msg("Loading already evaluated MCL results...")
  load(cfg$k_gc_da_sg)
  cat(" done\n")
}else{
  msg("Evaluating MCL communities...\n")
  k_gc_da_sg <- lapply(k_gc, function(X) {
    mclapply(X$da$com,
               evaluate_da_communities,
               da = X$da,
               mc.cores = cfg$dt_cores,
               #max.vector.size = 3e+09,
               mc.cleanup = TRUE,
               mc.silent = TRUE) %>%
      bind_rows() %>% inner_join(X$da, by = "com")
  })
  names(k_gc_da_sg) <- inflation_list
  msg("Evaluating MCL communities... done\n")

  k_gc_da_sg_file <- file.path(results, paste0("k_gc_da_sg_", time_string, ".Rda"))
  msg(paste0("Saving evaluated MCL communities results in ", k_gc_da_sg_file, "..."))
  save(k_gc_da_sg, file = k_gc_da_sg_file, compress = FALSE)
  cat(" done\n")
}

# Entropy
if (cfg$k_gc_entropy != ""){
  msg("Loading already calculated entropy results...")
  load(cfg$k_gc_entropy)
  cat(" done\n")
}else{
  msg("Calculating MCL communities entropy... (Patience, it takes a while)\n")
  k_hh_g_da <- k_hh_g %>%
    activate(nodes) %>%
    as_tibble() %>%
    separate_rows(archit, sep = "__") %>%
    as.data.table()

  k_gc_entropy <- lapply(k_gc_da_sg, function(X) {
    parallel::mclapply(X$com,
                          FUN = evaluate_component_entropy,
                          g = k_hh_g_da,
                          df = X,
                          mc.cores = cfg$entropy_cores,
                          #max.vector.size = 3e+09,
                          #ignore.interactive = TRUE,
                          mc.cleanup = TRUE,
                          mc.silent = TRUE) %>% bind_rows()
  })
  names(k_gc_entropy) <- inflation_list
  msg("Calculating MCL communities entropy... done\n")

  k_gc_entropy_file <- file.path(results, paste0("k_gc_entropy_", time_string, ".Rda"))
  msg(paste0("Saving evaluated MCL communities results in ", k_gc_entropy_file, "..."))
  save(k_gc_entropy, file = k_gc_entropy_file, compress = FALSE)
  cat(" done\n")
}

msg("Calculating MCL entropy summaries...")
k_gc_entropy_summary <- map_df(k_gc_entropy, function(X){
  X %>% mutate(e_0 = ifelse(entropy == 0, "eq", "gt")) %>%
    group_by(e_0) %>%
    count() %>%
    ungroup() %>%
    mutate(p = n/sum(n))
}, .id = "inflation")
cat(" done\n")

# Calculate inter and intra score modes
msg("Calculating inter/intra MCL communities score-per-column mode values...\n")

k_hh_g_dt <- k_hh_g %>%
  igraph::as_data_frame(what="edges") %>%
  rename(cl_name1 = from, cl_name2 = to, score_col = weight) %>%
  as.data.table()

k_orig_hh_mode <- parallel::mclapply(k_g_mcl_list, function(X){
  k_hh_g_dt %>%
    dt_left_join(X %>% dt_select(vertex, com) %>% dt_mutate(com = as.character(com), vertex = as.character(vertex)) %>% rename(cl_name1 = vertex) %>% unique() %>% as.data.table, by = "cl_name1") %>%
    dt_left_join(X %>% dt_select(vertex, com) %>% dt_mutate(com = as.character(com), vertex = as.character(vertex)) %>% rename(cl_name2 = vertex) %>% unique() %>% as.data.table, by = "cl_name2") %>%
    dt_filter(com.x == com.y) %>% mutate(com = com.x) %>% group_by(com) %>% summarise(mode = estimate_mode(score_col))
}, mc.cores = cfg$max_gc_jobs, #max.vector.size = 3e+09,
 mc.cleanup = TRUE, mc.silent = TRUE)
  #ignore.interactive = TRUE)

com_orig_intra_score <- map_df(k_orig_hh_mode, function(X){tibble(mode = estimate_mode(X$mode))}, .id = "inflation")
com_orig_inter_score <- map_df(k_gc, function(X){tibble(mode = estimate_mode(E(X$graph)$weight))}, .id = "inflation")

k_orig_hh_mode_summary <- com_orig_intra_score %>%
  mutate(class = "intra") %>%
  bind_rows(com_orig_inter_score %>% mutate(class = "inter"))
msg("Calculating inter/intra MCL communities score-per-column mode values... done\n\n")

# Get non-redundant set of DAs
# library(zoo)
# k_da_c <- map(p_doms %>% filter(class != "mono") %>% .$exp_rep %>% unique %>% strsplit("\\|"),
#               function(x) {
#                 g <- rollapply(data = x, 2, by=1, c) %>%
#                   graph_from_data_frame(directed = TRUE) %>% as_tbl_graph()
#               })


# Prefilter combinations that are too distant
# Get list of domains used during the MCL clustering
d <- p_doms %>% select(exp_rep) %>% unique() %>% as.data.table()

if (cfg$da_dist != ""){
  msg("Loading already calculated domain architectures distances...")
  lo_env$da_dist <- data.table::fread(input = cfg$da_dist, header = TRUE, showProgress = FALSE,
                                      nThread = cfg$dt_cores, sep = ",",  colClasses=c("character", "character", "numeric"))
  cat(" done\n")
}else{
  msg("Generating a non-redundant set of domain architectures...\n")
  n_comp <- (nrow(d)*(nrow(d) - 1))/2

  msg(paste("Calculating", scales::comma(n_comp),"pairwise domain architectures", "distances using", red("cosine"), "distance with a", red("q-gram of 3...")))
  lo_env$da_dist <- stringdist::stringdistmatrix(d$exp_rep, useNames = TRUE, method = "cosine", q = 3, nthread = cfg$dt_cores) %>%
    melt(as.matrix(da_dist), varnames = c("item1", "item2")) %>% rename(distance=value) %>% filter(distance!=0)
  cat(" done\n")

  da_dist_file <- file.path(results, paste0("da_dist_", time_string, ".tsv"))
  msg(paste0("Saving pairwise domain architectures distances in ", da_dist_file, "..."))
  fwrite(lo_env$da_dist, file = da_dist_file, col.names = TRUE, showProgress = FALSE)
  cat(" done\n")
}

msg("Filtering pairwise distances < 0.9...")
#da_dist_0.9 <- lo_env$da_dist %>% dt_filter(distance < 0.9) %>% dt_arrange(-distance)
da_dist_0.9 <- lo_env$da_dist %>% dt_filter(distance < 0.9 | is.na(distance)) %>% dt_arrange(-distance)
d[, n := seq_len(nrow(d))]
d.2 <- da_dist_0.9 %>%
  as_tibble() %>%
  mutate(idx1 = plyr::mapvalues(item1, from = d$exp_rep, to = d$n, warn_missing = FALSE),
         idx2 = plyr::mapvalues(item2, from = d$exp_rep, to = d$n, warn_missing = FALSE)) %>%
  as.data.table()

cat(" done\n")

#iterations <- 1000
# pb <- txtProgressBar(max = nproc, style = 3)
# progress <- function(n) setTxtProgressBar(pb, n)
# opts <- list(progress = progress)
# d.2 <- d.2 %>% as.data.table()

# pb <- txtProgressBar(min = 0, max = nrow(d.2), style = 3)
# k_das_refinement <- d.2[, c("n1", "n2", "V1_t", "V2_t", "DET") := {setTxtProgressBar(pb, .GRP);
#   n1 = str_count(item1, pattern = "\\|");
#   n2 = str_count(item2, pattern = "\\|");
#   V1_t = ifelse(str_count(item1, pattern = "\\|") <= str_count(item2, pattern = "\\|"), as.character(item1), as.character(item2));
#   V2_t = ifelse(str_count(item1, pattern = "\\|") > str_count(item2, pattern = "\\|"), as.character(item1), as.character(item2));
#   DET = ifelse(grepl(V1_t, V2_t, fixed = TRUE), TRUE, FALSE); list(n1, n2, V1_t, V2_t, DET)}, seq_len(nrow(d.2))]
# close(pb)


if (cfg$k_das_refinement != ""){
  msg("Loading nested domain architectures results...")
  load(cfg$k_das_refinement)
  cat(" done\n")
}else{
  # msg("Starting SNOW cluster...")
  # cl <- makeCluster(cfg$das_refinement_cores)
  # registerDoSNOW(cl)
  # cat( " done\n")

  msg("Starting nested domain architectures identification\n")
  # pb <- txtProgressBar(max = cfg$das_refinement_cores, style = 3)
  # progress <- function(n) setTxtProgressBar(pb, n)
  # opts <- list(progress = progress)
  #
  # k_das_refinement <- foreach(i=isplitRows(d.2, chunks = cfg$das_refinement_cores), .combine = rbindlist, .packages = c("tidyverse", "igraph", "tidygraph", "purrr", "stringr"), .options.snow = opts) %dopar% {
  #   purrr::pmap_dfr(i %>% as_tibble(),
  #                   ~tibble(idx1 = ..4,
  #                           idx2 = ..5,
  #                           item1 = ..1,
  #                           item2 = ..2,
  #                           n1 = str_count(item1, pattern = "\\|"),
  #                           n2 = str_count(item2, pattern = "\\|"),
  #                           V1_t = ifelse(str_count(item1, pattern = "\\|") <= str_count(item2, pattern = "\\|"), as.character(item1), as.character(item2)),
  #                           V2_t = ifelse(str_count(item1, pattern = "\\|") > str_count(item2, pattern = "\\|"), as.character(item1), as.character(item2)),
  #                           DET = ifelse(grepl(V1_t, V2_t, fixed = TRUE), TRUE, FALSE)))
  # }
  #
  # close(pb)
  # stopCluster(cl)

  pb <- txtProgressBar(max = nrow(d.2), style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  k_das_refinement <- d.2[, c("n1", "n2", "V1_t", "V2_t", "DET") := {setTxtProgressBar(pb, .GRP);
    n1 = str_count(item1, pattern = "\\|");
    n2 = str_count(item2, pattern = "\\|");
    V1_t = ifelse(str_count(item1, pattern = "\\|") <= str_count(item2, pattern = "\\|"), as.character(item1), as.character(item2));
    V2_t = ifelse(str_count(item1, pattern = "\\|") > str_count(item2, pattern = "\\|"), as.character(item1), as.character(item2));
    DET = ifelse(grepl(V1_t, V2_t, fixed = TRUE), TRUE, FALSE); list(n1, n2, V1_t, V2_t, DET)}, seq_len(nrow(d.2))]
  close(pb)

  msg("Starting nested domain architectures identification... done\n")
  k_das_refinement_file <- file.path(results, paste0("k_das_refinement_", time_string, ".Rda"))
  msg(paste0("Saving nested domain architectures results in ", k_das_refinement_file, "..."))
  save(k_das_refinement, file = k_das_refinement_file, compress = FALSE)
  cat(" done\n")
}

msg("Refining domain architectures...")
# Domains that are part of a larger domain
sub_da <- k_das_refinement  %>% dt_filter(DET == TRUE) %>% dt_select(V1_t) %>% unique()

# Large domains
large_da <- k_das_refinement %>% dt_filter(!(V2_t %in% sub_da$V1_t)) %>% dt_select(V2_t) %>% unique()

# get completeness
completness <- p_doms %>%
  group_by(archit) %>%
  summarise(complete = sum(complete),
            partial = sum(partial),
            N = complete + partial) %>%
  mutate(complete = complete/N,
         partial = partial/N)

# DAs that are part of a larger DA but at least they are seen in more than 75% of
# complete ORFs
sub_da <- k_das_refinement %>%
  dt_filter(DET == TRUE) %>%
  dt_select(V1_t) %>% dplyr::rename(archit = V1_t) %>% unique() %>%
  inner_join(completness, by = "archit") %>% arrange(desc(complete)) %>% filter(complete < 0.75) %>% as_tibble()

# create list with minimun set of DAs
# all DAs
final_da <- p_doms %>%
  select(exp_rep) %>%
  unique %>%
  filter(!(exp_rep %in% sub_da$archit))

cat(paste0(scales::comma(nrow(final_da)), " domain architectures found\n"))

# Summarize results -------------------------------------------------------
msg("Summarising MCL results...")
k_partition_stats <- map_df(k_g_mcl_list, function(X) {
  tibble(com_orig_n = length(unique(X$com)),
         com_orig_1mem = X %>% group_by(com) %>% count() %>% filter(n == 1) %>% nrow())
}, .id = "inflation") %>%
  inner_join(com_orig_intra_score %>% dplyr::rename(com_orig_intra_score = mode), by = "inflation")


k_partition_stats <- k_partition_stats %>%
  inner_join(map_df(k_gc_da_sg, function(X) {
    tibble(com_orig_1comp = X %>% filter(n_comp == 1) %>% nrow())
  }, .id = "inflation"), by = "inflation") %>%
  mutate(com_orig_1com_prop = com_orig_1comp/com_orig_n,
         com_orig_1mem_prop = com_orig_1mem/com_orig_n)


# Plot results
# ggpubr::ggarrange(k_orig_hh_mode_summary %>%
#                     select(-mode) %>%
#                     mutate(inflation = fct_rev(inflation)) %>%
#                     ggplot(aes(x=inter, xend=intra, y=inflation, group=inflation)) +
#                     geom_dumbbell(size=1, color="#e3e2e1",
#                                   colour_x = "#5b8124", colour_xend = "#bad744",
#                                   dot_guide=TRUE, dot_guide_size=0.25) +
#                     theme_light(),
#                   k_partition_stats %>%
#                     mutate(inflation = fct_rev(inflation)) %>%
#                     ggplot(aes(inflation, com_orig_n)) +
#                     geom_lollipop() +
#                     ggpubr::rotate() +
#                     theme_light(),
#                   k_partition_stats %>%
#                     mutate(inflation = fct_rev(inflation)) %>%
#                     ggplot(aes(inflation, com_orig_1mem_prop)) +
#                     geom_lollipop() +
#                     ggpubr::rotate() +
#                     theme_light(),
#                   k_partition_stats %>%
#                     mutate(inflation = fct_rev(inflation)) %>%
#                     ggplot(aes(inflation, com_orig_1com_prop)) +
#                     geom_lollipop() +
#                     ggpubr::rotate() +
#                     theme_light(),
#                   k_hh_gc_com_entropy_summary %>%
#                     filter(e_0 == "eq") %>%
#                     mutate(inflation = fct_rev(inflation)) %>%
#                     ggplot(aes(inflation, p)) +
#                     geom_lollipop() +
#                     ggpubr::rotate() +
#                     theme_light(),
#                   nrow = 1, ncol = 5, align = "hv"
# )


# Make a decision for the best inflation value ----------------------------

narchs <- nrow(final_da)
k_partition_stats_eval <- k_partition_stats %>%
  inner_join(k_gc_entropy_summary %>% filter(e_0 == "eq"), by = "inflation") %>%
  mutate(com_orig_1mem_prop = 1 - com_orig_1mem_prop,
         com_orig_n = 1/abs(com_orig_n - narchs),
         com_orig_n = com_orig_n/max(com_orig_n),
         com_orig_intra_score = com_orig_intra_score/max(com_orig_intra_score)) %>%
  select(inflation, com_orig_1com_prop, com_orig_1mem_prop, p, com_orig_intra_score, com_orig_n) %>%
  rename(grp = inflation) %>%
  gather(var, value, -grp) %>%
  group_by(grp) %>%
  mutate(nextval = lead(value, default = value[1]),
         angle = (1/5)*(2*pi),
         area = value*nextval*sin(angle)/2) %>%
  mutate(total = sum(area)) %>%
  ungroup()

cat(" done\n")

coord_radar <- function (theta = "x", start = 0, direction = 1)
{
  theta <- match.arg(theta, c("x", "y"))
  r <- if (theta == "x")
    "y"
  else "x"
  ggproto("CoordRadar", CoordPolar, theta = theta, r = r, start = start,
          direction = sign(direction),
          is_linear = function(coord) TRUE)
}


k_partition_stats_eval_plot <- k_partition_stats_eval %>%
  mutate(var1 = plyr::mapvalues(var, from = c("com_orig_1com_prop", "com_orig_1mem_prop", "p", "com_orig_intra_score", "com_orig_n"),
                                to = 1:5), var1 = as.numeric(var1)) %>% ungroup() %>%
  ggplot(aes(var1, value)) +
  #geom_polygon() +
  geom_point() +
  geom_polygon(, fill = "grey", size = 0.2, color = "black", alpha = 0.5) +
  geom_text(aes(0,0, label = round(total, 2)), color = "red") +
  facet_wrap(~grp ) +
  scale_y_continuous("", limits = c(0, 1), expand = c(0,0)) +
  scale_x_continuous("", breaks = 1:5, expand = c(0,0)) +
  theme_bw() +
  coord_radar()


k_partition_stats_eval_file <- file.path(results, paste0("k_partition_stats_eval_", time_string, ".tsv"))
msg("Saving MCL results stats...")
write_tsv(k_partition_stats_eval, k_partition_stats_eval_file, col_names = TRUE)
cat(" done\n")

k_partition_stats_eval_plot_file <- file.path(results, paste0("k_partition_stats_eval_plot_", time_string, ".pdf"))
k_partition_stats_eval_data_file <- file.path(results, paste0("k_partition_stats_eval_data_", time_string, ".tsv"))
msg("Saving MCL results radar plot...\n")
ggsave(filename = k_partition_stats_eval_plot_file, plot = k_partition_stats_eval_plot)
k_partition_stats_eval %>%
  mutate(var1 = plyr::mapvalues(var, from = c("com_orig_1com_prop", "com_orig_1mem_prop", "p", "com_orig_intra_score", "com_orig_n"),
                                to = 1:5), var1 = as.numeric(var1)) %>% ungroup() %>%
  write_tsv(k_partition_stats_eval_data_file, col_names = TRUE)
msg("Saving MCL results radar plot... done\n")


best_inflation <- k_partition_stats_eval %>%
  ungroup() %>%
  select(grp, total) %>%
  unique() %>%
  top_n(n = 1, wt = total) %>%
  .$grp

if (length(best_inflation >1 )) {
    best_inflation <- max(best_inflation)
}
msg(paste0("\nOptimal inflation value identified at ", red(best_inflation), "\n\n"))


# Add missing clusters to MCL communities ----------------------------------

# Add the missing nodes to existing communities
# Try to find the best hits for the missing nodes
# Based on coverage/probability/score per position


msg(paste0("Selecting MCL ", best_inflation, " value..."))
k_missing_dt <- lo_env$k_hhblits_missing %>% ungroup()
k_mcl_coms <- k_g_mcl_list[[best_inflation]]
k_max_com <- max(k_mcl_coms$com)
cat(" done\n")


k_mids_l <- k_missing_ids %>% length()
if (k_mids_l > 0){
  msg(paste0("Trying to assign ", scales::comma(k_mids_l), " missing clusters to existing communities...\n"))

  # missing_ids -> 7047

  # Try to identify these cluster that can have some homology to existing
  # We only take the queries for missing ids
  # We check that the target is in the MCL communities
  k_missing_d <- k_missing_dt %>%
    as.data.table() %>%
    dt_filter(cl_name1 %in% k_missing_ids) %>%
    dt_filter(cl_name2 %in% k_mcl_coms$vertex) %>%
    dt_left_join(k_mcl_coms %>% mutate(vertex = as.character(vertex)) %>% rename(cl_name2 = vertex), by = "cl_name2")

  msg("First pass: Using more relaxed HHBLITS filtering thresholds (prob > 50 and cov >= 40)...")
  # MCL clusters with more relaxed filters (p50 and cov >= 40)
  # we just keep the best hit
  k_missing_d_pass <- k_missing_d %>%
    dt_filter(probability >= 50, (q_cov >= 0.4 | t_cov >= 0.4)) %>%
    as_tibble() %>%
    group_by(cl_name1) %>%
    arrange(-probability, -q_cov, -t_cov) %>%
    #mutate(com1 = majority_vote(com)$majority,
    mutate(nhit = row_number()) %>%
    arrange(cl_name1, nhit) %>%
    ungroup() %>%
    #filter(probability >= 30, nhit <= 3) %>%
    filter(nhit == 1) %>%
    select(cl_name1, cl_name2, score_col, com) %>%
    rename(weight = score_col)
  cat(" done\n")
  # Find the ones we couldn't classify
  # 299
  k_to_assign <- setdiff(k_missing_ids, k_missing_d_pass$cl_name1 %>% unique())
  k_assigned <- setdiff(k_missing_ids, k_to_assign)

  # We need to check if any of the remaining no assigned clusters have any homology
  # to the just classified using more relaxed parameters

  msg("Second pass: Find secondary relationships between the newly assigned clusters and missing clusters...")
  k_missing_d_pass_1 <- bind_rows(lo_env$k_hhblits_missing %>%
                                    dt_inner_join(k_missing_d_pass %>% select(cl_name1, com), by = "cl_name1") %>%
                                    dt_filter(cl_name2 %in% k_to_assign) %>% dt_mutate(cl_name = cl_name2, com_f = com),
                                  lo_env$k_hhblits_missing %>%
                                    dt_inner_join(k_missing_d_pass %>% select(cl_name2, com), by = "cl_name2") %>%
                                    dt_filter(cl_name1 %in% k_to_assign) %>% dt_mutate(cl_name = cl_name1, com_f = com)) %>%
    dt_filter(probability > 50, (q_cov > 0.4 | t_cov > 0.4)) %>%
    group_by(cl_name) %>%
    arrange(-probability, -q_cov, -t_cov) %>%
    #mutate(com1 = majority_vote(com)$majority,
    mutate(nhit = row_number()) %>%
    arrange(cl_name, nhit) %>%
    ungroup() %>%
    #filter(probability >= 30, nhit <= 3) %>%
    filter(nhit == 1) %>%
    dt_select(cl_name, com_f) %>%
    dplyr::rename(com = com_f) %>%
    as_tibble()
  cat(" done\n")

  # The ones that cannot be assigned, we will try to find any existing connections between no classified
  # and run MCL with the best inflation value and create new clusters

  msg("Third pass: Non assigned clusters will create new MCL communities...")
  k_missing_d_pass_2 <- bind_rows(lo_env$k_hhblits_missing %>%
                                    dt_filter(cl_name1 %in% k_to_assign, cl_name2 %in% k_to_assign) %>%
                                    dt_filter(probability > 50, (q_cov > 0.4 | t_cov > 0.4)) %>%
                                    dt_inner_join(k_missing_d_pass_1 %>% select(cl_name, com) %>% rename(cl_name1 = cl_name), by = "cl_name1") %>%
                                    dt_mutate(cl_name = cl_name1, com_f = com),
                                  lo_env$k_hhblits_missing %>%
                                    dt_filter(cl_name1 %in% k_to_assign, cl_name2 %in% k_to_assign) %>%
                                    dt_filter(probability > 50, (q_cov > 0.4 | t_cov > 0.4)) %>%
                                    dt_inner_join(k_missing_d_pass_1 %>% select(cl_name, com) %>% rename(cl_name2 = cl_name), by = "cl_name2") %>%
                                    dt_mutate(cl_name = cl_name2, com_f = com)) %>%
    group_by(cl_name) %>%
    arrange(-probability, -q_cov, -t_cov) %>%
    #mutate(com1 = majority_vote(com)$majority,
    mutate(nhit = row_number()) %>%
    arrange(cl_name, nhit) %>%
    ungroup() %>%
    #filter(probability >= 30, nhit <= 3) %>%
    filter(nhit == 1) %>%
    dt_select(cl_name, com_f) %>%
    dplyr::rename(com = com_f) %>%
    as_tibble()
  cat(" done\n")

  # msg(paste0("We have been able to assign ",
  #            green(scales::comma(length(assigned))),
  #            " missing clusters. Still missing ",
  #            red(scales::comma(length(to_assign))),
  #            " clusters\n"
  # ))

  msg("Collecting and aggregating all communities...")
  k_missing_ids_1 <- setdiff(k_missing_ids, c(k_missing_d_pass$cl_name1 %>% unique, k_missing_d_pass_1$cl_name, k_missing_d_pass_2$cl_name))

  k_missing_d_pass_3 <- tibble(cl_name = k_missing_ids_1) %>% mutate(com = k_max_com + row_number())

  k_communities <-  bind_rows(k_missing_d_pass %>% select(cl_name1, com) %>% rename(cl_name = cl_name1),
                             k_missing_d_pass_1,
                             k_missing_d_pass_2,
                             k_missing_d_pass_3 %>% mutate(cl_name = as.character(cl_name)),
                             k_mcl_coms %>% rename(cl_name = vertex) %>% mutate(cl_name = as.character(cl_name)) %>% select(cl_name, com)) %>% unique()
  k_n_comp <- k_communities$com %>% unique() %>% length()
  k_n_clus <-  k_communities$cl_name %>% unique() %>% length()
  cat(" done\n")
}else{
  msg(paste0("Trying to assign ",
             scales::comma(k_mids_l),
             " missing clusters to existing communities...",
             red("Skipped"), "\n"))
  msg("All clusters already assigned to a component\n")
  k_communities <- k_mcl_coms %>% rename(cl_name = vertex) %>% mutate(cl_name = as.character(cl_name)) %>% select(cl_name, com) %>% unique()
  k_n_comp <- k_communities$com %>% unique() %>% length()
  k_n_clus <-  k_communities$cl_name %>% unique() %>% length()
}

k_correct <- length(cl_cat %>% dt_filter(category == "K") %>% .$cl_name) == k_n_clus
if(k_correct){
  msg(paste0("All clusters assigned: ", green(k_correct), "\n"))
}else{
  msg(paste0("All clusters assigned: ", red(k_correct), "\n"))
  #msg("Something went wrong... quitting\n!")
  #quit()
  k_max_com <- max(k_communities$com)
  missing_final <- cl_cat %>% dt_filter(category == "K") %>% mutate(cl_name=as.character(cl_name)) %>%
   filter(!cl_name %in% k_communities$cl_name) %>% mutate(com=k_max_com + row_number())
  k_communities <- bind_rows(k_communities,missing_final)
  k_n_comp <- k_communities$com %>% unique() %>% length()
  k_n_clus <-  k_communities$cl_name %>% unique() %>% length()
}


msg(paste0("We have been obtained ", scales::comma(k_n_comp), " from ", scales::comma(k_n_clus), " clusters\n"))

# Do we have all communities
#nrow(k_comp) == nrow(k_communities)

k_comp_file <- file.path(results, paste0("k_communities_", time_string, ".tsv"))
msg(paste0("Saving communities results to ", k_comp_file, "..."))
readr::write_tsv(k_communities, path = k_comp_file, col_names = FALSE)
cat(" done\n\n")

msg(paste0("KNOWN component inference ", green("DONE"), "\n\n"))


###########################################################################
# Get communities for the KWP ----------------------------------------------
###########################################################################
msg(paste0("Starting KWP component inference\n\n"))

msg("Reading KWP HHBLITS all-vs-all results...")
lo_env$kwp_hhblits_all <- fread(input = cfg$kwp_hhblits_results,
                                header = FALSE, verbose = cfg$dt_verbose, nThread = cfg$dt_cores)
names(lo_env$kwp_hhblits_all) <- c("cl_name1", "cl_name2", "probability", "e-value", "Score",
                                   "Cols", "q_start", "q_stop", "t_start", "t_stop", "q_len", "t_len", "q_cov",
                                   "t_cov")
cat(paste(" Read", green(scales::comma(nrow(lo_env$kwp_hhblits_all))), "entries\n"))

msg("Filtering KWP HHBLITS all-vs-all results with P>=50 and cov > 0.6...")
lo_env$kwp_hhblits <- lo_env$kwp_hhblits_all %>% dt_filter(probability > 50, q_cov >= 0.6, t_cov >= 0.6)
lo_env$kwp_hhblits <- lo_env$kwp_hhblits %>% dt_filter(cl_name1 != cl_name2)

# How many clusters do we have (removed self-hits)
# 522,765 p50;c0.6
kwp_unique_hhblits_cl <- c(lo_env$kwp_hhblits$cl_name1, lo_env$kwp_hhblits$cl_name2) %>% unique()
cat(paste0(" Found ", scales::comma(length(kwp_unique_hhblits_cl)), " clusters\n"))


# Missing clusters
# Num: 27,244
msg("Identifying missing clusters...")
kwp_missing_ids <- setdiff(cl_cat %>% dt_filter(category == "KWP") %>% .$cl_name, kwp_unique_hhblits_cl)
lo_env$kwp_hhblits_missing <- lo_env$kwp_hhblits_all %>%
  dt_mutate(cl_name1 = as.character(cl_name1), cl_name2 = as.character(cl_name2)) %>%
  dt_filter(cl_name1 %in% kwp_missing_ids | cl_name2 %in% kwp_missing_ids) %>%
  as_tibble() %>% mutate(score_col = Score/Cols)
cat(paste0("Found ", scales::comma(length(kwp_missing_ids))), "clusters missing\n")

if(dim(lo_env$kwp_hhblits)[1]!= 0){
  # Let's create a graph
  msg("Inferring graph from HHBLITS results...")
  lo_env$kwp_hhb_bh <- lo_env$kwp_hhblits %>%
    as_tibble() %>%
    mutate(cl_name1 = as.character(cl_name1),
         cl_name2 = as.character(cl_name2)) %>%
         mutate(score_col = Score/Cols)

  kwp_hh_g <- lo_env$kwp_hhb_bh %>%
    select(cl_name1, cl_name2, score_col) %>%
    rename(weight = score_col) %>%
  igraph::graph_from_data_frame(directed = FALSE) %>%
  igraph::simplify(remove.multiple = TRUE, remove.loops = TRUE, edge.attr.comb = list("max")) %>%
  as_tbl_graph()

cat(paste0(" Graph contains ", scales::comma(vcount(kwp_hh_g)), " vertices and ", scales::comma(ecount(kwp_hh_g)), " edges\n"))

if (cfg$kwp_g_mcl_list != ""){
  msg("Loading already computed MCL clusters...")
  load(cfg$kwp_g_mcl_list)
  inflation_list <- seq(cfg$mcl_inflation_min, cfg$mcl_inflation_max, cfg$mcl_inflation_step)
  cat(" done\n")
}else{
  msg("Preparing graph for MCL clustering...")
  gx <- kwp_hh_g

  w <- E(gx)$weight
  w <- w - min(w) + 0.001

  E(gx)$weight <- w

  msg(paste0("Running MCL clustering with inflation values from ", cfg$mcl_inflation_min, " to ", cfg$mcl_inflation_max, " with ", cfg$mcl_inflation_step, " steps...\n"))
  inflation_list <- seq(cfg$mcl_inflation_min, cfg$mcl_inflation_max, cfg$mcl_inflation_step)
  mcl_bin <- cfg$mcl_bin
  kwp_g_mcl_list <- parallel::mclapply(inflation_list,
                                          optimal_mcl,
                                          G = kwp_hh_g,
                                          mcl_cores = cfg$mcl_cores,
                                          #ignore.interactive = TRUE,
                                          Gx = gx,
                                          #max.vector.size = 1e+07,
                                          mc.cores = cfg$mcl_jobs,
                                          mc.cleanup = TRUE,
                                          mc.silent = TRUE,
                                          tmp_dir = tmp)

  #k_g_mcl_list <- lapply(inflation_list, optimal_mcl, G = k_hh_g, Gx = gx, mcl_cores = cfg$mcl_cores, tmp_dir = tmp)
  names(kwp_g_mcl_list) <- inflation_list

  kwp_g_mcl_list_file <- file.path(results, paste0("kwp_g_mcl_list_", time_string, ".Rda"))
  msg(paste0("Saving computed MCL communities in ", kwp_g_mcl_list_file, "..."))
  save(kwp_g_mcl_list, file = kwp_g_mcl_list_file, compress = FALSE)
  cat(" done\n")


}

# We run mcl with different inflation parameters:
if (cfg$kwp_gc != ""){
  msg("Loading already contracted MCL communities...")
  load(cfg$kwp_gc)
  cat(" done\n")
}else{
  msg("Contracting MCL communities...\n")
  kwp_gc <- parallel::mclapply(kwp_g_mcl_list,
                                  contract_graphs,
                                  G = kwp_hh_g,
                                  #max.vector.size = 3e+09,
                                  mc.cores = cfg$max_gc_jobs,
                                  mc.cleanup = TRUE,
                                  mc.silent = TRUE)
  names(kwp_gc) <- inflation_list
  msg("Contracting MCL communities... done\n")

  kwp_gc_file <- file.path(results, paste0("kwp_gc_", time_string, ".Rda"))
  msg(paste0("Saving contracted MCL communities results in ", kwp_gc_file, "..."))
  save(kwp_gc, file = kwp_gc_file, compress = FALSE)
  cat(" done\n")
}

msg("Computing modularity for each MCL community result...")
kwp_hh_g_modularity <- map_df(kwp_g_mcl_list, function(X){
  vnames <- V(kwp_hh_g)$name
  tibble(modularity = modularity(kwp_hh_g, X[match(vnames,X$vertex),]$com))
}, .id = "inflation")
cat(" done\n")



# Inter and intra score modes
# Calculate nter and intra score modes
msg("Calculating inter/intra MCL communities score-per-column mode values...\n")

kwp_hh_g_dt <- kwp_hh_g %>%
  igraph::as_data_frame(what="edges") %>%
  rename(cl_name1 = from, cl_name2 = to, score_col = weight) %>%
  as.data.table()

kwp_orig_hh_mode <- parallel::mclapply(kwp_g_mcl_list, function(X){
  kwp_hh_g_dt %>%
    dt_left_join(X %>% dt_select(vertex, com) %>% dt_mutate(com = as.character(com), vertex = as.character(vertex)) %>% rename(cl_name1 = vertex) %>% unique() %>% as.data.table, by = "cl_name1") %>%
    dt_left_join(X %>% dt_select(vertex, com) %>% dt_mutate(com = as.character(com), vertex = as.character(vertex)) %>% rename(cl_name2 = vertex) %>% unique() %>% as.data.table, by = "cl_name2") %>%
    dt_filter(com.x == com.y) %>% mutate(com = com.x) %>% group_by(com) %>% summarise(mode = estimate_mode(score_col))
}, mc.cores = cfg$max_gc_jobs,
  #max.vector.size = 3e+09,
  mc.cleanup = TRUE, mc.silent = TRUE)

com_orig_intra_score <- map_df(kwp_orig_hh_mode, function(X){tibble(mode = estimate_mode(X$mode))}, .id = "inflation")
com_orig_inter_score <- map_df(kwp_gc, function(X){tibble(mode = estimate_mode(E(X$graph)$weight))}, .id = "inflation")

kwp_orig_hh_mode_summary <- com_orig_intra_score %>%
  mutate(class = "intra") %>%
  bind_rows(com_orig_inter_score %>% mutate(class = "inter"))
msg("Calculating inter/intra MCL communities score-per-column mode values... done\n")

msg("Summarising MCL results...")
kwp_partition_stats <- map_df(kwp_g_mcl_list, function(X) {
  tibble(com_orig_n = length(unique(X$com)),
         com_orig_1mem = X %>% group_by(com) %>% count() %>% filter(n == 1) %>% nrow())
}, .id = "inflation") %>%
  inner_join(com_orig_intra_score %>% dplyr::rename(com_orig_intra_score = mode), by = "inflation")
cat(" done\n")


kwp_partition_stats_file <- file.path(results, paste0("kwp_partition_stats_", time_string, ".tsv"))
msg("Saving MCL results stats...")
write_tsv(kwp_partition_stats, kwp_partition_stats_file, col_names = TRUE)
cat(" done\n\n")


# Add missing clusters to MCL communities ----------------------------------

# Add the missing nodes to existing communities
# Try to find the best hits for the missing nodes
# Based on coverage/probability/score per position


msg(paste0("Selecting MCL ", best_inflation, " value..."))
kwp_mcl_coms <- kwp_g_mcl_list[[best_inflation]]
kwp_max_com <- max(kwp_mcl_coms$com)
cat(" done\n")


kwp_mids_l <- kwp_missing_ids %>% length()

if (kwp_mids_l > 0){
  msg(paste0("Trying to assign ", scales::comma(kwp_mids_l), " missing clusters to existing communities...\n"))

  kwp_missing_dt <- lo_env$kwp_hhblits_missing %>% ungroup()

  # missing_ids -> 7047

  # Try to identify these cluster that can have some homology to existing
  # We only take the queries for missing ids
  # We check that the target is in the MCL communities
  kwp_missing_d <- kwp_missing_dt %>%
    as.data.table() %>%
    dt_filter(cl_name1 %in% kwp_missing_ids) %>%
    dt_filter(cl_name2 %in% kwp_mcl_coms$vertex) %>%
    dt_left_join(kwp_mcl_coms %>% mutate(vertex = as.character(vertex)) %>% rename(cl_name2 = vertex), by = "cl_name2")

  msg("First pass: Using more relaxed HHBLITS filtering thresholds (prob > 50 and cov >= 40)...")
  # MCL clusters with more relaxed filters (p50 and cov >= 40)
  # we just keep the best hit
  kwp_missing_d_pass <- kwp_missing_d %>%
    dt_filter(probability >= 50, (q_cov >= 0.4 | t_cov >= 0.4)) %>%
    as_tibble() %>%
    group_by(cl_name1) %>%
    arrange(-probability, -q_cov, -t_cov) %>%
    #mutate(com1 = majority_vote(com)$majority,
    mutate(nhit = row_number()) %>%
    arrange(cl_name1, nhit) %>%
    ungroup() %>%
    #filter(probability >= 30, nhit <= 3) %>%
    filter(nhit == 1) %>%
    select(cl_name1, cl_name2, score_col, com) %>%
    rename(weight = score_col)
  cat(" done\n")
  # Find the ones we couldn't classify
  # 299
  kwp_to_assign <- setdiff(kwp_missing_ids, kwp_missing_d_pass$cl_name1 %>% unique())
  kwp_assigned <- setdiff(kwp_missing_ids, kwp_to_assign)

  # We need to check if any of the remaining no assigned clusters have any homology
  # to the just classified using more relaxed parameters

  msg("Second pass: Find secondary relationships between the newly assigned clusters and missing clusters...")
  kwp_missing_d_pass_1 <- bind_rows(lo_env$kwp_hhblits_missing %>%
                                      dt_inner_join(kwp_missing_d_pass %>% select(cl_name1, com), by = "cl_name1") %>%
                                      dt_filter(cl_name2 %in% kwp_to_assign) %>% dt_mutate(cl_name = cl_name2, com_f = com),
                                    lo_env$kwp_hhblits_missing %>%
                                      dt_inner_join(kwp_missing_d_pass %>% select(cl_name2, com), by = "cl_name2") %>%
                                      dt_filter(cl_name1 %in% kwp_to_assign) %>% dt_mutate(cl_name = cl_name1, com_f = com)) %>%
    dt_filter(probability > 50, (q_cov > 0.4 | t_cov > 0.4)) %>%
    group_by(cl_name) %>%
    arrange(-probability, -q_cov, -t_cov) %>%
    #mutate(com1 = majority_vote(com)$majority,
    mutate(nhit = row_number()) %>%
    arrange(cl_name, nhit) %>%
    ungroup() %>%
    #filter(probability >= 30, nhit <= 3) %>%
    filter(nhit == 1) %>%
    dt_select(cl_name, com_f) %>%
    dplyr::rename(com = com_f) %>%
    as_tibble()
  cat(" done\n")

  # The ones that cannot be assigned, we will try to find any existing connections between no classified
  # and run MCL with the best inflation value and create new clusters

  msg("Third pass: Non assigned clusters will create new MCL communities...")
  kwp_missing_d_pass_2 <- bind_rows(lo_env$kwp_hhblits_missing %>%
                                      dt_filter(cl_name1 %in% kwp_to_assign, cl_name2 %in% kwp_to_assign) %>%
                                      dt_filter(probability > 50, (q_cov > 0.4 | t_cov > 0.4)) %>%
                                      dt_inner_join(kwp_missing_d_pass_1 %>% select(cl_name, com) %>% rename(cl_name1 = cl_name), by = "cl_name1") %>%
                                      dt_mutate(cl_name = cl_name1, com_f = com),
                                    lo_env$kwp_hhblits_missing %>%
                                      dt_filter(cl_name1 %in% kwp_to_assign, cl_name2 %in% kwp_to_assign) %>%
                                      dt_filter(probability > 50, (q_cov > 0.4 | t_cov > 0.4)) %>%
                                      dt_inner_join(kwp_missing_d_pass_1 %>% select(cl_name, com) %>% rename(cl_name2 = cl_name), by = "cl_name2") %>%
                                      dt_mutate(cl_name = cl_name2, com_f = com)) %>%
    group_by(cl_name) %>%
    arrange(-probability, -q_cov, -t_cov) %>%
    #mutate(com1 = majority_vote(com)$majority,
    mutate(nhit = row_number()) %>%
    arrange(cl_name, nhit) %>%
    ungroup() %>%
    #filter(probability >= 30, nhit <= 3) %>%
    filter(nhit == 1) %>%
    dt_select(cl_name, com_f) %>%
    dplyr::rename(com = com_f) %>%
    as_tibble()
  cat(" done\n")

  msg("Collecting and aggregating all communities...")
  kwp_missing_ids_1 <- setdiff(kwp_missing_ids, c(kwp_missing_d_pass$cl_name1 %>% unique, kwp_missing_d_pass_1$cl_name, kwp_missing_d_pass_2$cl_name))

  kwp_missing_d_pass_3 <- tibble(cl_name = kwp_missing_ids_1) %>% mutate(com = kwp_max_com + row_number())

  kwp_communities <-  bind_rows(kwp_missing_d_pass %>% select(cl_name1, com) %>% rename(cl_name = cl_name1),
                               kwp_missing_d_pass_1,
                               kwp_missing_d_pass_2,
                               kwp_missing_d_pass_3 %>% mutate(cl_name = as.character(cl_name)),
                               kwp_mcl_coms %>% rename(cl_name = vertex) %>% mutate(cl_name = as.character(cl_name)) %>% select(cl_name, com)) %>% unique()
  kwp_n_comp <- kwp_communities$com %>% unique() %>% length()
  kwp_n_clus <-  kwp_communities$cl_name %>% unique() %>% length()
  cat(" done\n")

}else{
  msg(paste0("Trying to assign ",
             scales::comma(kwp_mids_l),
             " missing clusters to existing communities...",
             red("Skipped"), "\n"))
  msg("All clusters already assigned to a component\n")
  kwp_communities <- kwp_mcl_coms %>% rename(cl_name = vertex) %>% mutate(cl_name = as.character(cl_name)) %>% select(cl_name, com) %>% unique()
  kwp_n_comp <- kwp_communities$com %>% unique() %>% length()
  kwp_n_clus <-  kwp_communities$cl_name %>% unique() %>% length()
}
} else {
  kwp_communities <- kwp_missing_ids %>% tibble::enframe(name = NULL) %>%
   rename(cl_name=value) %>% mutate(com=row_number(), cl_name=as.character(cl_name))
  kwp_n_comp <- kwp_communities$com %>% unique() %>% length()
  kwp_n_clus <-  kwp_communities$cl_name %>% unique() %>% length()
}
msg(paste0("We have been obtained ", scales::comma(kwp_n_comp), " from ", scales::comma(kwp_n_clus), " clusters\n"))

# Do we have all communities
kwp_correct <- length(cl_cat %>% dt_filter(category == "KWP") %>% .$cl_name) == kwp_n_clus
if(kwp_correct){
  msg(paste0("All clusters assigned: ", green(kwp_correct), "\n"))
}else{
  msg(paste0("All clusters assigned: ", red(kwp_correct), "\n"))
  msg("Something went wrong... quitting\n!")
  quit()
}

kwp_comp_file <- file.path(results, paste0("kwp_communities_", time_string, ".tsv"))
msg(paste0("Saving communities results to ", kwp_comp_file, "..."))
readr::write_tsv(kwp_communities, path = kwp_comp_file, col_names = FALSE)
cat(" done\n\n")

msg(paste0("KWP component inference ", green("DONE"), "\n\n"))

###########################################################################
# Get communities for the GU ----------------------------------------------
###########################################################################
msg(paste0("Starting GU component inference\n\n"))

msg("Reading GU HHBLITS all-vs-all results...")
lo_env$gu_hhblits_all <- fread(input = cfg$gu_hhblits_results,
                               header = FALSE, verbose = cfg$dt_verbose, nThread = cfg$dt_cores)
names(lo_env$gu_hhblits_all) <- c("cl_name1", "cl_name2", "probability", "e-value", "Score",
                                  "Cols", "q_start", "q_stop", "t_start", "t_stop", "q_len", "t_len", "q_cov",
                                  "t_cov")
cat(paste(" Read", green(scales::comma(nrow(lo_env$gu_hhblits_all))), "entries\n"))

msg("Filtering GU HHBLITS all-vs-all results with P>=50 and cov > 0.6...")
lo_env$gu_hhblits <- lo_env$gu_hhblits_all %>% dt_filter(probability > 50, q_cov >= 0.6, t_cov >= 0.6)
lo_env$gu_hhblits <- lo_env$gu_hhblits %>% dt_filter(cl_name1 != cl_name2)

# How many clusters do we have (removed self-hits)
# 522,765 p50;c0.6
gu_unique_hhblits_cl <- c(lo_env$gu_hhblits$cl_name1, lo_env$gu_hhblits$cl_name2) %>% unique()
cat(paste0(" Found ", scales::comma(length(gu_unique_hhblits_cl)), " clusters\n"))


# Missing clusters
# Num: 27,244
msg("Identifying missing clusters...")
gu_missing_ids <- setdiff(cl_cat %>% dt_filter(category == "GU") %>% .$cl_name, gu_unique_hhblits_cl)
lo_env$gu_hhblits_missing <- lo_env$gu_hhblits_all %>%
  dt_mutate(cl_name1 = as.character(cl_name1), cl_name2 = as.character(cl_name2)) %>%
  dt_filter(cl_name1 %in% gu_missing_ids | cl_name2 %in% gu_missing_ids) %>%
  as_tibble() %>% mutate(score_col = Score/Cols)
cat(paste0("Found ", scales::comma(length(gu_missing_ids))), "clusters missing\n")

if(dim(lo_env$gu_hhblits)[1]!= 0){
  # Let's create a graph
  msg("Inferring graph from HHBLITS results...")
  lo_env$gu_hhb_bh <- lo_env$gu_hhblits %>%
    as_tibble() %>%
    mutate(cl_name1 = as.character(cl_name1),
         cl_name2 = as.character(cl_name2)) %>%
    mutate(score_col = Score/Cols)

  gu_hh_g <- lo_env$gu_hhb_bh %>%
  select(cl_name1, cl_name2, score_col) %>%
  rename(weight = score_col) %>%
  igraph::graph_from_data_frame(directed = FALSE) %>%
  igraph::simplify(remove.multiple = TRUE, remove.loops = TRUE, edge.attr.comb = list("max")) %>%
  as_tbl_graph()

cat(paste0(" Graph contains ", scales::comma(vcount(gu_hh_g)), " vertices and ", scales::comma(ecount(gu_hh_g)), " edges\n"))

if (cfg$gu_g_mcl_list != ""){
  msg("Loading already computed MCL clusters...")
  load(cfg$gu_g_mcl_list)
  inflation_list <- seq(cfg$mcl_inflation_min, cfg$mcl_inflation_max, cfg$mcl_inflation_step)
  cat(" done\n")
}else{
  msg("Preparing graph for MCL clustering...")
  gx <- gu_hh_g

  w <- E(gx)$weight
  w <- w - min(w) + 0.001

  E(gx)$weight <- w

  cat(" done\n")

  msg(paste0("Running MCL clustering with inflation values from ", cfg$mcl_inflation_min, " to ", cfg$mcl_inflation_max, " with ", cfg$mcl_inflation_step, " steps...\n"))
  inflation_list <- seq(cfg$mcl_inflation_min, cfg$mcl_inflation_max, cfg$mcl_inflation_step)
  mcl_bin <- cfg$mcl_bin
  gu_g_mcl_list <- parallel::mclapply(inflation_list,
                                         optimal_mcl,
                                         G = gu_hh_g,
                                         mcl_cores = cfg$mcl_cores,
                                         #ignore.interactive = TRUE,
                                         Gx = gx,
                                         #max.vector.size = 1e+07,
                                         mc.cores = cfg$mcl_jobs,
                                         mc.cleanup = TRUE,
                                         mc.silent = TRUE,
                                         tmp_dir = tmp)

  #k_g_mcl_list <- lapply(inflation_list, optimal_mcl, G = k_hh_g, Gx = gx, mcl_cores = cfg$mcl_cores, tmp_dir = tmp)
  names(gu_g_mcl_list) <- inflation_list

  gu_g_mcl_list_file <- file.path(results, paste0("gu_g_mcl_list_", time_string, ".Rda"))
  msg(paste0("Saving computed MCL communities in ", gu_g_mcl_list_file, "..."))
  save(gu_g_mcl_list, file = gu_g_mcl_list_file, compress = FALSE)
  cat(" done\n")


}

# We run mcl with different inflation parameters:
if (cfg$gu_gc != ""){
  msg("Loading already contracted MCL communities...")
  load(cfg$gu_gc)
  cat(" done\n")
}else{
  msg("Contracting MCL communities...\n")
  gu_gc <- parallel::mclapply(gu_g_mcl_list,
                                 contract_graphs,
                                 G = gu_hh_g,
                                 #max.vector.size = 3e+09,
                                 mc.cores = cfg$max_gc_jobs,
                                 mc.cleanup = TRUE,
                                 mc.silent = TRUE)
  names(gu_gc) <- inflation_list
  msg("Contracting MCL communities... done\n")

  gu_gc_file <- file.path(results, paste0("gu_gc_", time_string, ".Rda"))
  msg(paste0("Saving contracted MCL communities results in ", gu_gc_file, "..."))
  save(gu_gc, file = gu_gc_file, compress = FALSE)
  cat(" done\n")
}

msg("Computing modularity for each MCL community result...")
gu_hh_g_modularity <- map_df(gu_g_mcl_list, function(X){
  vnames <- V(gu_hh_g)$name
  tibble(modularity = modularity(gu_hh_g, X[match(vnames,X$vertex),]$com))
}, .id = "inflation")
cat(" done\n")



# Inter and intra score modes
# Calculate nter and intra score modes
msg("Calculating inter/intra MCL communities score-per-column mode values...\n")

gu_hh_g_dt <- gu_hh_g %>%
  igraph::as_data_frame(what="edges") %>%
  rename(cl_name1 = from, cl_name2 = to, score_col = weight) %>%
  as.data.table()

gu_orig_hh_mode <- parallel::mclapply(gu_g_mcl_list, function(X){
  gu_hh_g_dt %>%
    dt_left_join(X %>% dt_select(vertex, com) %>% dt_mutate(com = as.character(com), vertex = as.character(vertex)) %>% rename(cl_name1 = vertex) %>% unique() %>% as.data.table, by = "cl_name1") %>%
    dt_left_join(X %>% dt_select(vertex, com) %>% dt_mutate(com = as.character(com), vertex = as.character(vertex)) %>% rename(cl_name2 = vertex) %>% unique() %>% as.data.table, by = "cl_name2") %>%
    dt_filter(com.x == com.y) %>% mutate(com = com.x) %>% group_by(com) %>% summarise(mode = estimate_mode(score_col))
}, mc.cores = cfg$max_gc_jobs, #max.vector.size = 3e+09,
 mc.cleanup = TRUE, mc.silent = TRUE)

com_orig_intra_score <- map_df(gu_orig_hh_mode, function(X){tibble(mode = estimate_mode(X$mode))}, .id = "inflation")
com_orig_inter_score <- map_df(gu_gc, function(X){tibble(mode = estimate_mode(E(X$graph)$weight))}, .id = "inflation")

gu_orig_hh_mode_summary <- com_orig_intra_score %>%
  mutate(class = "intra") %>%
  bind_rows(com_orig_inter_score %>% mutate(class = "inter"))
msg("Calculating inter/intra MCL communities score-per-column mode values... done\n\n")

msg("Summarising MCL results...")
gu_partition_stats <- map_df(gu_g_mcl_list, function(X) {
  tibble(com_orig_n = length(unique(X$com)),
         com_orig_1mem = X %>% group_by(com) %>% count() %>% filter(n == 1) %>% nrow())
}, .id = "inflation") %>%
  inner_join(com_orig_intra_score %>% dplyr::rename(com_orig_intra_score = mode), by = "inflation")
cat(" done\n")


gu_partition_stats_file <- file.path(results, paste0("gu_partition_stats_", time_string, ".tsv"))
msg("Saving MCL results stats...")
write_tsv(gu_partition_stats, gu_partition_stats_file, col_names = TRUE)
cat(" done\n")


# Add missing clusters to MCL communities ----------------------------------

# Add the missing nodes to existing communities
# Try to find the best hits for the missing nodes
# Based on coverage/probability/score per position


msg(paste0("Selecting MCL ", best_inflation, " value..."))
gu_mcl_coms <- gu_g_mcl_list[[best_inflation]]
gu_max_com <- max(gu_mcl_coms$com)
cat(" done\n")


gu_mids_l <- gu_missing_ids %>% length()

if (gu_mids_l > 0){
  msg(paste0("Trying to assign ", scales::comma(gu_mids_l), " missing clusters to existing communities...\n"))
  gu_missing_dt <- lo_env$gu_hhblits_missing %>% ungroup()

  # missing_ids -> 7047

  # Try to identify these cluster that can have some homology to existing
  # We only take the queries for missing ids
  # We check that the target is in the MCL communities
  gu_missing_d <- gu_missing_dt %>%
    as.data.table() %>%
    dt_filter(cl_name1 %in% gu_missing_ids) %>%
    dt_filter(cl_name2 %in% gu_mcl_coms$vertex) %>%
    dt_left_join(gu_mcl_coms %>% mutate(vertex = as.character(vertex)) %>% rename(cl_name2 = vertex), by = "cl_name2")

  msg("First pass: Using more relaxed HHBLITS filtering thresholds (prob > 50 and cov >= 40)...")
  # MCL clusters with more relaxed filters (p50 and cov >= 40)
  # we just keep the best hit
  gu_missing_d_pass <- gu_missing_d %>%
    dt_filter(probability >= 50, (q_cov >= 0.4 | t_cov >= 0.4)) %>%
    as_tibble() %>%
    group_by(cl_name1) %>%
    arrange(-probability, -q_cov, -t_cov) %>%
    #mutate(com1 = majority_vote(com)$majority,
    mutate(nhit = row_number()) %>%
    arrange(cl_name1, nhit) %>%
    ungroup() %>%
    #filter(probability >= 30, nhit <= 3) %>%
    filter(nhit == 1) %>%
    select(cl_name1, cl_name2, score_col, com) %>%
    rename(weight = score_col)
  cat(" done\n")
  # Find the ones we couldn't classify
  # 299
  gu_to_assign <- setdiff(gu_missing_ids, gu_missing_d_pass$cl_name1 %>% unique())
  gu_assigned <- setdiff(gu_missing_ids, gu_to_assign)

  # We need to check if any of the remaining no assigned clusters have any homology
  # to the just classified using more relaxed parameters

  msg("Second pass: Find secondary relationships between the newly assigned clusters and missing clusters...")
  gu_missing_d_pass_1 <- bind_rows(lo_env$gu_hhblits_missing %>%
                                     dt_inner_join(gu_missing_d_pass %>% select(cl_name1, com), by = "cl_name1") %>%
                                     dt_filter(cl_name2 %in% gu_to_assign) %>% dt_mutate(cl_name = cl_name2, com_f = com),
                                   lo_env$gu_hhblits_missing %>%
                                     dt_inner_join(gu_missing_d_pass %>% select(cl_name2, com), by = "cl_name2") %>%
                                     dt_filter(cl_name1 %in% gu_to_assign) %>% dt_mutate(cl_name = cl_name1, com_f = com)) %>%
    dt_filter(probability > 50, (q_cov > 0.4 | t_cov > 0.4)) %>%
    group_by(cl_name) %>%
    arrange(-probability, -q_cov, -t_cov) %>%
    #mutate(com1 = majority_vote(com)$majority,
    mutate(nhit = row_number()) %>%
    arrange(cl_name, nhit) %>%
    ungroup() %>%
    #filter(probability >= 30, nhit <= 3) %>%
    filter(nhit == 1) %>%
    dt_select(cl_name, com_f) %>%
    dplyr::rename(com = com_f) %>%
    as_tibble()
  cat(" done\n")

  # The ones that cannot be assigned, we will try to find any existing connections between no classified
  # and run MCL with the best inflation value and create new clusters

  msg("Third pass: Non assigned clusters will create new MCL communities...")
  gu_missing_d_pass_2 <- bind_rows(lo_env$gu_hhblits_missing %>%
                                     dt_filter(cl_name1 %in% gu_to_assign, cl_name2 %in% gu_to_assign) %>%
                                     dt_filter(probability > 50, (q_cov > 0.4 | t_cov > 0.4)) %>%
                                     dt_inner_join(gu_missing_d_pass_1 %>% select(cl_name, com) %>% rename(cl_name1 = cl_name), by = "cl_name1") %>%
                                     dt_mutate(cl_name = cl_name1, com_f = com),
                                   lo_env$gu_hhblits_missing %>%
                                     dt_filter(cl_name1 %in% gu_to_assign, cl_name2 %in% gu_to_assign) %>%
                                     dt_filter(probability > 50, (q_cov > 0.4 | t_cov > 0.4)) %>%
                                     dt_inner_join(gu_missing_d_pass_1 %>% select(cl_name, com) %>% rename(cl_name2 = cl_name), by = "cl_name2") %>%
                                     dt_mutate(cl_name = cl_name2, com_f = com)) %>%
    group_by(cl_name) %>%
    arrange(-probability, -q_cov, -t_cov) %>%
    #mutate(com1 = majority_vote(com)$majority,
    mutate(nhit = row_number()) %>%
    arrange(cl_name, nhit) %>%
    ungroup() %>%
    #filter(probability >= 30, nhit <= 3) %>%
    filter(nhit == 1) %>%
    dt_select(cl_name, com_f) %>%
    dplyr::rename(com = com_f) %>%
    as_tibble()
  cat(" done\n")

  msg("Collecting and aggregating all communities...")
  gu_missing_ids_1 <- setdiff(gu_missing_ids, c(gu_missing_d_pass$cl_name1 %>% unique, gu_missing_d_pass_1$cl_name, gu_missing_d_pass_2$cl_name))

  gu_missing_d_pass_3 <- tibble(cl_name = gu_missing_ids_1) %>% mutate(com = gu_max_com + row_number())

  gu_communities <-  bind_rows(gu_missing_d_pass %>% select(cl_name1, com) %>% rename(cl_name = cl_name1),
                              gu_missing_d_pass_1,
                              gu_missing_d_pass_2,
                              gu_missing_d_pass_3 %>% mutate(cl_name = as.character(cl_name)),
                              gu_mcl_coms %>% rename(cl_name = vertex) %>% mutate(cl_name = as.character(cl_name)) %>% select(cl_name, com)) %>% unique()
  gu_n_comp <- gu_communities$com %>% unique() %>% length()
  gu_n_clus <-  gu_communities$cl_name %>% unique() %>% length()
  cat(" done\n")

}else{
  msg(paste0("Trying to assign ", scales::comma(gu_mids_l), " missing clusters to existing communities...", red("Skipped"), "\n"))
  msg("All clusters already assigned to a component\n")
  gu_communities <- gu_mcl_coms %>% rename(cl_name = vertex) %>% mutate(cl_name = as.character(cl_name)) %>% select(cl_name, com) %>% unique()
  gu_n_comp <- gu_communities$com %>% unique() %>% length()
  gu_n_clus <-  gu_communities$cl_name %>% unique() %>% length()
  cat(" done\n")
}
} else {
  gu_communities <- gu_missing_ids %>% tibble::enframe(name = NULL) %>%
   rename(cl_name=value) %>% mutate(com=row_number(), cl_name=as.character(cl_name))
  gu_n_comp <- gu_communities$com %>% unique() %>% length()
  gu_n_clus <-  gu_communities$cl_name %>% unique() %>% length()
}

msg(paste0("We have been obtained ", scales::comma(gu_n_comp), " from ", scales::comma(gu_n_clus), " clusters\n"))

# Do we have all communities
gu_correct <- length(cl_cat %>% dt_filter(category == "GU") %>% .$cl_name) == gu_n_clus
if(gu_correct){
  msg(paste0("All clusters assigned: ", green(gu_correct), "\n"))
}else{
  msg(paste0("All clusters assigned: ", red(gu_correct), "\n"))
  msg("Something went wrong... quitting\n!\n")
  quit()
}

gu_comp_file <- file.path(results, paste0("gu_communities_", time_string, ".tsv"))
msg(paste0("Saving communities results to ", gu_comp_file, "..."))
readr::write_tsv(gu_communities, path = gu_comp_file, col_names = FALSE)
cat(" done\n\n")

msg(paste0("GU component inference ", green("DONE"), "\n\n"))


###########################################################################
# Get communities for the EU ----------------------------------------------
###########################################################################
if(file.exists(cfg$eu_hhblits_results)){

  msg(paste0("Starting EU component inference\n\n"))

  msg("Reading EU HHBLITS all-vs-all results...")
  lo_env$eu_hhblits_all <- fread(input = cfg$eu_hhblits_results,
    header = FALSE, verbose = cfg$dt_verbose, nThread = cfg$dt_cores)
    names(lo_env$eu_hhblits_all) <- c("cl_name1", "cl_name2", "probability", "e-value", "Score",
    "Cols", "q_start", "q_stop", "t_start", "t_stop", "q_len", "t_len", "q_cov",
    "t_cov")
    cat(paste(" Read", green(scales::comma(nrow(lo_env$eu_hhblits_all))), "entries\n"))

    msg("Filtering EU HHBLITS all-vs-all results with P>=50 and cov > 0.6...")
    lo_env$eu_hhblits <- lo_env$eu_hhblits_all %>% dt_filter(probability > 50, q_cov >= 0.6, t_cov >= 0.6)
    lo_env$eu_hhblits <- lo_env$eu_hhblits %>% dt_filter(cl_name1 != cl_name2)

    # How many clusters do we have (removed self-hits)
    # 522,765 p50;c0.6
    eu_unique_hhblits_cl <- c(lo_env$eu_hhblits$cl_name1, lo_env$eu_hhblits$cl_name2) %>% unique()
    cat(paste0(" Found ", scales::comma(length(eu_unique_hhblits_cl)), " clusters\n"))


    # Missing clusters
    # Num: 27,244
    msg("Identifying missing clusters...")
    eu_missing_ids <- setdiff(cl_cat %>% dt_filter(category == "EU") %>% .$cl_name, eu_unique_hhblits_cl)
    lo_env$eu_hhblits_missing <- lo_env$eu_hhblits_all %>%
    dt_mutate(cl_name1 = as.character(cl_name1), cl_name2 = as.character(cl_name2)) %>%
    dt_filter(cl_name1 %in% eu_missing_ids | cl_name2 %in% eu_missing_ids) %>%
    as_tibble() %>% mutate(score_col = Score/Cols)
    cat(paste0("Found ", scales::comma(length(eu_missing_ids))), "clusters missing\n")

    if(dim(lo_env$eu_hhblits)[1]!= 0){
      # Let's create a graph
      msg("Inferring graph from HHBLITS results...")
      lo_env$eu_hhb_bh <- lo_env$eu_hhblits %>%
      as_tibble() %>%
      mutate(cl_name1 = as.character(cl_name1),
      cl_name2 = as.character(cl_name2)) %>%
      mutate(score_col = Score/Cols)

      eu_hh_g <- lo_env$eu_hhb_bh %>%
      select(cl_name1, cl_name2, score_col) %>%
      rename(weight = score_col) %>%
      igraph::graph_from_data_frame(directed = FALSE) %>%
      igraph::simplify(remove.multiple = TRUE, remove.loops = TRUE, edge.attr.comb = list("max")) %>%
      as_tbl_graph()

      cat(paste0(" Graph contains ", scales::comma(vcount(eu_hh_g)), " vertices and ", scales::comma(ecount(eu_hh_g)), " edges\n"))

      if (cfg$eu_g_mcl_list != ""){
        msg("Loading already computed MCL clusters...")
        load(cfg$eu_g_mcl_list)
        inflation_list <- seq(cfg$mcl_inflation_min, cfg$mcl_inflation_max, cfg$mcl_inflation_step)
        cat(" done\n")
      }else{
        msg("Preparing graph for MCL clustering...")
        gx <- eu_hh_g

        w <- E(gx)$weight
        w <- w - min(w) + 0.001

        E(gx)$weight <- w
        cat(" done\n")

        msg(paste0("Running MCL clustering with inflation values from ", cfg$mcl_inflation_min, " to ", cfg$mcl_inflation_max, " with ", cfg$mcl_inflation_step, " steps...\n"))
        inflation_list <- seq(cfg$mcl_inflation_min, cfg$mcl_inflation_max, cfg$mcl_inflation_step)
        mcl_bin <- cfg$mcl_bin
        eu_g_mcl_list <- parallel::mclapply(inflation_list,
          optimal_mcl,
          G = eu_hh_g,
          mcl_cores = cfg$mcl_cores,
          #ignore.interactive = TRUE,
          Gx = gx,
          #max.vector.size = 1e+07,
          mc.cores = cfg$mcl_jobs,
          mc.cleanup = TRUE,
          mc.silent = TRUE,
          tmp_dir = tmp)

          #k_g_mcl_list <- lapply(inflation_list, optimal_mcl, G = k_hh_g, Gx = gx, mcl_cores = cfg$mcl_cores, tmp_dir = tmp)
          names(eu_g_mcl_list) <- inflation_list

          eu_g_mcl_list_file <- file.path(results, paste0("eu_g_mcl_list_", time_string, ".Rda"))
          msg(paste0("Saving computed MCL communities in ", eu_g_mcl_list_file, "..."))
          save(eu_g_mcl_list, file = eu_g_mcl_list_file, compress = FALSE)
          cat(" done\n")


        }

        # We run mcl with different inflation parameters:
        if (cfg$eu_gc != ""){
          msg("Loading already contracted MCL communities...")
          load(cfg$eu_gc)
          cat(" done\n")
        }else{
          msg("Contracting MCL communities...\n")
          eu_gc <- parallel::mclapply(eu_g_mcl_list,
            contract_graphs,
            G = eu_hh_g,
            #max.vector.size = 3e+09,
            mc.cores = cfg$max_gc_jobs,
            mc.cleanup = TRUE,
            mc.silent = TRUE)
            names(eu_gc) <- inflation_list
            msg("Contracting MCL communities... done\n")

            eu_gc_file <- file.path(results, paste0("eu_gc_", time_string, ".Rda"))
            msg(paste0("Saving contracted MCL communities results in ", eu_gc_file, "..."))
            save(eu_gc, file = eu_gc_file, compress = FALSE)
            cat(" done\n")
          }

          msg("Computing modularity for each MCL community result...")
          eu_hh_g_modularity <- map_df(eu_g_mcl_list, function(X){
            vnames <- V(eu_hh_g)$name
            tibble(modularity = modularity(eu_hh_g, X[match(vnames,X$vertex),]$com))
          }, .id = "inflation")
          cat(" done\n")



          # Inter and intra score modes
          # Calculate nter and intra score modes
          msg("Calculating inter/intra MCL communities score-per-column mode values...\n")

          eu_hh_g_dt <- eu_hh_g %>%
          igraph::as_data_frame(what="edges") %>%
          rename(cl_name1 = from, cl_name2 = to, score_col = weight) %>%
          as.data.table()

          eu_orig_hh_mode <- parallel::mclapply(eu_g_mcl_list, function(X){
            eu_hh_g_dt %>%
            dt_left_join(X %>% dt_select(vertex, com) %>% dt_mutate(com = as.character(com), vertex = as.character(vertex)) %>% rename(cl_name1 = vertex) %>% unique() %>% as.data.table, by = "cl_name1") %>%
            dt_left_join(X %>% dt_select(vertex, com) %>% dt_mutate(com = as.character(com), vertex = as.character(vertex)) %>% rename(cl_name2 = vertex) %>% unique() %>% as.data.table, by = "cl_name2") %>%
            dt_filter(com.x == com.y) %>% mutate(com = com.x) %>% group_by(com) %>% summarise(mode = estimate_mode(score_col))
          }, mc.cores = cfg$max_gc_jobs, #max.vector.size = 3e+09,
          mc.cleanup = TRUE, mc.silent = TRUE)

          com_orig_intra_score <- map_df(eu_orig_hh_mode, function(X){tibble(mode = estimate_mode(X$mode))}, .id = "inflation")
          com_orig_inter_score <- map_df(eu_gc, function(X){tibble(mode = estimate_mode(E(X$graph)$weight))}, .id = "inflation")

          eu_orig_hh_mode_summary <- com_orig_intra_score %>%
          mutate(class = "intra") %>%
          bind_rows(com_orig_inter_score %>% mutate(class = "inter"))
          msg("Calculating inter/intra MCL communities score-per-column mode values... done\n")

          msg("Summarising MCL results...")
          eu_partition_stats <- map_df(eu_g_mcl_list, function(X) {
            tibble(com_orig_n = length(unique(X$com)),
            com_orig_1mem = X %>% group_by(com) %>% count() %>% filter(n == 1) %>% nrow())
          }, .id = "inflation") %>%
          inner_join(com_orig_intra_score %>% dplyr::rename(com_orig_intra_score = mode), by = "inflation")
          cat(" done\n")


          eu_partition_stats_file <- file.path(results, paste0("eu_partition_stats_", time_string, ".tsv"))
          msg("Saving MCL results stats...")
          write_tsv(eu_partition_stats, eu_partition_stats_file, col_names = TRUE)
          cat(" done\n\n")


          # Add missing clusters to MCL communities ----------------------------------

          # Add the missing nodes to existing communities
          # Try to find the best hits for the missing nodes
          # Based on coverage/probability/score per position


          msg(paste0("Selecting MCL ", best_inflation, " value..."))
          eu_mcl_coms <- eu_g_mcl_list[[best_inflation]]
          eu_max_com <- max(eu_mcl_coms$com)
          cat(" done\n")


          eu_mids_l <- eu_missing_ids %>% length()

          if (eu_mids_l > 0){
            msg(paste0("Trying to assign ", scales::comma(eu_mids_l), " missing clusters to existing communities...\n"))
            eu_missing_dt <- lo_env$eu_hhblits_missing %>% ungroup()

            # missing_ids -> 7047

            # Try to identify these cluster that can have some homology to existing
            # We only take the queries for missing ids
            # We check that the target is in the MCL communities
            eu_missing_d <- eu_missing_dt %>%
            as.data.table() %>%
            dt_filter(cl_name1 %in% eu_missing_ids) %>%
            dt_filter(cl_name2 %in% eu_mcl_coms$vertex) %>%
            dt_left_join(eu_mcl_coms %>% mutate(vertex = as.character(vertex)) %>% rename(cl_name2 = vertex), by = "cl_name2")

            msg("First pass: Using more relaxed HHBLITS filtering thresholds (prob > 50 and cov >= 40)...")
            # MCL clusters with more relaxed filters (p50 and cov >= 40)
            # we just keep the best hit
            eu_missing_d_pass <- eu_missing_d %>%
            dt_filter(probability >= 50, (q_cov >= 0.4 | t_cov >= 0.4)) %>%
            as_tibble() %>%
            group_by(cl_name1) %>%
            arrange(-probability, -q_cov, -t_cov) %>%
            #mutate(com1 = majority_vote(com)$majority,
            mutate(nhit = row_number()) %>%
            arrange(cl_name1, nhit) %>%
            ungroup() %>%
            #filter(probability >= 30, nhit <= 3) %>%
            filter(nhit == 1) %>%
            select(cl_name1, cl_name2, score_col, com) %>%
            rename(weight = score_col)
            cat(" done\n")
            # Find the ones we couldn't classify
            # 299
            eu_to_assign <- setdiff(eu_missing_ids, eu_missing_d_pass$cl_name1 %>% unique())
            eu_assigned <- setdiff(eu_missing_ids, eu_to_assign)

            # We need to check if any of the remaining no assigned clusters have any homology
            # to the just classified using more relaxed parameters

            msg("Second pass: Find secondary relationships between the newly assigned clusters and missing clusters...")
            eu_missing_d_pass_1 <- bind_rows(lo_env$eu_hhblits_missing %>%
              dt_inner_join(eu_missing_d_pass %>% select(cl_name1, com), by = "cl_name1") %>%
              dt_filter(cl_name2 %in% eu_to_assign) %>% dt_mutate(cl_name = cl_name2, com_f = com),
              lo_env$eu_hhblits_missing %>%
              dt_inner_join(eu_missing_d_pass %>% select(cl_name2, com), by = "cl_name2") %>%
              dt_filter(cl_name1 %in% eu_to_assign) %>% dt_mutate(cl_name = cl_name1, com_f = com)) %>%
              dt_filter(probability > 50, (q_cov > 0.4 | t_cov > 0.4)) %>%
              group_by(cl_name) %>%
              arrange(-probability, -q_cov, -t_cov) %>%
              #mutate(com1 = majority_vote(com)$majority,
              mutate(nhit = row_number()) %>%
              arrange(cl_name, nhit) %>%
              ungroup() %>%
              #filter(probability >= 30, nhit <= 3) %>%
              filter(nhit == 1) %>%
              dt_select(cl_name, com_f) %>%
              dplyr::rename(com = com_f) %>%
              as_tibble()
              cat(" done\n")

              # The ones that cannot be assigned, we will try to find any existing connections between no classified
              # and run MCL with the best inflation value and create new clusters

              msg("Third pass: Non assigned clusters will create new MCL communities...")
              eu_missing_d_pass_2 <- bind_rows(lo_env$eu_hhblits_missing %>%
                dt_filter(cl_name1 %in% eu_to_assign, cl_name2 %in% eu_to_assign) %>%
                dt_filter(probability > 50, (q_cov > 0.4 | t_cov > 0.4)) %>%
                dt_inner_join(eu_missing_d_pass_1 %>% select(cl_name, com) %>% rename(cl_name1 = cl_name), by = "cl_name1") %>%
                dt_mutate(cl_name = cl_name1, com_f = com),
                lo_env$eu_hhblits_missing %>%
                dt_filter(cl_name1 %in% eu_to_assign, cl_name2 %in% eu_to_assign) %>%
                dt_filter(probability > 50, (q_cov > 0.4 | t_cov > 0.4)) %>%
                dt_inner_join(eu_missing_d_pass_1 %>% select(cl_name, com) %>% rename(cl_name2 = cl_name), by = "cl_name2") %>%
                dt_mutate(cl_name = cl_name2, com_f = com)) %>%
                group_by(cl_name) %>%
                arrange(-probability, -q_cov, -t_cov) %>%
                #mutate(com1 = majority_vote(com)$majority,
                mutate(nhit = row_number()) %>%
                arrange(cl_name, nhit) %>%
                ungroup() %>%
                #filter(probability >= 30, nhit <= 3) %>%
                filter(nhit == 1) %>%
                dt_select(cl_name, com_f) %>%
                dplyr::rename(com = com_f) %>%
                as_tibble()
                cat(" done\n")

                msg("Collecting and aggregating all communities...")
                eu_missing_ids_1 <- setdiff(eu_missing_ids, c(eu_missing_d_pass$cl_name1 %>% unique, eu_missing_d_pass_1$cl_name, eu_missing_d_pass_2$cl_name))

                eu_missing_d_pass_3 <- tibble(cl_name = eu_missing_ids_1) %>% mutate(com = eu_max_com + row_number())

                eu_communities <-  bind_rows(eu_missing_d_pass %>% select(cl_name1, com) %>% rename(cl_name = cl_name1),
                eu_missing_d_pass_1,
                eu_missing_d_pass_2,
                eu_missing_d_pass_3 %>% mutate(cl_name = as.character(cl_name)),
                eu_mcl_coms %>% rename(cl_name = vertex) %>% mutate(cl_name = as.character(cl_name)) %>% select(cl_name, com)) %>% unique()
                eu_n_comp <- eu_communities$com %>% unique() %>% length()
                eu_n_clus <-  eu_communities$cl_name %>% unique() %>% length()
                cat(" done\n")

              }else{
                msg(paste0("Trying to assign ",
                scales::comma(eu_mids_l),
                " missing clusters to existing communities...",
                red("Skipped"), "\n"))
                msg("All clusters already assigned to a component\n")
                eu_communities <- eu_mcl_coms %>% rename(cl_name = vertex) %>% mutate(cl_name = as.character(cl_name)) %>% select(cl_name, com) %>% unique()
                eu_n_comp <- eu_communities$com %>% unique() %>% length()
                eu_n_clus <-  eu_communities$cl_name %>% unique() %>% length()
              }
            } else {
              eu_communities <- eu_missing_ids %>% tibble::enframe(name = NULL) %>%
              rename(cl_name=value) %>% mutate(com=row_number(), cl_name=as.character(cl_name))
              eu_n_comp <- eu_communities$com %>% unique() %>% length()
              eu_n_clus <-  eu_communities$cl_name %>% unique() %>% length()
            }

            msg(paste0("We have been obtained ", scales::comma(eu_n_comp), " from ", scales::comma(eu_n_clus), " clusters\n"))

            # Do we have all communities
            eu_correct <- length(cl_cat %>% dt_filter(category == "EU") %>% .$cl_name) == eu_n_clus
            if(eu_correct){
              msg(paste0("All clusters assigned:", green(eu_correct), "\n"))
            }else{
              msg(paste0("All clusters assigned:", red(eu_correct), "\n"))
              msg("Something went wrong... quitting\n!")
              quit()
            }

            eu_comp_file <- file.path(results, paste0("eu_communities_", time_string, ".tsv"))
            msg(paste0("Saving communities results to ", eu_comp_file, "..."))
            readr::write_tsv(eu_communities, path = eu_comp_file, col_names = FALSE)
            cat(" done\n\n")

            msg(paste0("EU component inference ", green("DONE"), "\n\n"))

            all_comp_file <- file.path(results, paste0("all_communities_", time_string, ".tsv"))
            msg(paste0("Saving all communities results to ", all_comp_file, "..."))
            bind_rows(k_communities %>% mutate(com = paste0("k_c_", com), category = "k"),
            kwp_communities %>% mutate(com = paste0("kwp_c_", com), category = "kwp"),
            gu_communities %>% mutate(com = paste0("gu_c_", com), category = "gu"),
            eu_communities %>% mutate(com = paste0("eu_c_", com), category = "eu")
          ) %>%
          write_tsv(path = all_comp_file, col_names = TRUE)
          cat(" done\n")

          all_comp_file1 <- file.path(dirname(results), paste0("cluster_communities.tsv"))
          msg(paste0("Saving all communities results to ", all_comp_file1, "..."))
          bind_rows(k_communities %>% mutate(com = paste0("k_c_", com), category = "k"),
          kwp_communities %>% mutate(com = paste0("kwp_c_", com), category = "kwp"),
          gu_communities %>% mutate(com = paste0("gu_c_", com), category = "gu"),
          eu_communities %>% mutate(com = paste0("eu_c_", com), category = "eu")
        ) %>%
        write_tsv(path = all_comp_file1, col_names = TRUE)
        cat(" done\n")
}

all_comp_file1 <- file.path(dirname(results), paste0("cluster_communities.tsv"))
msg(paste0("Saving all communities results to ", all_comp_file1, "..."))
bind_rows(k_communities %>% mutate(com = paste0("k_c_", com), category = "k"),
kwp_communities %>% mutate(com = paste0("kwp_c_", com), category = "kwp"),
gu_communities %>% mutate(com = paste0("gu_c_", com), category = "gu")) %>%
write_tsv(path = all_comp_file1, col_names = TRUE)
cat(" done\n")
