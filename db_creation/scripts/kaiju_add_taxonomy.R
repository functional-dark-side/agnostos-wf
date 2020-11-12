#!/usr/bin/env Rscript
# Check if basic packages are installed -----------------------------------

is.installed <- function(pkg){
  is.element(pkg, installed.packages()[,1])
}

if (!is.installed("crayon") || !is.installed("config") || !is.installed("taxonomizr")){
  cat("We will try to install the packages crayon, optparse and config... (this will be only be done once)\n")
  Sys.sleep(5)
  if (!is.installed("crayon")){
    suppressMessages(install.packages("crayon", repos = "http://cran.us.r-project.org"))
  }
  if (!is.installed("config")){
    suppressMessages(install.packages("config", repos = "http://cran.us.r-project.org"))
  }
  if (!is.installed("taxonomizr")){
    suppressMessages(install.packages("taxonomizr", repos = "http://cran.us.r-project.org"))
  }
}

suppressMessages(library(crayon))
suppressMessages(library(optparse))


# Check if packaged are installed -----------------------------------------

cat("\nChecking if all packages are installed...\n\n")

needed <- c("tidyverse", "data.table","maditr", "optparse", "taxonomizr")

missing_package <- FALSE
# For loop to run through each of the packages
for (p in 1:length(needed)){
  if(is.installed(needed[p])){
    cat(sprintf("%-12s: %s", needed[p], green("Installed\n")))
  }else{
    cat(sprintf("%-12s: %s", needed[p], red("Not installed\n")))
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
  make_option(c("-r", "--results"), type="character", default=NULL,
              help="Kraken2/Kaiju results path", metavar="character"),
  make_option(c("-o", "--output"), type="character", default=NULL,
              help="Kraken2 results with taxonomy", metavar="character"),
  make_option(c("-n", "--names"), type="character", default=NULL,
              help="GTDB-NCBI-like names.dmp", metavar="character"),
  make_option(c("-d", "--nodes"), type="character", default=NULL,
              help="GTDB-NCBI-like nodes.dmp", metavar="character"),
  make_option(c("-p", "--threads"), type="numeric", default=1,
              help="Threads", metavar="numeric"),
  make_option(c("-t", "--tmp"), type="character", default=NULL,
              help="TMP path", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$results) | is.null(opt$output) | is.null(opt$names) | is.null(opt$nodes)){
  print_help(opt_parser)
  stop("You need to provide the path to the Kraken2/Kaiju results, Kraken2/Kaiju ouput and the dmp files.\n", call.=FALSE)
}

setDTthreads(opt$threads)
options(datatable.verbose = FALSE)

TMP <- file.path(opt$tmp)
if(!dir.exists(path.expand(TMP))) dir.create(path.expand(TMP))
#unixtools::set.tempdir(path.expand(TMP))

# Load taxonomy data
msg("Loadind .dmp files...")
taxaNodes <- read.nodes.sql(opt$nodes, sqlFile = file.path(TMP, "nodesNames.sqlite"), overwrite=TRUE)
taxaNames <- read.names.sql(opt$names, sqlFile = file.path(TMP, "nodesNames.sqlite"), overwrite=TRUE)
cat(" done\n")
# Read results
msg("Reading kraken2/kaiju results...")
results <- fread(opt$results, header = FALSE, showProgress = FALSE, nThread = opt$threads, fill = TRUE, sep = "\t")
cat(" done\n")

msg("Adding taxonomy...")
classified <- results %>%
  filter(V1 == 'C')
classified <- cbind(classified, getTaxonomy(classified$V3,taxaNodes,taxaNames))
cat(" done\n")

msg("Writing results...")
fwrite(x = classified, file = opt$output, sep = "\t", col.names = FALSE)
cat(" done\n")
