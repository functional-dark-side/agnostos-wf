#!/usr/bin/env Rscript

library(tidyverse)
library(data.table)
library(optparse)


# Script command line options ---------------------------------------------

option_list <- list(
  make_option(c("-n", "--name"),
    type = "character", default = NULL,
    help = "Data name", metavar = "character"
  ),
  make_option(c("-c", "--comm"),
    type = "character", default = NULL,
    help = "New cluster communities", metavar = "character"
  ),
  make_option(c("-o", "--ocomm"),
    type = "character", default = NULL,
    help = "Original cluster communities", metavar = "character"
  ),
  make_option(c("-i", "--icomm"),
    type = "character", default = NULL,
    help = "Integrated cluster communities", metavar = "character"
  ),
  make_option(c("-s", "--shared"),
    type = "character", default = NULL,
    help = "Are shared cluster revalidated?", metavar = "character"
  )
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$name) | is.null(opt$comm) |
  is.null(opt$ocomm) | is.null(opt$icomm) |
  is.null(opt$shared)) {
  print_help(opt_parser)
  stop("You need to provide the path to the new and original cluster communities, and output paths\n", call. = FALSE)
}

options(datatable.verbose = FALSE)

data_name=as.character(opt$name)

comm <- fread(opt$comm,header=T,sep="\t") %>%
  mutate(cl_name=as.character(cl_name)) %>%
  mutate(com=paste0(com,"_",data_name))

ocomm <- fread(opt$ocomm, header=T, sep="\t") %>%
  mutate(cl_name=as.character(cl_name))

shared=as.character(opt$shared)

if(shared == "true"){
  revalid_com <- ocomm %>% filter(cl_name %in% comm$cl_name) %>% .$com %>% unique()
  ocomm_filt <- ocomm %>% filter(!com %in% revalid_com)
  old_cl_new_comm <- ocomm %>% filter(com %in% revalid_com) %>%
    left_join(comm, by="cl_name") %>%
    group_by(com.x) %>%
    tidyr::fill(category.y, .direction="downup") %>%
    mutate(com=ifelse(category.y==category.x & is.na(com.y), NA,
                      ifelse(category.y==category.x & !is.na(com.y), com.y,
                             ifelse(category.y!=category.x & is.na(com.y), com.x, com.y)))) %>%
    tidyr::fill(com, .direction="downup") %>% ungroup() %>%
    select(cl_name,com) %>% mutate(category=gsub('_.*$','',com))

  icomm <- comm %>% filter(!cl_name %in% old_cl_new_comm$cl_name) %>%
    rbind(ocomm_filt) %>% rbind(old_cl_new_comm)
}else{
  icomm <- ocomm %>% filter(!cl_name %in% comm$cl_name) %>%
    rbind(comm)
}

write.table(icomm, opt$icomm, row.names=F, sep="\t", quote=F)
