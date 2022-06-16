#!/usr/bin/env Rscript

# Check if basic packages are installed -----------------------------------

is.installed <- function(pkg){
  is.element(pkg, installed.packages()[,1])
}

if (!is.installed("proxy") || !is.installed("textreuse") || !is.installed("vegan")){
  cat("We will try to install the packages... (this will be only be done once)\n")
  Sys.sleep(5)
  if (!is.installed("proxy")){
    suppressMessages(install.packages("proxy", repos = "http://cran.us.r-project.org"))
  }
  if (!is.installed("textreuse")){
    suppressMessages(install.packages("textreuse", repos = "http://cran.us.r-project.org"))
  }
  if (!is.installed("vegan")){
    suppressMessages(install.packages("vegan", repos = "http://cran.us.r-project.org"))
  }
}

library(tidyverse)
library(proxy)
library(stringr)
library(data.table)
library(textreuse)
library(vegan)
library(parallel)
library(optparse)

# Script command line options ---------------------------------------------

option_list <- list(
  make_option(c("-i", "--input"), type="character", default=NULL,
              help="Annotated clusters path", metavar="character"),
  make_option(c("-p", "--pfam_terms"), type="character", default=NULL,
              help="Pfam shared terms path", metavar="character"),
  make_option(c("-o", "--output"), type="character", default=NULL,
              help="Functional validation results", metavar="character"),
  make_option(c("-f", "--functions"), type="character", default=NULL,
              help="Script with extra functions path", metavar="character"),
  make_option(c("-t", "--threads"), type="numeric", default=1,
              help="Number of threads", metavar="numeric")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$input) | is.null(opt$output) |
    is.null(opt$pfam_terms) | is.null(opt$functions) |
    is.null(opt$threads)){
  print_help(opt_parser)
  stop("You need to provide the path to the annotated clusters and pfam terms, and to the output file\n", call.=FALSE)
}

setDTthreads(opt$threads)
options(datatable.verbose = FALSE)
options(readr.num_columns = 0)

#Validation of clusters with pfam annotations
# Read annotated clusters
cluster_annot <- read_tsv(opt$input, col_names = F) %>%
  setNames(c("cl_name","rep","memb","annot","acc","clan","partial")) %>%
  dplyr::select(rep,memb,partial,annot,clan)

# Pfam terminal (C, N) domains of same proteins
pfam_shared_term <- read.table(opt$pfam_terms,
                               stringsAsFactors = F, header = T)

# Cluster lists
cluster.list <- split(cluster_annot, list(cluster_annot$rep), drop=TRUE)

rm(cluster_annot)
gc()

# Main function
shingl.jacc <- function(clstr, pfam){
  # functional annotation data per cluster
  test_multi <- as.data.frame(clstr) %>% setNames(c("rep","memb","partial","annot","clan"))
  # cluster size
  size <- dim(test_multi)[1]
  # remove not annotated memeber
  clstrnoNA <- test_multi %>% tidyr::drop_na()
  # number of annotated members (cluster size)
  annotated <- dim(clstrnoNA)[1]
  prop.annot <- annotated/size
  #select members and annotations (pfam domains and clans)
  m1 <- clstrnoNA %>% dplyr::select(memb,annot,clan)
  ma <- max(str_count(m1$annot, "\\|")) + 1 # max number of multiple domains on the same orfs
  mc <- max(str_count(m1$clan, "\\|")) + 1 # max number of multiple clans on the same orfs

  # Homog_pf: all members annotated to the same Pfam domain
  if(length(unique(m1$annot))==1){
      m1 <- m1 %>% dplyr::select(memb,annot) %>% mutate(pres=1) %>% distinct()
      ds <- Jaccard(m1)
      if(dim(m1)[1]==1){  #Only one annotated member
        median <- 1
      }else {
        median <- median(as.matrix(ds)[lower.tri(as.matrix(ds))])
      }
      rep <- clstrnoNA$rep[1]
      median.sc <- median*prop.annot
      type="Homog_pf"
      prop.type = prop.annot
      partial.prop <- dim(merge(m1, clstrnoNA, by="memb") %>% filter(partial!="00"))[1]/annotated
      res <- data.frame(rep=rep,
                        jacc_median_raw=median,
                        jacc_median_sc=median.sc,
                        type=type,
                        prop_type=prop.type,
                        prop_partial=partial.prop,
                        stringsAsFactors =F)
  }
  # Homog_clan: all members annotated to the same Pfam clan
  else if(length(unique(m1$clan))==1 & any(m1$clan!="no_clan")==T){
      m1 <- m1 %>% dplyr::select(memb,clan) %>% mutate(pres=1) %>% distinct()
      ds <- Jaccard(m1)
      median <- median(as.matrix(ds)[lower.tri(as.matrix(ds))])
      rep <- clstrnoNA$rep[1]
      median.sc <- median*prop.annot
      type="Homog_clan"
      prop.type = prop.annot
      partial.prop <- dim(merge(m1, clstrnoNA,by="memb") %>% filter(partial!="00"))[1]/annotated
      res <- data.frame(rep=rep,
                        jacc_median_raw=median,
                        jacc_median_sc=median.sc,
                        type=type,
                        prop_type=prop.type,
                        prop_partial=partial.prop,
                        stringsAsFactors =F)
      }
  # Not homog, check for Pfam different terminal domains of the same protein
  else{
        pfam_term <- pfam_shared_term
        # if the cluster contains any member with multiple annotation
        # split the annotation in multiple rows
        if(ma>1){
          multi_annot <- strsplit(m1$annot, split = "\\|")
          m1 <- data.frame(memb = rep(m1$memb, sapply(multi_annot, length)), annot = unlist(multi_annot), clan=rep(m1$clan, sapply(multi_annot, length)), stringsAsFactors = F)
        }
        # annotation metches with the list of Pfam terminal domains of same proteins
        term <- m1 %>% left_join(pfam_term,by=c("annot"="pfam")) %>% filter(annot != cname)
        # if there is a correspondance, replace terminal-domain names with the common ones
        if(dim(term)[1]>0){
        m1 <- m1 %>% mutate(annot=plyr::mapvalues(as.vector(.$annot), from = pfam_term$pfam, to = pfam_term$cname))
        }
        # Homog_pf_term: using the common name, do we have homogeneous annotations?
        if(length(unique(m1$annot))==1){
            m1 <- m1 %>% dplyr::select(memb,annot) %>% mutate(pres=1) %>% distinct
            ds <- Jaccard(m1)
            median <- median(as.matrix(ds)[lower.tri(as.matrix(ds))])
            rep <- clstrnoNA$rep[1]
            median.sc <- median*prop.annot
            type="Homog_pf_term"
            prop.type = prop.annot
            partial.prop <- dim(merge(m1, clstrnoNA,by="memb") %>% filter(partial!="00"))[1]/annotated
            res <- data.frame(rep=rep,
                              jacc_median_raw=median,
                              jacc_median_sc=median.sc,
                              type=type,
                              prop_type=prop.type,
                              prop_partial=partial.prop,
                              stringsAsFactors =F)
        }
        # Still not homogeneous annotations..
        else{
          # Only mono annotations
            if(ma==1 & mc==1){
              m1.a <- m1 %>% dplyr::select(memb,annot) %>% mutate(pres=1) %>% distinct
              ds <- Jaccard(m1.a)
              median.a <- median(as.matrix(ds)[lower.tri(as.matrix(ds))])
              rep <- clstrnoNA$rep[1]
              median.sc.a <- median.a*prop.annot
              m1.c <- m1 %>% dplyr::select(memb,clan) %>% mutate(pres=1) %>% distinct
              ds <- Jaccard(m1.c)
              median.c <- median(as.matrix(ds)[lower.tri(as.matrix(ds))])
              rep <- clstrnoNA$rep[1]
              median.sc.c <- median.c*prop.annot
               if(median.sc.a >= median.sc.c){
                type="Mono_pf"
                prop.type = prop.annot
                partial.prop <- dim(merge(m1, clstrnoNA,by="memb") %>% filter(partial!="00"))[1]/annotated
                res <- data.frame(rep=rep,
                                  jacc_median_raw=median.a,
                                  jacc_median_sc=median.sc.a,
                                  type=type,
                                  prop_type=prop.type,
                                  prop_partial=partial.prop,
                                  stringsAsFactors =F)
               }else{
                type="Mono_clan"
                prop.type = prop.annot
                partial.prop <- dim(merge(m1, clstrnoNA,by="memb") %>% filter(partial!="00"))[1]/annotated
                res <- data.frame(rep=rep,
                                  jacc_median_raw=median.c,
                                  jacc_median_sc=median.sc.c,
                                  type=type,
                                  prop_type=prop.type,
                                  prop_partial=partial.prop,
                                  stringsAsFactors =F)
               }
              }else{ # The cluster presents also multiple annotations
                # split multiple annotations (in many columns), fill empty cells with NAs
                m1.c <- m1 %>% dplyr::select(memb,clan) %>% mutate(pres=1) %>% distinct
                m2.c <- cbind(m1.c,str_split_fixed(m1.c$clan, "\\|", mc))
                m2.c[,2] <- NULL
                empty_as_na <- function(x){
                  if("factor" %in% class(x)) x <- as.character(x) ## since ifelse wont work with factors
                  ifelse(as.character(x)!="", x, NA)
                }
                m2.c <- m2.c %>% mutate_all(funs(empty_as_na))
                #evaluate 2-grams (or K-shingling with k=2)
                #remove doc with one word-domain (n_domain==1), keep them separated as Singl_pf/clan
                res.c <- MultiAnnot(m2.c,"clan", clstrnoNA=clstrnoNA, size=size)
                #with annot
                m1.a <- m1 %>% dplyr::select(memb,annot) %>% group_by(memb) %>%
                  mutate(annot=paste(annot,collapse = "|")) %>% mutate(pres=1) %>%
                  distinct %>% ungroup()
                m2.a <- cbind(m1.a,str_split_fixed(m1.a$annot, "\\|", ma))
                m2.a[,2] <- NULL
                empty_as_na <- function(x){
                  if("factor" %in% class(x)) x <- as.character(x) ## since ifelse wont work with factors
                  ifelse(as.character(x)!="", x, NA)
                }
                m2.a <- m2.a %>% mutate_all(funs(empty_as_na))
                #evaluate 2-grams (or K-shingling with k=2)
                #remove doc with one word (n_domain==1)..keep them separated
                res.a <- MultiAnnot(m2.c,"pf",clstrnoNA=clstrnoNA, size=size)
                # Choose the best results between domain-annot and clans
                if(any(!grepl("no_clan",m1$clan))==T && any(as.numeric(res.a$jacc_median_sc)< as.numeric(res.c$jacc_median_sc))==T){
                  res <- res.c
                }else{
                  res <- res.a
                }
              }
        }
  }
  return(res)
}

source(opt$functions)
cores <- (as.numeric(opt$threads) - 2)
res.list <- mclapply(cluster.list, shingl.jacc, pfam=pfam_shared_term, mc.cores = cores)
results <- plyr::ldply(res.list, data.frame)
res.parsed.1 <- results %>% select(rep,jacc_median_raw,jacc_median_sc,type,prop_type,prop_partial) %>% group_by(rep) %>%
  mutate(count = n_distinct(type)) %>% filter(count==1)
res.parsed.2 <- results %>% select(rep,jacc_median_raw,jacc_median_sc,type,prop_type,prop_partial) %>% group_by(rep) %>%
  mutate(count = n_distinct(type)) %>% filter(count==2) %>%
  filter(prop_type==max(prop_type)) %>%
  filter(prop_partial==min(prop_partial)) %>%
  filter(jacc_median_sc==max(jacc_median_sc)) %>%
  filter(jacc_median_raw==max(jacc_median_raw)) %>% group_by(rep) %>% slice(1)

shingl_jacc.res <- rbind(res.parsed.1,res.parsed.2)
shingl_jacc.res <- shingl_jacc.res %>% select(-count) %>%
        mutate(annot_categ=ifelse(grepl('Homolog',type), "HA",
                ifelse(grepl('Mono',type),"MoDA","MuDA")))

rm(res.list,results,res.parsed.1,res.parsed.2)
gc()

# Write results
write.table(shingl_jacc.res, opt$output,col.names=T,row.names=F,sep="\t",quote=F)
