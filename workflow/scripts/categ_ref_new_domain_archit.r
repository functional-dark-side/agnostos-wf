#!/usr/bin/env Rscript

library(data.table)
library(tidyverse)
library(stringr)

args <- commandArgs(TRUE)

print(args)

annot_kept <- fread(args[1], stringsAsFactors = F, header = F) %>%
    select(V1, V4) %>%
    setNames(c("cl_name", "annot"))

annot_pf <- annot_kept %>%
    separate_rows(annot, sep = "\\|") %>%
    group_by(cl_name) %>%
    mutate(type = ifelse(all(!grepl("^DUF", annot)), "PF", "DUF")) %>%
    group_by(cl_name, type) %>%
    mutate(original = paste(annot, collapse = "|")) %>%
    distinct()

annot_K <- annot_kept %>%
    filter(cl_name %in% (annot_pf %>% filter(type == "PF") %>%
        .$cl_name))

write.table(annot_K, paste(args[2], "_new_k_ids_annot.tsv", sep = ""),
    col.names = F, row.names = F, quote = F, sep = "\t"
)

annot_GU <- annot_kept %>%
    filter(!cl_name %in% annot_K$cl_name)

write.table(annot_GU, paste(args[2], "_new_gu_ids_annot.tsv", sep = ""),
    col.names = F, row.names = F, quote = F, sep = "\t"
)