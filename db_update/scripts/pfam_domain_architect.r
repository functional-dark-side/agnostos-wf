library(data.table)
library(tidyverse)
library(stringr)

args <- commandArgs(TRUE)

print(args)

annot_kept  <- fread(paste("zcat",args[1],sep=" "),
                    stringsAsFactors = F, header = F) %>%
  dplyr::select(V1,V2,V3,V4,V5,V6,V7) %>%
  setNames(c("cl_name","repres","rep_annot","memb","memb_annot","clan","partial"))

annot_kept1 <- annot_kept %>% mutate(class=ifelse(repres==memb,"rep","memb")) %>% select(-repres, -memb)
#select repres annot and info
annot_kept_rep <- annot_kept1 %>% filter(class=="rep") %>% select(cl_name,rep_annot,clan,partial)
#summary of the cluster
prop_partial <- annot_kept %>% select(cl_name,partial) %>% mutate(par=ifelse(partial=="00",0,1), count=1) %>%
  group_by(cl_name) %>% mutate(prop_partial=sum(par)/sum(count)) %>% distinct(cl_name, prop_partial)

annot_summ <- annot_kept %>% select(cl_name,memb_annot,clan) %>%
  drop_na() %>% group_by(cl_name,memb_annot,clan) %>% count()

annot_duf <- annot_kept %>% select(cl_name,memb_annot) %>% drop_na() %>% mutate(categ=ifelse(grepl('DUF',memb_annot),"duf","pf")) %>%
  group_by(cl_name,categ) %>% count() %>% ungroup() %>% group_by(cl_name) %>% mutate(p_categ=n/sum(n)) %>%
  distinct(cl_name,categ,p_categ) %>% top_n(1,p_categ) %>% filter(p_categ>0.5 | p_categ==0.5 & categ=="pf")
# Collect annot of repres and of members
annot_summ_rep_all <- merge(annot_duf, annot_kept_rep, by="cl_name")
annot_summ_rep_all <- merge(annot_summ_rep_all, prop_partial, by="cl_name")
annot_summ_rep_all <- merge(annot_summ_rep_all, annot_summ, by="cl_name")
#good_summ_rep_all$V3[is.na(good_summ_rep_all$V3)] <- as.character(good_summ_rep_all$V4[is.na(good_summ_rep_all$V3)])
annot_summ_rep_all <- annot_summ_rep_all %>% group_by(cl_name) %>% mutate(p_annot=n/sum(n)) %>% add_count(cl_name)
#select only the most abundant annotation for each cluster
annot_summ_rep <- top_n(annot_summ_rep_all, 1,p_annot)
# if repres without annot replace it with the most abundant one in the cluster
annot_summ_rep$rep_annot[is.na(annot_summ_rep$rep_annot)] <- as.character(annot_summ_rep$memb_annot[is.na(annot_summ_rep$rep_annot)])
# repres and most abundant annot are the same
equal <- annot_summ_rep %>% filter(memb_annot==rep_annot) %>%
  select(-n) %>% group_by(cl_name) %>% slice(1) %>%
  setNames(c("cl_name","pf_categ","perc_pf_categ","rep_annot","rep_clan","rep_partial","prop_partial","memb_annot","memb_clan",
             "prop_prev_annot"))

#they are different --> choose the repres one or the most abundant??
not_equal <- annot_summ_rep %>% filter(!cl_name %in% equal$cl_name)%>%
  select(-n) %>% group_by(cl_name) %>% slice(1) %>%
  setNames(c("cl_name","pf_categ","perc_pf_categ","rep_annot","rep_clan","rep_partial","prop_partial","memb_annot","memb_clan",
             "prop_prev_annot"))
not_equal <- not_equal %>% group_by(cl_name) %>% mutate(cl_annot=ifelse(rep_clan==memb_clan || rep_partial=="00",rep_annot,memb_annot)) %>% slice(1)

cl_pfam_arch <- rbind(equal %>% rename(cl_annot=memb_annot) %>% distinct(cl_name,cl_annot,pf_categ,perc_pf_categ),
                      not_equal %>% select(cl_name,cl_annot,pf_categ,perc_pf_categ))

write.table(cl_pfam_arch, paste(args[2],"pfam_domain_architect.tsv",sep="/"), col.names = T, row.names = F, quote = F, sep = "\t")

kept_duf <- cl_pfam_arch %>% filter(pf_categ=="duf" & perc_pf_categ==1)
write.table(kept_duf, paste(args[2],"kept_DUFs.tsv",sep="/"), col.names = T, row.names = F, quote = F, sep = "\t")
kept_pf <- cl_pfam_arch %>% anti_join(kept_duf,by="cl_name")
write.table(kept_pf, paste(args[2],"kept_PF.tsv",sep="/"), col.names = T, row.names = F, quote = F, sep = "\t")

