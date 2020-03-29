# Shingling function, specify k-shingles
Shingling <- function(document, k) {
  shingles <- character( length = length(document) - k + 1 )

  for( i in 1:( length(document) - k + 1 ) ) {
    shingles[i] <- paste( document[ i:(i + k - 1) ], collapse = " " )
  }

  return( unique(shingles) )
}

# Jaccard similarity
JaccardSimilarity <- function(x, y) {
  non_zero <- which(x | y)
  set_intersect <- sum( x[non_zero] & y[non_zero] )
  set_union <- length(non_zero)
  return(set_intersect / set_union)
}

# Jaccard function in case Shingling is not needed, specify type of annotations: Pfam domains of clans
Jaccard <- function(doc_annot){
  doc_annot <- setNames(doc_annot , c("memb","annot","pres"))
  J <- spread(doc_annot,annot,pres, fill=0)
  rownames(J) <- J$memb
  J <- J[-1]
  pr_DB$set_entry( FUN = JaccardSimilarity, names = c("JaccardSimilarity") )
  # jaccard similarity distance matrix
  ds <- dist(J, method = "JaccardSimilarity" )
  # delete the new entry
  pr_DB$delete_entry("JaccardSimilarity")
  return(ds)
}

# Cluster with multi-domain annotations
MultiAnnot <- function(annot_data, annot, clstrnoNA, size){
  # In case there are some member with mono-domain annotations (no shingling)
  singl <- annot_data %>% filter(is.na(annot_data[,4])==T) %>% mutate(pres=1) %>% select(memb,`1`,pres) %>%
    setNames(c("memb","annot","pres"))
  if(dim(singl)[1] > 0){
    n.singl <- dim(singl)[1]
    prop.singl <- n.singl/size
    ds <- Jaccard(singl)
    if(n.singl==1){
      median <- 1
    } else {
      median <- median(as.matrix(ds)[lower.tri(as.matrix(ds))])
    }
    rep <- clstrnoNA$rep[1]
    median.sc <- median*prop.singl
    type=paste("Singl_",annot,sep = "")
    prop.type=prop.singl
    partial.prop <- dim(merge(singl, clstrnoNA,by="memb") %>% filter(partial!="00"))[1]/n.singl
    res.s <- data.frame(rep=rep,
                        jacc_median_raw=median,
                        jacc_median_sc=median.sc,
                        type=type,
                        prop_type=prop.type,
                        prop_partial=partial.prop,
                        stringsAsFactors =F)
    # Multi-domains annotations
    multi <- annot_data %>% filter(is.na(annot_data[,4])==F)
    m <- dim(multi)[2] - 1
    multi.1 <- multi %>% select(-1,-2)
    multi.1 <- split(multi.1, seq(nrow(multi.1)))
    # "shingle" our example document, with k = n-grams, we choose a k=2
    multi.sh <- lapply(multi.1, function(x){
      Shingling(x, k = 2) })
    # unique shingles sets across all documents
    doc_dict <- unlist(multi.sh) %>% unique()
    # "characteristic" matrix
    M <- lapply(multi.sh, function(set, dict) {
      as.integer(dict %in% set)
    }, dict = doc_dict) %>% data.frame()
    # set the names for both rows and columns
    setnames(M, paste( "doc", 1:length(multi.sh), sep = "_" ) )
    rownames(M) <- doc_dict
    # how similar is two given document, jaccard similarity
    # create a new entry in the registry
    pr_DB$set_entry( FUN = JaccardSimilarity, names = c("JaccardSimilarity") )
    # jaccard similarity distance matrix
    dm <- dist(t(M), method = "JaccardSimilarity" )
    # delete the new entry
    pr_DB$delete_entry("JaccardSimilarity")

    n.multi <- dim(multi)[1]
    prop.multi <- n.multi/size
    if(n.multi==1){
      median <- 1
    } else {
      median <- median(as.matrix(dm)[lower.tri(as.matrix(dm))])
    }
    rep <- clstrnoNA$rep[1]
    median.sc <- median*prop.multi
    type=paste("Multi_",annot,sep = "")
    prop.type=prop.multi
    partial.prop <- dim(merge(multi, clstrnoNA,by="memb") %>% filter(partial!="00"))[1]/n.multi
    res.m <- data.frame(rep=rep,
                        jacc_median_raw=median,
                        jacc_median_sc=median.sc,
                        type=type,
                        prop_type=prop.type,
                        prop_partial=partial.prop,
                        stringsAsFactors =F)
    res <- rbind(res.s,res.m)
  }else{    # only Multi domains
    #evaluate 2-grams (or K-shingling with k=2)
    multi <- annot_data #filter(is.na(m2[,3])==F) %>% no necessary
    m <- dim(multi)[2] - 1
    multi.1 <- multi %>% select(-1)
    multi.1 <- split(multi.1, seq(nrow(multi.1)))
    # "shingle" our example document, with k = n-grams, we choose a k=2
    multi.sh <- lapply(multi.1, function(x){
      Shingling(x, k = 2) })

    # unique shingles sets across all documents
    doc_dict <- unlist(multi.sh) %>% unique()

    # "characteristic" matrix
    M <- lapply(multi.sh, function(set, dict) {
      as.integer(dict %in% set)
    }, dict = doc_dict) %>% data.frame()

    # set the names for both rows and columns
    setnames( M, paste( "doc", 1:length(multi.sh), sep = "_" ) )
    rownames(M) <- doc_dict
    # how similar is two given document, jaccard similarity
    # create a new entry in the registry
    pr_DB$set_entry( FUN = JaccardSimilarity, names = c("JaccardSimilarity") )
    # jaccard similarity distance matrix
    dm <- dist(t(M), method = "JaccardSimilarity" )
    # delete the new entry
    pr_DB$delete_entry("JaccardSimilarity")

    n.multi <- dim(multi)[1]
    prop.multi <- n.multi/size
    if(n.multi==1){
      median <- 1
    } else {
      median <- median(as.matrix(dm)[lower.tri(as.matrix(dm))])
    }
    type <- paste("Multi_",annot,sep = "")
    prop.type=prop.multi
    rep <- clstrnoNA$rep[1]
    median.sc <- median*prop.multi
    partial.prop <- dim(merge(multi, clstrnoNA,by="memb") %>% filter(partial!="00"))[1]/n.multi
    res <- data.frame(rep=rep,
                        jacc_median_raw=median,
                        jacc_median_sc=median.sc,
                        type=type,
                        prop_type=prop.type,
                        prop_partial=partial.prop,
                        stringsAsFactors =F)
  }
  return(res)
}

