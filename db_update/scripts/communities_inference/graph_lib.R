require(tidyverse)
require(igraph)


# Run MCL clustering ------------------------------------------------------

run_mcl <- function(X, graph = graph, options = NULL, mcl_bin = mcl_bin, tmp_dir = "/tmp") {
  if (is.null(options)) {
    mcl_options <- "--abc -I 2.0"
  } else {
    mcl_options <- options
  }
  R <- as.character(stringi::stri_rand_strings(1, 10))
  out_file <- file.path(tmp_dir, paste(X, "-mcl_", R, ".out", sep = ""))
  g_name <- file.path(tmp_dir, paste(X, "_", R, ".tsv", sep = ""))
  mcl_command <- paste(mcl_bin, g_name, mcl_options, "-o", out_file, "> /dev/null", sep = " ")

  # convert names to numeric
  graph <- graph %>%
    igraph::as_data_frame(what = "edges") %>%
    as_tibble()

  data.table::fwrite(graph %>% data.table::as.data.table(),
                     file = g_name, col.names = F, showProgress = T, sep = " ")

  # V(graph)$id <- V(graph)$name write.graph(graph, g_name, format = 'ncol')

  system(mcl_command, ignore.stdout = TRUE, ignore.stderr = TRUE)

  # map_file <- paste(out_folder, '/', X, '.map', sep = '')
  mcl_results <- parse_mcl(out_file) #%>% inner_join(graph %>% activate(nodes) %>%
  unlink(c(g_name, out_file))
  mcl_results
}



# Parse MCL results -------------------------------------------------------

parse_mcl <- function(X) {
  modules <- fread(input = X, sep = "\t", header = FALSE, fill = TRUE) %>% as_tibble()
  cat(paste("Parsing MCL output.", nrow(modules), "communities detected.\n"))

  com <- modules %>%
    mutate(com = row_number()) %>%
    gather(nr, vertex, -com) %>%
    filter(!is.na(vertex)) %>%
    select(-nr) %>%
    arrange(vertex)
  com
}



# Funtions for processing the graphs --------------------------------------



calculate_intra_score_m <- function(G = G) {
  cat("Preparing data...\n\tGetting nodes...")
  nodes <- G %>% as_data_frame(what = "vertices") %>% as_tibble()
  cat(" done\n\tGetting edges...")
  edges <- G %>% as_data_frame(what = "edges") %>% as_tibble()
  cat(" done\n\tJoining nodes with edges...")
  g_df <- edges %>%
    left_join(nodes, by = c(from = "name")) %>%
    left_join(nodes, by = c(to = "name"))
  cat(" done\nEstimating intra community score modes...")
  score_m <- g_df %>%
    filter(com.x == com.y) %>%
    select(com.x, weight) %>%
    group_by(com.x) %>%
    summarise(mode = estimate_mode(weight)) %>% ungroup() %>% dplyr::rename(com = com.x)
  cat(" done\n")
  score_m
}

contract_graphs <- function(X, G = G) {
  msg_sub("Filtering graph...")
  G <- G %>% activate(nodes) %>% inner_join(X %>% rename(name = vertex) %>%
                                              mutate(name = as.character(name)), by = c(name = "name")) #%>% select(-class, -complete, -partial)
  cat("done.\n")
  msg_sub("Contracting vertices...")
  attr <- as.numeric(as.factor(get.vertex.attribute(G, "com")))
  gc <- contract.vertices(G, attr, vertex.attr.comb = list(name = function(x) paste(x, collapse = "|"),
                                                           com = function(x) unique(x),
                                                           archit = function(x) paste(unique(x), collapse = "__"),
                                                           original = function(x) paste(unique(x), collapse = "__")))
  cat("done.\n")
  msg_sub("Simplifying edges...")
  gc <- igraph::simplify(gc, edge.attr.comb = list(weight = function(x) {
    estimate_mode(x)
  }), remove.loops = TRUE, )
  cat("done...\n")
  msg_sub("Cleaning up...")
  rm(G)
  gc()
  cat("done.\n")
  list(graph = gc %>% as_tbl_graph(), da = gc %>% as_tbl_graph() %>% as_tibble())
}

#G <- k_hh_g_big
#X <- k_hh_gc_da_com[[1]]
contract_graphs_refined <- function(X, G = G) {
  G <- G %>% activate(nodes) %>% inner_join(X %>% select(name, com) %>% unique(), by = c(name = "name"))
  G <- G %>% activate(nodes) %>% select(name, com, archit, original)
  attr <- as.numeric(as.factor(get.vertex.attribute(G, "com")))
  gc <- contract.vertices(G, attr,
                          vertex.attr.comb = list(name = function(x) paste(x, collapse = "|"),
                                                  com = function(x) unique(x),
                                                  archit = function(x) paste(unique(x), collapse = "__"),
                                                  original = function(x) paste(unique(x), collapse = "__")))

  gc <- igraph::simplify(gc, edge.attr.comb = list(weight = function(x) {
    estimate_mode(x)
  }), remove.loops = TRUE)
  rm(G)
  gc()

  list(graph = gc %>% as_tbl_graph(), da = gc %>% as_tbl_graph() %>% as_tibble())
}

estimate_mode <- function(x) {
  if (length(x) == 2) {
    mean(x)
  } else if (length(x) == 1){
    x
  }else if (length(x) == 0){
    NA
  }else {
    d <- density(x)
    d$x[which.max(d$y)]
  }
}


# What does mean?

#X <- k_gc[[1]]
# PDDEXK|UvrD
evaluate_da_communities <- function(x, da) {
  c <- map_df(strsplit(da %>% filter(com == x) %>% .$archit %>% strsplit("__") %>%
                         unlist(recursive = F), "\\|"), function(x) {
                           if (length(x) > 1) {
                             combn(x, 2, simplify = T) %>% t %>% as_tibble()
                           } else {
                             tibble(V1 = x, V2 = x)
                           }
                         }) %>% unique() %>% graph_from_data_frame(directed = F) %>% igraph::simplify() %>% igraph::decompose()
  tibble(com = x, n_comp = length(c), comps = list(c)) %>% unique()
}


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



optimal_mcl <- function(I = I, G = G, Gx = Gx, mcl_cores = 1, tmp_dir = "/tmp") {
  options <- paste("-te", mcl_cores, "--abc -I", I, sep = " ")
  I <- as.character(I)
  run_mcl(X = "g", graph = Gx, options = options, mcl_bin = mcl_bin, tmp_dir = tmp_dir)
}


# profvis(refine_communities(X$com[[1]], g = k_hh_g_big_da, df = X))


refine_communities <- function(c, g, df) {

  # Check if clusters shoud be together in the MCL component
  # Check if clusters shoud be split in the MCL component


  # # df <- X
  # c <- 78
  # df <- k_hh_gc_da_sg[[19]]
  #c <- 5
  #  g <- k_hh_g_big_da

  Y <- df %>% filter(com == c)

  #if (Y$n_comp > 1) {

  # Get communities from the MCL cluster and add their ids
  comps <- map_df(unlist(Y$comps, recursive = F), function(x) as_tbl_graph(x) %>% mutate(degree = centrality_degree()) %>%
                    as_tibble, .id = "id")


  # Aggregate the names of each component in a list
  comps_re <- comps %>% group_by(id) %>% summarise(comps = list(name))

  # Get names of the component members
  c_names <- Y %>% select(name) %>% separate_rows(name, sep = "\\|") %>% unique()

  # Get the archits of each member
  g <- g %>%
    filter(name %in% c_names$name) %>%
    inner_join(p_doms %>% dplyr::rename(name = cl_name) %>% select(name, original)) %>%
    select(name, original) %>% dplyr::rename(o_archit = original)


  # Get the multi domain clusters
  multi_dom <- g %>%
    filter(grepl("\\|", o_archit))

  # Get the single domain clusters
  single_dom <- g %>%
    filter(!(grepl("\\|", o_archit)))

  # For multi domain clusters
  if (nrow(multi_dom) > 0) {
    # Simplify names and get component id
    multi_dom <- multi_dom %>%
      mutate(narchit = o_archit) %>%
      separate_rows(narchit, sep = "\\|") %>%
      mutate(cname = plyr::mapvalues(narchit,
                                     p_doms$original,
                                     p_doms$archit, warn_missing = FALSE)) %>% #filter(name == "25187266")
      inner_join(comps %>% rename(cname = name)) %>% ungroup() %>% #filter(name == "25187266")
      # dplyr::group_by(name) %>% #slice(which.max(degree)) %>%
      select(id, name, narchit) %>% dplyr::rename(o_archit = narchit) #%>% filter(name == "25187266")

    # Do the same as above for singletons
    nodeid <- map_df(comps_re$id, function(x) {
      single_dom %>%
        mutate(cname = plyr::mapvalues(o_archit,
                                       p_doms$original,
                                       p_doms$archit, warn_missing = FALSE)) %>%
        filter(cname %in% unlist(comps_re[x, ]))
    }, .id = "id") %>% select(-cname) %>%
      bind_rows(multi_dom) %>% arrange(name)

  }else{
    nodeid <- map_df(comps_re$id, function(x) {
      single_dom %>%
        mutate(cname = plyr::mapvalues(o_archit,
                                       p_doms$original,
                                       p_doms$archit, warn_missing = FALSE)) %>%
        filter(cname %in% unlist(comps_re[x, ]))
    }, .id = "id") %>% select(-cname)
  }

  # nodeid <- map_df(comps_re$id, function(x) {
  #   g %>% filter(archit %in% unlist(comps_re[x, ]))
  # }, .id = "id") %>% bind_rows(multi_dom)

  # Add clan to pfams
  # p_clan_term <- p_doms %>%
  #   separate_rows(archit, sep = "\\|") %>%
  #   dplyr::rename(pfam = archit) %>%
  #   left_join(p_clan) %>%
  #   select(pfam, clan) %>%
  #   unique() %>%
  #   filter(!is.na(clan))
  #

  # Add clan to the table
  nodeid <- nodeid %>%
    left_join(p_clan %>% dplyr::rename(o_archit = pfam)) %>%
    mutate(cname = plyr::mapvalues(o_archit,
                                   p_doms$original,
                                   p_doms$archit, warn_missing = FALSE))
  # Find the consensus clan for each component
  nodeid_clan <- nodeid %>%
    group_by(id) %>%
    summarise(clan_mv = ifelse(length(unique(clan)) < 2 && is.na(unique(clan)), paste0("NOCLAN_", id), majority_vote(clan)$majority)) %>%
    mutate(group_no = as.integer(factor(clan_mv)))

  # Rename communities to reflect the new communities
  nodeid_orig <- nodeid %>%
    left_join(nodeid_clan, by = c(id = "id")) %>%
    mutate(com = paste(c, group_no, sep = "_"))

  #term_group <- nodeid_orig %>% group_by(cname) %>% summarise(new_group = majority_vote(com)$majority) %>%
  #    filter(!is.na(cname))
  term_group <- nodeid_orig %>%
    group_by(com) %>%
    count() %>%
    ungroup() %>%
    left_join(nodeid_orig, by = "com") %>%
    select(cname, com, n) %>%
    filter(!is.na(cname)) %>%
    unique() %>%
    group_by(cname) %>%
    mutate(N = n()) %>%
    filter(N > 1) %>%
    top_n(n = 1, wt = n)

  to_filt <- nodeid_orig %>%
    select(-com) %>%
    inner_join(term_group) %>%
    select(name, com) %>% unique()

  if (nrow(term_group) > 0){
    nodeid_orig <- nodeid_orig %>%
      mutate(com = ifelse(name %in% to_filt$name,
                          plyr::mapvalues(name, from = to_filt$name, to = to_filt$com),
                          com)
      )

  }else{
    nodeid_orig <- nodeid_orig
  }
  nodeid_orig #%>%
  #mutate(com = ifelse(length(unique(com)) < 2 && is.na(unique(com)), com, majority_vote(com)$majority))

  # } else {
  #
  #   comps <- map_df(unlist(Y$comps, recursive = F), function(x) as_tbl_graph(x) %>%
  #                     as_tibble, .id = "id")
  #
  #   comps_re <- comps %>% group_by(id) %>% summarise(comps = list(name))
  #
  #   c_names <- Y %>% select(name) %>% separate_rows(name, sep = "\\|")
  #
  #   g <- g %>% filter(name %in% c_names$name)
  #
  #   multi_dom <- g %>% filter(grepl("\\|", archit))
  #
  #   if (nrow(multi_dom) > 0) {
  #     multi_dom <- multi_dom %>%
  #       mutate(narchit = archit) %>%
  #       separate_rows(narchit, sep = "\\|") %>%
  #       inner_join(comps %>% rename(narchit = name)) %>% ungroup() %>%
  #       dplyr::group_by(name) %>% slice(which.max(degree)) %>%
  #       select(id, name, narchit) %>% rename(archit = narchit)
  #
  #     nodeid <- map_df(comps_re$id, function(x) {
  #       g %>% filter(archit %in% unlist(comps_re[x, ]))
  #     }, .id = "id") %>% bind_rows(multi_dom)
  #   }else{
  #     nodeid <- map_df(comps_re$id, function(x) {
  #       g %>% filter(archit %in% unlist(comps_re[x, ]))
  #     }, .id = "id")
  #   }
  #
  #
  #   nodeid <- nodeid %>% dplyr::rename(pfam = archit) %>% mutate(clan = plyr::mapvalues(pfam,
  #                                                                                       p_clan$pfam,
  #                                                                                       p_clan$clan, warn_missing = FALSE),
  #                                                                cname = plyr::mapvalues(pfam,
  #                                                                                        p_term$pfam,
  #                                                                                        p_term$cname, warn_missing = FALSE))
  #   # get clan information
  #
  #   nodeid_clan <- nodeid %>%
  #     group_by(id) %>%
  #     summarise(clan_mv = ifelse(length(unique(clan)) < 2 && is.na(unique(clan)), paste0("NOCLAN_", id), majority_vote(clan)$majority)) %>%
  #     mutate(group_no = as.integer(factor(clan_mv)))
  #
  #   nodeid_orig <- nodeid %>% left_join(nodeid_clan, by = c(id = "id")) %>%
  #     mutate(com = paste(c, 1, sep = "_")) %>% select(-cname) %>%
  #     left_join(p_term, by = "pfam")
  #
  #   term_group <- nodeid_orig %>%
  #     group_by(com) %>%
  #     count() %>%
  #     left_join(nodeid_orig, by = "com") %>%
  #     select(cname, com, n) %>%
  #     filter(!is.na(cname)) %>%
  #     unique() %>%
  #     group_by(cname) %>%
  #     mutate(N = n()) %>%
  #     filter(N > 1) %>%
  #     top_n(n = 1, wt = n)
  #
  #   if (nrow(term_group) > 0){
  #     nodeid_orig %>%
  #       mutate(com = ifelse(cname %in% term_group$cname,
  #                           plyr::mapvalues(cname,
  #                                           from = term_group$cname,
  #                                           to = term_group$com), com
  #       )
  #       )
  #   }else{
  #     nodeid_orig
  #   }
  #
  #
  # }
}

evaluate_component_entropy <- function(c, g, df) {

  # Check if clusters shoud be together in the MCL component
  # Check if clusters shoud be split in the MCL component


  # # df <- X
  # c <- 78
  # df <- k_hh_gc_da_sg[[1]]
  #c <- 5
  #  g <- k_hh_g_da

  Y <- df %>% filter(com == c)

  #if (Y$n_comp > 1) {

  # Get communities from the MCL cluster and add their ids
  comps <- map_df(unlist(Y$comps, recursive = F), function(x) as_tbl_graph(x) %>% mutate(degree = centrality_degree()) %>%
                    as_tibble, .id = "id")


  # Aggregate the names of each component in a list
  comps_re <- comps %>% group_by(id) %>% summarise(comps = list(name))

  # Get names of the component members
  c_names <- Y %>% select(name) %>% separate_rows(name, sep = "\\|") %>% unique()

  # Get the archits of each member
  g <- g %>%
    dt_filter(name %in% c_names$name) %>%
    #dt_inner_join(p_doms %>% dplyr::rename(name = cl_name) %>% select(name, original)) %>%
    dt_select(name, original) %>% dplyr::rename(o_archit = original)


  # Get the multi domain clusters
  multi_dom <- g %>%
    dt_filter(grepl("\\|", o_archit))

  # Get the single domain clusters
  single_dom <- g %>%
    dt_filter(!(grepl("\\|", o_archit)))

  p_doms_m <- p_doms %>% dt_filter(original %in% (multi_dom$o_archit %>% unique()))
  p_doms_s <- p_doms %>% dt_filter(original %in% (single_dom$o_archit %>% unique()))
  # For multi domain clusters
  if (nrow(multi_dom) > 0) {
    # Simplify names and get component id
    multi_dom <- multi_dom %>%
      mutate(narchit = o_archit) %>%
      separate_rows(narchit, sep = "\\|") %>%
      mutate(cname = plyr::mapvalues(narchit,
                                     p_doms_m$original,
                                     p_doms_m$archit, warn_missing = FALSE)) %>% #filter(name == "25187266")
      inner_join(comps %>% rename(cname = name)) %>%
      ungroup() %>% #filter(name == "25187266")
      # dplyr::group_by(name) %>% #slice(which.max(degree)) %>%
      select(id, name, narchit) %>%
      dplyr::rename(o_archit = narchit) #%>% filter(name == "25187266")

    # Do the same as above for singletons
    nodeid <- map_df(comps_re$id, function(x) {
      single_dom %>%
        mutate(cname = plyr::mapvalues(o_archit,
                                       p_doms_s$original,
                                       p_doms_s$archit, warn_missing = FALSE)) %>%
        filter(cname %in% unlist(comps_re[x, ]))
    }, .id = "id") %>% select(-cname) %>%
      bind_rows(multi_dom)

  }else{
    nodeid <- map_df(comps_re$id, function(x) {
      single_dom %>%
        mutate(cname = plyr::mapvalues(o_archit,
                                       p_doms_s$original,
                                       p_doms_s$archit, warn_missing = FALSE)) %>%
        filter(cname %in% unlist(comps_re[x, ]))
    }, .id = "id") %>% select(-cname)
  }

  # nodeid <- map_df(comps_re$id, function(x) {
  #   g %>% filter(archit %in% unlist(comps_re[x, ]))
  # }, .id = "id") %>% bind_rows(multi_dom)

  # Add clan to pfams
  # p_clan_term <- p_doms %>%
  #   separate_rows(archit, sep = "\\|") %>%
  #   dplyr::rename(pfam = archit) %>%
  #   left_join(p_clan) %>%
  #   select(pfam, clan) %>%
  #   unique() %>%
  #   filter(!is.na(clan))
  #

  # Add clan to the table
  p_doms_a <- p_doms %>% dt_filter(original %in% (nodeid$o_archit %>% unique()))
  nodeid <- nodeid %>%
    left_join(p_clan %>% dplyr::rename(o_archit = pfam)) %>%
    mutate(cname = plyr::mapvalues(o_archit,
                                   p_doms$original,
                                   p_doms$archit, warn_missing = FALSE))
  # Find the consensus clan for each component
  cl <- nodeid %>%
    group_by(id) %>%
    summarise(clan_mv = ifelse(length(unique(clan)) < 2 && is.na(unique(clan)), paste0("NOCLAN_", id), majority_vote(clan)$majority))
  freqs <- table(cl$clan_mv)/length(cl$clan_mv)

  tibble(com = c, entropy = -sum(freqs * log2(freqs)), n_comp = Y$n_comp)
}


# Progress bar ------------------------------------------------------------
