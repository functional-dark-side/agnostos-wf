## Workflow output data structure

### Creation-mode

The results of the creation-mode workflow are structured in folders corresponding to the workflow steps.

-   gene_prediction/
    -   orf_seqs.fasta - gene prediction sequence amino acid in fasta format.
    -   orf_partial_info.tsv - gene prediction sequence headers and Prodigal partiality information.
-   pfam_annotation/
    -   pfam_annot_parsed.tsv - parsed and filtered hmmsearch output in long format.
-   mmseqs_clustering/
    -   seqDB - MMseqs2 sequence database.
    -   cluDB - MMseqs2 cluster database.
    -   clu_seqDB - MMseqs2 cluster sequence database.
    -   cluDB.tsv - clusters in long format (representative - memeber "newline" member).
    -   cluDB_name_index.txt - cluster index of MMseqs2-ids and cluster-names.
    -   cluDB_name_rep_size.tsv - cluster names - representastive - cluster size.
    -   cluDB_info.tsv - summary table of the clustering process: cluster name - representative - memeber - length - cluster size.
    -   cluDB_no_singletons.tsv - clusters with more that one gene.
    -   cluDB_singletons.tsv - clusters with only one gene.
-   spurious_shadow/
    -   spurious_shadow_info.tsv - table containing information about quality of all genes.
-   annot_and_clust/
    -   pfam_name_acc_clan_multi.tsv.gz - all Pfam annotatated genes. Fields: gene - pfam_name - pfam_acc - pfam_clan. The multi-domain annotations are separated by a pipe "|" ("domainA|domainB").
    -   annotated_clusters.tsv - table with Pfam annotated clusters (long format).
    -   not_annotated_clusters.tsv - table with not annotated clusters (long format).
    -   singletons_pfam_annot.tsv - Pfam annotated singletons.
-   validation/
    -   functional_val_results.tsv - general/parsed cluster functional validation results.
    -   compositional_validation_filt_stats.tsv - cluster sequence composition stats.
    -   compositional_validation_rejected_orfs.tsv - list of non-homologous genes.
    -   compositional_validation_results.tsv - general/parsed cluster compositional validation results.
    The following folders contain the results of the compositional validation for each clusters, separated in folders named accordingly to the cluster names:
    -   alignments/
    -   SSN/
    -   validation_plots_for_R.rda
    -   good_clusters.tsv - set of good clusters main information.
    -   validation_results.tsv - validation results in tab separated format.
    -   validation_results_stats.tsv - validation results main cluster and gene statistics.
-   cluster_refinement/
    -   cluster_orfs_to_remove.tsv - list of genes x cluster to remove from the good clusters.
    -   refined_clusters.tsv (and refined_clusterDB) - refined cluster database.
    -   refined_annotated_clusters.tsv (and refined_annotated_clusterDB) - annotated clusters subset.
    -   refined_not_annotated_clusters.tsv (and refined_not_annotated_clusterDB) - not annotated clusters subset.
-   cluster_classification/
    -   noannot_vs_uniref90.tsv - MMseqs2 search vs UniRef90 results (blast tab format).
    -   uniref-nohits_vs_NR.tsv - MMseqs2 search vs NCBI nr results (blast tab format).
    -   cluster_pfam_domain_architectures.tsv - tabble with cluster consensus Pfam domain architectures.
    -   k/kwp/gu/eu_ids.txt - lists of cluster ids/names of the different categories (pre-refinement).
    -   k/kwp/gu_annotations.tsv - cluster annotations (pre-refinement).
-   cluster_categories/
    -   refined_clusterDB
    -   k/kwp/gu/eu_ids.txt - lists of cluster ids/names of the different categories (refined).
    -   K/KWP/GU_annotations.tsv - cluster annotations (refined).
    -   cluster_ids_categ.tsv - summary table of refined cluster names and categories.
    -   cluster_ids_categ_orfs.tsv.gz - summary table of refined cluster names, categories and genes.
    -   eu_hhbl_parsed.tsv - EU refinement hhblits results.
    -   eu_hhbl_new_gu_ids.txt
    -   eu_hhbl_new_kwp_ids.txt
    -   kwp_hhbl_name_acc_clan_multi.tsv - KWP refinement hhblits results.
    -   kwp_hhbl_new_gu_ids_annot.tsv
    -   kwp_hhbl_new_k_ids_annot.tsv
-   cluster_category_DB/
    -   k/kwp/gu/eu_orfs.txt - category gene headers.
    -   k/kwp/gu/eu_clu_orfs.fasta - category genes amino acid sequences.
    -   k/kwp/gu/eu_clseqdb.index - category cluster sequence MMseqs2 database.
    -   clu_hhm_db - cluster HMM profiles HH-suite database
    The following databases are the HH-SUITE DBs of the different cluster categories:
    -   k/kwp/gu/eu_aln
    -   k/kwp/gu/eu_a3m_db
    -   k/kwp/gu/eu_cons
    -   k/kwp/gu/eu_cs219.ffdata
    -   k/kwp/gu/eu_hhm_db.index
-   cluster_category_stats/
    -   cluster_kaiju_taxonomy.tsv
    -   cluster_category_dpd_perc.tsv - category level of darkness and disorder.
    -   cluster_dpd_perc.tsv - cluster level of darkness and disorder.
    -   cluster_category_completeness.tsv - percentage of complete gene x cluster.
    -   HQ_clusters.tsv - set of high quality clusters (clusters with high percentage of complete genes).
    -   cluster_category_summary_stats.tsv - summary table containing various information about the clusters.
    -   only_category_summary_stats.tsv - summary table containing various information about the cluster categories.
-   cluster_communities/
    -   k/kwp/gu/eu_hhblits.tsv - all-vs-all category hhblits raw results.
    -   YYYY-MM-DD-XXXXXX/ - folder containing all the community inference result files.
    -   cluster_communities.tsv - summary table containing the correspondence cluster-community.
-   report/
    -   workflow_report.html

##### Folder containing the files necessary for the DB_update module:
-   clusterDB_results/
    -   cluDB_name_origin_size.tsv - table with cluster names, their origin and their size.
    -   cluster_ids_categ.tsv - table with refined cluster names and categories.
    -   cluster_ids_categ_genes.tsv.gz  - table with refined cluster names, categories and genes.
    -   cluster_communities.tsv - summary table containing the correspondence cluster-community.
    -   cluster_category_summary_stats.tsv - summary table containing various information about the clusters.
    -   pfam_name_acc_clan_multi.tsv.gz - all genes Pfam annotations.
    -   K/KWP/GU/EU_annotations.tsv.gz - category annotations.
    -   orf_partial_info.tsv.gz - list of gene headers and their level of completness, based on the Prodigal prediction results.
    -   HQ_clusters.tsv - set of high quality clusters (clusters with high percentage of complete genes).
    -   spurious_shadow_info.tsv.gz - summary table with gene quality information.
    -   mmseqs_profiles/
            - clu_hhm_db - cluster HHM profiles in HH-suite format
            - clu_hmm_db - cluster HMM profiles MMseqs2 database (for profile searches).

    Plus, the mmseqs_clustering/ folder has to be copied here as well, including the cluDB, the seqDB and the file cluDB_name_rep_size.tsv

##### General summary table for the gene cluster DB
    -  clusterDB_results/DB_genes_summary_info.tsv
      -  gene_callers_id
      -  cl_name
      -  cl_size
      -  category
      -  is.HQ
      -  community
      -  pfam

### Update-mode

The results structure is the same as the creation-mode one (the clusters processed through the workflow steps are those not found in the original DB), plus a folder containing the cluster-update summary results (derived from the merging of the new with the original clusterDB):

-   integrated_cluster_DB/
    -   cluDB_name_origin_size.tsv - table with cluster names, their origin and their size.
    -   cluster_ids_categ.tsv - table with refined cluster names and categories.
    -   cluster_ids_categ_genes.tsv.gz  - table with refined cluster names, categories and genes.
    -   cluster_communities.tsv - summary table containing the correspondence cluster-community.
    -   cluster_category_summary_stats.tsv - summary table containing various information about the clusters.
    -   singletons_gene_cl_categories.tsv - table with singleton genes, cluster names and categories (if the configurtion "singl" was set to "true").
    -   pfam_name_acc_clan_multi.tsv.gz - all genes Pfam annotations.
    -   K/KWP/GU/EU_annotations.tsv.gz - category annotations.
    -   orf_partial_info.tsv.gz - list of gene headers and their level of completness, based on the Prodigal prediction results.
    -   HQ_clusters.tsv - set of high quality clusters (clusters with high percentage of complete genes).
    -   spurious_shadow_info.tsv.gz - summary table with gene quality information.
    -   mmseqs_profiles/
        -   clu_hmm_db - cluster HMM profiles MMseqs2 database (for profile searches).


Tables summarising the results and eventually the cluster contextual data

-   output_tables/
    -   contig_genes.tsv - genome - contig - gene ids
    -   DB_cluster_annotations.tsv - summary of K, KWP and GU annotations per cluster
    -   DB_genes_clusters_communities.tsv - gene - cluster - category - community
    -   DB_genes_summary_info_red.tsv - reduced set of information per cluster
    -   DB_genes_summary_info_exp.tsv - expanded set of information per cluster
    Contextual data (if the pre-existing DB is or originates from the agnostosDB)
        Fields:
        -  gene_callers_ids - gene identifier
        -  cl_name - cluster identifier
        -  contig - contig identifier
        -  gene_x_contig - numner of genes per contig
        -  db - database of origin (agnostosDB, new, etc..)
        -  cl_size - number of genes per cluster
        -  category - AGNOSTOS GC category
        -  pfam - GC pfam domain annotation, in the form of domain architectures (domains separated by "|")
        -  is.HQ - Logic (TRUE if the GC is high-quality)
        -  is.LS - Logic (TRUE if the GC is lineage-specific)
        -  lowest_rank - lineage-secific taxonomic rank
        -  lowest_level - lineage-secific taxonomic level
        -  niche_breadth_sign - Levin's niche breadth distribution

    -   DB_clusters_niche_breadth.tsv - clusters with significant niche breadth values in metagenomes
    -   DB_lineage_specific_clusters.tsv - linege-specific clusters within the GTDB phylogeney
    -   DB_mutant_phenotype_clusters.tsv - clusters with mutant phenotype (Proce et al. 2018)
    -   DB_clusters_in_metagenomes.tsv
    -   DB_clusters_in_gtdb_genomes.tsv

The cluster-update results in the form of MMseqs2 databases are stored in the "mmseqs_clustering/" folder.

### Profile-search

The profile search output consist in one main file containing the search results, and two additional files which are generated only if the "gene-info" file, containing the gene-to-contig correspondance is specified in input.

The output files:

-   "your_name_search_res_best-hits.tsv": best-hits with categories (gene_callers_id-cl_name-category-evalue).

-   "your_name_search_res_summary-categ.tsv": proportion of different categories per contig.

-   "your_name_search_res_summary-classes.tsv": proportion of classes per contig, where the classes are defined grouping the categories into "unknown" (EUs and GUs) and "known" (Ks and KWPs)
