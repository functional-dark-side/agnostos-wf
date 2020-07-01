## AgnostosDB dataset structure

#### agnostosDB/

-   cluster_ids_categ.tsv.gz - Table with cluster ids and categories.
-   cluster_ids_categ_genes.tsv.gz - Table with cluster ids and categories and genes (ORFs-headers).
-   cluster_db_size_categ_origin.tsv.gz - Table with cluster ids the sequence original database/provenience, the cluster size and the cluster categories (ORFs-headers).
-   cluDB_name_origin_size.tsv.gz - Table with cluster ids, the database/provenience of origin, and the cluster size (number of orfs)
-   cluster_category_summary_stats.tsv - A summary table with the refined cluster statistics. It includes cluster size, gene length stats, cluster gene completeness stats, HQ clusters, cluster level of darkness and disorder, and taxonomic info (prevalent taxonomy and taxonomic entropy).
-   orf_partial_info.tsv.gz - list of gene headers and their level of completness, based on the Prodigal prediction results.
-   HQ_clusters.tsv.gz - set of high quality (HQ) clusters (clusters with high percentage of complete ORFs). Fields: cluster name, category.
-   cluster_communities.tsv.gz - Table containing the correspondence between cluster and cluster communitites. Fields: cluster name, community name, category.
-   pfam_name_acc_clan_multi.tsv.gz - All genes Pfam multi-domain annotations. Fields: orf, Pfam_name, Pfam_accession, Pfam_clan.
-   cluster_pfam_domain_architectures.tsv.gz - Cluster consensus Pfam domain architectures.
-   K/KWP/GU_annotations.tsv.gz - Clusters annotations (Pfam, Uniref90, NCBI nr and Uniclust).
-   spurious_shadow_info.tsv.gz - Table containing information about spurious and shadow ORFs for each cluster. Fields: orf, length, cl_name, size, proportion of shadow per cluster, is.shadow, is.spurious.

<br>

-   **agnostosDB_dbf02445-20200519_mmseqs_clustering/**
        -   seqDB (seqDB.dbtype, seqDB.index, seqDB.lookup, seqDB_h, seqDB_h.index) - MG+GTDB MMseqs sequence database.
        -   cluDB (cluDB.index) -  MG+GTDB MMseqs cluster database.
        -   cluDB_name_rep_size.tsv - MG+GTDB cluster name, representative and size.
-   **agnostosDB_dbf02445-20200519_mmseqs_profiles/**
        -   clu_hmm_db (clu_hmm_db.index, clu_hmm_db.dbtype, clu_hmm_db_h, clu_hmm_db_h.dbtype, clu_hmm_db_h.index) - Cluster HMMs profile database in MMseqs format. This is the target database to use to perform profile searches against the MG+GTDB clusterDB.
-   **agnostosDB_dbf02445-20200519_mmseqs_cluseqdb/**
        -   clu_seqDB (clu_seqDB.index, clu_seqDB.dbtype) - Refined cluster sequence database, MMseqs format.
-   **agnostosDB_dbf02445-20200519_hh-suite-db/**
        -  Cluster database in hh-suite format (check hh-suite documentation for more information).
        -   clu_aln - (cluster multiple sequence alignments (MSAs), retrieved with FAMSA)
        -   clu_consensus - cluster consensus sequences
        -   clu_a3m - reformatted cluster MSAs in A3M format
        -   clu_hmm - cluster HHM-formatted HMMs profiles
        -   clu_cs219 -  cluster column-state sequences for prefiltering

<br>
Additional metadata (not required for the workflow)

-   **agnostosDB_dbf02445-20200519_original-data/**
    -   project]_metag_contig_genes.tsv.gz - metagenomic contigs and predicted genes per metagenome.
    -   GTDB_geneome_contig_genes.tsv.gz - genomic contigs and predicted genes per genome.
-   **agnostosDB_dbf02445-20200519_environmental/**
    -   mg_cluster_sample_norfs_coverage.tsv.gz - Metagenomic cluster contextual data table, with cluster name, sample_ID, number of ORFs and coverage for each sample.
    -   niche_breadth/
        -   gCL_nb_all.Rda, gClCo_nb_all.Rda - R object containing two lists with all niche breadth values.
        -   gCl_nb_all_mv.Rda, gClCo_nb_all_mv.Rda - R object containing two tables with the mean values.
-   **agnostosDB_dbf02445-20200519_phylogenetic/**
    -   gtdb_cluster_genome_norfs.tsv.gz - Genomic clusters contextual data table, with cluster name, genome accession (GTDB) and number of ORFs
    -   mg_gtdb_lineage_specific_clusters.tsv.gz - Lineage specific clusters.
    -   mg_gtdb_lineage_specific_communities.tsv.gz - Lineage specific clusters communitites.

#### AgnostosDB numbers overview

<img alt="mg_gtdb_numbers.png" src="assets/mg_gtdb_numbers.png" width="80%" height="" >
