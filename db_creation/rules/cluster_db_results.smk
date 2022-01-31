rule cluster_db_results:
    input:
        iclu_cat = config["rdir"] + "/cluster_categories/cluster_ids_categ.tsv",
        iclu_com = config["rdir"] + "/cluster_communities/cluster_communities.tsv",
        iclu_stats = config["rdir"] + "/cluster_category_stats/cluster_category_summary_stats.tsv",
        ihq_clu = config["rdir"] + "/cluster_category_stats/HQ_clusters.tsv",
        iclu_hmm = config["rdir"] + "/cluster_category_DB/clu_hmm_db"
    threads: 7
    params:
        rdir = config["rdir"] + "/clusterDB_results",
        parser = "/scripts/output_tables.r",
        data_name = config["data_name"],
        or_clu_seq = config["rdir"] + "/cluster_categories/refined_clusterDB",
        or_clu_gene = config["rdir"] + "/cluster_categories/cluster_ids_categ_genes.tsv.gz",
        singl = config["singl"],
        s_categ = config["rdir"] + "/cluster_classification/singleton_gene_cl_categories.tsv",
        clu_info = config["rdir"] + "/mmseqs_clustering/cluDB_info.tsv",
        or_clu_origin = config["rdir"] + "/mmseqs_clustering/cluDB_name_rep_size.tsv",
        or_sp_sh = config["rdir"] + "/spurious_shadow/spurious_shadow_info.tsv",
        or_multi_annot = config["rdir"] + "/annot_and_clust/pfam_name_acc_clan_multi.tsv",
        or_partial = config["rdir"] + "/combine_samples/orf_partial_info.tsv",
        sp_sh = config["rdir"] + "/clusterDB_results/spurious_shadow_info.tsv.gz",
        multi_annot = config["rdir"] + "/clusterDB_results/pfam_name_acc_clan_multi.tsv.gz",
        partial = config["rdir"] + "/clusterDB_results/orf_partial_info.tsv.gz",
        clu_gene = config["rdir"] + "/clusterDB_results/cluster_ids_categ_genes.tsv.gz",
        clu_origin = config["rdir"] + "/clusterDB_results/cluDB_name_origin_size.tsv",
        clu_hmm = config["rdir"] + "/clusterDB_results/mmseqs-profiles/",
        clu_seq = config["rdir"] + "/clusterDB_results/mmseqs-cluseqdb/clu_seqDB",
        singl_cl_gene_categ = config["rdir"] + "/clusterDB_results/singleton_cl_ids_categ_genes.tsv.gz",
        clu_out_tbl = config["rdir"] + "/clusterDB_results/DB_genes_summary_info.tsv"
    log:
        out = "logs/cludb_stdout.log",
        err = "logs/cludb_stderr.err"
    benchmark:
        "benchmarks/clusterDB_results/cludb.tsv"
    output:
        clu_cat = config["rdir"] + "/clusterDB_results/cluster_ids_categ.tsv.gz",
        clu_com = config["rdir"] + "/clusterDB_results/cluster_communities.tsv.gz",
        clu_stats = config["rdir"] + "/clusterDB_results/cluster_category_summary_stats.tsv.gz",
        hq_clu = config["rdir"] + "/clusterDB_results/HQ_clusters.tsv.gz"
    shell:
        """

        set -x
        set -e

        # Summary table with cluster db origin (original/shared/new)
        awk -v N={params.data_name} '{{print $1,N,$3}}' {params.or_clu_origin} > {params.clu_origin}
        sed -i 's/ /\\t/g' {params.clu_origin}

        # All gene headers and partiality information
        gzip -c {params.or_partial}  > {params.partial}

        # Spurious and shadow genes information:
        gzip -c {params.or_sp_sh} > {params.sp_sh}

        # All gene Pfam annotations:
        gzip -c {params.or_multi_annot} > {params.multi_annot}

        # All cluster category annotation files
        ODIR=$(dirname {input.iclu_cat})

        # Combine with new ones
        gzip -c ${{ODIR}}/K_annotations.tsv > {params.rdir}/K_annotations.tsv.gz
        gzip -c ${{ODIR}}/KWP_annotations.tsv > {params.rdir}/KWP_annotations.tsv.gz
        gzip -c ${{ODIR}}/GU_annotations.tsv > {params.rdir}/GU_annotations.tsv.gz

        # Integrated set of cluster categories
        gzip -c {input.iclu_cat} > {output.clu_cat}

        # and the cluster genes
        if [ {params.singl} == true ]; then
            awk -vOFS="\\t" '{{print $2,$3,$1}}' {params.s_categ} | gzip -c > {params.singl_cl_gene_categ}
        fi
        cp {params.or_clu_gene} {params.clu_gene}

        # Integrated set of cluster communities
        gzip -c {input.iclu_com} > {output.clu_com}

        # New integarted cluster HMMs DB (for MMseqs profile searches)
        mkdir -p {params.clu_hmm}
        cp -r {input.iclu_hmm}* {params.clu_hmm}/

        # New integarted cluster sequence DB (MMseqs sequence database with clustered index)
        seqdir=$(dirname {params.clu_seq})
        mkdir -p ${{seqdir}}
        cp {params.or_clu_seq} {params.clu_seq}
        cp {params.or_clu_seq}.index {params.clu_seq}.index
        cp {params.or_clu_seq}.dbtype {params.clu_seq}.dbtype

        # Integrated set of high quality (HQ) clusters
        gzip -c {input.ihq_clu} > {output.hq_clu}

        ./{params.parser} --clu_or {params.clu_origin} \
                          --cat {output.clu_cat} \
                          --clu_info {params.clu_info} \
                          --comm {output.clu_com} \
                          --hq_clu {output.hq_clu} \
                          --k_annot {params.multi_annot} \
                          --is_singl {params.singl} \
                          --s_cat {params.s_categ} \
                          --res {params.clu_out_tbl} \
                          --threads {threads} 2>{log.err} 1>{log.out}

        gzip {params.clu_origin}

        # Integrated cluster summary information
        gzip -c {input.iclu_stats}  > {output.clu_stats}

        """

rule cludb_res_done:
    input:
        clu_cat = config["rdir"] + "/clusterDB_results/cluster_ids_categ.tsv.gz",
        clu_com = config["rdir"] + "/clusterDB_results/cluster_communities.tsv.gz",
        clu_stats = config["rdir"] + "/clusterDB_results/cluster_category_summary_stats.tsv.gz",
        hq_clu = config["rdir"] + "/clusterDB_results/HQ_clusters.tsv.gz"
    output:
        cludb_done = touch(config["rdir"] + "/clusterDB_results/cludb_res.done")
    run:
        shell("echo 'DB DONE'")
