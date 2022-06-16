rule cluster_db_results:
    input:
        iclu_cat = config["rdir"] + "/cluster_categories/cluster_ids_categ.tsv",
        iclu_com = config["rdir"] + "/cluster_communities/cluster_communities.tsv",
        iclu_stats = config["rdir"] + "/cluster_category_stats/cluster_category_summary_stats.tsv",
        ihq_clu = config["rdir"] + "/cluster_category_stats/HQ_clusters.tsv",
        iclu_hhm = config["rdir"] + "/cluster_category_DB/clu_hhm_db"
    threads: 7
    conda:
        config["conda_env"]
    params:
        mmseqs_bin = config["mmseqs_bin"],
        data_name = config["data_name"],
        singl = config["singl"],
        parser = config["wdir"] +  "/scripts/creation_output_tables.r",
        rdir = config["rdir"] + "/clusterDB_results",
        or_clu_gene = config["rdir"] + "/cluster_categories/cluster_ids_categ_genes.tsv.gz",
        s_categ = config["rdir"] + "/cluster_classification/singleton_gene_cl_categories.tsv",
        clu_info = config["rdir"] + "/mmseqs_clustering/cluDB_info.tsv",
        or_clu_origin = config["rdir"] + "/mmseqs_clustering/cluDB_name_rep_size.tsv",
        or_sp_sh = config["rdir"] + "/spurious_shadow/spurious_shadow_info.tsv",
        or_multi_annot = config["rdir"] + "/annot_and_clust/pfam_name_acc_clan_multi.tsv.gz",
        or_partial = config["rdir"] + "/gene_prediction/orf_partial_info.tsv",
        sp_sh = config["rdir"] + "/clusterDB_results/spurious_shadow_info.tsv.gz",
        multi_annot = config["rdir"] + "/clusterDB_results/pfam_name_acc_clan_multi.tsv.gz",
        partial = config["rdir"] + "/clusterDB_results/orf_partial_info.tsv.gz",
        clu_gene = config["rdir"] + "/clusterDB_results/cluster_ids_categ_genes.tsv.gz",
        clu_origin = config["rdir"] + "/clusterDB_results/cluDB_name_origin_size.tsv",
        clu_hhm = config["rdir"] + "/clusterDB_results/mmseqs-profiles/",
        clu_seq = config["rdir"] + "/clusterDB_results/mmseqs-cluseqdb/clu_seqDB",
        singl_cl_gene_categ = config["rdir"] + "/clusterDB_results/singleton_gene_cl_categories.tsv.gz"
    log:
        out = "logs/cludb_stdout.log",
        err = "logs/cludb_stderr.err"
    benchmark:
        "benchmarks/clu_db_res.tsv"
    output:
        clu_cat = config["rdir"] + "/clusterDB_results/cluster_ids_categ.tsv.gz",
        clu_com = config["rdir"] + "/clusterDB_results/cluster_communities.tsv.gz",
        clu_stats = config["rdir"] + "/clusterDB_results/cluster_category_summary_stats.tsv.gz",
        clu_out_tbl = config["rdir"] + "/clusterDB_results/DB_genes_summary_info.tsv",
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
        if [ ! -f {params.multi_annot} ]; then
            mv {params.or_multi_annot} {params.multi_annot}
        fi

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

        # New integarted cluster HHMs DB and mmseqs profiles
        mkdir -p {params.clu_hhm}
        cp -r {input.iclu_hhm}* {params.clu_hhm}/
        {params.mmseqs_bin} convertprofiledb {input.iclu_hhm} {params.clu_hhm}/clu_hmm_db \
                                             --threads {threads} 2>{log.err}

        # Integrated set of high quality (HQ) clusters
        gzip -c {input.ihq_clu} > {output.hq_clu}

        {params.parser} --clu_or {params.clu_origin} \
                          --cat {output.clu_cat} \
                          --clu_info {params.clu_info} \
                          --comm {output.clu_com} \
                          --hq_clu {output.hq_clu} \
                          --k_annot {params.multi_annot} \
                          --is_singl {params.singl} \
                          --s_cat {params.s_categ} \
                          --res {output.clu_out_tbl} \
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
        clu_out_tbl = config["rdir"] + "/clusterDB_results/DB_genes_summary_info.tsv",
        hq_clu = config["rdir"] + "/clusterDB_results/HQ_clusters.tsv.gz"
    output:
        cludb_done = touch(config["rdir"] + "/clusterDB_results/cludb_res.done")
    run:
        shell("echo 'INTEGRATED DB DONE'")
