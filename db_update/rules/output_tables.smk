rule output_tables:
    input:
        genes = config["rdir"] + "/gene_prediction/orf_seqs.txt",
        cat = config["rdir"] + "/integrated_cluster_DB/cluster_ids_categ.tsv.gz"
    threads: 28
    conda:
        config["conda_env"]
    params:
        contig = config["rdir"] + "/output_tables/contig_genes.tsv",
        new_data_name = config["new_data_name"],
        orig_db = config["ordir"],
        clu_info = config["rdir"] + "/mmseqs_clustering/cluDB_info.tsv",
        clu_origin = config["rdir"] + "/integrated_cluster_DB/cluDB_name_origin_size.tsv.gz",
        comm = config["rdir"] + "/integrated_cluster_DB/cluster_communities.tsv.gz",
        hq_clu = config["rdir"] + "/integrated_cluster_DB/HQ_clusters.tsv.gz",
        k_annot = config["rdir"] + "/integrated_cluster_DB/K_annotations.tsv.gz",
        kwp_annot = config["rdir"] + "/integrated_cluster_DB/KWP_annotations.tsv.gz",
        gu_annot = config["rdir"] + "/integrated_cluster_DB/GU_annotations.tsv.gz",
        singl = config["singl"],
        singl_cat = config["rdir"] + "/integrated_cluster_DB/singleton_cl_ids_categ_genes.tsv.gz",
        parser = "/scripts/output_tables.r"
    output:
        res = config["rdir"] + "/output_tables/DB_genes_summary_info_exp.tsv",
        annot = config["rdir"] + "/output_tables/DB_cluster_annotations.tsv"
    log:
        out = "logs/tbl_stdout.log",
        err = "logs/tbl_stderr.err"
    benchmark:
        "benchmarks/out_tbl.tsv"
    shell:
        """

        set -x
        set -e

        DB=$(basename {params.orig_db})

        if [[ ${{DB}} == "agnostosDB" ]]; then
            # Download agnostosDB ecological analysis data
            wget https://ndownloader.figshare.com/files/23066879 -O ecological.tar.gz
            tar -xzvf ecological.tar.gz

            # Download agnostosDB phylogentic analysis data
            wget https://ndownloader.figshare.com/files/23066864 -O phylogenetic.tar.gz
            tar -xzvf phylogenetic.tar.gz

            # Download agnostosDB experimental analysis data
            wget https://ndownloader.figshare.com/files/23066864 -O phylogenetic.tar.gz
            tar -xzvf experimental.tar.gz
        fi

        # Run R script to retrieve a general table containing a summary of the integration
        # add contextual data if the original DB is the agnostosDB
        awk -vOFS='\\t' '{{split($1,a,"_\\\+|_-"); print a[1],$1}}' {input.genes} > {params.contig}

        ./{params.parser} --clu_or {params.clu_origin} \
                          --contig {params.contig} \
                          --cat {input.cat} \
                          --clu_info {params.clu_info} \
                          --name {params.new_data_name} \
                          --comm {params.comm} \
                          --hq_clu {params.hq_clu} \
                          --k_annot {params.k_annot} \
                          --kwp_annot {params.kwp_annot} \
                          --gu_annot {params.gu_annot} \
                          --orig_db {params.orig_db} \
                          --is_singl {params.singl} \
                          --s_cat {params.singl_cat} \
                          --threads {threads} 2>{log.err} 1>{log.out}

        """

rule output_tables_done:
    input:
        res = config["rdir"] + "/output_tables/DB_genes_summary_info_exp.tsv",
        annot = config["rdir"] + "/output_tables/DB_cluster_annotations.tsv"
    output:
        res_done = touch(config["rdir"] + "/output_tables/res.done")
    run:
        shell("echo 'Parsing of results DONE'")
