rule output_tables:
    input:
        genes = config["rdir"] + "gene_prediction/orf_seqs.txt",
        clu_origin = config["rdir"] + "integrated_cluster_DB/cluDB_name_origin_size.tsv.gz",
        cat = config["rdir"] + "integrated_cluster_DB/cluster_ids_categ.tsv.gz"
    threads: 28
    conda:
        config["conda_env"]
    params:
        contig = config["rdir"] + "/output_tables/contig_genes.tsv",
        new_data_name = config["new_data_name"],
        clu_info = config["rdir"] + "mmseqs_clustering/cluDB_info.tsv",
        comm = config["rdir"] + "integrated_cluster_DB/cluster_communities.tsv.gz",
        hq_clu = config["rdir"] + "integrated_cluster_DB/HQ_clusters.tsv.gz",
        k_annot = config["rdir"] + "integrated_cluster_DB/K_annotations.tsv.gz",
        kwp_annot = config["rdir"] + "integrated_cluster_DB/KWP_annotations.tsv.gz",
        gu_annot = config["rdir"] + "integrated_cluster_DB/GU_annotations.tsv.gz",
        parser = "/scripts/output_tables.r"
    output:
        res = config["rdir"] + "/output_tables/DB_genes_summary_info_exp.tsv",
        annot = config["rdir"] + "/output_tables/DB_cluster_annotations.tsv",
        metag = config["rdir"] + "output_tables/DB_clusters_in_metagenomes.tsv"
    log:
        out = "logs/tbl_stdout.log",
        err = "logs/tl_stderr.err"
    benchmark:
        "benchmarks/out_tbl.tsv"
    shell:
        """

        set -x
        set -e

        # Download agnostosDB environmental analysis data
        wget https://ndownloader.figshare.com/files/23066879 -O environmental.tar.gz
        tar -xzvf environmental.tar.gz

        # Download agnostosDB phylogentic analysis data
        wget https://ndownloader.figshare.com/files/23066864 -O phylogenetic.tar.gz
        tar -xzvf phylogenetic.tar.gz

        # Run R script to retrieve a general table containing a summary of the integration and contextual data
        awk '{split($1,a,"_\\+|_-"); print a[1]"\t"$1}' {input.genes} > {params.contig}
        # For Anvi'o data
        # awk 'NR>1{print $2"\t"$1}' anvio_gene_calls.tsv > {params.contig}

        ./{params.parser} --clu_or {input.clu_or} \
                          --contig {params.contig} \
                          --cat {input.cat} \
                          --clu_info {params.cval} \
                          --name {params.new_data_name} \
                          --comm {params.comm} \
                          --hq_clu {params.hq_clu} \
                          --k_annot {params.k_annot} \
                          --kwp_annot {params.kwp_annot} \
                          --gu_annot {params.gu_annots} \
                          --threads {threads} 2>{log.err} 1>{log.out}

        """

rule output_tables_done:
    input:
        res = config["rdir"] + "/output_tables/DB_genes_summary_info_exp.tsv",
        annot = config["rdir"] + "/output_tables/DB_cluster_annotations.tsv",
        metag = config["rdir"] + "output_tables/DB_clusters_in_metagenomes.tsv"
    output:
        res_done = touch(config["rdir"] + "/output_tables/res.done")
    run:
        shell("echo 'Parsing of results DONE'")
