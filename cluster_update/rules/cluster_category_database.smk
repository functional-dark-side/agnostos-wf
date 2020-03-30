rule cluster_category_database:
    input:
        k = config["rdir"] + "/cluster_categories/k_ids.txt",
        kwp = config["rdir"] + "/cluster_categories/kwp_ids.txt",
        gu = config["rdir"] + "/cluster_categories/gu_ids.txt",
        eu = config["rdir"] + "/cluster_categories/eu_ids.txt",
        index = config["rdir"] + "/mmseqs_clustering/new_cluDB_name_index.txt"
    threads: 28
    params:
        mmseqs_bin = config["mmseqs_bin"],
        mpi_runner = config["mpi_runner"],
        famsa_bin = config["famsa_bin"],
        categ_db = "scripts/category_db_files.sh",
        hhsuite = config["hhsuite"],
        reformat = "scripts/reformat_file.sh",
        consensus = "scripts/consensus.sh",
        hhmake = "scripts/hhmake.sh",
        step = "category database",
        outdir = config["rdir"] + "/cluster_category_DB",
        idir = config["rdir"] + "/cluster_categories",
        cluseqdb = config["rdir"] + "/cluster_refinement/refined_clusterDB",
        clu_hhm = config["rdir"] + "/cluster_category_DB/clu_hhm_db"
    log:
        out = "logs/catdb_stdout.log",
        err = "logs/catdb_stderr.err"
    benchmark:
        "benchmarks/cluster_category_DB/catdb.tsv"
    output:
        eu_cons = config["rdir"] + "/cluster_category_DB/eu_cons.index",
        eu_hmm = config["rdir"] + "/cluster_category_DB/eu_hhm_db.index",
        k_cons = config["rdir"] + "/cluster_category_DB/k_cons.index",
        k_hmm = config["rdir"] + "/cluster_category_DB/k_hhm_db.index",
        gu_cons = config["rdir"] + "/cluster_category_DB/gu_cons.index",
        gu_hmm = config["rdir"] + "/cluster_category_DB/gu_hhm_db.index",
        kwp_cons = config["rdir"] + "/cluster_category_DB/kwp_cons.index",
        kwp_hmm = config["rdir"] + "/cluster_category_DB/kwp_hhm_db.index",
        clu_hmm = config["rdir"] + "/cluster_category_DB/clu_hmm_db"
    shell:
        """
        set -e
        set -x

        export OMPI_MCA_btl=^openib
        export OMP_NUM_THREADS={threads}
        export OMP_PROC_BIND=FALSE

        {params.categ_db} --cluseq_db {params.cluseqdb} \
                          --step "{params.step}" \
                          --index {input.index} \
                          --mpi_runner "{params.mpi_runner}" \
                          --mmseqs {params.mmseqs_bin} \
                          --aln {params.famsa_bin} \
                          --hhsuite {params.hhsuite} \
                          --reformat {params.reformat} \
                          --consensus {params.consensus} \
                          --hhmake {params.hhmake} \
                          --outdir {params.outdir} \
                          --idir {params.idir} \
                          --clu_hhm {params.clu_hhm} \
                          --threads {threads} 2>{log.err} 1>{log.out}

        # Create a comprehensive HMMs cluster DB in MMseqs format to perform profile searches
        {params.mmseqs_bin} convertprofiledb {params.clu_hhm} {output.clu_hmm} \
            --threads {threads} 2>{log.err}

        rm {params.clu_hhm} {params.clu_hhm}.index {params.clu_hhm}.dbtype

        """

rule cluster_categ_db_done:
    input:
        eu_cons = config["rdir"] + "/cluster_category_DB/eu_cons.index",
        eu_hmm = config["rdir"] + "/cluster_category_DB/eu_hhm_db.index",
        k_cons = config["rdir"] + "/cluster_category_DB/k_cons.index",
        k_hmm = config["rdir"] + "/cluster_category_DB/k_hhm_db.index",
        gu_cons = config["rdir"] + "/cluster_category_DB/gu_cons.index",
        gu_hmm = config["rdir"] + "/cluster_category_DB/gu_hhm_db.index",
        kwp_cons = config["rdir"] + "/cluster_category_DB/kwp_cons.index",
        kwp_hmm = config["rdir"] + "/cluster_category_DB/kwp_hhm_db.index",
        clu_hmm = config["rdir"] + "/cluster_category_DB/clu_hmm_db"
    output:
        cat_db_done = touch(config["rdir"] + "/cluster_category_DB/cat_db.done")
    run:
        shell("echo 'CATEGORY DB DONE'")
