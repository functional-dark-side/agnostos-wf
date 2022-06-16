rule cluster_category_database:
    input:
        k = config["rdir"] + "/cluster_categories/k_ids.txt",
        kwp = config["rdir"] + "/cluster_categories/kwp_ids.txt",
        gu = config["rdir"] + "/cluster_categories/gu_ids.txt",
        eu = config["rdir"] + "/cluster_categories/eu_ids.txt"
    threads: 28
    params:
        mmseqs_bin = config["mmseqs_bin"],
        mpi_runner = config["mpi_runner"],
        famsa_bin = config["famsa_bin"],
        hhsuite = config["hhsuite"],
        module = config["module"],
        step = "category database",
        categ_db = "scripts/category_db_files.sh",
        reformat = "scripts/reformat_file.sh",
        consensus = "scripts/consensus.sh",
        hhmake = "scripts/hhmake.sh",
        cluseqdb = config["rdir"] + "/cluster_refinement/refined_clusterDB",
        idir = config["rdir"] + "/cluster_categories",
        outdir = config["rdir"] + "/cluster_category_DB"
    log:
        out = "logs/catdb_stdout.log",
        err = "logs/catdb_stderr.err"
    benchmark:
        "benchmarks/cat_db.tsv"
    output:
        k_cons = config["rdir"] + "/cluster_category_DB/k_cons.index",
        k_hmm = config["rdir"] + "/cluster_category_DB/k_hhm_db.index",
        clu_hhm = config["rdir"] + "/cluster_category_DB/clu_hhm_db"
    shell:
        """
        set -e
        set -x

        export OMPI_MCA_btl=^openib
        export OMP_NUM_THREADS={threads}
        export OMP_PROC_BIND=FALSE

        {params.categ_db} --cluseq_db {params.cluseqdb} \
                          --step "{params.step}" \
                          --mpi_runner "{params.mpi_runner}" \
                          --mmseqs {params.mmseqs_bin} \
                          --aln {params.famsa_bin} \
                          --hhsuite {params.hhsuite} \
                          --reformat {params.reformat} \
                          --consensus {params.consensus} \
                          --hhmake {params.hhmake} \
                          --outdir {params.outdir} \
                          --idir {params.idir} \
                          --clu_hhm {output.clu_hhm} \
                          --threads {threads} 2>{log.err} 1>{log.out}

        """

rule cluster_categ_db_done:
    input:
        k_cons = config["rdir"] + "/cluster_category_DB/k_cons.index",
        k_hmm = config["rdir"] + "/cluster_category_DB/k_hhm_db.index",
        clu_hhm = config["rdir"] + "/cluster_category_DB/clu_hhm_db"
    output:
        cat_db_done = touch(config["rdir"] + "/cluster_category_DB/cat_db.done")
    run:
        shell("echo 'CATEGORY DB DONE'")
