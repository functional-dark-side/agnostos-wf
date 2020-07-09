rule cluster_compositional_validation:
    input:
        fval = config["rdir"] + \
            "/functional_validation/functional_val_results.tsv"
    params:
        mpi_runner = config["mpi_runner"],
        mmseqs_bin = config["mmseqs_bin"],
        famsa_bin = config["famsa_bin"],
        odseq_bin = config["odseq_bin"],
        parasail_bin = config["parasail_bin"],
        seqtk_bin = config["seqtk_bin"],
        datamash = config["datamash_bin"],
        parallel_bin = config["parallel_bin"],
        threads_collect = config["threads_collect"],
        igraph_lib = config["igraph_lib"],
        cvals = config["wdir"] + "/scripts/compositional_validation.sh",
        stats = config["wdir"] + "/scripts/get_stats.r",
        isconnect = config["wdir"] + "/scripts/is_connected",
        filterg = config["wdir"] + "/scripts/filter_graph",
        collect = config["wdir"] + "/scripts/collect_val_results.sh",
        clseqdb = config["rdir"] + "/mmseqs_clustering/new_clu_seqDB",
        outdb = config["rdir"] + "/compositional_validation/comp_valDB",
        cl_index = config["rdir"] + "/mmseqs_clustering/new_cluDB_name_index.txt",
        stat_dir = config["rdir"] + "/compositional_validation/stats",
        outdir = config["rdir"] + "/compositional_validation"
    conda:
        config["conda_env"]
    threads: 7
    output:
        cl_cval = config["rdir"] + \
            "/compositional_validation/compositional_validation_results.tsv",
        cval_rej = config["rdir"] + \
            "/compositional_validation/compositional_validation_rejected_orfs.tsv"
    log:
        out = "logs/cval_stdout.log",
        err = "logs/cval_stderr.err"
    benchmark:
        "benchmarks/compositional_validation/cval.tsv"
    shell:
        """
        set -e
        set -x

        export OMPI_MCA_btl=^openib
        export OMP_NUM_THREADS={threads}
        export OMP_PROC_BIND=FALSE

        {params.igraph_lib}

        # Run compositional validation in mpi-mode
        {params.mpi_runner} {params.mmseqs_bin} apply {params.clseqdb} {params.outdb} --threads {threads} \
            -- {params.cvals} --derep {params.mmseqs_bin} \
                       --msa {params.famsa_bin} \
                       --msaeval {params.odseq_bin} \
                       --ssn {params.parasail_bin} \
                       --gfilter {params.filterg} \
                       --gconnect {params.isconnect} \
                       --seqp {params.seqtk_bin} \
                       --datap {params.datamash} \
                       --stats {params.stats} \
                       --index {params.cl_index} \
                       --outdir {params.outdir} \
                       --threads {threads} 2>{log.err} 1>{log.out}

        # Collect results:
        # collect cluster main compositional validation stats and cluster rejected (bad-aligned) ORFs
        {params.collect} {params.stat_dir} {output.cl_cval} {output.cval_rej} {params.parallel_bin} {params.threads_collect} 2>{log.err} 1>{log.out}

        rm -rf {params.outdb} {params.outdb}.index {params.outdb}.dbtype

        """


rule cluster_cvalidation_done:
    input:
        cl_cval = config["rdir"] + \
            "/compositional_validation/compositional_validation_results.tsv",
        cval_rej = config["rdir"] + \
            "/compositional_validation/compositional_validation_rejected_orfs.tsv"
    output:
        cval_done = touch(
            config["rdir"] + "/compositional_validation/cval.done")
    run:
        shell("echo 'COMPOSITIONAL VALIDATION DONE'")
