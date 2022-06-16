rule cluster_compositional_validation:
    input:
        fval = config["rdir"] + \
            "/validation/functional_val_results.tsv"
    params:
        mpi_runner = config["mpi_runner"],
        mmseqs_bin = config["mmseqs_bin"],
        famsa_bin = config["famsa_bin"],
        odseq_bin = config["odseq_bin"],
        parasail_bin = config["parasail_bin"],
        seqtk_bin = config["seqtk_bin"],
        datamash = config["datamash_bin"],
        parallel_bin = config["parallel_bin"],
        igraph_lib = config["igraph_lib"],
        parasail_lib = config["parasail_lib"],
        threads_collect = config["threads_collect"],
        module = config["module"],
        cvals = config["wdir"] + "/scripts/compositional_validation.sh",
        stats = config["wdir"] + "/scripts/get_stats.r",
        isconnect = config["wdir"] + "/scripts/is_connected",
        filterg = config["wdir"] + "/scripts/filter_graph",
        collect = config["wdir"] + "/scripts/collect_val_results.sh",
        outdir = config["rdir"] + "/validation",
        stat_dir = config["rdir"] + "/validation/stats",
        clseqdb = config["rdir"] + "/mmseqs_clustering/clu_seqDB",
        new_clseqdb = config["rdir"] + "/mmseqs_clustering/new_clu_seqDB",
        outdb = config["rdir"] + "/validation/comp_valDB"
    conda:
        config["conda_env"]
    threads: 7
    priority: 50
    output:
        cl_cval = config["rdir"] + \
            "/validation/compositional_validation_results.tsv",
        cval_rej = config["rdir"] + \
            "/validation/compositional_validation_rejected_orfs.tsv"
    log:
        out = "logs/cval_stdout.log",
        err = "logs/cval_stderr.err"
    benchmark:
        "benchmarks/cval.tsv"
    shell:
        """
        set -e
        set -x

        export OMPI_MCA_btl=^openib
        export OMP_NUM_THREADS={threads}
        export OMP_PROC_BIND=FALSE

        {params.igraph_lib}
        {params.parasail_lib}

        # Run compositional validation in mpi-mode
        if [[ {params.module} == "creation" ]]; then
            DB={params.clseqdb}
        else
            DB={params.new_clseqdb}
        fi

        {params.mpi_runner} {params.mmseqs_bin} apply ${{DB}} {params.outdb} \
            -- {params.cvals} --derep {params.mmseqs_bin} \
                       --msa {params.famsa_bin} \
                       --msaeval {params.odseq_bin} \
                       --ssn {params.parasail_bin} \
                       --gfilter {params.filterg} \
                       --gconnect {params.isconnect} \
                       --seqp {params.seqtk_bin} \
                       --datap {params.datamash} \
                       --stats {params.stats} \
                       --outdir {params.outdir} \
                       --slots {threads} \
                       --threads {threads} 2>{log.err} 1>{log.out}

        # Collect results:
        # collect cluster main compositional validation stats and cluster rejected (bad-aligned) ORFs
        {params.collect} {params.stat_dir} {output.cl_cval} {output.cval_rej} {params.parallel_bin} {params.threads_collect} 2>>{log.err} 1>>{log.out}

        rm -rf {params.outdb} {params.outdb}.index {params.outdb}.dbtype
        rm -rf {params.outdir}/stats {params.outdir}/log_vals {params.outdir}/rejected

        """


rule cluster_cvalidation_done:
    input:
        cl_cval = config["rdir"] + \
            "/validation/compositional_validation_results.tsv",
        cval_rej = config["rdir"] + \
            "/validation/compositional_validation_rejected_orfs.tsv"
    output:
        cval_done = touch(
            config["rdir"] + "/validation/cval.done")
    run:
        shell("echo 'COMPOSITIONAL VALIDATION DONE'")
