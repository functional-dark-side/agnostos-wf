rule mmseqs_cluster_update:
    input:
        or_clu = config["rdir"] + "/mmseqs_clustering/cluDB",
        new_seqs = "data/...."
    params:
        mmseqs_bin = config["mmseqs_bin"],
        mmseqs_tmp = config["mmseqs_tmp"],
        mmseqs_local_tmp = config["mmseqs_local_tmp"],
        mmseqs_cov = 0.8,
        mmseqs_id = 0.3,
        mmseqs_cov_mode = 0,
        mmseqs_ens = 5,
        mmseqs_split_mem = config["mmseqs_split_mem"],
        or_seqdb = config["rdir"] + "/mmseqs_clustering/seqDB",
        new_seqdb = config["rdir"] + "/mmseqs_cluster_update/new_seqDB",
        conc_seqdb = config["rdir"] + "/mmseqs_cluster_update/concat_seqDB",
        updt_seqdb = config["rdir"] + "/mmseqs_cluster_update/updated_seqDB",
        mmseqs_split = config["mmseqs_split"],
        mmseqs_mpi_runner = config["mpi_runner"]
    threads: 28
    output:
        updt_cludb = config["rdir"] + "/mmseqs_cluster_update/updated_cluDB"
    log:
        out = "logs/mmseqs_clustering_stdout.log",
        err = "logs/mmseqs_clustering_stderr.err"
    benchmark:
        "benchmarks/mmseqs_clustering/clu.tsv"
    shell:
        """
        set -x
        set -e

        export OMPI_MCA_btl=^openib
        export OMP_NUM_THREADS={threads}
        export OMP_PROC_BIND=FALSE

        {params.mmseqs_bin} createdb {input.new_seqs} {params.new_seqdb} 2>{log.err} 1>{log.out}

        # Concat new and old seq DBs
        {params.mmseqs_bin} concatdbs {params.or_seqdb} {params.new_seqdb} \
            {params.conc_seqdb} \
            --threads {threads} #--preserve-keys

        RUNNER="{params.mmseqs_mpi_runner}" {params.mmseqs_bin} clusterupdate \
          {params.or_seqdb}
          {params.conc_seqdb} \
          {params.or_cludb} \
          {params.updt_seqDB} \
          {params.updt_cludb} \
          {params.mmseqs_tmp} \
          --threads {threads} \
          -c {params.mmseqs_cov} \
          --cov-mode {params.mmseqs_cov_mode} \
          --min-seq-id {params.mmseqs_id} \
          -s {params.mmseqs_ens} \
          --split-memory-limit {params.mmseqs_split_mem} \
          --split {params.mmseqs_split} 2>>{log.err} 1>>{log.out}

        """

rule mmseqs_update_done:
    input:
        updt_cludb = config["rdir"] + "/mmseqs_cluster_update/updated_cluDB"
    output:
        updt_done = touch(
            config["rdir"] + "/mmseqs_cluster_update/updt_done.tsv")
    run:
        shell("echo 'MMSEQS2 CLUSTERING DONE'")
