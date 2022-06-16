rule mmseqs_clustering:
    input:
        orfs = config["rdir"] + "/gene_prediction/orf_seqs.fasta"
    params:
        mmseqs_bin = config["mmseqs_bin"],
        mmseqs_mpi_runner = config["mpi_runner"],
        mmseqs_local_tmp = config["mmseqs_local_tmp"],
        mmseqs_cov = 0.8,
        mmseqs_id = 0.3,
        mmseqs_cov_mode = 0,
        mmseqs_ens = 5,
        mmseqs_tmp = config["rdir"] + "/mmseqs_clustering/tmp",
        seqdb = config["rdir"] + "/mmseqs_clustering/seqDB",
        cludb = config["rdir"] + "/mmseqs_clustering/cluDB"

    threads: 28
    priority: 50
    output:
        clu = config["rdir"] + "/mmseqs_clustering/cluDB.tsv"
    log:
        out = "logs/mmseqs_clustering_stdout.log",
        err = "logs/mmseqs_clustering_stderr.err"
    benchmark:
        "benchmarks/clustering.tsv"
    shell:
        """
        set -x
        set -e

        export OMPI_MCA_btl=^openib
        export OMP_NUM_THREADS={threads}
        export OMP_PROC_BIND=FALSE

        {params.mmseqs_bin} createdb {input.orfs} {params.seqdb} 2>{log.err} 1>{log.out}

          {params.mmseqs_bin} cluster \
          {params.seqdb} \
          {params.cludb} \
          {params.mmseqs_tmp} \
          --local-tmp {params.mmseqs_local_tmp} \
          --threads {threads} \
          -c {params.mmseqs_cov} \
          --cov-mode {params.mmseqs_cov_mode} \
          --min-seq-id {params.mmseqs_id} \
          -s {params.mmseqs_ens} \
          --mpi-runner "{params.mmseqs_mpi_runner}" 2>>{log.err} 1>>{log.out}

        {params.mmseqs_bin} createtsv {params.seqdb} {params.seqdb} {params.cludb} {output.clu}
        """


rule mmseqs_clustering_done:
    input:
        clu = config["rdir"] + "/mmseqs_clustering/cluDB.tsv"
    output:
        cls_done = touch(config["rdir"] + "/mmseqs_clustering/clu.done")
    run:
        shell("echo 'MMSEQS2 CLUSTERING DONE'")
