rule mmseqs_clustering_update:
    input:
        new_seqs = config["rdir"] + "/gene_prediction/orf_seqs.fasta",
        annot = config["rdir"] + "/pfam_annotation/pfam_annot_parsed.tsv"
    params:
        mmseqs_bin = config["mmseqs_bin"],
        mmseqs_tmp = config["mmseqs_tmp"],
        mmseqs_local_tmp = config["mmseqs_local_tmp"],
        mmseqs_cov = 0.8,
        mmseqs_id = 0.3,
        mmseqs_cov_mode = 0,
        mmseqs_ens = 5,
        mmseqs_split_mem = config["mmseqs_split_mem"],
        or_dir = config["ordir"],
        or_cludb = config["ordir"] + "/mmseqs_clustering/cluDB",
        or_seqdb = config["ordir"] + "/mmseqs_clustering/seqDB",
        new_seqdb = config["rdir"] + "/mmseqs_clustering/new_seqDB",
        conc_seqdb = config["rdir"] + "/mmseqs_clustering/concat_seqDB",
        updt_seqdb = config["rdir"] + "/mmseqs_clustering/seqDB",
        updt_cludb = config["rdir"] + "/mmseqs_clustering/cluDB",
        mmseqs_split = config["mmseqs_split"],
        mmseqs_mpi_runner = config["mpi_runner"]
    threads: 28
    output:
        updt_clu = config["rdir"] + "/mmseqs_clustering/cluDB.tsv"
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

        # Create directory for original data to be integrated
        mkdir -p {params.or_dir}

        # Download original mmseqs clustering that you want to use for the integration
        if [[ ! -s {params.or_cludb}.index ]]; then
            wget https://ndownloader.figshare.com/files/23066651 -O {params.or_dir}/mmseqs_clustering.tar.gz
            tar -C {params.or_dir} -zxvf {params.or_dir}/mmseqs_clustering.tar.gz
            rm {params.or_dir}/mmseqs_clustering.tar.gz
        fi
        {params.mmseqs_bin} createdb {input.new_seqs} {params.new_seqdb} 2>{log.err} 1>{log.out}

        # Create symbolic link between original cluDB and seqDB (required by mmeseqs to run the update)

        ln -sf {params.or_seqdb}_h {params.or_cludb}_h
        ln -sf {params.or_seqdb}_h.index {params.or_cludb}_h.index
        ln -sf {params.or_seqdb}_h.dbtype {params.or_cludb}_h.dbtype
        ln -sf {params.or_seqdb}.lookup {params.or_cludb}.lookup

        # Concat new and old seq DBs
        {params.mmseqs_bin} concatdbs {params.or_seqdb} {params.new_seqdb} \
            {params.conc_seqdb} \
            --threads 1 #--preserve-keys 2>>{log.err} 1>>{log.out}
        {params.mmseqs_bin} concatdbs {params.or_seqdb}_h {params.new_seqdb}_h \
            {params.conc_seqdb}_h \
            --threads 1  2>>{log.err} 1>>{log.out}

        # Update the clusterDB:
        {params.mmseqs_bin} clusterupdate \
          {params.or_seqdb} \
          {params.conc_seqdb} \
          {params.or_cludb} \
          {params.updt_seqdb} \
          {params.updt_cludb} \
          {params.mmseqs_tmp} \
          --local-tmp {params.mmseqs_local_tmp} \
          --mpi-runner "{params.mmseqs_mpi_runner}" \
          --threads {threads} \
          -c {params.mmseqs_cov} \
          --cov-mode {params.mmseqs_cov_mode} \
          --min-seq-id {params.mmseqs_id} \
          -s {params.mmseqs_ens} 2>>{log.err} 1>>{log.out}

        {params.mmseqs_bin} createtsv \
          {params.updt_seqdb} \
          {params.updt_seqdb} \
          {params.updt_cludb} \
          {output.updt_clu} \
          --threads {threads} 2>>{log.err} 1>>{log.out}

        # Release some space removing not necessary files
        rm {params.or_cludb} {params.or_cludb}.index {params.or_cludb}_h {params.or_cludb}_h.index {params.or_cludb}.dbtype {params.or_cludb}_h.dbtype {params.or_cludb}.lookup
        rm {params.or_seqdb} {params.or_seqdb}.index {params.or_seqdb}_h {params.or_seqdb}_h.index {params.or_seqdb}.dbtype {params.or_seqdb}_h.dbtype {params.or_seqdb}.lookup

        """

rule mmseqs_cluster_update_done:
    input:
        updt_cludb = config["rdir"] + "/mmseqs_clustering/cluDB.tsv"
    output:
        updt_done = touch(
            config["rdir"] + "/mmseqs_clustering/updt.done")
    run:
        shell("echo 'MMSEQS2 CLUSTERING DONE'")
