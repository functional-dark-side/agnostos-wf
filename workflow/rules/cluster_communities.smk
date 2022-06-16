rule cluster_communities_dbs:
    input:
      k= config["rdir"] + "/cluster_category_DB/k_hhm_db.index",
      stats = config["rdir"] + "/cluster_category_stats/cluster_category_summary_stats.tsv"
    threads: 7
    conda:
      config["conda_env"]
    params:
      mmseqs_bin = config["mmseqs_bin"],
      mpi_runner = config["mpi_runner"],
      vmtouch = config["vmtouch"],
      hhblits = config["hhblits_bin_mpi"],
      hhparse = "scripts/categ_hhparser.sh",
      hhb_tmp_db = config["rdir"] + "/cluster_communities/hhbl_tmp",
    output:
      comm=config["rdir"] + "/cluster_communities/k_hhblits.tsv"
    log:
      out="logs/comm_stdout.log",
      err="logs/comm_stderr.err"
    benchmark:
      "benchmarks/comm_db.tsv"
    shell:
      """
      set -e
      set -x

      export OMPI_MCA_btl=^openib
      export OMP_NUM_THREADS={threads}
      export OMP_PROC_BIND=FALSE

      # HHblits all-clusters vs all-clusters for each category
      CATEG=$(echo -e "eu\ngu\nkwp\nk")

      IN=$(dirname {input.k})
      OUT=$(dirname {output.comm})

      for categ in $CATEG; do

        RES=${{OUT}}/${{categ}}_hhblits.tsv

        if [[ ! -s ${{RES}} ]]; then
            if [[ ! -s {params.hhb_tmp_db}.index  ]]; then
                {params.vmtouch} -f ${{IN}}/${{categ}}*
                {params.mpi_runner} {params.hhblits} -i ${{IN}}/${{categ}}_hhm_db \
                                                 -o {params.hhb_tmp_db} \
                                                 -n 2 -cpu 1 -v 0 \
                                                 -d ${{IN}}/${{categ}}
                mv {params.hhb_tmp_db}.ffdata {params.hhb_tmp_db}
                mv {params.hhb_tmp_db}.ffindex {params.hhb_tmp_db}.index
            fi
            {params.mpi_runner} {params.mmseqs_bin} apply {params.hhb_tmp_db} {params.hhb_tmp_db}.parsed \
                --threads 1 \
                -- {params.hhparse} 2>{log.err}

                sed -e 's/\\x0//g' {params.hhb_tmp_db}.parsed > ${{RES}} 2>{log.err}

                rm -rf {params.hhb_tmp_db} {params.hhb_tmp_db}.index {params.hhb_tmp_db}.dbtype {params.hhb_tmp_db}.parsed* {params.hhb_tmp_db}.ff*
        fi

      done
      """

rule cluster_communities_inference:
    input:
      k = config["rdir"] + "/cluster_communities/k_hhblits.tsv",
      stats = config["rdir"] + "/cluster_category_stats/cluster_category_summary_stats.tsv"
    threads: 7
    conda:
      config["conda_env"]
    params:
      get_comm = "scripts/communities_inference/get_communities.R",
      pfam_clan = config["pfam_clan"],
      db_mode = config["db_mode"],
      comm_config = "config/config_communities.yaml"
    output:
      comm=config["rdir"] + "/cluster_communities/cluster_communities.tsv"
    log:
      out="logs/comm_inf_stdout.log",
      err="logs/comm_inf_stderr.err"
    benchmark:
      "benchmarks/comm_inf.tsv"
    shell:
      """
      OUT=$(dirname {output.comm})
      
      if [ ! -s {params.pfam_clan} ]; then
        echo "Dowloading Pfam-A clan information"
        wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam34.0/Pfam-A.clans.tsv.gz -O {params.pfam_clan}
      fi

      # Start cluster community inference
      ./{params.get_comm} -c {params.comm_config} 1>{log.out} 2>{log.err}

      rm -rf ${{OUT}}/tmp

      if [[ {params.db_mode} == "memory" ]]; then
         rm {params.pfam_clan}
      fi

      """

rule cluster_comm_done:
    input:
      comm=config["rdir"] + "/cluster_communities/cluster_communities.tsv"
    output:
      comm_done = touch(config["rdir"] + "/cluster_communities/comm.done")
    run:
      shell("echo 'COMMUNITIES INFERENCE DONE'")
