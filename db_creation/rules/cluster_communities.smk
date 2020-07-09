rule cluster_communities:
    input:
      k= config["rdir"] + "/cluster_category_DB/k_hhm_db.index",
      stats = config["rdir"] + "/cluster_category_stats/cluster_category_summary_stats.tsv"
    threads: 7
    conda:
        "../envs/workflow.yml"
    params:
      mmseqs_bin = config["mmseqs_bin"],
      mpi_runner = config["mpi_runner"],
      hhblits = config["hhblits_bin"],
      hhparse = "scripts/categ_hhparser.sh",
      hhb_tmp_db = config["rdir"] + "/cluster_communities/hhbl_tmp",
      get_comm = "scripts/communities_inference/get_communities.R",
      comm_config = "config/config_communities.yml"
    output:
      comm=config["rdir"] + "/cluster_communities/cluster_communities.tsv"
    log:
      out="logs/comm_stdout.log",
      err="logs/comm_stderr.err"
    benchmark:
      "benchmarks/cluster_communities/comm.tsv"
    shell:
      """
      # HHblits all-clusters vs all-clusters for each category
      CATEG=$(echo -e "eu\ngu\nkwp\nk")

      IN=$(dirname {input.k})
      OUT=$(dirname {output.comm})

      for categ in $CATEG; do

        RES=${{OUT}}/${{categ}}_hhblits.tsv

        if [[ ! -s ${{RES}} ]]; then
            {params.mpi_runner} {params.mmseqs_bin} apply ${{IN}}/${{categ}}_hhm_db {params.hhb_tmp_db} \
                --threads {threads} \
                -- {params.hhparse} {params.hhblits} ${{IN}}/${{categ}} ${{OUT}}/hhr 2>{log.err}

                sed -e 's/\\x0//g' {params.hhb_tmp_db} > ${{RES}} 2>{log.err}

                rm -rf {params.hhb_tmp_db} {params.hhb_tmp_db}.index {params.hhb_tmp_db}.dbtype ${{OUT}}/hhr
        fi

      done

      # Start cluster community inference
      ./{params.get_comm} -c {params.comm_config} 1>{log.out} 2>{log.err}
      
      rm -rf ${{OUT}}/tmp
      """

rule cluster_comm_done:
    input:
      comm=config["rdir"] + "/cluster_communities/cluster_communities.tsv"
    output:
      comm_done = touch(config["rdir"] + "/cluster_communities/comm.done")
    run:
      shell("echo 'COMMUNITIES INFERENCE DONE'")
