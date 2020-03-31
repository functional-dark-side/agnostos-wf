rule cluster_validation_results:
      input:
        cval = config["rdir"] + "/compositional_validation/compositional_validation_results.tsv",
        fval = config["rdir"] + "/functional_validation/functional_val_results.tsv"
      threads: 28
      params:
        cl_annot = config["rdir"] + "/annot_and_clust/annotated_clusters.tsv",
        cl_noannot = config["rdir"] + "/annot_and_clust/not_annotated_clusters.tsv",
        val_res = "scripts/validation_summary.r",
        valdb = config["rdir"] + "/validation/validation_results.sqlite3",
        val_annot = config["rdir"] + "/validation/validation_annot_noannot.tsv",
        val_stats = config["rdir"] + "/validation/validation_results_stats.tsv",
        fval_plots = config["rdir"] + "/validation/functional_validation_plots.rda",
        cval_plots = config["rdir"] + "/validation/compositional_validation_plots.rda",
        tmp = config["rdir"] + "/validation/tmp"
      output:
        val_res = config["rdir"] + "/validation/validation_results.tsv",
        good = config["rdir"] + "/validation/good_clusters.tsv"
      log:
        out="logs/val_stdout.log",
        err="logs/val_stderr.err"
      benchmark:
        "benchmarks/validation/val.tsv"
      shell:
        """
        set -e
        set -x

        #Retrieve old cluster representatives information
        # Not annotated clusters
        #
        join -11 -21 <(awk '{{print $1,$2}}' {input.cval} | sort -k1,1 --parallel={threads}) \
          <(awk '!seen[$1]++{{print $1,$2,"noannot",$3}}' {params.cl_noannot} | sort -k1,1 --parallel={threads}) > {params.val_annot}
        # Annotated clusters
        join -11 -21 <(awk '{{print $1,$2}}' {input.cval} | sort -k1,1 --parallel={threads}) \
          <(awk '!seen[$1]++{{print $1,$2,"annot",$3}}' {params.cl_annot} | sort -k1,1 --parallel={threads} ) >> {params.val_annot}

        awk -vOFS='\\t' '{{print $1,$2,$3,$5,$4}}'  {params.val_annot} \
            > {params.tmp} && mv {params.tmp} {params.val_annot}

        # Combine with functional validation results
        # Results in SQlite as table of database and plots(as R objects)
         [[ -f {params.valdb} ]] && rm {params.valdb}

         ./{params.val_res} --valdb {params.valdb} \
                            --fval_res {input.fval} \
                            --cval_res {input.cval} \
                            --val_annot {params.val_annot} \
                            --val_res {output.val_res} \
                            --val_stats {params.val_stats} \
                            --good {output.good} \
                            --fplots {params.fval_plots} \
                            --cplots {params.cval_plots} 2>{log.err} 1>{log.out}
        """

rule cluster_validation_done:
    input:
        val_res = config["rdir"] + "/validation/validation_results.tsv",
        good = config["rdir"] + "/validation/good_clusters.tsv"
    output:
        val_done = touch(config["rdir"] + "/validation/val.done")
    run:
        shell("echo 'CLUSTER VALIDATION DONE'")
