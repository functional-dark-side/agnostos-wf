rule cluster_functional_validation:
    input:
        cl_annot = config["rdir"] + "/annot_and_clust/annotated_clusters.tsv"
    threads: 28
    conda:
        config["conda_env"]
    params:
        funct_valr = "scripts/functional_validation.r",
        funct_val_fun = "scripts/funct_val_functions.r",
        pfam_shared_terms = config["pfam_shared_terms"],
        db_mode = config["db_mode"]
    output:
        fval_res = config["rdir"] + \
            "/validation/functional_val_results.tsv"
    log:
        out = "logs/fval_stdout.log",
        err = "logs/fval_stderr.err"
    benchmark:
        "benchmarks/fval.tsv"
    shell:
        """
        set -e
        set -x

        # Pfam list common domain terms
        if [ ! -s {params.pfam_shared_terms} ]; then
            echo "Dowloading Pfam list of shared domain names"
            wget https://figshare.com/ndownloader/files/31127782 -O {params.pfam_shared_terms}
        fi

        ./{params.funct_valr} --input {input.cl_annot} \
                              --pfam_terms {params.pfam_shared_terms} \
                              --output {output.fval_res} \
                              --functions {params.funct_val_fun} \
                              --threads {threads} 2>{log.err} 1>{log.out}

        if [[ {params.db_mode} == "memory" ]]; then
            rm {params.pfam_shared_terms}
        fi

        """

rule cluster_fvalidation_done:
    input:
        fval_res = config["rdir"] + \
            "/validation/functional_val_results.tsv"
    output:
        fval_done = touch(config['rdir'] + "/validation/fval.done")
    run:
        shell("echo 'FUNCTIONAL VALIDATION DONE'")
