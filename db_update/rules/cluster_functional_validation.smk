rule cluster_functional_validation:
    input:
        cl_annot = config["rdir"] + "/annot_and_clust/annotated_clusters.tsv",
        sp_sh = config["rdir"] + "/spurious_shadow/spurious_shadow_info.tsv"
    threads: 28
    conda:
        "../envs/workflow.yml"
    params:
        funct_valr = "scripts/functional_validation.r",
        funct_val_fun = "scripts/funct_val_functions.r",
        pfam_shared_terms = config["pfam_shared_terms"]
    output:
        fval_res = config["rdir"] + \
            "/functional_validation/functional_val_results.tsv"
    log:
        out = "logs/fval_stdout.log",
        err = "logs/fval_stderr.err"
    benchmark:
        "benchmarks/functional_validation/fval.tsv"
    shell:
        """
        set -e
        set -x

        ./{params.funct_valr} --input {input.cl_annot} \
                              --pfam_terms {params.pfam_shared_terms} \
                              --output {output.fval_res} \
                              --functions {params.funct_val_fun} \
                              --threads {threads} 2>{log.err} 1>{log.out}

        """

rule cluster_fvalidation_done:
    input:
        fval_res = config["rdir"] + \
            "/functional_validation/functional_val_results.tsv"
    output:
        fval_done = touch(config['rdir'] + "/functional_validation/fval.done")
    run:
        shell("echo 'FUNCTIONAL VALIDATION DONE'")
