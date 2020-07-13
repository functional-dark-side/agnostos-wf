rule workflow_report:
    input:
        cat_db = config["rdir"] + "/cluster_category_DB/k_hhm_db.index",
        clu_cat = config["rdir"] + "/clusterDB_results/cluster_ids_categ.tsv.gz",
        comm = config["rdir"] + "/cluster_communities/cluster_communities.tsv"
    threads: 28
    conda:
        config["conda_env"]
    params:
        basedir = config["rdir"],
        outdir = config["rdir"] + "/report/",
        input_data_dir = config["data"],
        report_maker = "scripts/report_maker.r",
        wf_report = "scripts/workflow_report.Rmd"
    output:
        report=config["rdir"] + "/report/workflow_report.html"
    log:
        out = "logs/report_stdout.log",
        err = "logs/report_stderr.err"
    benchmark:
        "benchmarks/wf-report.tsv"
    shell:
        """
        Rscript --vanilla {params.report_maker} --basedir {params.basedir} \
                                                --outdir  {params.outdir} \
                                                --input {params.input_data_dir} \
                                                --wf_report {params.wf_report} \
                                                --output {output.report} 2>{log.err} 1>{log.out}

        """
