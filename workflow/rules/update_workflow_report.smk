rule update_workflow_report:
    input:
        iclu_com = config["rdir"] + "/integrated_cluster_DB/cluster_communities.tsv.gz",
        comm = config["rdir"] + "/cluster_communities/cluster_communities.tsv",
        annot = config["rdir"] + "/output_tables/DB_cluster_annotations.tsv"
    threads: 28
    conda:
        config["conda_env"]
    params:
        basedir = config["rdir"],
        outdir = config["rdir"] + "/report/",
        input_data = config["data"],
        name_data = config["data_name"],
        stage = config["data_stage"],
        report_maker = "scripts/report_maker.r",
        wf_report = "scripts/update_workflow_report.Rmd"
    output:
        report=config["rdir"] + "/report/workflow_report.html"
    log:
        out = "logs/report_stdout.log",
        err = "logs/report_stderr.err"
    benchmark:
        "benchmarks/updt_wf-report.tsv"
    shell:
        """
        Rscript --vanilla {params.report_maker} --basedir {params.basedir} \
                                                --outdir  {params.outdir} \
                                                --name {params.name_data} \
                                                --input {params.input_data} \
                                                --stage {params.stage} \
                                                --wf_report {params.wf_report} \
                                                --output {output.report} 2>{log.err} 1>{log.out}

        """
