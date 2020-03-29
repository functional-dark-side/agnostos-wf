rule pfam_annotation_results:
    input:
        annot=config["rdir"] + "/pfam_annotation/pfam_annotations.tsv"
    params:
        evalue=1e-05,
        coverage=0.4
    output:
        pf_annot=config["rdir"] + "/pfam_annotation/pfam_annot_parsed.tsv"
    log:
        err="logs/annot_stderr.err"
    benchmark:
        "benchmarks/pfam_annotation/pfannot_res.tsv"
    shell:
        "scripts/hmmsearch_res_parser.sh {input.annot} {params.evalue} {params.coverage} > {output.pf_annot} 2>>{log.err}"

rule pfam_annotation_res_done:
    input:
        annot = config["rdir"] + "/pfam_annotation/pfam_annot_parsed.tsv"
    output:
        annot_done = touch(config["rdir"] + "/pfam_annotation/annot_res.done")
    run:
        shell("echo 'ANNOTATION RES DONE'")
