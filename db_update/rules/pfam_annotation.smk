rule pfam_annotation:
    input:
        gp=config["rdir"] + "/gene_prediction/orf_seqs.fasta"
    params:
        hmmer_bin = config["hmmer_bin"],
        mpi_runner = config["mpi_runner"],
        pfamdb = config["pfam_db"],
        hmmout = config["rdir"] + "/pfam_annotation/hmmsearch_pfam_annot.out",
        hmmlog = config["rdir"] + "/pfam_annotation/hmmsearch_pfam_annot.log"
    conda:
        "../envs/workflow.yml"
    output:
        annot=config["rdir"] + "/pfam_annotation/pfam_annotations.tsv"
    threads: 28
    log:
        out="logs/pfannot_stdout.log",
        err="logs/pfannot_stderr.err"
    benchmark:
        "benchmarks/pfam_annotation/pfannot.tsv"
    shell:
        """
        set -x
        set -e

        export OMPI_MCA_btl=^openib
        export OMP_NUM_THREADS={threads}
        export OMP_PROC_BIND=FALSE

        NPFAM=$(grep -c '^NAME' {params.pfamdb})
        NSEQS=$(grep -c '^>' {input.gp})
        N=$(($NSEQS * $NPFAM))

        # Run hmmsearch (MPI-mode)
        {params.mpi_runner} {params.hmmer_bin} --mpi --cut_ga -Z "${{N}}" --domtblout {params.hmmout} -o {params.hmmlog} {params.pfamdb} {input.gp} 2>{log.err} 1>{log.out}

        # Collect the results
        grep -v '^#' {params.hmmout} > {output.annot}

        """

rule pfam_annotation_done:
    input:
        annot = config["rdir"] + "/pfam_annotation/pfam_annotations.tsv"
    output:
        annot_done = touch(config["rdir"] + "/pfam_annotation/annotation.done")
    run:
        shell("echo 'ANNOTATION DONE'")
