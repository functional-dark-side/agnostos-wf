rule pfam_annotation:
    input:
        gp=config["rdir"] + "/gene_prediction/orf_seqs.fasta"
    params:
        hmmer_bin = config["hmmer_bin"],
        mpi_runner = config["mpi_runner"],
        pfamdb = config["pfam_db"],
        db_mode = config["db_mode"],
        parser = config["wdir"] + "/scripts/hmmsearch_res_parser.sh",
        hmmout = config["rdir"] + "/pfam_annotation/hmmsearch_pfam_annot.out",
        hmmlog = config["rdir"] + "/pfam_annotation/hmmsearch_pfam_annot.log",
        annot=config["rdir"] + "/pfam_annotation/pfam_annotations.tsv",
        evalue=1e-05,
        coverage=0.4
    conda:
        config["conda_env"]
    output:
        pf_annot=config["rdir"] + "/pfam_annotation/pfam_annot_parsed.tsv"
    threads: 28
    log:
        out="logs/pfannot_stdout.log",
        err="logs/pfannot_stderr.err"
    benchmark:
        "benchmarks/pfam_annot.tsv"
    shell:
        """
        set -x
        set -e

        export OMPI_MCA_btl=^openib
        export OMP_NUM_THREADS={threads}
        export OMP_PROC_BIND=FALSE

        # Pfam database
        if [ ! -s {params.pfamdb} ]; then
            DB=$(dirname {params.pfamdb})
            mkdir -p ${{DB}}
            echo "Dowloading Pfam-A database"
            wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam34.0/Pfam-A.hmm.gz -O {params.pfamdb}.gz
            gunzip {params.pfamdb}.gz
        fi

        NPFAM=$(grep -c '^NAME' {params.pfamdb})
        NSEQS=$(grep -c '^>' {input.gp})
        N=$(($NSEQS * $NPFAM))

        if [ ! -s {params.annot} ]; then
            # Run hmmsearch (MPI-mode)
            {params.mpi_runner} {params.hmmer_bin} --mpi --cut_ga -Z "${{N}}" --domtblout {params.hmmout} -o {params.hmmlog} {params.pfamdb} {input.gp} 2>{log.err} 1>{log.out}

            # Collect the results
            grep -v '^#' {params.hmmout} > {params.annot}
        fi
        
        # Parse the results
        {params.parser} {params.annot} {params.evalue} {params.coverage} > {output.pf_annot} 2>>{log.err}

        # remove the search log files
        rm {params.hmmout} {params.hmmlog}
        # remove the unparsed pfam_annotations
        rm {params.annot}

        if [[ {params.db_mode} == "memory" ]]; then
            rm {params.pfamdb}
        fi

        """

rule pfam_annotation_done:
    input:
        pf_annot=config["rdir"] + "/pfam_annotation/pfam_annot_parsed.tsv"
    output:
        annot_done = touch(config["rdir"] + "/pfam_annotation/annotation.done")
    run:
        shell("echo 'ANNOTATION DONE'")
