rule gene_prediction:
    input:
        contigs = config['data']
    params:
        stage = config["data_stage"],
        prodigal_mode = config["prodigal_mode"],
        prodigal_bin = config["prodigal_bin"],
        data_partial = config["data_partial"],
        rename_orfs = "scripts/rename_orfs.awk",
        partial_info = "scripts/get_orf_partial_info.awk",
        gff_output = config["rdir"] + "/gene_prediction/orfs_info.gff",
        tmp = config["rdir"] + "/gene_prediction/tmpl"
    conda:
        config["conda_env"]
    output:
        fa = config["rdir"] + "/gene_prediction/orf_seqs.fasta",
        partial = config["rdir"] + "/gene_prediction/orf_partial_info.tsv"
    log:
        out = "logs/gene_stdout.log",
        err = "logs/gene_stderr.err"
    benchmark:
        "benchmarks/gene_prediction.tsv"
    shell:
        """
        set -x
        set -e

        if [[ {params.stage} = "contigs" ]]; then

            {params.prodigal_bin} -i {input.contigs} -a {output.fa} -m -p {params.prodigal_mode} -f gff  -o {params.gff_output} -q 2>{log.err} 1>{log.out}

            awk -f {params.rename_orfs} {output.fa} > {params.tmp} && mv {params.tmp} {output.fa}

            awk -f {params.partial_info} {params.gff_output} > {output.partial}

        elif [[ {params.stage} = "genes" ]]; then

            ln -sf {input.contigs} {output.fa}

            ln -sf {params.data_partial} {output.partial}

        elif [[ {params.stage} = "anvio_genes" ]]; then

            ln -sf {input.contigs} {output.fa}

            awk -vOFS="\\t" 'NR>1{{if($6==0) print $1,"00"; else print $1,"11";}}' {params.data_partial} > {output.partial}

        fi

        """

rule gene_prediction_done:
    input:
        fa = config["rdir"] + "/gene_prediction/orf_seqs.fasta",
        partial = config["rdir"] + "/gene_prediction/orf_partial_info.tsv"
    output:
        gp_done = touch(config["rdir"] + "/gene_prediction/gp.done")
    run:
        shell("echo 'GP DONE'")
