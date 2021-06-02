rule gene_prediction:
    input:
        contigs = config['data'] + "/{smp}_contigs.fasta"
    params:
        prodigal_mode = config["prodigal_mode"],
        prodigal_bin = config["prodigal_bin"],
        rename_orfs = "scripts/rename_orfs.awk",
        partial_info = "scripts/get_orf_partial_info.awk",
        gff_output = config["rdir"] + "/gene_prediction/{smp}_orfs_info.gff",
        tmp = config["rdir"] + "/gene_prediction/{smp}_tmpl",
        outdir = config["rdir"] + "/gene_prediction"
    conda:
        config["conda_env"]
    output:
        fa = config["rdir"] + "/gene_prediction/{smp}_orfs.fasta",
        orfs = config["rdir"] + "/gene_prediction/{smp}_orfs.txt",
        partial = config["rdir"] + "/gene_prediction/{smp}_partial_info.tsv"
    log:
        out = "logs/{smp}_gene_stdout.log",
        err = "logs/{smp}_gene_stderr.err"
    benchmark:
        "benchmarks/gene_prediction/{smp}.gp.tsv"
    shell:
        """
        set -x
        set -e

        {params.prodigal_bin} -i {input.contigs} -a {output.fa} -m -p {params.prodigal_mode} -f gff  -o {params.gff_output} -q 2>{log.err} 1>{log.out}

        awk -f {params.rename_orfs} {output.fa} > {params.tmp} && mv {params.tmp} {output.fa}

        grep '^>' {output.fa} | sed "s/^>//" > {output.orfs}

        awk -f {params.partial_info} {params.gff_output} > {output.partial}

        """

rule gene_prediction_done:
    input:
        fa = config["rdir"] + "/gene_prediction/{smp}_orfs.fasta",
        orfs = config["rdir"] + "/gene_prediction/{smp}_orfs.txt",
        partial = config["rdir"] + "/gene_prediction/{smp}_partial_info.tsv"
    output:
        gp_done = touch(config["rdir"] + "/gene_prediction/{smp}.gp.done")
    run:
        shell("echo 'GP DONE'")
