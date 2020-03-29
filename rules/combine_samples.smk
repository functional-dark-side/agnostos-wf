def files(wildcards, input):
    return " ".join(input)

rule combine_orf_fasta:
    input:
        orfsfa = expand(config["rdir"] + "/gene_prediction/{smp}_orfs.fasta", smp=SAMPLES)
    params:
        files = files
    output:
        all_orfsfa = config["rdir"] + "/combine_samples/orf_seqs.fasta"
    log:
        out = "logs/combining_stdout.log",
        err = "logs/combining_stderr.err"
    benchmark:
        "benchmarks/combine_samples/comb1.tsv"
    shell:
        """
        # Combine predicted ORFs fasta sequences from different samples
        cat {params.files} > {output.all_orfsfa}
        """

rule combine_orf_headers:
    input:
        orfstxt = expand(config["rdir"] + "/gene_prediction/{smp}_orfs.txt", smp=SAMPLES)
    params:
        files = files
    output:
        all_orfstxt = config["rdir"] + "/combine_samples/orf_seqs.txt"
    log:
        out = "logs/combining_stdout.log",
        err = "logs/combining_stderr.err"
    benchmark:
        "benchmarks/combine_samples/comb2.tsv"
    shell:
        """
        # Combine predicted ORFs headers from different samples
        cat {params.files} > {output.all_orfstxt}
        """

rule combine_partial_info:
    input:
        orfstxt = expand(config["rdir"] + "/gene_prediction/{smp}_partial_info.tsv", smp=SAMPLES)
    params:
        files = files
    output:
        all_partial = config["rdir"] + "/combine_samples/orf_partial_info.tsv"
    log:
        out = "logs/combining_stdout.log",
        err = "logs/combining_stderr.err"
    benchmark:
        "benchmarks/combine_samples/comb3.tsv"
    shell:
        """
        # Combine predicted ORFs partiality information from different samples
        cat {params.files} > {output.all_partial}
        """

rule combine_orf_annot:
    input:
        annot = expand(config["rdir"] + "/pfam_annotation/{smp}_pfam_parsed.tsv", smp=SAMPLES)
    params:
        files = files
    output:
        all_annot = config["rdir"] + "/combine_samples/pfam_annot_parsed.tsv"
    log:
        out = "logs/combining_stdout.log",
        err = "logs/combining_stderr.err"
    benchmark:
        "benchmarks/combine_samples/comb4.tsv"
    shell:
        """
        # Combine Pfam annotation results
        cat {params.files} > {output.all_annot}
        """

rule combine_samples_done:
    input:
        all_annot = config["rdir"] + "/combine_samples/pfam_annot_parsed.tsv",
        all_partial = config["rdir"] + "/combine_samples/orf_partial_info.tsv",
        all_orfstxt = config["rdir"] + "/combine_samples/orf_seqs.txt",
        all_orfsfa = config["rdir"] + "/combine_samples/orf_seqs.fasta"

    output:
        annot_done = touch(config["rdir"] + "/combine_samples/comb.done")
    run:
            shell("echo 'COMBINATION DONE'")
