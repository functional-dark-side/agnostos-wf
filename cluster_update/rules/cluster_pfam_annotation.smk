rule cluster_pfam_annotation:
    input:
        clu=config["rdir"] + "/mmseqs_clustering/cluDB_no_singletons.tsv",
        annot=config["rdir"] + "/pfam_annotation/pfam_annot_parsed.tsv",
        partial = config["rdir"] + "/gene_prediction/orf_partial_info.tsv"
    priority: 1
    threads: 28
    params:
        cl_annotr = "scripts/clu_annot.r",
        pfam_clan = config["pfam_clan"],
        singl = config["rdir"] + "/mmseqs_clustering/cluDB_singletons.tsv",
        s_annot = config["rdir"] + "/annot_and_clust/singletons_pfam_annot.tsv",
        or_multi_annot = config["ordir"] + "/pfam_name_acc_clan_multi.tsv",
        multi_annot = config["rdir"] + "/annot_and_clust/pfam_name_acc_clan_multi.tsv",
        tmp = config["rdir"] + "/annot_and_clust/tmp",
        concat = "scripts/concat_multi_annot.awk"
    output:
        cl_annot=config["rdir"] + "/annot_and_clust/annotated_clusters.tsv",
        cl_noannot=config["rdir"] + "/annot_and_clust/not_annotated_clusters.tsv"
    log:
        out="logs/clu_annot_stdout.log",
        err="logs/clu_annot_stderr.err"
    benchmark:
        "benchmarks/cluster_pfam_annotation/clu_annot.tsv"
    shell:
        """
        set -e
        set -x

        ## 1. Pfam-clan info and multi-domain format for the ORF Pfam annotations

        # Add Pfam clan info to the Pfam annotation results
        join -11 -21 <( awk '{{print $1,$3,$8,$9}}' {input.annot} | sort -k1,1) \
        <( zcat {params.pfam_clan} | awk -vFS='\\t' -vOFS='\\t' '{{print $4,$1,$2}}' \
        | awk -vFS='\\t' -vOFS='\\t' '{{for(i=1; i<=NF; i++) if($i ~ /^ *$/) $i = "no_clan"}};1' \
        | sort -k1,1) > {params.multi_annot}

        # Reorder columns
        # the previous output is: pfam_name - orf - start - end - pfam_acc - pfam_clan
        # we want to get: orf - pfam_name - pfam_acc - pfam_clan - start - end
        awk -vOFS='\\t' '{{print $2,$1,$5,$6,$3,$4}}' {params.multi_annot} > {params.tmp} && mv {params.tmp} {params.multi_annot}

        # Multiple annotations on the same line, separated by “|” (the annotations were ordered first by alignment position)
        sort -k1,1 -k5,6g {params.multi_annot} | \
            awk -vFS='\\t' -vOFS='\\t' '{{print $1,$2,$3,$4}}' | sort -k1,1 | \
            awk -f {params.concat} | sed 's/ /\t/g' > {params.tmp}

        mv {params.tmp} {params.multi_annot}

        cat {params.or_multi_annot} >> {params.multi_annot}

        ## 2. Cluster annotations

        # The r script "clu_annot.r" distribute the Pfam annotation in the clusters,
        # creating two sets: "annotated_clusters" and "not_annotated_clusters"
        ./{params.cl_annotr} --pfam_annot {params.multi_annot} \
                             --clusters {input.clu} \
                             --partial {input.partial} \
                             --output_annot {output.cl_annot} \
                             --output_noannot {output.cl_noannot}

        ## 3. Singleton annotations

        join -12 -21 <(sort -k2,2 {params.singl} ) \
        <(sort -k1,1 {params.multi_annot}) > {params.s_annot}

        """

rule cluster_pfam_annotation_done:
    input:
        cl_annot=config["rdir"] + "/annot_and_clust/annotated_clusters.tsv",
        cl_noannot=config["rdir"] + "/annot_and_clust/not_annotated_clusters.tsv"
    output:
        annot_done = touch(config["rdir"] + "/annot_and_clust/cl_annot.done")
    run:
        shell("echo 'CLUSTER PFAM ANNOTATION DONE'")
