rule cluster_classification:
    input:
        ref_noannot = config["rdir"] + \
            "/cluster_refinement/refined_not_annotated_clusters.tsv",
        ref_annot = config["rdir"] + \
            "/cluster_refinement/refined_annotated_clusters.tsv"
    threads: 28
    conda:
       config["conda_env"]
    params:
        wd = config["wdir"],
        mmseqs_bin = config["mmseqs_bin"],
        ffindex_apply = config["ffindex_apply"],
        mpi_runner = config["mpi_runner"],
        filterbyname = config["filterbyname"],
        famsa_bin = config["famsa_bin"],
        consensus = "scripts/consensus.sh",
        hhcons_bin = config["hhcons_bin"],
        mmseqs_local_tmp = config["mmseqs_local_tmp"],
        mmseqs_search = "scripts/mmseqs_double_search.sh",
        evalue_filt = "scripts/evalue_filter.awk",
        ripgrep = "rg",
        seqkit_bin = config["seqkit_bin"],
        dom_arch = config["rdir"] + \
            "/cluster_classification/cluster_pfam_domain_architectures.tsv",
        da_r = "scripts/cluster_pfam_domain_architectures.r",
        patterns = "scripts/hypothetical_grep.tsv",
        uniref_fasta = config["uniref90_db"],
        uniref_prot = config["uniref90_prot"],
        uniref_res = config["rdir"] + \
            "/cluster_classification/noannot_vs_uniref90.tsv",
        nr_fasta = config["nr_db"],
        nr_prot = config["nr_prot"],
        nr_res = config["rdir"] + \
            "/cluster_classification/uniref-nohits_vs_NR.tsv",
        evalue_thr = 0.6,
        hypo_thr = 1.0,
        name_index = config["rdir"] + "/mmseqs_clustering/cluDB_name_index.txt",
        good_cl = config["rdir"] + "/validation/good_clusters.tsv",
        refdb_noannot = config["rdir"] + \
            "/cluster_refinement/refined_not_annotated_clusterDB",
        refindex_noannot = config["rdir"] + \
            "/cluster_refinement/refined_not_annotated_clusterDB.index",
        ref_cons = config["rdir"] + \
            "/cluster_classification/refined_not_annotated_cluster_cons.fasta",
        clu_annot = config["rdir"] + "/annot_and_clust/annotated_clusters.tsv",
        pfam = config["pfam_shared_terms"],
        db_aln = config["rdir"] + \
            "/cluster_classification/ref_noannot_clu_aln",
        db_cons = config["rdir"] + \
            "/cluster_classification/ref_noannot_clu_cons",
        outdir = config["rdir"] + "/cluster_classification"
    output:
        k = config["rdir"] + "/cluster_classification/k_ids.txt",
        kwp = config["rdir"] + "/cluster_classification/kwp_ids.txt",
        gu = config["rdir"] + "/cluster_classification/gu_ids.txt",
        eu = config["rdir"] + "/cluster_classification/eu_ids.txt"
    log:
        out = "logs/class_stdout.log",
        err = "logs/class_stderr.err"
    benchmark:
        "benchmarks/cluster_classification/class.tsv"
    shell:
        """
        set -e
        set -x

        export OMPI_MCA_btl=^openib
        export OMP_NUM_THREADS={threads}
        export OMP_PROC_BIND=FALSE

        # Extract consensus sequnces
        # Create alignments and retrieve the consensus sequences
        if [ ! -f {params.db_aln} ]; then
            {params.mpi_runner} {params.mmseqs_bin} apply {params.refdb_noannot} {params.db_aln} --threads {threads} \
                -- {params.famsa_bin} STDIN STDOUT 2> /dev/null
        fi

        if [ ! -f {params.db_cons} ]; then
            {params.mpi_runner} {params.mmseqs_bin} apply {params.db_aln} {params.db_cons} --threads {threads} \
                -- {params.consensus} {params.hhcons_bin}
        fi

        # Extract the consensus sequences as fasta
        sed -e 's/\\x0//' {params.db_cons} > {params.ref_cons}

        # Rename them with using the real cluster-names
        join -11 -21 <(sort -k1,1 {params.name_index}) \
        <({params.seqkit_bin} fx2tab {params.ref_cons} | sort -k1,1 ) > {params.outdir}/cons_tmpl
        # Export the renamed (second field correspond to the cluster name) consensus sequences in fasta format
        awk -vOFS='\\t' '{{print $2,$3}}' {params.outdir}/cons_tmpl | \
         {params.seqkit_bin} tab2fx -w 0 > {params.ref_cons}

        rm {params.db_aln}* {params.db_cons}* {params.outdir}/cons_tmpl

        # Search the not annotated cluster consensus sequences against the UniRef90 database
        {params.mmseqs_search} --search {params.mmseqs_bin} \
                               --mpi_runner "{params.mpi_runner}" \
                               --ltmp {params.mmseqs_local_tmp} \
                               --cons {params.ref_cons} \
                               --db_fasta {params.uniref_fasta} \
                               --db_info {params.uniref_prot} \
                               --evalue_filter {params.evalue_filt} \
                               --evalue_threshold {params.evalue_thr} \
                               --hypo_threshold {params.hypo_thr} \
                               --hypo_patterns {params.patterns} \
                               --grep {params.ripgrep} \
                               --output {params.uniref_res} \
                               --outdir {params.outdir} \
                               --threads {threads} 2>{log.err} 1>{log.out}

        # The output of the search are 2 files: the MMseqs2 output in tabular format: "noannot_vs_uniref90.tsv"
        # the list of hypothetical proteins: "noannot_vs_uniref90_hypo.txt"
        # Parse results and search nohits against NCBI nr
        awk '!seen[$1]++{{print $1}}' {params.uniref_res} > {params.outdir}/noannot_vs_uniref90_hits.txt

        {params.filterbyname} in={params.ref_cons} \
                              out={params.outdir}/noannot_vs_uniref90_nohits.fasta \
                              names={params.outdir}/noannot_vs_uniref90_hits.txt \
                              include=f ignorejunk overwrite=true 2>>{log.err} 1>>{log.out}
        ####
        # Search the non-matching consensus sequences against the NCBI nr database

        {params.mmseqs_search} --search {params.mmseqs_bin} \
                                --mpi_runner "{params.mpi_runner}" \
                                --ltmp {params.mmseqs_local_tmp} \
                                --cons {params.outdir}/noannot_vs_uniref90_nohits.fasta \
                                --db_fasta {params.nr_fasta} \
                                --db_info {params.nr_prot} \
                                --evalue_filter {params.evalue_filt} \
                                --evalue_threshold {params.evalue_thr} \
                                --hypo_threshold {params.hypo_thr} \
                                --hypo_patterns {params.patterns} \
                                --grep {params.ripgrep} \
                                --output {params.nr_res} \
                                --outdir {params.outdir} \
                                --threads {threads} 2>{log.err} 1>{log.out}

        # The output of the search are 2 files: the MMseqs2 output in tabular format: "uniref-nohits_vs_NR.tsv"
        # the list of hypothetical proteins: "uniref-nohits_vs_NR_hypo.txt"

        # Parse results and define the first categories
        awk '!seen[$1]++{{print $1}}' {params.nr_res} \
        > {params.outdir}/uniref-nohits_vs_NR_hits.txt

        {params.filterbyname} in={params.outdir}/noannot_vs_uniref90_nohits.fasta \
                        out={params.outdir}/uniref-nohits_vs_NR_nohits.fasta \
                        names={params.outdir}/uniref-nohits_vs_NR_hits.txt \
                        include=f ignorejunk overwrite=true 2>>{log.err} 1>>{log.out}

        # Environmental unknowns (EUs)
        grep '^>' {params.outdir}/uniref-nohits_vs_NR_nohits.fasta | sed 's/^>//' > {output.eu}

        # Knowns without Pfam (KWPs)
        # Not-hypothetical hits from the search against UniRef90
        join -11 -21 -v1 <(sort {params.outdir}/noannot_vs_uniref90_hits.txt) \
          <(sort {params.outdir}/noannot_vs_uniref90_hypo.txt) > {params.outdir}/noannot_vs_uniref90_char.txt
        # Not-hypothetical hits from the search against NCBI nr
        join -11 -21 -v1 <(sort {params.outdir}/uniref-nohits_vs_NR_hits.txt) \
          <(sort {params.outdir}/uniref-nohits_vs_NR_hypo.txt) > {params.outdir}/uniref-nohits_vs_NR_char.txt

        cat {params.outdir}/noannot_vs_uniref90_char.txt \
        {params.outdir}/uniref-nohits_vs_NR_char.txt > {output.kwp}

        # Add annotation info:
        cat {params.outdir}/noannot_vs_uniref90_E{params.evalue_thr}_prot.tsv \
            {params.outdir}/uniref-nohits_vs_NR_E{params.evalue_thr}_prot.tsv > {params.outdir}/noannot_uniref_nr_annotations.tsv
        ## KWP annotations
        join -11 -21 <(sort -k1,1 {output.kwp} ) \
            <(sort -k1,1 {params.outdir}/noannot_uniref_nr_annotations.tsv) > {params.outdir}/kwp_annotations.tsv

        # Knowns (Ks) and Genomic unknowns (GUs)

        # Load the annotated cluster file in R to retrieve a consensus pfam domain architecture
        {params.da_r} --ref_annot {input.ref_annot} \
                            --clu_annot {params.clu_annot} \
                            --good_cl {params.good_cl} \
                            --pfam_terms {params.pfam} \
                            --dom_arch {params.dom_arch} \
                            --threads {threads} 2>>{log.err} 1>>{log.out}

        # Extract Pfam domain functional information (PF or DUF) from the DA file: {params.dom_arch}

        # Add the clusters annotated to pfam domains of unknown function (DUFs) to the GUs set
        cat <(awk 'NR>1 && $9=="DUF"{{print $1}}' {params.dom_arch}) \
         {params.outdir}/noannot_vs_uniref90_hypo.txt \
         {params.outdir}/uniref-nohits_vs_NR_hypo.txt > {output.gu}

        ## GU annotations
        join -11 -21 <(sort -k1,1 {output.gu} ) \
            <(sort -k1,1 {params.outdir}/noannot_uniref_nr_annotations.tsv) > {params.outdir}/gu_annotations.tsv
        awk 'NR>1 && $9=="DUF"{{print $1,"PFAM","0.0",$5}}' {params.dom_arch} >> {params.outdir}/gu_annotations.tsv

        # Retrieve the Knowns cluster set
        awk 'NR>1 && $9!="DUF"{{print $1}}' {params.dom_arch} > {output.k}

        ## K annotations
        awk -vOFS='\\t' 'NR>1 && $9!="DUF"{{print $1,"PFAM","0.0",$5}}' {params.dom_arch} >> {params.outdir}/k_annotations.tsv

        """

rule cluster_classific_done:
    input:
        k = config["rdir"] + "/cluster_classification/k_ids.txt",
        kwp = config["rdir"] + "/cluster_classification/kwp_ids.txt",
        gu = config["rdir"] + "/cluster_classification/gu_ids.txt",
        eu = config["rdir"] + "/cluster_classification/eu_ids.txt"
    output:
        ref_done = touch(config["rdir"] + "/cluster_classification/class.done")
    run:
        shell("echo 'CLUSTER CLASSIFICATION DONE'")
