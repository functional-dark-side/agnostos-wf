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
        vmtouch = config["vmtouch"],
        mmseqs_bin = config["mmseqs_bin"],
        local_tmp = config["mmseqs_local_tmp"],
        mpi_runner = config["mpi_runner"],
        filterbyname = config["filterbyname"],
        famsa_bin = config["famsa_bin"],
        hhcons_bin = config["hhcons_bin"],
        seqkit_bin = config["seqkit_bin"],
        singl = config["singl"],
        uniref_db = config["uniref90_db"],
        uniref_prot = config["uniref90_prot"],
        nr_db = config["nr_db"],
        nr_prot = config["nr_prot"],
        pfam = config["pfam_shared_terms"],
        db_mode = config["db_mode"],
        evalue_thr = 0.6,
        hypo_thr = 1.0,
        consensus = "scripts/consensus.sh",
        mmseqs_search = "scripts/mmseqs_double_search.sh",
        evalue_filt = "scripts/evalue_filter.awk",
        ripgrep = "rg",
        patterns = "scripts/hypothetical_grep.tsv",
        da_r = "scripts/cluster_pfam_domain_architectures.r",
        good_cl = config["rdir"] + "/validation/good_clusters.tsv",
        refdb_noannot = config["rdir"] + \
            "/cluster_refinement/refined_not_annotated_clusterDB",
        genes = config["rdir"] + "/mmseqs_clustering/clu_seqDB",
        s_annot = config["rdir"] + "/annot_and_clust/singletons_pfam_annot.tsv",
        s_all = config["rdir"] + "/mmseqs_clustering/cluDB_singletons.tsv",
        clu_annot = config["rdir"] + "/annot_and_clust/annotated_clusters.tsv",
        outdir = config["rdir"] + "/cluster_classification",
        db_aln = config["rdir"] + \
            "/cluster_classification/ref_noannot_clu_aln",
        db_cons = config["rdir"] + \
            "/cluster_classification/ref_noannot_clu_cons",
        ref_cons = config["rdir"] + \
            "/cluster_classification/refined_not_annotated_cluster_cons.fasta",
        uniref_res = config["rdir"] + \
            "/cluster_classification/noannot_vs_uniref90.tsv",
        nr_res = config["rdir"] + \
            "/cluster_classification/uniref-nohits_vs_NR.tsv",
        dom_arch = config["rdir"] + \
            "/cluster_classification/cluster_pfam_domain_architectures.tsv",
        s_categ = config["rdir"] + "/cluster_classification/singleton_gene_cl_categories.tsv"
    output:
        k = config["rdir"] + "/cluster_classification/k_ids.txt",
        kwp = config["rdir"] + "/cluster_classification/kwp_ids.txt",
        gu = config["rdir"] + "/cluster_classification/gu_ids.txt",
        eu = config["rdir"] + "/cluster_classification/eu_ids.txt"
    log:
        out = "logs/class_stdout.log",
        err = "logs/class_stderr.err"
    benchmark:
        "benchmarks/clu_class.tsv"
    shell:
        """
        set -e
        set -x

        export OMPI_MCA_btl=^openib
        export OMP_NUM_THREADS={threads}
        export OMP_PROC_BIND=FALSE

        # Extract consensus sequnces
        if [ ! -f {params.ref_cons} ]; then
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
            sed -e 's/\\x0//g' {params.db_cons} > {params.ref_cons}

            rm {params.db_aln}* {params.db_cons}*

            # Add singletons not annotated here to ref_cons if config["singl"]=="true"
            # keep the name ref_cons for the two searches (easier)
            if [ {params.singl} == "true" ]; then
                join -11 -21 -v1 <(awk '{{print $2}}' {params.s_all} | sort -k1,1) \
                <(awk '{{print $1}}' {params.s_annot} | sort -k1,1) > {params.outdir}/singletons_not_annot.txt
                sed -e 's/\\x0//g' {params.genes} > {params.outdir}/genes.fasta
                {params.filterbyname} in={params.outdir}/genes.fasta \
                                    out={params.outdir}/singletons_not_annot.fasta \
                                    names={params.outdir}/singletons_not_annot.txt \
                                    include=t ignorejunk overwrite="true" 2>>{log.err} 1>>{log.out}
                cat {params.outdir}/singletons_not_annot.fasta >> {params.ref_cons}
                rm {params.outdir}/genes.fasta
            fi
        fi

        # Search the not annotated cluster consensus sequences against the UniRef90 database

        if [ ! -s {params.uniref_prot} ]; then
            echo "Dowloading UniRef90 DB using mmseqs"
            {params.mmseqs_bin} databases UniRef90 {params.uniref_db} {params.local_tmp} --threads {threads} --remove-tmp-files
            # Create the protein description file:
            echo "extracting protein information"
            sed 's/\\x0//g' {params.uniref_db}_h | sed 's/ /__/g' | sed 's/__/\\t/' | sed 's/__/_/g' | gzip > {params.uniref_prot}
        fi

        if [ ! -f {params.outdir}/noannot_vs_uniref90_hits.txt ]; then
            {params.vmtouch} -f {params.uniref_db}*
            {params.mmseqs_search} --search {params.mmseqs_bin} \
                                    --mpi_runner "{params.mpi_runner}" \
                                    --ltmp {params.local_tmp} \
                                    --cons {params.ref_cons} \
                                    --db_target {params.uniref_db} \
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
                              include=f ignorejunk overwrite="true" 2>>{log.err} 1>>{log.out}

            if [[ {params.db_mode} == "memory" ]]; then
                rm {params.uniref_db}* {params.uniref_prot}
            fi

        fi
        ####
        # Search the non-matching consensus sequences against the NCBI nr database

        if [ ! -s {params.nr_prot} ]; then
            echo "Dowloading NR DB using mmseqs"
            {params.mmseqs_bin} databases NR {params.nr_db} {params.local_tmp} --threads {threads} --remove-tmp-files
            # Create the protein description file:
            echo "extracting protein information"
            sed 's/\\x0//g' {params.nr_db}_h | sed 's/ /__/g' | sed 's/__/\\t/' | sed 's/__/_/g' | gzip > {params.nr_prot}
        fi

        if [ ! -f {params.outdir}/uniref-nohits_vs_NR_hits.txt ]; then
            {params.vmtouch} -f {params.nr_db}*
            {params.mmseqs_search} --search {params.mmseqs_bin} \
                                    --mpi_runner "{params.mpi_runner}" \
                                    --ltmp {params.local_tmp} \
                                    --cons {params.outdir}/noannot_vs_uniref90_nohits.fasta \
                                    --db_target {params.nr_db} \
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
                        include=f ignorejunk overwrite="true" 2>>{log.err} 1>>{log.out}

            if [[ {params.db_mode} == "memory" ]]; then
                rm {params.nr_db}* {params.nr_prot}
            fi

        fi

        # Environmental unknowns (EUs)
        if [ {params.singl} == "true" ]; then
          grep '^>' {params.outdir}/uniref-nohits_vs_NR_nohits.fasta | sed 's/^>//' > {params.outdir}/all_eu
          # cluster EU
          join -11 -21 -v1 <(sort -k1,1 {params.outdir}/all_eu) \
          <(sort -k1,1 {params.outdir}/singletons_not_annot.txt) > {output.eu}
          # singleton EU
          join -11 -21 <(sort -k1,1 {params.outdir}/all_eu) \
          <(sort -k1,1 {params.outdir}/singletons_not_annot.txt) > {params.outdir}/singl_eu
        else
            grep '^>' {params.outdir}/uniref-nohits_vs_NR_nohits.fasta | sed 's/^>//' > {output.eu}
        fi

        # Knowns without Pfam (KWPs)
        # Not-hypothetical hits from the search against UniRef90
        join -11 -21 -v1 <(sort {params.outdir}/noannot_vs_uniref90_hits.txt) \
          <(sort {params.outdir}/noannot_vs_uniref90_hypo.txt) > {params.outdir}/noannot_vs_uniref90_char.txt
        # Not-hypothetical hits from the search against NCBI nr
        join -11 -21 -v1 <(sort {params.outdir}/uniref-nohits_vs_NR_hits.txt) \
          <(sort {params.outdir}/uniref-nohits_vs_NR_hypo.txt) > {params.outdir}/uniref-nohits_vs_NR_char.txt

        if [ {params.singl} == "true" ]; then
              cat {params.outdir}/noannot_vs_uniref90_char.txt \
              {params.outdir}/uniref-nohits_vs_NR_char.txt > {params.outdir}/all_kwp
              # cluster kWP
              join -11 -21 -v1 <(sort -k1,1 {params.outdir}/all_kwp) \
              <(sort -k1,1 {params.outdir}/singletons_not_annot.txt) > {output.kwp}
              # singleton KWP
              join -11 -21 <(sort -k1,1 {params.outdir}/all_kwp) \
              <(sort -k1,1 {params.outdir}/singletons_not_annot.txt) > {params.outdir}/singl_kwp
        else
              cat {params.outdir}/noannot_vs_uniref90_char.txt \
               {params.outdir}/uniref-nohits_vs_NR_char.txt > {output.kwp}
        fi

        # Add annotation info:
        cat {params.outdir}/noannot_vs_uniref90_E{params.evalue_thr}_prot.tsv \
            {params.outdir}/uniref-nohits_vs_NR_E{params.evalue_thr}_prot.tsv > {params.outdir}/noannot_uniref_nr_annotations.tsv
        ## KWP annotations
        join -11 -21 <(sort -k1,1 {output.kwp} ) \
            <(sort -k1,1 {params.outdir}/noannot_uniref_nr_annotations.tsv) > {params.outdir}/kwp_annotations.tsv

        # Knowns (Ks) and Genomic unknowns (GUs)

        if [ ! -s {params.pfam} ]; then
            echo "Dowloading Pfam list of shared domain names"
            wget https://figshare.com/ndownloader/files/31127782 -O {params.pfam}
        fi

        # Load the annotated cluster file in R to retrieve a consensus pfam domain architecture
        {params.da_r} --ref_annot {input.ref_annot} \
                            --clu_annot {params.clu_annot} \
                            --good_cl {params.good_cl} \
                            --pfam_terms {params.pfam} \
                            --dom_arch {params.dom_arch} \
                            --threads {threads} 2>>{log.err} 1>>{log.out}

        if [[ {params.db_mode} == "memory" ]]; then
            rm {params.pfam}
        fi

        # Extract Pfam domain functional information (PF or DUF) from the DA file: {params.dom_arch}

        # Add the clusters annotated to pfam domains of unknown function (DUFs) to the GUs set
        if [ {params.singl} == "true" ]; then
            cat <(awk 'NR>1 && $9=="DUF"{{print $1}}' {params.dom_arch}) \
                <(awk '$3~/^DUF/{{print $1}}' {params.s_annot}) \
                {params.outdir}/noannot_vs_uniref90_hypo.txt \
                {params.outdir}/uniref-nohits_vs_NR_hypo.txt > {params.outdir}/all_gu
             # cluster GU
             join -11 -21 -v1 <(sort -k1,1 {params.outdir}/all_gu) \
             <(awk '{{print $2}}' {params.s_all} | sort -k1,1 ) > {output.gu}
             # singleton GU
             join -11 -21 <(sort -k1,1 {params.outdir}/all_gu) \
             <(awk '{{print $2}}' {params.s_all} | sort -k1,1) > {params.outdir}/singl_gu
        else
            cat <(awk 'NR>1 && $9=="DUF"{{print $1}}' {params.dom_arch}) \
                {params.outdir}/noannot_vs_uniref90_hypo.txt \
                {params.outdir}/uniref-nohits_vs_NR_hypo.txt > {output.gu}
        fi

        ## GU annotations
        join -11 -21 <(sort -k1,1 {output.gu} ) \
            <(sort -k1,1 {params.outdir}/noannot_uniref_nr_annotations.tsv) > {params.outdir}/gu_annotations.tsv
        awk 'NR>1 && $9=="DUF"{{print $1,"PFAM","0.0",$6}}' {params.dom_arch} >> {params.outdir}/gu_annotations.tsv

        #rm {params.outdir}/noannot_uniref_nr_annotations.tsv

        # Retrieve the Knowns cluster set
        if [ {params.singl} == "true" ]; then
            cat <(awk 'NR>1 && $9!="DUF"{{print $1}}' {params.dom_arch}) \
                <(awk '$3!~/^DUF/{{print $1}}' {params.s_annot}) > {params.outdir}/all_k
             # cluster K
             join -11 -21 -v1 <(sort -k1,1 {params.outdir}/all_k) \
             <(awk '{{print $1}}' {params.s_annot} | sort -k1,1 ) > {output.k}
             # singleton K
             join -11 -21 <(sort -k1,1 {params.outdir}/all_k) \
             <(awk '{{print $1}}' {params.s_annot} | sort -k1,1) > {params.outdir}/singl_k
             # singleton categories
             cat <(awk '{{print $1,"EU"}}' {params.outdir}/singl_eu) \
                 <(awk '{{print $1,"GU"}}' {params.outdir}/singl_gu) \
                 <(awk '{{print $1,"KWP"}}' {params.outdir}/singl_kwp) \
                 <(awk '{{print $1,"K"}}' {params.outdir}/singl_k) > {params.outdir}/singl_gene_categ
            # singletons gene cluster categories
            join -12 -21 <(sort -k2,2 {params.s_all}) \
            <(sort -k1,1 {params.outdir}/singl_gene_categ) |\
             sed 's/ /\\t/g' > {params.s_categ}
            rm {params.outdir}/singl_* {params.outdir}/singletons_not_annot* {params.outdir}/all_*
        else
            awk 'NR>1 && $9!="DUF"{{print $1}}' {params.dom_arch} > {output.k}
        fi

        ## K annotations
        awk -vOFS='\\t' 'NR>1 && $9!="DUF"{{print $1,"PFAM","0.0",$6}}' {params.dom_arch} >> {params.outdir}/k_annotations.tsv

        # Clear results
        #rm -rf {params.outdir}/noannot_vs_uniref90_* {params.outdir}/uniref-nohits_vs_NR_*
        #rm -rf {params.refdb_noannot}*

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
