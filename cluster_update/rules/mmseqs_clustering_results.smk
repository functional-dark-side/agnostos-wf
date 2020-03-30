rule mmseqs_clustering_results:
    input:
        clu = config["rdir"] + "/mmseqs_clustering/cluDB.tsv",
        all_orfsfa = config["rdir"] + "/gene_prediction/orf_seqs.fasta",
    threads: config['threads_res']
    params:
        mmseqs_bin = config["mmseqs_bin"],
        mmseqs_tmp = config["mmseqs_tmp"],
        local_tmp = config["mmseqs_local_tmp"],
        mpi_runner = config["mpi_runner"],
        seqtk_bin = config["seqtk_bin"],
        seqdb = config["rdir"] + "/mmseqs_clustering/seqDB",
        cludb = config["rdir"] + "/mmseqs_clustering/cluDB",
        cluseqdb = config["rdir"] + "/mmseqs_clustering/clu_seqDB",
        clean_seqs = config["wdir"] + "/scripts/clean_aa_seqs.sh",
        new_cluseqdb = config["rdir"] + "/mmseqs_clustering/new_clu_seqDB",
        wide = config["rdir"] + "/mmseqs_clustering/cluDB_wide.tsv",
        length = config["rdir"] + "/mmseqs_clustering/cluDB_length.tsv",
        orig_info1 = config["ordir"] + "/mmseqs_clustering/cluDB_name_rep_size.tsv",
        new_info1 = config["rdir"] + "/mmseqs_clustering/cluDB_new_name_rep_size.tsv",
        info1 = config["rdir"] + "/mmseqs_clustering/cluDB_all_name_rep_size.tsv",
        info2 = config["rdir"] + "/mmseqs_clustering/cluDB_name_rep_size.tsv",
        tmp = config["rdir"] + "/mmseqs_clustering/cluDB_tmpl",
        tmp1 = config["rdir"] + "/mmseqs_clustering/cluDB_tmpl1",
        awk_wide = config["wdir"] + "/scripts/convert_long_wide.awk",
        awk_row_col = config["wdir"] + "/scripts/row_columns.awk",
        rename_or = config["rdir"] + "/mmseqs_clustering/clu_rename_or",
        rename_new = config["rdir"] + "/mmseqs_clustering/clu_rename_new",
        or_new_info = config["rdir"] + "/mmseqs_clustering/cluDB_name_rep_size_old_new.tsv",
        naming = config["wdir"] + "/scripts/cluster_naming.sh",
        namedb = config["rdir"] + "/mmseqs_clustering/cluDB_ids.db",
        shared = config["rdir"] + "/mmseqs_clustering/cluDB_shared_name_rep_size.tsv",
        original = config["rdir"] + "/mmseqs_clustering/cluDB_original_name_rep_size.tsv",
        index = config["rdir"] + "/mmseqs_clustering/cluDB_name_index.txt"
    conda:
        "../envs/workflow.yml"
    log:
        out = "logs/mmseqs_clustering_res_stdout.log",
        err = "logs/mmseqs_clustering_res_stderr.err"
    benchmark:
        "benchmarks/mmseqs_clustering/clu_res.tsv"
    output:
        clusters = config["rdir"] + \
            "/mmseqs_clustering/cluDB_no_singletons.tsv",
        singl = config["rdir"] + "/mmseqs_clustering/cluDB_singletons.tsv",
        clu_info = config["rdir"] + "/mmseqs_clustering/cluDB_info.tsv",
        new_index = config["rdir"] + "/mmseqs_clustering/new_cluDB_name_index.txt"
    shell:
        """
        set -e
        set -x

        # Create sequences database for the clustering
        # This DB can be accessed with ffindex, to extractt separated fasta files for each cluster and perform operations on them
        # Ex: ffindex_apply_mpi cluseqdb cluseqdb.index -- your_program/script
        if [[ ! -s {params.cluseqdb}_orig.index ]]; then
            {params.mmseqs_bin} createseqfiledb {params.seqdb} {params.cludb} {params.cluseqdb}_orig 2>{log.err}
        fi

        #Check if sequences contain "*", and if yes, clean them

        match=$( grep "\*" {params.seqdb} | head -n1 )

        if [ -n "$match" ]; then
            {params.mpi_runner} {params.mmseqs_bin} apply \
                {params.cluseqdb}_orig \
                {params.cluseqdb}  \
                --threads {threads} \
                -- {params.clean_seqs} 2>>{log.err}
        fi

        # To convert this cluster results tab separated file in wide format (repres member member member ..)
        awk -f {params.awk_wide} {input.clu} > {params.wide} 2>>{log.err}

        # Cluster name rep and number of ORFs (size)
        awk -vFS='\\t' -vOFS='\\t' '{{print NR,$1,NF-1}}' {params.wide} > {params.info1}

        # Join with the original clustering using the representative sequences
        join -12 -22 <(sort -k2,2 {params.orig_info1}) \
          <(sort -k2,2 {params.info1}) > {params.rename_or} 2>>{log.err}

        join -12 -22 -v2 <(sort -k2,2 {params.orig_info1}) \
          <(sort -k2,2 {params.info1}) > {params.rename_new} 2>>{log.err}

        # Store the information about changes in the original Clusters
        sort -k2,2n {params.rename_or} | awk -vOFS='\\t' '{{print $2,$1,$3,$5}}' > {params.or_new_info}
        # Order and rejoin old and new cluster ids/NAMES (the new ids will start from the last of the original clustering)
        OR=$(wc -l {params.or_new_info} | cut -d ' ' -f1)

        awk -vOR=$OR -vOFS='\\t' '{{print NR+OR,$1,$3}}' {params.rename_new} > {params.new_info1} 2>>{log.err}

        cat <(awk -vOFS='\\t' '{{print $1,$2,$4}}' {params.or_new_info} ) \
          {params.new_info1} > {params.info1}

        # Cluster naming
        # The official cluster names are going to be based on the line number of the wide formatted file
        # We are also going to produce a correspondence file to access the clusters in the MMseqs2 indices
        # Retrieving MMseqs-entry-name and sequence header from the clu_seqDB
        if [ -s {params.cluseqdb} ]; then
            {params.mpi_runner} {params.mmseqs_bin} apply {params.cluseqdb} {params.namedb} --threads {threads} -- {params.naming}
        else
            {params.mpi_runner} {params.mmseqs_bin} apply {params.cluseqdb}_orig {params.namedb} --threads {threads} -- {params.naming}
        fi
        # Join the sequences with the representatives: this way we combine mmseqs-name with cluster-name (row number)
        join -12 -22 <(sed -e 's/\\x0//' {params.namedb} | sort -k2,2 --parallel={threads} -T {params.local_tmp}) \
        <(sort -k2,2 --parallel={threads} -T {params.mmseqs_tmp} {params.info1} ) > {params.index} 2>>{log.err}

        awk -vOFS='\\t' '{{print $2,$3}}' {params.index} > {params.tmp} && mv {params.tmp} {params.index}

        # Join with the long format cluster file
        join -11 -22 <(sort --parallel={threads} -T {params.local_tmp} -k1,1 {input.clu} ) \
            <(sort --parallel={threads} -k2,2 {params.info1} ) > {params.tmp} 2>>{log.err}

        # Add length info
        sed -e 's/\\x0//g' {params.cluseqdb} | {params.seqtk_bin} comp | cut -f1,2 > {params.length}

        join -11 -22 <(sort -k1,1 --parallel={threads} -T {params.local_tmp} {params.length}) \
         <(sort -k2,2 --parallel={threads} -T {params.mmseqs_tmp} {params.tmp}) > {output.clu_info} 2>>{log.err}

        # Reorder fields (cl_name rep orf length size)
        sort -k4,4n {output.clu_info} | awk -vOFS='\\t' '{{print $4,$3,$1,$2,$5}}' > {params.tmp}

        cp {params.tmp} {output.clu_info}

        ## Retrieve the different cluster sets:
        # 1. Only-original clusters: cluDB_original_name_rep_size.tsv

        awk -vOFS='\\t' '$3==$4{{print $1,$2,$3}}' {params.or_new_info} > {params.original}

        # 2. Shared original-new clusters: cluDB_shared_name_rep_size.tsv

        awk -vOFS='\\t' '$3!=$4{{print $1,$2,$4}}' {params.or_new_info} > {params.shared}

        # 3. Only-new clusters: cluDB_new_name_rep_size.tsv
        # we already have this set: {params.new_info1}

        ## We are going to re-validate and classified the clusters in new_cluDB and
        # the shared_cluDB previously singletons now with size > 1

        cat {params.new_info1} \
        <(awk -vFS='\\t' -vOFS='\\t' '$3==1 && $4>1 {{print $1,$2,$4}}' {params.or_new_info} ) > {params.info2}

        # Subset cluster-index-name:

        join -11 -22 <(awk '{{print $1}}' {params.info2} | sort -k1,1 ) \
        <(sort -k2,2 {params.index}) > {params.tmp}

        awk -vOFS='\\t' '{{print $2,$1}}' {params.tmp} > {output.new_index}

        # And then subset the cluseqDB:
        if [ -s {params.cluseqdb} ]; then
            {params.mmseqs_bin} createsubdb <(awk '{{print $1}}' {output.new_index}) \
                                            {params.cluseqdb} \
                                            {params.new_cluseqdb} 2>>{log.err}
        else
            {params.mmseqs_bin} createsubdb <(awk '{{print $1}}' {output.new_index}) \
                                            {params.cluseqdb}_orig \
                                            {params.new_cluseqdb} 2>>{log.err}
        fi

        # Subset the clusters based on their size
        # Clusters with more than 1 member
        awk -vOFS='\\t' '$3>=2{{print $1,$2}}' {params.info2} > {params.tmp1}
        join -12 -21 <(sort -k2,2 {params.tmp1}) \
        <( awk -vFS='\\t' -vOFS='\\t' '{{print $2,$3}}' {output.clu_info} | sort -k1,1) > {params.tmp} 2>>{log.err}

        awk -vOFS='\\t' '{{print $2,$1,$3}}' {params.tmp} > {output.clusters}

        # Singletons
        awk -vOFS='\\t' '$3=="1"{{print $1,$2}}' {params.info2} > {output.singl}

        rm {params.tmp} {params.tmp1} {params.namedb} {params.namedb}.index {params.namedb}.dbtype

        """

rule mmseqs_clustering_results_done:
    input:
        clusters = config["rdir"] + \
            "/mmseqs_clustering/cluDB_no_singletons.tsv",
        singl = config["rdir"] + "/mmseqs_clustering/cluDB_singletons.tsv",
        clu_info = config["rdir"] + "/mmseqs_clustering/cluDB_info.tsv",
        new_index = config["rdir"] + "/mmseqs_clustering/new_cluDB_name_index.txt"
    output:
        cls_done = touch(config["rdir"] + "/mmseqs_clustering/clu_res.done")
    run:
        shell("echo 'CLUSTERING RESULTS DONE'")
