rule mmseqs_clustering_updt_results:
    input:
        clu = config["rdir"] + "/mmseqs_clustering/cluDB.tsv"
    threads: config['threads_res']
    params:
        mmseqs_bin = config["mmseqs_bin"],
        mmseqs_tmp = config["mmseqs_tmp"],
        local_tmp = config["mmseqs_local_tmp"],
        mpi_runner = config["mpi_runner"],
        seqtk_bin = config["seqtk_bin"],
        eval_shared = config["eval_shared"],
        shared_type = config["shared_type"],
        clean_seqs = config["wdir"] + "/scripts/clean_aa_seqs.sh",
        awk_wide = config["wdir"] + "/scripts/convert_long_wide.awk",
        awk_row_col = config["wdir"] + "/scripts/row_columns.awk",
        naming = config["wdir"] + "/scripts/cluster_naming.sh",
        concat_clu = config["wdir"] + "/scripts/concatenate_cludb.sh",
        seqdb = config["rdir"] + "/mmseqs_clustering/seqDB",
        cludb = config["rdir"] + "/mmseqs_clustering/cluDB",
        cluseqdb = config["rdir"] + "/mmseqs_clustering/clu_seqDB",
        new_cluseqdb = config["rdir"] + "/mmseqs_clustering/new_clu_seqDB",
        wide = config["rdir"] + "/mmseqs_clustering/cluDB_wide.tsv",
        length = config["rdir"] + "/mmseqs_clustering/cluDB_length.tsv",
        orig_info1 = config["ordir"] + "/mmseqs_clustering/cluDB_name_rep_size.tsv",
        new_info1 = config["rdir"] + "/mmseqs_clustering/cluDB_new_name_rep_size.tsv",
        info1 = config["rdir"] + "/mmseqs_clustering/cluDB_all_name_rep_size.tsv",
        info2 = config["rdir"] + "/mmseqs_clustering/cluDB_name_rep_size.tsv",
        tmp = config["rdir"] + "/mmseqs_clustering/cluDB_tmpl",
        tmp1 = config["rdir"] + "/mmseqs_clustering/cluDB_tmpl1",
        rename_or = config["rdir"] + "/mmseqs_clustering/clu_rename_or",
        rename_new = config["rdir"] + "/mmseqs_clustering/clu_rename_new",
        or_new_info = config["rdir"] + "/mmseqs_clustering/cluDB_name_rep_size_old_new.tsv",
        or_clu_cat = config["ordir"] + "/cluster_ids_categ.tsv.gz",
        namedb = config["rdir"] + "/mmseqs_clustering/cluDB_ids.db",
        shared = config["rdir"] + "/mmseqs_clustering/cluDB_shared_name_rep_size.tsv",
        original = config["rdir"] + "/mmseqs_clustering/cluDB_original_name_rep_size.tsv",
        index = config["rdir"] + "/mmseqs_clustering/cluDB_name_index.txt",
        clu_info = config["rdir"] + "/mmseqs_clustering/cluDB_info.tsv"
    conda:
        config["conda_env"]
    log:
        out = "logs/mmseqs_clustering_res_stdout.log",
        err = "logs/mmseqs_clustering_res_stderr.err"
    benchmark:
        "benchmarks/clu_updt_res.tsv"
    output:
        clusters = config["rdir"] + \
            "/mmseqs_clustering/cluDB_no_singletons.tsv",
        singl = config["rdir"] + "/mmseqs_clustering/cluDB_singletons.tsv"
    shell:
        """
        set -e
        set -x

        mkdir -p {params.mmseqs_tmp}

        # Create sequences database for the clustering
        # This DB can be accessed with ffindex, to extractt separated fasta files for each cluster and perform operations on them
        # Ex: ffindex_apply_mpi cluseqdb cluseqdb.index -- your_program/script
        if [[ ! -s {params.cluseqdb}.index ]]; then

            {params.mmseqs_bin} createseqfiledb {params.seqdb} {params.cludb} {params.cluseqdb} --threads {threads} 2>{log.err}

            #Check if sequences contain "*", and if yes, clean them using the following command
            #{params.mpi_runner} {params.mmseqs_bin} apply \
            #    {params.cluseqdb} \
            #    {params.cluseqdb}_clean  \
            #    --threads {threads} \
            #    -- {params.clean_seqs} 2>>{log.err}
        fi

        if [[ -s {params.cludb}.0 ]]; then
            {params.concat_clu} {params.cludb} {threads}
        fi

        if [[ -s {params.cluseqdb}.0 ]]; then
            {params.concat_clu} {params.cluseqdb} {threads}
        fi

        # To convert this cluster results tab separated file in wide format (repres member member member ..)
        awk -f {params.awk_wide} <(awk '!seen[$2]++' {input.clu}) > {params.wide} 2>>{log.err}

        # Cluster name rep and number of ORFs (size)
        awk -vFS='\\t' -vOFS='\\t' '{{print NR,$1,NF-1}}' {params.wide} > {params.info1}


        # Join with the original clustering using the representative sequences
        join -12 -22 <(sort -k2,2 {params.orig_info1}) \
          <(sort -k2,2 {params.info1}) > {params.rename_or} 2>>{log.err}

        join -12 -22 -v2 <(sort -k2,2 {params.orig_info1}) \
          <(sort -k2,2 {params.info1}) > {params.rename_new} 2>>{log.err}

        # Store the information about changes in the original Clusters
        sort -k2,2n {params.rename_or} | awk -vOFS='\\t' '{{print $2,$1,$3,$5}}' > {params.or_new_info}
        # Order and re-join old and new cluster ids/NAMES (the new ids will start from the last of the original clustering)
        OR=$(wc -l {params.or_new_info} | cut -d ' ' -f1)

        awk -vOR=$OR -vOFS='\\t' '!seen[$0]++{{print NR+OR,$1,$3}}' {params.rename_new} > {params.new_info1} 2>>{log.err}

        cat <(awk -vOFS='\\t' '!seen[$0]++{{print $1,$2,$4}}' {params.or_new_info} ) \
          {params.new_info1} > {params.info1}

        # Cluster naming
        # The official cluster names are going to be based on the line number of the wide formatted file
        # We are also going to produce a correspondence file to access the clusters in the MMseqs2 indices
        # Retrieving MMseqs-entry-name and sequence header from the clu_seqDB
        if [[ ! -s {params.index} ]]; then
            {params.mpi_runner} {params.mmseqs_bin} apply {params.cluseqdb} {params.namedb} --threads {threads} -- {params.naming}

            # Join the sequences with the representatives: this way we combine mmseqs-name with cluster-name (row number)
            join -12 -22 <(sed -e 's/\\x0//g' {params.namedb} | sort -k2,2 --parallel={threads} -T {params.local_tmp}) \
                <(sort -k2,2 --parallel={threads} -T {params.mmseqs_tmp} {params.info1} ) > {params.index}
            awk -vOFS='\\t' '!seen[$0]++{{print $2,$3}}' {params.index} > {params.tmp} && mv {params.tmp} {params.index}
        fi

        # Rename the cluster sequence DB with the cluster names that are goiung to be used in the next steps
        {params.mmseqs_bin} renamedbkeys {params.index} {params.cluseqdb} {params.cluseqdb}.new

        mv {params.cluseqdb}.new {params.cluseqdb}
        mv {params.cluseqdb}.new.index {params.cluseqdb}.index
        mv {params.cluseqdb}.new.dbtype {params.cluseqdb}.dbtype

        # Join with the long format cluster file
        join -11 -22 <(sort --parallel={threads} -T {params.local_tmp} -k1,1 {input.clu}) \
            <(sort --parallel={threads} -T {params.local_tmp} -k2,2 {params.info1} ) > {params.tmp1}

        # Add length info
        sed -e 's/\\x0//g' {params.cluseqdb} | {params.seqtk_bin} comp | cut -f1,2 > {params.length}

        join -11 -22 <(sort -k1,1 --parallel={threads} -T {params.local_tmp} {params.length}) \
         <(sort -k2,2 --parallel={threads} -T {params.mmseqs_tmp} {params.tmp1}) > {params.tmp}

        # Reorder fields (cl_name rep orf length size)
        sort -k4,4n --parallel={threads} -T {params.local_tmp} {params.tmp} | \
         awk -vOFS='\\t' '{{print $4,$3,$1,$2,$5}}' > {params.clu_info}

        ## Retrieve the different cluster sets:
        # 1. Only-original clusters: cluDB_original_name_rep_size.tsv

        awk -vOFS='\\t' '$3==$4{{print $1,$2,$3}}' {params.or_new_info} | awk '!seen[$0]++' > {params.original}

        # 2. Shared original-new clusters: cluDB_shared_name_rep_size.tsv

        awk -vOFS='\\t' '$3!=$4{{print $1,$2,$4}}' {params.or_new_info} | awk '!seen[$0]++' > {params.shared}

        # 3. Only-new clusters: cluDB_new_name_rep_size.tsv
        # we already have this set: {params.new_info1}

        ## We are going to re-validate and classified the clusters in new_cluDB and

        if [[ {params.eval_shared} == "true" ]]; then

            if [[ {params.shared_type} == "all" ]]; then
                # Evaluate all the shared GCs (GCs that got integrated with new genes)
                cat {params.new_info1} \
                    <(awk -vFS='\\t' -vOFS='\\t' '!seen[$0]++{{print $1,$2,$3}}' {params.shared} ) > {params.info2}
            fi

            if [[ {params.shared_type} == "discarded" ]]; then
                # Add the previously discarded shared GCs
                # download if not there the original cluster - category table
                if [[ ! -s {params.or_clu_cat} ]]; then
                    wget https://ndownloader.figshare.com/files/23067140 -O {params.or_clu_cat}
                fi
                # anti join shared GCs with the original cluster-category to get the discarded GCs
                join -11 -21 -v1 <(sort -k1,1 --parallel={threads} -T {params.local_tmp} {params.shared}) \
                    <(zcat {params.or_clu_cat} | awk '{{print $1}}' |\
                     sort -k1,1 --parallel={threads} -T {params.local_tmp}) > {params.shared}.temp

                cat {params.new_info1} \
                    <(awk -vFS='\\t' -vOFS='\\t' '!seen[$0]++{{print $1,$2,$3}}' {params.shared}.temp ) > {params.info2}
                rm {params.shared}.temp
            fi

            if [[ {params.shared_type} == "increase" ]]; then
                # Or the shared GCs with more than % new genes
                # Ex: increase of 30% in number of sequences
                awk '{{print $1,$2,$4,($4-$3)/$3}}' {params.or_new_info} > {params.or_new_info}.temp
                cat {params.new_info1} \
                    <(awk -vFS='\\t' -vOFS='\\t' '!seen[$0]++ && $4>=0.3{{print $1,$2,$3}}' {params.or_new_info}.temp ) > {params.info2}
                rm {params.or_new_info}.temp
            fi
        else
            # Evaluate also the previously singletons now integrated in GCs with new genes
            cat {params.new_info1} \
                <(awk -vFS='\\t' -vOFS='\\t' '$3==1 && $4>1 && !seen[$0]++{{print $1,$2,$4}}' {params.or_new_info} ) > {params.info2}
        fi

        # Subset the cluseqDB:
        if [[ ! -s {params.new_cluseqdb} ]]; then
            {params.mmseqs_bin} createsubdb <(awk '{{print $1}}' {params.info2}) \
                                            {params.cluseqdb} \
                                            {params.new_cluseqdb} 2>>{log.err}
        fi

        # Subset the clusters based on their size
        # Clusters with more than 1 member
        awk -vOFS='\\t' '$3>=2{{print $1,$2}}' {params.info2} > {params.tmp1}
        join -12 -21 <(sort -k2,2 {params.tmp1}) \
        <(sort -k1,1 --parallel={threads} -T {params.local_tmp} <(awk '!seen[$2]++' {input.clu} )) > {params.tmp}

        awk -vOFS='\\t' '!seen[$0]++{{print $2,$1,$3}}' {params.tmp} > {output.clusters}

        # Singletons
        awk -vOFS='\\t' '$3=="1"{{print $1,$2}}' {params.info2} | awk '!seen[$0]++' > {output.singl}

        rm {params.tmp} {params.tmp1} {params.namedb} {params.namedb}.index {params.namedb}.dbtype
        rm {params.wide} {params.length} {params.rename_or} {params.rename_new}

        """

rule mmseqs_clustering_results_done:
    input:
        clusters = config["rdir"] + \
            "/mmseqs_clustering/cluDB_no_singletons.tsv",
        singl = config["rdir"] + "/mmseqs_clustering/cluDB_singletons.tsv"
    output:
        cls_done = touch(config["rdir"] + "/mmseqs_clustering/clu_res.done")
    run:
        shell("echo 'CLUSTERING RESULTS DONE'")
