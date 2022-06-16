rule mmseqs_clustering_results:
    input:
        clu = config["rdir"] + "/mmseqs_clustering/cluDB.tsv"
    threads:
        config['threads_default']
    priority: 40
    params:
        mmseqs_bin = config["mmseqs_bin"],
        seqtk_bin = config["seqtk_bin"],
        mpi_runner = config["mpi_runner"],
        local_tmp = config["mmseqs_tmp"],
        awk_wide = config["wdir"] + "/scripts/convert_long_wide.awk",
        awk_row_col = config["wdir"] + "/scripts/row_columns.awk",
        concat_clu = config["wdir"] + "/scripts/concatenate_cludb.sh",
        naming = config["wdir"] + "/scripts/cluster_naming.sh",
        seqdb = config["rdir"] + "/mmseqs_clustering/seqDB",
        cludb = config["rdir"] + "/mmseqs_clustering/cluDB",
        cluseqdb = config["rdir"] + "/mmseqs_clustering/clu_seqDB",
        wide = config["rdir"] + "/mmseqs_clustering/cluDB_wide.tsv",
        length = config["rdir"] + "/mmseqs_clustering/cluDB_length.tsv",
        info1 = config["rdir"] + "/mmseqs_clustering/cluDB_name_rep_size.tsv",
        tmp = config["rdir"] + "/mmseqs_clustering/cluDB_tmpl",
        tmp1 = config["rdir"] + "/mmseqs_clustering/cluDB_tmpl1",
        namedb = config["rdir"] + "/mmseqs_clustering/cluDB_ids.db"
    conda:
        config["conda_env"]
    log:
        out = "logs/mmseqs_clustering_res_stdout.log",
        err = "logs/mmseqs_clustering_res_stderr.err"
    benchmark:
        "benchmarks/clu_res.tsv"
    output:
        clusters = config["rdir"] + \
            "/mmseqs_clustering/cluDB_no_singletons.tsv",
        singl = config["rdir"] + "/mmseqs_clustering/cluDB_singletons.tsv",
        clu_info = config["rdir"] + "/mmseqs_clustering/cluDB_info.tsv",
        index = config["rdir"] + "/mmseqs_clustering/cluDB_name_index.txt"
    shell:
        """
        set -e
        set -x

        # Create sequences database for the clustering
        {params.mmseqs_bin} createseqfiledb {params.seqdb} {params.cludb} {params.cluseqdb} --threads {threads} 2>{log.err}

        # To convert this cluster results tab separated file in wide format (repres member member member ..)
        awk -f {params.awk_wide} {input.clu} > {params.wide}

        # Cluster name rep and number of ORFs (size)
        awk -vFS='\\t' -vOFS='\\t' '{{print NR,$1,NF-1}}' {params.wide} > {params.info1}

        if [[ -s {params.cludb}.0 ]]; then
            {params.concat_clu} {params.cludb} {threads}
        fi

        if [[ -s {params.cluseqdb}.0 ]]; then
            {params.concat_clu} {params.cluseqdb} {threads}
        fi

        # Cluster naming
        # The official cluster names are going to be based on the line number of the wide formatted file
        # We are also going to produce a correspondence file to access the clusters in the MMseqs2 indices
        # Retrieving MMseqs-entry-name and sequence header from the clu_seqDB

        {params.mpi_runner} {params.mmseqs_bin} apply {params.cluseqdb} {params.namedb} -- {params.naming} --threads {threads}

        # Join the sequences with the representatives: this way we combine mmseqs-name with cluster-name (row number)
        join -12 -22 <(sed -e 's/\\x0//g' {params.namedb} | sort -k2,2 --parallel={threads} -T {params.local_tmp}) \
        <(sort -k2,2 --parallel={threads} -T {params.local_tmp} {params.info1} ) > {output.index}.tmp

        awk -vOFS='\\t' '{{print $2,$3}}' {output.index}.tmp > {output.index}

        # Rename the cluster sequence DB with the cluster names that are goiung to be used in the next steps
        {params.mmseqs_bin} renamedbkeys {output.index} {params.cluseqdb} {params.cluseqdb}.new

        mv {params.cluseqdb}.new {params.cluseqdb}
        mv {params.cluseqdb}.new.index {params.cluseqdb}.index
        mv {params.cluseqdb}.new.dbtype {params.cluseqdb}.dbtype

        # Join with the long format cluster file
        join -11 -22 <(sort --parallel={threads} -T {params.local_tmp} -k1,1 {input.clu} ) \
            <(sort --parallel={threads} -T {params.local_tmp} -k2,2 {params.info1} ) > {params.tmp}

        # Add length info
        sed -e 's/\\x0//g' {params.cluseqdb} | {params.seqtk_bin} comp | cut -f1,2 > {params.length}

        join -11 -22 <(sort -k1,1 --parallel={threads} -T {params.local_tmp} {params.length}) \
        <(sort -k2,2 --parallel={threads} -T {params.local_tmp} {params.tmp}) > {output.clu_info}

        # Reorder fields (cl_name rep orf length size)
        sort -k4,4n --parallel={threads} -T {params.local_tmp} {output.clu_info} | \
         awk -vOFS='\\t' '{{print $4,$3,$1,$2,$5}}' > {params.tmp}

        cp {params.tmp} {output.clu_info}

        # Subset the clusters based on their size

        # Clusters with more than 1 member
        awk -vOFS='\\t' '$3>=2{{print $1,$2}}' {params.info1} > {params.tmp1}
        join -12 -21 <(sort -k2,2 --parallel={threads} -T {params.local_tmp} {params.tmp1}) \
         <(sort -k1,1 --parallel={threads} -T {params.local_tmp} {input.clu} ) > {params.tmp}

        awk -vOFS='\\t' '{{print $2,$1,$3}}' {params.tmp} > {output.clusters}

        # Singletons
        awk -vOFS='\\t' '$3=="1"{{print $1,$2}}' {params.info1} > {output.singl}

        rm {params.tmp} {params.tmp1} {params.namedb} {params.namedb}.index {params.namedb}.dbtype {output.index}.tmp

        rm {params.wide} {params.length}

        """

rule mmseqs_clustering_results_done:
    input:
        clusters = config["rdir"] + \
            "/mmseqs_clustering/cluDB_no_singletons.tsv",
        singl = config["rdir"] + "/mmseqs_clustering/cluDB_singletons.tsv",
        clu_info = config["rdir"] + "/mmseqs_clustering/cluDB_info.tsv",
        index = config["rdir"] + "/mmseqs_clustering/cluDB_name_index.txt"
    output:
        cls_done = touch(config["rdir"] + "/mmseqs_clustering/clu_res.done")
    run:
        shell("echo 'CLUSTERING RESULTS DONE'")
