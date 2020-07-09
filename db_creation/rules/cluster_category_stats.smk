rule cluster_category_stats:
    input:
        eu_db = config["rdir"] + "/cluster_category_DB/eu_cons.index",
        cl_cat = config["rdir"] + "/cluster_categories/cluster_ids_categ.tsv",
    threads: 28
    params:
        mmseqs_bin = config["mmseqs_bin"],
        mpi_runner = config["mpi_runner"],
        taxdb = config["taxdb"],
        DPD = config["DPD"],
        dpd_info = config["dpd_info"],
        mmseqs_tax = "scripts/mmseqs_taxonomy.sh",
        stats = "scripts/cluster_category_stats.r",
        ref = config["rdir"] + "/cluster_refinement/refined_clusters.tsv",
        refdb = config["rdir"] + "/cluster_refinement/refined_clusterDB",
        cl_cat_orfs = config["rdir"] + "/cluster_categories/cluster_ids_categ_orfs.tsv",
        tax_dir = config["rdir"] + "/cluster_category_stats/taxonomy",
        tax = config["rdir"] + "/cluster_category_stats/taxonomy/cluster_category_taxonomy.tsv",
        dark = config["rdir"] + "/cluster_category_stats/darkness/cluster_category_darkness.tsv",
        dark_dir = config["rdir"] + "/cluster_category_stats/darkness",
        compl = config["rdir"] + "/cluster_category_stats/cluster_category_completeness.tsv",
        outdir = config["rdir"] + "/cluster_category_stats"
    conda:
        "../envs/workflow.yml"
    output:
        HQ_clusters = config["rdir"] + "/cluster_category_stats/HQ_clusters.tsv",
        cat_stats = config["rdir"] + "/cluster_category_stats/cluster_category_summary_stats.tsv"
    log:
        out="logs/stats_stdout.log",
        err="logs/stats_stderr.err"
    benchmark:
        "benchmarks/cluster_category_stats/cat_stats.tsv"
    shell:
        """
        ## Cluster mmseqs2 taxonomy
        mkdir -p {params.tax_dir}
        ./{params.mmseqs_tax} --search {params.mmseqs_bin} \
                              --input {params.refdb} \
                              --taxdb {params.taxdb} \
                              --cl_info {params.cl_cat_orfs}  \
                              --output {params.tax} \
                              --outdir {params.outdir} \
                              --mpi_runner "{params.mpi_runner}" \
                              --threads {threads} 2>{log.err} 1>{log.out}

        ## Cluster level of darkness
        mkdir -p {params.dark_dir}
        # Extract all sequences from the refined database set:
        sed -e 's/\\x0//g' {params.refdb} | gzip > {params.dark_dir}/refined_cl_orfs.fasta.gz

        # Create MMseqs2 databases
        {params.mmseqs_bin} createdb {params.dark_dir}/refined_cl_orfs.fasta.gz {params.dark_dir}/refined_cl_orfs_db
        {params.mmseqs_bin} createdb {params.DPD} {params.dark_dir}/dpd_db
        # Search
        {params.mmseqs_bin} search {params.dark_dir}/refined_cl_orfs_db {params.dark_dir}/dpd_db \
          {params.dark_dir}/refined_cl_orfs_dpd_db {params.dark_dir}/tmp \
          --threads {threads} --max-seqs 300 \
          -e 1e-20 --cov-mode 0 -c 0.6 --mpi-runner "{params.mpi_runner}"

        {params.mmseqs_bin} convertalis {params.dark_dir}/refined_cl_orfs_db {params.dark_dir}/dpd_db \
          {params.dark_dir}/refined_cl_orfs_dpd_db {params.dark_dir}/refined_cl_orfs_dpd.tsv \
          --format-mode 2 --threads {threads} \
          --format-output 'query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,qcov,tcov'

        rm {params.dark_dir}/refined_cl_orfs_db* {params.dark_dir}/dpd_db* {params.dark_dir}/refined_cl_orfs_dpd_db* {params.dark_dir}/tmp

        # Extract best-hits
        export LANG=C; export LC_ALL=C; sort -k1,1 -k11,11g -k13,13gr -k14,14gr {params.dark_dir}/refined_cl_orfs_dpd.tsv | \
          sort -u -k1,1 --merge > {params.dark_dir}/refined_cl_orfs_dpd_bh.tsv

        # Join with cluster categories
        join -11 -23 <(awk '{{print $1,$2}}' {params.dark_dir}/refined_cl_orfs_dpd_bh.tsv | sort -k1,1) \
          <(sort -k3,3 {params.cl_cat_orfs}) > {params.dark}

        sed -i 's/ /\\t/g' {params.dark}

        ## Cluster general stats


        ./{params.stats} --ref_clu {params.ref} \
                       --clu_categ {input.cl_cat} \
                       --clu_tax {params.tax} \
                       --clu_dark {params.dark} \
                       --dpd_info {params.dpd_info} \
                       --compl {params.compl} \
                       --hq_clu {output.HQ_clusters} \
                       --summ_stats {output.cat_stats} \
                       --output {params.outdir} 2>>{log.err} 1>>{log.out}

        """

rule cluster_categ_stats_done:
  input:
    HQ_clusters = config["rdir"] + "/cluster_category_stats/HQ_clusters.tsv",
    cat_stats = config["rdir"] + "/cluster_category_stats/cluster_category_summary_stats.tsv"
  output:
    cat_stats_done = touch(config["rdir"] + "/cluster_category_stats/cat_stats.done")
  run:
      shell("echo 'COMMUNITIES INFERENCE DONE'")
