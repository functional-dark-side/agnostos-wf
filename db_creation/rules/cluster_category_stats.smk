rule cluster_category_stats:
    input:
        k_db = config["rdir"] + "/cluster_category_DB/k_cons.index",
        cl_cat = config["rdir"] + "/cluster_categories/cluster_ids_categ.tsv",
    threads: 28
    params:
        mmseqs_bin = config["mmseqs_bin"],
        kaiju_bin = config["kaiju_bin"],
        mpi_runner = config["mpi_runner"],
        tmpl = config["mmseqs_local_tmp"],
        vmtouch = config["vmtouch"],
        taxdb = config["taxdb"],
        gtdb = config["gtdb_tax"],
        DPD = config["DPD"],
        dpd_info = config["dpd_info"],
        mmseqs_tax = "scripts/mmseqs_taxonomy.sh",
        kaiju_tax = "scripts/kaiju_taxonomy.sh",
        kaiju_parse = "scripts/kaiju_add_taxonomy.R",
        stats = "scripts/cluster_category_stats.r",
        ref = config["rdir"] + "/cluster_refinement/refined_clusters.tsv",
        refdb = config["rdir"] + "/cluster_refinement/refined_clusterDB",
        cl_cat_genes = config["rdir"] + "/cluster_categories/cluster_ids_categ_genes.tsv.gz",
        tax_dir = config["rdir"] + "/cluster_category_stats/taxonomy",
        tax = config["rdir"] + "/cluster_category_stats/taxonomy/cluster_mmseqs_taxonomy.tsv",
        kaiju_res = config["rdir"] + "/cluster_category_stats/taxonomy/cluster_kaiju_taxonomy.tsv",
        dark = config["rdir"] + "/cluster_category_stats/darkness/cluster_category_darkness.tsv",
        dark_dir = config["rdir"] + "/cluster_category_stats/darkness",
        compl = config["rdir"] + "/cluster_category_stats/cluster_category_completeness.tsv",
        outdir = config["rdir"] + "/cluster_category_stats"
    conda:
        config["conda_env"]
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
        ## Cluster mmseqs2 taxonomy with UniProtKB
        mkdir -p {params.tax_dir}
        if [[ ! -s {params.tax} ]]; then
            {params.vmtouch} -f {params.taxdb}
            ./{params.mmseqs_tax} --search {params.mmseqs_bin} \
                              --input {params.refdb} \
                              --taxdb {params.taxdb} \
                              --cl_info {params.cl_cat_genes}  \
                              --output {params.tax} \
                              --outdir {params.outdir} \
                              --mpi_runner "{params.mpi_runner}" \
                              --threads {threads} 2>{log.err} 1>{log.out}
        fi
        ## Cluster Kaiju taxonomy with GTDB r89
        if [[ ! -s {params.kaiju_res} ]]; then
            {params.vmtouch} -f {params.gtdb}
            ./{params.kaiju_tax} --search {params.kaiju_bin} \
                              --input {params.refdb} \
                              --taxdb {params.gtdb} \
                              --parsing {params.kaiju_parse} \
                              --output {params.kaiju_res} \
                              --tmpl {params.tmpl} \
                              --threads {threads} 2>>{log.err} 1>>{log.out}
        fi
        ## Cluster level of darkness
        mkdir -p {params.dark_dir}
        if [[ ! -s {params.dark} ]]; then
          # Extract all sequences from the refined database set:
          sed -e 's/\\x0//g' {params.refdb} | gzip > {params.dark_dir}/refined_cl_genes.fasta.gz

          # Create MMseqs2 databases
          {params.mmseqs_bin} createdb {params.dark_dir}/refined_cl_genes.fasta.gz {params.dark_dir}/refined_cl_genes_db
          {params.mmseqs_bin} createdb {params.DPD} {params.dark_dir}/dpd_db
          # Search
          {params.mmseqs_bin} search {params.dark_dir}/refined_cl_genes_db {params.dark_dir}/dpd_db \
            {params.dark_dir}/refined_cl_genes_dpd_db {params.dark_dir}/tmp \
            --threads {threads} --max-seqs 300 \
            -e 1e-20 --cov-mode 0 -c 0.6 --mpi-runner "{params.mpi_runner}"

          {params.mmseqs_bin} convertalis {params.dark_dir}/refined_cl_genes_db {params.dark_dir}/dpd_db \
            {params.dark_dir}/refined_cl_genes_dpd_db {params.dark_dir}/refined_cl_genes_dpd.tsv \
            --format-mode 2 --threads {threads} \
            --format-output 'query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,qcov,tcov'

          rm -rf {params.dark_dir}/refined_cl_orfs.fasta.gz {params.dark_dir}/refined_cl_genes_db* {params.dark_dir}/dpd_db* {params.dark_dir}/refined_cl_genes_dpd_db* {params.dark_dir}/tmp

          # Extract best-hits
          export LANG=C; export LC_ALL=C; sort -k1,1 -k11,11g -k13,13gr -k14,14gr {params.dark_dir}/refined_cl_genes_dpd.tsv | \
            sort -u -k1,1 --merge > {params.dark_dir}/refined_cl_genes_dpd_bh.tsv

          # Join with cluster categories
          join -11 -23 <(awk '{{print $1,$2}}' {params.dark_dir}/refined_cl_genes_dpd_bh.tsv | sort -k1,1) \
            <(sort -k3,3 <(zcat {params.cl_cat_genes})) > {params.dark}

          sed -i 's/ /\\t/g' {params.dark}
        fi
        ## Cluster general stats

        ./{params.stats} --ref_clu {params.ref} \
                       --clu_categ {input.cl_cat} \
                       --mmseqs_tax {params.tax} \
                       --kaiju_tax {params.kaiju_res} \
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
