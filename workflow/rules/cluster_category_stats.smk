rule cluster_category_stats:
    input:
        k_db = config["rdir"] + "/cluster_category_DB/k_cons.index",
        cl_cat = config["rdir"] + "/cluster_categories/cluster_ids_categ.tsv",
    threads: 28
    params:
        workdir = config["wdir"],
        mmseqs_bin = config["mmseqs_bin"],
        kaiju_bin = config["kaiju_bin"],
        eggnog_bin = config["eggnog_bin"],
        mpi_runner = config["mpi_runner"],
        tmpl = config["mmseqs_local_tmp"],
        vmtouch = config["vmtouch"],
        conf_tax_db = config["conf_tax_db"],
        gtdb = config["gtdb_tax"],
        nr = config["nr_tax"],
        DPD = config["DPD"],
        dpd_info = config["dpd_info"],
        aa_mutants = config["aa_mutants"],
        mutantDB = config["mutantDB"],
        db_mode = config["db_mode"],
        kaiju_tax = "scripts/kaiju_taxonomy.sh",
        kaiju_parse = "scripts/kaiju_add_taxonomy.R",
        stats = "scripts/cluster_category_stats.r",
        ref = config["rdir"] + "/cluster_refinement/refined_clusters.tsv",
        refdb = config["rdir"] + "/cluster_refinement/refined_clusterDB",
        cl_cat_genes = config["rdir"] + "/cluster_categories/cluster_ids_categ_genes.tsv.gz",
        outdir = config["rdir"] + "/cluster_category_stats",
        kaiju_res = config["rdir"] + "/cluster_category_stats/cluster_kaiju_taxonomy.tsv",
        mutant_dir = config["rdir"] + "/cluster_category_stats/mutant_phenotypes",
        mutants = config["rdir"] + "/cluster_category_stats/mutant_phenotypes/cluster_mutants.tsv",
        dark = config["rdir"] + "/cluster_category_stats/darkness/cluster_category_darkness.tsv",
        dark_dir = config["rdir"] + "/cluster_category_stats/darkness",
        eggnog_dir = config["rdir"] + "/cluster_category_stats/eggnog",
        compl = config["rdir"] + "/cluster_category_stats/cluster_category_completeness.tsv"
    conda:
        config["conda_env"]
    output:
        HQ_clusters = config["rdir"] + "/cluster_category_stats/HQ_clusters.tsv",
        cat_stats = config["rdir"] + "/cluster_category_stats/cluster_category_summary_stats.tsv"
    log:
        out="logs/stats_stdout.log",
        err="logs/stats_stderr.err"
    benchmark:
        "benchmarks/cat_stats.tsv"
    shell:
        """

        if [ {params.conf_tax_db} == "gtdb" ]; then
            TAXDB={params.gtdb}
            if [ ! -s {params.gtdb} ]; then
                echo "Dowloading GTDB-r89 kaiju DB"
                DB=$(dirname {params.mutantDB})
                wget https://ndownloader.figshare.com/files/24745184 -O ${{DB}}/gtdb.tar.gz
                tar xzvf ${{DB}}/gtdb.tar.gz --directory ${{DB}}
                rm ${{DB}}/gtdb.tar.gz
            fi
        else
            TAXDB={params.nr}

            if [[ ! -s {params.nr} ]]; then
                DIR=$(basedir ${{TAXDB}})
                cd ${{DIR}}
                {params.kaiju_bin}-makedb -s nr_euk
                cd {params.workdir}
            fi

        fi

        ## Cluster Kaiju taxonomy with GTDB r89
        if [[ ! -s {params.kaiju_res} ]]; then
            {params.vmtouch} -f ${{TAXDB}}
            ./{params.kaiju_tax} --search {params.kaiju_bin} \
                              --input {params.refdb} \
                              --taxdb ${{TAXDB}} \
                              --parsing {params.kaiju_parse} \
                              --output {params.kaiju_res} \
                              --tmpl {params.tmpl} \
                              --threads {threads} 2>{log.err} 1>{log.out}
        fi

        if [[ {params.db_mode} == "memory" ]]; then
            rm ${{TAXDB}}
        fi

        # Extract all sequences from the refined database set:
        sed -e 's/\\x0//g' {params.refdb} > {params.outdir}/refined_cl_genes.fasta

        # Create MMseqs2 databases for the refined genes
        {params.mmseqs_bin} createdb {params.outdir}/refined_cl_genes.fasta {params.outdir}/refined_cl_genes_db

        ## Cluster level of darkness
        mkdir -p {params.dark_dir}

        if [ ! -s {params.DPD} ]; then
            echo "Dowloading DPD"
            wget https://ndownloader.figshare.com/files/23756312 -O {params.DPD}
            wget https://ndownloader.figshare.com/files/23756306 -O {params.dpd_info}
        fi

        if [[ ! -s {params.dark} ]]; then
          {params.mmseqs_bin} createdb {params.DPD} {params.dark_dir}/dpd_db
          # Search
          {params.mmseqs_bin} search {params.outdir}/refined_cl_genes_db {params.dark_dir}/dpd_db \
            {params.dark_dir}/refined_cl_genes_dpd_db {params.dark_dir}/tmp \
            --threads {threads} --max-seqs 300 \
            -e 1e-20 --cov-mode 0 -c 0.6 --mpi-runner "{params.mpi_runner}"

          {params.mmseqs_bin} convertalis {params.outdir}/refined_cl_genes_db {params.dark_dir}/dpd_db \
            {params.dark_dir}/refined_cl_genes_dpd_db {params.dark_dir}/refined_cl_genes_dpd.tsv \
            --format-mode 2 --threads {threads} \
            --format-output 'query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,qcov,tcov'

          rm -rf {params.dark_dir}/dpd_db* {params.dark_dir}/refined_cl_genes_dpd_db* {params.dark_dir}/tmp

          # Extract best-hits
          export LANG=C; export LC_ALL=C; \
          sort -k1,1 -k11,11g -k13,13gr -k14,14gr --parallel {threads} -T {params.tmpl} \
          {params.dark_dir}/refined_cl_genes_dpd.tsv | \
            sort -u -k1,1 --parallel {threads} -T {params.tmpl} --merge > {params.dark_dir}/refined_cl_genes_dpd_bh.tsv

          # Join with cluster categories
          join -11 -23 <(awk '{{print $1,$2}}' {params.dark_dir}/refined_cl_genes_dpd_bh.tsv | \
          sort -k1,1 --parallel {threads} -T {params.tmpl}) \
            <(sort -k3,3 --parallel {threads} -T {params.tmpl} <(zcat {params.cl_cat_genes})) > {params.dark}

          sed -i 's/ /\\t/g' {params.dark}
        fi

        ## Cluster vs Price et al. mutants

        mkdir -p {params.mutant_dir}

        # Mutant phenotypes (Price et al. 2018)
        if [ ! -s {params.aa_mutants} ]; then
            ## Amino acid sequences
            wget https://fit.genomics.lbl.gov/cgi_data/aaseqs -O {params.aa_mutants}
        fi
        if [ ! -s {params.mutantDB} ]; then
            ## Contextual data
            wget https://fit.genomics.lbl.gov/cgi_data/feba.db -O {params.mutantDB}
        fi

        if [[ ! -s {params.mutants} ]]; then
            {params.mmseqs_bin} createdb {params.aa_mutants} {params.mutant_dir}/mutant_db
            # Search
            {params.mmseqs_bin} search {params.outdir}/refined_cl_genes_db {params.mutant_dir}/mutant_db \
                {params.mutant_dir}/refined_cl_genes_mutant_db {params.mutant_dir}/tmp \
                --threads {threads} --max-seqs 300 \
                -e 1e-20 --cov-mode 0 -c 0.6 --mpi-runner "{params.mpi_runner}"

            {params.mmseqs_bin} convertalis {params.outdir}/refined_cl_genes_db {params.mutant_dir}/mutant_db \
                {params.mutant_dir}/refined_cl_genes_mutant_db {params.mutant_dir}/refined_cl_genes_mutants.tsv \
                --format-mode 2 --threads {threads} \
                --format-output 'query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,qcov,tcov'

            rm -rf {params.outdir}/refined_cl_genes_db* {params.mutant_dir}/tmp
            rm -rf {params.mutant_dir}/mutant_db* {params.mutant_dir}/refined_cl_genes_mutant_db*

            # Extract best-hits
            export LANG=C; export LC_ALL=C;\
             sort -k1,1 -k11,11g -k13,13gr -k14,14gr --parallel {threads} -T {params.tmpl}\
             {params.mutant_dir}/refined_cl_genes_mutants.tsv | \
                sort -u -k1,1 --merge --parallel {threads} -T {params.tmpl} > {params.mutant_dir}/refined_cl_genes_mutants_bh.tsv

            # Join with cluster categories
            join -11 -23 <(awk '{{print $1,$2}}' {params.mutant_dir}/refined_cl_genes_mutants_bh.tsv |\
             sort -k1,1 --parallel {threads} -T {params.tmpl}) \
                <(sort -k3,3 --parallel {threads} -T {params.tmpl}\
                 <(zcat {params.cl_cat_genes})) > {params.mutants}

            sed -i 's/ /\\t/g' {params.mutants}
        fi

        ## EggNOG annotations
        mkdir -p {params.eggnog_dir}
        if [[ ! -s {params.eggnog_dir}/cluster_eggnogs.tsv ]]; then

            pip install biopython

            {params.eggnog_bin} -m diamond --no_annot --no_file_comment --cpu {threads} \
                                -i {params.outdir}/refined_cl_genes.fasta \
                                --output {params.eggnog_dir}/refined_nogs --override

            NOG=$(dirname {params.eggnog_bin})
            scp ${{NOG}}/data/* /dev/shm/
            {params.eggnog_bin} --annotate_hits_table {params.eggnog_dir}/refined_nogs.emapper.seed_orthologs \
                                --no_file_comments -o {params.eggnog_dir}/refined \
                                --cpu {threads} --override --data_dir /dev/shm
            rm /dev/shm/*

            awk -vFS='\\t' -vOFS='\\t' '{{print $1,$7,$8}}' {params.eggnog_dir}/refined.emapper.annotations |\
                sed 's/ /_/g' > {params.eggnog_dir}/refined_eggnogs.tsv

            # Join with cluster categories
            join -11 -23 <(awk '{{print $1,$2,$3}}' {params.eggnog_dir}/refined_eggnogs.tsv |\
             sort -k1,1 --parallel {threads} -T {params.tmpl}) \
                <(sort -k3,3 --parallel {threads} -T {params.tmpl}\
                 <(zcat {params.cl_cat_genes})) > {params.eggnog_dir}/cluster_eggnogs.tsv

            sed -i 's/ /\\t/g' {params.eggnog_dir}/cluster_eggnogs.tsv
            rm  {params.eggnog_dir}/refined_nogs*  {params.eggnog_dir}/refined_eggnogs.tsv
        fi

        rm {params.outdir}/refined_cl_genes.fasta

        ## Cluster general stats

        ./{params.stats} --ref_clu {params.ref} \
                       --clu_categ {input.cl_cat} \
                       --kaiju_tax {params.kaiju_res} \
                       --clu_dark {params.dark} \
                       --dpd_info {params.dpd_info} \
                       --compl {params.compl} \
                       --hq_clu {output.HQ_clusters} \
                       --mutantDB {params.mutantDB} \
                       --mutants {params.mutants} \
                       --eggnog {params.eggnog_dir}/cluster_eggnogs.tsv \
                       --summ_stats {output.cat_stats} \
                       --output {params.outdir} 2>>{log.err} 1>>{log.out}

        if [[ {params.db_mode} == "memory" ]]; then
            rm {params.dpd_info} {params.DPD}
            rm {params.aa_mutants} {params.mutantDB}
        fi

        """

rule cluster_categ_stats_done:
  input:
    HQ_clusters = config["rdir"] + "/cluster_category_stats/HQ_clusters.tsv",
    cat_stats = config["rdir"] + "/cluster_category_stats/cluster_category_summary_stats.tsv"
  output:
    cat_stats_done = touch(config["rdir"] + "/cluster_category_stats/cat_stats.done")
  run:
      shell("echo 'COMMUNITIES INFERENCE DONE'")
