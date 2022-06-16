rule spurious_shadow:
    input:
        fasta = config["rdir"] + "/gene_prediction/orf_seqs.fasta",
        clusters = config["rdir"] + "/mmseqs_clustering/cluDB_no_singletons.tsv"
    threads: 16
    priority: 30
    conda:
        config["conda_env"]
    params:
        hmmer_bin = config["hmmer_bin"],
        mpi_runner = config["mpi_runner"],
        hmmpress = config["hmmpress_bin"],
        antifamdb = config["antifam_db"],
        local_tmp = config["mmseqs_local_tmp"],
        module = config["module"],
        eval_shared = config["eval_shared"],
        db_mode = config["db_mode"],
        shadowr = "scripts/shadow_orfs.r",
        clu_info = config["rdir"] + "/mmseqs_clustering/cluDB_info.tsv",
        hmmout = config["rdir"] + "/spurious_shadow/hmmsearch_antifam_sp.out",
        hmmlog = config["rdir"] + "/spurious_shadow/hmmsearch_antifam_sp.log",
        partial = config["rdir"] + "/gene_prediction/orf_partial_info.tsv",
        spur = config["rdir"] + "/spurious_shadow/spurious_orfs.tsv",
        all_shad = config["rdir"] + "/spurious_shadow/all_shadow_orfs.tsv",
        only_shad = config["rdir"] + "/spurious_shadow/shadow_orfs_info.tsv",
        or_sp_sh = config["ordir"] + "/spurious_shadow_info.tsv.gz",
        tmp1 = config["rdir"] + "/spurious_shadow/tmpl1",
        tmp2 = config["rdir"] + "/spurious_shadow/tmpl2"
    output:
        sp_sh = config["rdir"] + "/spurious_shadow/spurious_shadow_info.tsv"
    log:
        out = "logs/spsh_stdout.log",
        err = "logs/spsh_stderr.err"
    benchmark:
        "benchmarks/spurious_shadow.tsv"
    shell:
        """
        set -x
        set -e

        export OMPI_MCA_btl=^openib
        export OMP_NUM_THREADS={threads}
        export OMP_PROC_BIND=FALSE

        # 1. Detection of spurious ORFs

        if [ ! -s {params.antifamdb} ]; then
            echo "Dowloading AntiFam database"
            DB=$(dirname {params.antifamdb})
            wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/AntiFam/current/Antifam.tar.gz -P ${{DB}}
            tar xvfz ${{DB}}/Antifam.tar.gz --directory ${{DB}}
            rm ${{DB}}/Antifam.tar.gz
            {params.hmmpress} {params.antifamdb}
        fi

        NANTI=$(grep -c '^NAME' {params.antifamdb})
        NSEQS=$(grep -c '^>' {input.fasta})
        N=$(($NSEQS * $NANTI))

        # Run hmmsearch (MPI mode)
        {params.mpi_runner} {params.hmmer_bin} --mpi --cut_ga -Z "${{N}}" --domtblout {params.hmmout} -o {params.hmmlog} {params.antifamdb} {input.fasta} 2>{log.err} 1>{log.out}

        if [[ {params.db_mode} == "memory" ]]; then
            DB=$(dirname {params.antifamdb})
            rm -rf ${{DB}}/AntiFam*
        fi

        # Parse the results
        grep -v '^#' {params.hmmout} > {params.spur}.tmp || true > {params.spur}.tmp 2>>{log.err}

        rm {params.hmmout} {params.hmmlog}

        if [[ -s {params.spur}.tmp ]]; then
            awk '$13<=1e-05 && !seen[$1]++{{print $1}}' {params.spur}.tmp > {params.spur} 2>>{log.err}
            rm {params.spur}.tmp
        else
            mv {params.spur}.tmp {params.spur}
        fi

        # 2. Detection of shadow ORFs
        ./{params.shadowr} --orfs {params.partial} \
                           --shadows {params.all_shad} \
                           --threads {threads} 2>>{log.err}

        ## 2.1 Parsing of results adding cluster information
        ## Add cluster info on the first column of ORFs
        join -12 -23 <(sort --parallel={threads} -k2,2 {params.all_shad}) \
                <(sort --parallel={threads} -T {params.local_tmp} -k3,3 {params.clu_info} ) > {params.tmp1}

        ## reorder columns...
        awk -vOFS='\\t' '{{print $6,$7,$8,$9,$1,$11,$12,$13,$3,$4,$5,$2,$10}}' {params.tmp1} > {params.tmp2} && mv {params.tmp2} {params.tmp1}

        ## Add cluster info on the second column of ORFs
        join -11 -23 <(sort --parallel={threads} -k1,1 {params.tmp1}) \
                <(sort --parallel={threads} -T {params.local_tmp} -k3,3 {params.clu_info} ) > {params.tmp2}

        ## reorder columns...
        awk -vOFS='\\t' '{{print $1,$14,$15,$16,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13}}' {params.tmp2} > {params.tmp1}

        ### 2.2 Filter to get only the real shadow ORFs
        ### print ORFs and "TRUE" or "BOTH" if they fall in the definition of shadows
        awk -vFS='\\t' -vOFS='\\t' '{{if($4>$11) print $8,"TRUE"; else if($11>$4) print $1,"TRUE"; else if($4==$11 && $3<$10) print1; else if($4==$11 && $3>$10) print $1,"TRUE";}}' {params.tmp1} > {params.tmp2}
        awk -vFS='\\t' -vOFS='\\t' '{{if($4==$11 && $3==$10) print $1,"BOTH";}}' {params.tmp1} >> {params.tmp2}
        awk -vFS='\\t' -vOFS='\\t' '{{if($4==$11 && $3==$10) print $8,"BOTH";}}' {params.tmp1} >> {params.tmp2}

        ### re-join with the cluster information table
        join -13 -21 -a1 -a2 <(sort --parallel={threads} -T {params.local_tmp} -k3,3 {params.clu_info} ) \
          <(sort --parallel={threads} -k1,1 <(awk '!seen[$0]++' {params.tmp2}) ) > {params.only_shad}

        ### Add a "FALSE" to the non-shadow ORFs
        awk -vOFS='\\t' '!seen[$0]++{{if($6=="TRUE" || $6=="BOTH") print $2,$1,$3,$4,$5,$6; else print $2,$1,$3,$4,$5,"FALSE";}}' {params.only_shad} > {params.tmp1}
        mv {params.tmp1} {params.only_shad}

        # 3. Combine spurious and shadow information
        join -12 -21 -a1 <(sort --parallel={threads} -k2,2 {params.only_shad}) \
            <(awk '{{print $1,"TRUE"}}' {params.spur} | sort --parallel={threads} -k1,1) > {params.tmp1}

        awk '{{if($7=="TRUE") print $1,$4,$2,$3,$5,$6,$7; else print $1,$4,$2,$3,$5,$6,"FALSE";}}' {params.tmp1} | \
            awk '{{print $1,$2,$3,$5,$6,$7}}' > {output.sp_sh}.tmp

        # calculate proportion of shadows per cluster
        # grep 'TRUE\|BOTH' {params.only_shad} | \
        awk '($6 == "TRUE" || $6 == ""BOTH) && (!seen[$0]++){{a[$1"\\t"$5]+=1}}END{{for(i in a) print i,a[i]}}' | awk '{{print $1,$2,$3/$2}}' > {params.tmp1}

        # join with the whole table on cluster ID
        join -13 -21 -a1 <(sort --parallel={threads} -T {params.local_tmp} -k3,3 {output.sp_sh}.tmp) \
            <(sort --parallel={threads} -k1,1 {params.tmp1}) > {params.tmp2}

        # Add 0 for the clusters with no shadows and
        # reorder column output table: orf - length - cl_name - cl_size - prop_shadow - is.shadow - is.spurious
        awk -vOFS='\\t' '{{ for(i=1; i<=7; i++) if($i ~ /^ *$/) $i = 0 }};1' {params.tmp2} | \
            awk -vOFS='\\t' '{{print $2,$3,$1,$4,$7,$5,$6}}' > {output.sp_sh}.tmp

        if [[ {params.module} == "update" ]]; then
            # Download original dataset spurious and shadow info and combine with the new dataset
            if [[ ! -s {params.or_sp_sh} ]]; then
                wget https://ndownloader.figshare.com/files/23067131 -O {params.or_sp_sh}
            fi
            if [[ {params.eval_shared} == "true" ]]; then
                join -11 -21 -v1 <(zcat {params.or_sp_sh} |\
                 sort -k1,1 --parallel {threads} -T {params.local_tmp}) \
                <(awk '{{print $1}}' {output.sp_sh}.tmp |\
                 sort -k1,1 --parallel {threads} -T {params.local_tmp} ) > {params.tmp1}
                 sed -i 's/ /\\t/g' {params.tmp1}
                 cat {output.sp_sh}.tmp {params.tmp1} > {output.sp_sh}
            else
                cat {output.sp_sh}.tmp <(zcat {params.or_sp_sh}) > {output.sp_sh}
            fi
            rm {output.sp_sh}.tmp
        else
            mv {output.sp_sh}.tmp {output.sp_sh}
        fi

        rm {params.tmp1} {params.tmp2} {params.all_shad} {params.spur} {params.only_shad}
        """

rule spurious_shadow_done:
    input:
        sp_sh = config["rdir"] + "/spurious_shadow/spurious_shadow_info.tsv"
    output:
        sp_sh_done = touch(config["rdir"] + "/spurious_shadow/spsh.done")
    run:
        shell("echo 'Spurious and Shadows detection DONE'")
