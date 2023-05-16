rule cluster_category_refinement:
    input:
        k = config["rdir"] + "/cluster_classification/k_ids.txt",
        kwp = config["rdir"] + "/cluster_classification/kwp_ids.txt",
        gu = config["rdir"] + "/cluster_classification/gu_ids.txt",
        eu = config["rdir"] + "/cluster_classification/eu_ids.txt"
    threads: config["threads_cat_ref"]
    conda:
        config["conda_env"]
    params:
        mmseqs_bin = config["mmseqs_bin"],
        mpi_runner = config["mpi_runner"],
        hhsuite = config["hhsuite"],
        famsa_bin = config["famsa_bin"],
        vmtouch = config["vmtouch"],
        hhblits_bin_mpi = config["hhblits_bin_mpi"],
        hhblits_prob = config["hhblits_prob"],
        uniclust_db = config["uniclust_db"],
        pfam_db = config["pfam_hh_db"],
        pfam_clan = config["pfam_clan"],
        hypo_filt = config["hypo_filt"],
        module = config["module"],
        db_mode = config["db_mode"],
        step_unk = "unknown refinement",
        step_k = "known refinement",
        categ_db = "scripts/category_db_files.sh",
        reformat = "scripts/reformat_file.sh",
        consensus = "scripts/consensus.sh",
        hhmake = "scripts/hhmake.sh",
        hhblits_search = "scripts/hhblits_search_parser.sh",
        hhpfam_parser = "scripts/hhpfam_parser.sh",
        concat = "scripts/concat_multi_annot.awk",
        new_da = "scripts/categ_ref_new_domain_archit.r",
        patterns = "scripts/hypothetical_grep.tsv",
        idir = config["rdir"] + "/cluster_classification",
        cluseqdb = config["rdir"] + "/cluster_refinement/refined_clusterDB",
        ref_clu = config["rdir"] + "/cluster_refinement/refined_clusters.tsv",
        categ_orfs = config["rdir"] + \
            "/cluster_categories/cluster_ids_categ_genes.tsv",
        outdir = config["rdir"] + "/cluster_categories",
        hmm_eu = config["rdir"] + "/cluster_categories/eu_hhm_db",
        hmm_kwp = config["rdir"] + "/cluster_categories/kwp_hhm_db",
        tmp_eu = config["rdir"] + "/cluster_categories/eu_hhbl",
        tmp_kwp = config["rdir"] + "/cluster_categories/kwp_hhbl"

    output:
        k = config["rdir"] + "/cluster_categories/k_ids.txt",
        kwp = config["rdir"] + "/cluster_categories/kwp_ids.txt",
        gu = config["rdir"] + "/cluster_categories/gu_ids.txt",
        eu = config["rdir"] + "/cluster_categories/eu_ids.txt",
        categ = config["rdir"] + "/cluster_categories/cluster_ids_categ.tsv"
    log:
        out = "logs/cat_ref_stdout.log",
        err = "logs/cat_ref_stderr.err"
    benchmark:
        "benchmarks/cat_ref.tsv"
    shell:
        """
        set -e
        set -x

        export OMPI_MCA_btl=^openib
        export OMP_NUM_THREADS={threads}
        export OMP_PROC_BIND=FALSE

        ### Environmental unknown refinement searching for remote homologies (in the Uniclust DB)
        if [[ ! -f {params.hmm_eu} ]]; then
            {params.categ_db} --cluseq_db {params.cluseqdb} \
                        --step "{params.step_unk}" \
                        --mpi_runner "{params.mpi_runner}" \
                        --mmseqs {params.mmseqs_bin} \
                        --aln {params.famsa_bin} \
                        --hhsuite {params.hhsuite} \
                        --reformat {params.reformat} \
                        --consensus {params.consensus} \
                        --hhmake {params.hhmake} \
                        --outdir {params.outdir} \
                        --idir {params.idir} \
                        --clu_hhm "none" \
                        --threads {threads} 2>{log.err} 1>{log.out}
        fi

        # Pfam HHBLITS DB
        if [ ! -s {params.pfam_db}_hhm.ffdata ]; then
            echo "Dowloading Pfam hh-suite database"
            DB=$(dirname {params.pfam_db})
            wget http://wwwuser.gwdg.de/~compbiol/data/hhsuite/databases/hhsuite_dbs/pfamA_34.0.tar.gz -O ${{DB}}/pfam.tar.gz
            tar xvfz ${{DB}}/pfam.tar.gz --directory ${{DB}}/
            rm ${{DB}}/pfam.tar.gz
        fi

        # Uniclust HHBLITS DB
        if [ ! -s {params.uniclust_db}_a3m.ffdata ]; then
            DB=$(dirname {params.uniclust_db})
            echo "Dowloading Uniclust hh-suite database"
            aria2c --file-allocation=none -d ${{DB}} -c -x 10 -s 10 http://wwwuser.gwdg.de/~compbiol/uniclust/2021_03/UniRef30_2021_03.tar.gz
            tar xvfz {params.uniclust_db}.tar.gz --directory ${{DB}}/
            rm {params.uniclust_db}.tar.gz
        fi

        # vmtouch the DBs
        {params.vmtouch} -f {params.uniclust_db}*
        {params.vmtouch} -f {params.pfam_db}*
        {params.vmtouch} -f {params.hmm_eu}

        # If EU GCs <= 10,000 search them vs Uniclust and parse results
        EU=$(wc -l {input.eu} | cut -d' ' -f1)
        if [[ ! -f {params.tmp_eu}_parsed.tsv ]]; then
            if [[ $EU -le 10000 ]]; then
                # Run directly EU vs Uniclust
                if [[ ! -f {params.tmp_eu}.index ]]; then
                    {params.mpi_runner} {params.hhblits_bin_mpi} -i {params.hmm_eu} \
                                                                -o {params.tmp_eu} \
                                                                -n 2 -cpu 1 -v 0 \
                                                                -d {params.uniclust_db} 2>>{log.err} 1>>{log.out}
                    mv {params.tmp_eu}.ffdata {params.tmp_eu}
                    mv {params.tmp_eu}.ffindex {params.tmp_eu}.index
                fi
                {params.mpi_runner} {params.mmseqs_bin} apply {params.tmp_eu} {params.tmp_eu}.parsed \
                    -- {params.hhblits_search} {params.outdir}/templ {threads} {params.hhblits_prob} --threads 1 2>>{log.err} 1>>{log.out}

                 sed -e 's/\\x0//g' {params.tmp_eu}.parsed > {params.tmp_eu}_parsed.tsv
                 rm -rf {params.outdir}/templ
            else
                # Run EU vs PFAM and then EU vs Uniclust
                ## EU vs Pfam
                if [[ ! -f {params.tmp_eu}.index ]]; then
                    {params.mpi_runner} {params.hhblits_bin_mpi} -i {params.hmm_eu} \
                                                                -o {params.tmp_eu} \
                                                                -n 2 -cpu 1 -v 0 \
                                                                -d {params.pfam_db} 2>>{log.err} 1>>{log.out}
                    mv {params.tmp_eu}.ffdata {params.tmp_eu}
                    mv {params.tmp_eu}.ffindex {params.tmp_eu}.index
                fi
                {params.mpi_runner} {params.mmseqs_bin} apply {params.tmp_eu} {params.tmp_eu}.parsed_pfam \
                    -- {params.hhblits_search} {params.outdir}/templ {threads} {params.hhblits_prob} --threads 1 2>>{log.err} 1>>{log.out}
                rm -rf {params.outdir}/templ
                sed -e 's/\\x0//g' {params.tmp_eu}.parsed_pfam > {params.tmp_eu}.tsv

                #Parsing hhblits result files and filtering for hits with coverage > 0.4, and removing overlapping pfam domains/matches
                ./{params.hhpfam_parser} {params.tmp_eu}.tsv {threads} > {params.tmp_eu}_filt.tsv

                if [ -s {params.tmp_eu}_filt.tsv ]; then
                    # join with the pfam names and clans
                    join -11 -21 <(awk '{{split($1,a,"."); print a[1],$3,$8,$9}}' {params.tmp_eu}_filt.tsv | sort -k1,1) \
                        <(gzip -dc {params.pfam_clan} |\
                        awk -vFS='\\t' -vOFS='\\t' '{{print $1,$2,$4}}' |\
                        awk -vFS='\\t' -vOFS='\\t' '{{for(i=1; i<=NF; i++) if($i ~ /^ *$/) $i = "no_clan"}};1' \
                        | sort -k1,1) > {params.tmp_eu}_name_acc_clan.tsv

                    # Multi domain format
                    awk '{{print $2,$3,$4,$5,$1,$6}}' {params.tmp_eu}_name_acc_clan.tsv |\
                        sort -k1,1 -k2,3g | \
                        awk -vOFS='\\t' '{{print $1,$4,$5,$6}}' | \
                        awk -f {params.concat} > {params.tmp_eu}_name_acc_clan_multi.tsv

                        rm {params.tmp_eu}_name_acc_clan.tsv

                    if [ -s {params.tmp_eu}_name_acc_clan_multi.tsv ]; then
                        # Divide the new hits with pfam into DUFs and not DUFs
                        ./{params.new_da} {params.tmp_eu}_name_acc_clan_multi.tsv {params.tmp_eu}

                        # New K clusters
                        awk '{{print $1}}' {params.tmp_eu}_new_k_ids_annot.tsv >> {output.k}
                        cat {input.k} >> {output.k}
                        # New GU clusters
                        awk '{{print $1}}' {params.tmp_eu}_new_gu_ids_annot.tsv >> {output.gu}
                        # Remaning EU clusters
                        join -11 -21 -v1 <(sort -k1,1 {input.eu}) \
                            <(awk '{{print $1}}' {params.tmp_eu}_name_acc_clan_multi.tsv | sort -k1,1) > {output.eu}.1
                        # Subset the EU hmm DB
                        {params.mmseqs_bin} createsubdb {output.eu}.1 {params.hmm_eu} {params.hmm_eu}.left
                        mv {params.hmm_eu}.left {params.hmm_eu}
                        mv {params.hmm_eu}.left.index {params.hmm_eu}.index
                    fi
                fi
                # Run remaining EU GCs vs Uniclust
                {params.vmtouch} -f {params.hmm_eu}
                if [[ ! -f {params.tmp_eu}_unicl.index ]]; then
                    {params.mpi_runner} {params.hhblits_bin_mpi} -i {params.hmm_eu} \
                                                                 -o {params.tmp_eu}_unicl \
                                                                 -n 2 -cpu 1 -v 0 \
                                                                 -d {params.uniclust_db} 2>>{log.err} 1>>{log.out}
                   mv {params.tmp_eu}_unicl.ffdata {params.tmp_eu}_unicl
                   mv {params.tmp_eu}_unicl.ffindex {params.tmp_eu}_unicl.index
                fi
                {params.mpi_runner} {params.mmseqs_bin} apply {params.tmp_eu}_unicl {params.tmp_eu}.parsed_unicl \
                    -- {params.hhblits_search} {params.outdir}/templ {threads} {params.hhblits_prob} --threads 1 2>>{log.err} 1>>{log.out}

                sed -e 's/\\x0//g' {params.tmp_eu}.parsed_unicl > {params.tmp_eu}_parsed.tsv
                rm -rf {params.outdir}/templ
            fi
        fi

        # Parse EU refinement results
        LC_ALL=C rg -j 4 -i -f {params.patterns} {params.tmp_eu}_parsed.tsv | \
            awk '{{print $0"\thypo"}}' > {params.tmp_eu}_hypo_char
        LC_ALL=C rg -j 4 -v -i -f {params.patterns} {params.tmp_eu}_parsed.tsv | \
            awk '{{print $0"\tchar"}}' >> {params.tmp_eu}_hypo_char

        sed -i 's/\\t\\t/\\t/g' {params.tmp_eu}_hypo_char

        if [ -s {params.tmp_eu}_hypo_char ]; then
            awk -vP={params.hypo_filt} -vFS='\\t' \
                '{{a[$2][$16]++}}END{{for (i in a) {{N=a[i]["hypo"]/(a[i]["hypo"]+a[i]["char"]); if (N >= P){{print i}}}}}}' {params.tmp_eu}_hypo_char \
                > {params.tmp_eu}_new_gu_ids.txt

            join -11 -21 -v1 <(awk '!seen[$2]++{{print $2}}' {params.tmp_eu}_parsed.tsv | sort -k1,1) \
                 <(sort -k1,1 {params.tmp_eu}_new_gu_ids.txt) > {params.tmp_eu}_new_kwp_ids.txt

            if [[ -f {output.eu}.1 ]]; then
                join -11 -21 -v1 <(sort -k1,1 {output.eu}.1) \
                    <(awk '!seen[$2]++{{print $2}}' {params.tmp_eu}_parsed.tsv | sort -k1,1) > {output.eu}
                rm {output.eu}.1
            else
                join -11 -21 -v1 <(sort -k1,1 {input.eu}) \
                    <(awk '!seen[$2]++{{print $2}}' {params.tmp_eu}_parsed.tsv | sort -k1,1) > {output.eu}
            fi

            cat {input.kwp} {params.tmp_eu}_new_kwp_ids.txt >> {output.kwp}

            cat {input.gu} {params.tmp_eu}_new_gu_ids.txt >> {output.gu}
            if [[ ! -s {output.k} ]]; then
                cp {input.k} {output.k}
            fi
        else
            echo "No EU was found to be a remote homolog of an already observed protein"
            cp {input.k} {output.k}
            cp {input.kwp} {output.kwp}
            cp {input.gu} {output.gu}
            cp {input.eu} {output.eu}
        fi

        # clean the EU results ...
        rm -rf {params.hmm_eu}*
        rm -rf {params.tmp_eu} {params.tmp_eu}.index {params.tmp_eu}.dbtype
        rm -rf {params.tmp_eu}_unicl {params.tmp_eu}_unicl.index {params.tmp_eu}_unicl.dbtype
        rm -rf {params.tmp_eu}.parsed_* {params.tmp_eu}.tsv {params.tmp_eu}_filt.tsv
        rm -rf {params.tmp_eu}_parsed.tsv {params.tmp_eu}_hypo_char

        #######
        ## Known without Pfam refinement, searching for remote homologies (in the Pfam DB)

        # Pfam clans
        if [ ! -s {params.pfam_clan} ]; then
          echo "Dowloading Pfam-A clan information"
          wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam34.0/Pfam-A.clans.tsv.gz -O {params.pfam_clan}
        fi

        if [[ ! -f {params.hmm_kwp} ]]; then
            {params.categ_db} --cluseq_db {params.cluseqdb} \
                          --step "{params.step_k}" \
                          --mpi_runner "{params.mpi_runner}" \
                          --mmseqs {params.mmseqs_bin} \
                          --aln {params.famsa_bin} \
                          --hhsuite {params.hhsuite} \
                          --reformat {params.reformat} \
                          --consensus {params.consensus} \
                          --hhmake {params.hhmake} \
                          --outdir {params.outdir} \
                          --idir {params.idir} \
                          --clu_hhm "none" \
                          --threads {threads} 2>>{log.err} 1>>{log.out}
        fi

        rm -f {params.hmm_kwp}.ff*
        ln -s {params.hmm_kwp} {params.hmm_kwp}.ffdata
        ln -s {params.hmm_kwp}.index {params.hmm_kwp}.ffindex

        if [[ ! -f {params.tmp_kwp}.index || ! -f {params.tmp_kwp}.tsv ]]; then
            {params.vmtouch} -f {params.hmm_kwp}
            if [[ ! -f {params.tmp_kwp}.index ]]; then
                {params.mpi_runner} {params.hhblits_bin_mpi} -i {params.hmm_kwp} \
                                                            -o {params.tmp_kwp} \
                                                            -n 2 -cpu 1 -v 0 \
                                                            -d {params.pfam_db}
               mv {params.tmp_kwp}.ffdata {params.tmp_kwp}
               mv {params.tmp_kwp}.ffindex {params.tmp_kwp}.index
            fi
            {params.mpi_runner} {params.mmseqs_bin} apply {params.tmp_kwp} {params.tmp_kwp}.parsed \
                    -- {params.hhblits_search} {params.outdir}/templ {threads} {params.hhblits_prob} --threads 1
            rm -rf  {params.outdir}/templ
            # Parsing hhr result files and filtering for hits with probability ≥ 90%
            sed -e 's/\\x0//g' {params.tmp_kwp}.parsed > {params.tmp_kwp}.tsv
        fi

        #Parsing hhblits result files and filtering for hits with coverage > 0.4, and removing overlapping pfam domains/matches
        ./{params.hhpfam_parser} {params.tmp_kwp}.tsv {threads} > {params.tmp_kwp}_filt.tsv

        # join with the pfam names and clans
        join -11 -21 <(awk '{{split($1,a,"."); print a[1],$3,$8,$9}}' {params.tmp_kwp}_filt.tsv | sort -k1,1) \
            <(gzip -dc {params.pfam_clan} |\
            awk -vFS='\\t' -vOFS='\\t' '{{print $1,$2,$4}}' |\
            awk -vFS='\\t' -vOFS='\\t' '{{for(i=1; i<=NF; i++) if($i ~ /^ *$/) $i = "no_clan"}};1' | sort -k1,1) > {params.tmp_kwp}_name_acc_clan.tsv

        # Multi domain format
        awk '{{print $2,$3,$4,$5,$1,$6}}' {params.tmp_kwp}_name_acc_clan.tsv |\
            sort -k1,1 -k2,3g | \
            awk -vOFS='\\t' '{{print $1,$4,$5,$6}}' | \
            awk -f {params.concat} > {params.tmp_kwp}_name_acc_clan_multi.tsv

        rm {params.tmp_kwp}_name_acc_clan.tsv

        if [ -s {params.tmp_kwp}_name_acc_clan_multi.tsv ]; then
            # Divide the new hits with pfam into DUFs and not DUFs
            ./{params.new_da} {params.tmp_kwp}_name_acc_clan_multi.tsv {params.tmp_kwp}

            # New K clusters
            awk '{{print $1}}' {params.tmp_kwp}_new_k_ids_annot.tsv >> {output.k}
            # New GU clusters
            awk '{{print $1}}' {params.tmp_kwp}_new_gu_ids_annot.tsv >> {output.gu}
            # New KWP clusters
            join -11 -21 -v1 <(sort -k1,1 {output.kwp}) \
                    <(awk '{{print $1}}' {params.tmp_kwp}_name_acc_clan_multi.tsv | sort -k1,1) > {params.tmp_kwp}_kwp
            mv {params.tmp_kwp}_kwp {output.kwp}
            # EU remain the same
        else
            # The categories mantain the same clusters
            echo "No KWP was found to be remote homolog of a Pfam domain"
        fi

        # clean the results ...
        rm -rf {params.hmm_kwp} {params.hmm_kwp}.index {params.hmm_kwp}.dbtype {params.hmm_kwp}.ff*
        rm -rf {params.tmp_kwp} {params.tmp_kwp}.index {params.tmp_kwp}.dbtype
        rm -rf {params.tmp_kwp}.parsed* {params.tmp_kwp}_filt.tsv {params.tmp_kwp}.tsv

        #Clean eventual duplicates
        awk '!seen[$0]++' {output.k} > {output.k}.tmp && mv {output.k}.tmp {output.k}
        awk '!seen[$0]++' {output.kwp} > {output.kwp}.tmp && mv {output.kwp}.tmp {output.kwp}
        awk '!seen[$0]++' {output.gu} > {output.gu}.tmp && mv {output.gu}.tmp {output.gu}
        awk '!seen[$0]++' {output.eu} > {output.eu}.tmp && mv {output.eu}.tmp {output.eu}

        # Create a file with cl_name and category
        cat <(awk -vOFS='\\t' '{{print $1,"K"}}' {output.k}) \
            <(awk -vOFS='\\t' '{{print $1,"KWP"}}' {output.kwp}) \
            <(awk -vOFS='\\t' '{{print $1,"GU"}}' {output.gu}) \
            <(awk -vOFS='\\t' '{{print $1,"EU"}}' {output.eu}) > {output.categ}

        # Add ORFs
        join -11 -21 <(sort -k1,1 {output.categ}) \
            <(awk '!seen[$1,$3]++{{print $1,$3}}' {params.ref_clu} | sort -k1,1) \
            > {params.categ_orfs}

        gzip {params.categ_orfs}

        # Gather cluster annotations obtained from the classification and the two refinement steps
        class=$(dirname {input.k})
        categ=$(dirname {output.k})

        # GU annotations:
        join -11 -21 <(sort -k1,1 {output.gu}) \
                <(awk '{{print $1,"Uniclust",$3,$4}}' {params.tmp_eu}_parsed.tsv | sort -k1,1) \
                > {params.outdir}/GU_annotations.tsv
        if [ -s {params.tmp_eu}_new_gu_ids_annot.tsv ]; then
            awk '{{print $1,"Pfam","0.0",$2}}' {params.tmp_eu}_new_gu_ids_annot.tsv >> {params.outdir}/GU_annotations.tsv
        fi
        awk '{{print $1,"Pfam","0.0",$2}}' {params.tmp_kwp}_new_gu_ids_annot.tsv >> {params.outdir}/GU_annotations.tsv
        cat {params.idir}/gu_annotations.tsv >> {params.outdir}/GU_annotations.tsv
        sed -i 's/ /\\t/g' {params.outdir}/GU_annotations.tsv

        # KWP annotations
        join -11 -21 <(sort -k1,1 {output.kwp}) \
            <(awk '{{print $1,"Uniclust",$3,$4}}' {params.tmp_eu}_parsed.tsv | sort -k1,1) \
            > {params.outdir}/KWP_annotations.tsv
        join -11 -21 <(sort -k1,1 {output.kwp}) \
                <(sort -k1,1 {params.idir}/kwp_annotations.tsv) >> {params.outdir}/KWP_annotations.tsv
        sed -i 's/ /\\t/g' {params.outdir}/KWP_annotations.tsv

        # K annotations
        if [ -s {params.tmp_eu}_new_k_ids_annot.tsv ]; then
            awk '{{print $1,"Pfam","0.0",$2}}' {params.tmp_eu}_new_k_ids_annot.tsv >> {params.outdir}/K_annotations.tsv
        fi
        cat {params.idir}/k_annotations.tsv \
                <(awk -vOFS='\\t' '{{print $1,"Pfam","0.0",$2}}' {params.tmp_kwp}_new_k_ids_annot.tsv) \
                    >> {params.outdir}/K_annotations.tsv
        sed -i 's/ /\\t/g' {params.outdir}/K_annotations.tsv


        if [[ {params.db_mode} == "memory" ]]; then
            rm {params.uniclust_db}* {params.pfam_db}*
        fi

      """

rule cluster_categ_ref_done:
    input:
        k = config["rdir"] + "/cluster_categories/k_ids.txt",
        kwp = config["rdir"] + "/cluster_categories/kwp_ids.txt",
        gu = config["rdir"] + "/cluster_categories/gu_ids.txt",
        eu = config["rdir"] + "/cluster_categories/eu_ids.txt",
        categ = config["rdir"] + "/cluster_categories/cluster_ids_categ.tsv"
    output:
        cat_ref_done = touch(
            config["rdir"] + "/cluster_categories/cat_ref.done")
    run:
        shell("echo 'CATEGORY REFINEMENT DONE'")
