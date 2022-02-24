# Installation and usage notes: setting up AGNOSTOS

The AGNOSTOS workflow **DB-creation** and **DB-update modules** were developed and tested in an HPC cluster setup with at least 4 nodes of 28 cores and 252 G of memory each, which uses [SLURM](https://slurm.schedmd.com/documentation.html) as Grid Batch Scheduler. To run one of these two modules follow the next stps:

1. Clone the repository: `git clone https://github.com/functional-dark-side/agnostos-wf` and `cd agnostos-wf/`

2. **Install packages not in Conda**. Please, check the installation script <installation_script.sh> (`sh installation_script.sh`) and in case you are missing some of the listed programs install them using the commands from the script.

  1. Note on MMseqs2: the program can be installed via conda, however the workflow was tested using the version "2f1db01c5109b07db23dc06df9d232e82b1b4b99" and newer releases could affect the workflow performance and results.

3. Check that you have the required external DBs listed in the [config.yaml](db_update/config/config.yaml) file (under "Databases"). In case you miss some of them, you can find the instructions for the download in the script <download_DBs.sh>. If you want to download all needed databases simply run `sh download_DBs.sh` (be patient this may take a while...).

4. **Check the configuration files (.yaml)** in the [config/](db_update/config) folders of the different modules. To change the programs (binaries) and output paths to your designated folders you can use the following commands:

(Example for the **DB-update module**)

```{bash}
# workflow directory
sed -i 's|vol/cloud/agnostos-wf/db_update|/your/wotkflow/path|g' db_update/config/config.yaml
sed -i 's|vol/cloud/agnostos-wf/db_update|/your/wotkflow/path|g' db_update/config/config_communities.yaml

# your data directory
sed -i 's|/vol/cloud/agnostos_test/db_update_data|/your/data/path|g' db_update/config/config.yaml

# your results directory
sed -i 's|/vol/cloud/agnostos_test/db_update|/your/results/path|g' db_update/config/config.yaml
sed -i 's|/vol/cloud/agnostos_test/db_update|/your/results/path|g' db_update/config/config_communities.yaml

# the directory of the existing GC database
sed -i 's|/vol/cloud/agnostosDB|/your/GC_DB/path|g' db_update/config/config.yaml

# the directory to the external databases
sed -i 's|/vol/cloud/agnostos-wf/databases|/your/external_database/path|g' db_update/config/config.yaml
sed -i 's|/vol/cloud/agnostos-wf/databases|/your/external_database/path|g' db_update/config/config_communities.yaml

# OPTIONAL: the directory to the binaries needed by the workflow,
# by default in the workflow folder under the directory bin/
sed -i 's|/vol/cloud/agnostos-wf/bin/|/your/binaries/path|g' db_update/config/config.yaml
```

Additionally you will have to specify if your data consists of contigs ("contigs"), self predicted genes sequences ("genes") or the gene prediction retrieved with [anvi'o](https://merenlab.org/software/anvio/help/7/programs/anvi-export-gene-calls/) ("anvio_genes"), and provide the name of the input files in the [config.yaml](db_update/config/config.yaml) file in the following entries:

```{yaml}
# Gene or contig file
new_data: "/your/data/path/your_genes.fasta"

# specify at which stage are your data, can be either "genes" or "contigs"
new_data_stage: "genes" #"contigs" or "genes" or "anvio_genes"

# If you already have the gene predictions, please provide path to gene completeness information
## In case your data comes from an anvi'o contigDB, please specify here the anvi'o gene_calls.tsv file,
## retrieved via "anvi-export-gene-calls -c CONTIGS.db -o anvio_gene_calls.tsv"
new_data_partial: "/vol/cloud/agnostos_test/db_update_data/new_genes_partial_info.tsv"
```

When everything is set you can run AGNOSTOS as follow:

(Example for the **DB-update module**)

```{bash}
cd db_update/
snakemake -s Snakefile --use-conda -j 100 --cluster-config config/cluster.yaml --cluster "sbatch --export=ALL -t {cluster.time} -c {threads} --ntasks-per-node {cluster.ntasks_per_node} --nodes {cluster.nodes} --cpus-per-task {cluster.cpus_per_task} --job-name {rulename}.{jobid} --partition {cluster.partition}" -R --until workflow_report
```

The **profile-search module** does not require an HPC environment and can be run on a local computer following the steps below and installing [MMseqs2](https://github.com/soedinglab/MMseqs2) in case you don't have it yet:

```{bash}
# download the AGNOSTOS seed database gene cluster profiles
wget https://ndownloader.figshare.com/files/23066963 -O mmseqs_profiles.tar.gz
tar -xzvf mmseqs_profiles.tar.gz

# download the AGNOSTOS seed database gene cluster categories
wget https://ndownloader.figshare.com/files/23067140 -O cluster_ids_categ.tsv.gz
gunzip cluster_ids_categ.tsv.gz

# Run the sequence-profile search
Profile_search/profile_search.sh --query your-genes.fasta --clu_hmm mmseqs_profiles/clu_hmm_db --clu_cat cluster_ids_categ.tsv --mmseqs /path/to/mmseqs --mpi FALSE --threads 8
```

NOTE: On MAC-OS-X you will probably need to install the gnu-getopt, which supports long options (--). For this you can use the command `conda install -c bioconda gnu-getopt` or with Homebrew as `brew install gnu-getopt`.
