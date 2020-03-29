## The Agnostos workflow

<div class="img_container" style="width:60%; margin:2em auto;">

<img src="/assets/unkn_fig1.png" width="80%">

</div>

#### Disclaimer: This is a work in progress!

The workflow is still under development. Please, conctact us in case you are interested in using it.

#### Snakemake workflow usage:

The "agnostos-wf" snakemake workflow was developed using/runs in the de.NBI Cloud.
We used a cluster setup with 10 nodes of 28 cores and 252 G of memory each.
The cluster was build using BiBiGrid and it is using SLURM as Grid Batch Scheduler.

##### Before running it:

1.  Run the installation script [installation_script.sh](util/installation_script.sh).

2.  Check the config files in the [config/](config) folder. Change the programs (binaries) and output paths to your designated folders.

3.  Check that you have the required external DBs listed in the [config.yaml](config/config.yaml) file (in the "databases/" folder). In case you miss some of them, you can find the instructions for the dowload in [download_DBs.sh](util/download_DBs.sh)

When everything is set...

##### Run the workflow:

**1.  Database creation-mode**: Start from a set of genomic/metagenomic contigs in fasta format and retrieve a database of categorised gene clusters and cluster communities.

```{bash}
snakemake --use-conda -j 100 --cluster-config config/cluster.yaml --cluster "sbatch --export=ALL -t {cluster.time} -c {threads} --ntasks-per-node {cluster.ntasks_per_node} --nodes {cluster.nodes} --cpus-per-task {cluster.cpus_per_task} --job-name {rulename}.{jobid} --partition {cluster.partition}" -R --until workflow_report
```

**2.  Update-mode**: Add your genomic/metagenomic contigs or genes to the MG+GTDB cluster database, dowloadable from [here](). A description of the MG+GTDB database files can be found in the [clusterDB_README.md](clusterDB_README.md).

-   The cluster-update workflow is in the [cluster_update/](cluster_update) folder. To run it, you just need to enter the folder, modify the [config.yaml](cluster_update/config/config.yaml) and [config_communities.yaml](cluster_update/config/config_communities.yaml)files specifying your input data and the output paths, and then run the command:

```{bash}
snakemake -s Snakefile --use-conda -j 100 --cluster-config config/cluster.yaml --cluster "sbatch --export=ALL -t {cluster.time} -c {threads} --ntasks-per-node {cluster.ntasks_per_node} --nodes {cluster.nodes} --cpus-per-task {cluster.cpus_per_task} --job-name {rulename}.{jobid} --partition {cluster.partition}" -R --until workflow_report
```

NB: Make sure to remove eventual "*" from your fasta sequences. Below you have an example on how to do it:

```{bash}
seqkit replace -s -p "\*" -r "" my_seqs.fasta > my_seqs_clean.fasta
```

**3.  Profile-search**: the profile-search vs the MG+GTDB HMM profile DB is not part of the Snakemake workflow. However, if you want to search your set of genes against our profiles you just need to dowload the profile DB from [here](), or the full DB, including contextual data, from [here]() (in this second case you can find the profile DB in "cluster_category_DB/clu_hmm_db").
The scripts can be found in the [Profile_search/](Profile_search) folder.
To run the search you just need the following command:

```{bash}
Profile_search/profile_search.sh --query your-genes.fasta --info your-genes_add_info.tsv --mmseqs /path/to/mmseqs --threads num-threads
```

The "--info" file is optional, and should be a table with the correspondence of the genes to the contigs and/or genomes/MAGs/samples. The format should be gene - genome (or sample_ID) etc.


##### Citation:

Vanni, C., Schechter M., Delmont, T.O., Buttigieg P.L., Acinas S., Barberan A., Casamayor E.O., Siren K., Steinegger M., Eren M. A., Glöckner F.O and Fernandez-Guerra A.. (2020). Light into the darkness: Unifying the known and unknown coding sequence space in microbiome analyses. In preparation.
