## The Agnostos workflow

<div class="img_container" style="width:60%; margin:2em auto;">

<img src="/assets/unkn_fig1.png" width="80%">

</div>

#### Disclaimer: This is a work in progress!

The workflow is still under development and it is subject to change. No guarentee is made regarding the functioning of the workflow and the accuracy of the results. Please, conctact us in case you are interested in using it.

### Snakemake workflow usage:

The "agnostos-wf" snakemake workflow was developed using/runs in the de.NBI Cloud.
We used a cluster setup with 10 nodes of 28 cores and 252 G of memory each.
The cluster was build using BiBiGrid and it is using SLURM as Grid Batch Scheduler.

#### Before running AGNOSTOS check the [usage file](AGNOSTOS_usage) to set it up on your computer!

### Run the workflow:

**1.  DB-creation module**: Start from a set of genomic/metagenomic contigs in fasta format and retrieve a database of categorised gene clusters and cluster communities.

```{bash}
cd db_creation/
snakemake --use-conda -j 100 --cluster-config config/cluster.yaml --cluster "sbatch --export=ALL -t {cluster.time} -c {threads} --ntasks-per-node {cluster.ntasks_per_node} --nodes {cluster.nodes} --cpus-per-task {cluster.cpus_per_task} --job-name {rulename}.{jobid} --partition {cluster.partition}" -R --until workflow_report
```

**2.  DB-Update module**: Add your genomic/metagenomic contigs or genes to the agnostosDB dataset is stored in Figshare (https://doi.org/10.6084/m9.figshare.12459056) and publicy available for download. In case you cannot download the whole dataset, seen to the large size of many of the files, the workflow will download the necessary files for each step and it will then remove them. A description of the agnostosDB files can be found in the [AgnostosDB_README.md](AgnostosDB_README.md).

-   The DB-update workflow is in the [db_update/](db_update) folder. To run it, you just need to enter the folder, modify the [config.yaml](db_update/config/config.yaml) and [config_communities.yml](db_update/config/config_communities.yml) files specifying your input data and the output paths (see [usage file](AGNOSTOS_usage)), and then run the command:

```{bash}
cd db_update/
snakemake -s Snakefile --use-conda -j 100 --cluster-config config/cluster.yaml --cluster "sbatch --export=ALL -t {cluster.time} -c {threads} --ntasks-per-node {cluster.ntasks_per_node} --nodes {cluster.nodes} --cpus-per-task {cluster.cpus_per_task} --job-name {rulename}.{jobid} --partition {cluster.partition}" -R --until workflow_report
```

NB:  If you want to run the DB-update module on the results of the DB-creation module, first copy or move the cluster database in the final DB-update results:

```{bash}
mv db_creation/mmseqs_clustering db_creation/clusterDB_results/
```

<br>

**Output**

The output of these 2 modules is described in the [Output_README.md](Output_README.md).

<br>

### Test dataset
Set of 10K contigs to test the DB-creation module and 5K contigs to test the DB-update module of the workflow.
The test-dataset can be downloaded from Figshare as follows:

```{bash}
mkdir -p agnostos_test
cd agnostos_test
wget https://ndownloader.figshare.com/files/25473332 -O db_creation_data.tar.gz
tar -xzvf db_creation_data.tar.gz

wget https://ndownloader.figshare.com/files/25473335 -O db_update_data.tar.gz
tar -xzvf db_update_data.tar.gz
```

A brief description of the dataset is available on Figshare.

<br>

**Profile-search**: the profile-search vs the AgnostosDB cluster HMM profiles database is not part of the Snakemake workflow. However, if you want to search your set of genes against our profiles you just need to dowload the AGNOSTOS [gene cluster profiles](https://ndownloader.figshare.com/files/23066963) and the [gene cluster categories](https://ndownloader.figshare.com/files/23067140) and make sure you have [MMseqs2](https://github.com/soedinglab/MMseqs2) installed.
The scripts can be found in the [Profile_search/](Profile_search) folder.
To run the search you just need the following command:

```{bash}
# download the AGNOSTOS seed database gene cluster profiles
wget https://ndownloader.figshare.com/files/23066963 -O mmseqs_profiles.tar.gz
tar -xzvf mmseqs_profiles.tar.gz

# download the AGNOSTOS seed database gene cluster categories
wget https://ndownloader.figshare.com/files/23067140 -O cluster_ids_categ.tsv.gz
gunzip cluster_ids_categ.tsv.gz

Profile_search/profile_search.sh --query your-genes.fasta --clu_hmm mmseqs_profiles/clu_hmm_db --clu_cat cluster_ids_categ.tsv --mmseqs /path/to/mmseqs --threads 8
```

As additional option you can specify an additional file using "--info". This file should be a table with the correspondence of the genes to the contigs and genomes/MAGs or samples. The format should be gene - contig - genome (or sample_ID) etc.

<br>

* * *

###### THE FUNCTIONAL DARK SIDE OF GENOMES AND METAGENOMES
To lern more about what we are doing check out our website  [dark.metagenomics.eu](https://dark.metagenomics.eu/).


* * *

#### Citation:

Vanni, C., Schechter, M., Acinas, S., Barberán, A., Buttigieg, P. L., Casamayor, E. O., Delmont, T. O., Duarte, C. M., Murat Eren, A., Finn, R., Kottmann, R., Mitchell, A., Sanchez, P., Siren, K., Steinegger, M., Glöckner, F. O., & Fernandez-Guerra, A. (2020). Light into the darkness: Unifying the known and unknown coding sequence space in microbiome analyses. In bioRxiv (p. 2020.06.30.180448). [https://doi.org/10.1101/2020.06.30.180448](https://doi.org/10.1101/2020.06.30.180448)
