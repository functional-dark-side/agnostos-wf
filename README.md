## AGNOSTOS - workflow

<div class="img_container" style="width:60%; margin:2em auto;">

<img src="/assets/unkn_fig1.png" width="80%">

</div>

#### Disclaimer: This is a work in progress!

The workflow is still under development and it is subject to change. No guarantee is made regarding the functioning of the workflow and the accuracy of the results. Please, contact us in case you are interested in using it.

### Snakemake workflow usage:

The "agnostos-wf" snakemake workflow was developed using/runs in the de.NBI Cloud.
We used a cluster setup with 10 nodes of 28 cores and 252 G of memory each.
The cluster was build using BiBiGrid and it is using SLURM as Grid Batch Scheduler.

#### Check our [GitHub wiki](https://github.com/functional-dark-side/agnostos-wf/wiki) or the [usage file](AGNOSTOS_usage.md) to set AGNOSTOS up on your computer and to run the different modules!


### Run the workflow:

#### WARNING: UPDATES!
1. The workflow general structure slightly changed since the first release: the DB-creation and DB-update modules are now gathered into one workflow folder and to switch between them you just need to specify the module in the [config.yaml](workflow/config.yaml) file.
2. It is now possible to classify the singletons into the 4 categories (K, KWP, GU and EU) by setting `singl: "true"` in the [config.yaml](workflow/config.yaml) file.
3. It is also possible to re-validate and classify the existing GCs that got integrated with new genes by setting
`eval_shared: "true"` in the [config.yaml](workflow/config.yaml) file. At the moment, all shared GCs is the default (`shared_type: "all"`), other possibilities are: "discarded" for reprocessing all the previously discarded GCs with new sequences and "increase" for reprocessing all the shared GCs with 30% new sequences.
4. The databases can be downloaded and then removed in each step, to reduced the amount of space needed by the workflow (`db_mode: "memory"`)


**1. The DB-creation module** starts from a set of genomic/metagenomic contigs in fasta format and retrieves a database of categorised gene clusters and cluster communities.

```{bash}
cd workflow/
snakemake --use-conda -j 100 --config module="creation" --cluster-config config/cluster.yaml --cluster "sbatch --export=ALL -t {cluster.time} -c {threads} --ntasks-per-node {cluster.ntasks_per_node} --nodes {cluster.nodes} --cpus-per-task {cluster.cpus_per_task} --job-name {rulename}.{jobid} --partition {cluster.partition}" -R --until creation_workflow_report
```

**2. The DB-Update module** is used to add your genomic/metagenomic contigs or genes to an existing gene cluster database such as the agnostosDB dataset, which is stored in Figshare (https://doi.org/10.6084/m9.figshare.12459056) and publicly available for download. In case you cannot download the whole dataset, seen to the large size of many of the files, the workflow will download the necessary files for each step and it will then remove them. A description of the agnostosDB files can be found in the [AgnostosDB_README.md](AgnostosDB_README.md).

- To run the DB-update module of the workflow, you just need to enter the folder, modify the [config.yaml](workflow/config/config.yaml) and [config_communities.yml](workflow/config/config_communities.yml) files specifying your input data and the output paths (see [usage file](AGNOSTOS_usage.md)), and then run the same command shown above, this time modifying the configuration parameter 'module' to "update":

```{bash}
cd workflow/
snakemake -s Snakefile --use-conda -j 100 --config module="update" --cluster-config config/cluster.yaml --cluster "sbatch --export=ALL -t {cluster.time} -c {threads} --ntasks-per-node {cluster.ntasks_per_node} --nodes {cluster.nodes} --cpus-per-task {cluster.cpus_per_task} --job-name {rulename}.{jobid} --partition {cluster.partition}" -R --until update_workflow_report
```

NB:  If you want to run the DB-update module on the results of the DB-creation module, first copy or move the cluster database in the final DB-creation result folder as follows:

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
wget https://figshare.com/ndownloader/files/31247614 -O db_creation_data.tar.gz
tar -xzvf db_creation_data.tar.gz

wget https://ndownloader.figshare.com/files/25473335 -O db_update_data.tar.gz
tar -xzvf db_update_data.tar.gz
```

A brief description of the dataset is available on Figshare.

<br>

**Profile-search**: the profile-search vs the AgnostosDB cluster HMM profiles database is not part of the Snakemake workflow. However, if you want to search your set of genes against our profiles you just need to download the AGNOSTOS [gene cluster profiles](https://figshare.com/ndownloader/files/30998305) and the [gene cluster categories](https://ndownloader.figshare.com/files/23067140) and make sure you have [MMseqs2](https://github.com/soedinglab/MMseqs2) installed.
The scripts can be found in the [Profile_search/](Profile_search) folder.
To run the search you just need the following command:

```{bash}
# download the AGNOSTOS seed database gene cluster profiles
wget https://figshare.com/ndownloader/files/30998305 -O mmseqs_profiles.tar.gz
tar -xzvf mmseqs_profiles.tar.gz

# download the AGNOSTOS seed database gene cluster categories
wget https://ndownloader.figshare.com/files/23067140 -O cluster_ids_categ.tsv.gz
gunzip cluster_ids_categ.tsv.gz

Profile_search/profile_search.sh --query your-genes.fasta --clu_hmm mmseqs_profiles/clu_hmm_db --clu_cat cluster_ids_categ.tsv --mmseqs /path/to/mmseqs --mpi FALSE --threads 8
```

As additional option you can specify an additional file using "--info". This file should be a table with the correspondence of the genes to the contigs and genomes/MAGs or samples. The format should be gene - contig - genome (or sample_ID) - optional-info.

<br>

* * *

###### THE FUNCTIONAL DARK SIDE OF GENOMES AND METAGENOMES
To learn more about what we are doing check out our website  [dark.metagenomics.eu](https://dark.metagenomics.eu/).


* * *

#### Citation:

Vanni, C., Schechter, M., Acinas, S., Barberán, A., Buttigieg, P. L., Casamayor, E. O., Delmont, T. O., Duarte, C. M., Murat Eren, A., Finn, R., Kottmann, R., Mitchell, A., Sanchez, P., Siren, K., Steinegger, M., Glöckner, F. O., & Fernandez-Guerra, A. Unifying the known and unknown microbial coding sequence space. eLife 2022. [https://doi.org/10.7554/eLife.67667](https://doi.org/10.7554/eLife.67667)
