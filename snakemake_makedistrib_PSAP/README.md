# Snakemake workflow to make PSAP null distributions 

## Introduction
This first workflow **snakemake_makedistrib_PSAP** allows the easy calculation of PSAP null distributions according to allele frequencies from a vcf file and a pathogenicity score (here the CADD score).
The input data is highly customizable and the use of Snakemake makes it possible to run the pipeline multiple times with different input parameters without running again the steps creating files that have already been created.

## Usage 

Once the conda environment with Snakemake installed is setup and activated, the user can modify the configuration file according to the desired input files and parameters, and run the pipeline.

### 1. Configure the pipeline

#### config.yaml file

**snakemake_makedistrib_PSAP/config/config.yaml** contains all pipeline parameters which can be tuned by the user:

|**Parameter** | **Description** |
|--------------|--------------|
|**snakemake_directory** |Path to `/snakemake_makedistrib_PSAP` directory on the user's machine|
|**genes** |bed file for the coordinates of genes' coding regions. The file for GRCh37 is provided with the pipeline and a new one for GRCh38 can be created using [BioMart](https://www.ensembl.org/info/data/biomart/index.html)|
|**cadd_regions** |bed file for the coordinates of CADD regions. The file for GRCh37 is provided with the pipeline or can be downloaded [here](https://lysine.univ-brest.fr/RAVA-FIRST/)|
|**coding_cadd_regions** |bed file for the coordinates of coding CADD regions. Created from the intersection of genes and CADD regions bed files|
|**cadd_file**|Score file used to calibrated PSAP null distributions. Format of a [CADD file](https://cadd.gs.washington.edu/download) currently supported. File needs to be indexed|
|**score_prefix** |Prefix describing the score file used|
|**allele_frequencies** |Allele frequencies file used to calibrated PSAP null distributions, in the format vcf.gz|
|**coverage** |bed file of regions well-covered regions in allele frequency database|
|**af_prefix** |Prefix describing the allele frequencies file used|
|**outdir** |Output directory|
|**outfile** |Prefix for output files, will be the name of the PSAP null distributions file|
|**units** |Unit of testing to construct PSAP null distributions, can take the values "gene", "cadd_region" or "coding_cadd_region"|

For the coverage file, the command line used to get the bed file provided for gnomAD genome V2 (good coverage = 90% of individuals at dp10) was:
```
zcat gnomad.genomes.r2.0.1.coverage.txt.gz | tail -n+2 | awk '{print $1"\t"($2-1)"\t"$2"\t"$7}' > gnomad.genome.r2.0.1.dp10.bed
```

The workflow can also be run for all three units at the same time by replacing `units = config["units"]` in the Snakefile by: `units = ["gene","coding_caddregion","caddregion"]`

> TIP: different analyses can be run using just one cloned repository. This is achieved by changing the outdir and outfile in the configuration file. Also different parameters values can be used in the different analyses.

### 2. Create the Conda environments

The Snakemake relies on conda environements which contain all the necessary dependencies to run the pipeline. The conda environements need to be installed first, which can take some time.
This step is easily done using the following command line:

```
snakemake --use-conda --conda-create-envs-only --conda-frontend conda -j 1
```

### 3. Run the pipeline.

Once the pipeline is configured and conda environments are created, the user just needs to run Snakemake pipeline to make PSAP null distributions:

```
snakemake --use-conda -j 22 
```

The mandatory arguments are:
* **--use-conda**: to use the conda environemnts created at the previous step
* **-j**: number of threads/jobs provided to snakemake. The pipeline splits each step by chromosome so running it on at least 22 threads is recommended. 

An optional argument **--dry-run** can be to used to make a dry run of the pipeline, check if there are no warnings and see all the files that will be created.
With the additionnal argument **--configfile=/path_to_configfile/name_config.yaml**, the user can use a configuration file in a different location instead of modifying the default config.yaml file in the snakemake directory.

> TIP: If the execution of a step of the pipeline fails or if you want to see how it progresses, you can check the `{outdir}/log` directory where a log file is created for each rule of the Snakemake by chromosome if applicable.
> A `{outdir}/benchmark` directory is also created automatically to monitor how much memory and cpu is used by each step of the pipeline.

## Pipeline steps

The steps of the pipeline are not meant to be used separately and will be run by Snakemake in the most efficient way depending on input and output files.

### 1. split_multiallelic_panel

The vcf file with allele frequencies from the reference database is split by chromosome and multi-allelic variants are split. Output located in `{outdir}/reference_panel`.

### 2. make_allele_frequencies_file

The allele frequency files created at the previous step are processed by default to keep only PASS SNVs, in well-covered regions of the database for the calculation of PSAP null distributions. 
If other filtering needs to be applied, or all filtering removed (no column `FILTER` in the input vcf), the following line can be modified in the Snakefile:
```
bcftools view -i 'FILTER=="PASS"' --regions-file {params.coverage} --types snps {input.file_allele_frequencies_split}
```
The vcf input file needs to have the columns `CHROM`, `POS`, `REF`, `ALT`, `CHROM`, `AF`, `AC` which will be kept in the output files. If any column has a different name or if another column needs to be used for allele frequencies,
the following line in the Snakefile can be modified:
```
bcftools query -f'%CHROM\t%POS\t%REF\t%ALT\t%AF\t%AC\n'
```
Output located in `{outdir}/reference_panel`.

### 3. make_notpassvariants_file

Variants that do not pass the filtering descrbed in the previous step will have to be removed from the calculation of PSAP null distributions downstream of the pipeline. Changes can be made to the default filtering as needed.
In this step, a file by chromosome is created with the coordinates of these variants. Output located in `{outdir}/reference_panel`.

### 4. list_units_wellcovered_panel

The bed files for genes, CADD regions and coding CADD regions are intersected with the bed file of good coverage in the allele frequencies database using [bedmap](https://bedops.readthedocs.io/en/latest/content/reference/statistics/bedmap.html).
A list of genes, CADD regions and coding CADD regions well-covered in at least 50% of their length in the database is computed. Final PSAP null distributions will only include this list of units. Output located in `{outdir}/reference_panel`.

### 5. make_input_table_forPSAP

In this step, information used to calculate PSAP null distributions is gathered in tables by chromosome using [bedtools](https://bedtools.readthedocs.io/en/latest/index.html). These tables include the chromosome, start, end, reference allele,
alternative allele, raw CADD score, PHRED CADD score, corresponding CADD region and gene if there is one of all possible SNVs of the genome in the well-covered regions by database. Output located in `{outdir}/score_tables`.

### 6. calculate_null_distributions_PSAP

This step carries out the main function of the pipeline, which is to calculate PSAP null distributions. Overall, it evaluates the probability of seeing a heterozygote or homozygote variant with such a maximal CADD score 
as high or higher in the chosen unit of testing (gene, CADD regions or coding CADD region) in the reference database. Variants from the step `make_notpassvariants_file` are removed from the analysis. The PHRED CADD score is used (column 7).
One file by chromosome and by model is created, each line corresponds to a unit of testing and columns to PSAP probabilities. Output located in `{outdir}/temp_lookup_bychr`.


### 7. make_final_null_distributions_PSAP

The null distributions by chromosome are gathered in a unique table for the autosomal dominant (het) model and autosomal recessive (hom model), for the units of testing well-covered in at least 50% of their length in the database. 
The final two output files created are `{outfile}_het.txt.gz` and `{outfile}_hom.txt.gz`. Output located in `{outdir}/final_lookuptables_PSAP`.

## Configuration of computation resources

This Snakemake pipeline is portable to cluster engines, see the [Snakemake cluster documentation](https://snakemake.readthedocs.io/en/stable/executing/cluster.html) for more detailed information. Taking the SLURM scheduler as an example, 
the pipeline can be run without any changes using the following command line:

```
snakemake --use-conda --slurm --default-resources slurm_account=<your SLURM account> slurm_partition=<your SLURM partition>
```

For a custom use of ressources and input parameters, a **profile** can be used. When the argument `--profile` is added, snakemake looks for a directory with the name of the given profile (here **slurm**) 
containing a `config.yaml` file. This file contains the parameter `cluster` which tells snakemake how to submit jobs to the cluster. In the case of slurm, the `sbatch` command is used with its arguments.
Other arguments include `jobs` which specifies the maximum number of jobs submitted at the same time, `default-resources` requested for each job and `resources`, which define the resource limits.
Where to save SLURM logs and what to call them is also specified. Note that this folder must already exist. All of these arguments can be modified in the file `slurm/config.yaml` by the user depending on 
the desired ressources to allocate to the pipeline. To run the pipeline according to a specific profile, the following command line is used:

```
snakemake --use-conda --profile slurm
```

Snakemake can also use generic cluster support like `qsub`, by giving an other argument to `--cluster` combined with ressources. The profile file can also be altered accordingly depending on user needs. The most simple command
line in that setting is:

```
snakemake --use-conda --cluster qsub --jobs 22
```

## Test data

This pipeline has been tested using allele frequencies from the genomes of the [gnomAD database](https://gnomad.broadinstitute.org/downloads/) v2 and [CADD v1.6](https://cadd.gs.washington.edu/download) both in GRCh37, 
with the default parameters of the pipeline. Bed files for genes and CADD regions are available with the pipeline as well. Both of these files are in open access, and can be used to test the pipeline and its functionalities.

___