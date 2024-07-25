# Snakemake workflow to apply PSAP null distributions 

## Introduction
This second workflow **snakemake_apply_PSAP** applies these null distributions to a vcf file of a patient or multiple patients or controls. 
The input data is highly customizable and the use of Snakemake makes it possible to run the pipeline multiple times with different input parameters without running again the steps creating files that have already been created.

## Usage 

Once the conda environmenet with Snakemake installed is setup and activated, the user can modify the configuration file according to the desired input files and parameters, and run the pipeline.
The VEP software is used for the vcf file annotation, without the need to install VEP on the machine thanks to the conda package manager. The necessary files for annotation are (cf next paragraph):
* VEP cache and FASTA (for the version 107 of VEP, which is the one currently used in the pipeline)
* CADD score files for SNVs and InDels

### 1. Configure the pipeline

#### config.yaml file

**snakemake_apply_PSAP/config/config.hg19.yaml** contains all pipeline parameters which are tuned by the user with example parameters for the GRCh37 assembly. **snakemake_apply_PSAP/config/config.hg38.yaml** contains example parameters to run the pipeline for data in GRCh38. Available parameters are as follows:

|**Parameter** | **Description** |
|--------------|--------------|
|**snakemake_directory**|Path to `/snakemake_apply_PSAP` directory on the user's machine|
|**vcf**|vcf file with individual data to score with PSAP, needs to be in format GRCh37 (no "chr" prefix for chromosomes) and bgzipped|
|**ped**|Corresponding PED file for the vcf file, tab-delimited, important columns are 2nd column = individual IDs, 5th column = sex (1=male, 2=female, other=unknown) and 6th column = status (1=unaffected, 2=affected)|
|**coverage**|bed file of regions well-covered regions in allele frequency database|
|**lookup_namefile**|Path and prefix of the lookup table for PSAP null distributions (lead to two files ending by "_het.txt.gz" and "_hom.txt.gz")|
|**variants_exclude**|List of variants to exclude from the vcf file, that were excluded from the calculation of null distributions. The file for gnomAD V2 is provided with the pipeline|
|**cadd_path**|Directory where the CADD files for annotation are located|
|**score_max**|Maximal possible value of the pathogenicity score (70 for CADD) used to calibrate null distributions|
|**genes**|bed file for the coordinates of genes coding regions. The file for GRCh37 is provided with the pipeline and a new one can be created using [BioMart](https://www.ensembl.org/info/data/biomart/index.html)|
|**cadd_regions**|bed file for the coordinates of CADD regions. The file for GRCh37 is provided with the pipeline or can be downloaded [here](https://lysine.univ-brest.fr/RAVA-FIRST/)|
|**vep_cache**|Directory of the VEP cache. Instructions on how to download the cache can be found on the [Ensembl website](https://www.ensembl.org/info/docs/tools/vep/script/vep_cache.html#cache) with the necessary files located [here](https://ftp.ensembl.org/pub/release-107/variation/indexed_vep_cache/) for VEP 107|
|**vep_cache_merged**|"TRUE" if VEP cache is merged, otherwise "FALSE"|
|**vep_fasta**|VEP FASTA file with its path which can be downloaded [here](https://ftp.ensembl.org/pub/release-107/fasta/) for VEP 107|
|**outdir**|Output directory|
|**outfile**|Prefix for output files, will be the name of the PSAP null distributions file|
|**unit**|Unit of testing to construct PSAP null distributions, can take the values "gene", "cadd_region" or "coding_cadd_region"|
|**cadd_version**|Version of CADD to be used for annotation|
|**assembly**|Both "GRCh37" and "GRCh38" options are currently supported|
|**compound_heterozygote_model**|TRUE or FALSE, if TRUE will calculate PSAP p-values for the compound heterozygote (CHET) model|
|**hemizygote_model**|TRUE or FALSE, if TRUE will calculate PSAP p-values for the hemizygote model for males depending on the sex indicated in the ped file (HEM PSAP null distributions pre-calculated for genes only)|
|**indel_file** File with InDels annotated by CADD website or "NA". Column names of InDel file must be "#Chrom Pos Ref Alt RawScore PHRED" (format of CADD 1.6 InDel file)|

> For the hemizygote model: this option will calculate hemizygote model for males and AD/AR/CHET models for females on chrX. 

> TIP: the user running the pipeline needs to have permission to write in the folder where the input vcf file is located. If not, the vcf file needs to already be indexed as the first step of the pipeline involves indexing the input vcf file. Suggested command for indexing a vcf file, using bcftools :

```
bcftools index {vcf}
```

> WARNING: The input vcf file for the PSAP pipeline must only contain biallelic variants (multiallelic variants need to be split in different lines).
This process can be done using a tool like [VCFprocessor](https://lysine.univ-brest.fr/vcfprocessor/functions.html#splitmultiallelic) or during a Quality Control of the vcf file using a tool 
like the [R package RAVAQ](https://gitlab.com/gmarenne/ravaq). We provide a script optional/qc_vcf.R to perform the recommended QC anf formatting for the vcf file. The script can be run using the command:

```
Rscript qc_vcf.R {vcf} {ped} {outfile} {outdir}
```

Users also need to check the format of the vcf file. For the GRCh37 assembly, no "chr" prefix is expected before the name of the chromosome. For the GRCh38 assembly, a "chr" prefix is expected before the name of the chromosome. If not the case, the user needs to format their vcf file accordingly, as the pipeline expects this format for annotation purposes. 

For the coverage file, the command line used to get the bed file provided for gnomAD genome V2 (good coverage = 90% of individuals at dp10) was:

```
zcat gnomad.genomes.r2.0.1.coverage.txt.gz | tail -n+2 | awk '{print $1"\t"($2-1)"\t"$2"\t"$7}' > gnomad.genome.r2.0.1.dp10.bed
```

The chosen lookup table for PSAP null distributions has to match the "unit" and "cadd_version" specified in the configuration file. Already calculated lookup tables include : 
*latest_gnomadgen_string_ensembl_cadd1.6_af_nosing_lookup : genes
*latest_gnomadgen_string_coding_caddregions_cadd1.6_af_nosing_lookup : coding CADD regions
*latest_gnomadgen_string_caddregions_cadd1.6_af_nosing_lookup : CADD regions 
Available PSAP null distributions have used gnomAD V2 genome for allele frequencies calibration (AF specifically) and CADD v1.6 for the pathogenicity score, all in GRCh37.

For the CADD files, the path and version need to lead to the following files:

```
{cadd_path}/CADD_v{cadd_version}/whole_genome_SNVs.tsv.gz
{cadd_path}/CADD_v{cadd_version}/InDels.tsv.gz
```
The formats of the [CADD files](https://cadd.gs.washington.edu/download) are currently supported. Files need to be indexed.

The file with InDels anotated by CADD is optional. If no file is provided, InDels will not be scored with PSAP and only SNVs will be kept for the analysis.
The script optional/make_indel_file_forcadd.sh can be used to generate input file with InDels to upload to the CADD website in the correct format. The script can be run using the following command 
(bcftools and tabix need to be installed in the environement for this command to work):
```
make_indel_file_forcadd.sh {vcf} {outfile} {outdir}
```

> TIP: different analyses can be run using the same repository. This is achieved by changing the outdir and outfile in the configuration file. Also different parameters values can be used in the different analyses.

### 2. Create the Conda environments

The Snakemake relies on conda environements which contain all the necessary dependencies to run the pipeline. The conda environements need to be installed first, which can take some time.
This step is easily done using the following command line:

```
snakemake --use-conda --conda-create-envs-only --conda-frontend conda -j 1 --configfile config/config.hg19.yaml
```

### 3. Run the pipeline.

Once the pipeline is configured and conda environments are created, the user just needs to run the Snakemake pipeline to score the vcf file with PSAP:

```
snakemake --use-conda --conda-frontend conda -j 20  --configfile config/config.hg19.yaml
```

The mandatory arguments are:
* **--use-conda**: to use the conda environments created at the previous step. The **--conda-frontend conda** specifies the use of conda to run the environments, instead of mamba.
* **-j**: number of threads/jobs provided to snakemake. The VEP annotation uses 20 threads if provided (not more), and the number of threads will condition the number of individuals analyzed at the same time. 

An optional argument **--dry-run** can be used to make a dry run of the pipeline, check if there are no warnings and see all the files that will be created.
With the additionnal argument **--configfile=/path_to_configfile/name_config.yaml**, the user can use a configuration file in a different location instead of modifying the default config.yaml file in the snakemake directory.

> TIP: If the execution of a step of the pipeline fails or if you want to see how it progresses, you can check the `{outdir}/log` directory where a log file is created for each rule of the Snakemake by individual if applicable.
> A `{outdir}/benchmark` directory is also created automatically to monitor how much memory and cpu is used by each step of the pipeline.

## Pipeline steps

The steps of the pipeline are not meant to be used separately and will be run by Snakemake in the most efficient way depending on input and output files.

### 1. filter_regions

The vcf file is filtered to exclude regions not well-covered in the database used to construct PSAP null distribution. This step ensures the reliability of PSAP results, so that variants analyzed in the vcf
were also present in the construction of PSAP null distributions. Output located in `{outdir}`.

### 2. write_column_names

Column names from the vcf file are written in a separate file. This file will then be used during PSAP annotation to check fields present in input files. Output located in `{outdir}`.

### 3. vep_annotation

During this step, the vcf file is annotated by the VEP software. Multiple annotations critical for downstream PSAP analysis are added, among which: the corresponding gene and CADD region, if applicable, and the CADD score of variants.
VEP cache and FASTA file need to be downloaded prior to the analysis, but the software itself is run through the conda environment. Annotation is slightly different depending on the unit of testing. Output located in `{outdir/annotated}`.

### 4. filter_variants

Variants that were filtered out in the database used to calculate PSAP null distributions are filtered out from the VEP output file in this step. Again, this ensures the reliability of PSAP results. Output located in `{outdir/annotated}`.

### 5. apply_PSAP_calculations

This step carries out the main function of the pipepline, which is the calculation of PSAP p-values for variants in the vcf input file. These calculations are run by individual. Preprocessing steps harmonize the data between the vcf file
and the VEP annotated file using the [gaston R package](https://cran.r-project.org/package=gaston). If an InDel CADD score file is provided, InDels are included in the analysis. Otherwise only autosomal SNVs are kept. Then, the variant with the maximal
CADD score by unit of testing (gene, CADD region or coding CADD region) is scored with the PSAP null distribtion for the autosomal dominant (heterozygote variants) or recessive model (homozygote variants). Additional models will be run if indicated in the configuration file, and can include the compound heterozygote and/or hemizygote models. Output located in `{outdir/annotated}`.

### 6. make_report_file

If affected individuals are present in the ped file (value 2 in the 6th column), this step makes a report file that merges the individual PSAP records. If a variant is present in multiple individuals, it will be a single record in the report file.
If there are control individuals in the vcf (value 1 in the 6th column of the ped file), a validation step is carried out. If the variant of the affected individual(s) is not seen in controls, the validation column has an "no_controls" value,
otherwise it has the value "violation". Output located in `{outdir}/annotated}`.

### 7. compress_output_files

Output files from PSAP are bgzipped to minimize the space taken by PSAP results, especially if a large number of individuals was analyzed. Report file is not compressed. Output located in `{outdir/annotated}`.

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
snakemake --use-conda --cluster qsub --jobs 20
```

## Test data

A test vcf file, ped file and InDel file are available in the `/example` folder. They correspond to the names of the files in the configuration file (path need to be altered). They can be used as an example to format the user's files 
or to test run the pipeline. The vcf file has been generated by sampling 10 individuals of European descent from the [1000 Genomes Project](https://www.internationalgenome.org/data), restricting to regions covered by the exome, 
and inserting known best-reviewed coding pathogenic variants (one in each individual exome) from the [ClinVar database](https://www.ncbi.nlm.nih.gov/clinvar/). The list of inserted variants and their genotypes can be found in 
the file `added_pathogenic_mutations_AD_AR_1kG.txt`. Running PSAP on these synthetic exomes shows how the method performs in a real-life setting.

___