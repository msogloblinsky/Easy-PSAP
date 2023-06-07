#!/usr/bin/env bash

mkdir -p "${snakemake_params[outdir]}"/annotated

INPUT_FILE="${snakemake_input[filtered_vcf]}"
VEP_CACHE_PATH="${snakemake_params[vep_cache]}"
VEP_CACHE_MERGED="${snakemake_params[vep_cache_merged]}"
VEP_FASTA="${snakemake_params[vep_fasta]}"
OUTDIR="${snakemake_params[outdir]}"
OUTFILE="${snakemake_params[outfile]}"
UNIT="${snakemake_params[unit]}"
CADD_PATH="${snakemake_params[cadd_path]}"
CADD_VERSION="${snakemake_params[cadd_version]}"
GENES="${snakemake_params[genes]}"
CADD_REGIONS="${snakemake_params[cadd_regions]}"
ASSEMBLY="${snakemake_params[assembly]}"
species="homo_sapiens"

if [ $VEP_CACHE_MERGED = "TRUE" ]; then cache="--cache --merged --offline --dir $VEP_CACHE_PATH" ; fi
if [ $VEP_CACHE_MERGED = "FALSE" ]; then cache="--cache --offline --dir $VEP_CACHE_PATH" ; fi

fasta="$VEP_FASTA"
annotations="--show_ref_allele --variant_class --symbol  --plugin CADD,$CADD_PATH/CADD_v$CADD_VERSION/whole_genome_SNVs.tsv.gz,$CADD_PATH/CADD_v$CADD_VERSION/InDels.tsv.gz --custom $GENES,Gene_Ensembl,bed,overlap,0 --custom $CADD_REGIONS,CADD_region,bed,overlap,0"
speed="--fork ${snakemake[threads]} --buffer_size ${snakemake_resources[mem]}" 
options="$cache $speed --species $species --assembly $ASSEMBLY --fasta $fasta $annotations --use_given_ref --tab --compress_output bgzip"

if [ $UNIT = "gene" ]; then vep $options --pick_allele_gene -i ${INPUT_FILE} -o ${OUTDIR}/annotated/${OUTFILE}.vep.tsv.gz ; fi 2> "${snakemake_log[0]}"
if [ $UNIT = "cadd_region" ] || [ $UNIT == "coding_cadd_region" ] ; then vep $options --pick_allele -i ${INPUT_FILE} -o ${OUTDIR}/annotated/${OUTFILE}.vep.tsv.gz ; fi 2> "${snakemake_log[0]}"
