#!/bin/bash

VCF=$1 
COVERAGE=$2
OUTDIR=$3
OUTFILE=$4

# FILTER REGIONS : GOOD QUALITY IN DATABASE +/- CODING AND SPLICING REGIONS
	bcftools index -f $VCF
	bcftools view  -R $COVERAGE $VCF -Oz > ${OUTDIR}/${OUTFILE}_filtered_coverage.vcf.gz
	bcftools sort ${OUTDIR}/${OUTFILE}_filtered_coverage.vcf.gz -Oz > ${OUTDIR}/${OUTFILE}_sortedfiltered_coverage.vcf.gz
	bcftools index ${OUTDIR}/${OUTFILE}_sortedfiltered_coverage.vcf.gz
	
