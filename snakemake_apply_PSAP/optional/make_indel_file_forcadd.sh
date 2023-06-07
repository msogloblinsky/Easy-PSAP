#!/bin/bash

VCF=$1 
OUTDIR=$2
OUTFILE=$3

bcftools view -G --no-header --types indels ${VCF} > ${OUTDIR}/${OUTFILE}.indels.vcf
split -l 80000 --numeric-suffixes ${OUTDIR}/${OUTFILE}.indels.vcf ${OUTDIR}/${OUTFILE}.indels --additional-suffix=.vcf
for file in ${OUTDIR}/${OUTFILE}.indels*.vcf ; do bgzip ${{file}} ; done