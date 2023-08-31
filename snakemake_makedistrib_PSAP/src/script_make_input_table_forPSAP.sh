#!/bin/bash

cadd_file=$1
genes=$2
cadd_regions=$3
coverage=$4
outdir=$5
score_prefix=$6
chr=$7
assembly=$8          
	
tabix ${cadd_file} ${chr} | bgzip > ${outdir}/score_tables/${score_prefix}_chr${chr}.tsv.gz

zcat ${outdir}/score_tables/${score_prefix}_chr${chr}.tsv.gz | awk 'BEGIN{OFS="\t"}{print $1,$2,$2,$3,$4,$5,$6}' |
bedtools intersect -a stdin -b ${cadd_regions} -wa -wb | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4,$5,$6,$7,$11 }' | 
bedtools intersect -a stdin -b ${genes} -loj | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4,$5,$6,$7,$8,$12 }' | 
bedtools intersect -a stdin -b ${coverage} |
bgzip > ${outdir}/score_tables/coverage_genes_caddregions_${score_prefix}_chr${chr}.tsv.gz


if [[ ${assembly} -eq GRCh38 ]]
then
	zcat ${outdir}/score_tables/${score_prefix}_chr${chr}.tsv.gz | awk 'BEGIN{OFS="\t"}{print "chr"$1,$2,$2,$3,$4,$5,$6}' |
	bedtools intersect -a stdin -b ${cadd_regions} -wa -wb | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4,$5,$6,$7,$11 }' | 
	bedtools intersect -a stdin -b ${genes} -loj | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4,$5,$6,$7,$8,$12 }' | 
	bedtools intersect -a stdin -b ${coverage} |
	bgzip > ${outdir}/score_tables/coverage_genes_caddregions_${score_prefix}_chr${chr}.tsv.gz
elif [[ ${assembly} -eq GRCh37 ]]
then
	zcat ${outdir}/score_tables/${score_prefix}_chr${chr}.tsv.gz | awk 'BEGIN{OFS="\t"}{print $1,$2,$2,$3,$4,$5,$6}' |
	bedtools intersect -a stdin -b ${cadd_regions} -wa -wb | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4,$5,$6,$7,$11 }' | 
	bedtools intersect -a stdin -b ${genes} -loj | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4,$5,$6,$7,$8,$12 }' | 
	bedtools intersect -a stdin -b ${coverage} |
	bgzip > ${outdir}/score_tables/coverage_genes_caddregions_${score_prefix}_chr${chr}.tsv.gz
fi