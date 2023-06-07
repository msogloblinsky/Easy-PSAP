#!/bin/Rscript
arg <- commandArgs(trailingOnly=T)
vcf <- arg[1]
ped <- arg[2]
outfile <- arg[3]
outdir <- arg[4]

suppressPackageStartupMessages(library(RAVAQ))

ped_file <- read.table(ped, stringsAsFactors=F)
if(length(unique(ped_file[,6])) > 1 ){groupIDlist= split( ped_file[,1] , f = ped_file[,6] ) }else{groupIDlist = ped_file[,1]}

ravaq( vcffile=vcf, groupIDlist=groupIDlist, outdir=outdir, outprefix=outfile,
MAX_AB_GENO_DEV = 0.25, MAX_ABHET_DEV = "disabled",  MIN_CALLRATE = "disabled", MIN_FISHER_CALLRATE = "disabled",
step=1 )

ravaq( vcffile=paste0(outdir,"/",outfile,".step1.vcf.gz"), groupIDlist=groupIDlist, outdir=outdir, outprefix=outfile,
MAX_AB_GENO_DEV = 0.25, MAX_ABHET_DEV = "disabled",  MIN_CALLRATE = "disabled", MIN_FISHER_CALLRATE = "disabled",
step=3 )
