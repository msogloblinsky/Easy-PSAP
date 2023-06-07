#!/bin/Rscript
arg <- commandArgs(trailingOnly=T)
suppressPackageStartupMessages(library(data.table))

outfile <- arg[1]
outdir <- arg[2]
variants_exclude <- arg[3]

annotations <-fread(paste(outdir,"/annotated/",outfile,".vep.tsv.gz",sep=""),sep="\t",stringsAsFactors=F,header=TRUE,quote = "", data.table=F)
exclude <- fread(variants_exclude, data.table=FALSE)
annotations$Key <- paste0(annotations$Location,"_",annotations$REF_ALLELE,">",annotations$Allele)
exclude$Key <- paste0(exclude$V1,":",exclude$V3,"_",exclude$V4,">",exclude$V5)

print(paste0("Total number of annotated variants : ",dim(annotations)[1] ))
print(paste0("Number of not PASS variants in gnomAD genome removed : ",length(which(annotations$Key %in% exclude$Key))))

annotations_filt <- annotations[-which(annotations$Key %in% exclude$Key),]
fwrite(annotations_filt[, -grep("Key", names(annotations_filt))], paste(outdir,"/annotated/filtered_",outfile,".vep.tsv.gz",sep=""),compress = "gzip", col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)