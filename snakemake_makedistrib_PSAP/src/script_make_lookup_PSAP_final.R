args <- commandArgs(trailingOnly=TRUE)

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(data.table))

list_units <- args[1]
units <- args[2]
outdir <- args[3]
outfile <- args[4]

#Make final lookup table
keep <- fread(list_units, data.table=F, nThread=40) 
keep <- keep %>% unlist %>% as.vector

## HET
lookup_het <- data.frame()
for(i in 1:22){
data.tmp <- fread(paste0(outdir,"/temp_lookup_bychr/",outfile,"_",units,"_het_chr",i,".txt"), data.table=F, nThread=40)
lookup_het <- rbind.data.frame(lookup_het,data.tmp)
}
lookup_het <- filter(lookup_het, V1 %in% keep)
lookup_het <- lookup_het %>% group_by(V1) %>% filter(row_number() == 1)
fwrite(lookup_het, paste0(outdir,"/final_lookuptables_PSAP/",outfile,"_",units,"_het.txt"), nThread=40, quote=F, sep="\t", row.names=F, col.names=F)


## HOM
lookup_hom <- data.frame()
for(i in 1:22){
data.tmp <- fread(paste0(outdir,"/temp_lookup_bychr/",outfile,"_",units,"_hom_chr",i,".txt"),, data.table=F, nThread=40)
lookup_hom <- rbind.data.frame(lookup_hom,data.tmp)
}
lookup_hom <- filter(lookup_hom, V1 %in% keep)
lookup_hom <- lookup_hom %>% group_by(V1) %>% filter(row_number() == 1)
fwrite(lookup_hom, paste0(outdir,"/final_lookuptables_PSAP/",outfile,"_",units,"_hom.txt"), nThread=40, quote=F, sep="\t", row.names=F, col.names=F)
