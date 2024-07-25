arg <- commandArgs(trailingOnly=TRUE)

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(splitstackshape))

outdir <- arg[1]
af_prefix <- arg[2]
units <- c("gene","coding_cadd_region","cadd_region")

for(i in 1:length(units)){

coverage_units <- read.table(paste0(outdir,"/reference_panel/coverage_",units[i],"_",af_prefix,".bed"))
suppressWarnings(coverage_units_split <- cSplit(coverage_units, "V5", sep=",", direction="long"))
coverage_units_split_group <- coverage_units_split %>% group_by(V5) %>% summarize(f = mean(V1))
list_units <- data.frame(unlist(coverage_units_split_group[which(coverage_units_split_group$f>0.5),"V5"]))
write.table(list_units, paste0(outdir,"/reference_panel/list_goodcoverage_50pct_",units[i],"_",af_prefix,".bed"), col.names=F, row.names=F, quote=F)
}
