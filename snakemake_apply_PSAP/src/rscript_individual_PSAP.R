#!/bin/Rscript
arg <- commandArgs(trailingOnly=T) 
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(gaston))
suppressPackageStartupMessages(library(splitstackshape))

outfile <- arg[1]
outdir <- arg[2]
indv.id <- arg[3]
dir <- arg[4]
ped <- read.table(arg[5],stringsAsFactors=F)
vcf.file <- arg[6]
lookup_namefile <- arg[7]
unit_testing <- arg[8]
merged <- arg[9]
indel.file <- ifelse(arg[10] == "NA", NA, arg[10])

if(file.exists(paste(outdir,"/annotated/",outfile,"_",indv.id,"_PSAP.txt",sep="")) == FALSE ){

gene_col <- "SYMBOL"
func_col <- "Consequence"
loc_col <- "VARIANT_CLASS"
aa_col <- "Amino_acids"
score <- "CADD_PHRED"
caddreg<- "CADD_region"
scale = seq(0,70,0.05)

unit_column <- ifelse(unit_testing == "gene", "Gene_Ensembl", "CADD_region")
lookup <-fread(file=paste(lookup_namefile,"_het.txt.gz",sep=""),stringsAsFactors=FALSE,header=FALSE,data.table=FALSE,select=1)
lookup.present <- as.vector(lookup[,1])

print("read in data")

annotations <-fread(paste(outdir,"/annotated/filtered_",outfile,".vep.tsv.gz",sep=""),sep="\t",stringsAsFactors=F,header=TRUE,quote = "", data.table=F)
vcf.header <- read.table(paste(outdir,"/",outfile,".header",sep=""), sep="\t",stringsAsFactors=F,comment.char="@")
stopifnot(indv.id %in% vcf.header) 

chr <- str_split_fixed(annotations$Location, ":", 2)
chrpos <- cbind(chr[,1],str_split_fixed(chr[,2], "-",2))
colnames(chrpos) <- c("Chr","Start","End")
annotations <- cbind.data.frame(annotations,chrpos)
annotations$Chr <- as.character(annotations$Chr) ; annotations$Start <- as.character(annotations$Start) ; annotations$End <- as.character(annotations$End)
 
options(scipen=999)

exome.raw<-read.vcf(vcf.file) ; exome.raw <- select.inds(exome.raw, id == indv.id)
exome.raw <- select.snps(exome.raw, N1 !=0 | N2 !=0)
exome.raw <- exome.raw@snps

#for deletions : add one bp, remove 1st character from ref and alt #for insertions remove 1st character from ref and alt #replace missing by -
insertions <- which(nchar(exome.raw[,"A1"]) < nchar(exome.raw[,"A2"]))
deletions <- which(nchar(exome.raw[,"A1"]) > nchar(exome.raw[,"A2"]))
exome.raw[insertions,"A2"] <- sub('.', '', exome.raw[insertions,"A2"]) ; exome.raw[insertions,"A1"] <- "-"
exome.raw[deletions,"pos"] <- exome.raw[deletions,"pos"] + 1 ; exome.raw[deletions,"A1"] <- sub('.', '', exome.raw[deletions,"A1"]) ; exome.raw[deletions,"A2"] <- "-"

annotations$Key <- paste(annotations$Chr,":",annotations$Start,"_",annotations$REF_ALLELE,">",annotations$Allele, sep = "")

exome.raw$pos <- as.character(exome.raw$pos)
exome.raw$ID <- paste(exome.raw$chr,":",exome.raw$pos,"_",exome.raw$A1,">",exome.raw$A2, sep = "")
exome.raw<-exome.raw[,c("ID","quality","filter","N1","N2")]

if(unit_testing == "gene"){
	if(merged == "TRUE"){annotations <- annotations[which(annotations$SOURCE == "Ensembl"),]}
	annotations.tmp <- suppressWarnings(cSplit(annotations, "Gene_Ensembl", ",", "long"))
	forannotations <- annotations.tmp[-which(annotations.tmp$SYMBOL != annotations.tmp$Gene_Ensembl | annotations.tmp$Gene_Ensembl == "-" ),]
}
if(unit_testing == "cadd_region"){
	forannotations <- suppressWarnings(cSplit(annotations, "CADD_region", ",", "long"))
}
if(unit_testing == "coding_cadd_region"){
	annotations.tmp <- suppressWarnings(cSplit(annotations, "CADD_region", ",", "long"))
	forannotations <- annotations.tmp[-which(annotations.tmp$Gene_Ensembl == "-"),]
}

tmp.annotations<-forannotations[which(forannotations$Key %in% exome.raw$ID),]
exome.raw.fin<-merge(tmp.annotations,exome.raw,by.x="Key",by.y="ID",all.x=T)

exome.raw.fin <- as.data.frame(exome.raw.fin)
exome.raw.fin <- exome.raw.fin[! duplicated(exome.raw.fin[,c("Key",unit_column)]),]

print("Extract genotypes")
exome.raw.fin[,indv.id] <- ifelse(exome.raw.fin$N2 == 1, "hom", "het")

print("clean data process")
exome.raw.fin[,score] <- as.numeric(exome.raw.fin[,score])
lookup.remove <- which(! exome.raw.fin[,unit_column] %in% lookup.present)

if(length(lookup.remove) != 0){
	exome <- exome.raw.fin[-unique(lookup.remove),]
}else{exome <- exome.raw.fin}

if (is.na(indel.file)==FALSE){
	indel.scores = read.table(indel.file, sep='\t', quote='', header=T, comment.char='',skip=1)

	if (any(! names(indel.scores)[1:3]==c("X.Chrom","Pos","Ref"))){
		stop("Error with CADD indel input file format\n Please see documentation for details\n Expecting second line to contain #Chrom Pos Ref Alt RawScore PHRED\n")
	}

	ichr<-grep("chr",indel.scores$X.Chrom)
	vchr<-grep("chr",exome$Chr)
	
	if (length(vchr)>0 & length(ichr)==0){
		indel.scores$X.Chrom<-paste("chr",indel.scores$X.Chrom,sep="")
	}
	
	adash<-grep("-",indel.scores$Alt)
	rdash<-grep("-",indel.scores$Ref)
	
	if (length(adash)==0 & length(rdash)==0){
		indel.scores$Ref <- as.character(indel.scores$Ref)
		indel.scores$Alt <- as.character(indel.scores$Alt)
		indel.scores$Pos <- as.integer(apply(indel.scores, 1, function(x) ifelse(nchar(x[3])>1,as.numeric(x[2])+1,x[2])))
		indel.scores$Ref <- apply(indel.scores, 1, function(x) ifelse(nchar(x[3])>1,substring(x[3],2),"-"))
		indel.scores$Alt <- apply(indel.scores, 1, function(x) ifelse(nchar(x[4])>1,substring(x[4],2),"-"))
		indel.scores <- unique(indel.scores[,c("X.Chrom","Pos","Ref","Alt","RawScore","PHRED")])
	}
	
	drops <- c("CADD_RAW","CADD_PHRED")
	indels = exome[which(exome[,loc_col] %in% c("deletion","insertion")), !names(exome) %in% drops]
	indels = merge(indels, indel.scores, by.x=c("Chr","Start","REF_ALLELE","Allele"), by.y=c("X.Chrom","Pos","Ref","Alt"), all.x=T, all.y=F)
	names(indels)[names(indels)=="RawScore"] <- "CADD_RAW"
	names(indels)[names(indels)=="PHRED"] <- score
	indels <- indels[,which(names(indels) %in% names(exome))]
	exome <- rbind(exome[which(! exome[,loc_col] %in% c("deletion","insertion")),], indels)
}

MAX <- 70
MIN <- 0
exome[,score][which(exome[,score] > MAX)] <- MAX
exome[,score][which(exome[,score] < MIN)] <- MIN
exome$score_norm <- (exome[,score]-MIN)/(MAX-MIN)

if(is.na(indel.file)){
exome <- exome[which(exome$VARIANT_CLASS == "SNV"),]
}

info <- exome[which(is.na(exome[,score]) == F) ,
c("Key","Chr","Start","Allele","REF_ALLELE","SYMBOL","Gene","Consequence","VARIANT_CLASS","IMPACT","Amino_acids","Protein_position","Gene_Ensembl","CADD_region","score_norm",score,indv.id)]

id.raw <- exome.raw.fin$Key
id.final <- info$Key
not_include <- unique(exome.raw.fin[which(! id.raw %in% id.final),
c("Key","Chr","Start","Allele","REF_ALLELE","SYMBOL","Gene","Consequence","VARIANT_CLASS","IMPACT","Amino_acids","Protein_position","Gene_Ensembl","CADD_region",score,indv.id)])

rm(list=c("exome","exome.raw.fin"))
tmp <- info[-which(names(info) == indv.id)]
tmp["Geno"] <- info[,which(names(info) == indv.id)]

print("data cleaned, beginning PSAP annotation")
out <- data.frame()
source(paste(dir,"/src/rscript_apply_PSAP.R",sep=""))
out <- out[which(!names(out) %in% c("i","j","score_norm"))]
out <- out[order(out$PSAP),]

print("writing file")
write.table(out,paste(outdir,"/annotated/",outfile,"_",indv.id,"_PSAP.txt",sep=""),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
}