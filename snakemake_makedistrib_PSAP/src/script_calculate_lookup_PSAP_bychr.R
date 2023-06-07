args <- commandArgs(trailingOnly=TRUE)

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(data.table))

score_tables <- args[1]
allele_frequencies_files <- args[2]
notpass_variants <- args[3]
units <- args[4]
outdir <- args[5]
outfile <- args[6]
chr <- as.numeric(args[7])

timeT = Sys.time()
options("digits" = 10)
setDTthreads(20)

get_cadd_het_final <- function(table_data){
#Group by unique CADD, with the probability of at least 1 or 2 heterygotes or 1 homozygote for duplicated CADDs
y <- table_data %>%                                     
  dplyr::group_by(cadd_phred) %>%
  dplyr::summarise(proba_het1 =  1 - prod(1-proba_het), 
  proba_het2 = ifelse(n() ==1, 0, 1 - prod(1-proba_het) - sum(proba_het*(prod(1-proba_het)/(1-proba_het))) ), 
  proba_hom1 = 1 - prod(1-proba_hom), n=n() ) %>% 
  as.data.frame()


#Ordering CADDs by descending order and calculating the proability for each to have no heterozygote or no homozygote at all
y2 <- y[order(-y$cadd_phred),]
y2$nohet <- 1- y2$proba_het1
y2$nohom <- 1- y2$proba_hom1


#Calculation of the probability to be the maximum CADD for the dominant model (a) and recessive model (d)
res1 <- lapply(1:nrow(y2),function(i){
						a = ifelse(i==1, y2$proba_het1[i],y2$proba_het1[i] * prod(y2$nohet[1:(i-1)]) ) ;
						d = ifelse(i==1, y2$proba_hom1[i],y2$proba_hom1[i] * prod(y2$nohom[1:(i-1)]) ) ;
						return(c(a,d))
						}
						)
res_mat <- matrix(unlist(res1),nrow=2)

#Distribution of maximum CADDs for the heterozygote model
scale = seq(0,70,0.05)
res_mat_het <- data.frame(t(rbind(res_mat[1,], t(y2$cadd_phred))))
colnames(res_mat_het) <- c("V1","V2")
res_mat_het$V1 <- res_mat_het$V1/(sum(res_mat_het$V1) )

df1 <- res_mat_het %>%
 dplyr::mutate(bin = cut(V2, breaks = scale, right=FALSE)) %>%
 dplyr::group_by(bin) %>%
 dplyr::summarize(prob = sum(V1) ) %>%
 dplyr::ungroup() %>%
 tidyr::complete(bin, fill = list(prob = 0)) 

as.data.frame(t(cumsum(df1$prob[1400:1])[1400:1]))
}


get_cadd_hom_final <- function(table_data){
y <- table_data %>%                                     
  dplyr::group_by(cadd_phred) %>%
  dplyr::summarise(proba_het1 =  1 - prod(1-proba_het), 
  proba_het2 = ifelse(n() ==1, 0, 1 - prod(1-proba_het) - sum(proba_het*(prod(1-proba_het)/(1-proba_het))) ), 
  proba_hom1 = 1 - prod(1-proba_hom), n=n() ) %>% 
  as.data.frame()

y2 <- y[order(-y$cadd_phred),]
y2$nohet <- 1- y2$proba_het1
y2$nohom <- 1- y2$proba_hom1

res1 <- lapply(1:nrow(y2),function(i){
						a = ifelse(i==1, y2$proba_het1[i],y2$proba_het1[i] * prod(y2$nohet[1:(i-1)]) ) ;
						d = ifelse(i==1, y2$proba_hom1[i],y2$proba_hom1[i] * prod(y2$nohom[1:(i-1)]) ) ;
						return(c(a,d))
						}
						)
res_mat <- matrix(unlist(res1),nrow=2)

#Distribution of maximum CADDs for the homozygote model
scale = seq(0,70,0.05)
res_mat_hom <- data.frame(t(rbind(res_mat[2,], t(y2$cadd_phred))))
colnames(res_mat_hom) <- c("V1","V2")
res_mat_hom$V1 <- res_mat_hom$V1/(sum(res_mat_hom$V1) )


df1 <- res_mat_hom %>%
 dplyr::mutate(bin = cut(V2, breaks = scale, right=FALSE)) %>%
 dplyr::group_by(bin) %>%
 dplyr::summarize(prob = sum(V1) ) %>%
 dplyr::ungroup() %>%
 tidyr::complete(bin, fill = list(prob = 0)) 

as.data.frame(t(cumsum(df1$prob[1400:1])[1400:1]))
}


print("Read in data")

timeF = Sys.time()
timeT = Sys.time()
data <- fread(score_tables, data.table=FALSE)
setnames(data, names(data), c("chr","start","end","ref","alt","cadd_raw","cadd_phred","cadd_region","gene"))

if(units == "gene"){data <- filter(data, gene != ".") ; column_testing="gene"
}else if(units == "coding_caddregion"){data <- filter(data, gene != ".") ; data <- data[,-9] ; data <- data[!duplicated(data),] ; column_testing="cadd_region"
}else if(units == "caddregion"){data <- data[,-9] ; data <- distinct_all(data, .keep_all = TRUE) ; column_testing="cadd_region"}

print(Sys.time()-timeT)
print("Filtering data")

timeT = Sys.time()
#Exclude variants in regions not well-covered : already done 
#Exclude not pass var
exclude <- fread(notpass_variants, data.table=FALSE)
names(exclude) <- c("chr","start","ref","alt")
data.exclude <- anti_join(data, exclude, by = c("chr","start","ref","alt")) 

#Get allele frequencies
refpanel_af <- fread(allele_frequencies_files, data.table=FALSE, nThread=40)
setnames(refpanel_af, names(refpanel_af), c("CHROM","POS","REF","ALT","AF","AC"))
print(Sys.time()-timeT)

#Calculation of probabilities for each alternative allele
timeT = Sys.time()
x <- left_join(data.exclude,refpanel_af, by=c("start" = "POS", "alt" = "ALT"))
x$proba <- as.numeric(ifelse(is.na(x$AF), 0, x$AF))

#Calculation of probabilities to be heterozygous or homozygous for each allele
x$proba_het <- 2*x$proba*(1-x$proba)
x$proba_hom <- x$proba**2
x$cadd_phred <- round(x$cadd_phred,3)
print(Sys.time()-timeT)

# Direct calculation of null distributions by chr
print("Calculating null distributions")
timeT = Sys.time()
cadd_het_final <- setDT(x)[,get_cadd_het_final(.SD),column_testing]
print(Sys.time()-timeT)
print("Model HET done")

timeT = Sys.time()
cadd_hom_final <- setDT(x)[,get_cadd_hom_final(.SD),column_testing]
print(Sys.time()-timeT)
print("Model HOM done")

timeT = Sys.time()
print("Writing output files")
fwrite(cadd_het_final, file = paste0(outdir,"/temp_lookup_bychr/",outfile,"_",units,"_het_chr",chr,".txt"), quote=F, sep="\t", row.names=F, col.names=F)
fwrite(cadd_hom_final, file = paste0(outdir,"/temp_lookup_bychr/",outfile,"_",units,"_hom_chr",chr,".txt"), quote=F, sep="\t", row.names=F, col.names=F)

print(paste0("Chromosome ",chr," done")) 
print(Sys.time()-timeF)