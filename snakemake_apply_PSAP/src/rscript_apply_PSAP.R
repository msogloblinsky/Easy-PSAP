## Calculate PSAP p-values for AD and AR models
tmp$Dz.Model <- NA
 
##For AD 
exome.AD <- subset(tmp, Geno == "het")
if (nrow(exome.AD)>0) {
	exome.AD$Dz.Model[which(exome.AD$Geno == "het")] <- "DOM-het"
	lookup <-fread(file=paste(lookup_namefile,"_het.txt.gz",sep=""),stringsAsFactors=FALSE,header=FALSE,data.table=FALSE)
	exome.AD$j<-findInterval(exome.AD[,score],scale)+1
	if(length(which(exome.AD$j == 1402)) > 0){ exome.AD[which(exome.AD$j == 1402),] <- 1401}
	exome.AD$i<-as.integer(factor(exome.AD[,unit_column],levels=lookup[,1]))
	exome.AD$PSAP<-unlist(apply(exome.AD,1,function(x,tab) {i<-as.numeric(x["i"]); j<-as.numeric(x["j"]); y<-tab[i,j]; return(y)},lookup))
	tmp.AD <- exome.AD %>% group_by(!!as.name(unit_column)) %>% filter(!!as.name(score)	== max(!!as.name(score))) %>% distinct(!!as.name(unit_column), .keep_all= TRUE) %>% arrange(!!as.name(unit_column)) %>% as.data.frame()
	out <- rbind(out, tmp.AD)	
	rm(list = c("exome.AD", "lookup"))
}
print("AD model complete")

##For AR 
exome.AR <- subset(tmp, Geno == "hom")
if (nrow(exome.AR)>0) {
	exome.AR$Dz.Model[which(exome.AR$Geno == "hom")] <- "REC-hom"
	lookup <-fread(file=paste(lookup_namefile,"_hom.txt.gz",sep=""),stringsAsFactors=FALSE,header=FALSE,data.table=FALSE)
	exome.AR$j<-findInterval(exome.AR[,score],scale)+1
	if(length(which(exome.AR$j == 1402)) > 0){ exome.AR[which(exome.AR$j == 1402),] <- 1401}
	exome.AR$i<-as.integer(factor(exome.AR[,unit_column],levels=lookup[,1]))
	exome.AR$PSAP<-unlist(apply(exome.AR,1,function(x,tab) {i<-as.numeric(x["i"]); j<-as.numeric(x["j"]); y<-tab[i,j]; return(y)},lookup))
	tmp.AR <- exome.AR %>% group_by(!!as.name(unit_column)) %>% filter(get(score)	== max(get(score))) %>% distinct(!!as.name(unit_column), .keep_all= TRUE) %>% arrange(!!as.name(unit_column)) %>% as.data.frame()
	out <- rbind(out, tmp.AR)
	rm(list = c("exome.AR", "lookup"))
}
print("REC-hom model complete") 



