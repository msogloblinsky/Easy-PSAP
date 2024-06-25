## Calculate PSAP p-values for AD and AR models
tmp$Dz.Model <- NA
 
##For AD 
exome.AD <- subset(tmp, Geno == "het")
if (nrow(exome.AD)>0) {
	exome.AD$Dz.Model[which(exome.AD$Geno == "het")] <- "DOM-het"
	if(hem_model){exome.AD$Dz.Model[which(exome.AD$Geno == "het_female")] <- "X-linked-het"}
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
	if(hem_model){exome.AR$Dz.Model[which(exome.AR$Geno == "hom_female")] <- "X-linked-hom"}
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

##For CHET
if(chet_model){
	exome.CHET <- subset(tmp, Geno == "het") 
	dups <- names(table(exome.CHET[,unit_column])[which(table(exome.CHET[,unit_column]) >= 2)])
	exome.CHET <- exome.CHET[exome.CHET[,unit_column] %in% dups,]
	if (nrow(exome.CHET) > 0) {
		exome.CHET$Geno[which(exome.CHET$Geno == "het")] <- "chet"
		exome.CHET$Dz.Model[which(exome.CHET$Geno == "chet")] = "REC-chet"
		if(hem_model){exome.CHET$Geno[which(exome.CHET$Geno == "het_female")] <- "chet_female"
		exome.CHET$Dz.Model[which(exome.CHET$Geno == "chet_female")] = "X-linked-chet"}
		lookup <-fread(file=paste(lookup_namefile,"_chet.txt.gz",sep=""),stringsAsFactors=FALSE,header=FALSE,data.table=FALSE)
		exome.CHET$j<-findInterval(exome.CHET[,score],scale)+1
		exome.CHET$i<-as.integer(factor(exome.CHET[,unit_column],levels=lookup[,1]))
		exome.CHET$PSAP<-unlist(apply(exome.CHET,1,function(x,tab) {i<-as.numeric(x["i"]); j<-as.numeric(x["j"]); y<-tab[i,j]; return(y)},lookup))
		tmp.CHET <- do.call(rbind,by(exome.CHET,exome.CHET[,unit_column],function(x,score) { x[,score]<-as.numeric(x[,score]); x<-x[order(x[,score],decreasing=T,na.last=NA),]; return(x[1:2,]) },score=score)) 
		units <- unique(tmp.CHET[,unit_column])
		for (unite in units){
			tmp.CHET$PSAP[which(tmp.CHET[,unit_column] == unite)] <- tmp.CHET$PSAP[which(tmp.CHET[,unit_column] == unite)][2]	
		}			
		out <- rbind(out, tmp.CHET)
		rm(list = c("exome.CHET", "lookup"))
	}
	print("REC-chet model complete")
}

# For male: additional mode HEM
if (gender == "male" & hem_model == TRUE){
	exome.HEM <- subset(tmp, Geno == "hem") # for HEM mode
	if (nrow(exome.HEM) > 0) {
		exome.HEM$Dz.Model[which(exome.HEM$Geno == "hem" & (exome.HEM$Chr == "X" | exome.HEM$Chr == "chrX"))] <- "X-linked-hem"
		lookup <- fread(file=paste(lookup_namefile,"_XY_hem_chrX.txt.gz",sep=""),stringsAsFactors=FALSE,header=FALSE,data.table=FALSE,select=1)
		exome.HEM$j<-findInterval(exome.HEM[,"score_norm"],scale)+1
		exome.HEM$i<-as.integer(factor(exome.HEM[,gene_col],levels=lookup[,1]))
		exome.HEM$PSAP<-unlist(apply(exome.HEM,1,function(x,tab) {i<-as.numeric(x["i"]); j<-as.numeric(x["j"]); y<-tab[i,j]; return(y)},lookup))
		tmp.HEM <- exome.HEM %>% group_by(!!as.name(unit_column)) %>% filter(!!as.name(score)	== max(!!as.name(score))) %>% distinct(!!as.name(unit_column), .keep_all= TRUE) %>% arrange(!!as.name(unit_column)) %>% as.data.frame()
		out <- rbind(out, tmp.HEM)	
	}
print("HEM model complete")
}

