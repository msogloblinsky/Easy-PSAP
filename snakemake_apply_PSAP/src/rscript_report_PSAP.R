# Make report file for AD and AR PSAP results for affected individuals
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dplyr))

arg = commandArgs(trailingOnly=T)
outfile = arg[1]
outdir = arg[3]
unit_testing = arg[4]
chet_model <- ifelse(arg[5] == "TRUE" | arg[5] == "True", TRUE, FALSE)
hem_model <- ifelse(arg[6] == "TRUE" | arg[6] == "True", TRUE, FALSE)
unit_column <- ifelse(unit_testing == "gene", "Gene_Ensembl", "CADD_region")

ped<-read.table(file=arg[2],header=F,stringsAsFactors=F, sep="\t")
af<-ped$V2[which(ped$V6==2)]
n.af<-length(af)
id<-af

if (sum(which(ped$V6==2)) == 0) {
stop("No affected individuals in the ped file\n Please see documentation for details\n Expecting value 2 in the 6th column for some or all individuals")
}

models = c("DOM-het","REC-hom")
if(chet_model) {models = c(models,"REC-chet")}
if(hem_model) {models = c(models,"X-linked-hom","X-linked-het","X-linked-hem")}
if(hem_model & chet_model) {models = c(models,"X-linked-chet")}
genos = c("het","hom")

print("analysing affected individuals")

af.dat = list()
for(i in 1:length(af)){
  dat<-fread(file=paste(outdir,"/annotated/",outfile,"_",af[i],"_PSAP.txt",sep=""),sep="\t",header=T,stringsAsFactors=F,check.names=F,quote="", data.table=F)
  dat$pid = af[i]
  print(af[i])
  
  for(m in 1:length(models)){
    tmp = filter(dat,Dz.Model == models[m])
    if(models[m] == "REC-chet" & nrow(tmp) > 0){
      a1 = filter(dat,Dz.Model == "DOM-het" & !!as.name(unit_column) %in% tmp[,unit_column])
      a1 = merge(a1[-which(names(a1) %in% c("PSAP","Dz.Model"))],tmp[c(unit_column,"PSAP","Dz.Model")])
      tmp = rbind(tmp,a1)
      tmp$chet.id=paste(tmp$Key,tmp$PSAP,sep=":")
    }      
    if(i == 1){
      af.dat[[m]] = tmp
    }else{
      af.dat[[m]] = rbind(af.dat[[m]],tmp)
    }
  }
}

candidates = data.frame()
for(m in 1:length(models)){
  print(models[m])
  if(models[m] == "REC-chet"){
	#tmp <- af.dat[[m]] %>% group_by(chet.id) %>% mutate(pid=paste(pid,collapse=","), n=n() ) %>% distinct(across(-c("pid","chet.id"))) %>% arrange(Key) %>% ungroup %>% select(-chet.id)
	tmp <- do.call(rbind,by(af.dat[[m]],af.dat[[m]]$chet.id, function(dat) return(data.frame(unique(dat[which(!names(dat) %in% c("pid","chet.id"))]),pid=paste(unique(dat$pid),collapse=","),n=length(unique(dat$pid)) ))))
  }else{
    tmp = do.call(rbind,by(af.dat[[m]],af.dat[[m]]$Key, function(dat) return(data.frame(unique(dat[which(names(dat) != "pid")]),pid=paste(unique(dat$pid),collapse=","),n=length(unique(dat$pid)) ))))
  }
  candidates = rbind(candidates,tmp)
}

if(length(which(ped$V6==1)) > 0){
	uf = ped$V2[which(ped$V6==1)]
	n.uf = length(uf)
	id<-c(af,uf)
	print("validating against unrelated individuals")
	uf.dat = data.frame()
	for(i in uf){
		dat<-fread(file=paste(outdir,"/annotated/",outfile,"_",i,"_PSAP.txt",sep=""),sep="\t",header=T,stringsAsFactors=F,check.names=F,data.table=F)
		uf.dat = rbind(uf.dat,dat)
	}
	if(unit_testing %in% c("cadd_region","coding_cadd_region")){
	uf.dat = unique(uf.dat[c("CADD_region","Chr","Start","REF_ALLELE","Allele","Dz.Model","Geno","Key","PSAP")]) }
	if(unit_testing == "gene"){
	uf.dat = unique(uf.dat[c("Gene_Ensembl","Chr","Start","REF_ALLELE","Allele","Dz.Model","Geno","Key","PSAP")]) }
  
  candidates$validation = "ok"
  candidates$validation[which(candidates$Dz.Model == "DOM-het" & candidates$Key %in% uf.dat$Key)] = "violation"
  candidates$validation[which(candidates$Dz.Model != "DOM-het" & candidates$Key %in% uf.dat$Key[which(uf.dat$Dz.Model != "DOM-het")])] = "violation"

if(chet_model){
	violations = filter(candidates,Dz.Model == "REC-chet" & validation == "violation")
        candidates$pid = as.character(candidates$pid)
        gene=unique(violations[,unit_column]) 
        candidates$pid_beforeupdateRECCHET = candidates$pid  
	for(g in gene){                                      
		ok.chets = which(candidates[,unit_column] == g & candidates$Dz.Model == "REC-chet" & candidates$validation == "ok")
		violated.pid = unique(unlist(strsplit(as.character(violations$pid[which(violations[,unit_column] == g)]),",",fixed=F)))
		for(p in violated.pid){  
			if(length(ok.chets) > 0){
				for(row in ok.chets){ 
					if(length(grep(",",candidates$pid[row])) > 0){
						pids = unlist(strsplit(candidates$pid[row],",",fixed=T))
						pids = pids[which(pids != p)]
						if(length(pids) > 0){ 
							candidates$pid[row] = paste(pids,collapse=",") 
						}else{
							candidates$pid[row] = NA
						}
					}else{
						pids = candidates$pid[row]
						pids = pids[which(pids != p)]     
						if(length(pids) > 0){ 
							candidates$pid[row] = pids 
						}else{
							candidates$pid[row] = NA
						}
					}
				}
			}
		}
	}
}else{
  candidates$validation = "no_controls"
}
}

print("validation complete")

validated=candidates 
rm(candidates)

print("generating report file")
if(unit_testing %in% c("cadd_region","coding_cadd_region")){
write.table(validated[order(validated$PSAP,validated$CADD_region,validated$Dz.Model),],file=paste(outdir,"/annotated/",outfile,".report.txt",sep=""),sep="\t",col.names=T,row.names=F,quote=F)}
if(unit_testing == "gene"){
write.table(validated[order(validated$PSAP,validated$Gene_Ensembl,validated$Dz.Model),],file=paste(outdir,"/annotated/",outfile,".report.txt",sep=""),sep="\t",col.names=T,row.names=F,quote=F) }


