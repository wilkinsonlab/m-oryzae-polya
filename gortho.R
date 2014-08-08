source("/home/marco/Downloads/gProfileR/R/gProfileR.R")
args <- commandArgs(TRUE)
t<-scan(args[1], what="character")
s<-as.character(args[2])
d<-as.character(args[3])
r<-gorth(t, source_organism=s, target_organism=d)
r<-r[order(r$alias.number),]
write.table(file=paste(args[1], "_orthologs_", d, sep=""), r$target.ensg, quote=F, row.names=F, col.names=F)
