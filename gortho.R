source("/home/marco/Downloads/gProfileR/R/gProfileR.R")
args <- commandArgs(TRUE)
t<-read.table(args[1])
r<-gorth(t[,1], source_organism="moryzae", target_organism="scerevisiae")
r<-r[order(r$alias.number),]
write.table(file=paste(args[1], "_yeast", sep=""), r$target.ensg, quote=F, row.names=F, col.names=F)

