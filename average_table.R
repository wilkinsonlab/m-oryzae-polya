args <- commandArgs(TRUE)
setwd(".")
file <- args[1]
m=read.csv(file, sep="\t", header=F,row.names=1)
m<-apply(m, 2, function(x) {cx <- cumsum(x); N <- length(x); (cx[20:N] - c(0, cx[1:(N-20)]))/20} )
write.table(file=paste(file, ".average", sep=""), m, sep="\t", row.names=T, col.names=F, quote=F)

