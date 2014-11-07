args <- commandArgs(TRUE)
setwd(".")
file <- args[1]
m=read.csv(file, sep="\t")
m<-apply(m, 2, function(x) {cx <- cumsum(x); N <- length(x); (cx[5:N] - c(0, cx[1:(N-5)]))/5} )
write.table(file=paste(file, ".average", sep=""), m, sep="\t", row.names=F, col.names=F)

