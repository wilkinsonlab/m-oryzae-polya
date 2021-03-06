args <- commandArgs(TRUE)
suppressMessages(library(DESeq))
setwd(".")
counts<-read.table(args[1], row.names=1)
design <- rep (c("a","b"),each=3)
cds  <-  newCountDataSet(counts, design)
cds  <-  estimateSizeFactors( cds)
cds <- estimateDispersions( cds, fitType="local"  )
res  <-  nbinomTest(  cds,  "a",  "b")
write.csv(res, file = args[2],  row.names=F, quote=FALSE)
png(file = paste("images/", args[1], "_pval.png", sep=""), width=600, height=600)
hist(res$pval, breaks=100, col="skyblue", border="slateblue", main="")
dev.off()
png(file = paste("images/", args[1], "_padj.png", sep=""), width=600, height=600)
hist(res$padj, breaks=100, col="skyblue", border="slateblue", main="")
dev.off()

