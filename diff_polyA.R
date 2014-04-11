args <- commandArgs(TRUE)
suppressMessages(library(DEXSeq))
suppressMessages(library(DESeq2))
library("parallel")
setwd(".")
sampleTable <- data.frame(
 row.names = c( "a1", "a2", "a3", "b1", "b2", "b3" ),
 countFile = c( "_out1", "_out2", "_out3", "_out4", "_out5", "_out6" ),
 condition = c( "a","a","a","b","b","b" ),
libType = c( "single-end", "single-end", "single-end", "single-end", "single-end", "single-end" ) )
ecs <- read.HTSeqCounts(as.character(sampleTable$countFile), sampleTable)
ecs <- estimateSizeFactors( ecs )
ecs <- estimateDispersions( ecs,minCount=1,nCores=7)
ecs <- fitDispersionFunction( ecs )
ecs <- testForDEU( ecs,nCores=7 )
ecs <- estimatelog2FoldChanges( ecs,nCores=7 )
res <- DEUresultTable(ecs)
write.csv(res, file = paste(args[1], ".csv", sep=""),  row.names=F, quote=FALSE)


png(file = paste("images/", args[1], "_MA.png", sep=""), width=600, height=600)
mar.default <- c(5,4,4,2) + 0.1
par(mar = mar.default + c(0, 1, 0, 0))
res_ <- res[,5:7]
res_[,1] <- res[,6]
res_[,2] <- res[,7]
res_[,3] <- res[,5]
res_[,3] <- res_[,3] < 0.05
res_[,3][is.na(res_)[,3]] <- FALSE
plotMA(res_, ylim=c(-4,4),main=paste(args[1], sep=""), cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
dev.off()
