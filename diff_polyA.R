args <- commandArgs(TRUE)
suppressMessages(library(DEXSeq))
suppressMessages(library(DESeq))
library("parallel")
setwd(".")
sampleTable <- data.frame(
 row.names = c( "a1", "a2", "a3", "b1", "b2", "b3" ),
 countFile = c( "_out1", "_out2", "_out3", "_out4", "_out5", "_out6" ),
 condition = c( "a","a","a","b","b","b" ),
libType = c( "single-end", "single-end", "single-end", "single-end", "single-end", "single-end" ) )
ecs <- read.HTSeqCounts(as.character(sampleTable$countFile), sampleTable)
ecs <- estimateSizeFactors( ecs )
ecs <- estimateDispersions( ecs,minCount=1,nCores=4)
ecs <- fitDispersionFunction( ecs )
ecs <- testForDEU( ecs,nCores=4 )
ecs <- estimatelog2FoldChanges( ecs,nCores=4 )
res <- DEUresultTable(ecs)
write.csv(res, file = paste(args[1], ".csv", sep=""),  row.names=F, quote=FALSE)
