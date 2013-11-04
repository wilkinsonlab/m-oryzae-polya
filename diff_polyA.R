args <- commandArgs(TRUE)
suppressMessages(library(DEXSeq))
setwd(".")
sampleTable <- data.frame(
 row.names = c( "a1", "a2", "a3", "b1", "b2", "b3" ),
 countFile = c( "_out1", "_out2", "_out3", "_out4", "_out5", "_out6" ),
 condition = c( "a","a","a","b","b","b" ),
libType = c( "single-end", "single-end", "single-end", "single-end", "single-end", "single-end" ) )
ecs <- read.HTSeqCounts(as.character(sampleTable$countFile), sampleTable)
ecs <- estimateSizeFactors( ecs )
ecs <- estimateDispersions( ecs,minCount=1 )
ecs <- fitDispersionFunction( ecs )
ecs <- testForDEU( ecs )
ecs <- estimatelog2FoldChanges( ecs )
res <- DEUresultTable(ecs)
write.csv(res, file = args[1],  row.names=F, quote=FALSE)
