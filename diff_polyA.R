args <- commandArgs(TRUE)
suppressMessages(library(DEXSeq))
suppressMessages(library(DESeq2))
library("BiocParallel")
setwd(".")

countFiles = c( "_out1", "_out2", "_out3", "_out4", "_out5", "_out6" )
sampleTable = data.frame(
row.names = c( "a1", "a2", "a3",
"b1", "b2", "b3" ),
condition = c("a", "a", "a",
"b", "b", "b"),
libType = c( "single-end", "single-end", "single-end",
"single-end", "single-end", "single-end" ) )

BPPARAM = MulticoreParam(workers=7)
dxd = DEXSeqDataSetFromHTSeq(countFiles, sampleData=sampleTable, design= ~ sample + exon + condition:exon )
dxd <- estimateSizeFactors( dxd )
dxd <- estimateDispersions( dxd, BPPARAM=BPPARAM)
dxd <- testForDEU( dxd, BPPARAM=BPPARAM )
dxd = estimateExonFoldChanges( dxd,  BPPARAM=BPPARAM)
res = DEXSeqResults( dxd )
write.csv(res, file = paste(args[1], ".csv", sep=""),  row.names=F, quote=FALSE)

