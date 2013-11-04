args <- commandArgs(TRUE)
suppressMessages(library("RColorBrewer"))
suppressMessages(library("gplots"))
suppressMessages(library(DESeq))
setwd(".")
a <- args[1]
b <- args[2]
c1<-read.table(paste(a, "-1.expr", sep=""),row.names=1) 
c2<-read.table(paste(a, "-2.expr", sep=""),row.names=1)
c3<-read.table(paste(a, "-3.expr", sep=""),row.names=1)
c4<-read.table(paste(b, "-1.expr", sep=""),row.names=1)
c5<-read.table(paste(b, "-2.expr", sep=""),row.names=1)
c6<-read.table(paste(b, "-3.expr", sep=""),row.names=1)
counts<-cbind(c1,c2,c3,c4,c5,c6)
colnames(counts)<-c(paste(a, "-1.expr", sep=""),paste(a, "-2.expr", sep=""),paste(a, "-3.expr", sep=""),paste(b, "-1.expr", sep=""),paste(b, "-2.expr", sep=""),paste(b, "-3.expr", sep=""))
design <- rep (c(a,b),each=3)
cds  <-  newCountDataSet(counts, design)
cds  <-  estimateSizeFactors( cds)
cds <- estimateDispersions( cds )
res  <-  nbinomTest(  cds,  a,  b)
resSig <- res[ res$padj < 0.05, ]


#summary = read.table("gene_summary.txt", sep="\t", header=TRUE)
#res = cbind(res, summary$Description)
#names(res)[names(res)=="summary$Description"] <- "desc"
write.csv(res, file = paste("diff_expr/", a, "_vs_", b, "_expr.csv", sep=""),  row.names=F, quote=FALSE)
write.table(counts, file = paste("diff_expr/", a, "_vs_", b, "_expr.count", sep=""),  sep="\t", col.names=F, quote=FALSE)

#plotDispEsts( cds )
#plotMA(res,col = ifelse(res$padj>=0.05, "black", "red3"))
png(file = paste("images/", a, "_vs_", b, "_pval.png", sep=""), width=600, height=600)
hist(res$pval, breaks=100, col="skyblue", border="slateblue", main="")
dev.off()
png(file = paste("images/", a, "_vs_", b, "_padj.png", sep=""), width=600, height=600)
hist(res$padj, breaks=100, col="skyblue", border="slateblue", main="")
dev.off()


cdsFullBlind = estimateDispersions( cds, method = "blind" )
vsdFull = varianceStabilizingTransformation( cdsFullBlind )

# heatmap samples
dists = dist( t( exprs(vsdFull) ) )
png(file = paste("images/", a, "_vs_", b, "_samples_heatmap.png", sep=""), width=600, height=600)
heatmap( as.matrix( dists ), symm=TRUE, cexCol=0.7, cexRow=0.7 )
dev.off()

# heatmap 30 most expressed genes
select = order(rowMeans(counts(cds)), decreasing=TRUE)[1:30]
hmcol = colorRampPalette(brewer.pal(9, "GnBu"))(100)
png(file = paste("images/", a, "_vs_", b, "_top30_heatmap.png", sep=""), width=600, height=600)
heatmap.2(exprs(vsdFull)[select,], col = hmcol, trace="none", margin=c(10, 6))
dev.off()

