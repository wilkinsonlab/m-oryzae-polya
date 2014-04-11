args <- commandArgs(TRUE)
suppressMessages(library("RColorBrewer"))
suppressMessages(library("gplots"))
suppressMessages(library(DESeq2))
directory<-"."
a <- args[1]
b <- args[2]
c1<-paste(a, "-1.expr", sep="")
c2<-paste(a, "-2.expr", sep="")
c3<-paste(a, "-3.expr", sep="")
c4<-paste(b, "-1.expr", sep="")
c5<-paste(b, "-2.expr", sep="")
c6<-paste(b, "-3.expr", sep="")
sampleFiles<-c(c1,c2,c3,c4,c5,c6)
sampleCondition<-c("a","a","a","b","b","b")
sampleTable<-data.frame(sampleName=sampleFiles, fileName=sampleFiles, condition=sampleCondition)
ddsHTSeq<-DESeqDataSetFromHTSeqCount(sampleTable=sampleTable, directory=directory, design=~condition)
colData(ddsHTSeq)$condition<-factor(colData(ddsHTSeq)$condition, levels=c("a","b"))
dds<-DESeq(ddsHTSeq)
res<-results(dds)
res<-data.frame(rownames(res),res)
colnames(res)[1] <- "gene"
write.csv(res, file = paste("diff_expr/", a, "_vs_", b, "_expr.csv", sep=""), row.names=F, quote=F)

png(file = paste("images/", a, "_vs_", b, "_MA.png", sep=""), width=600, height=600)
mar.default <- c(5,4,4,2) + 0.1
par(mar = mar.default + c(0, 1, 0, 0))
plotMA(dds,ylim=c(-10,10),main=paste(a, " -> ", b, sep=""), cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
dev.off()

png(file = paste("images/", a, "_vs_", b, "_DE.png", sep=""), width=600, height=600)
plotDispEsts(dds)
dev.off()

rld <- rlogTransformation(dds, blind=TRUE)
vsd <- varianceStabilizingTransformation(dds, blind=TRUE)


# heatmap samples
distsRL <- dist(t(assay(rld)))
mat <- as.matrix(distsRL)
png(file = paste("images/", a, "_vs_", b, "_samples_heatmap.png", sep=""), width=600, height=600)
heatmap.2(mat, trace="none", margin=c(13, 13))
dev.off()

# heatmap 30 most expressed genes
select = order(rowMeans(counts(dds)), decreasing=TRUE)[1:30]
png(file = paste("images/", a, "_vs_", b, "_top30_heatmap.png", sep=""), width=600, height=600)
heatmap.2(assay(vsd)[select,],  trace="none", margin=c(10, 6))
dev.off()

