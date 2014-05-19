args <- commandArgs(TRUE)
suppressMessages(library("RColorBrewer"))
suppressMessages(library("gplots"))
suppressMessages(library(DESeq2))
suppressMessages(library( "genefilter" ))

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
plotMA(dds,alpha=.05,ylim=c(-6,6),main=paste(a, " -> ", b, sep=""))
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
topVarGenes <- head( order( rowVars( assay(rld) ), decreasing=TRUE ), 30 )
png(file = paste("images/", a, "_vs_", b, "_top30_heatmap.png", sep=""), width=1024, height=768)
assay_table<-assay(rld)[ topVarGenes, ]
summary_table <- read.table("gene_summary.txt", sep="\t", quote="\"", row.names=1)
merge_table<-merge(assay_table,summary_table,by=0)
merge_table$Row.names<-paste(merge_table$Row.names,merge_table$V2)
merge_table$V2=NULL
rownames(merge_table) <- merge_table[,1]
merge_table$Row.names<-NULL
rownames(merge_table)<-sub(" \\[Source.*", "", rownames(merge_table))
heatmap.2( as.matrix(merge_table), scale="row",trace="none", dendrogram="column",col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255),margin=c(10,25))
dev.off()


