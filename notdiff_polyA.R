args <- commandArgs(TRUE)
suppressMessages(library(DESeq2))
directory<-"."
c1<-"_out1" 
c2<-"_out2"
c3<-"_out3"
c4<-"_out4"
c5<-"_out5"
c6<-"_out6"
sampleFiles<-c(c1,c2,c3,c4,c5,c6)
sampleCondition<-c("a","a","a","b","b","b")
sampleTable<-data.frame(sampleName=sampleFiles, fileName=sampleFiles, condition=sampleCondition)
ddsHTSeq<-DESeqDataSetFromHTSeqCount(sampleTable=sampleTable, directory=directory, design=~condition)
colData(ddsHTSeq)$condition<-factor(colData(ddsHTSeq)$condition, levels=c("a","b"))
dds<-DESeq(ddsHTSeq, fitType="local")
res<-results(dds)
res<-data.frame(rownames(res),res)
colnames(res)[1] <- "gene"
write.csv(res, file = args[1],  row.names=F, quote=FALSE)

