args <- commandArgs(TRUE)
t <- read.table("_t")
t$V3 = (2**t$V1*t$V2*2)/(2**t$V1+1)
t$V4 = (t$V2 * 2) - t$V3
t$V5 = t$V3 - t$V4
png(file = paste("images/", args[1], "_fold.png", sep=""), width=600, height=600)
plot(t$V5,t$V1,ylim=c(-4,4),xlim=c(-50,50),cex=0.8, xlab="expression difference", ylab="log2foldChange")
#plot(t$V2,t$V1,log="x",ylim=c(-4,4),cex=0.8, xlab="baseMean", ylab="log2foldChange")
abline(h=0)
dev.off()
