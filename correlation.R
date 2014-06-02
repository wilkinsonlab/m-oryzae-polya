args <- commandArgs(TRUE)
setwd(".")
a <- args[1]
c1<-read.table(paste(a, "-1.expr", sep=""),row.names=1) 
c2<-read.table(paste(a, "-2.expr", sep=""),row.names=1)
c3<-read.table(paste(a, "-3.expr", sep=""),row.names=1)

mar.default <- c(5,4,4,2) + 0.1
png(file = paste("images/", a, "_1_vs_2_expr_cor.png", sep=""), width=600, height=600)
par(mar = mar.default + c(0, 1, 0, 0))
plot(log(c1$V2, 2), log(c2$V2, 2), main=a, xlab="log2 counts in replicate 1", ylab="log2 counts in replicate 2", cex.lab=1.5, cex.main=1.5)
text(13, 1, paste("Spearman cor ", round(cor(c1,c2,method="spearman"), 3)) ,cex = 1.5)
dev.off()
png(file = paste("images/", a, "_1_vs_3_expr_cor.png", sep=""), width=600, height=600)
par(mar = mar.default + c(0, 1, 0, 0))
plot(log(c1$V2, 2), log(c3$V2, 2), main=a, xlab="log2 counts in replicate 1", ylab="log2 counts in replicate 3", cex.lab=1.5, cex.main=1.5)
text(13, 1, paste("Spearman cor ", round(cor(c1,c3,method="spearman"), 3)) ,cex = 1.5)
dev.off()
png(file = paste("images/", a, "_2_vs_3_expr_cor.png", sep=""), width=600, height=600)
par(mar = mar.default + c(0, 1, 0, 0))
plot(log(c2$V2, 2), log(c3$V2, 2), main=a, xlab="log2 counts in replicate 2", ylab="log2 counts in replicate 3", cex.lab=1.5, cex.main=1.5)
text(13, 1, paste("Spearman cor ", round(cor(c2,c3,method="spearman"), 3)) ,cex = 1.5)
dev.off()

cat("1 vs 2,", paste(cor(c1,c2,method="spearman"), "\n"))
cat("1 vs 3,", paste(cor(c1,c3,method="spearman"), "\n"))
cat("2 vs 3,", paste(cor(c2,c3,method="spearman"), "\n"))


p1<-read.table(paste(a, "-1.polyA_all_m", sep="")) 
p2<-read.table(paste(a, "-2.polyA_all_m", sep=""))
p3<-read.table(paste(a, "-3.polyA_all_m", sep=""))
p1<-data.frame(paste(p1$V5, p1$V2, sep=":"), p1$V1)
p2<-data.frame(paste(p2$V5, p2$V2, sep=":"), p2$V1)
p3<-data.frame(paste(p3$V5, p3$V2, sep=":"), p3$V1)
colnames(p1) <- c("polyA", "val")
colnames(p2) <- c("polyA", "val")
colnames(p3) <- c("polyA", "val")
m1<-merge(p1, p2,by = "polyA", all.x=TRUE, all.y=TRUE)
m2<-merge(p1, p3,by = "polyA", all.x=TRUE, all.y=TRUE)
m3<-merge(p2, p3,by = "polyA", all.x=TRUE, all.y=TRUE)
m1[is.na(m1)] <- 0
m2[is.na(m2)] <- 0
m3[is.na(m3)] <- 0

p1_a<-read.table(paste(a, "-1.polyA", sep=""))
p2_a<-read.table(paste(a, "-2.polyA", sep=""))
p3_a<-read.table(paste(a, "-3.polyA", sep=""))
p1_a<-data.frame(paste(p1_a$V5, p1_a$V2, sep=":"), p1_a$V1)
p2_a<-data.frame(paste(p2_a$V5, p2_a$V2, sep=":"), p2_a$V1)
p3_a<-data.frame(paste(p3_a$V5, p3_a$V2, sep=":"), p3_a$V1)
colnames(p1_a) <- c("polyA", "val")
colnames(p2_a) <- c("polyA", "val")
colnames(p3_a) <- c("polyA", "val")
m1_a<-merge(p1_a, p2_a,by = "polyA", all.x=TRUE, all.y=TRUE)
m2_a<-merge(p1_a, p3_a,by = "polyA", all.x=TRUE, all.y=TRUE)
m3_a<-merge(p2_a, p3_a,by = "polyA", all.x=TRUE, all.y=TRUE)
m1_a[is.na(m1_a)] <- 0
m2_a[is.na(m2_a)] <- 0
m3_a[is.na(m3_a)] <- 0


png(file = paste("images/", a, "_1_vs_2_polyA_cor.png", sep=""), width=600, height=600)
par(mar = mar.default + c(0, 1, 0, 0))
plot(log(m1_a$val.x, 2), log(m1_a$val.y, 2), main=a, xlab="log2 counts in replicate 1", ylab="log2 counts in replicate 2", cex.main=1.5, cex.lab=1.5,  col="black")
points(log(m1$val.x, 2), log(m1$val.y, 2), col="red")
#text(10, 3, paste("Spearman cor ", round(cor(m1$val.x, m1$val.y,method="spearman"), 3)) ,cex = 1.2)
dev.off()
png(file = paste("images/", a, "_1_vs_3_polyA_cor.png", sep=""), width=600, height=600)
par(mar = mar.default + c(0, 1, 0, 0))
plot(log(m2_a$val.x, 2), log(m2_a$val.y, 2), main=a, xlab="log2 counts in replicate 1", ylab="log2 counts in replicate 3", cex.main=1.5, cex.lab=1.5, col="black")
points(log(m2$val.x, 2), log(m2$val.y, 2), col="red")
#text(10, 3, paste("Spearman cor ", round(cor(m2$val.x, m2$val.y,method="spearman"), 3)) ,cex = 1.2)
dev.off()
png(file = paste("images/", a, "_2_vs_3_polyA_cor.png", sep=""), width=600, height=600)
par(mar = mar.default + c(0, 1, 0, 0))	
plot(log(m3_a$val.x, 2), log(m3_a$val.y, 2), main=a, xlab="log2 counts in replicate 2", ylab="log2 counts in replicate 3", cex.main=1.5, cex.lab=1.5, col="black")
points(log(m3$val.x, 2), log(m3$val.y, 2), col="red")
#text(10, 3, paste("Spearman cor ", round(cor(m3$val.x, m3$val.y,method="spearman"), 3)) ,cex = 1.2)
dev.off()

cat("1 vs 2,", paste(cor(m1_a$val.x, m1_a$val.y,method="spearman"), ",", cor(m1$val.x, m1$val.y,method="spearman"), "\n"))
cat("1 vs 3,", paste(cor(m2_a$val.x, m2_a$val.y,method="spearman"), ",", cor(m2$val.x, m2$val.y,method="spearman"), "\n"))
cat("2 vs 3,", paste(cor(m3_a$val.x, m3_a$val.y,method="spearman"), ",", cor(m3$val.x, m3$val.y,method="spearman"), "\n"))


