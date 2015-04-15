args <- commandArgs(TRUE)

directory<-"."


a1 = read.table(file=args[1], row.names=1)
a2 = read.table(file=args[2], row.names=1)
a3 = read.table(file=args[3], row.names=1)
b1 = read.table(file=args[4], row.names=1)
b2 = read.table(file=args[5], row.names=1)
b3 = read.table(file=args[6], row.names=1)
a1 = a1 / t(replicate(nrow(a1), colSums(a1)))
a2 = a2 / t(replicate(nrow(a2), colSums(a2)))
a3 = a3 / t(replicate(nrow(a3), colSums(a3)))
b1 = b1 / t(replicate(nrow(b1), colSums(b1)))
b2 = b2 / t(replicate(nrow(b2), colSums(b2)))
b3 = b3 / t(replicate(nrow(b3), colSums(b3)))

WT_mean=rowMeans(cbind(a1, a2, a3))
mutant_mean=rowMeans(cbind(b1, b2, b3))
diff=mutant_mean-WT_mean
stdev=sd(diff)
zscore=diff/stdev
pval=2*pnorm(-abs(zscore))
padj=p.adjust(pval, method="BH")
write.csv(file=args[7],cbind(a1, a2, a3, b1, b2, b3, zscore, pval, padj))
