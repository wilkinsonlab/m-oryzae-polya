args <- commandArgs(TRUE)
setwd(".")
a <- as.numeric(args[1])
b <- as.numeric(args[2])
c <- as.numeric(args[3])
d <- as.numeric(args[4])
res <- paste( a, b, c, d, fisher.test(rbind(c(a,b),c(c,d)), alternative="g")$p.value, sep='\t')
cat(res)
