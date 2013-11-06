args <- commandArgs(TRUE)
res  = read.table(args[1])  
png(file = paste("images/", args[1], "_usage_ratio.png", sep=""), width=600, height=600)
plot(res, log="xy")
dev.off()

