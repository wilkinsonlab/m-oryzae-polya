args <- commandArgs(TRUE)
setwd(".")
filename <- args[1]
table <- read.csv(filename, sep="\t")
q_values <- p.adjust(table$p_value, method="BH")
new <- cbind(table[1:6], q_values, table[7:8])
write.table(new, file=filename,sep="\t", quote=FALSE, row.names=FALSE)
