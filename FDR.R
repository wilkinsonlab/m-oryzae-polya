args <- commandArgs(TRUE)
setwd(".")
filename <- args[1]
table <- read.csv(filename, sep="\t")
q_values <- p.adjust(table$p_value, method="BH")
table <- cbind(table, q_values)
table<-table[,c("term", "a", "b", "c", "d", "p_value", "q_values", "description", "genes")]
write.table(table, file=filename,sep="\t", quote=FALSE, row.names=FALSE)
