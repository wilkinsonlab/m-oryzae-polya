source("/home/marco/Downloads/gProfileR/R/gProfileR.R")
args <- commandArgs(TRUE)
<<<<<<< HEAD
query<-read.table(args[1])
back<-read.table(args[2])
res<-gprofiler(query, organism="moryzae", ordered_query=args[3], custom_bg=c(back))
write.table(res, file=paste(file=args[1], "_go_enrich.tsv", sep=""), sep="\t")
=======
t<-read.table(args[1])
r<-gprofiler(t, organism="moryzae", ordered_query=args[2])
write.csv(r, file=paste(file=args[1], "_go_enrich.csv", sep=""), sep="\t")
>>>>>>> aba5408e32a4b6d8c1600d920cdd5e933769fe6f
