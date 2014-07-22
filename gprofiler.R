source("/home/marco/Downloads/gProfileR/R/gProfileR.R")
args <- commandArgs(TRUE)
species<-args[1]
query<-scan(args[2], what="character")
back<-scan(args[3], what="character")
order<-args[4]
res<-gprofiler(query, organism=species, custom_bg=back, ordered_query=order)

write.table(res, file=paste(file=args[2], "_go_enrich.tsv", sep=""), sep="\t")
