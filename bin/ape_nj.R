library(ape)
args<-commandArgs(trailingOnly = TRUE)
m<-read.table(args[1],sep=",",row.names=1)
tree <- nj(as.dist(m))
out_name <- args[2]
write.tree(tree,out_name)