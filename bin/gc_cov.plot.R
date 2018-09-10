
args<-commandArgs(trailingOnly = TRUE)

tigInfo<-read.table(args[1],sep="\t")
unclass<-read.table(args[2],sep="\t")
colnames(tigInfo)<-c("contig","prot_len","coverage","gc","pid","sp_os","lineage","evalue","parse_lineage")
colnames(unclass)<-c("contig","coverage","gc","taxonomy")
gc_window<-args[3]
cov_window<-args[4]
output<-args[5]

gc_window = gsub("[()]","",gc_window)
cov_window = gsub("[()]","",cov_window)

gc_window = strsplit(gc_window,split=',')
cov_window = strsplit(cov_window,split=',')

choices = c(rgb(0,0,1,alpha=0.6),rgb(0.54,0,0,alpha=0.6))
tigInfo$color = "black"
try(tigInfo[tigInfo$evalue == 0.0,]$evalue<-2.225074e-308, silent = TRUE)
tigInfo$log_evalue<-log(tigInfo$evalue)
tigInfo[tigInfo$parse_lineage == "target",]$color<-choices[1]
tigInfo[tigInfo$parse_lineage == "nontarget",]$color<-choices[2]
unclass$color<-rgb(0,0,0,alpha=0.45)
rect_color = rgb(0,1,0,alpha=0.5)
xlim = c(0,1050)


leg<-matrix(nrow=4,ncol=3)
leg[,1]<-c("Target","Non-target","Unclassified","Window")
leg[,2]<-c(choices,rgb(0,0,0,alpha=0.45),"black")
leg[,3]<-c(19,19,19,22)

pdf(output)
plot(tigInfo$coverage,tigInfo$gc,cex=0.25,col=tigInfo$color,xlab="Coverage",ylab="GC content",pch=19,xlim=xlim)
rect(as.numeric(cov_window[[1]][1]),as.numeric(gc_window[[1]][1]),
     as.numeric(cov_window[[1]][2]),as.numeric(gc_window[[1]][2]),col=rect_color)
legend("topright",leg[,1],col=leg[,2],pch=c(19,19,19,22),pt.bg=rect_color,pt.cex=2,cex=0.9)

plot(tigInfo$coverage,tigInfo$gc,cex=0.25,col=tigInfo$color,xlab="Coverage",ylab="GC content",pch=19,xlim=xlim)
points(unclass$coverage,unclass$gc,cex=0.25,col=unclass$color,pch=19)
rect(as.numeric(cov_window[[1]][1]),as.numeric(gc_window[[1]][1]),
     as.numeric(cov_window[[1]][2]),as.numeric(gc_window[[1]][2]),col=rect_color)
legend("topright",leg[,1],col=leg[,2],pch=c(19,19,19,22),pt.bg=rect_color,pt.cex=2,cex=0.9)

plot(tigInfo$coverage,tigInfo$gc,cex=0.25,col=tigInfo$color,xlab="Coverage",ylab="GC content",pch=19,xlim=xlim)
legend("topright",leg[,1],col=leg[,2],pch=c(19,19,19,22),pt.bg=rect_color,pt.cex=2,cex=0.9)

plot(unclass$coverage,unclass$gc,cex=0.25,col=unclass$color,pch=19,xlab="Coverage",ylab="GC content",xlim=xlim)
legend("topright",leg[,1],col=leg[,2],pch=c(19,19,19,22),pt.bg=rect_color,pt.cex=2,cex=0.9)

dev.off()
