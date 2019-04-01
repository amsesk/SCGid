library(phytools)
library(RColorBrewer)
library(dplyr)

color_pal = sample(brewer.pal(9,"Set1"))

## Test Case ##
#tree = "/Users/kevinamses/Documents/Flux_Downlaods/dissertation/rhop_neo/scgid/RHalf-1f_scgid_output/rscu/RHalf-1f_rscu_nj.tre"
#annot = "/Users/kevinamses/Documents/Flux_Downlaods/dissertation/rhop_neo/scgid/RHalf-1f_scgid_output/rscu/RHalf-1f_rscuTree_annot.csv"

args<-commandArgs(trailingOnly = TRUE)
tree = args[1]
annot = args[2]
output = args[3]

tree = read.tree(tree)
colnames=c("contig","taxlvl","tar_ntar","selection")
annot = read.table(annot, sep=',', col.names = colnames)

edgecol = rep("black",tree$Nnode+length(tree$tip.label))
p=1
taxon_groups = as.character(unique(annot$taxlvl))
annot %>% group_by(taxlvl) %>% summarise(n=n())
taxon_groups = taxon_groups[-match('unclassified', taxon_groups)]
for (t in taxon_groups) {
  taxon_entire_clades = c()
  for (n in getDescendants(tree, length(tree$tip.label)+1)) {
    nd = getDescendants(tree, n)
    node_desc_tips = nd[nd<(length(tree$tip.label)+1)]
    if (all(annot[node_desc_tips,]$taxlvl == t)) {
      taxon_entire_clades = c(taxon_entire_clades, n)
    }
  }
  edgecol[tree$edge[,2] %in% taxon_entire_clades] = color_pal[p]
  p=p+1
}
pdf(output)
pt = plot(tree, type="fan",no.margin=TRUE, show.tip.label=TRUE, edge.width = 1.0,
          cex=0.3, edge.color=edgecol)
legend('topright', c(as.character(taxon_groups),'unclassified'), 
       col = c(color_pal[1:length(taxon_groups)],'black'), cex=1, pch=20)
pt = plot(tree, type="fan",no.margin=TRUE, show.tip.label=FALSE, edge.width = 1.0,
          cex=0.3, edge.color=edgecol)
legend('topright', c(as.character(taxon_groups),'unclassified'), 
       col = c(color_pal[1:length(taxon_groups)],'black'), cex=1, pch=20)

trainset_clades = c()
for (n in getDescendants(tree, length(tree$tip.label)+1)) {
  nd = getDescendants(tree,n)
  node_desc_tips = nd[nd<(length(tree$tip.label)+1)]
  if (all(annot[node_desc_tips,]$selection == "trainset")) {
    trainset_clades = c(trainset_clades,n)
  }
}
edgecol[tree$edge[,2] %in% trainset_clades] = color_pal[p]

pt = plot(tree, type="fan",no.margin=FALSE, show.tip.label=TRUE, edge.width = 1.0,
          cex=0.3, edge.color=edgecol)
legend('topright', c(as.character(taxon_groups),'unclassified','trainset'), 
       col = c(color_pal[1:length(taxon_groups)],'black',color_pal[p]), cex=1, pch=20)
pt = plot(tree, type="fan",no.margin=FALSE, show.tip.label=FALSE, edge.width = 1.0,
          cex=0.3, edge.color=edgecol)
legend('topright', c(as.character(taxon_groups),'unclassified','trainset'), 
       col = c(color_pal[1:length(taxon_groups)],'black',color_pal[p]), cex=1, pch=20)
dev.off()
#nodelabels(cex=1.0,frame = "circle")
  