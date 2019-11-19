library(phytools)
library(RColorBrewer)
library(dplyr)

color_pal = sample(brewer.pal(9,"Set1"))

#color_pal = brewer.pal(9,"Set1")[c(2:9,1)]
#color_pal = c("red","blue","pink","purple")

## Test Case ##
#tree = "/Users/kevinamses/Documents/Flux_Downlaods/dissertation/rhop_neo/scgid/RHalf-1f_scgid_output/rscu/RHalf-1f_rscu_nj.tre"
#annot = "/Users/kevinamses/Documents/Flux_Downlaods/dissertation/rhop_neo/scgid/RHalf-1f_scgid_output/rscu/RHalf-1f_rscuTree_annot.csv"

tree = "/home/aimzez/work/jason/NRRL_3617/scgid_spades/r1_scgid_output/rscu/r1_rscu_nj.tre"
annot = "/home/aimzez/work/jason/NRRL_3617/scgid_spades/r1_scgid_output/rscu/r1_rscuTree_annot.csv"

args<-commandArgs(trailingOnly = TRUE)
tree = args[1]
annot = args[2]
output = args[3]

tree = read.tree(tree)
colnames=c("contig","taxlvl","tar_ntar","selection")
annot = read.delim(annot, sep=',', col.names = colnames, header = FALSE)
annot = annot %>% mutate(taxlvl = as.character(taxlvl))

edgecol = rep("black",tree$Nnode+length(tree$tip.label))
p=1
#taxon_groups = as.character(unique(annot$taxlvl))
#annot %>% group_by(taxlvl) %>% summarise(n=n())

#only important taxon groups
imp<-annot %>% group_by(taxlvl) %>% summarise(n=n()) %>% 
  filter(n>length(annot$taxlvl)/100)
if (dim(imp)[1]+1 > length(color_pal)) {
  bigger_pal = colorRampPalette(color_pal)
  color_pal = bigger_pal(dim(imp)[1]+1)
}
try(annot[!annot$taxlvl %in% imp$taxlvl,]$taxlvl<-"Other")

taxon_groups<-annot %>% group_by(taxlvl) %>% summarise(n=n())
if (!is.na(match('unclassified', taxon_groups$taxlvl))) {
  taxon_groups = taxon_groups[-match('unclassified', taxon_groups$taxlvl),]
}

for (t in taxon_groups$taxlvl) {
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
pdf(output,width = 20, height=20)
#tiff(output)
pt = plot(tree, type="fan", no.margin=TRUE, show.tip.label=TRUE, edge.width = 1.0,
          cex=0.3, edge.color=edgecol)
legend('topright', c(as.character(taxon_groups$taxlvl),'unclassified'), 
       col = c(color_pal[1:length(taxon_groups$taxlvl)],'black'), cex=1, pch=20)
pt = plot(tree, type="fan",no.margin=TRUE, show.tip.label=FALSE, edge.width = 1.0,
          cex=0.3, edge.color=edgecol)
legend('topright', c(as.character(taxon_groups$taxlvl),'unclassified'), 
       col = c(color_pal[1:length(taxon_groups$taxlvl)],'black'), cex=1, pch=20)

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
legend('topright', c(as.character(taxon_groups$taxlvl),'unclassified','trainset'), 
       col = c(color_pal[1:length(taxon_groups$taxlvl)],'black',color_pal[p]), cex=1, pch=20)
pt = plot(tree, type="fan",no.margin=FALSE, show.tip.label=FALSE, edge.width = 1.0,
          cex=0.3, edge.color=edgecol)
legend('topright', c(as.character(taxon_groups$taxlvl),'unclassified','trainset'), 
       col = c(color_pal[1:length(taxon_groups$taxlvl)],'black',color_pal[p]), cex=1, pch=20)
dev.off()
#nodelabels(cex=1.0,frame = "circle")
  