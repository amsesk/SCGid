library(dplyr)
library(reshape2)
library(ggplot2)
library(RColorBrewer)

getTaxLevel <- function (row) {
  spl<-as.list(strsplit(row, ", ")[[1]])
  return(spl[2])
}

### Test case
#args = c("~/Documents/Flux_Downlaods/dissertation/rhop_neo/scgid/RHalf-1f_scgid_output/blob/RHalf-1f_info_table.tsv",
#         "~/Documents/Flux_Downlaods/dissertation/rhop_neo/scgid/RHalf-1f_scgid_output/blob/RHalf-1f_unclassified_info_table.tsv",
#         "(0.41,0.75","(0.0468227,8.0)")

### Pagey Test case
args = c("~/Documents/Flux_Downlaods/dissertation/scgid/stylopage_scgid_output/blob/stylopage_info_table.tsv",
         "~/Documents/Flux_Downlaods/dissertation/scgid/stylopage_scgid_output/blob/stylopage_unclassified_info_table.tsv",
         "(0.2559,0.4698)","(29.8954,1370.0834)")
### Read in data from <scgid gc-cov>

#tigInfo<-read.table(args[1],sep="\t",
#                    col.names = c("contigid","plen","coverage","gc","pid","hitsp","hitlin","evalue","parse_lin"))
#unclass<-read.table(args[2],sep="\t", col.names=c("contigid","coverage","gc","parse_lin"))
#unclass$color = "black"
#unclass$coverage = log(unclass$coverage)

#gc_window<-args[3]
#cov_window<-args[4]

args<-commandArgs(trailingOnly = TRUE)

tigInfo<-read.table(args[1],sep="\t")
unclass<-read.table(args[2],sep="\t")
colnames(tigInfo)<-c("contigid","plen","coverage","gc","pid","hitsp","hitlin","evalue","parse_lin")
colnames(unclass)<-c("contigid","coverage","gc","parse_lin")
gc_window<-args[3]
cov_window<-args[4]
output<-args[5]

lin<-as.character(tigInfo$hitlin)
lin<-sub("\\[","",lin)
lin<-sub("\\]","",lin)
linapp<-unlist(lapply(lin, getTaxLevel))
tigInfo$hitlin_level<-linapp

gc_window = gsub("[()]","",gc_window)
cov_window = gsub("[()]","",cov_window)
gc_window = as.numeric(unlist(strsplit(gc_window,split=',')))
cov_window = log(as.numeric(unlist(strsplit(cov_window,split=','))))
window = as.data.frame(cbind(cov_window, gc_window))
window$col<-"Best Window"


imp<-tigInfo %>% group_by(hitlin_level) %>% summarise (n = n()) %>% 
  filter(n>dim(tigInfo)[1]/100)

try(tigInfo[!tigInfo$hitlin_level %in% imp$hitlin_level,]$hitlin_level <- "Other", silent = TRUE)
try(tigInfo[tigInfo$evalue == 0,]$evalue<-2.225074e-308, silent = TRUE)
tigInfo<-tigInfo %>% mutate(log_evalue = log(evalue))


plt<-tigInfo %>% select(contigid,coverage,gc,parse_lin,log_evalue,hitlin_level)
plt$coverage<-log(plt$coverage)
#plt<-melt(plt, id.vars=c("contigid","coverage","gc"))
#cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
#tntPalette <- c("darkred","blue")

unclass$coverage = log(unclass$coverage)
unclass$included<-"Unclassified_Dumped"
unclass[((unclass$coverage>=cov_window[1] & unclass$coverage<=cov_window[2]) & 
           (unclass$gc>=gc_window[1] & unclass$gc<=gc_window[2])),]$included<-"Unclassified_Kept"

myDark2<-brewer.pal(dim(imp)[1]+1, "Dark2")
myDark2_blk<-c(myDark2, "black")
myDark2_blkblu<-c(myDark2_blk,"blue")

pdf(output)

my_theme_bw<-function(base_size=12, base_family="") {
    theme_bw (base_size=base_size, base_family=base_family) +
    theme (
      axis.title.y = element_text(size=34, margin = margin(t=0, r=10, b=0, l=0)),
      axis.title.x = element_text(size=34, margin = margin(t=10, r=0, b=0, l=0)),
      axis.text = element_text(size=30, color="black"),
      legend.text = element_text(size=20),
      legend.title = element_text(size=20),
      plot.margin = unit(c(1,1,1,1),"cm"),
      aspect.ratio = 1.0
    )
}
my_theme_bw_noY<-function(base_size=12, base_family="") {
  theme_bw (base_size=base_size, base_family=base_family) +
    theme (
      axis.title.y = element_blank(),
      axis.title.x = element_text(size=34, margin = margin(t=10, r=0, b=0, l=0)),
      axis.text.y = element_blank(),
      axis.text.x = element_text(size=30, color="black"),
      legend.text = element_text(size=20),
      legend.title = element_text(size=20),
      plot.margin = unit(c(1,1,1,1),"cm"),
      aspect.ratio = 1.0
    )
}
bw_bare<-function(base_size=12, base_family="") {
  theme_bw (base_size=base_size, base_family=base_family) +
    theme (
      axis.title.y = element_blank(),
      axis.title.x = element_blank(),
      axis.text.y = element_blank(),
      axis.text.x = element_blank(),
      legend.text = element_text(size=20),
      legend.title = element_text(size=20),
      #plot.margin = unit(c(1,1,1,1),"cm"),
      aspect.ratio = 1.0
    )
}

### Just annotated points
blob0d = ggplot(plt, aes(coverage, gc, color=hitlin_level, size=log_evalue/-150)) + 
  geom_point(alpha=0.20) +
  xlab("ln(Coverage)")+
  ylab("GC content") +
  scale_colour_manual(values=myDark2) +
  guides(
    color = FALSE,
    size = FALSE
    #color = guide_legend(override.aes = list(alpha=1, size=5), title="Hit Taxonomy"), 
    #size = guide_legend(override.aes=list(alpha=1), title="Scaled Hit Significance")
    ) +

  bw_bare()
  #my_theme_bw()
  #my_theme_bw_noY()

blob1d = ggplot(plt, aes(coverage, gc, color=hitlin_level, size=log_evalue/-150)) + 
  geom_point(alpha=0.20) +
  xlab("ln(Coverage)")+
  ylab("GC content") +
  scale_colour_manual(values=myDark2) +
  guides(
    color = FALSE,
    size = FALSE
    #color = guide_legend(override.aes = list(alpha=1, size=5), title="Hit Taxonomy"), 
    #size = guide_legend(override.aes=list(alpha=1), title="Scaled Hit Significance")
  ) +
  #scale_y_continuous( expand=c(0,0)) +
  #scale_x_continuous( expand=c(0,0)) +
  ### Window annotations
  geom_vline(xintercept=log(326.08817752423016), color="black",linetype=2) +
  
  geom_vline(xintercept = cov_window[1], color="black") +
  geom_vline(xintercept = cov_window[2], color="black") +
  annotate("rect", xmin=cov_window[1], ymin=0.255566, xmax=cov_window[2], 
           ymax=0.585654, fill="black", alpha=0.2) +
  bw_bare()

blob2d = ggplot(plt, aes(coverage, gc, color=hitlin_level, size=log_evalue/-150)) + 
  geom_point(alpha=0.20) +
  xlab("ln(Coverage)")+
  ylab("GC content") +
  scale_colour_manual(values=myDark2) +
  guides(
    color = FALSE,
    size = FALSE
    #color = guide_legend(override.aes = list(alpha=1, size=5), title="Hit Taxonomy"), 
    #size = guide_legend(override.aes=list(alpha=1), title="Scaled Hit Significance")
  ) +
  #scale_y_continuous( expand=c(0,0)) +
  #scale_x_continuous( expand=c(0,0)) +
  ### Window annotations
  geom_hline(yintercept=0.39517770523381646, color="black", linetype=4) +
  geom_vline(xintercept=log(326.08817752423016), color="black",linetype=2) +
  
  geom_hline(yintercept = gc_window[1], color="black") +
  geom_hline(yintercept = gc_window[2], color="black") +
  annotate("rect", xmin=0, ymin=gc_window[1], xmax=7.223048, 
                            ymax=gc_window[2], fill="black", alpha=0.20) +
  
  geom_vline(xintercept = cov_window[1], color="black") +
  geom_vline(xintercept = cov_window[2], color="black") +
  annotate("rect", xmin=cov_window[1], ymin=0.255566, xmax=cov_window[2], 
           ymax=0.585654, fill="black", alpha=0.2) +
  bw_bare()
#my_theme_bw()
#my_theme_bw_noY()

render=ggplot_build(blob)
render$layout$panel_params[[1]]$x.range
render$layout$panel_params[[1]]$y.range

## annotated points and unlclass with window and color_changing_points
done = ggplot(plt, aes(coverage, gc, color=hitlin_level, size=log_evalue/-150)) + 
  geom_point(alpha=0.20) +
  xlab("log(Coverage)")+
  ylab("GC content") +
  scale_colour_manual(values=myDark2_blkblu) +
  scale_fill_manual(values=c("green")) +
  xlim(0.2731,7.2230) +
  ylim(0.2556,0.5857) +
  guides(
    color = guide_legend(override.aes = list(alpha=1, size=5), title="Hit Taxonomy"), 
    #color = FALSE,
    size = guide_legend(override.aes=list(alpha=1), title="Scaled Hit Significance"),
    #size = FALSE,
    fill = guide_legend(title="Windows")
    #fill = FALSE
  ) +
  geom_hline(yintercept=0.39517770523381646, color="black", linetype=4) +
  geom_vline(xintercept=log(326.08817752423016), color="black",linetype=2) +
  
  geom_hline(yintercept = gc_window[1], color="black") +
  geom_hline(yintercept = gc_window[2], color="black") +
  annotate("rect", xmin=0.2731, ymin=gc_window[1], xmax=7.2230, 
           ymax=gc_window[2], fill="black", alpha=0.20) +
  
  geom_vline(xintercept = cov_window[1], color="black") +
  geom_vline(xintercept = cov_window[2], color="black") +
  annotate("rect", xmin=cov_window[1], ymin=0.2556, xmax=cov_window[2], 
           ymax=0.5857, fill="black", alpha=0.2) +
  
  geom_point(data=unclass, aes(coverage, gc, size=0.15, color=included), alpha=1.0, inherit.aes = FALSE) +
  geom_rect(data=window, aes(xmin=cov_window[1], ymin=gc_window[1], 
                             xmax=cov_window[2], ymax=gc_window[2], fill=col), color="black", alpha=0.2, inherit.aes = FALSE) +
  #bw_bare()
  my_theme_bw()

#### Line graphs
library(reshape2)
library(gridExtra)
library(ggpubr)
library(grid)

points = read.delim("~/Documents/Papers/scgid/figures/gct_window_calc/cov2_pts.csv", sep=',')
gc_points = read.delim("/Users/kevinamses/Documents/Papers/scgid/figures/gct_window_calc/gc2_pts.csv",sep=',')
colnames(gc_points)<-c("Nontarget","Target","Trade-off","step","direction")

points_melt = melt(points, id.vars=c('direction','step'))
gc_points_melt = melt(gc_points, id.vars=c('direction','step'))

stepwise_col = c("red","blue","purple")

cov_steps = ggplot(points_melt, aes(log(step),value, color=variable)) +
  geom_line(size=0.8) +
  guides(
    color = FALSE
  ) +
  geom_vline(xintercept=log(326.08817752423016), color="black",linetype=2) +
  geom_vline(xintercept = cov_window[1], color="black") +
  geom_vline(xintercept = cov_window[2], color="black") +
  annotate("rect", xmin=cov_window[1], ymin=0, xmax=cov_window[2], 
           ymax=1, fill="black", alpha=0.2) +
  scale_color_manual(values=stepwise_col) +
  xlim(0.2731,7.2230) +
  ylab("value") +
  xlab("ln(Coverage)") +
  bw_bare()

gc_steps = ggplot(gc_points_melt, aes(step, value, color=variable)) +
  geom_line(size=0.8) +
  guides(
    color = guide_legend(title="Window Statistics", size=5)
  ) +
  geom_vline(xintercept=0.39517770523381646, color="black", linetype=4) +
  geom_vline(xintercept = gc_window[1], color="black") +
  geom_vline(xintercept = gc_window[2], color="black") +
  
  scale_color_manual(values=stepwise_col) +
  annotate("rect", ymin=0, xmin=gc_window[1], ymax=1, 
           xmax=gc_window[2], fill="black", alpha=0.20) +
  xlim(0.2556,0.5857) +
  ylab("value") +
  xlab("GC") +
  my_theme_bw()
  #my_theme_bw()
gc_steps+coord_flip()


blank <- grid.rect(gp=gpar(col="white"))

grid.arrange(blob0d, blob1d, blob2d, gc_steps+coord_flip(), blank, cov_steps, done, blank, 
             nrow=2, ncol=4)
