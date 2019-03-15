library(dplyr)
library(reshape2)
library(ggplot2)

getTaxLevel <- function (row) {
  spl<-as.list(strsplit(row, ", ")[[1]])
  return(spl[2])
}

### Test case
#args = c("~/Documents/Flux_Downlaods/dissertation/rhop_neo/scgid/RHalf-1f_scgid_output/blob/RHalf-1f_info_table.tsv",
#         "~/Documents/Flux_Downlaods/dissertation/rhop_neo/scgid/RHalf-1f_scgid_output/blob/RHalf-1f_unclassified_info_table.tsv",
#         "(0.41,0.75","(0.0468227,8.0)")
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


imp<-tigInfo %>% group_by(hitlin_level) %>% summarise (n = n()) %>% 
  filter(n>dim(tigInfo)[1]/1000)

try(tigInfo[!tigInfo$hitlin_level %in% imp$hitlin_level,]$hitlin_level <- "Other", silent = TRUE)
try(tigInfo[tigInfo$evalue == 0,]$evalue<-2.225074e-308, silent = TRUE)
tigInfo<-tigInfo %>% mutate(log_evalue = log(evalue))


plt<-tigInfo %>% select(contigid,coverage,gc,parse_lin,log_evalue,hitlin_level)
plt$coverage<-log(plt$coverage)
#plt<-melt(plt, id.vars=c("contigid","coverage","gc"))
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
tntPalette <- c("darkred","blue")

unclass$color = "black"
unclass$coverage = log(unclass$coverage)

pdf(output)

### Just annotated points
ggplot(plt, aes(coverage, gc, color=hitlin_level, size=log_evalue/-150)) + 
  geom_point(alpha=0.20) +
  xlab("log(Coverage)")+
  ylab("GC content") +
  #scale_colour_manual(values=cbPalette) +
  scale_color_brewer(type="qual", palette = "Dark2") +
  guides(
    color = guide_legend(override.aes = list(alpha=1), title="Hit Taxonomy"), 
    size = guide_legend(override.aes=list(alpha=1), title="Scaled Hit Significance")
    ) +
  #geom_point(data=unclass, aes(log(coverage),gc, size=0.1), alpha=0.10, inherit.aes = FALSE) +
  guides (
    alpha = FALSE
  ) +
  #geom_rect(data=window, aes(xmin=cov_window[1], ymin=gc_window[1], 
                             #xmax=cov_window[2], ymax=gc_window[2]), fill = "green", color="black", alpha=0.2, inherit.aes = FALSE) +
  theme_bw()

## annotated points and unclassified
ggplot(plt, aes(coverage, gc, color=hitlin_level, size=log_evalue/-150)) + 
  geom_point(alpha=0.20) +
  xlab("log(Coverage)")+
  ylab("GC content") +
  scale_colour_manual(values=cbPalette) +
  guides(
    color = guide_legend(override.aes = list(alpha=1), title="Hit Taxonomy"), 
    size = guide_legend(override.aes=list(alpha=1), title="Scaled Hit Significance")
  ) +
  geom_point(data=unclass, aes(coverage,gc, size=0.1, color=parse_lin), alpha=1.0, inherit.aes = FALSE) +
  guides (
    alpha = FALSE
  ) +
  #geom_rect(data=window, aes(xmin=cov_window[1], ymin=gc_window[1], 
                             #xmax=cov_window[2], ymax=gc_window[2]), fill = "green", color="black", alpha=0.2, inherit.aes = FALSE) +
  theme_bw()

## annotated points with window
ggplot(plt, aes(coverage, gc, color=hitlin_level, size=log_evalue/-150)) + 
  geom_point(alpha=0.20) +
  xlab("log(Coverage)")+
  ylab("GC content") +
  scale_colour_manual(values=cbPalette) +
  guides(
    color = guide_legend(override.aes = list(alpha=1), title="Hit Taxonomy"), 
    size = guide_legend(override.aes=list(alpha=1), title="Scaled Hit Significance")
  ) +
  geom_point(data=unclass, aes(coverage,gc, size=0.1), alpha=1.0, inherit.aes = FALSE) +
  guides (
    alpha = FALSE
  ) +
  geom_rect(data=window, aes(xmin=cov_window[1], ymin=gc_window[1], 
                             xmax=cov_window[2], ymax=gc_window[2]), fill = "green", color="black", alpha=0.2, inherit.aes = FALSE) +
  theme_bw()

## annotated points and unlclass with window and color_changing_points
unclass[((unclass$coverage>=cov_window[1] & unclass$coverage<=cov_window[2]) & 
          (unclass$gc>=gc_window[1] & unclass$gc<=gc_window[2])),]$color<-"blue"
ggplot(plt, aes(coverage, gc, color=hitlin_level, size=log_evalue/-150)) + 
  geom_point(alpha=0.20) +
  xlab("log(Coverage)")+
  ylab("GC content") +
  scale_colour_manual(values=cbPalette) +
  guides(
    color = guide_legend(override.aes = list(alpha=1), title="Hit Taxonomy"), 
    size = guide_legend(override.aes=list(alpha=1), title="Scaled Hit Significance")
  ) +
  geom_point(data=unclass, aes(coverage, gc, size=0.1), color=unclass$color, alpha=1.0, inherit.aes = FALSE) +
  guides (
    alpha = FALSE
  ) +
  geom_rect(data=window, aes(xmin=cov_window[1], ymin=gc_window[1], 
                             xmax=cov_window[2], ymax=gc_window[2]), fill = "green", alpha=0.2, inherit.aes = FALSE) +
  #xlim(min(unclass$coverage)*1.10, max(unclass$coverage)*1.10) +
  #ylim(min(unclass$gc)*0.9, max(unclass$gc)*1.10) +
  theme_bw()

## Specified Target/Nontarget instead of taxonomy
plt_tnt = plt
ggplot(plt_tnt, aes(coverage, gc, color=parse_lin, size=log_evalue/-150)) + 
  geom_point(alpha=0.20) +
  xlab("log(Coverage)")+
  ylab("GC content") +
  scale_colour_manual(values=tntPalette) +
  guides(
    color = guide_legend(override.aes = list(alpha=1), title="Hit Taxonomy"), 
    size = guide_legend(override.aes=list(alpha=1), title="Scaled Hit Significance")
  ) +
  #geom_point(data=unclass, aes(log(coverage),gc, size=0.1), alpha=0.10, inherit.aes = FALSE) +
  guides (
    alpha = FALSE
  ) +
  #geom_rect(data=window, aes(xmin=cov_window[1], ymin=gc_window[1], 
  #xmax=cov_window[2], ymax=gc_window[2]), fill = "green", color="black", alpha=0.2, inherit.aes = FALSE) +
  theme_bw()

## Specified Target/Nontarget with window
ggplot(plt_tnt, aes(coverage, gc, color=parse_lin, size=log_evalue/-150)) + 
  geom_point(alpha=0.20) +
  xlab("log(Coverage)")+
  ylab("GC content") +
  scale_colour_manual(values=tntPalette) +
  guides(
    color = guide_legend(override.aes = list(alpha=1), title="Hit Taxonomy"), 
    size = guide_legend(override.aes=list(alpha=1), title="Scaled Hit Significance")
  ) +
  #geom_point(data=unclass, aes(log(coverage),gc, size=0.1), alpha=0.10, inherit.aes = FALSE) +
  guides (
    alpha = FALSE
  ) +
  geom_rect(data=window, aes(xmin=cov_window[1], ymin=gc_window[1], 
    xmax=cov_window[2], ymax=gc_window[2]), fill = "green", color="black", alpha=0.2, inherit.aes = FALSE) +
  theme_bw()
dev.off()
 