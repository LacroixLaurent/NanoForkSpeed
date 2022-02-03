# script FigS4
## 202120202
suppressMessages(library(tidyverse))
theme_set(theme_bw())
library(patchwork)
library(ggdist)

mypal <- c(paletteer::paletteer_d("ggthemes::Classic_20"),"grey40")
`%+%` <- paste0
setwd("/Users/ll/work/Ori/NFS_paper/")
#path_figures <- "/Users/ll/work/RStudioProjects/NanoForkSpeed/Figures_Article/figures/"
#pathdata <- "/Users/ll/work/RStudioProjects/NanoForkSpeed/Figures_Article/data/"
pathdata <- "/Users/ll/work/Ori/NFS_paper/GitHub_upload/data/"
path_figures <- "/Users/ll/work/Ori/NFS_paper/GitHub_upload/figures/"


toplot <- read_tsv(paste0(pathdata,"FigureS4_data.tsv.gz")) %>% mutate(B_pulse=factor(as.character(B_pulse),levels=unique(sort(B_pulse))))
totext <- toplot %>% group_by(B_pulse) %>% summarise(n=n()) %>% ungroup
tomea <- toplot %>% group_by(B_pulse) %>% summarise(mea=round(mean(speed))) %>% ungroup

p1 <- ggplot(toplot,aes(x=B_pulse,y=speed))+
	stat_slab(aes(fill=B_pulse),col=NA,alpha=0.6,scale=0.8,normalize="groups")+
	stat_pointinterval(.width=c(.5,.95),col="grey40")+
	geom_point(data=tomea,aes(x=B_pulse,y=mea),col="red",shape=95,size=10,show.legend=F)+
	coord_cartesian(ylim=c(0,5000))+
	geom_text(data=totext,aes(x=B_pulse,y=0,label=n),fontface="italic") +
	geom_text(data=tomea,aes(x=B_pulse,y=4800,label=mea),col="red") +
	scale_fill_manual("[BrdU] pulse",values=mypal)+
	ylab("Speed (bp/min)")+
	theme(axis.text.x = element_text(angle = 45,hjust=1),plot.tag=element_text(face="bold"))+
	labs(tag="b")

tomea2 <- toplot %>% group_by(B_pulse) %>% summarise(mea=round(mean(dY),3)) %>% ungroup

p2 <- ggplot(toplot,aes(x=B_pulse,y=dY))+
	stat_slab(aes(fill=B_pulse),col=NA,alpha=0.6,scale=0.8,normalize="groups")+
	stat_pointinterval(.width=c(.5,.95),col="grey40")+
	geom_point(data=tomea2,aes(x=B_pulse,y=mea),col="red",shape=95,size=10,show.legend=F)+
	coord_cartesian(ylim=c(0,1))+
	geom_text(data=totext,aes(x=B_pulse,y=0,label=n),fontface="italic") +
	geom_text(data=tomea2,aes(x=B_pulse,y=1,label=mea),col="red") +
	scale_fill_manual("[BrdU] pulse",values=mypal)+
	theme(axis.text.x = element_text(angle = 45,hjust=1),plot.tag=element_text(face="bold"))+
	labs(tag="a")+
	ylab("BrdU signal amplitude")
	
p0 <- p2 / p1 & theme(legend.position = "right",axis.text.x = element_blank(),axis.title.x = element_blank(),axis.ticks.x = element_blank())
p0 + plot_layout(guides = "collect")
ggsave(paste0(path_figures,"FigureS4.pdf"),h=7,w=6)

