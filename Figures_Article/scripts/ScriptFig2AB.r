### Script for figures 2AB
suppressMessages(library(tidyverse))
theme_set(theme_bw())
library(patchwork)
library(ggdist)
mypal <- c(paletteer::paletteer_d("ggthemes::Classic_20"),"grey40")
`%+%` <- paste0
path_figures <- "./Figures_Article/figures/"
pathdata <- "./Figures_Article/data/"


## Figure2A
toplot <- read_tsv(paste0(pathdata,"Figure2A_data.tsv.gz"))
toplot$name4 <- factor(toplot$name4,levels=c(paste0("BT1_run",1:21),"BT1_run21b"))
totext <- toplot %>% group_by(name4) %>% summarise(n=n()) %>% ungroup
tomea <- toplot %>% group_by(name4) %>% summarise(mea=round(mean(speed))) %>% ungroup

speedmea0 <- tomea %>% pull(mea) %>% mean
# 2128

p1 <- ggplot(toplot,aes(x=name4,y=speed))+
	geom_hline(aes(yintercept=speedmea0),lty=2)+
	stat_slab(aes(fill=nanotype),col=NA,alpha=0.6,scale=0.8)+
	stat_pointinterval(.width=c(.5,.95),col="grey40")+
	geom_point(data=tomea,aes(x=name4,y=mea),col="red",shape=95,size=10)+
	geom_text(data=totext,aes(x=name4,y=0,label=n),fontface="italic") +
	geom_text(data=tomea,aes(x=name4,y=3800,label=mea),col="red") +
	scale_fill_discrete("30°C",type=mypal[c(3,4,7)])+
	scale_colour_discrete("30°C",type=mypal[c(3,4,7)])+
	coord_cartesian(ylim=c(0,4000))+
	ylab("Speed (bp/min)")+
	labs(tag="a")+
	theme(axis.text.x = element_text(angle = 45,hjust=1),axis.title.x=element_blank(),plot.tag=element_text(face="bold"))

## Figure 2B
toplot2 <- read_tsv(paste0(pathdata,"Figure2B_data.tsv.gz")) %>%
	mutate(name4=factor(name4))

totext2 <- toplot2 %>% group_by(name4) %>% summarise(n=n()) %>% ungroup
tomea2 <- toplot2 %>% group_by(name4) %>% summarise(mea=round(mean(speed))) %>% ungroup

p2 <- ggplot(toplot2,aes(x=name4,y=speed))+
	stat_slab(aes(fill=nanotype),col=NA,alpha=0.6,scale=0.8)+
	stat_pointinterval(.width=c(.5,.95),col="grey40")+
	geom_point(data=tomea2,aes(x=name4,y=mea),col="red",shape=95,size=10)+
	geom_text(data=totext2,aes(x=name4,y=0,label=n),fontface="italic") +
	geom_text(data=tomea2,aes(x=name4,y=3800,label=mea),col="red") +
	scale_fill_discrete("25°C",type=mypal[c(4)])+
	scale_colour_discrete("25°C",type=mypal[c(4)])+
	coord_cartesian(ylim=c(0,4000))+
	ylab("Speed (bp/min)")+
	labs(tag="b")+
	theme(plot.tag=element_text(face="bold"),axis.text.x = element_text(angle = 45,hjust=1),axis.title.x=element_blank())

p0 <- (p1|p2) + plot_layout(guides = "keep", ncol=2,byrow=F,widths = c(23,3))
ggsave(paste0(path_figures,"Figure2AB.pdf"),h=5,w=18,p0)

