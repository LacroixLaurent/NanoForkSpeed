### script fig S10
### 20220202

suppressMessages(library(tidyverse))
library(patchwork)
library(ggdist)

theme_set(theme_bw())
mypal <- c(paletteer::paletteer_d("ggthemes::Classic_20"),"grey40")
`%+%` <- paste0
setwd("/Users/ll/work/Ori/NFS_paper/")
#path_figures <- "/Users/ll/work/RStudioProjects/NanoForkSpeed/Figures_Article/figures/"
#pathdata <- "/Users/ll/work/RStudioProjects/NanoForkSpeed/Figures_Article/data/"
pathdata <- "/Users/ll/work/Ori/NFS_paper/GitHub_upload/data/"
path_figures <- "/Users/ll/work/Ori/NFS_paper/GitHub_upload/figures/"
pval <- read_tsv(paste0(pathdata,"pval_table.tsv.gz"))


### CEN

toplot <- read_tsv(paste0(pathdata,"FigureS10A_data.tsv.gz"))
toplot$feat_name <- factor(toplot$feat_name,levels=c(paste0("CEN",1:16),"other"))
toplot2 <- toplot %>% group_by(feat_name) %>% summarise(n=n(),m=round(mean(mea)),sd=sd(mea,na.rm=T)) %>% ungroup
topval <- left_join(toplot2 %>% select(feat_name),pval %>% filter(Figure=="S10A") %>% select(feat_name=type,signif,pval_adj))

fs10a <- ggplot(toplot,aes(x=feat_name,y=mea))+
	stat_dots(col="grey40",alpha=0.6,size=1,shape=16,side="both")+
	geom_point(data=toplot2,aes(x=feat_name,y=m,col=feat_name),shape=95,alpha=0.6,size=10,show.legend=F)+
	geom_errorbar(data=toplot2,aes(x=feat_name,y=m,ymin=m-sd,ymax=m+sd,col=feat_name),width=0.2,alpha=0.6)+
	geom_text(data=toplot2,aes(x=feat_name,y=3900,label=m),col="red") +
	geom_text(data=topval,aes(x=feat_name,y=3600,label=signif),col="grey20",size=6) +
	geom_text(data=toplot2,aes(x=feat_name,y=0,label=n),fontface="italic") +
	coord_cartesian(ylim=c(0,4000))+
	scale_colour_manual("",values=mypal[c(3:18,1)])+
	ylab("Speed (bp/min)")+
	xlab("Centromeres")+
	theme(axis.ticks.x = element_blank(),axis.text.x = element_text(angle=45,hjust=1,colour=mypal[c(3:18,1)]),legend.position = "none",plot.tag=element_text(face="bold"),legend.key.size = unit(0.3, 'cm'))+
	labs(tag="a")


### TEL
toplot <- read_tsv(paste0(pathdata,"FigureS10B_data.tsv.gz"))
toplot$feat_name <- factor(toplot$feat_name,levels=unique(toplot$feat_name)[c(1:15,17,16)])
toplot2 <- toplot %>% group_by(feat_name) %>% summarise(n=n(),m=round(mean(mea)),sd=sd(mea,na.rm=T)) %>% ungroup
topval <- left_join(toplot2 %>% select(feat_name),pval %>% filter(Figure=="S10B") %>% select(feat_name=type,signif,pval_adj))

fs10b <- ggplot(toplot,aes(x=feat_name,y=mea))+
	stat_dots(col="grey40",alpha=0.6,size=1,shape=16,side="both")+
	geom_point(data=toplot2,aes(x=feat_name,y=m,col=feat_name),shape=95,alpha=0.6,size=10,show.legend=F)+
	geom_errorbar(data=toplot2,aes(x=feat_name,y=m,ymin=m-sd,ymax=m+sd,col=feat_name),width=0.2,alpha=0.6)+
	geom_text(data=toplot2,aes(x=feat_name,y=3900,label=m),col="red") +
	geom_text(data=topval,aes(x=feat_name,y=3600,label=signif),col="grey20",size=6) +
	geom_text(data=toplot2,aes(x=feat_name,y=0,label=n),fontface="italic") +
	coord_cartesian(ylim=c(0,4000))+
	scale_colour_manual("",values=mypal[c(3:18,1)])+
	ylab("Speed (bp/min)")+
	xlab("Telomeres")+
	theme(axis.ticks.x = element_blank(),axis.text.x = element_text(angle=45,hjust=1,colour=mypal[c(3:18,1)]),legend.position = "none",plot.tag=element_text(face="bold"),legend.key.size = unit(0.3, 'cm'))+
	labs(tag="b")

layout <- "
A
B
"

fs10 <- wrap_plots(fs10a,fs10b,design=layout)
ggsave(paste0(path_figures,"FigureS10.pdf"),h=8,w=9,fs10)

