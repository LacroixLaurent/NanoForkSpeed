## Script Fig4A-B
suppressMessages(library(tidyverse))
library(patchwork)
library(ggdist)

theme_set(theme_bw())
mypal <- c(paletteer::paletteer_d("ggthemes::Classic_20"),"grey40")
`%+%` <- paste0
path_figures <- "/Users/ll/work/RStudioProjects/NanoForkSpeed/Figures_Article/figures/"
pathdata <- "/Users/ll/work/RStudioProjects/NanoForkSpeed/Figures_Article/data/"
pval <- read_tsv(paste0(pathdata,"pval_table.tsv.gz"))

### Figure4A HU
# selecting forks
toplot <- read_tsv(paste0(pathdata,"Figure4A_data.tsv.gz")) %>% mutate(HU=factor(as.character(HU),levels=c(0,1,2.5,5,10,25,50,100)))
toplot2 <- toplot %>% group_by(HU) %>% summarise(m=round(mean(mea)),sd=sd(mea,na.rm=T),n=n()) %>% ungroup
topval <- left_join(tibble(HU=factor(c(0,1,2.5,5,10,25,50,100))),pval %>% filter(Figure=="4A") %>% select(HU=type,signif,pval_adj))
f4a <- ggplot(toplot,aes(x=HU,y=mea))+
	stat_dots(col="grey40",alpha=0.6,size=1,shape=16,side="both")+
	geom_point(data=toplot2,aes(x=HU,y=m,col=HU),shape=95,alpha=0.6,size=10,show.legend=F)+
	geom_errorbar(data=toplot2,aes(x=HU,y=m,ymin=m-sd,ymax=m+sd,col=HU),width=0.2,alpha=0.6)+
	geom_text(data=toplot2,aes(x=HU,y=0,label=n),fontface="italic") +
	geom_text(data=toplot2,aes(x=HU,y=3400,label=m),col="red") +
	geom_text(data=topval,aes(x=HU,y=3100,label=signif),col="grey20",size=6) +
	scale_colour_discrete("HU (mM)",type=mypal)+
	coord_cartesian(ylim=c(0,3500))+
	theme(axis.ticks.x = element_blank(),axis.text.x = element_text(colour=mypal),legend.position = "none",plot.tag=element_text(face="bold"))+
	xlab("HU (mM)")+
	ylab("Speed (bp/min)")+
	labs(tag="a")


### Figure4B mutants
# selecting forks
toplot <- read_tsv(paste0(pathdata,"Figure4B_data.tsv.gz")) %>% mutate(mutant=factor(mutant,levels=c("WT","rtt109","sml1","csm3","tof1","mrc1")))
toplot2 <- toplot %>% group_by(mutant) %>% summarise(m=round(mean(mea)),sd=sd(mea,na.rm=T),n=n()) %>% ungroup
topval <- left_join(tibble(mutant=c("WT","rtt109","sml1","csm3","tof1","mrc1")),pval %>% filter(Figure=="4B") %>% select(mutant=type,signif,pval_adj))
f4b <- ggplot(toplot,aes(x=mutant,y=mea))+
	stat_dots(col="grey40",alpha=0.6,size=1,shape=16,side="both")+
	geom_point(data=toplot2,aes(x=mutant,y=m,col=mutant),shape=95,alpha=0.6,size=10,show.legend=F)+
	geom_errorbar(data=toplot2,aes(x=mutant,y=m,ymin=m-sd,ymax=m+sd,col=mutant),width=0.2,alpha=0.6)+
	geom_text(data=toplot2,aes(x=mutant,y=0,label=n),fontface="italic") +
	geom_text(data=toplot2,aes(x=mutant,y=3400,label=m),col="red") +
	geom_text(data=topval,aes(x=mutant,y=3100,label=signif),col="grey20",size=6) +
	coord_cartesian(ylim=c(0,3500))+
	scale_colour_manual("Strain",breaks=c("WT","rtt109","sml1","csm3","tof1","mrc1"),label=expression("BT1",paste("BT1 rtt109",Delta),paste("BT1 sml1",Delta),paste("BT1 csm3",Delta),paste("BT1 tof1",Delta),paste("BT1 mrc1",Delta)),values=mypal[c(1,3:7)])+
	guides(col=guide_legend(label.hjust=0))+
	scale_x_discrete(labels=expression("BT1",paste("BT1 rtt109",Delta),paste("BT1 sml1",Delta),paste("BT1 csm3",Delta),paste("BT1 tof1",Delta),paste("BT1 mrc1",Delta)))+
	theme(axis.ticks.x = element_blank(),axis.text.x = element_text(angle=45,hjust=1,colour=mypal[c(1,3:7)]),legend.position = "none",plot.tag=element_text(face="bold"))+
	xlab("Strain")+
	ylab("Speed (bp/min)")+
	labs(tag="b")

p0 <- (f4a|f4b) + plot_layout( ncol=2, byrow=F,widths=c(9,7))
ggsave(paste0(path_figures,"Figure4AB.pdf"),h=4,w=12,p0)
