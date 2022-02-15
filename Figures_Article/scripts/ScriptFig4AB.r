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
toplot <- read_tsv(paste0(pathdata,"Figure4A_data.tsv.gz")) %>% mutate(HU=factor(as.character(HU),levels=c(0,1,2.5,5,10,25,50,100)))
toplot2 <- toplot %>% group_by(HU) %>% summarise(m=round(mean(mea)),sd=sd(mea,na.rm=T),n=n()) %>% ungroup
topval <- left_join(tibble(HU=factor(c(0,1,2.5,5,10,25,50,100))),pval %>% filter(Figure=="4A") %>% select(HU=type,signif,pval_adj))
f4a <- ggplot(toplot,aes(x=HU,y=mea))+
	geom_point(data=toplot2,aes(x=HU,y=m),col="red",shape=95,size=10,show.legend=F)+
	geom_errorbar(data=toplot2,aes(x=HU,y=m,ymin=m-sd,ymax=m+sd),col="red",width=0.2)+
	stat_dots(col="black",shape=16,side="both",dotsize=1.5,alpha=0.4)+
	geom_text(data=toplot2,aes(x=HU,y=0,label=n),col="black",fontface="italic",size=3) +
	geom_text(data=toplot2,aes(x=HU,y=3500,label=m),col="red",size=3) +
	geom_text(data=topval,aes(x=HU,y=3200,label=signif),size=5) +
	scale_colour_discrete("HU (mM)",type=mypal)+
	coord_cartesian(ylim=c(0,3500))+
	theme(legend.position = "none",plot.tag=element_text(face="bold"))+
	xlab("HU (mM)")+
	ylab("Speed (bp/min)")+
	labs(tag="a")

### Figure4B mutants
toplot <- read_tsv(paste0(pathdata,"Figure4B_data.tsv.gz")) %>% mutate(mutant=factor(mutant,levels=c("WT","rtt109","sml1","csm3","tof1","mrc1")))
toplot2 <- toplot %>% group_by(mutant) %>% summarise(m=round(mean(mea)),sd=sd(mea,na.rm=T),n=n()) %>% ungroup
topval <- left_join(tibble(mutant=c("WT","rtt109","sml1","csm3","tof1","mrc1")),pval %>% filter(Figure=="4B") %>% select(mutant=type,signif,pval_adj))
f4b <- ggplot(toplot,aes(x=mutant,y=mea))+
	geom_point(data=toplot2,aes(x=mutant,y=m),col="red",shape=95,size=10,show.legend=F)+
	geom_errorbar(data=toplot2,aes(x=mutant,y=m,ymin=m-sd,ymax=m+sd),col="red",width=0.2)+
	stat_dots(col="black",shape=16,side="both",dotsize=1.5,alpha=0.4)+
	geom_text(data=toplot2,aes(x=mutant,y=0,label=n),col="black",fontface="italic",size=3) +
	geom_text(data=toplot2,aes(x=mutant,y=3500,label=m),col="red",size=3) +
	geom_text(data=topval,aes(x=mutant,y=3200,label=signif),size=5) +
	coord_cartesian(ylim=c(0,3500))+
	scale_colour_manual("Strain",breaks=c("WT","rtt109","sml1","csm3","tof1","mrc1"),label=expression("BT1",paste("BT1 rtt109",Delta),paste("BT1 sml1",Delta),paste("BT1 csm3",Delta),paste("BT1 tof1",Delta),paste("BT1 mrc1",Delta)),values=mypal[c(1,3:7)])+
	guides(col=guide_legend(label.hjust=0))+
	scale_x_discrete(labels=expression("BT1",paste("BT1 rtt109",Delta),paste("BT1 sml1",Delta),paste("BT1 csm3",Delta),paste("BT1 tof1",Delta),paste("BT1 mrc1",Delta)))+
	theme(axis.text.x = element_text(angle=45,hjust=1),legend.position = "none",plot.tag=element_text(face="bold"))+
	xlab("Strain")+
	ylab("Speed (bp/min )")+
	labs(tag="b")


p0 <- (f4a|f4b) + plot_layout( ncol=2, byrow=F,widths=c(9,7))
ggsave(paste0(path_figures,"Figure4AB.pdf"),h=6,w=7,p0)
