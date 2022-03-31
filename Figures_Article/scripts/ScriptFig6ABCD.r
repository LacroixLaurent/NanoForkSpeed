### Script Fig6ABCD
suppressMessages(library(tidyverse))
library(patchwork)
library(ggdist)

theme_set(theme_bw())
mypal <- c(paletteer::paletteer_d("ggthemes::Classic_20"),"grey40")
`%+%` <- paste0
path_figures <- "./Figures_Article/figures/"
pathdata <- "./Figures_Article/data/"
pval <- read_tsv(paste0(pathdata,"pval_table.tsv.gz")) %>% mutate(type=case_when(type=="rRNA"~"rDNA",T~type))

### Figure6A speed feature
toplot <- read_tsv(paste0(pathdata,"Figure6A_data.tsv.gz"))
toplot$feat <- factor(toplot$feat,levels=c("centromere","telomere","rDNA","tRNA","other"))

toplot <- toplot %>% mutate(feat=fct_recode(feat,"Centromeres"="centromere","Telomeres"="telomere","rDNA"="rDNA","tRNA genes"="tRNA","Rest of the\n genome"="other"))
toplot2 <- toplot %>% group_by(feat) %>% summarise(m=round(mean(mea)),sd=sd(mea,na.rm=T),n=n()) %>% ungroup
topval <- left_join(toplot2 %>% select(feat),pval %>% filter(Figure=="6A") %>% select(feat=type,signif,pval_adj,signif1,signif2) %>% mutate(feat=fct_recode(feat,"Centromeres"="centromere","Telomeres"="telomere","rDNA"="rDNA","tRNA genes"="tRNA")))


f6a <- ggplot(toplot,aes(x=feat,y=mea))+
	geom_point(data=toplot2,aes(x=feat,y=m),col="red",shape=95,size=10,show.legend=F)+
	geom_errorbar(data=toplot2,aes(x=feat,y=m,ymin=m-sd,ymax=m+sd),col="red",width=0.2)+
	stat_dots(col="black",shape=16,side="both",dotsize=1.5,alpha=0.4)+
	geom_text(data=toplot2,aes(x=feat,y=0,label=n),col="black",fontface="italic",size=3) +
	geom_text(data=toplot2,aes(x=feat,y=3500,label=m),col="red",size=3) +
	geom_text(data=topval,aes(x=feat,y=3200,label=signif1),size=5) +
	geom_text(data=topval,aes(x=feat,y=3250,label=signif2),size=3) +
	scale_colour_manual("Features",values=mypal[c(3,5,7,9,1)])+
	coord_cartesian(ylim=c(0,3500))+
	theme(axis.text.x = element_text(angle=45,hjust=1),legend.position = "none",plot.tag=element_text(face="bold"))+
	xlab("Genomic features")+
	ylab("Speed (bp/min)")+
	labs(tag="a")


### Figure6b
toplot <- read_tsv(paste0(pathdata,"Figure6B_data.tsv.gz")) %>%
	mutate(type=fct_recode(type,"Co-directional"="CD","Head-on"="HO"))
toplot2 <- toplot %>% group_by(type) %>% summarise(m=round(mean(mea)),sd=sd(mea,na.rm=T),n=n()) %>% ungroup
f6b <- ggplot(toplot,aes(x=type,y=mea))+
	geom_point(data=toplot2,aes(x=type,y=m),col="red",shape=95,size=10,show.legend=F)+
	geom_errorbar(data=toplot2,aes(x=type,y=m,ymin=m-sd,ymax=m+sd),col="red",width=0.2)+
	stat_dots(col="black",shape=16,side="both",dotsize=1.8,alpha=0.4)+
	geom_text(data=toplot2,aes(x=type,y=3500,label=m),col="red",size=3) +
	annotate(geom="text",x=1.5,y=3200,label="*",col="black",size=5) +
	geom_text(data=toplot2,aes(x=type,y=0,label=n),col="black",fontface="italic",size=3) +
	scale_colour_manual("tRNA",values=mypal[c(1,7)],labels=c("Co-Directional","Head-On"))+
	coord_cartesian(ylim=c(0,3500))+
	theme(axis.text.x = element_text(angle=45,hjust=1),legend.position = "none",plot.tag=element_text(face="bold"))+
	xlab("tRNA genes")+
	ylab("Speed (bp/min)")+
	labs(tag="b")

### Figure6c lead/lag
toplot <- read_tsv(paste0(pathdata,"Figure6C_data.tsv.gz"))
toplot$type <- factor(toplot$type, levels=c("leading","lagging")) %>% fct_recode(.,"Leading"="leading","Lagging"="lagging")
toplot2 <- toplot %>% group_by(type) %>% summarise(m=round(mean(mea)),sd=sd(mea,na.rm=T),n=n()) %>% ungroup

f6c <- ggplot(toplot,aes(x=type,y=mea))+
	geom_point(data=toplot2,aes(x=type,y=m),col="red",shape=95,size=10,show.legend=F)+
	geom_errorbar(data=toplot2,aes(x=type,y=m,ymin=m-sd,ymax=m+sd),col="red",width=0.2)+
	stat_dots(col="black",shape=16,side="both",dotsize=2.7,alpha=0.4)+
	geom_text(data=toplot2,aes(x=type,y=3500,label=m),col="red",size=3) +
	geom_text(data=toplot2,aes(x=type,y=0,label=n),col="black",fontface="italic",size=3) +
	annotate(geom="text",x=1.5,y=3200,label="*",col="black",size=5) +
	scale_colour_manual("",values=mypal[c(1,7)])+
	coord_cartesian(ylim=c(0,3500))+
	theme(axis.text.x = element_text(angle=45,hjust=1),legend.position = "none",plot.tag=element_text(face="bold"))+
	xlab("Strand")+
	ylab("Speed (bp/min)")+
	labs(tag="c")


### Figure6d speed timing
toplot2 <- read_tsv(paste0(pathdata,"Figure6D_data.tsv.gz"))
library(ggpubr)
# rescale from 0 to 1
toplot2 %>% mutate(RT=(RTsc3-max(RTsc3,na.rm=T))/(min(RTsc3,na.rm=T)-max(RTsc3,na.rm=T))) -> toplot2

f6d <- ggplot(toplot2,aes(x=RT,y=speed))+
	geom_hex()+
	stat_cor(data=toplot2,aes(x=RT,y=speed),label.y = 5000,method="spearman",cor.coef.name="rho",col="grey20",fontface="italic",size=4)+
	scale_fill_continuous(type = "viridis",direction=-1,option="magma")+
#	scale_x_reverse()+
	ylim(c(0,5000))+
	ylab("Speed (bp/min)")+
	xlab("Replication timing (Early to Late)")+
	labs(tag="d")

layout <- "
AAAAABBCC
DDDDDDDDD"

p0 <- f6a+f6b+f6c+f6d+plot_layout(design = layout)
ggsave(paste0(path_figures,"Figure6.pdf"),h=9,w=7,p0)

