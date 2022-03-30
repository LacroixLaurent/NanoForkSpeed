### script fig S12
suppressMessages(library(tidyverse))
library(patchwork)
library(ggdist)

theme_set(theme_bw())
mypal <- c(paletteer::paletteer_d("ggthemes::Classic_20"),"grey40")
`%+%` <- paste0
path_figures <- "/Users/ll/work/RStudioProjects/NanoForkSpeed/Figures_Article/figures/"
pathdata <- "/Users/ll/work/RStudioProjects/NanoForkSpeed/Figures_Article/data/"
pval <- read_tsv(paste0(pathdata,"pval_table.tsv.gz"))

speedm <- 2130

# tRNA
toplot <- read_tsv(paste0(pathdata,"FigureS12_data.tsv.gz"))
tRNA_tib3 <- toplot %>% group_by(feat_name,type) %>% summarise(m=round(mean(mea)),sd=sd(mea,na.rm=T),n=n()) %>% ungroup
tRNA_tib3$feat_name <- fct_reorder(tRNA_tib3$feat_name, tRNA_tib3$m)

tRNA2plot <- left_join(tRNA_tib3,pval %>% filter(Figure=="S11") %>% select(feat_name=type,signif,pval_adj,signif1,signif2)) %>% filter(n>1)
tRNA2plot$feat_name <- fct_reorder(tRNA2plot$feat_name, tRNA2plot$m)

tRNA2plot$qtile <- cut(tRNA2plot$m,quantile(tRNA2plot$m,probs=seq(0,1,1/12)),include.lowest=T,labels=c(paste0("q",1:12)))

toplot1 <- tRNA2plot %>% filter(qtile=="q1")
p1 <- ggplot(toplot1)+
	coord_cartesian(ylim=c(0,5000))+
	geom_hline(aes(yintercept=speedm),lty=2)+
	geom_point(aes(x=feat_name,y=m,col=type),shape=95,size=10,show.legend=F)+
	geom_errorbar(aes(x=feat_name,y=m,ymin=m-sd,ymax=m+sd,col=type),width=0.2)+
	geom_text(aes(x=feat_name,y=0,label=n),fontface="italic") +
	geom_text(aes(x=feat_name,y=4800,label=round(m)),col="red") +
	geom_text(aes(x=feat_name,y=4250,label=signif1),size=5) +
	geom_text(aes(x=feat_name,y=4400,label=signif2),size=3) +
	scale_colour_manual("",values=mypal[c(1,7)])+
	theme(axis.text.x = element_text(angle = 45,hjust=1),axis.title.x = element_blank(),axis.ticks.x = element_blank())+
	ylab("Speed (bp/min)")+
	ggtitle("tRNA_q1")

toplot2 <- tRNA2plot %>% filter(qtile=="q2")
p2 <- ggplot(toplot2)+
	coord_cartesian(ylim=c(0,5000))+
	geom_hline(aes(yintercept=speedm),lty=2)+
	geom_point(aes(x=feat_name,y=m,col=type),shape=95,size=10,show.legend=F)+
	geom_errorbar(aes(x=feat_name,y=m,ymin=m-sd,ymax=m+sd,col=type),width=0.2)+
	geom_text(aes(x=feat_name,y=0,label=n),fontface="italic") +
	geom_text(aes(x=feat_name,y=4800,label=round(m)),col="red") +
	geom_text(aes(x=feat_name,y=4250,label=signif1),size=5) +
	geom_text(aes(x=feat_name,y=4400,label=signif2),size=3) +
	scale_colour_manual("",values=mypal[c(1,7)])+
	theme(axis.text.x = element_text(angle = 45,hjust=1),axis.title.x = element_blank(),axis.ticks.x = element_blank())+
	ylab("Speed (bp/min)")+
	ggtitle("tRNA_q2")

toplot3 <- tRNA2plot %>% filter(qtile=="q3")
p3 <- ggplot(toplot3)+
	coord_cartesian(ylim=c(0,5000))+
	geom_hline(aes(yintercept=speedm),lty=2)+
	geom_point(aes(x=feat_name,y=m,col=type),shape=95,size=10,show.legend=F)+
	geom_errorbar(aes(x=feat_name,y=m,ymin=m-sd,ymax=m+sd,col=type),width=0.2)+
	geom_text(aes(x=feat_name,y=0,label=n),fontface="italic") +
	geom_text(aes(x=feat_name,y=4800,label=round(m)),col="red") +
	geom_text(aes(x=feat_name,y=4250,label=signif1),size=5) +
	geom_text(aes(x=feat_name,y=4400,label=signif2),size=3) +
	scale_colour_manual("",values=mypal[c(1,7)])+
	theme(axis.text.x = element_text(angle = 45,hjust=1),axis.title.x = element_blank(),axis.ticks.x = element_blank())+
	ylab("Speed (bp/min)")+
	ggtitle("tRNA_q3")

toplot4 <- tRNA2plot %>% filter(qtile=="q4")
p4 <- ggplot(toplot4)+
	coord_cartesian(ylim=c(0,5000))+
	geom_hline(aes(yintercept=speedm),lty=2)+
	geom_point(aes(x=feat_name,y=m,col=type),shape=95,size=10,show.legend=F)+
	geom_errorbar(aes(x=feat_name,y=m,ymin=m-sd,ymax=m+sd,col=type),width=0.2)+
	geom_text(aes(x=feat_name,y=0,label=n),fontface="italic") +
	geom_text(aes(x=feat_name,y=4800,label=round(m)),col="red") +
	geom_text(aes(x=feat_name,y=4250,label=signif1),size=5) +
	geom_text(aes(x=feat_name,y=4400,label=signif2),size=3) +
	scale_colour_manual("",values=mypal[c(1,7)])+
	theme(axis.text.x = element_text(angle = 45,hjust=1),axis.title.x = element_blank(),axis.ticks.x = element_blank())+
	ylab("Speed (bp/min)")+
	ggtitle("tRNA_q4")


p0 <- p1/p2/p3/p4 & theme(legend.position = "bottom",axis.title.x=element_blank())
p0 + plot_layout(guides = "collect")
ggsave(paste0(path_figures,"FigureS12A.pdf"),h=12,w=12)

toplot1 <- tRNA2plot %>% filter(qtile=="q5")
p1 <- ggplot(toplot1)+
	coord_cartesian(ylim=c(0,5000))+
	geom_hline(aes(yintercept=speedm),lty=2)+
	geom_point(aes(x=feat_name,y=m,col=type),shape=95,size=10,show.legend=F)+
	geom_errorbar(aes(x=feat_name,y=m,ymin=m-sd,ymax=m+sd,col=type),width=0.2)+
	geom_text(aes(x=feat_name,y=0,label=n),fontface="italic") +
	geom_text(aes(x=feat_name,y=4800,label=round(m)),col="red") +
	geom_text(aes(x=feat_name,y=4250,label=signif1),size=5) +
	geom_text(aes(x=feat_name,y=4400,label=signif2),size=3) +
	scale_colour_manual("",values=mypal[c(1,7)])+
	theme(axis.text.x = element_text(angle = 45,hjust=1),axis.title.x = element_blank(),axis.ticks.x = element_blank())+
	ylab("Speed (bp/min)")+
	ggtitle("tRNA_q5")

toplot2 <- tRNA2plot %>% filter(qtile=="q6")
p2 <- ggplot(toplot2)+
	coord_cartesian(ylim=c(0,5000))+
	geom_hline(aes(yintercept=speedm),lty=2)+
	geom_point(aes(x=feat_name,y=m,col=type),shape=95,size=10,show.legend=F)+
	geom_errorbar(aes(x=feat_name,y=m,ymin=m-sd,ymax=m+sd,col=type),width=0.2)+
	geom_text(aes(x=feat_name,y=0,label=n),fontface="italic") +
	geom_text(aes(x=feat_name,y=4800,label=round(m)),col="red") +
	geom_text(aes(x=feat_name,y=4250,label=signif1),size=5) +
	geom_text(aes(x=feat_name,y=4400,label=signif2),size=3) +
	scale_colour_manual("",values=mypal[c(1,7)])+
	theme(axis.text.x = element_text(angle = 45,hjust=1),axis.title.x = element_blank(),axis.ticks.x = element_blank())+
	ylab("Speed (bp/min)")+
	ggtitle("tRNA_q6")

toplot3 <- tRNA2plot %>% filter(qtile=="q7")
p3 <- ggplot(toplot3)+
	coord_cartesian(ylim=c(0,5000))+
	geom_hline(aes(yintercept=speedm),lty=2)+
	geom_point(aes(x=feat_name,y=m,col=type),shape=95,size=10,show.legend=F)+
	geom_errorbar(aes(x=feat_name,y=m,ymin=m-sd,ymax=m+sd,col=type),width=0.2)+
	geom_text(aes(x=feat_name,y=0,label=n),fontface="italic") +
	geom_text(aes(x=feat_name,y=4800,label=round(m)),col="red") +
	geom_text(aes(x=feat_name,y=4250,label=signif1),size=5) +
	geom_text(aes(x=feat_name,y=4400,label=signif2),size=3) +
	scale_colour_manual("",values=mypal[c(1,7)])+
	theme(axis.text.x = element_text(angle = 45,hjust=1),axis.title.x = element_blank(),axis.ticks.x = element_blank())+
	ylab("Speed (bp/min)")+
	ggtitle("tRNA_q7")

toplot4 <- tRNA2plot %>% filter(qtile=="q8")
p4 <- ggplot(toplot4)+
	coord_cartesian(ylim=c(0,5000))+
	geom_hline(aes(yintercept=speedm),lty=2)+
	geom_point(aes(x=feat_name,y=m,col=type),shape=95,size=10,show.legend=F)+
	geom_errorbar(aes(x=feat_name,y=m,ymin=m-sd,ymax=m+sd,col=type),width=0.2)+
	geom_text(aes(x=feat_name,y=0,label=n),fontface="italic") +
	geom_text(aes(x=feat_name,y=4800,label=round(m)),col="red") +
	geom_text(aes(x=feat_name,y=4250,label=signif1),size=5) +
	geom_text(aes(x=feat_name,y=4400,label=signif2),size=3) +
	scale_colour_manual("",values=mypal[c(1,7)])+
	theme(axis.text.x = element_text(angle = 45,hjust=1),axis.title.x = element_blank(),axis.ticks.x = element_blank())+
	ylab("Speed (bp/min)")+
	ggtitle("tRNA_q8")


p0 <- p1/p2/p3/p4 & theme(legend.position = "bottom",axis.title.x=element_blank())
p0 + plot_layout(guides = "collect")
ggsave(paste0(path_figures,"FigureS12B.pdf"),h=12,w=12)

toplot1 <- tRNA2plot %>% filter(qtile=="q9")
p1 <- ggplot(toplot1)+
	coord_cartesian(ylim=c(0,5000))+
	geom_hline(aes(yintercept=speedm),lty=2)+
	geom_point(aes(x=feat_name,y=m,col=type),shape=95,size=10,show.legend=F)+
	geom_errorbar(aes(x=feat_name,y=m,ymin=m-sd,ymax=m+sd,col=type),width=0.2)+
	geom_text(aes(x=feat_name,y=0,label=n),fontface="italic") +
	geom_text(aes(x=feat_name,y=4800,label=round(m)),col="red") +
	geom_text(aes(x=feat_name,y=4250,label=signif1),size=5) +
	geom_text(aes(x=feat_name,y=4400,label=signif2),size=3) +
	scale_colour_manual("",values=mypal[c(1,7)])+
	theme(axis.text.x = element_text(angle = 45,hjust=1),axis.title.x = element_blank(),axis.ticks.x = element_blank())+
	ylab("Speed (bp/min)")+
	ggtitle("tRNA_q9")

toplot2 <- tRNA2plot %>% filter(qtile=="q10")
p2 <- ggplot(toplot2)+
	coord_cartesian(ylim=c(0,5000))+
	geom_hline(aes(yintercept=speedm),lty=2)+
	geom_point(aes(x=feat_name,y=m,col=type),shape=95,size=10,show.legend=F)+
	geom_errorbar(aes(x=feat_name,y=m,ymin=m-sd,ymax=m+sd,col=type),width=0.2)+
	geom_text(aes(x=feat_name,y=0,label=n),fontface="italic") +
	geom_text(aes(x=feat_name,y=4800,label=round(m)),col="red") +
	geom_text(aes(x=feat_name,y=4250,label=signif1),size=5) +
	geom_text(aes(x=feat_name,y=4400,label=signif2),size=3) +
	scale_colour_manual("",values=mypal[c(1,7)])+
	theme(axis.text.x = element_text(angle = 45,hjust=1),axis.title.x = element_blank(),axis.ticks.x = element_blank())+
	ylab("Speed (bp/min)")+
	ggtitle("tRNA_q10")

toplot3 <- tRNA2plot %>% filter(qtile=="q11")
p3 <- ggplot(toplot3)+
	coord_cartesian(ylim=c(0,5000))+
	geom_hline(aes(yintercept=speedm),lty=2)+
	geom_point(aes(x=feat_name,y=m,col=type),shape=95,size=10,show.legend=F)+
	geom_errorbar(aes(x=feat_name,y=m,ymin=m-sd,ymax=m+sd,col=type),width=0.2)+
	geom_text(aes(x=feat_name,y=0,label=n),fontface="italic") +
	geom_text(aes(x=feat_name,y=4800,label=round(m)),col="red") +
	geom_text(aes(x=feat_name,y=4250,label=signif1),size=5) +
	geom_text(aes(x=feat_name,y=4400,label=signif2),size=3) +
	scale_colour_manual("",values=mypal[c(1,7)])+
	theme(axis.text.x = element_text(angle = 45,hjust=1),axis.title.x = element_blank(),axis.ticks.x = element_blank())+
	ylab("Speed (bp/min)")+
	ggtitle("tRNA_q11")

toplot4 <- tRNA2plot %>% filter(qtile=="q12")
p4 <- ggplot(toplot4)+
	coord_cartesian(ylim=c(0,5000))+
	geom_hline(aes(yintercept=speedm),lty=2)+
	geom_point(aes(x=feat_name,y=m,col=type),shape=95,size=10,show.legend=F)+
	geom_errorbar(aes(x=feat_name,y=m,ymin=m-sd,ymax=m+sd,col=type),width=0.2)+
	geom_text(aes(x=feat_name,y=0,label=n),fontface="italic") +
	geom_text(aes(x=feat_name,y=4800,label=round(m)),col="red") +
	geom_text(aes(x=feat_name,y=4250,label=signif1),size=5) +
	geom_text(aes(x=feat_name,y=4400,label=signif2),size=3) +
	scale_colour_manual("",values=mypal[c(1,7)])+
	theme(axis.text.x = element_text(angle = 45,hjust=1),axis.title.x = element_blank(),axis.ticks.x = element_blank())+
	ylab("Speed (bp/min)")+
	ggtitle("tRNA_q12")


p0 <- p1/p2/p3/p4 & theme(legend.position = "bottom",axis.title.x=element_blank())
p0 + plot_layout(guides = "collect")
ggsave(paste0(path_figures,"FigureS12C.pdf"),h=12,w=12)

