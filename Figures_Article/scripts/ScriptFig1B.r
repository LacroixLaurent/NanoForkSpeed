# script figure 1B
suppressMessages(library(tidyverse))
theme_set(theme_bw())
library(patchwork)
library(ggrastr)

mypal <- c(paletteer::paletteer_d("ggthemes::Classic_20"),"grey40")
`%+%` <- paste0
#setwd("/Users/ll/work/Ori/NFS_paper/")
path_figures <- "/Users/ll/work/RStudioProjects/NanoForkSpeed/Figures_Article/figures/"
pathdata <- "/Users/ll/work/RStudioProjects/NanoForkSpeed/Figures_Article/data/"
toplot <- readRDS(paste0(pathdata,"Figure1B_data.rds"))
b2a.thr=0.02
pl <- list()
for (i in 1:nrow(toplot))
{
test <- toplot %>% dplyr::slice(i)
pl[[i]] <- ggplot(test$signalr[[1]]) +
	rasterise(geom_point(data=test$signalr[[1]] ,aes(x=positions,y=Bprob,col="data.raw"),size=0.2,alpha=0.5,shape=16),dev="cairo",dpi=300)+
	geom_line(aes(x=positions,y=signal,col="data.smoothed"))+
	geom_line(data=test$RDP[[1]],aes(x=x,y=y,col="RDP_segment"))+
	geom_text(data=test$sl2[[1]],
		aes(x=sl.x,y=-0.02,col="RDP_segment",label=sl.pat2,fontface="bold",vjust=1),
		show.legend = F,size=2)+
	geom_hline(yintercept=b2a.thr,linetype="dashed",alpha=0.3) +
	geom_vline(data=test$forks[[1]],aes(xintercept=X0),linetype="longdash",alpha=0.3)+
	geom_vline(data=test$forks[[1]],aes(xintercept=X1),linetype="twodash",alpha=0.3)+
	geom_label(data=test$forks[[1]],
		aes(x=X0,y=0.95,fontface="bold",col="Fork_pulse",label="X0"),
		size=1.8, show.legend = F,label.padding=unit(0.3,"mm"))+
	geom_label(data=test$forks[[1]],
		aes(x=X1,y=0.95,fontface="bold",col="Fork_pulse",label="X1"),
		size=1.8, show.legend = F,label.padding=unit(0.3,"mm"))+
	geom_segment(data=test$forks[[1]],
		aes(x=X0,xend=X1,y=(0.8),yend=(0.8),col="Fork_pulse"),
		size=0.5,arrow=arrow(length = unit(0.1,"cm")), show.legend = F)+
	geom_text(data=test$forks[[1]],
		aes(x=(X0+X1)/2,y=(0.7),fontface="bold",col="Fork_pulse",label=speed),
		size=2, show.legend = F)+
	xlab(paste0(test$chrom," (kb)"))+
	scale_x_continuous(
		labels=scales::unit_format(big.mark ="",suffix="",scale=1e-3,sep=""),
		expand=c(0,0))+
	ylab("BrdU signal")+
	guides(col = guide_legend(title = "Legend",override.aes = list(lwd = 1,labels="")))+
	theme(legend.position = "right",legend.title=element_blank())+
	scale_color_manual("",labels = c("Raw data","Smoothed signal","RDP segment","Fork pulse"),breaks = c("data.raw","data.smoothed","RDP_segment","Fork_pulse"),values = mypal[c(2,1,3,5,6)])+
#	annotate(geom="text",x=(test$end-4e3),y=0.95,label=i,fontface="bold.italic",size=4)+
	coord_cartesian(ylim=c(-0.07,1))
}


p0 <- (pl[[1]]+
#labs(tag="b")+theme(plot.tag=element_text(face="bold"))+
pl[[2]]+theme(axis.title.y=element_blank(),axis.text.y=element_blank()))/
(pl[[3]]+
pl[[4]]+theme(axis.title.y=element_blank(),axis.text.y=element_blank()))/
(pl[[5]]+
pl[[6]]+theme(axis.title.y=element_blank(),axis.text.y=element_blank()))/
(pl[[7]]+
pl[[8]]+theme(axis.title.y=element_blank(),axis.text.y=element_blank()))/
(pl[[9]]+
pl[[10]]+theme(axis.title.y=element_blank(),axis.text.y=element_blank())) &
 theme(legend.position = "bottom")
p0 + plot_layout(guides = "collect")
ggsave(paste0(path_figures,"Figure1B.pdf"),h=8,w=7,device=cairo_pdf)
ggsave(paste0(path_figures,"Figure1B.png"),h=8,w=7)

p01 <- p0 + plot_layout(guides = "collect")
(plot_spacer()+p01) + plot_layout(heights=c(1,5))
ggsave(paste0(path_figures,"Figure1Bscale.pdf"),h=9,w=7,device=cairo_pdf)
ggsave(paste0(path_figures,"Figure1Bscale.png"),h=9,w=7)

