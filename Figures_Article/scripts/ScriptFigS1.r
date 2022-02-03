### Script FigS1
suppressMessages(library(tidyverse))
library(patchwork)
library(ggpubr)

theme_set(theme_bw())
mypal <- c(paletteer::paletteer_d("ggthemes::Classic_20"),"grey40")
`%+%` <- paste0
path_figures <- "/Users/ll/work/RStudioProjects/NanoForkSpeed/Figures_Article/figures/"
pathdata <- "/Users/ll/work/RStudioProjects/NanoForkSpeed/Figures_Article/data/"

pBc <- c(0,10,20,30,40,50,60,70,80,90,100)/100

toplot <- read_tsv(paste0(pathdata,"FigureS1A_data.tsv.gz"))
toplot$type <- factor(toplot$type, levels=c("Megalodon_BrdU","RepNano_CNN","RepNano_TM","DNAscent_v2","DNAscent_v1"))
id.mse <- function(pB,B.detected) {mean((pB-B.detected)^2)}
tomse <- split(toplot,toplot$type)
res_mse <- tibble(mse=sapply(tomse, function(x) id.mse(x$pBMS,x$B.mean)) %>% round(.,4),type=c("Megalodon_BrdU","RepNano_CNN","RepNano_TM","DNAscent_v2","DNAscent_v1"))


fs1a <- ggplot(toplot,aes(x=pBMS,y=B.mean,col=type))+
	geom_point(size=2,show.legend=F)+
	coord_cartesian(ylim=c(0,1),xlim=c(0,1))+
	geom_smooth(method="lm",se=F)+
	geom_abline(slope=1,intercept=0,linetype="dashed")+
	stat_regline_equation(aes(label =  paste(..eq.label.., ..rr.label.., sep = "~~~~")),label.y=c(1,0.9,0.8,0.7,0.6),show.legend=F)+
	scale_color_manual(breaks = c("Megalodon_BrdU","RepNano_CNN","RepNano_TM","DNAscent_v2","DNAscent_v1"),values = mypal[c(1,3,5,7,9)])+
	geom_text(data=res_mse,aes(x=0,y=c(0.95,0.85,0.75,0.65,0.55),label=paste0("mse=",mse),col=type),hjust=0,show.legend=F)+
	xlab("BrdU ratio (MS)")+
	ylab("BrdU ratio (nanopore)")+theme(plot.tag=element_text(face="bold"))



bin1k_meg3 <- read_tsv(paste0(pathdata,"FigureS1B_data.tsv.gz"))
fs1b <- ggplot(bin1k_meg3)+
	geom_freqpoly(aes(x=Bmean,col=factor(pB,levels=pBc),y=..ndensity..),binwidth=0.02)+
	coord_cartesian(xlim=c(0,1))+
	scale_x_continuous(breaks=c(0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1))+
	paletteer::scale_color_paletteer_d("ggthemes::Classic_20")+
	xlab("Mean BrdU ratio per kb")+
	ylab("Normalized Density")+
	ggtitle("Megalodon_BrdU")+theme(plot.tag=element_text(face="bold"))

bin1k_DS1 <- read_tsv(paste0(pathdata,"FigureS1C_data.tsv.gz"))
fs1c <- ggplot(bin1k_DS1)+
	geom_freqpoly(aes(x=Bmean,col=factor(pB,levels=pBc),y=..ndensity..),binwidth=0.02)+
	coord_cartesian(xlim=c(0,1))+
	scale_x_continuous(breaks=c(0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1))+
	paletteer::scale_color_paletteer_d("ggthemes::Classic_20")+
	xlab("Mean BrdU ratio per kb")+
	ylab("Normalized Density")+
	ggtitle("DNAscent_v1")+theme(plot.tag=element_text(face="bold"))

bin1k_DS2 <- read_tsv(paste0(pathdata,"FigureS1D_data.tsv.gz"))
fs1d <- ggplot(bin1k_DS2)+
	geom_freqpoly(aes(x=Bmean,col=factor(pB,levels=pBc),y=..ndensity..),binwidth=0.02)+
	coord_cartesian(xlim=c(0,1))+
	scale_x_continuous(breaks=c(0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1))+
	paletteer::scale_color_paletteer_d("ggthemes::Classic_20")+
	xlab("Mean BrdU ratio per kb")+
	ylab("Normalized Density")+
	ggtitle("DNAscent_v2")+theme(plot.tag=element_text(face="bold"))

bin1k_TM <- read_tsv(paste0(pathdata,"FigureS1E_data.tsv.gz"))
fs1e <- ggplot(bin1k_TM)+
	geom_freqpoly(aes(x=Bmean,col=factor(pB,levels=pBc),y=..ndensity..),binwidth=0.02)+
	coord_cartesian(xlim=c(0,1))+
	scale_x_continuous(breaks=c(0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1))+
	paletteer::scale_color_paletteer_d("ggthemes::Classic_20")+
	xlab("Mean BrdU ratio per kb")+
	ylab("Normalized Density")+
	ggtitle("RepNano_TM")+theme(plot.tag=element_text(face="bold"))

bin1k_CNN5 <- read_tsv(paste0(pathdata,"FigureS1F_data.tsv.gz"))
fs1f <- ggplot(bin1k_CNN5)+
	geom_freqpoly(aes(x=Bmean,col=factor(pB,levels=pBc),y=..ndensity..),binwidth=0.02)+
	coord_cartesian(xlim=c(0,1))+
	scale_x_continuous(breaks=c(0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1))+
	paletteer::scale_color_paletteer_d("ggthemes::Classic_20")+
	xlab("Mean BrdU ratio per kb")+
	ylab("Normalized Density")+
	ggtitle("RepNano_CNN")+theme(plot.tag=element_text(face="bold"))


fs1aa <- fs1a+ theme(legend.position = c(0.8, 0.2),legend.title=element_blank(),legend.background=element_rect(fill = NA),legend.key.size = unit(0.3, 'cm'))

combined <- fs1aa +fs1b+ fs1c+ fs1d+ fs1e+ fs1f & guides(colour=guide_legend(title="BrdU\nratio\n(b to f)"))
combined[[1]] <- combined[[1]] + plot_layout(guides="keep")
combined + plot_layout(ncol=2,guides="collect")+plot_annotation(tag_levels = 'a')

ggsave(paste0(path_figures,"FigureS1.pdf"),h=14,w=12)

