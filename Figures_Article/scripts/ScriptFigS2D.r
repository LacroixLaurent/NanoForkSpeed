### Script FigS2D
suppressMessages(library(tidyverse))
library(patchwork)
library(ggrastr)
theme_set(theme_bw())

mypal <- c(paletteer::paletteer_d("ggthemes::Classic_20"),"grey40")
`%+%` <- paste0
path_figures <- "/Users/ll/work/RStudioProjects/NanoForkSpeed/Figures_Article/figures/"
pathdata <- "/Users/ll/work/RStudioProjects/NanoForkSpeed/Figures_Article/data/"

### figS2D
toplot <- readRDS(paste0(pathdata,"FigureS2D_BT3_data.rds"))
pl_BT3 <- list()
for (i in 1:nrow(toplot))
{
test <- toplot %>% dplyr::slice(i)
pl_BT3[[i]] <- ggplot(test$signalr[[1]]) +
	rasterise(geom_point(data=test$signalr[[1]] ,aes(x=positions,y=Bprob,col="data.raw"),size=0.2,alpha=0.5,shape=16),dev="cairo",dpi=300)+
geom_line(aes(x=positions,y=signal,col="data.smoothed"))+
xlab(paste0(test$chrom," (kb)"))+
	scale_x_continuous(
		labels=scales::unit_format(big.mark ="",suffix="",scale=1e-3,sep="",accuracy=1)#,
		#expand=c(0,0)
		)+
ylab("BrdU signal")+
guides(col = guide_legend(title = "Legend",override.aes = list(lwd = 1,labels="")))+
	theme(legend.position = "right",legend.title=element_blank())+
	scale_color_manual(labels = c("Raw data","Smoothed signal"),values = mypal[c(2,1)])+
	coord_cartesian(ylim=c(0,1))
}

toplot <- readRDS(paste0(pathdata,"FigureS2D_BT2_data.rds"))
pl_BT2 <- list()
for (i in 1:nrow(toplot))
{
test <- toplot %>% dplyr::slice(i)
pl_BT2[[i]] <- ggplot(test$signalr[[1]]) +
	rasterise(geom_point(data=test$signalr[[1]] ,aes(x=positions,y=Bprob,col="data.raw"),size=0.2,alpha=0.5,shape=16),dev="cairo",dpi=300)+
geom_line(aes(x=positions,y=signal,col="data.smoothed"))+
xlab(paste0(test$chrom," (kb)"))+
	scale_x_continuous(
		labels=scales::unit_format(big.mark ="",suffix="",scale=1e-3,sep="",accuracy=1)#,
		#expand=c(0,0)
		)+
ylab("BrdU signal")+
guides(col = guide_legend(title = "Legend",override.aes = list(lwd = 1,labels="")))+
theme(legend.position = "right",legend.title=element_blank())+
	scale_color_manual(labels = c("Raw data","Smoothed signal"),values = mypal[c(2,1)])+
	coord_cartesian(ylim=c(0,1))
}

toplot <- readRDS(paste0(pathdata,"FigureS2D_BT1_data.rds"))
pl_BT1 <- list()
for (i in 1:nrow(toplot))
{
test <- toplot %>% dplyr::slice(i)
pl_BT1[[i]] <- ggplot(test$signalr[[1]]) +
	rasterise(geom_point(data=test$signalr[[1]] ,aes(x=positions,y=Bprob,col="data.raw"),size=0.2,alpha=0.5,shape=16),dev="cairo",dpi=300)+
geom_line(aes(x=positions,y=signal,col="data.smoothed"))+
xlab(paste0(test$chrom," (kb)"))+
	scale_x_continuous(
		labels=scales::unit_format(big.mark ="",suffix="",scale=1e-3,sep="",accuracy=1)#,
		#expand=c(0,0)
		)+
ylab("BrdU signal")+
guides(col = guide_legend(title = "Legend",override.aes = list(lwd = 1,labels="")))+
	theme(legend.position = "right",legend.title=element_blank())+
	scale_color_manual(labels = c("Raw data","Smoothed signal"),values = mypal[c(2,1)])+
	coord_cartesian(ylim=c(0,1))
}

toplot <- readRDS(paste0(pathdata,"FigureS2D_MCM_data.rds"))
pl_MCM <- list()
for (i in 1:nrow(toplot))
{
test <- toplot %>% dplyr::slice(i)
pl_MCM[[i]] <- ggplot(test$signalr[[1]]) +
	rasterise(geom_point(data=test$signalr[[1]] ,aes(x=positions,y=Bprob,col="data.raw"),size=0.2,alpha=0.5,shape=16),dev="cairo",dpi=300)+
geom_line(aes(x=positions,y=signal,col="data.smoothed"))+
xlab(paste0(test$chrom," (kb)"))+
	scale_x_continuous(
		labels=scales::unit_format(big.mark ="",suffix="",scale=1e-3,sep="",accuracy=1)#,
		#expand=c(0,0)
		)+
ylab("BrdU signal")+
guides(col = guide_legend(title = "Legend",override.aes = list(lwd = 1,labels="")))+
	theme(legend.position = "right",legend.title=element_blank())+
	scale_color_manual(labels = c("Raw data","Smoothed signal"),values = mypal[c(2,1)])+
coord_cartesian(ylim=c(0,1))
}


p0 <- pl_MCM[[1]]+pl_MCM[[2]]+pl_MCM[[3]]+pl_MCM[[4]]+pl_MCM[[5]]+pl_BT1[[1]]+pl_BT1[[2]]+pl_BT1[[3]]+pl_BT1[[4]]+pl_BT1[[5]]+pl_BT2[[1]]+pl_BT2[[2]]+pl_BT2[[3]]+pl_BT2[[4]]+pl_BT2[[5]]+pl_BT3[[1]]+pl_BT3[[2]]+pl_BT3[[3]]+pl_BT3[[4]]+pl_BT3[[5]] & theme(legend.position="bottom")
p0[[1]] <- p0[[1]]+ggtitle("MCM869")#+
#labs(tag="d")+theme(plot.tag=element_text(face="bold"))
p0[[6]] <- p0[[6]]+ggtitle("BT1")
p0[[11]] <- p0[[11]]+ggtitle("BT2")
p0[[16]] <- p0[[16]]+ggtitle("BT3")

p0 +plot_layout(byrow = FALSE,ncol=4,guides="collect")
#ggsave(paste0(path_figures,"FigureS2D.png"),h=9,w=14)
ggsave(paste0(path_figures,"FigureS2D.pdf"),h=9,w=14,device=cairo_pdf)

