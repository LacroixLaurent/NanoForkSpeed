### Script FigS6
suppressMessages(library(tidyverse))
library(patchwork)
library(ggprism)

theme_set(theme_bw())
mypal <- c(paletteer::paletteer_d("ggthemes::Classic_20"),"grey40")
`%+%` <- paste0
path_figures <- "/Users/ll/work/RStudioProjects/NanoForkSpeed/Figures_Article/figures/"
pathdata <- "/Users/ll/work/RStudioProjects/NanoForkSpeed/Figures_Article/data/"

### Manip
toplot0 <- readRDS(paste0(pathdata,"FigureS6E_data.rds"))

pl_BT1 <- list()
for (i in 1:nrow(toplot0))
{
test <- toplot0 %>% dplyr::slice(i)
pl_BT1[[i]] <- ggplot(test$signalr[[1]]) +
#geom_point(aes(x=positions,y=Bprob,col="data.raw"),size=0.2,alpha=0.3)+
geom_line(aes(x=positions,y=signal,col="data.smoothed"))+
xlab(paste0(test$chrom," (kb)"))+
	scale_x_continuous(
		labels=scales::unit_format(big.mark ="",suffix="",scale=1e-3,sep="",accuracy=1),
		expand=c(0,0))+
ylab("BrdU signal")+
guides(col = guide_legend(title = "Legend",override.aes = list(lwd = 1,labels="")))+
theme(legend.position = "right")+
scale_color_manual(breaks = c("data.raw","data.smoothed"),values = mypal[c(2,1)])+
coord_cartesian(ylim=c(0,1))
}

# multi
toplot1 <- readRDS(paste0(pathdata,"FigureS6D_data.rds"))

pl_Multi <- list()
for (i in 1:nrow(toplot1))
{
test <- toplot1 %>% dplyr::slice(i)
pl_Multi[[i]] <- ggplot(test$signalr[[1]]) +
geom_line(aes(x=positions,y=signal,col="data.smoothed"))+
xlab("Coordinates (kb)")+
	scale_x_continuous(
		labels=scales::unit_format(big.mark ="",suffix="",scale=1e-3,sep="",accuracy=1),
		expand=c(0,0))+
ylab("BrdU signal")+
guides(col = guide_legend(title = "Legend",override.aes = list(lwd = 1,labels="")))+
theme(legend.position = "right")+
scale_color_manual(breaks = c("data.raw","data.smoothed"),values = mypal[c(2,1)])+
coord_cartesian(ylim=c(0,1))
}

# multiGT
toplot2 <- readRDS(paste0(pathdata,"FigureS6C_data.rds"))

pl_MultiGT <- list()
for (i in 1:nrow(toplot2))
{
test <- toplot2 %>% dplyr::slice(i)
pl_MultiGT[[i]] <- ggplot(test$signalr[[1]]) +
geom_line(aes(x=positions,y=signal,col="data.smoothed"))+
xlab("Coordinates (kb)")+
	scale_x_continuous(
		labels=scales::unit_format(big.mark ="",suffix="",scale=1e-3,sep="",accuracy=1),
		expand=c(0,0))+
ylab("BrdU signal")+
guides(col = guide_legend(title = "Legend",override.aes = list(lwd = 1,labels="")))+
theme(legend.position = "right")+
scale_color_manual(breaks = c("data.raw","data.smoothed"),values = mypal[c(2,1)])+
coord_cartesian(ylim=c(0,1))
}

# mono
toplot3 <- readRDS(paste0(pathdata,"FigureS6B_data.rds"))

pl_Mono <- list()
for (i in 1:nrow(toplot3))
{
test <- toplot3 %>% dplyr::slice(i)
pl_Mono[[i]] <- ggplot(test$signalr[[1]]) +
geom_line(aes(x=positions,y=signal,col="data.smoothed"))+
xlab("Coordinates (kb)")+
	scale_x_continuous(
		labels=scales::unit_format(big.mark ="",suffix="",scale=1e-3,sep="",accuracy=1),
		expand=c(0,0))+
ylab("BrdU signal")+
guides(col = guide_legend(title = "Legend",override.aes = list(lwd = 1,labels="")))+
theme(legend.position = "right")+
scale_color_manual(breaks = c("data.raw","data.smoothed"),values = mypal[c(2,1)])+
coord_cartesian(ylim=c(0,1))
}

# monoGT
toplot4 <- readRDS(paste0(pathdata,"FigureS6A_data.rds"))

pl_MonoGT <- list()
for (i in 1:nrow(toplot4))
{
test <- toplot4 %>% dplyr::slice(i)
pl_MonoGT[[i]] <- ggplot(test$signalr[[1]]) +
geom_line(aes(x=positions,y=signal,col="data.smoothed"))+
xlab("Coordinates (kb)")+
	scale_x_continuous(
		labels=scales::unit_format(big.mark ="",suffix="",scale=1e-3,sep="",accuracy=1),
		expand=c(0,0))+
ylab("BrdU signal")+
guides(col = guide_legend(title = "Legend",override.aes = list(lwd = 1,labels="")))+
theme(legend.position = "right")+
scale_color_manual(breaks = c("data.raw","data.smoothed"),values = mypal[c(2,1)])+
coord_cartesian(ylim=c(0,1))
}


p0 <- pl_MonoGT[[1]]+pl_MonoGT[[2]]+pl_MonoGT[[3]]+pl_MonoGT[[4]]+pl_MonoGT[[5]]+pl_Mono[[1]]+pl_Mono[[2]]+pl_Mono[[3]]+pl_Mono[[5]]+pl_Mono[[4]]+pl_MultiGT[[1]]+pl_MultiGT[[2]]+pl_MultiGT[[3]]+pl_MultiGT[[4]]+pl_MultiGT[[5]]+pl_Multi[[1]]+pl_Multi[[2]]+pl_Multi[[3]]+pl_Multi[[4]]+pl_Multi[[5]]+pl_BT1[[1]]+pl_BT1[[2]]+pl_BT1[[3]]+pl_BT1[[4]]+pl_BT1[[5]] & theme(legend.position="none")
p0[[1]] <- p0[[1]]+ggtitle("Single forks without noise")+labs(tag="a")+theme(plot.tag=element_text(face="bold"))
p0[[6]] <- p0[[6]]+ggtitle("Single forks with noise")+labs(tag="b")+theme(plot.tag=element_text(face="bold"))
p0[[11]] <- p0[[11]]+ggtitle("Multiple forks without noise")+labs(tag="c")+theme(plot.tag=element_text(face="bold"))
p0[[16]] <- p0[[16]]+ggtitle("Multiple forks with noise")+labs(tag="d")+theme(plot.tag=element_text(face="bold"))
p0[[21]] <- p0[[21]]+ggtitle("BT1")+labs(tag="e")+theme(plot.tag=element_text(face="bold"))

p0 +plot_layout(byrow = FALSE,ncol=5,guides="collect")
ggsave(paste0(path_figures,"FigureS6A-E.pdf"),h=9,w=18)

### figure simu
toplot <- read_tsv(paste0(pathdata,"FigureS6F_data.tsv.gz"))

fs6f <- ggplot(toplot,aes(x=x,y=value,col=name))+
	geom_line()+
	scale_colour_manual("",labels=c("Experimental data","Simulated data"),values=mypal[c(9,7)])+
	xlab("Position (nt)")+
	ylab("BrdU signal")+
	labs(tag="f")+
	theme(plot.tag=element_text(face="bold"))

toplot <- read_tsv(paste0(pathdata,"FigureS6G_data.tsv.gz"))

fs6g <- ggplot(toplot,aes(x=distance,y=value,col=name))+
	geom_line()+
	scale_colour_manual("",labels=c("Experimental data","Simulated data"),values=mypal[c(9,7)])+
	xlab("Distance (nt)")+
	ylab("Autocorrelation")+
	labs(tag="g")+
	theme(plot.tag=element_text(face="bold"))

fs6f+fs6g +plot_layout(guides="collect")
ggsave(paste0(path_figures,"FigureS6FG.pdf"),h=4,w=9)
