### Script for figures 3 and S7
suppressMessages(library(tidyverse))
theme_set(theme_bw())
library(patchwork)
library(ggdist)
mypal <- c(paletteer::paletteer_d("ggthemes::Classic_20"),"grey40")
`%+%` <- paste0
path_figures <- "/Users/ll/work/RStudioProjects/NanoForkSpeed/Figures_Article/figures/"
pathdata <- "/Users/ll/work/RStudioProjects/NanoForkSpeed/Figures_Article/data/"

## fig3A
toplotA <- read_tsv(paste0(pathdata,"Figure3A_data.tsv.gz"))%>% mutate(type=fct_relevel(type,c("Single tracks without noise","Single tracks with noise","Multiple tracks without noise","Multiple tracks with noise")))
totextA <- toplotA %>% group_by(type) %>% summarise(n=n())
tomedA <- toplotA %>% group_by(type) %>% summarise(med=round(median(speed_error)))

f3a <- ggplot(toplotA,aes(x=type,y=speed_error))+
	stat_slab(aes(fill=type),col=NA,alpha=0.6,scale=0.8,normalize="groups")+
	stat_pointinterval(.width=c(.5,.95),col="grey40")+
	ylim(c(-1500,1500))+
	ylab("Speed error (bp/min)")+
	geom_text(data=totextA,aes(x=type,y=-1500,label=n),fontface="italic") +
	geom_text(data=tomedA,aes(x=type,y=1200,label=med),col="red")+
	scale_fill_manual("",values=mypal[c(2,1,8,7)])+
	theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())+
	labs(tag="a")+
	theme(plot.tag=element_text(face="bold"))
toplotA %>% group_by(type) %>% summarise(iqr=IQR(speed_error))
#  type                            iqr
#1 Single tracks without noise    123.
#2 Single tracks with noise       371.
#3 Multiple tracks without noise  186
#4 Multiple tracks with noise     438.

### Figure3B
toplotB <- read_tsv(paste0(pathdata,"Figure3B_data.tsv.gz"))
set.seed(123)
g4 <- tibble(g=rnorm(n=100000,mean=2486,sd=150))
f3b <- ggplot(toplotB)+
	geom_density(aes(x=speed0,col="Deconvolved distribution"))+
	geom_density(data=g4,aes(x=g,y=..density..*0.63,col="Fitted mode\n(m=2486,sd=150,p=63%)"))+
	coord_cartesian(xlim=c(0,4000))+
	scale_colour_manual("",values=mypal[c(6,5)],breaks=c("Deconvolved distribution","Fitted mode\n(m=2486,sd=150,p=63%)"))+
	xlab("Speed")+
	ylab("Density")+
	labs(tag="b")+
	theme(plot.tag=element_text(face="bold"))+
	xlab("Speed (bp/min)")


## fig3C
toplotC <- read_tsv(paste0(pathdata,"Figure3C_data.tsv.gz"))%>% mutate(type=fct_relevel(type,c("Experimental data","Multiple tracks with noise","Multiple tracks without noise","True Speed")))

f3c <- ggplot(toplotC)+
	stat_bin(aes(x=speed,y=..density..,col=type),geom="step",binwidth=200,position="identity")+
	stat_bin(data=toplotC %>% filter(type %in% c("Multiple tracks with noise")),aes(x=speed,y=..density..,col=type),geom="step",binwidth=200,position="identity")+
	stat_bin(data=toplotC %>% filter(type %in% c("Experimental data")),aes(x=speed,y=..density..,col=type),geom="step",binwidth=200,position="identity")+
	coord_cartesian(xlim=c(0,4000))+
	scale_colour_manual("",values=mypal[c(9,7,8,6)])+
	labs(tag="c")+
	theme(plot.tag=element_text(face="bold"))+
	ylab("Density")+
	xlab("Speed (bp/min)")


p0 <- f3a+f3b+f3c+ plot_layout(ncol=1)
ggsave(paste0(path_figures,"Figure3.pdf"),h=9,w=6,p0)

### fig s7
toplotS7A <- read_tsv(paste0(pathdata,"FigureS7A_data.tsv.gz"))
totext <- toplotS7A %>% group_by(gp) %>% summarise(n=n(),med=round(median(speed_error)),med2=round(median(speed_error2),2))

fs7a <- ggplot(toplotS7A)+
	geom_boxplot(aes(x=gp,y=speed_error,group=gp),outlier.shape=NA)+
	geom_hline(aes(yintercept=-250),linetype="dotted")+
	geom_hline(aes(yintercept=250),linetype="dotted")+
	geom_hline(aes(yintercept=0),linetype="dashed")+
	geom_text(data=totext,aes(x=gp,y=-1300,label=n),fontface="italic",size=3) +
	geom_text(data=totext,aes(x=gp,y=1000,label=med),col="red",size=3) +
	geom_text(data=totext,aes(x=gp,y=800,label=paste0("(",med2,")")),col="red",fontface="italic",size=3) +
	coord_cartesian(ylim=c(-1500,1000),xlim=c(0,3000))+
	labs(tag="a")+
	theme(plot.tag=element_text(face="bold"))+
	theme(axis.title.x=element_blank(),axis.text.x=element_blank())+
	ggtitle("Single tracks without noise")+
	ylab("Speed error (bp/min)")+
	xlab("True Speed (bp/min)")

toplotS7B <- read_tsv(paste0(pathdata,"FigureS7B_data.tsv.gz"))
totext <- %>% group_by(gp) %>% summarise(n=n(),med=round(median(speed_error)),med2=round(median(speed_error2),2))

fs7b <- ggplot(toplotS7B)+
	geom_boxplot(aes(x=gp,y=speed_error,group=gp),outlier.shape=NA)+
	geom_hline(aes(yintercept=-250),linetype="dotted")+
	geom_hline(aes(yintercept=250),linetype="dotted")+
	geom_hline(aes(yintercept=0),linetype="dashed")+
	geom_text(data=totext,aes(x=gp,y=-1300,label=n),fontface="italic",size=3) +
	geom_text(data=totext,aes(x=gp,y=1000,label=med),col="red",size=3) +
	geom_text(data=totext,aes(x=gp,y=800,label=paste0("(",med2,")")),col="red",fontface="italic",size=3) +
	coord_cartesian(ylim=c(-1500,1000),xlim=c(0,3000))+
	labs(tag="b")+
	theme(plot.tag=element_text(face="bold"))+
	theme(axis.title.x=element_blank(),axis.text.x=element_blank())+
	ggtitle("Single tracks with noise")+
	ylab("Speed error (bp/min)")+
	xlab("True Speed (bp/min)")

toplotS7C <- read_tsv(paste0(pathdata,"FigureS7C_data.tsv.gz"))
totext <- toplotS7C %>% group_by(gp) %>% summarise(n=n(),med=round(median(speed_error)),med2=round(median(speed_error2),2))

fs7c <- ggplot(toplotS7C)+
	geom_boxplot(aes(x=gp,y=speed_error,group=gp),outlier.shape=NA)+
	geom_hline(aes(yintercept=-250),linetype="dotted")+
	geom_hline(aes(yintercept=250),linetype="dotted")+
	geom_hline(aes(yintercept=0),linetype="dashed")+
	geom_text(data=totext,aes(x=gp,y=-1300,label=n),fontface="italic",size=3) +
	geom_text(data=totext,aes(x=gp,y=1000,label=med),col="red",size=3) +
	geom_text(data=totext,aes(x=gp,y=800,label=paste0("(",med2,")")),col="red",fontface="italic",size=3) +
	coord_cartesian(ylim=c(-1500,1000),xlim=c(0,3000))+
	labs(tag="c")+
	theme(plot.tag=element_text(face="bold"))+
	theme(axis.title.x=element_blank(),axis.text.x=element_blank())+
	ggtitle("Multiple tracks without noise")+
	ylab("Speed error (bp/min)")+
	xlab("True Speed (bp/min)")

toplotS7D <- read_tsv(paste0(pathdata,"FigureS7D_data.tsv.gz"))
totext <- toplotS7D %>% group_by(gp) %>% summarise(n=n(),med=round(median(speed_error)),med2=round(median(speed_error2),2))

fs7d <- ggplot(toplotS7D)+
	geom_boxplot(aes(x=gp,y=speed_error,group=gp),outlier.shape=NA)+
	geom_hline(aes(yintercept=-250),linetype="dotted")+
	geom_hline(aes(yintercept=250),linetype="dotted")+
	geom_hline(aes(yintercept=0),linetype="dashed")+
	geom_text(data=totext,aes(x=gp,y=-1300,label=n),fontface="italic",size=3) +
	geom_text(data=totext,aes(x=gp,y=1000,label=med),col="red",size=3) +
	geom_text(data=totext,aes(x=gp,y=800,label=paste0("(",med2,")")),col="red",fontface="italic",size=3) +
	coord_cartesian(ylim=c(-1500,1000),xlim=c(0,3000))+
	labs(tag="d")+
	theme(plot.tag=element_text(face="bold"))+
	ggtitle("Multiple tracks with noise")+
	ylab("Speed error (bp/min)")+
	xlab("True Speed (bp/min)")+
	scale_x_continuous(labels=c("0"="0-100","1000"="1000-1100","2000"="2000-2100","3000"="3000-3100"))

p0 <- fs7a+fs7b+fs7c+fs7d + plot_layout(ncol=1)

ggsave(paste0(path_figures,"FigureS7.pdf"),h=10,w=16,p0)

