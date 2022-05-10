# Script for figures S2A and S2B
#load required library
suppressMessages(library(tidyverse))

path_figures <- "./Figures_Article/figures/"
pathdata <- "./Figures_Article/data/"


# Import data for growth curve
data_gc<-read.csv(paste0(pathdata,"FigureS2A_data.csv"),sep=",",header=TRUE)

# Transform Strain by levels for plotting
data_gc%>%mutate(Strain=factor(Strain,levels = c("WT","BT1","MCM869")))->data_gc

# Plot growth curve saved as a pdf file
pdf(file = paste0(path_figures,"FigureS2A.pdf"), width = 8, height = 4)
ggplot(data_gc)+
geom_line(aes(x=Time,y=log(OD600,base=10),group=Strain,color=Strain),size=0.5)+
geom_point(aes(x=Time,y=log(OD600,base=10),group=Strain,color=Strain),size=1)+
scale_color_manual(values=c("chartreuse4","blue","firebrick1"))+
xlab("Time (minute)")+
ylab("LOG10(OD600)")+
theme_bw()
dev.off()

# Import data for doubling time
data_dt<-read.csv(paste0(pathdata,"FigureS2B_data.csv"),sep=",",header=TRUE)

# Transform Strain by levels for plotting
data_dt%>%mutate(Strain=factor(Strain,levels = c("WT","BT1","MCM869")))->data_dt

# Compute median, min and max doubling time for each strain
data_dt%>%group_by(Strain)%>%
summarise(med_dt=median(Doubling_time),max_dt=max(Doubling_time),min_dt=min(Doubling_time))->stat_dt

# Plot doubling time stats saved as a pdf file
pdf(file = paste0(path_figures,"FigureS2B.pdf"), width = 8, height = 4)
ggplot(data_dt)+
geom_point(aes(x=Strain,y=Doubling_time,color=Strain),size=1.5)+
ylim(0,150)+
ylab("Doubling time (minute)")+
xlab("Strain")+
geom_errorbar(data=stat_dt,aes(ymin=min_dt,ymax=max_dt,x=Strain),width = 0.15,alpha=0.9)+
geom_point(data=stat_dt,aes(x=Strain,y=med_dt,group=Strain),shape=3,size=2,alpha=0.9)+
scale_color_manual(values=c( "chartreuse4","blue","firebrick1"))+
theme_bw()
dev.off()


