## script figure S3
suppressMessages(library(GenomicRanges))
seqinf <- readRDS("/Users/ll/work/Ori/seqinfS288CrDNA.rds")
suppressMessages(library(rtracklayer))
suppressMessages(library(tidyverse))
library(patchwork)
library(ggprism)
library(ggcorrplot)

theme_set(theme_bw())
mypal <- c(paletteer::paletteer_d("ggthemes::Classic_20"),"grey40")
`%+%` <- paste0
path_figures <- "./Figures_Article/figures/"
pathdata <- "./Figures_Article/data/"

## Figure S3B
toplot <- read_tsv(paste0(pathdata,"FigureS3B_data.tsv.gz"))
mcornt <- as.data.frame(toplot)
rownames(mcornt) <- mcornt$X0
mcornt <- mcornt[,1:3]

pE <- ggcorrplot(mcornt,lab=T,lab_size=8,type="upper",outline.color = NA)+
	scale_fill_gradient2(limit=c(0,1),low="blue",mid="white",high="red",midpoint=0.5)+
	labs(tag="b")+
	theme(plot.tag=element_text(face="bold"),panel.grid=element_blank())
ggsave(paste0(path_figures,"FigureS3B.pdf"),h=6,w=6,pE)


# Figure S3A
ARS <- import(paste0(pathdata,"ARS_newman.bed"))
chrom <- as(seqinf,"GRanges")[-17]

i=15
	ROI <- chrom[i]
	ARS_ROI <- as_tibble(ARS[overlapsAny(ARS,ROI)]) %>% dplyr::rename(featname=name)
bs=100


pl_ARS <- ggplot(ARS_ROI)+
	geom_rect(aes(xmin=start,xmax=end,ymin=0,ymax=1),col=mypal[13],fill=mypal[13])+
	xlab(paste0(as.character(seqnames(ROI))," (kb)"))+
	ylab("ORI")+
	scale_x_continuous(
		guide = "prism_minor",
		labels=scales::unit_format(big.mark ="",suffix="",scale=1e-3,sep=""),
		limits=c(1,seqlengths(seqinf)[i]),
		breaks=seq(0,seqlengths(seqinf)[i],100000),
		minor_breaks=seq(0,seqlengths(seqinf)[i],20000),
		expand=c(0,0))+
	theme(axis.ticks.y = element_blank(),axis.text.y = element_blank(),panel.grid.major.y=element_blank(),panel.grid.minor.y=element_blank())

rfd_nfs <- import(paste0(pathdata,"BT1_NFS_RFD_bs0.001k_lr0_wiNA.bw"))
toplot_rfd <- tibble(pos=start(ROI):end(ROI), value=as.numeric(unlist(coverage(rfd_nfs,weight=rfd_nfs$score)[ROI]))) %>% mutate(type=ifelse(sign(value)>=0,"RFD_R","RFD_L"))

toplot_rfd2 <- toplot_rfd %>% mutate(pos = floor(pos/bs)*bs+1) %>%
		group_by(pos) %>%
		summarise(value = mean(value,na.rm=T), .groups = "drop") %>% mutate(type=ifelse(sign(value)>=0,"RFD_R","RFD_L"))

pl_RFD_NFS <- ggplot(toplot_rfd2)+
	geom_point(aes(x=pos,y=value,col=type),size=0.1)+
	scale_x_continuous(
		guide = "prism_minor",
		labels=scales::unit_format(big.mark ="",suffix="",scale=1e-3,sep=""),
		limits=c(1,seqlengths(seqinf)[i]),
		breaks=seq(0,seqlengths(seqinf)[i],100000),
		minor_breaks=seq(0,seqlengths(seqinf)[i],20000),
		expand=c(0,0))+
	coord_cartesian(ylim=c(-1,1))+
	xlab(paste0(as.character(seqnames(ROI))," (kb)"))+
	ylab("RFD")+
	labs(tag="a")+theme(plot.tag=element_text(face="bold"))+
	ggtitle("NFS (BT1)")+
	scale_color_manual("RFD",breaks=c("RFD_R","RFD_L"),values=c("red","blue"))+
	guides(colour = guide_legend(override.aes = list(size=2)))

rfd_OK <- import(paste0(pathdata,"MCM869_OKseq_RFD_rm.bw"))
toplot_OK <- tibble(pos=start(ROI):end(ROI),
value=as.numeric(unlist(coverage(rfd_OK,weight=rfd_OK$score)[ROI]))) %>% mutate(type=ifelse(sign(value)>=0,"RFD_R","RFD_L"))
toplot_OK2 <- toplot_OK %>% mutate(pos = floor(pos/bs)*bs+1) %>%
		group_by(pos) %>%
		summarise(value = mean(value,na.rm=T), .groups = "drop") %>% mutate(type=ifelse(sign(value)>=0,"RFD_R","RFD_L"))

pl_RFD_OK <- ggplot(toplot_OK2)+
	geom_point(aes(x=pos,y=value,col=type),size=0.1)+
	scale_x_continuous(
		guide = "prism_minor",
		labels=scales::unit_format(big.mark ="",suffix="",scale=1e-3,sep=""),
		limits=c(1,seqlengths(seqinf)[i]),
		breaks=seq(0,seqlengths(seqinf)[i],100000),
		minor_breaks=seq(0,seqlengths(seqinf)[i],20000),
		expand=c(0,0))+
	coord_cartesian(ylim=c(-1,1))+
	xlab(paste0(as.character(seqnames(ROI))," (kb)"))+
	ylab("RFD")+
	ggtitle("OK-seq (MCM869)")+
	scale_color_manual("RFD",breaks=c("RFD_R","RFD_L"),values=c("red","blue"))+
	guides(colour = guide_legend(override.aes = list(size=2)))

rfd_FS <- import(paste0(pathdata,"MCM869_ForkSeq_RFD_bs0.001k_lr0_wiNA.bw"))
toplot_FS <- tibble(pos=start(ROI):end(ROI),
value=as.numeric(unlist(coverage(rfd_FS,weight=rfd_FS$score)[ROI]))) %>% mutate(type=ifelse(sign(value)>=0,"RFD_R","RFD_L"))
toplot_FS2 <- toplot_FS %>% mutate(pos = floor(pos/bs)*bs+1) %>%
		group_by(pos) %>%
		summarise(value = mean(value,na.rm=T), .groups = "drop") %>% mutate(type=ifelse(sign(value)>=0,"RFD_R","RFD_L"))

pl_RFD_FS <- ggplot(toplot_FS2)+
	geom_point(aes(x=pos,y=value,col=type),size=0.1)+
	scale_x_continuous(
		guide = "prism_minor",
		labels=scales::unit_format(big.mark ="",suffix="",scale=1e-3,sep=""),
		limits=c(1,seqlengths(seqinf)[i]),
		breaks=seq(0,seqlengths(seqinf)[i],100000),
		minor_breaks=seq(0,seqlengths(seqinf)[i],20000),
		expand=c(0,0))+
	coord_cartesian(ylim=c(-1,1))+
	xlab(paste0(as.character(seqnames(ROI))," (kb)"))+
	ylab("RFD")+
	ggtitle("FORK-seq (MCM869)")+
	scale_color_manual("RFD",breaks=c("RFD_R","RFD_L"),values=c("red","blue"))+
	guides(colour = guide_legend(override.aes = list(size=2)))

pl_RFD <- pl_RFD_NFS/pl_RFD_FS/pl_RFD_OK & theme(legend.position = "right",axis.text.x = element_blank(),axis.title.x=element_blank())
pA <- pl_RFD/pl_ARS+plot_layout(ncol = 1, heights = c(3,3,3,1),guides = "collect")
ggsave(paste0(path_figures,"FigureS3A.pdf"),w=7,h=4,pA)

### fig S3C
library(Hmisc)

toplot <- read_tsv(paste0(pathdata,"FigureS3C_data.tsv.gz"))

pdf(paste0(path_figures,"FigureS3C.pdf"),height=5,width=6)
Ecdf(toplot %>% filter(type=="ini.all") %>% pull(dtac),col=mypal[7],q=0.5,xlim=c(0,30000),xlab="Distance to nearest ORI centre",ylab="ECDF",lwd=1,subtitles=F)
Ecdf(toplot %>% filter(type=="rd.ini.all") %>% pull(dtac),col=mypal[8],q=0.5,add=T,lwd=1,lty=2,subtitles=F)
Ecdf(toplot %>% filter(type=="ter.all") %>% pull(dtac),col=mypal[1],q=0.5,add=T,lwd=1,subtitles=F)
Ecdf(toplot %>% filter(type=="rd.ter.all") %>% pull(dtac),col=mypal[2],q=0.5,add=T,lwd=1,lty=3,subtitles=F)
legend("bottomright",legend=c(paste0("NFS_ini (n=",length(toplot %>% filter(type=="ini.all") %>% pull(dtac)),", med=",round(median(toplot %>% filter(type=="ini.all") %>% pull(dtac))),"bp)"),paste0("NFS_ini shuffled (n=",length(toplot %>% filter(type=="rd.ini.all") %>% pull(dtac)),", med=",round(median(toplot %>% filter(type=="rd.ini.all") %>% pull(dtac))),"bp)"),paste0("NFS_ter (n=",length(toplot %>% filter(type=="ter.all") %>% pull(dtac)),", med=",round(median(toplot %>% filter(type=="ter.all") %>% pull(dtac))),"bp)"),paste0("NFS_ter shuffled (n=",length(toplot %>% filter(type=="rd.ter.all") %>% pull(dtac)),", med=",round(median(toplot %>% filter(type=="rd.ter.all") %>% pull(dtac))),"bp)")),bty='n',text.col=mypal[c(7,8,1,2)],cex=0.75)
mtext(expression(bold("c")),side=3,line=1,at=-5000,cex=1.4)
dev.off()

