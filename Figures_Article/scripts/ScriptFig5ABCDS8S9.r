### Script Fig5, S8 and S9
# LL 20220202
## Figure 5ABCD
suppressMessages(library(GenomicRanges))
suppressMessages(library(rtracklayer))
suppressMessages(library(tidyverse))
library(patchwork)
library(ggprism)

theme_set(theme_bw())
mypal <- c(paletteer::paletteer_d("ggthemes::Classic_20"),"grey40")
`%+%` <- paste0
setwd("/Users/ll/work/Ori/NFS_paper/")
#path_figures <- "/Users/ll/work/RStudioProjects/NanoForkSpeed/Figures_Article/figures/"
#pathdata <- "/Users/ll/work/RStudioProjects/NanoForkSpeed/Figures_Article/data/"
pathdata <- "/Users/ll/work/Ori/NFS_paper/GitHub_upload/data/"
path_figures <- "/Users/ll/work/Ori/NFS_paper/GitHub_upload/figures/"

seqinf <- readRDS(paste0(pathdata,"seqinfS288CrDNA.rds"))
root_title <- "BT1_wt"
chrom <- as(seqinf,"GRanges")[-c(17,18)]

### binned speed
speed_leg <- factor(c("speed","shuffled"),levels=c("speed","shuffled"))
speed_pal <- mypal[c(1,7)]
names(speed_pal) <- speed_leg

### 20k speed
bs0=100
bs=20000

bin_med <- import(paste0(pathdata,root_title,"_speedmed_bin",bs/1000,"k_center.bw"))
bin_CIinf <- import(paste0(pathdata,root_title,"_speedmedCIinf_bin",bs/1000,"k_center.bw"))
bin_CIsup <- import(paste0(pathdata,root_title,"_speedmedCIsup_bin",bs/1000,"k_center.bw"))
bin_shmed <- import(paste0(pathdata,root_title,"_shuffled_medspeedmed_bin",bs/1000,"k_center.bw"))
bin_shCIinf <- import(paste0(pathdata,root_title,"_shuffled_q01speedmed_bin",bs/1000,"k_center.bw"))
bin_shCIsup <- import(paste0(pathdata,root_title,"_shuffled_q99speedmed_bin",bs/1000,"k_center.bw"))
bin_cov <- import(paste0(pathdata,root_title,"_forcov_bin",bs/1000,"k_center.bw"))
medcov <- median(bin_cov$score,na.rm=T)
speedmedgen <- median(bin_med$score,na.rm=T)

pl20k <- lapply(seq_along(chrom), function(i) {
	ROI <- chrom[i]
toplot_bin <- tibble(pos=start(ROI):end(ROI),
speed=as.numeric(unlist(coverage(bin_med,weight=bin_med$score)[ROI])),
speedINF=as.numeric(unlist(coverage(bin_CIinf,weight=bin_CIinf$score)[ROI])),
speedSUP=as.numeric(unlist(coverage(bin_CIsup,weight=bin_CIsup$score)[ROI])),
speed_shuffled=as.numeric(unlist(coverage(bin_shmed,weight=bin_shmed$score)[ROI])),
speedINFsh=as.numeric(unlist(coverage(bin_shCIinf,weight=bin_shCIinf$score)[ROI])),
speedSUPsh=as.numeric(unlist(coverage(bin_shCIsup,weight=bin_shCIsup$score)[ROI])))
toplot_bin[toplot_bin==0] <- NA

toplot_bin2 <- toplot_bin %>% mutate(pos = round(pos/bs0)*bs0+1) %>%
		group_by(pos) %>%
		summarise(speed = mean(speed,na.rm=T), speedINF = mean(speedINF,na.rm=T),speedSUP = mean(speedSUP,na.rm=T),speed_shuffled = mean(speed_shuffled,na.rm=T),speedINFsh = mean(speedINFsh,na.rm=T),speedSUPsh = mean(speedSUPsh,na.rm=T),.groups = "drop")

toplot_cov <- tibble(pos=start(ROI):end(ROI),
coverage=as.numeric(unlist(coverage(bin_cov,weight=bin_cov$score)[ROI])))
toplot_cov[toplot_cov==0] <- NA
toplot_cov2 <- toplot_cov %>% 
		mutate(pos = round(pos/bs0)*bs0+1) %>%
		group_by(pos) %>%
		summarise(coverage = mean(coverage,na.rm=T), .groups = "drop")

p1 <- ggplot(toplot_bin2)+
	geom_line(aes(x=pos,y=speed,col="speed"))+
	geom_ribbon(aes(x=pos,ymin=speedINF,ymax=speedSUP,fill="speed"),alpha=0.3,col=NA)+
	geom_line(aes(x=pos,y=speed_shuffled,col="shuffled"))+
	geom_ribbon(aes(x=pos,ymin=speedINFsh,ymax=speedSUPsh,fill="shuffled"),alpha=0.3,col=NA)+
	geom_hline(yintercept=speedmedgen,linetype=2,alpha=0.5)+
	xlab(paste0(as.character(seqnames(ROI))," (kb)"))+
	ylab("Speed (bp/min)")+
	scale_x_continuous(
		guide = "prism_minor",
		labels=scales::unit_format(big.mark ="",suffix="",scale=1e-3,sep=""),
		limits=c(1,seqlengths(seqinf)[i]),
		breaks=seq(0,seqlengths(seqinf)[i],100000),
		minor_breaks=seq(0,seqlengths(seqinf)[i],50000),
		expand=c(0,0))+
	scale_color_manual("",values = speed_pal)+
	scale_fill_manual("",values = speed_pal)+
	theme(legend.key = element_rect(colour = "black"))+
	coord_cartesian(ylim=c(500,3500),expand=F)
p1cov <- ggplot(toplot_cov2)+
	geom_hline(yintercept=medcov,linetype=2,alpha=0.5)+
	geom_line(aes(x=pos,y=coverage),col=speed_pal[1])+
	scale_x_continuous(
		guide = "prism_minor",
		labels=scales::unit_format(big.mark ="",suffix="",scale=1e-3,sep=""),
		limits=c(1,seqlengths(seqinf)[i]),
		breaks=seq(0,seqlengths(seqinf)[i],100000),
		minor_breaks=seq(0,seqlengths(seqinf)[i],50000),
		expand=c(0,0))+
	xlab(paste0(as.character(seqnames(ROI))," (kb)"))+
	coord_cartesian(ylim=c(0,500),expand=F)

return(list(p1,p1cov))
})



###MWW
wtl.1k <- import(paste0(pathdata,root_title,"_wl_adj_bin1k_center.bw"))
wtl.2k <- import(paste0(pathdata,root_title,"_wl_adj_bin2k_center.bw"))
wtl.3k <- import(paste0(pathdata,root_title,"_wl_adj_bin3k_center.bw"))
wtl.4k <- import(paste0(pathdata,root_title,"_wl_adj_bin4k_center.bw"))
wtl.5k <- import(paste0(pathdata,root_title,"_wl_adj_bin5k_center.bw"))
wtl.6k <- import(paste0(pathdata,root_title,"_wl_adj_bin6k_center.bw"))
wtl.7k <- import(paste0(pathdata,root_title,"_wl_adj_bin7k_center.bw"))
wtl.8k <- import(paste0(pathdata,root_title,"_wl_adj_bin8k_center.bw"))
wtl.9k <- import(paste0(pathdata,root_title,"_wl_adj_bin9k_center.bw"))
wtl.10k <- import(paste0(pathdata,root_title,"_wl_adj_bin10k_center.bw"))
wtl.15k <- import(paste0(pathdata,root_title,"_wl_adj_bin15k_center.bw"))
wtl.20k <- import(paste0(pathdata,root_title,"_wl_adj_bin20k_center.bw"))

wtg.1k <- import(paste0(pathdata,root_title,"_wg_adj_bin1k_center.bw"))
wtg.2k <- import(paste0(pathdata,root_title,"_wg_adj_bin2k_center.bw"))
wtg.3k <- import(paste0(pathdata,root_title,"_wg_adj_bin3k_center.bw"))
wtg.4k <- import(paste0(pathdata,root_title,"_wg_adj_bin4k_center.bw"))
wtg.5k <- import(paste0(pathdata,root_title,"_wg_adj_bin5k_center.bw"))
wtg.6k <- import(paste0(pathdata,root_title,"_wg_adj_bin6k_center.bw"))
wtg.7k <- import(paste0(pathdata,root_title,"_wg_adj_bin7k_center.bw"))
wtg.8k <- import(paste0(pathdata,root_title,"_wg_adj_bin8k_center.bw"))
wtg.9k <- import(paste0(pathdata,root_title,"_wg_adj_bin9k_center.bw"))
wtg.10k <- import(paste0(pathdata,root_title,"_wg_adj_bin10k_center.bw"))
wtg.15k <- import(paste0(pathdata,root_title,"_wg_adj_bin15k_center.bw"))
wtg.20k <- import(paste0(pathdata,root_title,"_wg_adj_bin20k_center.bw"))

bs0 <- 1000

pl_MWW3 <- lapply(seq_along(chrom), function(i) {
	ROI <- chrom[i]
toplot_wtl <- bind_rows(
tibble(pos=start(ROI):end(ROI),wt=as.numeric(unlist(coverage(wtl.1k,weight=wtl.1k$score)[ROI])),type="slow1k",type2="1k"),
tibble(pos=start(ROI):end(ROI),wt=as.numeric(unlist(coverage(wtl.2k,weight=wtl.2k$score)[ROI])),type="slow2k",type2="2k"),
tibble(pos=start(ROI):end(ROI),wt=as.numeric(unlist(coverage(wtl.3k,weight=wtl.3k$score)[ROI])),type="slow3k",type2="3k"),
tibble(pos=start(ROI):end(ROI),wt=as.numeric(unlist(coverage(wtl.4k,weight=wtl.4k$score)[ROI])),type="slow4k",type2="4k"),
tibble(pos=start(ROI):end(ROI),wt=as.numeric(unlist(coverage(wtl.5k,weight=wtl.5k$score)[ROI])),type="slow5k",type2="5k"),
tibble(pos=start(ROI):end(ROI),wt=as.numeric(unlist(coverage(wtl.6k,weight=wtl.6k$score)[ROI])),type="slow6k",type2="6k"),
tibble(pos=start(ROI):end(ROI),wt=as.numeric(unlist(coverage(wtl.7k,weight=wtl.7k$score)[ROI])),type="slow7k",type2="7k"),
tibble(pos=start(ROI):end(ROI),wt=as.numeric(unlist(coverage(wtl.8k,weight=wtl.8k$score)[ROI])),type="slow8k",type2="8k"),
tibble(pos=start(ROI):end(ROI),wt=as.numeric(unlist(coverage(wtl.9k,weight=wtl.9k$score)[ROI])),type="slow9k",type2="9k"),
tibble(pos=start(ROI):end(ROI),wt=as.numeric(unlist(coverage(wtl.10k,weight=wtl.10k$score)[ROI])),type="slow10k",type2="10k"),
tibble(pos=start(ROI):end(ROI),wt=as.numeric(unlist(coverage(wtl.15k,weight=wtl.15k$score)[ROI])),type="slow15k",type2="15k"),
tibble(pos=start(ROI):end(ROI),wt=as.numeric(unlist(coverage(wtl.20k,weight=wtl.20k$score)[ROI])),type="slow20k",type2="20k")
)
toplot_wtl$type2 <- factor(toplot_wtl$type2, levels=c("1k","2k","3k","4k","5k","6k","7k","8k","9k","10k","15k","20k"))

toplot_wtl2 <- toplot_wtl %>%
		mutate(pos = floor((pos-1)/bs0)*bs0+1) %>%
		group_by(type2,pos) %>%
		summarise(wt = mean(wt,na.rm=T), .groups = "drop")


toplot_wtg <- bind_rows(
tibble(pos=start(ROI):end(ROI),wt=as.numeric(unlist(coverage(wtg.1k,weight=wtg.1k$score)[ROI])),type="fast1k",type2="1k"),
tibble(pos=start(ROI):end(ROI),wt=as.numeric(unlist(coverage(wtg.2k,weight=wtg.2k$score)[ROI])),type="fast2k",type2="2k"),
tibble(pos=start(ROI):end(ROI),wt=as.numeric(unlist(coverage(wtg.3k,weight=wtg.3k$score)[ROI])),type="fast3k",type2="3k"),
tibble(pos=start(ROI):end(ROI),wt=as.numeric(unlist(coverage(wtg.4k,weight=wtg.4k$score)[ROI])),type="fast4k",type2="4k"),
tibble(pos=start(ROI):end(ROI),wt=as.numeric(unlist(coverage(wtg.5k,weight=wtg.5k$score)[ROI])),type="fast5k",type2="5k"),
tibble(pos=start(ROI):end(ROI),wt=as.numeric(unlist(coverage(wtg.6k,weight=wtg.6k$score)[ROI])),type="fast6k",type2="6k"),
tibble(pos=start(ROI):end(ROI),wt=as.numeric(unlist(coverage(wtg.7k,weight=wtg.7k$score)[ROI])),type="fast7k",type2="7k"),
tibble(pos=start(ROI):end(ROI),wt=as.numeric(unlist(coverage(wtg.8k,weight=wtg.8k$score)[ROI])),type="fast8k",type2="8k"),
tibble(pos=start(ROI):end(ROI),wt=as.numeric(unlist(coverage(wtg.9k,weight=wtg.9k$score)[ROI])),type="fast9k",type2="9k"),
tibble(pos=start(ROI):end(ROI),wt=as.numeric(unlist(coverage(wtg.10k,weight=wtg.10k$score)[ROI])),type="fast10k",type2="10k"),
tibble(pos=start(ROI):end(ROI),wt=as.numeric(unlist(coverage(wtg.15k,weight=wtg.15k$score)[ROI])),type="fast15k",type2="15k"),
tibble(pos=start(ROI):end(ROI),wt=as.numeric(unlist(coverage(wtg.20k,weight=wtg.20k$score)[ROI])),type="fast20k",type2="20k")
)
toplot_wtg$type2 <- factor(toplot_wtg$type2, levels=c("1k","2k","3k","4k","5k","6k","7k","8k","9k","10k","15k","20k"))

toplot_wtg2 <- toplot_wtg %>%
		mutate(pos = floor((pos-1)/bs0)*bs0+1) %>%
		group_by(type2,pos) %>%
		summarise(wt = mean(wt,na.rm=T), .groups = "drop")
		
## group and fix test threshold to 1e-2
MWW.th <- 1e-2
toplot_wt2 <- full_join(toplot_wtg2,toplot_wtl2,by=c("type2","pos"),suffix=c(".f",".s")) %>%
		mutate(wt2=map2_chr(wt.f,wt.s, function(x,y) case_when((x<MWW.th)~"faster",(y<MWW.th)~"slower",(x>=MWW.th & y>=MWW.th)~"n.s.",T~"NA")) %>% as_factor)

custom_breaks <- seq(start(ROI)-1, end(ROI), 10000)

p2 <- ggplot(toplot_wt2)+geom_tile(aes(x=pos+500,y=type2,fill=wt2))+
	scale_fill_manual("",values=c("purple","green3","white","grey65"),breaks = c("slower", "faster", "n.s.", "NA"),labels = c("slower", "faster", "n.s.", "N/A"),na.value = "grey65")+
	xlab(paste0(as.character(seqnames(ROI))," (kb)"))+
	ylab("MWW.test.adj")+
	coord_cartesian(expand=F)+
	scale_x_continuous(
		guide = "prism_minor",
		labels=scales::unit_format(big.mark ="",suffix="",scale=1e-3,sep=""),
		limits=c(1,seqlengths(seqinf)[i]),
		breaks=seq(0,seqlengths(seqinf)[i],100000),
		minor_breaks=seq(0,seqlengths(seqinf)[i],50000),
		expand=c(0,0))+
	theme(axis.ticks.y = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.key = element_rect(colour = "black"))

return(p2)
})


### creat geno_track, edit it for i=12 (rDNA), i=3 (HMLR)
feat.list <- GRangesList(
	import(paste0(pathdata,"sc3_CEN.bed")),
	import(paste0(pathdata,"sc3_TRNA.bed")),
	import(paste0(pathdata,"sc3_RRNA_locus.bed")),
	import(paste0(pathdata,"sc3_TEL.bed")),
	import(paste0(pathdata,"sc3_HMLR.bed")),
	import(paste0(pathdata,"ARS_newman.bed"))
)
names(feat.list) <- c("CEN","tRNA","rDNA","TEL","HML/HMR","ORI")
#feat.list <- endoapply(feat.list,NewSeqinfo,seqin=seqinf)
feat.list <- lapply(seq_along(feat.list), function(x) {feat.list[[x]]$type=names(feat.list)[x];return(feat.list[[x]])})
feat <- do.call(c,feat.list)

### genomic track chrXII
geno_leg <- factor(c("CEN","tRNA","rDNA","ORI","TEL"),levels=c("CEN","tRNA","rDNA","ORI","TEL"))
geno_pal <- mypal[c(5,3,9,13,19)]
names(geno_pal) <- geno_leg
i=12
	ROI <- chrom[i]
	featROI <- as_tibble(feat[overlapsAny(feat,ROI)]) %>% rename(featname=name)
pl_geno_12 <- ggplot(featROI)+
	geom_rect(data=featROI %>% filter(type %in% c("rDNA","tRNA")),aes(xmin=start,xmax=end,ymin=2,ymax=3,col=type,fill=type))+
	geom_rect(data=featROI%>% filter(type %in% c("CEN","TEL")),aes(xmin=start,xmax=end,ymin=1,ymax=2,col=type,fill=type))+
	geom_rect(data=featROI%>% filter(type %in% c("ORI")),aes(xmin=start,xmax=end,ymin=0,ymax=1,col=type,fill=type))+
	scale_color_manual("",values = geno_pal)+
	scale_fill_manual("",values = geno_pal)+
	xlab(paste0(as.character(seqnames(ROI))," (kb)"))+
	scale_x_continuous(
		guide = "prism_minor",
		labels=scales::unit_format(big.mark ="",suffix="",scale=1e-3,sep=""),
		limits=c(1,seqlengths(seqinf)[i]),
		breaks=seq(0,seqlengths(seqinf)[i],100000),
		minor_breaks=seq(0,seqlengths(seqinf)[i],50000),
		expand=c(0,0))+
	theme(axis.ticks.y = element_blank(),axis.text.y = element_blank(),panel.grid.major.y=element_blank(),panel.grid.minor.y=element_blank(),legend.key = element_rect(colour = "black"))

### genomic track chrIII
geno_leg <- factor(c("CEN","tRNA","HML/HMR","ORI","TEL"),levels=c("CEN","tRNA","HML/HMR","ORI","TEL"))
geno_pal <- mypal[c(5,3,17,13,19)]
names(geno_pal) <- geno_leg
i=3
	ROI <- chrom[i]
	featROI <- as_tibble(feat[overlapsAny(feat,ROI)]) %>% rename(featname=name)
pl_geno_3 <-ggplot(featROI)+
	geom_rect(data=featROI %>% filter(type %in% c("tRNA")),aes(xmin=start,xmax=end,ymin=2,ymax=3,col=type,fill=type))+
	geom_rect(data=featROI%>% filter(type %in% c("CEN","TEL","HML/HMR")),aes(xmin=start,xmax=end,ymin=1,ymax=2,col=type,fill=type))+
	geom_rect(data=featROI%>% filter(type %in% c("ORI")),aes(xmin=start,xmax=end,ymin=0,ymax=1,col=type,fill=type))+
	scale_color_manual("",values = geno_pal)+
	scale_fill_manual("",values = geno_pal)+
	xlab(paste0(as.character(seqnames(ROI))," (kb)"))+
	scale_x_continuous(
		guide = "prism_minor",
		labels=scales::unit_format(big.mark ="",suffix="",scale=1e-3,sep=""),
		limits=c(1,seqlengths(seqinf)[i]),
		breaks=seq(0,seqlengths(seqinf)[i],100000),
		minor_breaks=seq(0,seqlengths(seqinf)[i],50000),
		expand=c(0,0))+
	theme(axis.ticks.y = element_blank(),axis.text.y = element_blank(),panel.grid.major.y=element_blank(),panel.grid.minor.y=element_blank(),legend.key = element_rect(colour = "black"))

### genomic track chrOTHER
geno_leg <- factor(c("CEN","tRNA","ORI","TEL"),levels=c("CEN","tRNA","ORI","TEL"))
geno_pal <- mypal[c(5,3,13,19)]
names(geno_pal) <- geno_leg

pl_geno_other <- lapply(seq_along(chrom), function(i) {
	ROI <- chrom[i]
	featROI <- as_tibble(feat[overlapsAny(feat,ROI)]) %>% rename(featname=name)
ggplot(featROI)+
	geom_rect(data=featROI %>% filter(type %in% c("tRNA")),aes(xmin=start,xmax=end,ymin=2,ymax=3,col=type,fill=type))+
	geom_rect(data=featROI%>% filter(type %in% c("CEN","TEL")),aes(xmin=start,xmax=end,ymin=1,ymax=2,col=type,fill=type))+
	geom_rect(data=featROI%>% filter(type %in% c("ORI")),aes(xmin=start,xmax=end,ymin=0,ymax=1,col=type,fill=type))+
	scale_color_manual("",values = geno_pal)+
	scale_fill_manual("",values = geno_pal)+
	xlab(paste0(as.character(seqnames(ROI))," (kb)"))+
	scale_x_continuous(
		guide = "prism_minor",
		labels=scales::unit_format(big.mark ="",suffix="",scale=1e-3,sep=""),
		limits=c(1,seqlengths(seqinf)[i]),
		breaks=seq(0,seqlengths(seqinf)[i],100000),
		minor_breaks=seq(0,seqlengths(seqinf)[i],50000),
		expand=c(0,0))+
	theme(axis.ticks.y = element_blank(),axis.text.y = element_blank(),panel.grid.major.y=element_blank(),panel.grid.minor.y=element_blank(),legend.key = element_rect(colour = "black"))
})

###
pp <- list()
for (i in 1:16) 
{
pl_speed <- pl20k[[i]][[1]]+theme(axis.title.x=element_blank(),legend.key.size = unit(0.3, 'cm'))+coord_cartesian(ylim=c(900,2600))+labs(tag=names(chrom)[i])+theme(plot.tag=element_text(face="bold"))
pl_gen <- pl_geno_other[[i]]+theme(legend.title=element_blank(),axis.title.x=element_blank(),axis.ticks.x = element_blank(),axis.text.x = element_blank(),legend.key.size = unit(0.3, 'cm'))+ylab("Features")
pl_test <- pl_MWW3[[i]]+theme(legend.title=element_blank(),axis.title.x=element_blank(),legend.key.size = unit(0.3, 'cm'))+ylab("Scale")
pl_cov <- pl20k[[i]][[2]]+ylab("Coverage")
pp[[i]] <- pl_speed/pl_test/pl_gen/pl_cov+plot_layout(ncol = 1, heights = c(3.8,1.4,0.7,0.7))
}
i=3
pl_speed <- pl20k[[i]][[1]]+theme(axis.title.x=element_blank(),legend.key.size = unit(0.3, 'cm'))+coord_cartesian(ylim=c(900,2600))+labs(tag=names(chrom)[i])+theme(plot.tag=element_text(face="bold"))
pl_gen <- pl_geno_3+theme(legend.title=element_blank(),axis.title.x=element_blank(),axis.ticks.x = element_blank(),axis.text.x = element_blank(),legend.key.size = unit(0.3, 'cm'))+ylab("Features")
pl_test <- pl_MWW3[[i]]+theme(legend.title=element_blank(),axis.title.x=element_blank(),legend.key.size = unit(0.3, 'cm'))+ylab("Scale")
pl_cov <- pl20k[[i]][[2]]+ylab("Coverage")
pp[[i]] <- pl_speed/pl_test/pl_gen/pl_cov+plot_layout(ncol = 1, heights = c(3.8,1.4,0.7,0.7))

i=12
pl_speed <- pl20k[[i]][[1]]+theme(axis.title.x=element_blank(),legend.key.size = unit(0.3, 'cm'))+coord_cartesian(ylim=c(900,2600))+labs(tag=names(chrom)[i])+theme(plot.tag=element_text(face="bold"))
pl_gen <- pl_geno_12+theme(legend.title=element_blank(),axis.title.x=element_blank(),axis.ticks.x = element_blank(),axis.text.x = element_blank(),legend.key.size = unit(0.3, 'cm'))+ylab("Features")
pl_test <- pl_MWW3[[i]]+theme(legend.title=element_blank(),axis.title.x=element_blank(),legend.key.size = unit(0.3, 'cm'))+ylab("Scale")
pl_cov <- pl20k[[i]][[2]]+ylab("Coverage")
pp[[i]] <- pl_speed/pl_test/pl_gen/pl_cov+plot_layout(ncol = 1, heights = c(3.8,1.4,0.7,0.7))

f5 <- (pp[[11]]|pp[[12]])/(pp[[13]]|pp[[14]]) & theme(legend.position="right")
ggsave(paste0(path_figures,"Figure5.pdf"),h=14,w=14,plot=f5)

fs8 <- (pp[[1]]|pp[[2]])/(pp[[3]]|pp[[4]]) & theme(legend.position="right")
ggsave(paste0(path_figures,"FigureS8_1.pdf"),h=14,w=14,plot=fs8)
fs8 <- (pp[[5]]|pp[[6]])/(pp[[7]]|pp[[8]]) & theme(legend.position="right")
ggsave(paste0(path_figures,"FigureS8_2.pdf"),h=14,w=14,plot=fs8)
fs8 <- (pp[[9]]|pp[[10]])/(pp[[15]]|pp[[16]]) & theme(legend.position="right")
ggsave(paste0(path_figures,"FigureS8_3.pdf"),h=14,w=14,plot=fs8)


### localized zoom for chrXI
### 1k speed
bs0=10
bs=1000

bin_med <- import(paste0(pathdata,root_title,"_speedmed_bin",bs/1000,"k_center.bw"))
bin_CIinf <- import(paste0(pathdata,root_title,"_speedmedCIinf_bin",bs/1000,"k_center.bw"))
bin_CIsup <- import(paste0(pathdata,root_title,"_speedmedCIsup_bin",bs/1000,"k_center.bw"))
bin_shmed <- import(paste0(pathdata,root_title,"_shuffled_medspeedmed_bin",bs/1000,"k_center.bw"))
bin_shCIinf <- import(paste0(pathdata,root_title,"_shuffled_q01speedmed_bin",bs/1000,"k_center.bw"))
bin_shCIsup <- import(paste0(pathdata,root_title,"_shuffled_q99speedmed_bin",bs/1000,"k_center.bw"))
bin_cov <- import(paste0(pathdata,root_title,"_forcov_bin",bs/1000,"k_center.bw"))
medcov <- median(bin_cov$score,na.rm=T)
pl1k <- lapply(seq_along(chrom), function(i) {
	ROI <- chrom[i]
toplot_bin <- tibble(pos=start(ROI):end(ROI),
speed=as.numeric(unlist(coverage(bin_med,weight=bin_med$score)[ROI])),
speedINF=as.numeric(unlist(coverage(bin_CIinf,weight=bin_CIinf$score)[ROI])),
speedSUP=as.numeric(unlist(coverage(bin_CIsup,weight=bin_CIsup$score)[ROI])),
speed_shuffled=as.numeric(unlist(coverage(bin_shmed,weight=bin_shmed$score)[ROI])),
speedINFsh=as.numeric(unlist(coverage(bin_shCIinf,weight=bin_shCIinf$score)[ROI])),
speedSUPsh=as.numeric(unlist(coverage(bin_shCIsup,weight=bin_shCIsup$score)[ROI])))
toplot_bin[toplot_bin==0] <- NA

toplot_bin2 <- toplot_bin %>% mutate(pos = round(pos/bs0)*bs0+1) %>%
		group_by(pos) %>%
		summarise(speed = mean(speed,na.rm=T), speedINF = mean(speedINF,na.rm=T),speedSUP = mean(speedSUP,na.rm=T),speed_shuffled = mean(speed_shuffled,na.rm=T),speedINFsh = mean(speedINFsh,na.rm=T),speedSUPsh = mean(speedSUPsh,na.rm=T),.groups = "drop")

toplot_cov <- tibble(pos=start(ROI):end(ROI),
coverage=as.numeric(unlist(coverage(bin_cov,weight=bin_cov$score)[ROI])))
toplot_cov[toplot_cov==0] <- NA
toplot_cov2 <- toplot_cov %>% 
		mutate(pos = (round(pos/bs0)+1)*bs0) %>%
		group_by(pos) %>%
		summarise(coverage = mean(coverage,na.rm=T), .groups = "drop")

p1 <- ggplot(toplot_bin2)+
	geom_line(aes(x=pos,y=speed,col="speed"))+
	geom_ribbon(aes(x=pos,ymin=speedINF,ymax=speedSUP,fill="speed"),alpha=0.3,col=NA)+
	geom_line(aes(x=pos,y=speed_shuffled,col="shuffled"))+
	geom_ribbon(aes(x=pos,ymin=speedINFsh,ymax=speedSUPsh,fill="shuffled"),alpha=0.3,col=NA)+
	geom_hline(yintercept=speedmedgen,linetype=2,alpha=0.5)+
	xlab(paste0(as.character(seqnames(ROI))," (kb)"))+
	ylab("Speed (bp/min)")+
	scale_x_continuous(
		guide = "prism_minor",
		labels=scales::unit_format(big.mark ="",suffix="",scale=1e-3,sep=""),
		limits=c(1,seqlengths(seqinf)[i]),
		breaks=seq(0,seqlengths(seqinf)[i],100000),
		minor_breaks=seq(0,seqlengths(seqinf)[i],50000),
		expand=c(0,0))+
	scale_color_manual("",values = speed_pal)+
	scale_fill_manual("",values = speed_pal)+
	theme(legend.key = element_rect(colour = "black"))+
	coord_cartesian(ylim=c(500,3500),expand=F)
p1cov <- ggplot(toplot_cov2)+
	geom_hline(yintercept=medcov,linetype=2,alpha=0.5)+
	geom_line(aes(x=pos,y=coverage),col=speed_pal[1])+
	scale_x_continuous(
		guide = "prism_minor",
		labels=scales::unit_format(big.mark ="",suffix="",scale=1e-3,sep=""),
		limits=c(1,seqlengths(seqinf)[i]),
		breaks=seq(0,seqlengths(seqinf)[i],100000),
		minor_breaks=seq(0,seqlengths(seqinf)[i],50000),
		expand=c(0,0))+
	xlab(paste0(as.character(seqnames(ROI))," (kb)"))+
	coord_cartesian(ylim=c(0,500),expand=F)

return(list(p1,p1cov))
})



library(ggrepel)
ROI_z <- GRanges(seqnames="chrXI",range=IRanges(645000,665000),strand="*",seqinfo=seqinf)

genes_tb <- as_tibble(import(paste0(pathdata,"Saccharomyces_cerevisiae.R64-1-1.104.gtf"))) %>% filter(seqnames !="Mito") %>% mutate(chrom="chr" %+% seqnames) %>% filter(type=="gene") %>% select(chrom,start,end,strand,gene_id,gene_name)
# merging gene_ID,common_name
genes_tb <- genes_tb %>% mutate(name=ifelse(is.na(gene_name),gene_id,gene_name))
gene_GR <- with(genes_tb, GRanges(seqnames=chrom,range=IRanges(start,end),strand=strand,name=name,seqinfo=seqinf))

geno_leg2 <- factor(c("geneR","geneL"),levels=c("geneR","geneL"))
geno_pal2 <- mypal[c(15,17)]
names(geno_pal2) <- geno_leg2



pl_geno2 <- lapply(seq_along(chrom), function(i) {
	ROI <- chrom[i]
	featROI <- as_tibble(feat[overlapsAny(feat,ROI)]) %>% rename(featname=name)
	geneROIR <- as_tibble(gene_GR[overlapsAny(gene_GR,ROI) & strand(gene_GR)=="+"]) %>% mutate(type="geneR")
	geneROIL <- as_tibble(gene_GR[overlapsAny(gene_GR,ROI) & strand(gene_GR)=="-"]) %>% mutate(type="geneL")
ggplot(featROI)+
	geom_rect(data=geneROIR,aes(xmin=start,xmax=end,ymin=1,ymax=2,col=type,fill=type))+
	geom_label_repel(data=geneROIR,aes(x=(start+end)/2,y=1.5,label=name,col=type),size=2,max_iter=10000,show.legend=F,force=10)+
	geom_rect(data=geneROIL,aes(xmin=start,xmax=end,ymin=0,ymax=1,col=type,fill=type))+
	geom_label_repel(data=geneROIL,aes(x=(start+end)/2,y=0.5,label=name,col=type),size=2,max_iter=10000,show.legend=F,force=10)+
	scale_color_manual("",values = geno_pal2)+
	scale_fill_manual("",values = geno_pal2)+
	xlab(paste0(as.character(seqnames(ROI))," (kb)"))+
	scale_x_continuous(
		guide = "prism_minor",
		labels=scales::unit_format(big.mark ="",suffix="",scale=1e-3,sep=""),
		limits=c(1,seqlengths(seqinf)[i]),
		breaks=seq(0,seqlengths(seqinf)[i],100000),
		minor_breaks=seq(0,seqlengths(seqinf)[i],50000),
		expand=c(0,0))+
	theme(axis.ticks.y = element_blank(),axis.text.y = element_blank(),panel.grid.major.y=element_blank(),panel.grid.minor.y=element_blank(),legend.key = element_rect(colour = "black"))
})

	
# zoom
roi <- as_tibble(ROI_z)
ch <- roi$seqnames
pl_speed <- pl1k[[roi$seqnames]][[1]]+
	scale_x_continuous(labels=scales::unit_format(big.mark ="",suffix="",scale=1e-3,sep="",accuracy=1),limits=c(roi$start,min(seqlengths(seqinf)[roi$seqnames],roi$end)))+
	theme(axis.title.x=element_blank(),legend.key.size = unit(0.3, 'cm'))
pl_gen <- pl_geno2[[roi$seqnames]]+
	scale_x_continuous(labels=scales::unit_format(big.mark ="",suffix="",scale=1e-3,sep="",accuracy=1),limits=c(roi$start,min(seqlengths(seqinf)[roi$seqnames],roi$end)),expand=c(0,0))+
	theme(legend.title=element_blank(),axis.title.x=element_blank(),axis.ticks.x = element_blank(),axis.text.x = element_blank(),legend.key.size = unit(0.3, 'cm'))+
	ylab("Genes")
pl_test <- pl_MWW3[[roi$seqnames]]+
	scale_x_continuous(labels=scales::unit_format(big.mark ="",suffix="",scale=1e-3,sep="",accuracy=1),limits=c(roi$start,min(seqlengths(seqinf)[roi$seqnames],roi$end)))+
	theme(legend.title=element_blank(),axis.title.x=element_blank(),legend.key.size = unit(0.3, 'cm'))+ylab("Scale")
pl_cov <- pl1k[[roi$seqnames]][[2]]+
	scale_x_continuous(labels=scales::unit_format(big.mark ="",suffix="",scale=1e-3,sep="",accuracy=1),limits=c(roi$start,min(seqlengths(seqinf)[roi$seqnames],roi$end)))+
	coord_cartesian(ylim=c(0,50),expand=F)+
	ylab("Coverage")
	
pl_speed/pl_test/pl_gen/pl_cov+plot_layout(ncol = 1, heights = c(2,1,0.5,0.5))
ggsave(paste0(path_figures,"FigureS9_.pdf"),h=7,w=7)
