# script figure S11

library(GenomicRanges)
library(rtracklayer)
suppressMessages(library(tidyverse))
theme_set(theme_bw())
library(patchwork)
library(ggprism)

mypal <- paletteer::paletteer_d("ggthemes::Classic_20")
`%+%` <- paste0

path_figures <- "./Figures_Article/figures/"
pathdata <- "./Figures_Article/data/"

seqinf <- readRDS(paste0(pathdata,"seqinfS288CrDNA.rds"))

ini2plot <- read_tsv(paste0(pathdata,"FigureS11_data.tsv.gz"))

ARS <- import(paste0(pathdata,"ARS_newman.bed"))
rDNA <- import(paste0(pathdata,"sc3_RRNA.bed"))
ars <- as_tibble(ARS) %>% mutate(chrom=seqnames)
rdna <- as_tibble(rDNA) %>% mutate(chrom=seqnames)

geno_leg <- c("Ini. segment","ORI","RFB","rDNA")
geno_pal <- mypal[c(2,13,5,9)]
names(geno_pal) <- geno_leg

p1 <- ggplot(ini2plot)+
	geom_segment(aes(x=x0,xend=x1,y=width,yend=width,col="Ini. segment"),alpha=0.6)+
	geom_segment(data=ars %>% filter(chrom=="chrXII"),aes(x=start,xend=end,y=0,yend=40000,col="ORI"),size=0.5)+
	geom_rect(data=rdna %>% dplyr::slice(1),aes(xmin=start,xmax=end,ymin=200,ymax=39800,col="RFB"),size=0.5)+
	geom_rect(data=rdna %>% dplyr::slice(2:5),aes(xmin=start,xmax=end,ymin=-500,ymax=-1000,fill="rDNA"),size=1,col=NA)+
	xlab("chrXII (kb)")+
	ylab("Width (bp)")+
	scale_x_continuous(guide = "prism_minor",
		labels=scales::unit_format(big.mark ="",suffix="",scale=1e-3,sep=""),limits=c(1,seqlengths(seqinf)[12]),
		breaks=seq(0,seqlengths(seqinf)[12],20000),
		minor_breaks=seq(0,seqlengths(seqinf)[12],5000),
		expand=c(0,0))+
	ylim(c(-1000,40000))+
	coord_cartesian(xlim=c(350000,550000))+
	scale_color_manual("",values = geno_pal)+
	scale_fill_manual("",values = geno_pal)+
	theme(panel.grid.minor.y=element_blank(),panel.grid.minor.x=element_blank(),panel.grid.major.y=element_blank(),panel.grid.major.x=element_blank(),legend.key.size = unit(0.3, 'cm'))

ggsave(paste0(path_figures,"FigureS11.pdf"),w=8,h=2,p1)
