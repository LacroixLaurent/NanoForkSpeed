### Yeast Feature annotation from UCSC track
### Scfeat downloaded from https://genome.ucsc.edu/cgi-bin/hgTables from the sgdOther tracks
library(tidyverse)
pathdata <- "/Users/ll/work/RStudioProjects/NanoForkSpeed/Figures_Article/data/"
allfeat <- read_tsv(paste0(pathdata,"ScFeat.bed"),col_names=c("seqnames","start","end","name","score","strand"))

tel <- allfeat %>% filter(str_detect(name,"TEL.{2}[R|L]$"))

cen <- allfeat %>% filter(str_detect(name,"CEN")) %>% group_by(name) %>% mutate(len=end-start) %>% filter(len==max(len)) %>% ungroup %>% select(-len)

trna <- allfeat %>% filter(seqnames!="chrM") %>% filter(str_detect(name,"^t.(...)")) %>% group_by(name) %>% mutate(len=end-start) %>% filter(len==max(len)) %>% ungroup %>% select(-len) %>% distinct

hmlr <- allfeat %>% filter(name %in% c("HML","HMR"))

rrna <- allfeat %>% filter(seqnames!="chrM") %>% filter(str_detect(name,"^RDN") | str_detect(name,"^ETS")| str_detect(name,"^ITS")| str_detect(name,"^NTS")) %>% group_by(name) %>% mutate(len=end-start) %>% filter(len==max(len)) %>% ungroup %>% select(-len) %>% distinct

library(rtracklayer)
seqinf <- readRDS(paste0(pathdata,"seqinfS288CrDNA.rds"))

TEL <- sort(with(tel,GRanges(seqnames=seqnames,ranges=IRanges(start,end),strand=strand,name=name,seqinfo=seqinf)),ignore.strand=T)
CEN <- sort(with(cen,GRanges(seqnames=seqnames,ranges=IRanges(start,end),strand=strand,name=name,seqinfo=seqinf)),ignore.strand=T)
HMLR <- sort(with(hmlr,GRanges(seqnames=seqnames,ranges=IRanges(start,end),strand=strand,name=name,seqinfo=seqinf)),ignore.strand=T)
TRNA <- sort(with(trna,GRanges(seqnames=seqnames,ranges=IRanges(start,end),strand=strand,name=name,seqinfo=seqinf)),ignore.strand=T)
RRNA <- sort(with(rrna,GRanges(seqnames=seqnames,ranges=IRanges(start,end),strand=strand,name=name,seqinfo=seqinf)),ignore.strand=T)
## to isolate only the 2 complete repeats
RRNA_locus <- reduce(RRNA,ignore.strand=T)[1]
export(RRNA_locus,con=paste0(pathdata,"sc3_RRNA_locus.bed"))
export(TEL,con=paste0(pathdata,"sc3_TEL.bed"))
export(CEN,con=paste0(pathdata,"sc3_CEN.bed"))
export(HMLR,con=paste0(pathdata,"sc3_HMLR.bed"))
export(TRNA,con=paste0(pathdata,"sc3_TRNA.bed"))

# keep only RDN37-1,2 and RDN5-1,2 and add RFB
RRNA2 <- RRNA[RRNA$name %in% c("RDN37-1","RDN37-2","RDN5-1","RDN5-2")]
RRNA3 <- c(GRanges("chrXII",IRanges(460504,460684),strand="*",seqinfo=seqinf,name="RFB"),RRNA2)
export(RRNA3,con=paste0(pathdata,"sc3_RRNA.bed"))

