PLS data generation
================

This is a typical example how we proceed to detect forks, initiations
and terminations from the reads base\_called with megalodon and our
model (script=basecalling\_sample.sh) then parsed with the
Parsing\_function4Megalodon.r (script=basecalling\_sample.r).

### Background threshold determination

The first step of the analysis requires to check the distribution of the
signal in order to set the b2a threshold, used to identify the
background of the signal.  
This script outputs a PDF file corresponding to the distribution of the
BrdU signal in 1kb bins for the reads filtered the presence of 3 and
more linearized segments.

``` r
suppressMessages(library(tidyverse))
`%+%`<- paste0
source("./helper_function.r")
source("./PLS_function.r")
expmeg <- "JP3A_Megalodon_00_smdata.rds"
pathdata <- "~/work/Ori/newRaw_nt/"
ex.name <- expmeg %>% str_remove("_Megalodon_00_smdata.rds")
EXP <- readRDS(pathdata %+% expmeg)  %>% filter(chrom!="chrM")
plot_signal(EXP,EXPname=ex.name,EXP_b2a.thr0=0.02,alldata=F,nreads=5000,saved=T,plotit=T)
```

It is also possible to output the distribution for all the reads (of at
least 5kb) but usually this distribution is less informative due to the
overwhelming presence of reads with very low signal. From the second
plot of the PDF file (a zoom of the first plot of th e PDF), it is
possible to draw a line between the 0 level signal (background) and the
informative signal which usually come in two wide distribution. One of
medium level mode corresponding to the DNA replicated during the chase
and one with higher level mode corresponding to the DNA replicated
during the pulse.  
The position of this line was the threshold we fixed for our further
analysis. This threshold was quite similar between all the experiments
we analyzed, thus a fixed threshold of 0.02 was chosen to simplify the
pipeline.  
We also limit to 5000 the number of reads used for this evaluation to
speed up the computation but all the reads can be used by setting
nreads=NA.

### PLS analysis

Once the user has checked that the distribution of the signal and the
position of the background threshold (ex.b2a in the following script)
are satisfactory, the proper forks detection can be performed.  
This procedure create an output file saved in the .rds format. The data
are organised as a list or 4 elements:  
1- a list of tibble containing:  
1.1- the reads filtered for their size and the presence of 3 and more
linearized segments.  
1.2- the reads with detected forks.  
2- a tibble of all the detected forks.  
3- a tibble containing the initiations (Init) and the terminations (Ter)
detected.  
4- a table summary of different metrics of the experiment.

``` r
library(GenomicRanges)
library(rtracklayer)
suppressMessages(library(kmlShape))
suppressMessages(library(tidyverse))
library(formattable)
`%+%`<- paste0
source("PLS_function.r")

seqinf <- Seqinfo(seqnames=c("chrI","chrII","chrIII","chrIV","chrV","chrVI","chrVII","chrVIII","chrIX","chrX","chrXI","chrXII","chrXIII","chrXIV","chrXV","chrXVI","chrM","rDNA-10R"),seqlengths=c(230218,813184,316620,1531933,576874,270161,1090940,562643,439888,745751,666816,1078177,924431,784333,1091291,948066,85779,113097), isCircular=c(rep(F,16),T,F),genome="S288CrDNA")

pathdata <- "~/work/Ori/newRaw_nt/"
expmeg <- "JP3A_Megalodon_00_smdata.rds"
ex.name <- expmeg %>% str_remove("_Megalodon_00_smdata.rds")
ex.pulse <- 2
ex.b2a <- 0.02
EXP <- readRDS(pathdata %+% expmeg)  %>% filter(chrom!="chrM")
EXP_PLS <- PLSmaster(EXP,pulse0=ex.pulse,PLS.save=T,EXPname=ex.name %+% "_nt",b2a=ex.b2a)
# toprint <- as_tibble(EXP_PLS[[4]])
# toprint2 <- toprint %>% select(-1) %>% mutate_all(.,as.numeric)
# toprint3 <- bind_cols(toprint[,1],toprint2%>% round(.,2) %>% mutate_all(formatC, digit=4))
# formattable(toprint3,align = c("l", rep("c", NCOL(toprint3) - 1)))
```

As parsed data file were big, fork detection is performed on the split
data and then these results are merged using the PLS\_merging function
from the PLS\_function.r file.  
The merged data file keep the same organization with a slightly
simplified summary table. Those merged file are then used to produced
the figures and data discussed in the manuscript.

``` r
PLS_merging ("./","./","JP3A",suff="merged",file_list0="JP3A_nt_PLS_data.rds")
res <- readRDS("JP3Amerged_PLS_data.rds")
toprint4 <- bind_cols(res[[4]][,1],res[[4]] %>% select(-1) %>% round(.,2) %>% mutate_all(format, digit=4))
formattable(toprint4,align = c("l", rep("c", NCOL(toprint3) - 1)))
```
