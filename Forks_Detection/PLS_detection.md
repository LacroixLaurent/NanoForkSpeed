PLS data generation
================

This is a typical example how we proceed to detect forks, initiations
and terminations from the reads base\_called with megalodon and our
model (script=*basecalling\_sample.sh*) then parsed with the
*Parsing\_function4Megalodon.r* (script=*basecalling\_sample.r*).

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
ex.name <- "Exp_Test"
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
are organized as a list or 4 elements:  
1- a list of tibble containing:  
1.1- the reads filtered for their size and the presence of 3 and more
linear segments.  
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
source("./PLS_function.r")

seqinf <- Seqinfo(seqnames=c("chrI","chrII","chrIII","chrIV","chrV","chrVI","chrVII","chrVIII","chrIX","chrX","chrXI","chrXII","chrXIII","chrXIV","chrXV","chrXVI","chrM","rDNA-10R"),seqlengths=c(230218,813184,316620,1531933,576874,270161,1090940,562643,439888,745751,666816,1078177,924431,784333,1091291,948066,85779,113097), isCircular=c(rep(F,16),T,F),genome="S288CrDNA")

pathdata <- "~/work/Ori/newRaw_nt/"
expmeg <- "JP3A_Megalodon_00_smdata.rds"
ex.name <- "Exp_Test"
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
from the *PLS\_function.r* file.  
The merged data file keep the same organization with a slightly
simplified summary table. Those merged file are then used to produce the
figures and data discussed in the manuscript [Theulot et al.,
2022](https://doi.org/XX.XXXXX/JOURNAL/REF).

``` r
PLS_merging ("./","./","Exp_Test",suff="_merged",file_list0="Exp_Test_nt_PLS_data.rds")
res <- readRDS("Exp_Test_merged_PLS_data.rds")
toprint4 <- bind_cols(res[[4]][,1],res[[4]] %>% select(-1) %>% round(.,2) %>% mutate_all(format, digit=4))
formattable(toprint4,align = c("l", rep("c", NCOL(toprint4) - 1)))
```

### Outpout data format explanation

<<<<<<< HEAD
##### PLSmaster output

1 PLS\_data  
1.1: PLS\_data$allRDP3

This tibble contains all the reads analyzed containing at least 3 linear
segments after “Piecewise Linear Simplification”.  
- read\_id= read identifier  
- chrom= mapped chromosome  
- start= start of the mapped read (1-based)  
- end= end of the mapped read  
- strand= strand of the mapped read  
- gap\_pos= position of gaps introduced during the alignment  
- signalr= BrdU signal along the mapped read (positions=chromosomal
coordinate, Bprob= raw BrdU probability from megodon witu our model),
signal=smoothed BrdU signal using first a 100nt rolling mean then a
2500nt rolling weighted mean with a gaussian weight function)  
- length= length of the read  
- smBmedy= median of the smoothed BrdU signal  
- Bmedy= median of the raw BrdU signal  
- RDP= Piecewise Linear Simplification of the smoothed signal using the
Ramer Douglas Peucker algorithm (x,y = positions of extremities of the
linear segments)  
- RDP.length= number segment+1  
- sl= slope results (myslope1 function)  
- sl2= slope results after letter affectation (myslope2 function)  
- pat= read slope pattern after concatenation of sl2$sl.pat2  
- patR= position of rightward forks in the read slope pattern  
- patL= position of leftward forks in the read slope pattern  
- fork.R= position of rightward forks indicating if the forks can be
also a termination (not\_ter), positions of the identifed start of the
pulse (X0,Y0), end of the pulse (X1,Y1) and end of the lastunamibiguous
chase segment (X2, which might not coincide with the end of the chase),
average speed during the pulse (speed), average signal slope during the
pulse (sl\_pulse) and the begining of the chase (sl\_chase) and BrdU
signal amplitude (d.Y, &gt;0 for for rightward forks)  
- fork.L= position of leftward forks indicating if the forks can be also
a termination (not\_ter), positions of the identifed start of the pulse
(X0,Y0), end of the pulse (X1,Y1) and end of the lastunamibiguous chase
segment (X2, which might not coincide with the end of the chase),
average speed during the pulse (speed), average signal slope during the
pulse (sl\_pulse) and the begining of the chase (sl\_chase) and BrdU
signal amplitude (d.Y, &lt;0 for for leftward forks)  
- forks= concatenation of fork.R and forkL  
- n.forks= number of forks detected in the read

1.2: PLS\_data$with\_forks

same structure as PLS\_data$allRDP3
=======
#### PLSmaster output

1: PLS\_data 1.1: PLS\_data$allRDP3 This tibble contains all the reads
analyzed containing at least 3 linear segments after “Piecewise Linear
Simplification”. - read\_id= read identifier - chrom= mapped chromosome
- start= start of the mapped read (1-based) - end= end of the mapped
read - strand= strand of the mapped read - gap\_pos= position of gaps
introduced during the alignment - signalr= BrdU signal along the mapped
read (positions=chromosomal coordinate, Bprob= raw BrdU probability from
megodon witu our model), signal=smoothed BrdU signal using first a 100nt
rolling mean then a 2500nt rolling weighted mean with a gaussian weight
function) - length= length of the read - smBmedy= - Bmedy= - RDP= -
RDP\_length= - sl= - sl2= - pat= - patR - patL - fork.R - fork.L -
forks= - n.forks=
>>>>>>> b75e6609e7d70a3ade7a8c9033776f127c0ace24

#### PLS\_merging output