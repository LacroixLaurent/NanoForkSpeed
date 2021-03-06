NFS data generation
================

This is a typical example of how we proceed to detect forks, initiations
and terminations from the reads base called with megalodon and our model
(script=*basecalling_sample.sh*) then parsed with the
*Parsing_function4Megalodon.r* (script=*basecalling_sample.r*).

### Background threshold determination

The first step of the analysis requires to check the distribution of the
signal in order to set the b2a threshold, used to identify the
background of the signal.  
This script outputs a PDF file corresponding to the distribution of the
BrdU signal in 1kb bins for the reads filtered the presence of 3 and
more linear segments.

``` r
suppressMessages(library(tidyverse))
`%+%`<- paste0
source("./helper_function.r")
source("./NFS_function.r")
pathdata <- "~/work/Ori/NFS_paper/Zenodo_upload/"
ex.name <- "BT1_run4"
expmeg <- ex.name %+% "_Megalodon_00_smdata.rds"

EXP <- readRDS(pathdata %+% expmeg)  %>% filter(chrom!="chrM")
plot_signal(EXP,EXPname=ex.name,EXP_b2a.thr0=0.02,alldata=F,nreads=5000,saved=T,plotit=F)
```

It is also possible to output the distribution for all the reads (of at
least 5kb) but usually this distribution is less informative due to the
overwhelming presence of reads with very low signal. From the second
plot of the PDF file (a zoom of the first plot of the PDF), it is
possible to draw a line between the 0 level signal (background) and the
informative signal which usually come in two wide distribution. One of
medium level mode corresponding to the DNA replicated during the chase
and one with higher level mode corresponding to the DNA replicated
during the pulse.  
The position of this line was the threshold we fixed for our further
analysis. This threshold was quite similar between all the experiments
we analysed, thus a fixed threshold of 0.02 was chosen to simplify the
pipeline.  
We also limit to 5000 the number of reads used for this evaluation to
speed up the computation but all the reads can be used by setting
nreads=NA.

### NFS analysis

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
suppressMessages(library(kmlShape))
suppressMessages(library(tidyverse))
`%+%`<- paste0
source("./NFS_function.r")

pathdata <- "~/work/Ori/NFS_paper/Zenodo_upload/"
ex.name <- "BT1_run4"
expmeg <- ex.name %+% "_Megalodon_00_smdata.rds"
ex.save <- str_remove(expmeg,"_smdata.rds")
ex.pulse <- 2
ex.b2a <- 0.02
EXP <- readRDS(pathdata %+% expmeg)  %>% filter(chrom!="chrM")
EXP_NFS <- NFSmaster(EXP,pulse0=ex.pulse,NFS.save=T,EXPname=ex.save,b2a=ex.b2a)
```

As parsed data files are big, fork detection is performed on the split
data and then these results are merged using the NFS_merging function
from the *NFS_function.r* file.  
The merged data file keep the same organization with a slightly
simplified summary table. Those merged file are then used to produce the
figures and data discussed in the manuscript [Theulot et al.,
2022](https://doi.org/10.1038/s41467-022-31012-0).

``` r
source("./NFS_function.r")
NFS_merging ("./","./","BT1_run4",suff="_merged",file_list0="BT1_run4_Megalodon_00_NFS_data.rds")
```

### Session Info

``` r
library("devtools")
library(magrittr)
session_info() %>% capture.output(file="forks_detection_session_info.txt")
```

### Outpout data format explanation

##### NFSmaster output

1- NFS_data  
1.1- NFS_data$allRDP3

This tibble contains all the reads analysed containing at least 3 linear
segments after ???Piecewise Linear Simplification???.  
- read_id= read identifier  
- chrom= mapped chromosome  
- start= start of the mapped read (1-based)  
- end= end of the mapped read  
- strand= strand of the mapped read  
- gap_pos= position of gaps introduced during the alignment  
- signalr= BrdU signal along the mapped read (positions=chromosomal
coordinate, Bprob= raw BrdU probability from megodon with our model),
signal=smoothed BrdU signal using first a 100nt rolling mean then a
2500nt rolling weighted mean with a gaussian weight function)  
- length= length of the read  
- smBmedy= median of the smoothed BrdU signal  
- Bmedy= median of the raw BrdU signal  
- RDP= Piecewise Linear Simplification of the smoothed signal using the
Ramer Douglas Peucker algorithm (x,y = positions of extremities of the
linear segments)  
- RDP.length= number of RDP anchoring points=number of segments+1  
- sl2= slope results after letter affectation  
- forks= position of the forks indicating the positions of the
identified start of the pulse (X0,Y0), end of the pulse (X1,Y1) and end
of the last un ambiguous chase segment (X2, which might not coincide
with the end of the chase), average speed during the pulse (speed),
average signal slope during the pulse (sl_pulse) and the beginning of
the chase (sl_chase) and BrdU signal amplitude (d.Y, \>0 for rightward
forks,\<0 for leftward forks )  
- n.forks= number of forks detected in the read

1.2- NFS_data$with_forks

same structure as NFS_data$allRDP3 after filtering out of the forks
overlapping an alignment gap.

2- forks  
- chrom= mapped chromosome  
- strand= strand of the mapped read  
- direc= direction of the forks (L=left, R=right)  
- speed= estimated average speed for the pulse duration  
- d.Y= amplitude of the forks (\>0 for R, \<0 for L)  
- type= leading or lagging  
- X0= start of BrdU incorporation (position of the B/P and N/B
transitions for rightward and leftward forks)  
- X1= start of the thymidine chase (position of the (P\|A)/N and
P/(A\|N) transitions for rightward and leftward forks, respectively)  
- X2= end of the last non ambiguous chase segment  
- read_id= read identifier  
- trac= raw BrdU signal from 1kb before X0 to 50kb after X0 by 100nt
bin  
- exp= name of the experiment or of subfile

3- initer  
- chrom= mapped chromosome  
- strand= strand of the mapped read  
- read_id= read identifier  
- x0= X0 of the left (respectively X1 and right) fork for initiation
(respectively termination)  
- x1= X0 of the right (respectively X1 and left) fork for initiation
(respectively termination)  
- center= center of the x0-x1 segment (center=(x0+x1)/2)  
- spL= speed of the left fork  
- spR= speed of the right fork  
- yR= amplitude of the left fork  
- yL= amplitude of the right fork  
- type = Initiation or Termination

4- stats  
- EXPname= name of the experiment or of subfile  
- n_reads= number of starting reads  
- n_reads(len\>minlen)= number of reads whose length is above the minlen
set in the NFSmaster function (default=5kb)  
- sumlength(\>minlen)= sum of length of the reads longer than minlen
(default=5kb)  
- read_med_len_with_forks= median of length for reads with forks  
- b2a.thr= threshold used to identify non-BrdU containing part of the
reads  
- B_median= median of the reads BrdU signal median  
- B_mad= median absolute deviation of the reads BrdU signal median  
- n_reads(RDP\>3)= number of reads with 3 or more segment after
piece-wise linear simplification with RDP  
- n_reads_w\_forks= number of reads with forks  
- n_forks= number of forks  
- n_forks(no_gap)= number of forks after filtering out forks overlapping
with alignment gaps (??1kb)  
- speed.med= median of estimated speeds  
- dY.median= median of forks amplitude  
- nb_init= number of initiations  
- nb_ter= number of terminations  
- forkdens= number of forks per Mb of evaluated reads
(n_forks/sumlength*1e6)  
- initdens= number of initiations per Mb of evaluated reads
(n_init/sumlength*1e6)  
- terdens= number of terminations per Mb of evaluated reads
(n_ter/sumlength\*1e6)

#### NFS_merged output

1- reads  
This tibble contains all the reads analysed containing at least one
fork.  
- read_id= read identifier  
- chrom= mapped chromosome  
- start= start of the mapped read (1-based)  
- end= end of the mapped read  
- strand= strand of the mapped read  
- gap_pos= position of gaps introduced during the alignment  
- signalr= BrdU signal along the mapped read (positions=chromosomal
coordinate, Bprob= raw BrdU probability from megodon with our model),
signal=smoothed BrdU signal using first a 100nt rolling mean then a
2500nt rolling weighted mean with a gaussian weight function)  
- length= length of the read  
- smBmedy= median of the smoothed BrdU signal  
- Bmedy= median of the raw BrdU signal  
- RDP= Piecewise Linear Simplification of the smoothed signal using the
Ramer Douglas Peucker algorithm (x,y = positions of extremities of the
linear segments)  
- sl2= slope results after letter affectation  
- forks= position of the forks indicating the positions of the
identified start of the pulse (X0,Y0), end of the pulse (X1,Y1) and end
of the last un ambiguous chase segment (X2, which might not coincide
with the end of the chase), average speed during the pulse (speed),
average signal slope during the pulse (sl_pulse) and the beginning of
the chase (sl_chase) and BrdU signal amplitude (d.Y, \>0 for rightward
forks,\<0 for leftward forks )  
- n.forks= number of forks detected in the read

2- forks  
- chrom= mapped chromosome  
- strand= strand of the mapped read  
- direc= direction of the forks (L=left, R=right)  
- speed= estimated average speed for the pulse duration  
- d.Y= amplitude of the forks (\>0 for R, \<0 for L)  
- type= leading or lagging  
- X0= start of BrdU incorporation (position of the B/P and N/B
transitions for rightward and leftward forks)  
- X1= start of the thymidine chase (position of the (P\|A)/N and
P/(A\|N) transitions for rightward and leftward forks, respectively)  
- X2= end of the last non ambiguous chase segment  
- read_id= read identifier  
- trac= raw BrdU signal from 1kb before X0 to 50kb after X0 by 100nt
bin  
- exp= name of the experiment or of subfile

3- initer  
- chrom= mapped chromosome  
- strand= strand of the mapped read  
- read_id= read identifier  
- x0= X0 of the left (respectively X1 and right) fork for initiation
(respectively termination)  
- x1= X0 of the right (respectively X1 and left) fork for initiation
(respectively termination)  
- center= center of the x0-x1 segment (center=(x0+x1)/2)  
- spL= speed of the left fork  
- spR= speed of the right fork  
- yR= amplitude of the left fork  
- yL= amplitude of the right fork  
- type= Initiation or Termination

4- stats  
- Expname= name of the experiment or of subfile  
- n_reads= number of starting reads  
- n_reads2= number of reads whose length is above the minlen set in the
NFSmaster function (default=5kb)  
- sumlength= sum of length of the reads longer than minlen
(default=5kb)  
- b2a.thr= threshold used to identify non-BrdU containing part of the
reads  
- n_reads_RDP3= number of reads with 3 or more segment after piece-wise
linear simplification with RDP  
- n_reads_w\_forks= number of reads with forks  
- n_forks= number of forks after filtering out forks overlapping with
alignment gaps (??1kb)  
- speed_med= median of estimated speeds  
- dY_med= median of forks amplitude  
- nb_init= number of initiations  
- nb_ter= number of terminations  
- fork_dens= number of forks per Mb of evaluated reads
(n_forks/sumlength*1e6)  
- init_dens= number of initiations per Mb of evaluated reads
(n_init/sumlength*1e6)  
- ter_dens= number of terminations per Mb of evaluated reads
(n_ter/sumlength\*1e6)
