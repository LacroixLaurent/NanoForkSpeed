## Rscript to parse megalodon data

nc=4L

Exp <- "BT1_run4"
setwd(paste0("./data/",Exp,"_mega"))
system("samtools view -H mod_mappings.bam > header")
system("samtools view mod_mappings.bam | split - mod_splitted_ -l 50000 -a 2 -d --filter='cat header - | samtools view -b - > $FILE.bam'")
# this use of cat works with the ubuntu version of cat but not with the Apple/FreeBSD version.
nbam <- length(dir()[grep("mod_splitted_",dir())])-1

`%+%` <- paste0
source("./Parsing_function4megalodon.r")

if (nbam>9) {
	for (i in 0:9) {
		bamf <- paste0("mod_splitted_0",i,".bam")
		outf <- paste0(Exp,"_Megalodon_0",i)
		mega_parsing(bam.in=bamf,out.file=outf,ncores=nc)
		gc()
	}
	for (i in 10:nbam) {
		bamf <- paste0("mod_splitted_",i,".bam")
		outf <- paste0(Exp,"_Megalodon_",i)
		mega_parsing(bam.in=bamf,out.file=outf,ncores=nc)
		gc()
	}
}else{
	for (i in 0:nbam) {
		bamf <- paste0("mod_splitted_0",i,".bam")
		outf <- paste0(Exp,"_Megalodon_0",i)
		mega_parsing(bam.in=bamf,out.file=outf,ncores=nc)
		gc()
	}
}

### cleaning
system("rm mod_splitted_*.bam")
system("rm header")

# This part of the script used our parsing function
# the nc parameter should be used to set the number of core available
# keeping in mind that R uses a lot of memory.
#Our computer is a Dell Precision 7920 Tower equipped with:
# Processor: 2x Intel Xeon Silver 4214
# RAM: 192Gb
# GPU: GeForce RTX 2080 Ti, 11Gb DDR5 (CUDA v11.2)
# On our computer, more than 62Go of RAM were required.
# the parsing script output both smoothed and raw BrdU probability data.
# The size of the output is big.
# Forks detection is performed on the split files

library("devtools")
library(magrittr)
session_info() %>% capture.output(file="base_calling_session_info.txt")
