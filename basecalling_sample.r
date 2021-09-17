## Rscript to parse megalodon data
# LL20210917
Exp <- "EXP"
setwd(paste0("/mnt/data2/Laurent/Mega/",Exp,"__mega"))
system("samtools view -H mod_mappings.bam > header")
system("samtools view mod_mappings.bam | split - mod_splitted_ -l 50000 -a 2 -d --filter='cat header - | samtools view -b - > $FILE.bam'")
# this use of cat works with the ubuntu version of cat but not with the Apple/FreeBSD version.
nbam <- length(dir()[grep("mod_splitted_",dir())])-1

`%+%` <- paste0
source("/mnt/data/Laurent/mega_parsing_function_v2.r")

if (nbam>9) {
	for (i in 0:9) {
		bamf <- paste0("mod_splitted_0",i,".bam")
		outf <- paste0(Exp,"_Megalodon_0",i)
		mega_parsing_v2(bam.in=bamf,out.file=outf,ncores=4,bs0=100)
		gc()
	}
	for (i in 10:nbam) {
		bamf <- paste0("mod_splitted_",i,".bam")
		outf <- paste0(Exp,"_Megalodon_",i)
		mega_parsing_v2(bam.in=bamf,out.file=outf,ncores=4,bs0=100)
		gc()
	}
}else{
	for (i in 0:nbam) {
		bamf <- paste0("mod_splitted_0",i,".bam")
		outf <- paste0(Exp,"_Megalodon_0",i)
		mega_parsing_v2(bam.in=bamf,out.file=outf,ncores=4,bs0=100)
		gc()
	}
}
toto <- dir()
file_list <- toto[grep("mindata.rds",toto)]
out <- do.call(bind_rows,lapply(file_list, readRDS)) %>% arrange(chrom,start)
saveRDS(out,file=paste0(Exp,"_Megalodon_bin100.rds"))
### cleaning
system(paste0("rm ",Exp,"_Megalodon_??_bin100_mindata.rds"))
system("rm mod_splitted_*.bam")
system("rm header")

# this part of the script used our parsing function
# the ncores parameter should be used to set the number of core available
# keeping in mind tha tR use a lot of memory.
# on our computer, more than 40Go of RAM were required.
# the parsing script can also output the non binned data,
# but the size of the output file will be in the order of Go instead of tens of Mo.
# the binned file generated is the one we used to look to forks
