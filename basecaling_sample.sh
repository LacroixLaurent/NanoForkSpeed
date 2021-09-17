# Sample of the basecalling procedure used with our BrdU trainde model
# This requires the ONT megalodon and guppy software
# We used guppy v4.4.1 GPU and megalogon v2.2.9
# a symbolic link for guppy should be present in the working directory
# ont-guppy -> /users/rce/hyrien/src/ont-guppy_4.4.1
ln -s ~/src/ont-guppy_4.4.1 ./ont-guppy
# the configuration file should point to the modified base model
# see BrdU_Configutation4megalodon.cfg
# (the model file is provided (model_adversial_bug.json)) (but too big for github?)
# the cfg file should be placed in the src/ont-guppy_4.4.1/data folder (not sure it is mandatory)
# a conda env was create using the yml file provided within miniconda3
conda env create -f megalodon.yml
# then within this env, we used megalodon with the following parameters
conda activate mega
REF=/mnt/data/Refgen/S288CwExtrarDNA_ROMAN.fa
# we used a modified version of the S288C reference genome (R64.2.1)
# to which we add an extra chromosome (rDNA-10R) containing 10 rDNA tandem repats
CFG=BrdU_Configuration4megalodon.cfg
INPUT=/mnt/data3/RawData/EXP_fast5
OUTPUT=/mnt/data2/Laurent/Mega/EXP_mega

nohup megalodon $INPUT --outputs mod_mappings --reference $REF  --output-dir $OUTPUT  --guppy-config $CFG --disable-mod-calibration --overwrite --device cuda:0 --processes 20 --guppy-timeout 90 --mod-min-prob 0  &
# please notice than to experiment with many long reads, we had to increase the guppy-timeout to 600
# also notice that this computer was equiped with a NVidia RTX2080Ti GPU

# this generate a megalon's style mod-mappings.bam file in the EXP_mega folder
# this bam file can be big and we had to split it to be able to process it with our parsing script
# the following step are done within R (v4.0.5) on the same computer used for the basecalling

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
