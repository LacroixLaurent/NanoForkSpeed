### Megalodon data import
### LL 20210124
mega_parsing_v2 <- function(bam.in,out.file,ncores=1L,bs0=100,save.full=F)
{
require(parallel)
require(dplyr)
require(stringr)
require(furrr)
# import from bam: read_id, flag, chrom,start,end,seq,pos,prob
data0 <- try(system(paste0("samtools view ",bam.in," | cut -f 1,2,3,4,9,10,13,14"),intern=T))
# split lines
data1 <- mclapply(data0,function(x) strsplit(x,"\t")[[1]],mc.cores=ncores)
# reshape data lines
data2 <- mclapply(data1, function(x) 
{
id <- x[1]
if (x[2]=="0") {strand="+"} else {strand="-"}
chrom <- x[3]
start <- as.numeric(x[4])
end <- start+as.numeric(x[5])
if (strand=="+") {seq <- x[6]} else {seq <- as.character(Biostrings::reverseComplement(Biostrings::DNAString(x[6])))}
pos <- x[7] %>% str_remove(";") %>% str_remove(fixed("Mm:Z:T+B,"))
pos1 <- as.numeric(c(str_split(pos,",")[[1]]))
prob <- x[8] %>% str_remove("Ml:B:C,")
prob1 <- as.numeric(c(str_split(prob,",")[[1]]))/255
res <- list(id,chrom,start,end,strand,seq,pos1,prob1)
names(res) <- c("read_id","chrom","start","end","strand","seq","pos","prob")
return(res)
},mc.cores=ncores)

data3 <- mclapply(data2, function(x)
{
Tpos <- str_locate_all(x$seq,"T")[[1]][,1]
pos_out <- Tpos[cumsum(x$pos+1)]
prob_out <- x$prob
x$signal0 <- tibble(pos_out,prob_out)
x$pos <- NULL
x$prob <- NULL
x$seq <- NULL
return(x)
},mc.cores=ncores)

data4 <- mclapply(data3, function(x)
{
y <- tibble(read_id=x$read_id,chrom=x$chrom,start=x$start,end=x$end,strand=x$strand,signal_raw=list(x$signal0))
},mc.cores=ncores)

data5 <- do.call(bind_rows,data4)

plan(multicore, workers = ncores)

data6 <- data5 %>% 
	mutate(newsignal = future_pmap(
		., function(start,end,strand,signal_raw,...) {
			if(strand=="+") {
				positions= start - 1 + signal_raw %>% pull(pos_out);
				Bprob=signal_raw %>% pull(prob_out)
			}
			if(strand=="-") {
				positions= end + 1 - signal_raw %>% pull(pos_out);
				Bprob=signal_raw %>% pull(prob_out)
			}
			out=tibble(positions,Bprob) %>% arrange(positions)
		}))

data7 <- data6 %>%
	mutate(signalb=future_map(newsignal,
		function(y) {
		bs=bs0
		y %>%
		mutate(positions = round(positions/bs)*bs) %>%
		group_by(positions) %>%
		summarise(signal = mean(Bprob,na.rm=T), .groups = "drop")
		}))
if (save.full)
{saveRDS(data7, file=paste0(out.file,"_fulldata.rds"))}

data8 <- data7 %>% select(-c("signal_raw","newsignal")) %>% arrange(chrom,start)
saveRDS(data8, file=paste0(out.file,"_bin",bs0,"_mindata.rds"))
}
