### Megalodon data import
### LL 20210919
mega_parsing <- function(bam.in,out.file,ncores=1L)
{
	require(parallel)
	require(dplyr)
	require(stringr)
	require(furrr)
	require(tidyr)
	options(future.rng.onMisuse = "ignore")
	suppressMessages(require(RcppRoll))
	min_len <- 5000
	w_pg <- dnorm(1:2500,1251,300)
	min_gap <- 100

	# import from bam: read_id, flag, chrom,start,end,seq,pos,prob
	data0 <- try(system(paste0("samtools view ",bam.in," | cut -f 1,2,3,4,9,10,13,14"),intern=T))
	# split lines
	data1 <- mclapply(data0,function(x) strsplit(x,"\t")[[1]],mc.cores=ncores)
	rm(data0)
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
	rm(data1)
	# translate T position to sequence position and affect B probability
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
	rm(data2)
	data4 <- mclapply(data3, function(x)
	{
		y <- tibble(read_id=x$read_id,chrom=x$chrom,start=x$start,end=x$end,strand=x$strand,signal_raw=list(x$signal0))
	},mc.cores=ncores)
	rm(data3)
	data5 <- do.call(bind_rows,data4) %>% filter((end-start)>min_len)
	rm(data4)

	plan(multicore, workers = ncores)

	data6 <- data5 %>%
		mutate(newsignal = future_pmap(
			., function(start,end,strand,signal_raw,...) {
				# translate read coordinates to genomic coordinates
				if(strand=="+") {
					positions= start - 1 + signal_raw %>% pull(pos_out);
					Bprob=signal_raw %>% pull(prob_out)
				}
				# translate read coordinates to genomic coordinates and reverse the signal if strand==-
				if(strand=="-") {
					positions= end + 1 - signal_raw %>% pull(pos_out);
					Bprob=signal_raw %>% pull(prob_out)
				}
				out=tibble(positions,Bprob) %>% arrange(positions)
			})) %>%
		# remove the raw signal
		select(-signal_raw)
	rm(data5)
	# detecting gaps
	data7 <- data6 %>%
		mutate(gap_pos=future_map(newsignal, function(x) {
			y<-x %>% drop_na()
			# measure of the discontinuities in binned positions
			y$dif <- c(diff(y$positions),1)
			# mapping of the gap that are longer than the min_gap threshold
			gap_start <- y[which(y$dif>min_gap),]$positions
			gap_end <- y[which(y$dif>min_gap)+1,]$positions
			if (length(gap_start)>0) {
				res=tibble(gap_start,gap_end)
			}else{
				res=tibble(gap_start=0,gap_end=0)
			}
			res$gap_width <- res$gap_end-res$gap_start+1
			return(res)
		}))
	rm(data6)
	# smoothing
	data8 <- data7 %>%
		# add NA to fill all the coordinates
		mutate(signal1=future_pmap(.,function(newsignal,start,end,...) {tibble(positions=start:end,y=NA) %>% left_join(newsignal, by="positions") %>% select(-y)})) %>%
		mutate(signalr=future_map(signal1, function(x) {
			# do a first smooth to mimick the binning
			sm0 <- roll_mean(x$Bprob,by=1,align="center",n=100,na.rm=T,fill=NA)
			# smooth the data with a gaussian function
			x$signal <- roll_mean(sm0,by=1,align="center",weights=w_pg,normalize=T,na.rm=T,fill=NA)
			# remove the NA
			y <- x %>% drop_na();
			return(y)})) %>%
		select(-c(newsignal,signal1)) %>% arrange(chrom,start)

	# save the data
	saveRDS(data8, file=paste0(out.file,"_smdata.rds"))
}
