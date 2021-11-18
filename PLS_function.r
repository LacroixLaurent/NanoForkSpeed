### Function used in Theulot et al 2021
### Forks detection, orientation and speed measurement
### as well as Initiation and Termination detection
### work with splitted files from the parsing procedure done by 50000 reads chunks
### smoothing is performed at the parsing level
#### LL 20211108


### the main function
PLS.detect <- function(tibble.test,RDP.eps=0.1,slope.thr=0.25,b2a.thr=0.02,pulse=2)
	# signal already smoothed during the parsing
{
	suppressMessages(require(kmlShape))
	suppressMessages(require(tidyverse))
	## internal helper functions
	# General use parameters
	RDP.spar=NA
	gap_ext <- 1000L
	min_gap <- 100L

	# Ramer-Douglas-Peucker simplification using a homemade smoothing function)
	myRDP <- function(x,...)
	{
		DouglasPeuckerEpsilon(x$positions,x$signal,RDP.eps,spar=RDP.spar)
	}

	# Evaluation of the slopes of the simplified segments
	myslope1 <- function(RDP,...)
	{
		sl.x <- (RDP$x[1:(length(RDP$x)-1)]+RDP$x[-1])/2
		sl.y <- diff(RDP$y,1)/diff(RDP$x,1)*1e5
		d.y <- diff(RDP$y,1)
		sl.y2 <- sl.y
		sl.y2[abs(sl.y)<slope.thr] <- 0
		y0 <- lag(RDP$y)[-1]
		y1 <- RDP$y[-1]
		res <- tibble(sl.x,sl.y,sl.y2,d.y,y0,y1)
		return(res)
	}

	# Slope classification into B (before the pulse),A (after the pulse),P,N
	myslope2 <- function(sl,...)
	{
		sl %>%
			mutate(sl.pat=map(sl.y2,sign2let)) %>%
			unnest(cols = c(sl.pat)) %>%
			mutate(sl.pat2=case_when((sl.pat=="0" & meanB<b2a.thr) ~ "B",(sl.pat=="0" & meanB>=b2a.thr) ~ "A",(meanB<b2a.thr) ~ "B",T ~ sl.pat))
	}

	# Letter attribution for P (positive) and N (negative) slopes
	sign2let <- function(x)
	{
		if (x==0) {return("0")}
		if (x>0) {return("P")}
		if (x<0) {return("N")}
	}

	# Forks pattern used for forks identification
	patRnew <- "BP(P|A)*N+"
	patRnew2 <- "BP(P|A)*"
	patLnew <- "P+(N|A)*NB"
	patLnew2 <- "(N|A)*NB"

	# Rightward progressing fork mapping function
	map.fR <- function(sl2test,patRtest,t.pulse=2)
	{
		if (nrow(patRtest)!=0)
		{
			patRtest %>%
				mutate(X0=sl2test$x1[st]) %>%
				mutate(X1=sl2test$x1[mid]) %>%
				mutate(X2=sl2test$x1[en]) %>%
				mutate(sl_pulse=map2_dbl(st,mid,function(x,y) mean(sl2test$sl.y2[(x+1):y]))) %>%
				mutate(sl_chase=map2_dbl(mid,en,function(x,y) mean(sl2test$sl.y2[(x+1):y]))) %>%
				mutate(d.Y=map2_dbl(st,mid,function(x,y) sum(sl2test$d.y[(x+1):y]))) %>%
				mutate(speed=map2_dbl(X0,X1, function(x,y) round(abs(y-x)/t.pulse))) %>%
				mutate(Y0=sl2test$y1[st]) %>%
				mutate(Y1=sl2test$y1[mid]) %>%
				select(-c(st,en,mid)) -> out
		}else{
			out <-tibble(X0=as.integer(NA),X1=as.integer(NA),X2=as.integer(NA),sl_pulse=as.double(NA),sl_chase=as.double(NA),d.Y=as.double(NA),speed=as.double(NA),Y0=as.integer(NA),Y1=as.integer(NA),not_ter=as.logical(NA))
		}
		return(out)
	}

	# Leftward progressing fork mapping function
	map.fL <- function(sl2test,patLtest,t.pulse=2)
	{
		if (nrow(patLtest)!=0)
		{
			patLtest %>%
				mutate(X0=sl2test$x0[st]) %>%
				mutate(X1=sl2test$x0[mid]) %>%
				mutate(X2=sl2test$x0[en]) %>%
				mutate(sl_pulse=map2_dbl(st,mid,function(x,y) mean(sl2test$sl.y2[(x-1):y]))) %>%
				mutate(sl_chase=map2_dbl(mid,en,function(x,y) mean(sl2test$sl.y2[(x-1):y]))) %>%
				mutate(d.Y=map2_dbl(st,mid,function(x,y) sum(sl2test$d.y[(x-1):y]))) %>%
				mutate(speed=map2_dbl(X0,X1, function(x,y) round(abs(y-x)/t.pulse))) %>%
				mutate(Y0=sl2test$y0[st]) %>%
				mutate(Y1=sl2test$y0[mid]) %>%

				select(-c(st,en,mid)) -> out
		}else{
			out <-tibble(X0=as.integer(NA),X1=as.integer(NA),X2=as.integer(NA),sl_pulse=as.double(NA),sl_chase=as.double(NA),d.Y=as.double(NA),speed=as.double(NA),Y0=as.integer(NA),Y1=as.integer(NA),not_ter=as.logical(NA))
		}
		return(out)
	}


	## Main  detection function
	test2 <- tibble.test %>%
		# RDP simplification
		mutate(RDP=map(signalr,myRDP)) %>%
		# selection of reads with 3 or more segment to look for forks
		mutate(RDP.length=map_int(RDP,function(x) nrow(x))) %>%
		filter(RDP.length>3) %>%
		# compute mean B signal for the RDP segments
		mutate(meanB= map2(signalr,RDP,function(tibble1,tibble2) {
			tibble2 %>%
				mutate(x_end = lead(x)) %>%   # shifted copy of the x
				filter(!is.na(x_end)) %>%     # remove last segment with only a start
				mutate(mean_Bscore = map2_dbl(x,x_end,function(x, x_end) {
					tibble1 %>%
						filter(positions>=x & positions<=x_end) %>%
						pull(signal) %>%
						mean(na.rm=T)
				}))})) %>%
		# computing slopes for the RDP segments
		mutate(slopes=map(RDP,myslope1)) %>%
		# add mean B signal to the slopes
		mutate(sl=map2(meanB,slopes, function(z,y) {add_column(y, meanB=z %>% pull(mean_Bscore), x0=z %>% pull(x),x1=z %>% pull(x_end))})) %>%
		# slope to letter transformation
		mutate(sl2=map(sl,myslope2)) %>% select(-c(slopes,meanB))

	# letters concatenation
	test2$pat <- sapply(test2$sl2, function(x) paste0(x$sl.pat2,collapse=""))
	# R_forks pattern search
	test2$patR <- lapply(test2$pat, function(x) (
		{
			A <- str_locate_all(x, patRnew)[[1]]
			B <-do.call(rbind,str_locate_all(str_extract_all(x, patRnew)[[1]],patRnew2))
			midR <-A[,1]+B[,2]-1
			not_ter <- sapply(A[,2], function(pat) str_sub(x,start=pat+1,end=pat+1)!="B")
			res <- as_tibble(cbind.data.frame(A,midR,not_ter))
			colnames(res) <- c("st","en","mid","not_ter")
			return(res)
		})
	)
	# L_forks pattern search
	test2$patL <- lapply(test2$pat, function(x) (
		{
			A <- str_locate_all(x, patLnew)[[1]]
			B <-do.call(rbind,str_locate_all(str_extract_all(x, patLnew)[[1]],patLnew2))
			midL <-A[,1]+B[,1]-1
			not_ter <- sapply(A[,1], function(pat) str_sub(x,start=pat-1,end=pat-1)!="B")
			res <- as_tibble(cbind.data.frame(A,midL,not_ter))
			colnames(res) <- c("en","st","mid","not_ter")
			return(res)
		})
	)

	test3 <- test2 %>%
		# add pulse duration in the forks (L&R) tibbles
		mutate(fork.R=map2(sl2,patR, function(x,y) map.fR(x,y,t.pulse=pulse))) %>%
		mutate(fork.L=map2(sl2,patL, function(x,y) map.fL(x,y,t.pulse=pulse))) %>%
		# merging R and L forks tibbles
		mutate(forks=map2(fork.R,fork.L, function(w,z) (bind_rows(w,z) %>%
																											# remove empty lines
																											filter(!is.na(X0)) %>%
																											# remove termination type forks
																											filter(not_ter) %>%
																											arrange(X0)))) %>%
		# compute number of forks per read
		mutate(n.forks=map_int(forks,nrow))
	test4 <- test3 %>%
		# remove reads without forks (after ter filtering)
		filter(n.forks>0) %>%
		# extend gap + and - gap_ext=1000nt to catch artefactual forks
		mutate(gap_pos=map(gap_pos, function(x) {
			res <- x %>%
				mutate(start2=case_when((gap_width>min_gap)~(gap_start-gap_ext),T~0)) %>%
				mutate(end2=case_when((gap_width>min_gap)~(gap_end+gap_ext),T~0))
			return(res)
		})) %>%
		# tag forks overlapping with gaps
		mutate(forks=pmap(., function(forks,gap_pos,...) {
			# look for every fork if it overlaps the extended gap
			gapovl <- sapply(1:nrow(forks), function(i)
			{
				if (gap_pos$gap_width[1]>100) {
					res <- sapply(1:nrow(gap_pos), function(j)
					{(gap_pos$start2[j]>=min(forks$X0[i],forks$X2[i]) & gap_pos$start2[j]<=max(forks$X0[i],forks$X2[i]) | (gap_pos$end2[j]>=min(forks$X0[i],forks$X2[i]) & gap_pos$end2[j]<=max(forks$X0[i],forks$X2[i])))}) %>% sum %>% as.logical
				}else{
					res <- F
				}
			})
			# remove forks overlapping with extended gaps
			forks2 <- forks %>% add_column(gapovl) %>% filter(!gapovl)
			return(forks2)
		})) %>%
		# count of forks without gap per read
		mutate(n.forks2=map_int(forks, nrow)) %>%
		# remove reads with no forks after this filtering
		filter(n.forks2>0)
	res <- list(test3,test4)
	names(res) <- c("allRDP3","with_forks")
	return(res)
}

###


PLSmaster <- function(EXP,RDP.eps0=0.1, slope.thr0=0.25,pulse0=2,PLS.save=T,EXPname="EXP",minlen=5000,b2a=0.02)
	# in the output forks, L=left, R=right (in OKseq data, left=+ and right=-)
	# EXP is the tibble of the reads coming from the basecalling procedure with
	# the Parsing_function4megalodon.r function
	# RDP.eps0 set the tolerance for the RDP simplification
	# slope.thr0 is the threshold to call flat RDP segments
	# pulse0 is the duration of the BrdU pulse (in minute)
	# PLS.save is used to switch on the automatic saving of the results
	# EXPname is used to name the saved file
	# minlen can be adjusted select only reads longer than a threshold
	# b2a is the threshold used to call if a flat segment corresponds to B or A
{
	EXP_PLSall <- EXP %>%
		# remove the prefix "read_" if present
		mutate(read_id=map_chr(read_id, function(x) str_remove(x,"read_"))) %>%
		mutate(length=end-start) %>%
		# filtering according to a minimal length
		filter(length>minlen) %>%
		# compute median if the smoothed signal for eahc read
		mutate(smBmedy=map_dbl(signalr,function(z) median(z$signal,na.rm=T))) %>%
		# compute median if the raw signal for eahc read
		mutate(Bmedy=map_dbl(signalr,function(z) median(z$Bprob,na.rm=T)))

	# compute the sum of length of all the reads with length>minlen
	EXPlen <- EXP_PLSall %>% pull(length) %>% sum
	names(EXPlen) <- paste0("sumlength(>",(minlen/1000),"kb)")

	# filter for reads with 3 or more RDP segments
	# compute some stats
	EXP_stat <- c(EXPname,nrow(EXP),nrow(EXP_PLSall))
	names(EXP_stat)=c("EXPname","n_reads",paste0("n_reads(len>",minlen/1000,"kb)"))

	# threshold analysis
	### obsolete for megalodon 3+ version but still useful to detect Exp with high background
	### done with the Raw signal
	EXP_Bmed <- median(EXP_PLSall$Bmedy)
	EXP_Bmad <- mad(EXP_PLSall$Bmedy)
	if (b2a=="auto") {EXP_b2a.thr0 <- EXP_Bmed+3*EXP_Bmad} else {EXP_b2a.thr0 <- b2a}
	EXP_b2a <- c(EXP_b2a.thr0,EXP_Bmed,EXP_Bmad)
	names(EXP_b2a) <- c("b2a.thr","B_median","B_mad")

	# forks detection
	EXP_PLS_det <- PLS.detect(EXP_PLSall,RDP.eps=RDP.eps0,slope.thr=slope.thr0,b2a.thr=EXP_b2a.thr0,pulse=pulse0)

	# more stats computed
	EXP_PLSstat <- c(sapply(EXP_PLS_det,function(x) nrow(x)),sum(EXP_PLS_det[[1]]$n.forks),sum(EXP_PLS_det[[2]]$n.forks2))
	names(EXP_PLSstat) <- c("n_reads(RDP>3)","n_reads_w_forks","n_forks","n_forks(no_gap)")

	# outputting forks
	### for those forks, start and end correspond to the X0 and X2 limits to improve the RFD coverage
	# leading and lagging strand mapping
	### leading= strand+ forkR or strand- forkL
	### lagging= strand+ forkL or strand- forkR
	# maximum extension to extract the raw signal
	trac.xmax <- 50000
	if (nrow(EXP_PLS_det[[2]])>0)
	{
		EXPforks <- EXP_PLS_det[[2]] %>%
			select(chrom,strand,forks,signalr,read_id) %>%
			unnest(cols = c(forks)) %>%
			# selection of start and end for RFD
			mutate(st=pmin(X0,X2)) %>%
			mutate(en=pmax(X0,X2)) %>%
			# leading/lagging strand mapping
			mutate(type=case_when((strand=="+" & d.Y>0) | (strand=="-" & d.Y<0) ~ "leading", T ~ "lagging"))%>%
			# fork direction mapping
			mutate(direc=case_when(sign(d.Y)==1 ~ "R", T ~ "L")) %>%
			# signal extraction for mean_trace plotting and binning by 100
			### This is done with the smoothed signal. could be done with the raw one
			mutate(trac=pmap(list(signalr,X0,d.Y), function(y,x0,d) {
				if (d>0)
				{
					out <- y %>%
						filter(positions>=x0-1000 & positions<x0+trac.xmax) %>%
						mutate(positions = round(positions/100)*100) %>%
						group_by(positions) %>%
						summarise(signal = mean(Bprob,na.rm=T), .groups = "drop")
				}else{
					out <- y %>%
						filter(positions<=x0+1000 & positions>x0-trac.xmax) %>%
						mutate(positions = round(positions/100)*100) %>%
						group_by(positions) %>%
						summarise(signal = mean(Bprob,na.rm=T), .groups = "drop") %>% arrange(desc(positions))}
				return(out)
			})) %>%
			select(chrom,strand,st,en,direc,speed,d.Y,type,X0,X1,X2,read_id,trac) %>%
			# compute length of the extracted trace (in 100nt bin units)
			#			mutate(len=map_int(trac,nrow)) %>%
			# associated EXPname to fork
			mutate(exp=EXPname)
	}else{
		EXPforks <- tibble(chrom=character(),strand=character(),st=integer(),en=integer(),direc=character(),speed=integer(),d.Y=double(),type=character(),X0=integer(),X1=integer(),X2=integer(),read_id=character(),trac=list(),len=integer(),exp=character())
	}
	#	 computing the length (in 100nt bin) of the longest trace for padding
	#	EXP_maxlenfork <- max(EXPforks$len)
	#	names(EXP_maxlenfork) <- "fork_maxlen"

	# median speed
	EXP_med_speed <- EXPforks %>% pull(speed) %>% median(na.rm=T)
	names(EXP_med_speed) <- "speed_med"
	# median amplitude (of the RDP jump)
	EXP_med_dY <- EXPforks %>% pull(d.Y) %>% abs %>% median(na.rm=T)
	names(EXP_med_dY) <- "dY_median"

	# computing median length of reads above the minlen,
	# reads with 3 or more RDP segments
	# and reads with detected forks
	EXP_med_read_len <- c(
		median(EXP_PLSall$length),
		median(EXP_PLS_det[[1]]$length),
		median(EXP_PLS_det[[2]]$length))
	names(EXP_med_read_len) <- c(paste0("read_med_len_all(>",minlen/1000,"kb)"),"read_med_len_RDP3","read_med_len_with_forks")

	### Detection of Initiation and Termination
	# extrapolated center(Ext.center) ExC=X0L+spL/(spL+spR)*X0L_X0R
	# for ini, x0=X0L and x1=X0R, for ter, x0=X1R and x1=X1L
	# mapping R and L forks, and concatenating of the R and L letter for pattern search
	initer2 <- EXP_PLS_det[[2]] %>% mutate(forks=map(forks, function(x) {y <- x %>% filter(!gapovl) %>% mutate(forkdir=case_when(d.Y>0 ~"R",T~"L"));return(y)})) %>% mutate(patIT=map(forks, function(x) paste0(x$forkdir,collapse="")))
	# mapping initiations with the pattern LR
	initer2$patI <- lapply(initer2$patIT, function(x) (do.call(rbind,str_locate_all(x,"LR"))[,1]))
	# mapping termination with the pattern RL
	initer2$patT <- lapply(initer2$patIT, function(x) (do.call(rbind,str_locate_all(x,"RL"))[,1]))
	# positions extraction for initiation
	initer3 <- initer2 %>%
		mutate(ini=map2(forks,patI, function(x,y) {
			if (length(y)!=0)
			{
				x0 <- sapply(y, function(z) x$X0[z])
				x1 <- sapply(y, function(z) x$X0[z+1])
				center <- round((x0+x1)/2)
				spL <- sapply(y, function(z) x$speed[z])
				yL <- sapply(y, function(z) -x$d.Y[z])
				spR <- sapply(y, function(z) x$speed[z+1])
				yR <- sapply(y, function(z) x$d.Y[z+1])
				Ext.center <- x0+spL/(spL+spR)*abs(x1-x0)
				tibble(x0,x1,center,Ext.center,spL,yL,spR,yR)
			}else{
				tibble()
			}})) %>%
		# positions extraction for termination
		mutate(ter=map2(forks,patT, function(x,y) {
			if (length(y)!=0)
			{
				x0=sapply(y, function(z) x$X1[z])
				x1=sapply(y, function(z) x$X1[z+1])
				center=round((x0+x1)/2)
				spR <- sapply(y, function(z) x$speed[z])
				yR <- sapply(y, function(z) x$d.Y[z])
				spL <- sapply(y, function(z) x$speed[z+1])
				yL <- sapply(y, function(z) -x$d.Y[z+1])
				Ext.center <- x0+spL/(spL+spR)*abs(x1-x0)
				tibble(x0,x1,center,Ext.center,spL,yL,spR,yR)
			}else{
				tibble()
			}})) %>%
		select(chrom,strand,read_id,ini,ter)
	# computing speed and dY ratio for initiations
	if (sum(lengths(initer3$ini))>0) {
		iniforks <- initer3 %>%
			select(chrom,strand,read_id,ini) %>%
			unnest(cols=c(ini)) %>%
			mutate(sp.ratio=exp(abs(log(spL/spR)))) %>%
			mutate(dY.ratio=exp(abs(log(yL/yR)))) %>%
			mutate(type="Init")
	}else{
		iniforks <- tibble(chrom=character(),strand=character(),read_id=character(),x0=integer(),x1=integer(),center=double(),spL=double(),yL=double(),spR=double(),yR=double(),sp.ratio=double(),dY.ratio=double(),type=character())
	}
	# computing speed and dY ratio for terminations
	if (sum(lengths(initer3$ter))>0) {
		terforks <- initer3 %>%
			select(chrom,strand,read_id,ter) %>%
			unnest(cols=c(ter)) %>%
			mutate(sp.ratio=exp(abs(log(spL/spR)))) %>%
			mutate(dY.ratio=exp(abs(log(yL/yR))))%>%
			mutate(type="Ter")
	}else{
		terforks <- tibble(chrom=character(),strand=character(),read_id=character(),x0=integer(),x1=integer(),center=double(),spL=double(),yL=double(),spR=double(),yR=double(),sp.ratio=double(),dY.ratio=double(),type=character())
	}

	# merging initiation and termination tibbles
	initer_res <- bind_rows(iniforks,terforks)
	### add leading/lagging ration for speed and dY
	initer_res <- initer_res %>%
		mutate(sp.ratio.LeadLag=pmap_dbl(.,function(strand,spR,spL,...) case_when(strand=="+"~(spR/spL),T~(spL/spR)))) %>%
		mutate(dY.ratio.LeadLag=pmap_dbl(.,function(strand,yR,yL,...) case_when(strand=="+"~(yR/yL),T~(yL/yR))))

	# compute stats for initiations and terminations
	initer_stat <- c(nrow(iniforks),median(iniforks$sp.ratio,na.rm=T),median(iniforks$dY.ratio,na.rm=T),nrow(terforks),median(terforks$sp.ratio,na.rm=T),median(terforks$dY.ratio,na.rm=T))
	names(initer_stat) <- c("nb_init","med(init_sp_ratio)","med(init_dY_ratio)","nb_ter","med(ter_sp_ratio)","med(ter_dY_ratio)")

	AllEXPstats <- as.data.frame(t(c(EXP_stat,EXPlen,EXP_med_read_len,EXP_b2a,EXP_PLSstat,EXP_med_speed,EXP_med_dY,initer_stat)))

	# compute forks, initiations and termination densities per Mb sequenced reads
	forkdens <- EXP_PLSstat[4]/EXPlen*1e6
	names(forkdens) <- "forks_per_Mb"
	initdens <- initer_stat[1]/EXPlen*1e6
	terdens <- initer_stat[4]/EXPlen*1e6
	names(initdens) <- "init_per_Mb"
	names(terdens) <- "ter_per_Mb"
	dens <- cbind(forkdens,initdens,terdens)
	rownames(dens) <- 1
	AllEXPstats2 <- cbind(AllEXPstats,dens)

	# exporting results as a list containing 4 places:
	# 1- reads with 3 or more PLS segments and read with forks after ter and gap filtering
	# 2- forks
	# 3- Initiations and Terminations
	# 4- stats
	res <- list(EXP_PLS_det,EXPforks,initer_res,AllEXPstats2)
	names(res) <- c("PLS_data","forks","initer","stats")

	# saving the data if required
	if (PLS.save==T) {saveRDS(res,file=paste0(EXPname,"_PLS_data.rds"))}

	return(res)
}

## A function to merge PLS results and extract stats

PLS_merging <- function(dir_in,dir_out,ExpName,suff="",file_list0=NA)
{
	require(tidyverse)
	if (is.na(file_list0)) {file_list <- dir(dir_in,pattern=paste0(ExpName,"_"))}else{file_list=file_list0}
	### BEWARE pattern can be misleading
	PLS_reads <- do.call(bind_rows,lapply(file_list, function(x) {
		readRDS(paste0(dir_in,x))[[1]][[2]]})) %>%
		# keeping only the essentiel
		select(-c(RDP.length,sl,patR,patL,fork.R,fork.L,n.forks)) %>%
		arrange(chrom,start)

	PLS_forks <- do.call(bind_rows,lapply(file_list, function(x) {
		readRDS(paste0(dir_in,x))[[2]]})) %>%
		arrange(chrom,st) %>%
		mutate(exp=ExpName)

	PLS_initer <- do.call(bind_rows,lapply(file_list, function(x) {
		readRDS(paste0(dir_in,x))[[3]]}))  %>%
		arrange(chrom,x0,type)
	# redo stats
	Stats_in <- do.call(bind_rows,lapply(file_list, function(x) {
		readRDS(paste0(dir_in,x))[[4]]}))
	PLS_stats <- tibble(Expname=ExpName,n_reads=sum(as.numeric(Stats_in[,2])),n_reads2=sum(as.numeric(Stats_in[,3])),sumlength=sum(as.numeric(Stats_in[,4])),b2a.thr=as.numeric(Stats_in[1,8]),n_reads_RDP3=sum(as.numeric(Stats_in[,11])),n_reads_forks=sum(as.numeric(Stats_in[,12])),n_forks=sum(as.numeric(Stats_in[,14])),speed_med=median(PLS_forks$speed,na.rm=T),dY_med=median(abs(PLS_forks$d.Y),na.rm=T),nb_init=sum(as.numeric(Stats_in[,17])),nb_ter=sum(as.numeric(Stats_in[,20])),fork_dens=sum(as.numeric(Stats_in[,14]))/sum(as.numeric(Stats_in[,4]))*1e6,init_dens=sum(as.numeric(Stats_in[,17]))/sum(as.numeric(Stats_in[,4]))*1e6,ter_dens=sum(as.numeric(Stats_in[,20]))/sum(as.numeric(Stats_in[,4]))*1e6)

	res <- list(PLS_reads,PLS_forks,PLS_initer,PLS_stats)
	names(res) <- c("reads","forks","initer","stats")
	saveRDS(res,file=paste0(dir_out,ExpName,suff,"_PLS_data.rds"))

}



