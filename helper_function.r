### Helper Function used in Theulot et al 2021

#### LL 20210917

###
simpleRFD <- function(gr,lr=1,na2zero=F,expor=F,outname='myRFDdata')
{
	### L=+ to keep compatibility with OK seq data
	require(GenomicRanges)
	require(rtracklayer)

	bs <- 1

	cv_L <- coverage(gr[strand(gr)=='+'])
	cv_R <- coverage(gr[strand(gr)=='-'])

	cv <- cv_L+cv_R
	RFD <- (cv_R-cv_L)/(cv_R+cv_L)
	lr_index <- which(cv<=lr)
	RFD2 <- RFD
	RFD2[lr_index] <- NA

	naname <- '_wiNA'
	if (na2zero)
	{
		RFD[is.na(RFD)] <- 0
		RFD2[is.na(RFD2)] <- 0
		naname <- '_noNA'
	}

	if (expor)
	{
		export(cv,con=paste0(outname,'_cov_tot_bs',bs/1000,'k_lr',lr,'.bw'))
		export(cv_L,con=paste0(outname,'_cov_2left_bs',bs/1000,'k_lr',lr,'.bw'))
		export(cv_R,con=paste0(outname,'_cov_2right_bs',bs/1000,'k_lr',lr,'.bw'))
		export(RFD2,con=paste0(outname,'_RFD_bs',bs/1000,'k_lr',lr,naname,'.bw'))
	}
	res <- list(cv,cv_L,cv_R,RFD,RFD2)
	names(res) <- c('cv','cv_L','cv_R','RFD','RFD2')
	return(res)
}

##

### a function to check correlation between RFD (or other coverage like type of data)
cor.rfd <- function(a,b,met='s')
{cor(as.numeric(unlist(a)[!is.na(unlist(a)) & !is.na(unlist(b))]),as.numeric(unlist(b)[!is.na(unlist(a)) & !is.na(unlist(b))]),method=met)}
##

### a function to plot forks with RDP informations
plotforks <- function(toto,b2a.thr=0.02,fileout,plot.smooth=T)
{

	suppressMessages(require(tidyverse))
	require(gridExtra)
	require(RcppRoll)
	theme_set(theme_bw())

	mypal=RColorBrewer::brewer.pal(12,"Paired")
	pl=list()
	for (i in 1:nrow(toto))
	{
		test <- toto %>% dplyr::slice(i)
		pl[[i]] <- ggplot(test$signalb[[1]]) +
			geom_text(data=test$sl2[[1]],aes(x=sl.x,y=0,col="RDP_seg_type",label=sl.pat2,fontface="bold"), show.legend = F)+
			geom_line(aes(x=positions,y=signal,col="data.raw"),linetype="dashed",alpha=0.8)+
			geom_line(data=test$RDP[[1]],aes(x=x,y=y,col="RDP_segment"))+
			geom_hline(yintercept=b2a.thr,linetype="dashed") +
			geom_segment(data=test$forks[[1]],aes(x=X1,xend=X2,y=(0.5+sign(d.Y)/40),yend=(0.5+sign(d.Y)/40),col="PLS_fork_chase"),arrow=arrow(length = unit(0.2,"cm")), show.legend = F)+
			geom_segment(data=test$forks[[1]],aes(x=X0,xend=X1,y=(0.5+sign(d.Y)/40),yend=(0.5+sign(d.Y)/40),col="PLS_fork_pulse"),arrow=arrow(length = unit(0.1,"cm")), show.legend = F)+
			geom_text(data=test$forks[[1]],aes(x=(X0+X1)/2,y=(0.8+sign(d.Y)/20),fontface="bold",col="PLS_speed",label=speed.rdp),size=2, show.legend = F)+
			xlab(paste(test$chrom,test$start,test$end,test$strand,test$read_id,sep="_"))+
			guides(col = guide_legend(title = "Legend",override.aes = list(lwd = 1,labels="")))+
			theme(legend.position = "right")+
			scale_color_manual(breaks = c("data.raw","data.smoothed","RDP_segment","RDP_seg_type","PLS_fork_pulse","PLS_fork_chase","PLS_speed"),values = mypal[c(1,2,4,3,6,5,8)])+
			coord_cartesian(ylim=c(0,1))

		if (plot.smooth) {
			w <- w_pg25
			pl[[i]] <- pl[[i]]+
				geom_line(aes(x=positions,y=roll_mean(signal,weights=w,normalize=T,na.rm=T,fill=NA),col="data.smoothed"))
		}
	}
	pdf(fileout,height=12)
	if (nrow(toto)>=5)
	{for (j in seq(1,(nrow(toto)),5)[1:(nrow(toto)%/%5)])
	{do.call(grid.arrange,c(pl[j:(j+4)],ncol=1))}
	}
	if (nrow(toto)%%5 >0)
	{
		j=tail(seq(1,(nrow(toto)),5),1)
		do.call(grid.arrange,c(pl[j:(j+nrow(toto)%%5-1)],ncol=1,nrow=5))
	}
	dev.off()
}

##

### plot the distribution of the signal to set the b2a.thr

plot_signal <- function(EXP,xmax=1,EXPname="EXP",bs=1000,minlen=5000,EXP_b2a.thr0=0.02,alldata=F,nreads=NA)
{
	suppressMessages(require(tidyverse))
	require(ggpubr)
	theme_set(theme_bw())


	if (!is.na(nreads) & nrow(EXP)>nreads)
	{
		EXP2 <- sample_n(EXP,nreads)
	}else{
		EXP2 <- EXP
	}

	EXP_PLSall <- EXP2 %>%
		mutate(read_id=map_chr(read_id, function(x) str_remove(x,"read_"))) %>%
		select(read_id,chrom,start,end,strand,signalb) %>%
		mutate(length=end-start) %>%
		filter(length>minlen) %>%
		mutate(signalb=map(signalb, function(x) {x %>% dplyr::rename(signal=2)})) %>%
		mutate(RDP=map(signalb,myRDP0,RDP.eps=0.15)) %>%
		mutate(RDP.length=map_int(RDP,function(x) nrow(x))) %>%
		mutate(Bmedy=map_dbl(signalb,function(z) median(z$signal)))
	EXP_PLS3 <- EXP_PLSall %>% filter(RDP.length>3)
	# all data
	if (alldata==T) {
		test0 <- EXP_PLSall %>%
			mutate(noise= map(signalb, function(y) {
				y %>%
					mutate(positions = round(positions/bs)*bs) %>%
					group_by(positions) %>%
					summarise(Bmean = mean(signal),.groups = 'drop')%>%
					select(Bmean)
			})) %>%
			select(noise) %>%
			unnest(cols=c(noise))
		signal_plot0 <- ggplot(test0)+
			geom_histogram(aes(x=Bmean),binwidth=0.002,alpha=0.3)+
			geom_vline(aes(xintercept=EXP_b2a.thr0))+
			coord_cartesian(xlim=c(0,xmax))+
			scale_x_continuous(paste0("mean B signal by ",bs/1000,"kb"), breaks=seq(0,xmax,0.1))

		signal_plot1 <- ggplot(test0 %>% filter(Bmean>0.002))+
			geom_histogram(aes(x=Bmean),binwidth=0.002,alpha=0.3)+
			geom_vline(aes(xintercept=EXP_b2a.thr0))+
			coord_cartesian(xlim=c(0.002,xmax/4))+
			scale_x_continuous(paste0("mean B signal by ",bs/1000,"kb"), breaks=seq(0,xmax/4,0.01))+
			theme(axis.text.x = element_text(angle = 45,hjust=1))
		ggarrange(signal_plot0,signal_plot1,nrow=2)
		ggsave(paste0(EXPname,"_all_1kbmeansignal.pdf"),h=8,w=6)
	}
	# RDP0>3 data
	test1 <- EXP_PLS3 %>%
		mutate(noise= map(signalb, function(y) {
			y %>%
				mutate(positions = round(positions/bs)*bs) %>%
				group_by(positions) %>%
				summarise(Bmean = mean(signal),.groups = 'drop')%>%
				select(Bmean)
		}
		)) %>%
		select(noise) %>% unnest(cols=c(noise))
	signal_plot2 <- ggplot(test1)+
		geom_histogram(aes(x=Bmean),binwidth=0.002,alpha=0.3)+
		geom_vline(aes(xintercept=EXP_b2a.thr0))+
		coord_cartesian(xlim=c(0,xmax))+
		scale_x_continuous(paste0("mean B signal by ",bs/1000,"kb"), breaks=seq(0,xmax,0.1))

	signal_plot3 <- ggplot(test1 %>% filter(Bmean>0.002))+
		geom_histogram(aes(x=Bmean),binwidth=0.002,alpha=0.3)+
		geom_vline(aes(xintercept=EXP_b2a.thr0))+
		coord_cartesian(xlim=c(0.002,xmax/4))+
		scale_x_continuous(paste0("mean B signal by ",bs/1000,"kb"), breaks=seq(0,xmax/4,0.01))+
		theme(axis.text.x = element_text(angle = 45,hjust=1))

	ggarrange(signal_plot2,signal_plot3,nrow=2)

	ggsave(paste0(EXPname,"_RDP3_1kbmeansignal.pdf"),h=8,w=6)

}

### my GR shuffling
shuffleGR4=function(seqinf=seqinfS288CrDNA,chrnb=16,inputGR=inputData,gap=Ngaps2)
{	require(GenomicRanges)
	seqname=seqnames(seqinf)

	hit <- inputGR[seqnames(inputGR)==seqname[chrnb]]
	gapchr=gap[seqnames(gap)==seqname[chrnb]]
	# altenative to deal with no gap
	if (length(gapchr)==0) {gapchr=GRanges(seqnames=seqname[chrnb],ranges=IRanges(start=1,width=1),seqinfo=seqinfo(inputGR))}
	ravail <- ranges(gaps(gapchr)[seqnames(gaps(gapchr))==seqname[chrnb] & strand(gaps(gapchr))=="*"])
	#		st_avail <- unlist(as.vector(ravail))
	# broken in BioC3.7, should come back in BioC3.8
	# Temporary fix
	st_avail <- IRanges:::unlist_as_integer(ravail)
	#
	st_rdgr <- sample(st_avail,length(hit))
	if (length(hit)==1)
	{
		wi_rdgr <- width(hit)
	}else{
		wi_rdgr <- sample(width(hit))
		#necessary if only one range sample(width()) choose a number
		#betwen in 1:width() rather than one width
	}
	ra_rdgr <- sort(IRanges(start=st_rdgr,width=wi_rdgr))
	rgap <- ranges(gapchr)
	#sum(overlapsAny(ra_rdgr,ranges(gapchr)))

	keep <- IRanges()
	ra_rdgr2 <- IRanges()
	while ((sum(overlapsAny(ra_rdgr,rgap))!=0) | (sum(overlapsAny(ra_rdgr2,keep))!=0))
	{
		keep <- ra_rdgr[overlapsAny(ra_rdgr,rgap)==0]
		hit2 <- ra_rdgr[overlapsAny(ra_rdgr,rgap)!=0]
		st_rdgr2 <- sample(st_avail,length(hit2))
		if (length(hit2)==1)
		{
			wi_rdgr2 <- width(hit2)
		}else{
			wi_rdgr2 <- sample(width(hit2))
		}
		ra_rdgr2 <- IRanges(start=st_rdgr2,width=wi_rdgr2)
		ra_rdgr <- c(keep,ra_rdgr2)
	}
	rdgr <- sort(GRanges(seqnames=Rle(rep(seqname[chrnb],length(hit))),ranges=ra_rdgr,strand=Rle(rep('*',length(hit))),seqinfo=seqinfo(inputGR)))
	return(rdgr)
}

# function to resample on a genome

shuffleGRgen <- function(dummy=1,seqinf2=seqinfS288CrDNA,inputGR2=inputData,gap2=Ngaps2,chrlist=1:chnb)
{
	rdlist=GRangesList()
	for (i in chrlist) {rdlist[[i]] <- shuffleGR4(seqinf=seqinf2,chrnb=i,inputGR=inputGR2,gap=gap2)}
	y<- do.call(c,rdlist)
	return(y)
}

# Gap annotation
findNgaps <- function(x)
	# x is a DNAString
{ y=Rle(strsplit(as.character(x),NULL)[[1]])
y2=ranges(Views(y,y=='N'))
return(y2)	# y2 is a list of IRanges
}

### a function to change seqinf of a GRanges
NewSeqinfo <- function(GR,seqin) {
	seqlevels(GR,pruning.mode="coarse") <- seqlevels(seqin)
	seqinfo(GR) <- seqin
	return(GR)
}

function.cluster2 <- function(input,clust.dist0=clust.dist)
	# input is a tibble of genomic objects with a width of 1
{
	require(tidygenomics)
	test <- genome_cluster(input, by=c("chrom", "center", "center"),max_distance=clust.dist0)
	test2 <- test %>%
		group_by(cluster_id) %>%
		nest %>%
		mutate(eff=map_int(data,nrow)) %>%
		mutate(cl.med=map_dbl(data,function(x) median(x$center))) %>%
		mutate(cl.start=map_dbl(data,function(x) min(x$center)))%>%
		mutate(cl.end=map_dbl(data,function(x) max(x$center)))%>%
		mutate(cl.mad=map_dbl(data,function(x) mad(x$center)))%>%
		unnest(cols=c(data)) %>%
		select(-c(strand,spL,yL,spR,yR,sp.ratio,dY.ratio,type,sp.ratio.LeadLag,dY.ratio.LeadLag,exp))
	test3 <- test2 %>% select(chrom,cl.med,cl.start,cl.end,eff,cl.mad,cluster_id) %>% group_by(cluster_id) %>% slice(1)
	res <- list(test2,test3)
	return(res)
}


### a function to compute mean trace with padding

padandtrace <- function(EXPforks,pad.max=500)
{

	suppressMessages(require(tidyverse))

	position <- (-10:(pad.max-1))*100
	EXPtrace0f <- EXPforks %>%
		select(trac,exp,len) %>%
		mutate(sig=map(trac, function(x) x%>% pull(signal))) %>%
		mutate(sig2=map(sig, function(x) {if(length(x)<pad.max){res=c(x,rep(NA,(pad.max+10-length(x))))}else{res=x[1:(pad.max+10)]};return(res)}))
	EXPtrace1mean <- colMeans(do.call(rbind,EXPtrace0f$sig2),na.rm=T)
	out <- suppressMessages(bind_cols(position,EXPtrace1mean))
	names(out) <- c("position","mean_trace")
	return(out)
}

compute_meantrace <- function(PLStib,trac.xmax=50000)
{
	toplot <- PLStib
	pad.max0 <- trac.xmax/100

	totrace <- split(toplot, toplot$exp)
	totrace2 <- lapply(totrace, function(x) padandtrace(x,pad.max0))
	totrace3 <- lapply(totrace2, function(x) do.call(bind_cols,x))
	totrace4 <- tibble(exp=factor(names(totrace3),levels=levels(toplot$exp)),trace=totrace3)
	return(totrace4)
}

plot_meantrace <- function(totrace4,explist,ymax0=0.8,xmax=50000,normalise=F,root.title="",pathout="",nameout="EXP_meantrace",expor=T)
{
	totrace1 <- totrace4 %>% filter(exp %in% explist)
	if (normalise) {
		totrace1 <- totrace1 %>%
			mutate(trace=map(trace, function(x) {x %>% mutate(mean_trace=mean_trace/max(mean_trace,na.rm=T))}))
		ymax=1
		lab.y="%BrdU (norm)"
	}else{
		ymax=ymax0
		lab.y="%BrdU"
	}

	totrace5 <- totrace1 %>% unnest(cols=c(trace))
	p1 <- ggplot(totrace5) +
		geom_line(aes(x=position, y=mean_trace,col=exp))+
		coord_cartesian(xlim=c(-1000,xmax),ylim=c(0,ymax))+
		paletteer::scale_color_paletteer_d("ggthemes::Classic_20")+
		ylab(lab.y)+
		xlab("Position (relative to X0)")+
		ggtitle(paste0(root_title,"_",nameout))
	if(expor)
	{
		ggsave(paste0(pathout,root_title,"_",nameout,".pdf"), h=6,w=8)
	}else{p1}
}

### plot read length
plot_readlength <- function(EXP,EXP_PLS,fileout=NA,ymax=150000) {
	if (is.na(fileout))
	{
		fileout <- EXP.PLS[[2]]$exp[1]
	}
	toplot <- bind_rows(
		tibble(len=EXP %>% mutate(length=end-start) %>% filter(length>5000) %>% pull(length),leg="All reads >5kb"),
		tibble(len=EXP_PLS[[1]][[1]] %>% pull(length),leg="All reads RDP3"),
		tibble(len=EXP_PLS[[1]][[2]] %>% pull(length),leg="Reads with forks"))
	totext <- toplot %>% group_by(leg) %>% summarise(n=n()) %>% ungroup
	tomed <- toplot %>% group_by(leg) %>% summarise(med=round(median(len),0)) %>% ungroup
	ggplot(toplot)+
		geom_violin(aes(y=len,fill=leg,x=leg),col=NA,scale="width")+
		coord_cartesian(ylim=c(0,ymax))+
		geom_boxplot(aes(y=len,x=leg),outlier.shape=NA,width=0.2)+
		geom_text(data=totext,aes(x=leg,y=0,label=n),fontface="italic") +
		geom_text(data=tomed,aes(x=leg,ymax-1000,label=med),col="red") +
		scale_fill_brewer(palette="Set1")+
		ggtitle(fileout)+
		xlab("Read categories")+
		ylab("Length")
	ggsave(paste0(fileout,"_readlength.pdf"),h=4,w=6)
}

