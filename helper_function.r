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
		pl[[i]] <- ggplot(test$signal[[1]]) +
			geom_text(data=test$sl2[[1]],aes(x=sl.x,y=0,col="RDP_seg_type",label=sl.pat2,fontface="bold"), show.legend = F)+
			geom_line(aes(x=positions,y=signal,col="data.raw"),linetype="dashed",alpha=0.8)+
			geom_line(data=test$RDP[[1]],aes(x=x,y=y,col="RDP_segment"))+
			geom_hline(yintercept=b2a.thr,linetype="dashed") +
			geom_segment(data=test$forks[[1]],aes(x=X1,xend=X2,y=(0.5+sign(as.sl)/40),yend=(0.5+sign(as.sl)/40),col="RDP_fork_chase"),arrow=arrow(length = unit(0.2,"cm")), show.legend = F)+
			geom_segment(data=test$forks[[1]],aes(x=X0,xend=X1,y=(0.5+sign(as.sl)/40),yend=(0.5+sign(as.sl)/40),col="RDP_fork_pulse"),arrow=arrow(length = unit(0.1,"cm")), show.legend = F)+
			geom_text(data=test$forks[[1]],aes(x=(X0+X1)/2,y=(0.8+sign(as.sl)/20),fontface="bold",col="RDP_speed",label=speed.rdp),size=2, show.legend = F)+
			xlab(paste(test$chrom,test$start,test$end,test$strand,test$read_id,sep="_"))+
			guides(col = guide_legend(title = "Legend",override.aes = list(lwd = 1,labels="")))+
			theme(legend.position = "right")+
			scale_color_manual(breaks = c("data.raw","data.smoothed","RDP_segment","RDP_seg_type","RDP_fork_pulse","RDP_fork_chase","RDP_speed"),values = mypal[c(1,2,4,3,6,5,8)])+
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
