#load packages  
library(fda.usc)
library(ggplot2)
library(reshape2)
library(fda)
library(gridExtra)
library(ggthemes)
library(MASS)
library(rstudioapi)
source("../functions/functions.R")

#set working directory 
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# load minimum data -----------------------------------------
data<-read.csv("../data/mintemp.csv")
names(data)

rawtimes<-data[c(1:12)] #data frame of times for each series
rawseries<-data[c(13:24)] #data frame of each series 
names(rawseries)<-c("EHS", "CHS", "WHS", "ES", "CS", "WS", "EC", "CC", "WC", "EH", "CH", "WH")

# html code for gray (for plotting)
gray<-"#bababa"

###### Data Pre-processing ######

# Smoothing 


# create basis function
range<-matrix(NA,nrow=12,ncol=2)
for (i in 1:12){
  range[i,]<-range(rawtimes[,i])
}
start<-max(range[,1])
end<-min(range[,2])


# Since some of the series have missing values, each series is smoothed seperately
fdata<-list(NULL)
fullseries<-matrix(NA, nrow=ncol(rawseries), ncol=100)

for (i in 1:ncol(rawseries)){
  #if there's missing entries in the series
  if (any(is.na(rawseries[,i]))){ 
    rm<-which(is.na(rawseries[,i])) #indices for missing values
    new_time<-rawtimes[-c(rm),i] #remove missing entries from time 
    new_series<-rawseries[-c(rm),i] #remove missing entries from series
    daybasis <- create.bspline.basis(range=range(new_time), norder=3, nbasis=21)
    smooth <- smooth.basis(new_time, new_series, daybasis) #smooth the series
    fdata[[i]]<-fdata(smooth$fd, argvals = seq(from=start, to=end, length.out=100)) #rediscretize data into fdata format
  }
  #if no missing entries 
  else{ 
    new_time<-rawtimes[,i] #time
    new_series<-rawseries[,i] #series 
    daybasis <- create.bspline.basis(range=range(new_time), norder=3, nbasis=21)
    smooth <- smooth.basis(new_time, new_series, daybasis) #smooth the series
    fdata[[i]]<-fdata(smooth$fd, argvals = seq(from=start, to=end, length.out=100)) #rediscretize data into fdata format
  }
  fullseries[i,]<-fdata[[i]]$data #save the re-discretized series post-smoothing 
}

head(fullseries)

# save times corresponding to re-discretized series 
argvals<-fdata[[1]]$argvals 

# set dimnames 
dimnames(fullseries)<-list(names(rawseries), argvals) 

# view the order of series
dimnames(fullseries)[[1]]

####### FANOVA Tests ######

### Test for difference b/w all treatments 
group = c(1,1,1,2,2,2,3,3,3,4,4,4) #group series according to group membership 
N <- 10000 # number of bootstrap resamples 
set.seed(30) 
res <-obtain_resamples(argvals, fullseries, N, group, group_names=c("HSR","S","C","H")) # obtain resample curves 
pval <-obtain_pvalue(res$resamples_mean, res$sample_mean, c(1,1,1,1)) #obtain pvalue 
pval$pvalue #pvalue = 0.0246 - moderate evidence

### Test for pairwise difference 
# HSR - S
treatment<-c(1,1,0,0)
pval <-obtain_pvalue(res$resamples_mean, res$sample_mean, treatment)
pval$pvalue #0.0432
tukey <- tukey_correction(res$resamples_mean, res$sample_mean, treatment)
tukey$pvalue #0.1548

# HSR - C
treatment<-c(1,0,1,0)
pval <-obtain_pvalue(res$resamples_mean, res$sample_mean, treatment)
pval$pvalue #0.0357
tukey <- tukey_correction(res$resamples_mean, res$sample_mean, treatment)
tukey$pvalue #0.0628

# HSR - H 
treatment<-c(1,0,0,1)
pval <-obtain_pvalue(res$resamples_mean, res$sample_mean, treatment)
pval$pvalue #0.2565
tukey <- tukey_correction(res$resamples_mean, res$sample_mean, treatment)
tukey$pvalue #0.4913

# S - C
treatment<-c(0,1,1,0)
pval <-obtain_pvalue(res$resamples_mean, res$sample_mean, treatment)
pval$pvalue #0.1939
tukey <- tukey_correction(res$resamples_mean, res$sample_mean, treatment)
tukey$pvalue #0.7023

# S - H
treatment<-c(0,1,0,1)
pval <-obtain_pvalue(res$resamples_mean, res$sample_mean, treatment)
pval$pvalue #0.0077
tukey <- tukey_correction(res$resamples_mean, res$sample_mean, treatment)
tukey$pvalue #0.2156

# C - H 
treatment<-c(0,0,1,1)
pval <-obtain_pvalue(res$resamples_mean, res$sample_mean, treatment)
pval$pvalue #0.0424
tukey <- tukey_correction(res$resamples_mean, res$sample_mean, treatment)
tukey$pvalue #0.2208

### Test for interaction b/w heating and snowremoval
con <- c(1, -1, 1, -1) #contrast: HSR - S + C -H
pval<-obtain_pvalue_contrast(res$resamples_mean, res$sample_mean, contrast=con) 
pval$pvalue #0.6114

# visualization 
# obtain dataframe of difference curves, sample difference curve, and time
dat<-vis_dataframe(argvals = argvals, resamples_mean = res$resamples_mean, sample_mean = res$sample_mean, contrast = con)
# convert to long data format for ggplot 
datlong<-melt(dat$dataframe, id.vars="argvals")
# plot 
int_plot<-ggplot()+geom_line(data=datlong, aes(x=as.Date(argvals, origin="2011-01-01 00:00:00"), y=value, group=variable), colour=gray, alpha=0.2)+
  geom_line(data=subset(datlong, variable == "res.1"|variable == "sample"),aes(x=as.Date(argvals, origin="2011-01-01 00:00:00"),y=value, col=variable))+
  scale_color_manual(labels = c("Resamples", "Sample"), values=c(gray,'black'))+
  #scale_color_manual(labels = c("Resamples",bquote(paste("Sample(",tilde(H),"S)"))), values=c('black', gray))+
  theme_classic() + theme(legend.title = element_blank())+
  labs(title="",x ="Time", y = "Min Temperature Interaction", col="")+timescale+mytheme+
  annotate(geom="text", x=as.Date("2011-09-27"), y=Inf,label=paste("global p-value=",pval$pvalue),hjust = 1, vjust = 1)

### Test for main effect (heating)
group = c(1,1,1,2,2,2,2,2,2,1,1,1) #group 1 = HSR, H; group 2 = S,C 
N <- 10000 # number of bootstrap resamples 
set.seed(30) 
res <-obtain_resamples(argvals, fullseries, N, group, group_names=c("H","NH"))
pval <-obtain_pvalue(res$resamples_mean, res$sample_mean, c(1,1)) #obtain pvalue 
pval$pvalue # 0.0047 

# visualization 
# obtain dataframe of difference curves, sample difference curve, and time
dat<-vis_dataframe(argvals = argvals, resamples_mean = res$resamples_mean, sample_mean = res$sample_mean, contrast = con)
# convert to long data format for ggplot 
datlong<-melt(dat$dataframe, id.vars="argvals")
# plot 
main_plot1<-ggplot()+geom_line(data=datlong, aes(x=as.Date(argvals, origin="2011-01-01 00:00:00"), y=value, group=variable), colour=gray, alpha=0.2)+
  geom_line(data=subset(datlong, variable == "res.1"|variable == "sample"),aes(x=as.Date(argvals, origin="2011-01-01 00:00:00"),y=value, col=variable))+
  scale_color_manual(labels = c("Resamples", "Sample"), values=c(gray,'black'))+
  #scale_color_manual(labels = c("Resamples",bquote(paste("Sample(",tilde(H),"S)"))), values=c('black', gray))+
  theme(legend.title = element_blank())+
  labs(title="",x ="Time", y = bquote(paste("Min Temperature Difference (",tilde(H)," - N", tilde(H), ")", sep="")), col="")+timescale+mytheme+
  annotate(geom="text", x=as.Date("2011-09-27"), y=Inf,label=paste("global p-value=",pval$pvalue), hjust = 1, vjust = 1)

# westfall young correction (to obtain significance over time)
# create split in domain 
#split<-20
#argvalsindex<-split(seq(1:length(argvals)), ceiling(seq_along(1:length(argvals))/5))
argvalsindex<-list(seq(from=1, to=5), seq(from=5, to =10), seq(from=10, to=15), seq(from=15, to=20),
                   seq(from=20, to=25), seq(from=25, to =30), seq(from=30, to=35), seq(from=35, to=40),
                   seq(from=40, to=45), seq(from=45, to =50), seq(from=50, to=55), seq(from=55, to=60),
                   seq(from=60, to=65), seq(from=65, to =70), seq(from=70, to=75), seq(from=75, to=80),
                   seq(from=80, to=85), seq(from=85, to =80), seq(from=90, to=95), seq(from=95, to=100))
wpval<-westfall_young_pvalue(argvals=argvals, resamples_mean=res$resamples_mean, sample_mean=res$sample_mean, group=group, split=argvalsindex)
wpval$raw_pvalue
wpval$corr_pvalue

#determine cutoff 
sigdomain1<-NULL
for (i in which(wpval$corr_pvalue<=0.05)){
  sig1<-argvalsindex[[i]]
  sigdomain1<-c(sigdomain1,sig1)
}

sigdomain2<-NULL
for (i in which(wpval$corr_pvalue<=0.10 & wpval$corr_pvalue>0.05)){
  sig2<-argvalsindex[[i]]
  sigdomain2<-c(sigdomain2,sig2)
}


sig1<-data.frame(time=as.Date(argvals[sigdomain1],origin="2011-01-01 00:00:00"),sig_diff=dat$dataframe$sample[sigdomain1])
sig2<-data.frame(time=as.Date(argvals[sigdomain2],origin="2011-01-01 00:00:00"),sig_diff=dat$dataframe$sample[sigdomain2])

corrected_main_plot1<-main_plot1+geom_line(data=sig1, aes(x=time, y=sig_diff), colour="red", alpha=1, size=1)+geom_line(data=sig2, aes(x=time, y=sig_diff), colour="yellow", alpha=1, size=1)


### Test for main effect (snow removal)
group = c(1,1,1,1,1,1,2,2,2,2,2,2) #group 1 = HSR, S; group 2 = H,C
N <- 10000 # number of bootstrap resamples 
set.seed(30) 
res <-obtain_resamples(argvals, fullseries, N, group, group_names=c("S","NS"))
pval <-obtain_pvalue(res$resamples_mean, res$sample_mean, c(1,1)) #obtain pvalue 
pval$pvalue # 0.207

# visualization 
# obtain dataframe of difference curves, sample difference curve, and time
dat<-vis_dataframe(argvals = argvals, resamples_mean = res$resamples_mean, sample_mean = res$sample_mean, contrast = con)

# convert to long data format for ggplot 
datlong<-melt(dat$dataframe, id.vars="argvals")
# plot 
main_plot2<-ggplot()+geom_line(data=datlong, aes(x=as.Date(argvals, origin="2011-01-01 00:00:00"), y=value, group=variable), colour=gray, alpha=0.2)+
  geom_line(data=subset(datlong, variable == "res.1"|variable == "sample"),aes(x=as.Date(argvals, origin="2011-01-01 00:00:00"),y=value, col=variable))+
  scale_color_manual(labels = c("Resamples", "Sample"), values=c(gray,'black'))+
  #scale_color_manual(labels = c("Resamples",bquote(paste("Sample(",tilde(H),"S)"))), values=c('black', gray))+
  theme(legend.title = element_blank())+
  labs(title="",x ="Time", y = bquote(paste("Min Temperature Difference (",tilde(S)," - N", tilde(S), ")", sep="")), col="")+timescale+mytheme+
  annotate(geom="text", x=as.Date("2011-09-27"), y=Inf,label=paste("global p-value=",pval$pvalue),hjust = 1, vjust = 1)

# westfall young correction (to obtain significance over time)
# create split in domain 
#split<-20
#argvalsindex<-split(seq(1:length(argvals)), ceiling(seq_along(1:length(argvals))/5))
argvalsindex<-list(seq(from=1, to=5), seq(from=5, to =10), seq(from=10, to=15), seq(from=15, to=20),
                   seq(from=20, to=25), seq(from=25, to =30), seq(from=30, to=35), seq(from=35, to=40),
                   seq(from=40, to=45), seq(from=45, to =50), seq(from=50, to=55), seq(from=55, to=60),
                   seq(from=60, to=65), seq(from=65, to =70), seq(from=70, to=75), seq(from=75, to=80),
                   seq(from=80, to=85), seq(from=85, to =80), seq(from=90, to=95), seq(from=95, to=100))
wpval<-westfall_young_pvalue(argvals=argvals, resamples_mean=res$resamples_mean, sample_mean=res$sample_mean, group=group, split=argvalsindex)
wpval$raw_pvalue
wpval$corr_pvalue

#determine cutoff 
sigdomain1<-NULL
for (i in which(wpval$corr_pvalue<=0.05)){
  sig1<-argvalsindex[[i]]
  sigdomain1<-c(sigdomain1,sig1)
}

sigdomain2<-NULL
for (i in which(wpval$corr_pvalue<=0.10 & wpval$corr_pvalue>0.05)){
  sig2<-argvalsindex[[i]]
  sigdomain2<-c(sigdomain2,sig2)
}


sig1<-data.frame(time=as.Date(argvals[sigdomain1],origin="2011-01-01 00:00:00"),sig_diff=dat$dataframe$sample[sigdomain1])
sig2<-data.frame(time=as.Date(argvals[sigdomain2],origin="2011-01-01 00:00:00"),sig_diff=dat$dataframe$sample[sigdomain2])

corrected_main_plot2<-main_plot2+geom_line(data=sig1, aes(x=time, y=sig_diff), colour="red", alpha=1, size=1)+geom_line(data=sig2, aes(x=time, y=sig_diff), colour="yellow", alpha=1, size=1)


