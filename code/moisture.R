#load packages  
library(fda.usc)
library(ggplot2)
library(reshape2)
library(fda)
library(gridExtra)
library(ggthemes)
library(MASS)
library(rstudioapi)

#set working directory 
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# read in fanova functions from file
source("../functions/functions.R")

# html code for gray (for plotting)
gray<-"#bababa"
timescale<-scale_x_date(date_breaks = "2 weeks", date_labels =  "%b %d")
mytheme=theme_classic()+theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1))

# read in data
data<-read.csv("../data/moisture.csv")
names(data)
time<-data[,1] #date
rawseries<-as.matrix(data[c(2:13)]) #full series 
str(data)
datalong<-melt(data, id.vars="day2011")
unique(datalong$variable)
datalong$treatment[datalong$variable=="EHS_25cm_M" | datalong$variable=="CHS_25cm_M" | datalong$variable=="WHS_25cm_M"]<-"HSR"
datalong$treatment[datalong$variable=="ES_25cm_M" | datalong$variable=="CS_25cm_M" | datalong$variable=="WS_25cm_M"]<-"SR"
datalong$treatment[datalong$variable=="EC_25cm_M" | datalong$variable=="CC_25cm_M" | datalong$variable=="WC_25cm_M"]<-"C"
datalong$treatment[datalong$variable=="EH_25cm_M" | datalong$variable=="CH_25cm_M" | datalong$variable=="WH_25cm_M"]<-"H"
ggplot(data=datalong, aes(x=day2011, y=value, group=variable, color = treatment))+geom_line()

###### Data Pre-processing ######

# Smoothing 
# breaks control placement of knots 
breaks<-c(seq(from=range(time)[1], to=190,length.out=5),
          seq(from=193, to=230, length.out=10),
          seq(from=230, to=range(time)[2], length.out=5))

# create basis function
daybasis <- create.bspline.basis(range(time), norder=3, nbasis=21, breaks=breaks) 

# Since some of the series have missing values, each series is smoothed seperately
fdata<-list(NULL)
fullseries<-matrix(NA, nrow=ncol(rawseries), ncol=100)

for (i in 1:ncol(rawseries)){
  #if there's missing entries in the series
  if (any(is.na(rawseries[,i]))){ 
    rm<-which(is.na(rawseries[,i])) #indices for missing values
    new_time<-time[-c(rm)] #remove missing entries from time 
    new_series<-rawseries[-c(rm),i] #remove missing entries from series
    smooth <- smooth.basis(new_time, new_series, daybasis) #smooth the series
    fdata[[i]]<-fdata(smooth$fd, argvals = seq(from=range(new_time)[1], to=range(new_time)[2], length.out=100)) #rediscretize data into fdata format
  }
  #if no missing entries 
  else{ 
    new_time<-time #time
    new_series<-rawseries[,i] #series 
    smooth <- smooth.basis(new_time, new_series, daybasis) #smooth the series
    fdata[[i]]<-fdata(smooth$fd, argvals = seq(from=range(new_time)[1], to=range(new_time)[2], length.out=100)) #rediscretize data into fdata format
  }
  fullseries[i,]<-fdata[[i]]$data #save the re-discretized series post-smoothing 
}

# save times corresponding to re-discretized series 
argvals<-fdata[[1]]$argvals 
# set dimnames 
dimnames(fullseries)<-list(names(data[,-1]), argvals) 

# view the order of series
dimnames(fullseries)[[1]]

####### FANOVA Tests ######

### Test for difference b/w all treatments 
group = c(1,1,1,2,2,2,3,3,3,4,4,4) #group series according to group membership 
N <- 10000 # number of bootstrap resamples 
set.seed(30) 
res <-obtain_resamples(argvals, fullseries, N, group, group_labels=c("HSR","S","C","H")) # obtain resample curves 
pval <-obtain_pvalue(res$resamples_mean, res$sample_mean, include_group=c(1,2,3,4)) #obtain pvalue 
pval$pvalue #pvalue = 0.0854 - moderate evidence

### Test for pairwise difference 
# C - H
include_group<-c(1,2) # which groups to include in test 
pval <-obtain_pvalue(res$resamples_mean, res$sample_mean, include_group) # un-corrected pvalue 
print(c(pval$group_labels, pval$pvalue)) #0.3621
tukey <- tukey_correction(res$resamples_mean, res$sample_mean, include_group) # tukey corrected pvalue 
print(c(tukey$group_labels, tukey$pvalue)) #0.7637

# obtain tukey pvals for all group pairs 
pairs<-combn(c(1,2,3,4), 2) # obtain indices for all group pairing 
pairwise_pvals<-matrix(NA, nrow=dim(pairs)[2], ncol=3) # empty matrix to store pvalues 

for (i in 1:nrow(pairwise_pvals)){
  include_group<-c(pairs[,i]) # indices for pairing
  tukey <- tukey_correction(res$resamples_mean, res$sample_mean, include_group) # obtain tukey correction
  pairwise_pvals[i,] <- c(tukey$group_labels,tukey$pvalue) # save group name & pvalue
}
pairwise_pvals<-data.frame(pairwise_pvals)
names(pairwise_pvals)<-c("group 1", "group 2", "p-values")
pairwise_pvals


### Test for interaction b/w heating and snowremoval
con <- c(1, -1, 1, -1) #contrast: HSR - S + C -H
pval<-obtain_pvalue_contrast(res$resamples_mean, res$sample_mean, contrast=con) 
pval$pvalue #0.9327

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
  labs(title="",x ="Time", y = "Soil Moisture Interaction", col="")+timescale+mytheme+
  annotate(geom="text", x=as.Date("2011-09-27"), y=Inf,label=paste("global p-value=",pval$pvalue),hjust = 1, vjust = 1)

### Test for main effect (heating)
group = c(1,1,1,2,2,2,2,2,2,1,1,1) #group 1 = HSR, H; group 2 = S,C 
N <- 10000 # number of bootstrap resamples 
set.seed(30) 
res <-obtain_resamples(argvals, fullseries, N, group, group_labels=c("H","NH"))

pval <-obtain_pvalue(res$resamples_mean, res$sample_mean, include_group = c(1,2)) #obtain pvalue 
pval$pvalue # 0.2962 

# visualization 
# obtain dataframe of difference curves, sample difference curve, and time
#dat<-vis_dataframe(argvals = as.Date(argvals, origin="2011-01-01 00:00:00"), resamples_mean = res$resamples_mean, sample_mean = res$sample_mean, contrast = con)
dat<-vis_dataframe(argvals = argvals, resamples_mean = res$resamples_mean, sample_mean = res$sample_mean, contrast = con)
# convert to long data format for ggplot 
datlong<-melt(dat$dataframe, id.vars="argvals")
# plot 
main_plot1<-ggplot()+geom_line(data=datlong, aes(x=as.Date(argvals, origin="2011-01-01 00:00:00"), y=value, group=variable), colour=gray, alpha=0.2)+
  geom_line(data=subset(datlong, variable == "res.1"|variable == "sample"),aes(x=as.Date(argvals, origin="2011-01-01 00:00:00"),y=value, col=variable))+
  scale_color_manual(labels = c("Resamples", "Sample"), values=c(gray,'black'))+
  #scale_color_manual(labels = c("Resamples",bquote(paste("Sample(",tilde(H),"S)"))), values=c('black', gray))+
  theme(legend.title = element_blank())+
  labs(title="",x ="Time", y = bquote(paste("Soil Moisture Difference (",tilde(H)," - N", tilde(H), ")", sep="")), col="")+timescale+mytheme+
  annotate(geom="text", x=as.Date("2011-09-27"), y=Inf,label=paste("global p-value=",pval$pvalue), hjust = 1, vjust = 1)

### Test for main effect (snow removal)
group = c(1,1,1,1,1,1,2,2,2,2,2,2) #group 1 = HSR, S; group 2 = H,C
N <- 10000 # number of bootstrap resamples 
set.seed(30) 
res <-obtain_resamples(argvals, fullseries, N, group, group_labels=c("S","NS"))
pval <-obtain_pvalue(res$resamples_mean, res$sample_mean, include_group=c(1,2)) #obtain pvalue 
pval$pvalue # 0.0142

# visualization 
# obtain dataframe of difference curves, sample difference curve, and time
#dat<-vis_dataframe(argvals = as.Date(argvals, origin="2011-01-01 00:00:00"), resamples_mean = res$resamples_mean, sample_mean = res$sample_mean, contrast = con)
dat<-vis_dataframe(argvals = argvals, resamples_mean = res$resamples_mean, sample_mean = res$sample_mean, contrast = con)

# convert to long data format for ggplot 
datlong<-melt(dat$dataframe, id.vars="argvals")
# plot 
main_plot2<-ggplot()+geom_line(data=datlong, aes(x=as.Date(argvals, origin="2011-01-01 00:00:00"), y=value, group=variable), colour=gray, alpha=0.2)+
  geom_line(data=subset(datlong, variable == "res.1"|variable == "sample"),aes(x=as.Date(argvals, origin="2011-01-01 00:00:00"),y=value, col=variable))+
  scale_color_manual(labels = c("Resamples", "Sample"), values=c(gray,'black'))+
  #scale_color_manual(labels = c("Resamples",bquote(paste("Sample(",tilde(H),"S)"))), values=c('black', gray))+
  theme(legend.title = element_blank())+
  labs(title="",x ="Time", y = bquote(paste("Soil Moisture Difference (",tilde(S)," - N", tilde(S), ")", sep="")), col="")+timescale+mytheme+
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
wpval<-westfall_young_pvalue(argvals=argvals, resamples_mean=res$resamples_mean, sample_mean=res$sample_mean, include_group=c(1,2), split=argvalsindex)
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


