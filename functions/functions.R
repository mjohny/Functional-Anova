#function to seperate and order full series by group
group_data<-function(fulldata, group, group_names){
  grouped_data<-NULL
  for (i in unique(group)){
    grouped_data[[i]]<-fulldata[which(group==i),]
  }
  result<-grouped_data
  if (missing(group_names)){
  names(result)<-paste("Group",sort(unique(group)), sep="")}
  else {names(result)<-group_names}
  return(result)
}

# function to output resample means
obtain_resamples<-function(argvals, fulldata, N, group, group_names){
  # seperate full series by group
  if (missing(group_names)){
    series<-group_data(fulldata, group)
  }
  else{
    series<-group_data(fulldata, group, group_names)
  }
  
  sample_mean <-matrix(NA, nrow = length(argvals), ncol=length(unique(group)))
  
  # obtain sample means 
  for (j in 1:length(series)){
    smooth<-fdata(series[[j]], argvals=argvals, rangeval<-rangeval)
    smooth_mean<-func.mean(smooth)
    sample_mean[,j]<-smooth_mean$data
  }
  dimnames(sample_mean)<-list(c(seq(1:length(argvals))),names(series))
  # define empty vectors
  K<-NULL
  resamples<-NULL
  resamples_mean<-NULL
  res_mean<-array(NA, dim=c(N,length(argvals),length(unique(group))))
  
  # simulate resample curves
  for (i in 1:N){
    for (j in 1:length(series)){
      K[[j]]<-mvrnorm(n=nrow(series[[j]]), mu=rep(0, length(argvals)), var(series[[j]]))
      resamples[[j]]<-fdata(K[[j]],argvals=argvals, rangeval<-rangeval)
      resamples_mean[[j]]<-func.mean(resamples[[j]])
      res_mean[i,,j]<-resamples_mean[[j]]$data
    }
  }
  dimnames(res_mean)<-list(c(1:N),c(seq(1:length(argvals))),names(series))
  result <- list(group_names=names(series),resamples_mean=res_mean, sample_mean=sample_mean)
  return(result)
}

obtain_pvalue<-function(resamples_mean, sample_mean, group){
  N<-dim(resamples_mean)[1]
  d<-rep(NA, N)
  index<-which(group>0)
  for (i in 1:N){
    curves<-t(resamples_mean[i,,index])
    dists<-metric.lp(curves)
    d[i]<-sum(dists[lower.tri(dists)]^2)
  }
  
  #obtain sample distance
  sample_dists<-metric.lp(t(sample_mean[,index]))
  stat<-sum(sample_dists[lower.tri(sample_dists)]^2)
  
  #pvalue
  pvalue<-length(d[d>=stat])/N

  # group names
  group_names<-dimnames(sample_mean)[[2]][index]
  
  result<-list(group_names=group_names,pvalue=pvalue)

  return(result)
}

obtain_pvalue_contrast<-function(resamples_mean, sample_mean, contrast){
  #obtain resample distances
  N<-dim(resamples_mean)[1]
  con<-contrast
  d<-rep(NA, N)
  res_con<-0
  con_names<-NULL
  for (i in 1:dim(resamples_mean)[3]){
    res_con<-res_con+resamples_mean[,,i]*con[i]
    con_names<-paste(con_names,con[i],"*",dimnames(resamples_mean)[[3]][i],",", sep="")
  }
  d<-sqrt(rowSums(res_con^2)) # bootstrap distances 
  
  #obtain sample distance
  sample_con<-0
  for (i in 1:dim(sample_mean)[2]){
    sample_con<-sample_con+sample_mean[,i]*con[i]
  }
  stat<-sqrt(sum(sample_con^2))
  
  #pvalue
  pvalue<-length(d[d>=stat])/N
  
  result<-list(contrast=con_names, pvalue=pvalue)
  return(result)
}

# function to obtain tukey p-value 
tukey_correction<-function(resamples_mean, sample_mean, group){
  #obtain resample distances
  N<-dim(resamples_mean)[1]
  index<-which(group>0)
  #sample_mean_group<-array(NA, dim=c(1,length(argvals),2))
  
  #obtain sample distance
  sample_dists<-metric.lp(t(sample_mean[,index]))
  stat<-sum(sample_dists[lower.tri(sample_dists)]^2)
  
  #obtain tukey distances
  maxindex<-rep(NA, N)
  maxdist<-rep(NA, N)
  # obtain every pair of groups
  ngroups<-dim(resamples_mean)[3]
  g<-c(1:ngroups)
  pairs<-combn(g,2)
  distvals<-rep(NA,dim(pairs)[2])
  curves_pair<-array(NA, dim=c(2,length(argvals),dim(pairs)[2]))
  maxresamples<-array(NA, dim=c(N,length(argvals),2))
  # obtain distance between each pair 
  for (i in 1:N){
    for (j in 1:dim(pairs)[2]){
      curves_pair[,,j]<-t(resamples_mean[i,,pairs[,j]])
      dists<-metric.lp(curves_pair[,,j])
      # save distance 
      distvals[j]<-sum(dists[lower.tri(dists)]^2)
    }
    # select index of pair corresponding to largest distance
    maxindex[i]<-which.max(distvals)
    # save max distance 
    maxdist[i]<-max(distvals)
    # max resamples pair
    maxresamples[i,,]<-t(curves_pair[,,j])
  }
  
  # tukey-pvalue
  pvalue<-length(maxdist[maxdist>=stat])/N
  
  # tukey resample curves
  resamples_tukey<-maxresamples
  dimnames(resamples_tukey)<-list(dimnames(resamples_mean)[[1]],dimnames(resamples_mean)[[2]],c("maxgroup1","maxgroup2"))
  # sample curves 
  group_names<-dimnames(sample_mean)[[2]][index]
  sample_mean<-sample_mean[,index]

   
  result<-list(group_names=group_names,pvalue=pvalue, resamples_tukey=resamples_tukey, sample_mean=sample_mean)
  
  return(result)
  
}

# function to obtain difference curves 
obtain_contrast<-function(resamples_mean, sample_mean, contrast){
  con <- contrast
  con2<- as.character(contrast)
  diff<-0
  con_names<-NULL
  for (i in 1:dim(resamples_mean)[3]){
    if (sign(con[i])==-1 & i>1){con_names<-substr(con_names, 1, nchar(con_names)-1)}
    else {con_names<-con_names}
    diff<-diff+resamples_mean[,,i]*con[i]
    if (con2[i]=="-1"){con2[i]<-"-"}
    if (con2[i]=="1"){con2[i]<-""}
    else {con2[i]<-con2[i]}
    con_names<-paste(con_names,con2[i],dimnames(resamples_mean)[[3]][i],"+", sep="")
  }
  con_names<-substr(con_names, 1, nchar(con_names)-1)
  resamples_diff<-diff
  
  sample_diff<-0
  for (i in 1:dim(sample_mean)[2]){
    sample_diff<-sample_diff+sample_mean[,i]*con[i]
  }
  sample_diff<-sample_diff
  
  result<-list(contrast=con_names, resamples_diff=resamples_diff, sample_diff=sample_diff)
  return(result)
}

# create vis dataframe
vis_dataframe <- function(argvals,resamples_mean,sample_mean,contrast){
  # obtain resample difference curves
  res_diff<-obtain_contrast(resamples_mean,sample_mean,contrast)
  res_curves<-res_diff$resamples_diff
  sample_curves<-res_diff$sample_diff
  
  # create data frame
  df<-data.frame(argvals=c(argvals),res=t(res_curves),sample=c(sample_curves))
  contrasts<-res_diff$contrast
  
  result<-list(dataframe=df, contrasts=contrasts)
  return(result)}

# westfall young correction 
westfall_young_pvalue<-function(argvals, resamples_mean, sample_mean, group, split){
  N<-dim(resamples_mean)[1] #number of bootstrap resamples 
  group_index<-unique(group[which(group>0)]) #group indexs involved in test
  p_raw<-rep(NA, length(split)) #empty vector to store raw pvalue for each split 
  split<-split
  # list of vectors to store distances for each time split 
  resdistall<-list(NULL) # list of empty vectors for storing bootstrap distances
  for (j in 1:length(argvalsindex)){
    resdistall[[j]]<-rep(NA, N)
  }
  
  # obtain bootstrap distance for each split
  for (i in 1:N){
    for (j in 1:length(split)){
      split_index<-split[[j]]
      curves<-t(resamples_mean[i,split_index,group_index])
      dists<-metric.lp(curves)
      resdistall[[j]][i]<-sum(dists[lower.tri(dists)]^2)
    }
  }
  
  # obtain sample dist for each split 
  for (j in 1:length(split)){
    split_index<-split[[j]]
    curves<-t(sample_mean[split_index,group_index])
    dists<-metric.lp(curves)
    stat<-sum(dists[lower.tri(dists)]^2)
    d<-resdistall[[j]]
    p_raw[j]<-length(d[d>=stat])/N
  }
  
  # obtain corrected pvalues 
  p_star<-list(NULL)
  for (i in 1:length(argvalsindex)){
    p_star[[i]]<-rep(NA, N)
  }
  
  for (i in 1:N){
  for (j in 1:length(argvalsindex)){
    stat<-resdistall[[j]][i]
    d<-resdistall[[j]][-i]
    p_star[[j]][i]<-length(d[d>=stat])/(N-1)
  }
  }
  
  #order the raw p-value in increasing order 
  p_index<-order(p_raw) 
  p_raw2<-p_raw[p_index] #ordered p-values 
  
  #arrange p_stars in same order as p_raw and compute successive minima 
  p_star2<-list(NULL)
  
  for (i in 1:length(p_index)){
    p_star2[[i]]<-p_star[[p_index[i]]]}
  
  # obtain corrected pvalues 
  qr<-matrix(NA, nrow=N, ncol=length(argvalsindex))
  
  #combine p_stars into matrix where col = domain, and rows = N ittertions
  p_star2_comb<-NULL
  p_star2_comb<-cbind(c(p_star2[[1]]))
  for (i in 2:length(p_star2)){
    p_star2_comb<-cbind(p_star2_comb, c(p_star2[[i]]))
  }
  
  # compute succesive minimum over the rows 
  for (i in 1:nrow(qr)){
    qr[i,1]<-min(p_star2_comb[i,])
    for (j in 2:ncol(qr)){
      qr[i,j]<-min(p_star2_comb[i,-c(1:j-1)])
    }
  }
  
  # corrected p-value 
  p_cor<-rep(NA,length(argvalsindex))
  for (i in 1:length(p_cor)){
    p_cor[i]<-length(qr[,i][qr[,i]<=p_raw2[i]])/N
  }
  
  # enforce monotonicity
  p_cor2<-p_cor
  for (i in 1:19){
    if(p_cor2[i+1]>=p_cor2[i]){
      p_cor2[i+1]<-p_cor2[i+1]}
    else {p_cor2[i+1]<-p_cor2[i]}
  }
  
  #reorder corrected pvalue to original order 
  p_cor3<-rep(NA,length(argvalsindex))
  for (i in 1:length(p_cor3)){
    p_cor3[i]<-p_cor2[order(p_index)[i]]
  }
  
  result<-list(raw_pvalue=p_raw, p_star=p_star, corr_pvalue=p_cor3)
  return(result)
}


