get.summary.stats=function(breakpt,dat,max.time,ncat.data,ndata.types){
  breakpt1=c(1,breakpt,max.time)
  n=length(breakpt1)
  
  #get summarized results
  res=list()
  for (j in 1:ndata.types){
    #initialize matrix
    res[[j]]=matrix(0,n-1,ncat.data[j])
    
    #get summary for each interval
    for (i in 2:n){
      ind=breakpt1[i-1]:(breakpt1[i]-1)
      tmp=dat[ind,j]
      tmp1=table(tmp)
      ind=as.numeric(names(tmp1))
      res[[j]][i-1,ind]=tmp1
    }
  }
  res
}
#---------------------------------------------
log.marg.likel=function(alpha,summary.stats,ncat.data,ndata.types){
  #get ratios
  probs=rep(NA,ndata.types)
  for (i in 1:ndata.types){
    lnum=rowSums(lgamma(alpha+summary.stats[[i]]))  
    lden=lgamma(ncat.data[i]*alpha+rowSums(summary.stats[[i]]))
    p2=sum(lnum)-sum(lden)  
    p1=nrow(summary.stats[[i]])*(lgamma(ncat.data[i]*alpha)-ncat.data[i]*lgamma(alpha))
    probs[i]=p1+p2
  }
  sum(probs)
}
#---------------------------------------------
samp.move=function(breakpt,max.time,dat,alpha,ncat.data,ndata.types){
  breakpt.old=breakpt
  p=length(breakpt)
  new.brk=sample(1:max.time,size=1)
  brk.augmented=sort(unique(c(breakpt.old,new.brk)))
  
  p0=1  
  rand1=runif(1)
  if (p == 1) {
    #birth
    if (rand1 < 1/2){
      breakpt.new=brk.augmented
      p0=2/3 #death prob 2 -> 1 is (1/3) and birth prob 1 -> 2 is 1/2. 
    }
    #swap
    if (rand1 > 1/2) breakpt.new=new.brk
  }
  if (p > 1) {
    #birth
    if (rand1 < 1/3) {
      breakpt.new=brk.augmented
    }
    #death
    if (rand1 > 1/3 & rand1 < 2/3) {
      ind=sample(1:length(breakpt.old),size=1)
      breakpt.new=breakpt.old[-ind]
      if (p==2) p0=3/2 #birth prob from 1 -> 2 is 1/2 and death prob from 2 -> 1 is 1/3
    }
    #swap
    if (rand1 > 2/3) {
      ind=sample(1:length(breakpt.old),size=1)
      breakpt.new=sort(unique(c(breakpt.old[-ind],new.brk)))
    }
  }
  
  #get sufficient statistics
  stats.old=get.summary.stats(breakpt=breakpt.old,max.time=max.time,dat=dat,
                              ncat.data=ncat.data,ndata.types=ndata.types)
  stats.new=get.summary.stats(breakpt=breakpt.new,max.time=max.time,dat=dat,
                              ncat.data=ncat.data,ndata.types=ndata.types)
  
  #get marginal loglikel
  pold=log.marg.likel(alpha=alpha,summary.stats=stats.old,
                      ncat.data=ncat.data,ndata.types=ndata.types)
  pnew=log.marg.likel(alpha=alpha,summary.stats=stats.new,
                      ncat.data=ncat.data,ndata.types=ndata.types)+log(p0)
  prob=exp(pnew-pold)
  rand2=runif(1)
  
  if (rand2<prob) return(breakpt.new)
  return(breakpt.old)
}