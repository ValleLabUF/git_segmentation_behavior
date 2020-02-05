rm(list=ls(all=TRUE))
library('MCMCpack')
set.seed(4)

nobs=2000
nseg=20
tmp=runif(nseg)
prob=tmp/sum(tmp); prob
partition=rmultinom(1,size=nobs,prob=prob)
seg.index=rep(1:nseg,times=partition)
nloc=100

#generate probabilities for each type of data
ncat.data=c(8,3,10,2,4,6) #number of bins for each data type
ndata.type=length(ncat.data)
probs=list()
for (i in 1:ndata.type){
  probs[[i]]=rdirichlet(nseg,alpha=rep(0.1,ncat.data[i]))
}

#generate observations
obs=matrix(NA,nobs,ndata.type)
for (i in 1:nobs){
  for (j in 1:ndata.type){
    tmp=rmultinom(1,size=1,prob=probs[[j]][seg.index[i],])
    obs[i,j]=which(tmp==1)
  }
}

#name these variables
nome=paste0('y',1:ndata.type)
colnames(obs)=nome

#plot these data
par(mfrow=c(ndata.type,1),mar=rep(1,4))
for (i in 1:ndata.type){
  plot(obs[,i])
  abline(v=cumsum(partition))
}

#add time
# time1=1:nobs
# obs1=cbind(obs,time1)

#export results
write.csv(obs,'fake data.csv',row.names=F)
