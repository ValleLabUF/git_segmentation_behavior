rm(list=ls(all=TRUE))
set.seed(1)

source('gibbs functions.R')
source('gibbs sampler main function.R')
dat=read.csv('fake data.csv',as.is=T)
ndata.types=ncol(dat)

#priors
alpha=0.01
ngibbs=5000

#run gibbs sampler
breakpt=time.segm.behavior(dat=dat,ngibbs=ngibbs)

#compare estimated breakpoints to true break points
length(breakpt)
abline(v=breakpt,lty=3,col='green')  
