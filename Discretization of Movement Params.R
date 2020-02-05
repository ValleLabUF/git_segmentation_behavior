# Generate 8 bins for relative turning angle (TA)

assign.rel_angle.bin=function(dat){
 angle.bin.lims=seq(from=-pi, to=pi, by=pi/4)
 dat$TA<- NA

  for(i in 1:length(angle.bin.lims)) {
   tmp=which(dat$rel.angle >= angle.bin.lims[i] & dat$rel.angle < angle.bin.lims[i+1])
   dat[tmp,"TA"]=i
  }
  dat
}

dat<- assign.rel_angle.bin(dat)

#View histogram
ggplot(dat[!is.na(dat$TA),], aes(factor(TA))) +
  geom_bar() +
  labs(x = "Bin #")






#Generate 6 bins for step length (SL); use 90th percentile as last cutoff

assign.dist.bin=function(dat){
 max.dist=max(dat$dist, na.rm = T) #using value from entire dataset, not specific time segment
 upper90.thresh=as.numeric(quantile(dat$dist, 0.90, na.rm=T)) #using value from entire dataset
 
 dist.bin.lims=seq(from=0, to=upper90.thresh, length.out = 6)
 dist.bin.lims=c(dist.bin.lims, max.dist)
 dat$SL<- NA
 #names(dat)[7]<- "dist"

  for(i in 1:length(dist.bin.lims)) {
    tmp=which(dat$dist >= dist.bin.lims[i] & dat$dist < dist.bin.lims[i+1])
    dat[tmp,"SL"]=i
  }
 tmp=which(dat$dist == max.dist)
 dat[tmp,"SL"]=6
 
  dat
}

dat<- assign.dist.bin(dat)

#View histogram
ggplot(dat[!is.na(dat$SL),], aes(factor(SL))) +
  geom_bar() +
  labs(x = "Bin #")





#Assignment of binary response variable for changes in turning angle autocorrelation (TAA)
  #Try both with the sign of the direction

#sign
chng.rel_angle.sign=function(dat){
  dat$TAA<- NA
  
  for(i in 1:(nrow(dat)-1)) {
    dat$TAA[i+1]<- ifelse(sign(dat$rel.angle[i]) == sign(dat$rel.angle[i+1]), 1, 0)
  }
  dat
}

dat<- chng.rel_angle.sign(dat)



#view time series
ggplot(dat[!is.na(dat$TAA),], aes(x=1:nrow(dat[!is.na(dat$TAA),]),
                                           y=factor(TAA))) +
  geom_point(size=2) +
  scale_y_discrete("Change in Sign of Relative Turning Angle", breaks=c(0,1),
                   labels=c("Change Direction","Maintain Direction")) +
  xlab("Time")


