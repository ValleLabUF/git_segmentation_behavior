assign.rel_angle.bin=function(dat, angle.bin.lims){
  dat$TA<- NA
  
  for(i in 1:length(angle.bin.lims)) {
    tmp=which(dat$rel.angle >= angle.bin.lims[i] & dat$rel.angle < angle.bin.lims[i+1])
    dat[tmp,"TA"]=i
  }
  dat
}
#---------------------------------------------
assign.dist.bin=function(dat, dist.bin.lims){
  
  dat$SL<- NA
  
  for(i in 1:length(dist.bin.lims)) {
    tmp=which(dat$dist >= dist.bin.lims[i] & dat$dist < dist.bin.lims[i+1])
    dat[tmp,"SL"]=i
  }
  tmp=which(dat$dist == dist.bin.lims[length(dist.bin.lims)])
  dat[tmp,"SL"]=length(dist.bin.lims) - 1
  
  dat
}
#---------------------------------------------
round_track_time=function(dat, int, tol) {  #replacement for sett0() when wanting to only round some of the times
  dat<- df.to.list(dat)
  for (i in 1:length(dat)) {
    tmp=matrix(NA,nrow(dat[[i]]),2)
    if (length(int) == 1) {  #when using only 1 time interval
      for (j in 1:nrow(dat[[i]])) {
        if (is.na(dat[[i]]$dt[j])) {
          tmp[j, 1:2]<- NA
        } else if (dat[[i]]$dt[j] >= (int - tol) & dat[[i]]$dt[j] <= (int + tol)) {
          tmp[j, 1:2]<- c(int, as.numeric(round(dat[[i]]$date[j], units = "hours")))
        } else {
          tmp[j, 1:2]<- c(dat[[i]]$dt[j], dat[[i]]$date[j])
        }
      }
    } else {  #when using more than one time interval
      for (j in 1:nrow(dat[[i]])) {
        if (is.na(dat[[i]]$dt[j])) {
          tmp[j, 1:2]<- NA
        } else if (sum(dat[[i]]$dt[j] >= (int - tol) & dat[[i]]$dt[j] <= (int + tol)) > 0) {
          ind<- which(dat[[i]]$dt[j] >= (int - tol) & dat[[i]]$dt[j] <= (int + tol))
          tmp[j, 1:2]<- c(int[ind], as.numeric(round(dat[[i]]$date[j], units = "hours")))
        } else {
          tmp[j, 1:2]<- c(dat[[i]]$dt[j], dat[[i]]$date[j])
        }
      }
    }
    dat[[i]]$dt<- tmp[,1]
    dat[[i]]$date<- tmp[,2] %>% as.POSIXct(origin = '1970-01-01', tz = "UTC")
  }
  dat<- map_dfr(dat, `[`)
  dat
}
#---------------------------------------------
#discretize step length and turning angle by ID
#only subset obs w/ 1 hr (3600 s) time step; results in loss of irregular obs
behav.prep=function(dat,tstep) {
  id<- unique(dat$id)
  
  cond<- matrix(0, length(dat.list), 1)
  for (i in 1:length(dat.list)) {  #don't include IDs w < 1 obs of dt == tstep
    cond[i,]<- if(nrow(dat.list[[i]][dat.list[[i]]$dt == tstep,]) > 1) {
      1
    } else {
      0
    }
  }
  tmp<- cbind(id,cond)
  ind<- which(tmp[,2] == 1)
  
  n=length(ind)
  behav.list<- vector("list", n)
  names(behav.list)<- id[ind]
  
  for (i in 1:n) {
    behav.list[[i]]<- dat[dat$id==id[ind[i]],] %>% mutate(obs = 1:nrow(.)) %>% filter(dt==tstep) %>%
                  mutate(time1 = 1:nrow(.)) 
  }
  behav.list
}
#---------------------------------------------
behav.seg.image=function(dat, nbins) {  #Transform single var vectors into pres/abs matrices for heatmap; nbins is vector of bins per param in order
  SL<- matrix(0, nrow(dat), nbins[1])
  for (i in 1:nrow(dat)){
    SL[i,dat$SL[i]]=1
  }
  
  TA<- matrix(0, nrow(dat), nbins[2])
  for (i in 1:nrow(dat)){
    TA[i,dat$TA[i]]=1
  }
  
  
  behav.list<- list(SL=SL,TA=TA)
  behav.list
}
#---------------------------------------
assign.time.seg=function(dat, brkpts){  #add tseg assignment to each obs
  tmp=which(unique(dat$id) == brkpts$id)
  breakpt<- brkpts[tmp,-1] %>% purrr::discard(is.na) %>% as.numeric(.[1,])
  breakpt1=c(0,breakpt,Inf)
  tmp1<- which(diff(breakpt1) < 1)  #check for impossible time units
  breakpt1[tmp1+1]<- breakpt1[(tmp1+1)] + 1  #fix impossible time units
  n=length(breakpt1)
  res=matrix(NA,nrow(dat),1)
  for (i in 2:n){
    ind=which(breakpt1[i-1]<=dat$time1 & dat$time1<breakpt1[i])
    res[ind,]=i-1
  }
  dat$tseg<- as.vector(res)
  dat
}
#------------------------------------------------
df.to.list=function(dat) {  #only for id as col in dat
  id<- unique(dat$id)
  n=length(id)
  dat.list<- vector("list", n)
  names(dat.list)<- id
  
  for (i in 1:length(id)) {
    dat.list[[i]]<- dat[dat$id==id[i],]
  }
  dat.list
}
#------------------------------------------------
traceplot=function(data, type, identity) {  #create traceplots for nbrks or LML for all IDs
  for (i in 1:length(identity)) {
    par(ask=TRUE)
    plot(x=1:ngibbs, y=data[i,-1], type = "l", xlab = "Iteration",
         ylab = ifelse(type == "nbrks", "# of Breakpoints",
                       ifelse(type == "LML","Log Marginal Likelihood",
                              stop("Need to select one of 'nbrks' or 'LML' for plotting"))),
         main = paste("ID",rownames(data)[i]))
  }
  on.exit(par(ask = FALSE))
}
#---------------------------------------------
getML=function(dat,nburn) {  #select ML value that is beyond burn-in phase
  if (which.max(dat[-1]) < nburn) {
    ML<- dat[-1] %>% order(decreasing = T) %>% subset(. > nburn) %>% first()
  } else {
    ML<- which.max(dat[-1])
  }
  return(ML)
}
#---------------------------------------------
getBreakpts=function(dat,ML,identity) {  #extract breakpoints of ML per ID
  tmp<- list()
  
  for(i in 1:length(dat)) {
    ind<- ML[i]
    tmp[[i]]<- dat[[i]][[ind]]
  }
  
  names(tmp)<- identity
  max.length<- max(sapply(tmp, length))
  tmp<- lapply(tmp, function(x) { c(x, rep(NA, max.length-length(x)))})
  tmp<- map_dfr(tmp, `[`) %>% t() %>% data.frame()
  tmp<- cbind(id = identity, tmp)
  names(tmp)<- c('id', paste0("Brk_",1:(ncol(tmp)-1)))
  
  tmp
}
#------------------------------------------------
plot.heatmap.behav=function(data, nbins, brkpts, dat.res) {
  
  behav.heat<- behav.seg.image(data, nbins)
  
  SL<- data.frame(behav.heat$SL)
  names(SL)<- 1:nbins[1]
  SL<- SL %>% gather(key, value) %>% mutate(time=rep(data$time1, times=nbins[1]),
                                            behav=rep("SL", nrow(data)*nbins[1]))
  
  TA<- data.frame(behav.heat$TA)
  names(TA)<- 1:nbins[2]
  TA<- TA %>% gather(key, value) %>% mutate(time=rep(data$time1, times=nbins[2]),
                                            behav=rep("TA", nrow(data)*nbins[2]))
  
  behav.heat_long<- rbind(SL,TA)
  behav.heat_long$value<- factor(behav.heat_long$value)
  levels(behav.heat_long$value)<- c("Unoccupied","Occupied")
  
  ind=which(unique(data$id) == brkpts$id)
  breakpt<- brkpts[ind,-1] %>% purrr::discard(is.na) %>% t() %>% data.frame()
  names(breakpt)<- "breaks"
  
  print(
  ggplot(behav.heat_long, aes(x=time, y=key, fill=value)) +
    geom_tile() +
    facet_wrap(~behav, scales = 'free', nrow = 2) +
    scale_fill_viridis_d('') +
    scale_y_discrete(expand = c(0,0)) +
    scale_x_continuous(expand = c(0,0)) +
    geom_vline(data = breakpt, aes(xintercept = breaks), color = viridis(n=9)[7],
               size = 0.4) +
    labs(x = "Observations", y = "Bin", title = paste("ID", unique(data$id))) +
    theme_bw() +
    theme(axis.title = element_text(size = 18), axis.text = element_text(size = 12),
          strip.text = element_text(size = 12, face = 'bold'),
          title = element_text(size = 20))
  )
}
#------------------------------------------------
plot.heatmap=function(data, nbins, brkpts, dat.res, type) {  #type can either be 'loc' or 'behav'
  
  if (type == "loc") {
    par(ask = TRUE)
    map(data, ~plot.heatmap.loc(., nbins = nbins, brkpts = brkpts, dat.res = dat.res))
    par(ask = FALSE)
  } else if (type == "behav") {
    par(ask = TRUE)
    map(data, ~plot.heatmap.behav(., nbins = nbins, brkpts = brkpts, dat.res = dat.res))
    par(ask = FALSE)
  } else {
    stop("Need to select type as either 'loc' or 'behav'")
  }
}
