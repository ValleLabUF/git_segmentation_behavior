
library(tidyverse)
library(tictoc)
library(furrr)
library(viridis)
library(lubridate)


source('gibbs functions2.R')
source('helper functions.R')
source('gibbs sampler2.R')

dat<- read.csv("Modified Armadillo Data.csv", header = T, sep = ",")
dat$date<- dat$date %>% as_datetime()

#if dt within 1 min of 5 min, round to 5 min
dat<- round_track_time(dat = dat, int = 300, tol = 1/60*3600)

dat.list<- df.to.list(dat=dat)

#filter data for dt of interest
behav.list<- behav.prep(dat=dat, tstep = 300)  #add move params and filter by 300 s interval

#add discretized mvmt params
angle.bin.lims=seq(from=-pi, to=pi, by=pi/4)
max.dist=max(dat[dat$dt == 300,]$dist, na.rm = T) #using value from entire dataset, not specific time segment
upper90.thresh=as.numeric(quantile(dat[dat$dt == 300,]$dist, 0.90, na.rm=T)) #using value from entire dataset of given tstep
dist.bin.lims=seq(from=0, to=upper90.thresh, length.out = 6)
# dist.bin.lims=as.numeric(quantile(dat[dat$dt == 300,]$dist, c(0,0.25,0.50, 0.75,0.90), na.rm=T))
dist.bin.lims=c(dist.bin.lims, max.dist)

for (i in 1:length(behav.list)) {
behav.list[[i]]<- behav.list[[i]] %>% assign.dist.bin(dist.bin.lims = dist.bin.lims,
                                                 max.dist = max.dist) %>%
                                 assign.rel_angle.bin(angle.bin.lims = angle.bin.lims)
}


behav.list<- behav.list[sapply(behav.list, nrow) > 2]  #remove IDs w/ fewer than 3 obs
behav.list2<- lapply(behav.list, function(x) subset(x, select = c(id, SL, TA)))  #retain id and parameters on which to segment
dat_red<- behav.list2 %>% map(., ~mutate(., tseg = 1:nrow(.))) %>% map_dfr(., `[`) %>% dplyr::select(c(id, tseg, SL, TA))

#################################
#### Run Gibbs Sampler by ID ####
#################################

ngibbs = 40000


## Run Gibbs sampler
plan(multisession)  #run all MCMC chains in parallel
                    #refer to future::plan() for more details

dat.res<- behavior_segment(dat = behav.list2, ngibbs = ngibbs, nbins = c(4,8))
###Takes 12.5 min to run 40000 iterations for 26 IDs


## Traceplots
#type is either 'nbrks' or 'LML' for y-axis label
identity<- names(behav.list)

traceplot(data = dat.res$nbrks, type = "nbrks", identity = identity)
traceplot(data = dat.res$LML, type = "LML", identity = identity)


##Determine maximum likelihood (ML) for selecting breakpoints
ML<- apply(dat.res$LML, 1, function(x) getML(dat = x, nburn = 500))
brkpts<- getBreakpts(dat = dat.res$brkpts, ML = ML, brk.cols = 99)  #brk.cols is max matrix cols


## Heatmaps
plot.heatmap(data = behav.list, nbins = c(4,8), brkpts = brkpts, dat.res = dat.res, type = "behav")



#########################################
#### Assign Behavioral Time Segments ####
#########################################

dat_out<- map(behav.list, assign.time.seg, brkpts = brkpts) %>% map_dfr(`[`)  #assign time seg and make as DF

setwd("~/Documents/Snail Kite Project/Data/R Scripts/git_LDA_behavior")
write.csv(dat_out, "Snail Kite Gridded Data_TOHO_behav.csv", row.names = F)


