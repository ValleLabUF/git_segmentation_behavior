
library(tidyverse)
library(tictoc)
library(furrr)
library(viridis)
library(lubridate)


source('gibbs functions.R')
source('helper functions.R')
source('gibbs sampler.R')

dat<- read.csv("Modified Snow Leopard Data.csv", header = T, sep = ",")
dat$date<- dat$date %>% as_datetime()

#if dt within 5 min of 3 h, round to 3 h
dat<- round_track_time(dat = dat, int = 10800, tol = 300)

dat.list<- df.to.list(dat=dat)

#filter data for dt of interest
behav.list<- behav.prep(dat=dat, tstep = 10800)  #add move params and filter by 10800 s interval

#define bin number and limits for step lengths and turning angles
angle.bin.lims=seq(from=-pi, to=pi, by=pi/4)  #8 bins

max.dist=max(dat[dat$dt == 10800,]$dist, na.rm = T)
dist.bin.lims=quantile(dat[dat$dt == 10800,]$dist, c(0,0.25,0.50,0.75,0.90), na.rm=T)
dist.bin.lims=c(dist.bin.lims, max.dist)  #5 bins

for (i in 1:length(behav.list)) {
behav.list[[i]]<- behav.list[[i]] %>% assign.dist.bin(dist.bin.lims = dist.bin.lims) %>%
                                 assign.rel_angle.bin(angle.bin.lims = angle.bin.lims)
}


behav.list<- behav.list[sapply(behav.list, nrow) > 2]  #remove IDs w/ fewer than 3 obs
behav.list2<- lapply(behav.list, function(x) subset(x, select = c(id, SL, TA)))  #retain id and parameters on which to segment


#################################
#### Run Gibbs Sampler by ID ####
#################################

ngibbs = 40000

#prior
alpha=1

## Run Gibbs sampler
plan(multisession)  #run all MCMC chains in parallel
                    #refer to future::plan() for more details

dat.res<- behavior_segment(dat = behav.list2, ngibbs = ngibbs, nbins = c(5,8), alpha = alpha)
###Takes 8 min to run 40000 iterations for 5 IDs


## Traceplots
#type is either 'nbrks' or 'LML' for y-axis label
identity<- names(behav.list)

traceplot(data = dat.res$nbrks, type = "nbrks", identity = identity)
traceplot(data = dat.res$LML, type = "LML", identity = identity)


##Determine maximum likelihood (ML) for selecting breakpoints
ML<- apply(dat.res$LML, 1, function(x) getML(dat = x, nburn = 500))
brkpts<- getBreakpts(dat = dat.res$brkpts, ML = ML, identity = identity)


## Heatmaps
plot.heatmap(data = behav.list, nbins = c(5,8), brkpts = brkpts, dat.res = dat.res,
             type = "behav")



#########################################
#### Assign Behavioral Time Segments ####
#########################################

dat_out<- map(behav.list, assign.time.seg, brkpts = brkpts) %>%
  map_dfr(`[`)  #assign time seg and make as DF

setwd("~/Documents/Snail Kite Project/Data/R Scripts/ValleLabUF/git_LDA_behavior")
write.csv(dat_out, "Snow Leopard Data_behav.csv", row.names = F)


