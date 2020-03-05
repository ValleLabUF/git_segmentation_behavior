
library(tidyverse)
library(tictoc)
library(furrr)
library(viridis)
library(lubridate)
library(ggforce)


source('gibbs functions.R')
source('helper functions.R')
source('gibbs sampler.R')

dat<- read.csv("Snail Kite Gridded Data_TOHO.csv", header = T, sep = ",")
dat$date<- dat$date %>% as_datetime()

#if dt within 5 min of 1 hr, round to 1 hr
dat<- round_track_time(dat = dat, int = 3600, tol = 5/60*3600)

dat.list<- df.to.list(dat=dat)

#filter data for tstep of interest
behav.list<- behav.prep(dat=dat, tstep = 3600)  #add move params and filter by 3600 s interval

#define bin number and limits for step lengths and turning angles
angle.bin.lims=seq(from=-pi, to=pi, by=pi/4)  #8 bins

max.dist=max(dat[dat$dt == 3600,]$dist, na.rm = T)
dist.bin.lims=quantile(dat[dat$dt == 3600,]$dist, c(0,0.25,0.50,0.75,0.90), na.rm=T)
dist.bin.lims=c(dist.bin.lims, max.dist)  #5 bins


#Viz limits on continuous vars
behav.df<- map_dfr(behav.list, `[`)

ggplot(behav.df, aes(x=dist/1000)) +
  geom_density(fill = "lightblue") +
  xlim(0,5) +  #limit to only 5 km
  geom_vline(xintercept = dist.bin.lims/1000, linetype = "dashed") +
  theme_bw() +
  theme(panel.grid = element_blank(), axis.title = element_text(size = 16),
        axis.text = element_text(size = 12)) +
  labs(x = "\nStep Length (km)", y = "Density\n")

ggplot(behav.df, aes(x=rel.angle)) +
  geom_density(fill = "indianred") +
  geom_vline(xintercept = angle.bin.lims, linetype = "dashed") +
  theme_bw() +
  theme(panel.grid = element_blank(), axis.title = element_text(size = 16),
        axis.text = element_text(size = 12)) +
  labs(x = "\nTurning Angle (rad)", y = "Density\n")





#assign bins to obs
for (i in 1:length(behav.list)) {
  behav.list[[i]]<- behav.list[[i]] %>%
    assign.dist.bin(dist.bin.lims = dist.bin.lims) %>%
    assign.rel_angle.bin(angle.bin.lims = angle.bin.lims)
}

behav.list<- behav.list[sapply(behav.list, nrow) > 2]  #remove IDs w/ fewer than 3 obs
behav.list2<- lapply(behav.list, function(x) subset(x, select = c(id, SL, TA)))  #retain id and parameters on which to segment



#Viz discretization of params
behav.df2<- map_dfr(behav.list2, `[`)
behav.df2<- behav.df2 %>% gather(key, value, -id)

param.prop<- behav.df2 %>%
  group_by(key, value) %>%
  summarise(n=n()) %>%
  mutate(prop=n/nrow(behav.df)) %>%
  ungroup()  #if don't ungroup after grouping, ggforce won't work

param.prop<- param.prop[-14,]
param.prop[1:5, "value"]<- ((diff(dist.bin.lims)/2) + dist.bin.lims[1:5])/1000
param.prop[6:13, "value"]<- (diff(angle.bin.lims)/2) + angle.bin.lims[1:8]


ggplot(data = param.prop %>% filter(key == "SL"), aes(value, prop)) +
  geom_bar(stat = "identity", width = (diff(dist.bin.lims)-10)/1000,
           fill = "lightblue", color = "black") +
  facet_zoom(xlim = c(0,3)) +
  labs(x = "Step Length (km)", y = "Proportion") +
  theme_bw() +
  theme(panel.grid = element_blank(), axis.title = element_text(size = 16),
        axis.text = element_text(size = 12))


ggplot(data = param.prop %>% filter(key == "TA"), aes(value, prop)) +
  geom_bar(stat = "identity", fill = "indianred", color = "black") +
  labs(x = "Turning Angle (radians)", y = "Proportion") +
  theme_bw() +
  scale_x_continuous(breaks = c(-pi, -pi/2, 0, pi/2, pi),
                     labels = expression(-pi, -pi/2, 0, pi/2, pi)) +
  theme(panel.grid = element_blank(), axis.title = element_text(size = 16),
        axis.text = element_text(size = 12))



#################################
#### Run Gibbs Sampler by ID ####
#################################

ngibbs = 40000


## Run Gibbs sampler
plan(multisession)  #run all MCMC chains in parallel
                    #refer to future::plan() for more details

#prior
alpha=1

dat.res<- behavior_segment(dat = behav.list2, ngibbs = ngibbs, nbins = c(5,8), alpha = alpha)
###Takes 25 min to run 40000 iterations for 26 IDs


## Traceplots
#type is either 'nbrks' or 'LML' for y-axis label
identity<- names(behav.list)

traceplot(data = dat.res$nbrks, type = "nbrks", identity = identity)
traceplot(data = dat.res$LML, type = "LML", identity = identity)


##Determine maximum likelihood (ML) for selecting breakpoints
ML<- apply(dat.res$LML, 1, function(x) getML(dat = x, nburn = 500))
brkpts<- getBreakpts(dat = dat.res$brkpts, ML = ML, identity = identity)


## Heatmaps
plot.heatmap(data = behav.list, nbins = c(5,8), brkpts = brkpts, dat.res = dat.res, type = "behav")



#########################################
#### Assign Behavioral Time Segments ####
#########################################

dat_out<- map(behav.list, assign.time.seg, brkpts = brkpts) %>% map_dfr(`[`)  #assign time seg and make as DF

setwd("~/Documents/Snail Kite Project/Data/R Scripts/ValleLabUF/git_LDA_behavior")
write.csv(dat_out, "Snail Kite Gridded Data_TOHO_behav2.csv", row.names = F)


