set.seed(1)

library(tidyverse)
library(tictoc)
library(furrr)
library(viridis)
library(lubridate)
library(ggforce)


source('gibbs functions.R')
source('helper functions.R')
source('gibbs sampler.R')

dat<- read.csv("Modified Armadillo Data.csv", header = T, sep = ",")
dat$date<- dat$date %>% as_datetime()

#if dt within 1 min of 5 min, round to 5 min
dat<- round_track_time(dat = dat, int = 300, tol = (1/60)*3600)

#introduce noise so that small positive step lengths added and large turning angles
tmp<- which(dat$InBurrow == 1)
dat$x[tmp]<- dat$x[tmp] + runif(length(tmp), -0.5,0.5)
dat$y[tmp]<- dat$y[tmp] + runif(length(tmp), -0.5,0.5)
dat<- prep.data(dat = dat, coord.names = c("x","y"))

# dat$LC<- as.numeric(dat$lc)  #convert string to numeric
dat$InBurrow<- dat$InBurrow + 1  #can't have 0s in data

dat.list<- df.to.list(dat=dat, ind = "id")

#filter data for dt of interest
behav.list<- behav.prep(dat=dat, tstep = 300, dat.list = dat.list)  #add move params and filter by 300 s interval

#define bin number and limits for step lengths and turning angles
angle.bin.lims=seq(from=-pi, to=pi, by=pi/4)  #8 bins

max.dist=max(dat[dat$dt == 300,]$step, na.rm = T)
dist.bin.lims=quantile(dat[dat$dt == 300,]$step, c(0,0.25,0.50,0.75,0.95), na.rm=T)
dist.bin.lims=c(dist.bin.lims, max.dist)  #5 bins


#Viz limits on continuous vars
behav.df<- map_dfr(behav.list, `[`)

ggplot(behav.df, aes(x=step/1000)) +
  geom_density(fill = "lightblue") +
  xlim(0,1) +  #limit to only 1 km
  geom_vline(xintercept = dist.bin.lims/1000, linetype = "dashed") +
  theme_bw() +
  theme(panel.grid = element_blank(), axis.title = element_text(size = 16),
        axis.text = element_text(size = 12)) +
  labs(x = "\nStep Length (km)", y = "Density\n")

ggplot(behav.df, aes(x=angle)) +
  geom_density(fill = "indianred") +
  geom_vline(xintercept = angle.bin.lims, linetype = "dashed") +
  theme_bw() +
  theme(panel.grid = element_blank(), axis.title = element_text(size = 16),
        axis.text = element_text(size = 12)) +
  labs(x = "\nTurning Angle (rad)", y = "Density\n")



#assign bins to obs
behav.list<- map(behav.list, discrete_move_par, lims = list(dist.bin.lims, angle.bin.lims),
                 varIn = c("step", "angle"), varOut = c("SL", "TA"))

#replace TA NAs when in burrow to bin 9
# behav.list<- map(behav.list,
#                  ~mutate_at(., "TA",
#                             function(x) ifelse(is.na(.$TA) & .$InBurrow == 1, 9, .$TA)))
behav.list<- behav.list[map_dbl(behav.list, nrow) > 2]  #remove IDs w/ fewer than 3 obs
behav.list2<- lapply(behav.list, function(x) subset(x, select = c(id, SL, TA, InBurrow)))  #retain id and parameters on which to segment



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
  # facet_zoom(xlim = c(0,0.1)) +
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

#prior
alpha=1

#subset data and create list of proposed breakpts by ID
# test<- behav.list2[c(1,2,6)]  #IDs tm13, tm14, tm24
breaks<- map(behav.list, find.breaks, "InBurrow")

## Run Gibbs sampler
plan(multisession)  #run all MCMC chains in parallel
                    #refer to future::plan() for more details

dat.res<- behavior_segment(data = behav.list2, ngibbs = ngibbs, nbins = c(5,8,2), alpha = alpha,
                           breakpt = breaks)
###Takes 82 min to run 40000 iterations for all IDs


## Traceplots
#type is either 'nbrks' or 'LML' for y-axis label
traceplot(data = dat.res$nbrks, type = "nbrks")
traceplot(data = dat.res$LML, type = "LML")


##Determine maximum likelihood (ML) for selecting breakpoints
ML<- getML(dat = dat.res$LML, nburn = 20000)
brkpts<- getBreakpts(dat = dat.res$brkpts, ML = ML)


## Heatmaps
plot.heatmap(data = behav.list2, nbins = c(5,8,2), brkpts = brkpts, dat.res = dat.res,
             type = "behav", title = TRUE, legend = TRUE)



#########################################
#### Assign Behavioral Time Segments ####
#########################################

dat_out<- map(behav.list, assign.time.seg, brkpts = brkpts) %>% map_dfr(`[`)  #assign time seg and make as DF

setwd("~/Documents/Snail Kite Project/Data/R Scripts/ValleLabUF/git_LDA_behavior")
# write.csv(dat_out, "Armadillo Data_behav.csv", row.names = F)



# #Remove all but first Burrow obs per period
# #add time indicator for later merging of behavior
# dat$time1<- 1:nrow(dat)
# 
# #filter for focal ID
# dat.list<- df.to.list(dat, "id")
# 
# #identify first burrow obs per period
# breaks<- map(dat.list, find.breaks, "InBurrow")
# breaks1<- map(breaks, function(x) x[seq(1, length(x), by=2)])
# 
# #remove all burrow locs besides first obs
# #also make sure burrow obs retained have dt == 300
# dat.list2<- map2(dat.list, breaks1, ~{.x %>% 
#     slice(sort(c(which(.$InBurrow == 0), .y))) %>% 
#     mutate_at("dt", function(x) ifelse(.$InBurrow == 1, 300, .$dt))}) %>% 
#   map(., ~round_track_time(., int = 300, tol = (1/60)*3600))
# 
# dat2<- bind_rows(dat.list2)




# #Extract LU/LC from shp to use for segmentation model
# st_crs(classes_DL_padrao.shp)<- "+init=epsg:32721"
# st_crs(classes_STseca.shp)<- "+init=epsg:32721"
# 
# #Using 'sf'
# dat.sf<- st_as_sf(dat, coords = c("x","y"), crs = "+init=epsg:32721")
# 
# 
# foo<- st_join(dat.sf %>% filter(region == "N"), classes_DL_padrao.shp["SPRCLASSE"],
#               join = st_within)
# 
# dat.sf.N<- dat.sf %>% filter(region == "N")
# dat.sf.N$lc<- unlist(over(as(dat.sf.N, "Spatial"),
#                           as(classes_DL_padrao.shp["SPRCLASSE"], "Spatial")))
# 
# 
# dat.sf.S<- dat.sf %>% filter(region == "S")
# dat.sf.S$lc<- unlist(over(as(dat.sf.S, "Spatial"),
#                           as(classes_STseca.shp["SPRCLASSE"], "Spatial")))
# 
# 
# dat.sf2<- rbind(dat.sf.N, dat.sf.S) %>% 
#   st_drop_geometry()
# #---------------------
# 
# #Using 'sp'
# dat.spdf<- SpatialPointsDataFrame(coords = dat[,c("x","y")], data = dat,
#                                   proj4string = CRS("+init=epsg:32721"))
# 
# dat.spdf.N<- subset(dat.spdf, dat.spdf$region == "N")
# dat.spdf.S<- subset(dat.spdf, dat.spdf$region == "S")
# 
# dat.spdf.N$lc<- unlist(over(as(dat.sf.N, "Spatial"),
#                             as(classes_DL_padrao.shp["SPRCLASSE"], "Spatial")))
# dat.spdf.S$lc<- unlist(over(as(dat.sf.S, "Spatial"),
#                             as(classes_STseca.shp["SPRCLASSE"], "Spatial")))
# dat<- rbind(dat.spdf.N@data, dat.spdf.S@data)
