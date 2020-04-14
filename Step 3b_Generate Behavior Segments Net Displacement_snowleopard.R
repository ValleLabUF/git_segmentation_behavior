library(tidyverse)
library(tictoc)
library(furrr)
library(viridis)
library(lubridate)
library(ggforce)
library(gridExtra)


source('gibbs functions.R')
source('helper functions.R')
source('gibbs sampler.R')

dat<- read.csv("Modified Snow Leopard Data.csv", header = T, sep = ",")
dat$date<- dat$date %>% as_datetime()
dat$displ<- dat$R2n %>% sqrt()/1000  #create Net Displacement var and change from m to km

#Inspect time series of and displacement
ggplot(dat, aes(x=date, y=displ)) +
  geom_path() +
  facet_wrap(~id, scales = "free") + #separate IDs
  theme_bw() +
  labs(x = "Time", y = expression("Net Displacement "(km)))

#Inspect histograms of displacement
ggplot(dat, aes(x=displ)) +
  geom_histogram(binwidth = 0.5, fill = "black") +
  facet_wrap(~id, scales = "fixed") +
  theme_bw() +
  labs(x = expression("Net Displacement "(km)), y = "Frequency")



dat.list<- df.to.list(dat=dat, ind="id") %>%
            map(~mutate(.x, time1 = 1:nrow(.x)))  #add observation col


##Discretize net displacement
displ.bin.lims<- seq(0, ceiling(max(dat$displ, na.rm = TRUE)), by = 0.5)  #create bin limits 
dat.list<- map(dat.list, discrete_move_par, lims = list(displ.bin.lims), varIn = "displ",
               varOut = "ND")
dat.list2<- map(dat.list, subset, select = c(id, ND))


#############################
#############################
#### Normal 2-pass model ####
#############################
#############################


#################################
#### Run Gibbs Sampler by ID ####
#################################

ngibbs = 5000

nbins = length(displ.bin.lims)

#prior
alpha=1

## Run Gibbs sampler
plan(multisession)  #run all MCMC chains in parallel
#refer to future::plan() for more details

dat.res<- behavior_segment(dat = dat.list2, ngibbs = ngibbs, nbins = nbins, alpha = alpha)



## Traceplots
#type is either 'nbrks' or 'LML' for y-axis label
identity<- names(dat.list)

traceplot(data = dat.res$nbrks, type = "nbrks", identity = identity)
traceplot(data = dat.res$LML, type = "LML", identity = identity)


##Determine maximum likelihood (ML) for selecting breakpoints
ML<- apply(dat.res$LML, 1, function(x) getML(dat = x, nburn = 500))
brkpts<- getBreakpts(dat = dat.res$brkpts, ML = ML, identity = identity)


brkpts_t<- brkpts %>%  #change from wide to long format
  pivot_longer(-id, values_to = "value", values_drop_na = TRUE)

ggplot(map_dfr(dat.list, `[`), aes(x=time1, y=displ)) +
  geom_path() +
  geom_vline(data = brkpts_t, aes(xintercept = value, group = id), color = "blue") +
  facet_wrap(~id, scales = "free_x") + #separate IDs
  theme_bw() +
  labs(x = "Time", y = expression("Net Displacement "(km))) +
  theme(panel.grid = element_blank(),
        strip.text = element_text(face = "bold", size = 10))



##2nd pass of model

#normalize displ_n var to min value of tseg to create new baseline before collating all data again by ID
dat.list3<- map(dat.list, assign.time.seg, brkpts = brkpts) %>% 
              map(~df.to.list(.x, ind = "tseg")) %>%
              modify_depth(2, ~mutate(.x, displ_n = displ - min(displ, na.rm = T))) %>%
              modify_depth(1, ~map_dfr(.x, `[`))

#Compare differences in 'displ' on original and normalized scales to inspect differences
p1<- ggplot(map_dfr(dat.list3, `[`), aes(x=time1, y=displ)) +
  geom_path() +
  geom_vline(data = brkpts_t, aes(xintercept = value, group = id), color = "blue") +
  facet_wrap(~id, scales = "free_x", ncol = 1) + #separate IDs
  theme_bw() +
  labs(x = "Time", y = expression("Net Displacement "(km)), title = "Standard") +
  theme(panel.grid = element_blank(),
        strip.text = element_text(face = "bold", size = 10),
        plot.title = element_text(hjust = 0.5, face = "bold"))

p2<- ggplot(map_dfr(dat.list3, `[`), aes(x=time1, y=displ_n)) +
  geom_path() +
  geom_vline(data = brkpts_t, aes(xintercept = value, group = id), color = "blue") +
  facet_wrap(~id, scales = "free_x", ncol = 1) + #separate IDs
  theme_bw() +
  labs(x = "Time", y = expression("Normalized Net Displacement "(km)), title = "Normalized") +
  theme(panel.grid = element_blank(),
        strip.text = element_text(face = "bold", size = 10),
        plot.title = element_text(hjust = 0.5, face = "bold"))

#plot
grid.arrange(p1, p2, ncol = 2)



#discretize displ_n for 2nd pass of model
dat_norm.list<- map(dat.list3, discrete_move_par, lims = list(displ.bin.lims), varIn = "displ_n",
                    varOut = "ND2")
dat_norm.list2<- map(dat_norm.list, subset, select = c(id, ND2))



#################################
#### Run Gibbs Sampler by ID ####
#################################

ngibbs = 5000

nbins = length(displ.bin.lims)

#prior
alpha=1

## Run Gibbs sampler
plan(multisession)  #run all MCMC chains in parallel
#refer to future::plan() for more details

dat.res2<- behavior_segment(dat = dat_norm.list2, ngibbs = ngibbs, nbins = nbins, alpha = alpha)



## Traceplots
#type is either 'nbrks' or 'LML' for y-axis label
identity<- names(dat.list)

traceplot(data = dat.res2$nbrks, type = "nbrks", identity = identity)
traceplot(data = dat.res2$LML, type = "LML", identity = identity)


##Determine maximum likelihood (ML) for selecting breakpoints
ML2<- apply(dat.res2$LML, 1, function(x) getML(dat = x, nburn = 500))
brkpts2<- getBreakpts(dat = dat.res2$brkpts, ML = ML2, identity = identity)


brkpts_t2<- brkpts2 %>%  #change from wide to long format
  pivot_longer(-id, values_to = "value", values_drop_na = TRUE)

ggplot(map_dfr(dat.list3, `[`), aes(x=time1, y=displ_n)) +
  geom_path() +
  geom_vline(data = brkpts_t, aes(xintercept = value, group = id), color = "blue") +
  geom_vline(data = brkpts_t2, aes(xintercept = value, group = id), color = "red",
             linetype = "dashed") +
  facet_wrap(~id, scales = "free_x") + #separate IDs
  theme_bw() +
  labs(x = "Time", y = expression("Net Displacement "(km))) +
  theme(panel.grid = element_blank(),
        strip.text = element_text(face = "bold", size = 10))






#############################################
#############################################
#### 2-pass model using different alphas ####
#############################################
#############################################


########################################
#### Run Gibbs Sampler w/ alpha = 2 ####
########################################

ngibbs = 5000

nbins = length(displ.bin.lims)

#prior
alpha=2

## Run Gibbs sampler
plan(multisession)  #run all MCMC chains in parallel
#refer to future::plan() for more details

dat.res3<- behavior_segment(dat = dat.list2, ngibbs = ngibbs, nbins = nbins, alpha = alpha)



## Traceplots
#type is either 'nbrks' or 'LML' for y-axis label
identity<- names(dat.list)

traceplot(data = dat.res3$nbrks, type = "nbrks", identity = identity)
traceplot(data = dat.res3$LML, type = "LML", identity = identity)


##Determine maximum likelihood (ML) for selecting breakpoints
ML3<- apply(dat.res3$LML, 1, function(x) getML(dat = x, nburn = 500))
brkpts3<- getBreakpts(dat = dat.res3$brkpts, ML = ML3, identity = identity)


brkpts_t3<- brkpts3 %>%  #change from wide to long format
  pivot_longer(-id, values_to = "value", values_drop_na = TRUE)

ggplot(map_dfr(dat.list, `[`), aes(x=time1, y=displ)) +
  geom_path() +
  geom_vline(data = brkpts_t3, aes(xintercept = value, group = id), color = "blue") +
  facet_wrap(~id, scales = "free_x") + #separate IDs
  theme_bw() +
  labs(x = "Time", y = expression("Net Displacement "(km))) +
  theme(panel.grid = element_blank(),
        strip.text = element_text(face = "bold", size = 10))






##########################################
#### Run Gibbs Sampler w/ alpha = 0.5 ####
##########################################

ngibbs = 5000

nbins = length(displ.bin.lims)

#prior
alpha=0.5

## Run Gibbs sampler
plan(multisession)  #run all MCMC chains in parallel
#refer to future::plan() for more details

dat.res4<- behavior_segment(dat = dat.list2, ngibbs = ngibbs, nbins = nbins, alpha = alpha)



## Traceplots
#type is either 'nbrks' or 'LML' for y-axis label
identity<- names(dat.list)

traceplot(data = dat.res4$nbrks, type = "nbrks", identity = identity)
traceplot(data = dat.res4$LML, type = "LML", identity = identity)


##Determine maximum likelihood (ML) for selecting breakpoints
ML4<- apply(dat.res4$LML, 1, function(x) getML(dat = x, nburn = 500))
brkpts4<- getBreakpts(dat = dat.res4$brkpts, ML = ML4, identity = identity)


brkpts_t4<- brkpts4 %>%  #change from wide to long format
  pivot_longer(-id, values_to = "value", values_drop_na = TRUE)

ggplot(map_dfr(dat.list, `[`), aes(x=time1, y=displ)) +
  geom_path() +
  geom_vline(data = brkpts_t3, aes(xintercept = value, group = id), color = "blue") +
  geom_vline(data = brkpts_t4, aes(xintercept = value, group = id), color = "red",
             linetype = "dashed") +
  facet_wrap(~id, scales = "free_x") + #separate IDs
  theme_bw() +
  labs(x = "Time", y = expression("Net Displacement "(km))) +
  theme(panel.grid = element_blank(),
        strip.text = element_text(face = "bold", size = 10))
