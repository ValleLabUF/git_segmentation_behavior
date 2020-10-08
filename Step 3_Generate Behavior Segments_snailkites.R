
### Updated Snail Kite Segmentation (bayesmove)

library(bayesmove)
library(tidyverse)
library(future)
library(lubridate)
library(ggforce)


dat<- read.csv("Snail Kite Gridded Data_TOHO.csv", header = T, sep = ",")
dat$date<- dat$date %>% as_datetime()
# dat$rel.angle<- abs(dat$rel.angle)  #take abs value of angle

#if dt within 5 min of 1 hr, round to 1 hr
dat<- round_track_time(dat = dat, id = "id", int = 3600, tol = 300, time.zone = "UTC")
dat.list<- df_to_list(dat = dat, ind = "id")

# Filter observations
behav.list<- filter_time(dat.list = dat.list, int = 3600)

#define bin number and limits for step lengths and turning angles
angle.bin.lims=seq(from=-pi, to=pi, by=pi/4)  #8 bins

dist.bin.lims=quantile(dat[dat$dt == 3600,]$dist,
                       c(0,0.25,0.50,0.75,0.90,1), na.rm=T) #5 bins


#Viz limits on continuous vars
behav.df<- map_dfr(behav.list, `[`)

ggplot(behav.df, aes(x=dist/1000)) +
  geom_density(fill = "lightblue") +
  # xlim(0,5) +  #limit to only 5 km
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
behav.list<- map(behav.list,
                 discrete_move_var,
                 lims = list(dist.bin.lims, angle.bin.lims),
                 varIn = c("dist", "rel.angle"),
                 varOut = c("SL", "TA"))

behav.list<- behav.list[sapply(behav.list, nrow) > 2]  #remove IDs w/ fewer than 3 obs
behav.list2<- map(behav.list,
                  subset,
                  select = c(id, SL, TA))



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

set.seed(1)

# Define model params
alpha<- 1
ngibbs<- 40000
nbins<- c(5,8)

## Run Gibbs sampler
plan(multisession)  #run all MCMC chains in parallel
                    #refer to future::plan() for more details
dat.res<- segment_behavior(data = behav.list2[21], ngibbs = ngibbs, nbins = nbins,
                           alpha = alpha)
future:::ClusterRegistry("stop")  #close all threads and memory used
###Takes 25 min to run 40000 iterations for 26 IDs



## Traceplots
traceplot(data = dat.res$nbrks, ngibbs = ngibbs, type = "nbrks")
traceplot(data = dat.res$LML, ngibbs = ngibbs, type = "LML")


##Determine maximum likelihood (ML) for selecting breakpoints
MAP.est<- get_MAP(dat = dat.res$LML, nburn = ngibbs/2)
brkpts<- get_breakpts(dat = dat.res$brkpts, MAP.est = MAP.est)


# Plot breakpoints over the data
behav.list<- lapply(behav.list, function(x) x %>% 
                      mutate(dist = dist/1000))
plot_breakpoints(data = lapply(behav.list[21], function(x) x[1:2500,]),
                 as_date = FALSE, var_names = c("dist","rel.angle"),
                 var_labels = c('Step Length (km)', 'Turning Angle (rad)'),
                 brkpts = brkpts[,1:27])

# plot_breakpoints(data = behav.list[21], as_date = FALSE, var_names = c("SL","TA"),
#                  var_labels = c('SL Bins', 'TA Bins'), brkpts = brkpts)

# ggsave("Figure 6a (segmentation heatmap SNIK_12).png", width = 7, height = 5, units = "in",
#        dpi = 330)


################################################
#### Compare MAP Estimate Against Posterior ####
################################################


### Possibly include as an additional function in bayesmove
# par(ask = T)
# for (i in 1:nrow(brkpts)) {
#   brks.df<- data.frame(brks = unlist(dat.res$brkpts[[i]]))
#   brks_map<- brkpts[i, -1] %>% 
#     purrr::discard(is.na) %>% 
#     t() %>% 
#     data.frame()
#   names(brks_map)<- "brks"
#   
#   print(
#   ggplot() +
#     geom_density(data = brks.df, aes(brks), fill = "grey50", alpha = 0.5, size = 1) +
#     geom_vline(data = brks_map, aes(xintercept = brks), color = "blue") +
#     labs(x = "Breakpoints", y = "Density", title = brkpts$id[i]) +
#     theme_bw() +
#     theme(axis.text = element_text(size = 12),
#           axis.title = element_text(size = 16))
#   )
# }
# par(ask = F)


brks.df<- list()  #create empty list to store posterior samples
brks_map<- list() #create empty list to store MAP estimates
nburn<- ngibbs/2

for (i in 1:nrow(brkpts)) {
  brks.df[[i]]<- data.frame(id = names(dat.res$brkpts)[i],
                            brks = unlist(dat.res$brkpts[[i]])[(nburn+1):ngibbs])
  
  tmp<- brkpts[i, -1] %>%
    purrr::discard(is.na) %>%
    t() %>%
    data.frame() %>% 
    mutate(id = brkpts$id[i])
  names(tmp)[1]<- "brks"
  brks_map[[i]]<- tmp
}

brks.df<- bind_rows(brks.df)
brks_map<- bind_rows(brks_map)


ggplot() +
  geom_density(data = brks.df, aes(brks), fill = "grey50", alpha = 0.5, size = 1, adjust = 1/2) +
  geom_vline(data = brks_map, aes(xintercept = brks), color = "blue") +
  labs(x = "Breakpoints", y = "Density") +
  theme_bw() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 16),
        strip.text = element_text(size = 10, face = "bold")) +
  facet_wrap(~id, scales = "free", ncol = 5)


# ggsave("Compare MAP vs Posterior Breakpoints.png", width = 13, height = 9, units = "in",
#        dpi = 300)  
