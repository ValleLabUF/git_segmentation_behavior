library(dplyr)
library(lubridate)
library(ggplot2)

source('helper functions.R')


dat<- read.csv("Snail Kite Gridded Data.csv", header = T, sep = ",")
names(dat)[7]<- "dist"  #change to generic form
dat$ESTtime<- as_datetime(dat$ESTtime)
dat$date<- date(dat$ESTtime)
behav.list<- behav.prep(dat=dat, tstep = 3600)  #add move params and filter by 3600 s interval

ggplot(behav.list$`1`[1:70,]) +  #48 days betweeb multi-observation days for ID 1
  geom_path(aes(utmlong, utmlat)) +
  geom_point(aes(utmlong, utmlat), size=2) +
  facet_wrap(~date, scales = "free") +
  theme(axis.text = element_blank())


ggplot(behav.list$`12`[1:77,]) +
  geom_path(aes(utmlong, utmlat)) +
  geom_point(aes(utmlong, utmlat), size=2) +
  facet_wrap(~date, scales = "fixed") +
  theme(axis.text = element_blank())


ggplot(behav.list$`19`[1:69,]) +
  geom_path(aes(utmlong, utmlat)) +
  geom_point(aes(utmlong, utmlat), size=2) +
  facet_wrap(~date, scales = "fixed") +
  theme(axis.text = element_blank())


ggplot(behav.list$`27`[1:70,]) +
  geom_path(aes(utmlong, utmlat)) +
  geom_point(aes(utmlong, utmlat), size=2) +
  facet_wrap(~date, scales = "fixed") +
  theme(axis.text = element_blank())
