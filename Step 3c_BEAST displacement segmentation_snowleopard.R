library(Rbeast)
library(tidyverse)

#Load data
dat<- read.csv("Example data.csv", as.is = T)
dat<- dat$x
plot(dat, type = "l")

#Run BEAST model
dat.res<- beast(data = dat)
plot(dat.res)

