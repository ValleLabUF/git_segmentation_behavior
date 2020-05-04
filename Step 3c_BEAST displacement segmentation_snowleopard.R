library(Rbeast)
library(tidyverse)

# #Vignette
# data(avhrr_YellowStone)
# plot(avhrr_YellowStone,type='l')
# result=beast(avhrr_YellowStone)
# plot(result)
# dev.off()

#Load data
dat<- read.csv("Modified Snow Leopard Data.csv", as.is = T)
plot(sqrt(dat[dat$id == "Pari",]$R2n)/1000, type = "l")
pari<- sqrt(dat[dat$id == "Pari",]$R2n)/1000

#works when using obs 1-500, but not 1-1000
foo<- beast(data = pari[400:574])
plot(foo)




# res<- list()
# id<- unique(dat$id)
# for (i in 1:length(id)) {
#   res[[i]]<- beast(data = sqrt(dat[dat$id == id[i],]$R2n)/1000,
#                    option = findfrequency(sqrt(dat[dat$id == id[i],]$R2n)/1000))
# }
# 
# freq<- findfrequency(sqrt(dat[dat$id == id[1],]$R2n)/1000)
# foo<- beast(data = sqrt(dat[dat$id == id[1],]$R2n)/1000)
