#prueba calibR
setwd("E:/Pepe/2022/Github/CalibR")
source("calib4.R")#Call up a function
#INPUT
#data is gonna calibrate
d2=read.csv("mejillones2.csv",sep=";",dec=".",header = TRUE)
calibR4(input=d2,sigma = 0.68,curve="marine20",show.table = F,show.plot =F,colour="no")