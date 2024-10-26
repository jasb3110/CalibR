#prueba calibR
setwd("D:/Pepe/2022/Github/CalibR")
source("calib2.R")#Call up a function
#INPUT
#data is gonna calibrate
d2=read.csv("mejillones2.csv",sep=";",dec=".",header = TRUE)
calibR2(input=d2,sigma = 0.68,curve="marine20",show.table = F,show.plot =F,colour="no")