# CalibR

## Contents
- [Proposal](#proposal)
- [R code](#r-code)
- [Reference](#reference)


## Proposal 
This script tries to replicate the radiocarbon calibration of software Calib 8.02 [(1)](#reference)for any Operative system in R programing. therefore, the Data input is the same as Calib 8.02 input. Also, it is able to plot and create output data in a file in the format .csv. Bellow, I attached this script and example of published data Guiñez et al, 2014 [(2)](#reference).

## R code
```markdown
#########################################################################
#datation
dct=" ~/directory.file/"
setwd(dct)#directory
library("IntCal")
library("ggplot2")
library("ggh4x")
library("beepr")

# outcome folders
outy="calout"
if (file.exists(outy)) {
  cat("The folder already exists")
} else {
  dir.create(outy)
}

#INPUT
#data is gonna calibrate
d=read.csv("mejillones.csv",sep=";",dec=".",header = TRUE)#example 

################################
dd<- d[-c(1),]

a=0.68# 0.68, 0.95, 0.99 >>>one, two, three sigma
rrr=NULL
rr=NULL
r=NULL
sss=NULL
vvv=NULL
dd$max=NULL
gg=NULL
curv="marine20"

rsv=NA
for(i in 1:length(dd$Sample)){
  rsv=dd$Delta.R[i]
 if(as.numeric(dd$X14C.BP[i])+as.numeric(dd$X14C.Age.SD[i])+as.numeric(dd$Delta.R[i])+as.numeric(dd$Delta.R.SD[i])>50000|
    as.numeric(dd$X14C.BP[i])-as.numeric(dd$Delta.R[i])<603){
  
   if(as.numeric(dd$X14C.BP[i])+as.numeric(dd$X14C.Age.SD[i])+as.numeric(dd$Delta.R[i])+as.numeric(dd$Delta.R.SD[i])>50000){
     rsv=49999-as.numeric(dd$X14C.BP[i])
    dd$mean[i]=NA
    dd$lower[i]=NA
    dd$upper[i]=NA
    dd$median[i]=NA
    dd$max[i]=NA
    dd$error[i]=NA
    print(paste0("Can not calibrate dates beyond calibration curve!: row ",i))
    next()
   }
    
   if(as.numeric(dd$X14C.BP[i])-as.numeric(dd$Delta.R[i])<603){
     rsv=as.numeric(dd$X14C.BP[i])-604
     
     dd$mean[i]=NA
     dd$lower[i]=NA
     dd$upper[i]=NA
     dd$median[i]=NA
     dd$max[i]=NA
     dd$error[i]=NA
     print(paste0("Can not calibrate dates beyond calibration curve!: row ",i))
     next()
   } 
   
  }else{
  assign(paste0("rrr",i),calibrate(age =as.numeric(dd$X14C.BP[i]), error =as.numeric(dd$X14C.Age.SD[i]),cc =curv,prob =a ,yr.steps =1,threshold = 5e-04,rounded=2,reservoir=c(as.numeric(rsv),as.numeric(dd$Delta.R.SD[i]))))
  assign(paste0("r",i),as.data.frame(get(paste0("rrr",i))[1])[1])
  assign(paste0("rr",i),as.data.frame(get(paste0("rrr",i))[1])[2])
  dd$mean[i]=sum(as.data.frame(get(paste0("rrr",i))[1])[1]*as.data.frame(get(paste0("rrr",i))[1])[2])
  
  if(dim(get(paste0("rrr",i))[[2]])[1]==1){
  dd$lower[i]=get(paste0("rrr",i))[[2]][1]
  dd$upper[i]=get(paste0("rrr",i))[[2]][2]
  dd$median[i]=dd$lower[i]*.5+dd$upper[i]*.5
  dd$percent[i]=get(paste0("rrr",i))[[2]][3]
  }else{
  gg=as.data.frame(get(paste0("rrr",i))[[2]])  
  dd$lower[i]=gg$from[which(gg$perc==max(gg$perc))]
  dd$upper[i]=gg$to[which(gg$perc==max(gg$perc))]
  dd$median[i]=dd$lower[i]*.5+dd$upper[i]*.5
  dd$percent[i]=gg$perc[which(gg$perc==max(gg$perc))]
  }
  if(length(which(get(paste0("rr",i))==max(get(paste0("rr",i)))))==1){
    dd$max[i]=get(paste0("r",i))[1][which(get(paste0("rr",i))[,1]==max(get(paste0("rr",i))[,1])),]
  }else{
    vvv=get(paste0("r",i))[which(get(paste0("rr",i))==max(get(paste0("rr",i)))),]
    sss= abs(vvv-dd$mean[i])
    dd$max[i]= vvv[which(sss==min(sss))]
  }
  dd$error[i]=max(abs(dd$lower[i]-dd$max[i]),abs(dd$max[i]-dd$upper[i]))
  }
}  

#plotting in low resolution 
#for( i in 1:length(dd$Sample)){  
#X11();plot(as.data.frame(get(paste0("rrr",i))[1])[[1]],as.data.frame(get(paste0("rrr",i))[1])[[2]] ,type="l",xlab="Cal BP",ylab="Density",main =dd$Sample[i])
#  abline(v=dd$mean[i],col="gray")#mean value
#  abline(v=dd$lower[i],col="blue")# lower value
#  abline(v=dd$upper[i],col="red")#upper value
#  abline(v=dd$median[i],col="green")#median value
#  abline(v=dd$max[i],col="black")#maximum probability value
#}  
###########################################
#create or open folder with specified name
file=paste0(dct,"/",outy)
setwd(file)
f=dd$Lab[i]
if (file.exists(f)) {
  cat("The folder already exists")
} else {
  dir.create(f)
}
#outcome folder
folder=paste0(file,"/",f)
setwd(folder)
###########################################

labely=expression(paste("Density (",10^-3,")"))
for( i in 1:length(dd$Sample)){
  dr=as.data.frame(get(paste0("rrr",i))[1])
  dr[[2]]=1000*dr[[2]]
  if(is.na(dd$mean[i])){
    print(paste0("Can not plotted dates beyond calibration curve!: Convencial age ",dd$X14C.BP[i],"\u00B1",dd$X14C.Age.SD[i]))
  }else{
  plot=ggplot(data=dr, aes(x=dr[[1]], y=dr[[2]])) + geom_line()+
    theme_classic()+
    
    geom_vline(xintercept =dd$mean[i] ,color = "gray", size=.5)+
    
    geom_segment(aes(y =0,
                     yend =dr$V2[which(dr$cal.BP==dd$lower[i])][1],
                     x=dd$lower[i],
                     xend=dd$lower[i]),color = "blue", size=.5)+
    
    geom_segment(aes(y =dr$V2[which(dr$cal.BP==dd$lower[i])][1],
                     yend =dr$V2[which(dr$cal.BP==dd$lower[i])][1],
                     x=max(dr$cal.BP),
                     xend=dd$lower[i]),color = "blue", size=.5)+
    annotate("text",x=max(dr$cal.BP)*.5+dd$lower[i]*.5,y=dr$V2[which(dr$cal.BP==dd$lower[i])][1]*.95,label="Lower",size = 3,col="blue")+
 
    geom_segment(aes(y =0,
                     yend =dr$V2[which(dr$cal.BP==dd$upper[i])][1],
                     x=dd$upper[i],
                     xend=dd$upper[i]),color = "red", size=.5)+
    
    geom_segment(aes(y =dr$V2[which(dr$cal.BP==dd$upper[i])][1],
                     yend =dr$V2[which(dr$cal.BP==dd$upper[i])][1],
                     x=min(dr$cal.BP),
                     xend=dd$upper[i]),color = "red", size=.5)+
    
    annotate("text",x=min(dr$cal.BP)*.5+dd$upper[i]*.5,y=dr$V2[which(dr$cal.BP==dd$upper[i])][1]*.95,label="Upper",size = 3,col="red")+
    
    geom_vline(xintercept =dd$median[i] ,color = "green", size=.5)+
    geom_vline(xintercept =dd$max[i] ,color = "gray30", size=.5)+
    annotate("text",x=dd$mean[i],y=quantile(dr$V2)[2],label="Mean",size = 3,col="gray40",angle='90')+
    annotate("text",x=dd$median[i],y=quantile(dr$V2)[3],label="Median",size = 3,col="#66CC99",angle='90')+
    annotate("text",x=dd$max[i],y=quantile(dr$V2)[4],label="Maximum",size = 3,col="black",angle='90')+
    
    scale_x_continuous(limits = c(min(dr$cal.BP),max(dr$cal.BP)),breaks =scales::pretty_breaks(n = 5),guide = "axis_minor")+
    scale_y_continuous(limits = c(0,max(dr$V2)),breaks =scales::pretty_breaks(n = 5),guide = "axis_minor")+
    annotate("text",x=quantile(dr[[1]])[4]*.5+quantile(dr[[1]])[5]*.5,y=quantile(dr[[2]])[5]*.99,label=paste0(dd$Lab[i]," sample on ",dd$Depth[i],"cm"), size = 5,col="black")+
    annotate("text",x=quantile(dr[[1]])[4]*.5+quantile(dr[[1]])[5]*.5,y=quantile(dr[[2]])[5]*.94,label=paste0("Cal.age: ",dd$max[i],"\u00B1",dd$error[i]), size = 5,col="black")+
    annotate("text",x=quantile(dr[[1]])[4]*.5+quantile(dr[[1]])[5]*.5,y=quantile(dr[[2]])[5]*.89,label=paste0("\u0394R = ",dd$Delta.R[i],"\u00B1",dd$Delta.R.SD[i]), size = 5,col="black")+
    annotate("text",x=quantile(dr[[1]])[4]*.5+quantile(dr[[1]])[5]*.5,y=quantile(dr[[2]])[5]*.84,label=paste0("Cal. curve: ",curv), size = 5,col="black")+
    labs(title=paste0("Relative probability of sample "),x ="Cal yr BP", y = labely)
    theme(xis.ticks.length=unit(0.25,"cm"),ggh4x.axis.ticks.length.minor = rel(0.5),axis.ticks = element_line(size = 2),ggh4x.axis.ticks.length.minor = rel(0.5),axis.text.x=element_text(size=11,colour = "black",face="bold",angle=45, hjust=1),axis.text.y=element_text(size=11,colour = "black",face="bold",hjust=1),
          axis.title=element_text(size=14,face="bold"),title = element_text(size=16,colour = "black",face="bold"))
  
    ggsave(paste0("Core ",dd$Lab[i]," sample on ",dd$Depth[i],"cm.png"), dpi = 900,   width = 250,
           height = 159,unit="mm",plot =plot)  
    }
}    
#View(dd)
#OUTCOME
write.csv(dd,paste0(dd$Lab[i],".calibrated.csv"),sep=",",dec=".",col.names = TRUE)
print(paste0("Calibration finished of ",dd$Lab[i],"!!!"))
beep(8)#mario bross sound
################################################################################
```

## Reference

Stuiver, M., & Reimer, P. J. (1993). EXTENDED 14C DATA BASE AND REVISED CALIB 3.014C AGE CALIBRATION PROGRAM. Radiocarbon, 35(1), 215–230. https://doi.org/10.14210/bjast.v17.n2.pNB5-8

Guiñez, M., Valdés, J., Sifeddine, A., Boussafir, M., & Dávila, P. M. (2014). Anchovy population and ocean-climatic fluctuations in the Humboldt Current System during the last 700 years and their implications. Palaeogeography, Palaeoclimatology, Palaeoecology, 415, 210–224. https://doi.org/10.1016/j.palaeo.2014.08.026
