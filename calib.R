calibR=function(input=input,sigma=c(.68,.95,.99,"1s","2s","3s","1sigma","2sigma","3sigma"),curve=c(1,2,3,"intcal20", "marine20","shcal20",
                                                                                                   "marine13", "shcal13",
                                                                                                  "nh1", "nh2", "nh3", "sh1-2", "sh3", "nh1_monthly", "nh1_monthly", "nh2_monthly",
                                                                                                  "nh3_monthly", "sh1-2_monthly", "sh3_monthly", "kure", "levinKromer",
                                                                                                  "santos"),show.table=c(TRUE,FALSE,T,F,1,0),show.plot=c("minimal","default","both","FALSE",T,F,0,1)){
##################################### 
  # License: GNU
  # Author: José Solís, March 2023
  # email: solisbenites.jose@gmail.com
####################################
  require("IntCal")
  require("ggplot2")
  require("ggh4x")
  require("gridExtra")
  require("devtools")
  require("magrittr")
  require("scales")
  require("beepr")
###################################
  begin=Sys.time()#Begining time
  dd<- input[-c(1),]
  cat(paste0("To begin to calibrate of ",dd$Lab[1], sep="\n\n"))
  
  namess=c("Lab","Sample","X14C.BP","X14C.Age.SD","Lab.Error",         
           "Age.Span"           ,"Uncorrected.14C"    ,"Uncorrected.14C.SD" ,"d13C"               ,"d13C.SD",           
           "Delta.R"            ,"Delta.R.SD"         ,"Marine.Carbon"      ,"Description"        ,"CalCurve",          
           "Depth")
  
  l=c("Lab code", "Sample code",expression(phantom()^14*C~"\n(yrs BP)"),           
      expression(phantom()^14*C~SD~"\n(yrs BP)"), "Lab Error", "Age Span",          
      "Uncorrected~delta^14*C","Uncorrected~SD~delta^14*C","delta^13*C",              
      "delta^13*C~SD" ,"Delta*R~(yrs)" ,            
      "Delta*R~SD~(yrs)","Marine carbon\n(%)","Description","Calibration~curve",          
      "Depth~(cm)")
  
  ccc=c(1,2,3,"intcal20", "marine20","shcal20",
        "marine13", "shcal13",
        "nh1", "nh2", "nh3", "sh1-2", "sh3", "nh1_monthly", "nh1_monthly", "nh2_monthly",
        "nh3_monthly", "sh1-2_monthly", "sh3_monthly", "kure", "levinKromer",
        "santos")
  
  sigmas=c("1s","2s","3s","1sigma","2sigma","3sigma")
  sigmass.number=c(.68,.95,.99,.68,.95,.99)
  sigma0=as.data.frame(cbind(sigmas,sigmass.number))
  dd=as.data.frame(dd)
  labely=expression(paste("Density (",10^-3,")"))
  tabless=c(TRUE,FALSE,T,F,1,0)
  plotss=c("minimal","default","both","FALSE",T,F,0,1)
  labely=expression(paste("Density (",10^-3,")"))
  dct=getwd()#name of directory
  setwd(dct)#directory

# outcome folders
outy="calout"
if (file.exists(outy)) {
  cat("Calout folder already exists", sep="\n\n")
} else {
  dir.create(outy)
}

if(ncol(dd)!=length(namess)){
stop("Input database haven´t numbers of columns")
}else{
for(k in 1:ncol(dd)){
  if(colnames(dd)[k]!=namess[k]){
  cat(paste0("Error in colname ",k))
  stop("Input database haven´t numbers of columns")  
  }else{
    if(length(sigma)!=1){
    stop("Sigma should be one item") 
    }else{
      if(!is.character(sigma)){
        if(!is.numeric(sigma)){
        stop("Sigma should be alowed value") 
        }else{
          if(sigma<0&sigma>1){
          stop("Sigma should be a value between 0 to 1")
          }else{
          sigma=sigma
          }
        }
      }else{
          if(length(which(sigma0$sigmas!=sigma))==length(sigma0$sigmas)){
          stop("Sigma should be alowed value")  
          }else{
        s=sigma0$sigmass.number[which(sigma0$sigmas==sigma)]
        sigma=s
        }
       }
      }
    }
  }
}

if(length(show.plot)!=1){
  stop("plot value should be one item")
}else{
  if(sum(plotss!=show.plot,na.rm=T)!=7){
    stop("plot value don´t allowed, just to be minimal, default, both, FALSE, T, F, 0, 1")
  }  
}  

if(length(curve)!=1){
    stop("just picked one calibration curve")
  }else{
  if(sum(ccc!=curve,na.rm=T)!=(length(ccc)-1)){
    stop(paste0("name of calibration curve do not find, just to be intcal20, marine20, shcal20, marine13, shcal13, nh1, nh2, nh3,          
    sh1-2, sh3, nh1_monthly, nh1_monthly, nh2_monthly, nh3_monthly, sh1-2_monthly, sh3_monthly,  
    kure, levinKromer, santos"))
  }
}

#######################################  
#Input variables
sigma=sigma# 0.68, 0.95, 0.99 >>>one, two, three sigma
rrr=NULL
rr=NULL
r=NULL
sss=NULL
vvv=NULL
gg=NULL
rsv=NULL
sdrsv=NULL
c14=NULL
sdc14=NULL
dd$mean=NA
dd$lower=NA
dd$upper=NA
dd$median=NA
dd$percent=NA
dd$max=NA
dd$error=NA

for(i in 1:length(dd$Sample)){
  rsv=as.numeric(dd$Delta.R[i])
  sdrsv=as.numeric(dd$Delta.R.SD[i])
  c14=as.numeric(dd$X14C.BP[i])
  sdc14=as.numeric(dd$X14C.Age.SD[i])
  if(sum(is.null(c14)==T,is.null(sdc14)==T,na.rm=T)>0){
    next()
  }else{
    if(sum(is.na(c14)==T,is.na(sdc14)==T)>0){
      next()
    }else{
      if(sum(curve=="marine20"&((c14-sdc14-rsv-sdrsv)<603),curve=="marine20"&((c14+sdc14+rsv+sdrsv)>50000),na.rm = T)==1){ 
        next()
      }else{
        if(is.na(rsv)|is.null(rsv)==1){
          rsv=0
        }
        if(is.na(sdrsv)|is.null(sdrsv)==1){
          sdrsv=0 
        }
        
        assign(paste0("rrr",i),calibrate(age=c14, error=sdc14, cc = curve, prob=sigma, yr.steps=1, threshold=5e-04, rounded=4, reservoir=c(rsv,sdrsv),legend.cex = 1))
        invisible(get(paste0("rrr",i)))
        assign(paste0("r",i),as.data.frame(get(paste0("rrr",i))[1])[1])
        assign(paste0("rr",i),as.data.frame(get(paste0("rrr",i))[1])[2])
        dd$mean[i]=sum(as.data.frame(get(paste0("rrr",i))[1])[1]*as.data.frame(get(paste0("rrr",i))[1])[2])
        
        if(dim(get(paste0("rrr",i))[[2]])[1]==1){
          dd$lower[i]=get(paste0("rrr",i))[[2]][1]
          dd$upper[i]=get(paste0("rrr",i))[[2]][2]
          dd$median[i]=round(dd$lower[i]*.5+dd$upper[i]*.5)
          dd$percent[i]=get(paste0("rrr",i))[[2]][3]
        }else{
          gg=as.data.frame(get(paste0("rrr",i))[[2]])  
          dd$lower[i]=gg$from[which(gg$perc==max(gg$perc))]
          dd$upper[i]=gg$to[which(gg$perc==max(gg$perc))]
          dd$median[i]=round(dd$lower[i]*.5+dd$upper[i]*.5)
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
  }
}
##########################################
#create or open folder with specified name
file=paste0(dct,"/",outy)
setwd(file)
f=dd$Lab[i]
if (file.exists(f)) {
  cat("The outcome folder already exists", sep="\n\n")
} else {
  dir.create(f)
}
#outcome folder
folder=paste0(file,"/",f)
setwd(folder)
###########################################
#Plotting
c14=NULL
sdc14=NULL
rsv=NULL
sdrsv=NULL

for( i in 1:length(dd$Sample)){
  rsv=as.numeric(dd$Delta.R[i])
  sdrsv=as.numeric(dd$Delta.R.SD[i])
  c14=as.numeric(dd$X14C.BP[i])
  sdc14=as.numeric(dd$X14C.Age.SD[i])
  
  if(sum(is.null(c14)==T,is.null(sdc14)==T,na.rm=T)>0){
    warning("Can not calibrate dates: null value", sep="\n\n")
    next()
  }else{
    if(sum(is.na(c14)==T,is.na(sdc14)==T)>0){
      warning("Can not calibrate dates: NA value", sep="\n\n")
      next()
    }else{
      if(sum(curve=="marine20"&((c14-sdc14-rsv-sdrsv)<603),curve=="marine20"&((c14+sdc14+rsv+sdrsv)>50000),na.rm = T)==1){ 
  
        cat(paste0("Can´t plotted date beyond ",curve," calibration curve!: Convencial age ",dd$X14C.BP[i],"\u00B1",dd$X14C.Age.SD[i]), sep="\n\n")
        next()
      }else{
        if(is.na(rsv)|is.null(rsv)==1){
          rsv=0
        }
        if(is.na(sdrsv)|is.null(sdrsv)==1){
          sdrsv=0 
        }
        png(paste0(dd$Sample[i],".plot.png"),width = 200, heigh = 200, units = 'mm', res =1200)
        calibrate(age=c14, error=sdc14, cc = curve, prob=sigma, yr.steps=1, threshold=5e-04, rounded=4, reservoir=c(rsv,sdrsv),legend.cex = 1)
        dev.off()
        
        if(sum(show.plot=="both",show.plot=="default",show.plot==T,show.plot==1,na.rm=T)==1){
          x11()
          plot.new()
          calibrate(age=c14, error=sdc14, cc = curve, prob=sigma, yr.steps=1, threshold=5e-04, rounded=4, reservoir=c(rsv,sdrsv),legend.cex = 1)
          dev.off()
          }
        
        dr=as.data.frame(get(paste0("rrr",i))[1])
        dr[[2]]=1000*dr[[2]]
        
        if(sum(is.na(dd$Depth[i]),is.null(dd$Depth[i]),dd$Depth[i]=="",na.rm = T)>0){
          label.name=paste0(dd$Lab[i],"-",dd$Sample[i]) 
        }else{
          label.name=paste0(dd$Lab[i],"-",dd$Sample[i]," sample on ",dd$Depth[i],"cm")  
        }

        if (is.numeric(curve)) {
        if (curve==1) 
          curve="intcal20"
        else if (curve==2) 
          curve="marine20"
        else if (curve==3) 
          curve="shcal20"
        }
        
        plotting=ggplot(data=dr, aes(x=dr[[1]], y=dr[[2]])) + geom_line()+
          theme_classic()+
          
          geom_segment(aes(y =0,
                           yend =dr$V2[which(dr$cal.BP==trunc(dd$mean[i]))][1],
                           x=dd$mean[i] ,
                           xend=dd$mean[i] ),color = "gray", size=.5)+
          
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
          
          
          geom_segment(aes(y =0,
                           yend =dr$V2[which(dr$cal.BP==trunc(dd$median[i]))][1],
                           x=dd$median[i] ,
                           xend=dd$median[i] ),color = "green", size=.5)+
          
          geom_segment(aes(y =0,
                           yend =dr$V2[which(dr$cal.BP==trunc(dd$max[i]))][1],
                           x=dd$max[i] ,
                           xend=dd$max[i] ),color = "gray30", size=.5)+
          
          annotate("text",x=dd$mean[i],y=quantile(dr$V2)[2],label="Mean",size = 3,col="gray40",angle='90')+
          annotate("text",x=dd$median[i],y=quantile(dr$V2)[3],label="Median",size = 3,col="#66CC99",angle='90')+
          annotate("text",x=dd$max[i],y=quantile(dr$V2)[4],label="Maximum",size = 3,col="black",angle='90')+
          
          scale_x_continuous(limits = c(min(dr$cal.BP)*.99,max(dr$cal.BP)*1.01),breaks =scales::pretty_breaks(n = 5),guide = "axis_minor")+
          scale_y_continuous(limits = c(0,max(dr$V2)*1.02),breaks =scales::pretty_breaks(n = 5),guide = "axis_minor")+
          annotate("text",x=quantile(dr[[1]])[4]*.506+quantile(dr[[1]])[5]*.505,y=quantile(dr[[2]])[5]*1.0*1.02,label=label.name, size = 4,col="black")+
          annotate("text",x=quantile(dr[[1]])[4]*.506+quantile(dr[[1]])[5]*.505,y=quantile(dr[[2]])[5]*.97*1.02,label=paste0("Cal.age: ",dd$max[i],"\u00B1",dd$error[i]), size = 4,col="black")+
          annotate("text",x=quantile(dr[[1]])[4]*.506+quantile(dr[[1]])[5]*.505,y=quantile(dr[[2]])[5]*.94*1.02,label=paste0("\u0394R = ",rsv,"\u00B1",sdrsv), size = 4,col="black")+
          annotate("text",x=quantile(dr[[1]])[4]*.506+quantile(dr[[1]])[5]*.505,y=quantile(dr[[2]])[5]*.91*1.02,label=paste0("Probability: ",trunc(100*dd$percent[i])/100,"%"), size = 4,col="black")+
          annotate("text",x=quantile(dr[[1]])[4]*.506+quantile(dr[[1]])[5]*.505,y=quantile(dr[[2]])[5]*.88*1.02,label=paste0("Cal. curve: ",curve), size = 4,col="black")+
          labs(title=paste0("Relative probability of sample "),x ="Cal yr BP", y = labely)
        theme(axis.ticks.length = unit(0.25,"cm"),ggh4x.axis.ticks.length.minor = rel(0.5),axis.ticks = element_line(size = 2),ggh4x.axis.ticks.length.minor = rel(0.5),axis.text.x=element_text(size=11,colour = "black",face="bold",angle=45, hjust=1),axis.text.y=element_text(size=11,colour = "black",face="bold",hjust=1),
              axis.title=element_text(size=14,face="bold"),title = element_text(size=16,colour = "black",face="bold"))
        
        ggsave(paste0(label.name,".png"), dpi = 900,   width = 250,
               height = 159,unit="mm",plot =plotting)
        
        if(sum(show.plot=="minimal",show.plot=="both",show.plot==T,show.plot==1,na.rm=T)==1){
          x11()
          plot.new()
          plotting
          dev.off()
        }else{
          warning("Error in plotting of Calibration curve")
        }
        }
      }    
    }
  }

################################################################################
#OUTCOME

dd$mean=round(dd$mean)
write.csv(dd,paste0(dd$`Lab`[i],".calibrated.csv"),sep=",",dec=".",col.names = TRUE)

#plot input table
d1<- input[-c(1),]
colnames(d1)=l
myt=ttheme_minimal(base_size = 12, base_colour = "black", base_family = "",
                   parse = T , padding = unit(c(2,2), "mm"),colhead=list(fg_params = list(parse=TRUE), fontface=4L,bg_params=list(fill="gray90")))

png(paste0(dd$`Lab`[1],".input.png"), width = 20+ncol(d1)*425/15, heigh = 20+100/19*nrow(d1), units = 'mm', res =1200)
grid.table(d1,rows = NULL,theme=myt)
dev.off()

colnames(dd)[1:length(l)]=l
#plot output table
dd$percent=round(dd$percent)
png(paste0(dd$`Lab code`[1],".output.png"), width = 20+ncol(dd)*525/22, heigh = 20+100/19*nrow(dd), units = 'mm', res =1200)
grid.table(dd,rows = NULL,theme=myt)
dev.off()

if(length(show.table)!=1){
  stop("table value should be one item")
}else{
if(sum(tabless!=show.table,na.rm=T)!=3){
  stop("Table value don´t allowed, just to be TRUE,FALSE,T,F,1 and 0")
}else{  
  if(sum(show.table==TRUE,show.table==T,show.table==1,na.rm=T)==3){
  x11();grid.table(d1,rows = NULL,theme=myt)
  x11();grid.table(dd,rows = NULL,theme=myt)
  }else{
    warning("Error in plotting of table")
  }
  }  
 }
setwd(dct)
end=Sys.time()#ending time
tem=end-begin
cat("Working time was estimated how", sep="\n\n")
print(tem)
cat(paste0("Calibration finished of ",dd$`Lab code`[i]," successfully!!", sep="\n\n"))
beep(8)#mario bros sound
}