source ("C:/users/hac809/desktop/Best_ensemble/Corr.R")
library(raster)
library(ncdf4)
library(plyr)
library(rworldmap)
library(dplyr)
library(zoo)
require(forecast)

#Global aggregation
Global<-list()
for(j in 1:3){
  Global[[j]]<-join(Prod_LPJdf[[j]][,-1],Area_df[[j]][,-1])
  colnames(Global[[j]])<-c("Lon","Lat","Prod","Year","Area")
  Global[[j]]<-na.omit(Global[[j]])
  Global[[j]]<-aggregate(Global[[j]],by=list(Global[[j]]$Year),FUN=sum)
  Global[[j]]$Yield<-Global[[j]]$Prod/Global[[j]]$Area
  }
FGlobal<- read.csv("C:/Users/hac809/Documents/Pughtam-cropozone/Global_evaluation_inputs/FAO/FAOGl.csv", h=T)
FGlobal<-subset(FGlobal,year>1980&year<2010)

#moving average detrend
detrend<-function(x){
  dt<-ma(x,order=5,centre=TRUE)
  dt_res<- x-dt
  dt_res<-subset(dt_res,!is.na(dt_res))
}
#GLOBAL CORR
Corr_Glob<-list()
for (i in 1:3){
  Fdty<-detrend(FGlobal[,3*i+1])
  Ldty<-detrend(Global[[i]][,7])
  Corr_Glob[[i]]<-cor.test(Fdty,Ldty)
}

