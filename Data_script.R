library('plot.matrix')
library(raster)
library(ncdf4)
library(tidyverse)
setwd("C:/Users/Hector/Documents/Pughtam-cropozone/Best Ensemble")
load("./rice_best_model_agmerra_hist_default_yield_ric_annual_1980_2010_TEST.R_shifted_ts_processed_dt.Rdata");Ricedef<-t(mapi2)
load("./maize_best_model_agmerra_hist_default_yield_mai_annual_1980_2010_TEST.R_shifted_ts_processed_dt.Rdata");Maizedef<-t(mapi2)
load("./wheat_best_model_agmerra_hist_default_yield_whe_annual_1980_2010_TEST.R_shifted_ts_processed_dt.Rdata");Wheatdef<-t(mapi2)
table(Ricedef)
table(Maizedef)
table(Wheatdef)
plot(Ricedef,col=topo.colors,na.cell=FALSE, border=NA)
plot(Maizedef,col=topo.colors,na.cell=FALSE, border=NA)
plot(Wheatdef,col=topo.colors,na.cell=FALSE, border=NA)
ensemb<-list(Wheatdef,Maizedef,Ricedef)
#reading the model outputs
crop<-c("Wheat","Maize","Rice")
ggcms <- c("pdssat","epic-boku","epic-iiasa","gepic","papsim","pegasus","lpj-guess",
           "lpjml","cgms-wofost","clm-crop","epic-tamu","orchidee-crop","pepic","prysbi2")
irrggcm_yi<-list()
nirggcm_yi<-list()
temp1<-list()
temp2<-list()
for (i in 1:length(crop)){
  for(j in 1:length(ggcms)){
    temp1[[j]]<-intersect(list.files(path=paste0("./GGCMI",crop[i]),pattern=ggcms[j]),list.files(path=paste0("./GGCMI",crop[i]),pattern="firr"))
    temp2[[j]]<-intersect(list.files(path=paste0("./GGCMI",crop[i]),pattern=ggcms[j]),list.files(path=paste0("./GGCMI",crop[i]),pattern="noirr"))
  }
  names(temp1)<-ggcms;names(temp2)<-ggcms
  temp1<-compact(temp1);temp2<-compact(temp2)
  for(j in 1:length(temp1)){temp1[[j]]<-brick(paste0("./GGCMI",crop[i],"/",temp1[[j]]))}
  for(j in 1:length(temp2)){temp2[[j]]<-brick(paste0("./GGCMI",crop[i],"/",temp2[[j]]))}
  irrggcm_yi[[i]]<-temp1
  nirggcm_yi[[i]]<-temp2
}
names(irrggcm_yi)<-crop;names(nirggcm_yi)<-crop

#extracting the extra files from prysbi2
prysbi2W<-irrggcm_yi[[1]][[13]]
prysbi2M<-irrggcm_yi[[2]][[13]]
prysbi2R<-irrggcm_yi[[3]][[11]]

for (i in 1:length(irrggcm_yi)){
irrggcm_yi[[i]]<-irrggcm_yi[[i]][1:length(irrggcm_yi[[i]])-1]
}

#calling the areas estimated from SPAM
cropAr<-c("SPAMest_Wheat_","SPAMest_Maize_","SPAMest_Rice_")
cropArirr<-c("SPAMest_Wheatirr_","SPAMest_Maizeirr_","SPAMest_Riceirr_")
ArNew<-list()
ArNewi<-list()
ArNewr<-list()
setwd("C:/Users/Hector/Documents/Pughtam-cropozone/Global_evaluation_outputs")
for (j in 1:length(cropAr)){
  ArNew[[j]]<-list.files(pattern=cropAr[j])
  ArNewi[[j]]<-list.files(pattern=cropArirr[j])
  ArNew[[j]]<-subset(brick(ArNew[[j]]),20:50)
  ArNewi[[j]]<-subset(brick(ArNewi[[j]]),20:50)*(ArNew[[j]]/ArNew[[j]])
  ArNewr[[j]]<-ArNew[[j]]-ArNewi[[j]]
}

#Production calculation by model
Wheatprod<-list()
Maizprod<-list()
Riceprod<-list()
for (j in 1 : length(irrggcm_yi[[1]])){Wheatprod[[j]]<-(nirggcm_yi[[1]][[j]]*ArNewr[[1]]+irrggcm_yi[[1]][[j]]*ArNewi[[1]])}
for (j in 1 : length(irrggcm_yi[[2]])){Maizprod[[j]]<-(nirggcm_yi[[2]][[j]]*ArNewr[[2]]+irrggcm_yi[[2]][[j]]*ArNewi[[2]])}
for (j in 1 : length(irrggcm_yi[[3]])){Riceprod[[j]]<-(nirggcm_yi[[3]][[j]]*ArNewr[[3]]+irrggcm_yi[[3]][[j]]*ArNewi[[3]])}  
Wheatprod[[13]]<-prysbi2W*ArNew[[1]];names(Wheatprod)<-c(names(irrggcm_yi[[1]]),"prysbi2")
Maizprod[[13]]<-prysbi2M*ArNew[[2]];names(Maizprod)<-c(names(irrggcm_yi[[2]]),"prysbi2")
Riceprod[[11]]<-prysbi2R*ArNew[[3]];names(Riceprod)<-c(names(irrggcm_yi[[3]]),"prysbi2")
crop_prod<-list(Wheatprod,Maizprod,Riceprod)

#Yield calculation by model
Wheatyi<-list()
Maizyi<-list()
Riceyi<-list()
for (j in 1 : length(Wheatprod)){Wheatyi[[j]]<-Wheatprod[[j]]/ArNew[[1]]};names(Wheatyi)<-names(Wheatprod)
for (j in 1 : length(Maizprod)){Maizyi[[j]]<-Maizprod[[j]]/ArNew[[2]]};names(Maizyi) <-names(Maizprod)
for (j in 1 : length(Riceprod)){Riceyi[[j]]<-Riceprod[[j]]/ArNew[[3]]};names(Riceyi) <-names(Riceprod)
crop_yi<-list(Wheatyi,Maizyi,Riceyi)


#Delete package to avoid conflick with extract of raster
#.rs.unloadPackage("tidyr")

#
ensembras<-list()
for (i in 1:3){
ensembras[[i]]<-raster(ensemb[[i]])
extent(ensembras[[i]])<-(c(-180, 180, -90, 90))
};names(ensembras)<-crop

ggcms_ <- c("pdssat","`epic-boku`","`epic-iiasa`","gepic","papsim","pegasus","`lpj-guess`",
           "lpjml","`cgms-wofost`","`clm-crop`","`epic-tamu`","`orchidee-crop`","pepic","prysbi2")
temp1<-list()
ensemYi<-list()
for (k in 1:3){
  for (i in 1:31){
  temp1[[i]]<-raster(ncols=720,nrows=360)
  for (j in 1:length(ggcms_)){
    temp1[[i]][Which(ensembras[[k]]==j,cells=TRUE)]<-eval(parse(text=paste0("crop_yi[[k]]$",ggcms_[j],"[[i]][Which(ensembras[[k]]==j,cells=TRUE)]")))
  }
    }
  ensemYi[[k]]<-stack(temp1)
}


#checkin the loop output (remember that seconf index in crop_yi is lower that 14 depending on the crop)
# check1<-extract(crop_yi[[3]][[11]][[1]],Which(ensembras[[3]]==14,cells=TRUE))
# check2<-extract(ensemYi[[3]][[1]],Which(ensembras[[3]]==14,cells=TRUE))
# check3<-check1-check2
# unique(check3)
