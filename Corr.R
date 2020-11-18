source ("C:/Users/hac809/Desktop/Best_ensemble/Data_script.R")
library(ncdf4)
library(raster)
library(rworldmap)
library(dplyr)
library(plyr)
library(zoo)
require(forecast)
#Calling the yield output of LPJ-GUESS
setwd("C:/Users/hac809/Documents/Pughtam-cropozone/Best_Ensemble")
yield=read.table('yield.out', h=TRUE)
countries<-as.data.frame(gridCountriesDegreesHalf)
names(countries)<-c("UN","Lon","Lat")
yield<-join(yield,countries)
#Selecting cultivar by country based on HDI (Human Development Index)
#MaizCult<-read.csv("C:/Users/hac809/Documents/Pughtam-cropozone/Global_evaluation_inputs/FAO/MaizCultBC.csv",h=T)
#MaizCult<-MaizCult[c(2,4)]
MaizCult<-read.csv("C:/Users/hac809/Documents/Pughtam-cropozone/Global_evaluation_inputs/FAO/FAO_MaizBC.csv",h=T)
#MaizCult$Cultivar<-ifelse(MaizCult$FYield>5,1,2)
MaizCult$Cultivar<-ifelse(MaizCult$FYield>0.055*(MaizCult$Year-1960)+4.15,1,2)
MaizCult<-MaizCult[c(1,3,7)]
MaizCult<-subset(MaizCult,!is.na(UN)&UN!=-99)
MaizCult<-distinct(MaizCult)
#yield<-yield[,-16]
yield<-join(yield,MaizCult)
yield2=yield
yield2=na.locf(yield2,fromLast = TRUE)
yield2$TeCo<-ifelse(yield2$Cultivar==2,yield2$TeCG,yield2$TeCS)
yield2$TeCoi<-ifelse(yield2$Cultivar==2,yield2$TeCGi,yield2$TeCSi)
yield2<-yield2[,c(1,2,3,15,4,5,17,8,9,10,18,13,14)]
#yield2$TeW<-pmax(yield2$TeWW,yield2$TeSW)
#yield2$TeWi<-pmax(yield2$TeWWi,yield2$TeSWi)


vect<-c(1980:2010)
li<-list()
for (i in 1:length(vect)){
  li[[i]]<-subset(yield2,Year==vect[i])
  li[[i]]<-li[[i]][c(-3)]
}

#crops<-names(li[[i]])[c(11,10,5,3,7,6,4)]
crops<-names(li[[i]])[c(8,4,9,5,10,6,11,12,7)]
bb <- extent(-180, 180, -90, 90)
rasterize<-function(coord,value){
  areas<-SpatialPointsDataFrame(coords=coord,data.frame(value),proj4string = CRS("+proj=longlat +datum=WGS84"))
  areas2<-SpatialPixelsDataFrame(areas,tolerance=0.0192077,areas@data) 
  areas3<-as(areas2,"SpatialGridDataFrame")
  Rstr<-raster(areas3)
}
rastcrop<-list()
for(j in 1:length(crops)){
  crop_r<-list()
  for (i in 1:length(vect)){
    yield_coords<-cbind(li[[i]]$Lon,li[[i]]$Lat)
    crop_r[[i]]<-rasterize(yield_coords,li[[i]][crops[j]])
  }
  names(crop_r)<-paste0("Yi_",vect)
  rastcrop[[j]]<-stack(crop_r)
  rastcrop[[j]]<-extend(rastcrop[[j]],bb)
}
names(rastcrop)<-paste0(crops)

#Calling Ray files, if not available see the end of this script to call the original files
crop<-c("Wheat","Maize","Rice")
cropR<-c("Ray_Yield_Wheat","Ray_Yield_Maize_","Ray_Yield_Rice_")
crop_rayR<-list()
Prod_ray<-list()
setwd("C:/Users/hac809/Documents/Pughtam-cropozone/Global_evaluation_outputs")
for (j in 1:length(crop)){
  crop_rayR[[j]]<-list.files(pattern=cropR[j])
  crop_rayR[[j]]<-subset(brick(crop_rayR[[j]]),11:41)
  Prod_ray[[j]]<-ArNew[[j]]*crop_rayR[[j]]
}

names(crop_rayR)<-paste0("Ray",crop)
names(Prod_ray)<-paste0("RayP",crop)


#multiplied by 10 to convert yield from kg/m2 to Ton/ha production(ha)
#Divided by 0.88, 0.87, 0.91 since FAO production is assumed to have 12% water content
#Wheatprod<-((rastcrop[[2]]*ArNewr[[1]]+rastcrop[[1]]*ArNewi[[1]])*10/0.88)
WWprod<-((rastcrop[[2]]*ArNewr[[1]]+rastcrop[[1]]*ArNewi[[1]])*10/0.88)
SWprod<-((rastcrop[[4]]*ArNewr[[1]]+rastcrop[[3]]*ArNewi[[1]])*10/0.88)
#Corr by GR for wheat
WWcorr<- stack(WWprod,crop_rayR[[1]])
WWcorr1<- calc(WWcorr, fun=function(x) cor(x[1:31], x[32:62], method='pearson'))
SWcorr<- stack(SWprod,crop_rayR[[1]])
SWcorr1<- calc(SWcorr, fun=function(x) cor(x[1:31], x[32:62], method='pearson'))
Scomp<-overlay(SWcorr1,WWcorr1,fun=function(x,y) {ifelse(x>y,1,0)})
Wcomp<-overlay(WWcorr1,SWcorr1,fun=function(x,y) {ifelse(x>y,1,0)})
Wheatprod<-WWprod*Wcomp+SWprod*Scomp
Maizprod<-((rastcrop[[6]]*ArNewr[[2]]+rastcrop[[5]]*ArNewi[[2]])*10/0.88)
Riceprod<-((rastcrop[[8]]*ArNew[[3]])*10/0.87)
ProdBG<-list(Wheatprod,Maizprod,Riceprod)

#Selecting year lag based on correlation
detrend<-function(x){
  dt<-ma(x,order=5,centre=TRUE)
  dt_res<-x-dt
  return(dt_res)
}
ProdBG1<-list()
YielBG<-list()
for (i in 1:length(ProdBG)){
  temp1<-calc(ProdBG[[i]],detrend)
  temp1<-temp1[[-c(1,2,nlayers(temp1)-1,nlayers(temp1))]]
  raytest<-calc(crop_rayR[[i]],detrend)
  raytest<-raytest[[-c(1,2,nlayers(raytest)-1,nlayers(raytest))]]
  temp1<-stack(temp1,raytest)
  Corr<-calc(temp1, fun=function(x) cor(x[1:27], x[28:54], method='pearson'))
  Corr1<-calc(temp1, fun=function(x) cor(x[2:27], x[28:53], method='pearson'))
  Corr2<-calc(temp1, fun=function(x) cor(x[1:26], x[29:54], method='pearson'))
  Map1<-overlay(Corr,Corr1,Corr2, fun=function(x,y,z) {ifelse(x>y&x>z,1,0)})
  Map2<-overlay(Corr1,Corr2,Corr, fun=function(x,y,z) {ifelse(x>y&x>z,1,0)})
  Map3<-overlay(Corr2,Corr,Corr1, fun=function(x,y,z) {ifelse(x>y&x>z,1,0)})
  ProdBG1[[i]]<-ProdBG[[i]][[c(2:30)]]*Map1+ProdBG[[i]][[c(3:31)]]*Map2+ProdBG[[i]][[c(1:29)]]*Map3
  YielBG[[i]]<-ProdBG1[[i]]/ArNew[[i]][[c(2:30)]]
}


#Transforming LPJ and Ensemble data to data frames
cut<-function(x){
  cu<-x[c(-1,-2,-(length(x$year)),-(length(x$year)-1)),]
  return(data.frame(cu))
}

Deraster<-function(x){
  df<-list()
  for(i in 1:nlayers(x)){
    df[[i]]<-as.data.frame(x[[i]],xy=TRUE)
    df[[i]]$year<-1980+i
  }
  df<-lapply(df, setNames, c("Lon","Lat","Yield","Years"))
  return(bind_rows(df, .id = "label"))
}
ensemYi1<-list();for (i in 1:3){ensemYi1[[i]]<-ensemYi[[i]][[c(2:30)]]}

yield_LPJdf<-lapply(YielBG,FUN=Deraster)
Prod_LPJdf<-lapply(ProdBG1,FUN=Deraster)
yield_Ensemdf<-lapply(ensemYi1,FUN=Deraster)
Prod_Ensemdf<-lapply(ensemPro,FUN=Deraster)
Area_df<-lapply(ArNew,FUN=Deraster);
Area_df<-lapply(Area_df,setNames,c("label","Lon","Lat","Area","Years"))
Prod_Ensemdf1<-list()

for (i in 1:3){
  colnames(Prod_Ensemdf[[i]])[4]<-"Prod"
  Prod_Ensemdf1[[i]]<-subset(Prod_Ensemdf[[i]],Years==2010)
  Prod_Ensemdf1[[i]]$Prod2<-Prod_Ensemdf1[[i]]$Prod/1000000
}

y_ext<-list()
yield_bothdf<-list()
yearsg_ma<-list()
for (i in 1:3){
  y_ext[[i]]<-cbind(yield_LPJdf[[i]][,-1],Ensem=yield_Ensemdf[[i]][,4])
  y_ext[[i]]$label<-rep(1:259200,29)
  yield_bothdf[[i]]<-na.omit(y_ext[[i]])
  yield_bothdf[[i]]<-yield_bothdf[[i]][yield_bothdf[[i]]$label %in% names(which(table(yield_bothdf[[i]]$label)>7)), ]
  yearsg_ma[[i]]<-yield_bothdf[[i]][,c(4,6)]
  colnames(yearsg_ma[[i]])<-c("year","label")
  yearsg_ma[[i]]<-ddply(yearsg_ma[[i]], .(label),cut)
  yield_bothdf[[i]]<-yield_bothdf[[i]][,c(6,4,1,2,3,5)]
}


#DETREND
dtcropbg<-list()
for (j in 1:3){
  dtbg<-list()
  for(i in 1:2){
    detrendbg<-function(x){
      dt<-ma(x[,i+4],order=5,centre=TRUE)
      dt_res<- x[,i+4]-dt
      dt_res<-subset(dt_res,!is.na(dt_res))
      return(data.frame(dt_res))
    }
    dtbg[[i]]<-ddply(yield_bothdf[[j]], .(label),detrendbg)
  }
  dtcropbg[[j]]<-join(yearsg_ma[[j]],do.call("cbind",dtbg)[c(-3)])
  colnames(dtcropbg[[j]])<-c("year","label","LPJ","Ensem")
}
names(dtcropbg)<-crop

#CORRELATION By GC
corrfunG<-function(x){
  COR=cor.test(x$LPJ,x$Ensem)
  n1<-length(x$LPJ)
  n2<-length(x$Ensem)
  return(data.frame(COR$estimate,COR$p.value))
}

corrbg<-list()
corrbg1<-list()
#corrbg2<-list()
#corrindex<-list()
Cropc<-list()
y_exts<-list()
for ( i in 1:3){
  corrbg[[i]]<-ddply(dtcropbg[[i]], .(label),corrfunG)
  corrbg1[[i]]<-corrbg[[i]]
  corrbg1[[i]]$COR.estimate<-ifelse(corrbg1[[i]]$COR.p.value>0.05,corrbg1[[i]]$COR.estimate==0,corrbg1[[i]]$COR.estimate)
  y_exts[[i]]<-subset(y_ext[[i]],Years==1981)
  for(j in 1:nrow(corrbg1[[i]])){
    corrbg1[[i]]$Lon[j]<-y_exts[[i]][which(y_exts[[i]]$label==corrbg1[[i]]$label[j]),1]
    corrbg1[[i]]$Lat[j]<-y_exts[[i]][which(y_exts[[i]]$label==corrbg1[[i]]$label[j]),2]
    corrbg[[i]]$Lon[j]<-y_exts[[i]][which(y_exts[[i]]$label==corrbg[[i]]$label[j]),1]
    corrbg[[i]]$Lat[j]<-y_exts[[i]][which(y_exts[[i]]$label==corrbg[[i]]$label[j]),2]
  }}

scorrbg1<-list()
linesc<-list()
X11(width=4,height=9)
par(mfrow=c(3,1),omi=c(0.2,0.1,0,0),mai=c(0.2,0.1,0.1,0.1))
for (i in 1:3){
  scorrbg1[[i]]<-corrbg1[[i]][,c(4,5,2)]
  coordinates(scorrbg1[[i]])<-c("Lon","Lat")
  gridded(scorrbg1[[i]])<-TRUE
  scorrbg1[[i]]<-raster(scorrbg1[[i]])
  scorrbg1[[i]]<-extend(scorrbg1[[i]],extent(-180, 180, -90, 90))
  scorrbg1[[i]]<-as(scorrbg1[[i]],"SpatialGridDataFrame")
  mapYcorGR<-mapGriddedData(scorrbg1[[i]],catMethod = c(-1,-0.8,-0.6,-0.4,-0.2,-0.00001,0.00001,0.2,0.4,0.6,0.8,1),
                            colourPalette = c("red4","red","hotpink1","plum1","pink","white","deepskyblue1","dodgerblue","steelblue3","blue","blue4"),
                            borderCol = "gray",oceanCol="azure2",xlim=c(-180,180),ylim=c(57,90),landCol="gray",addLegend=FALSE)
  if(i==3){do.call(addMapLegend,c(mapYcorGR, legendLabels="all", legendWidth=1.1,digits=1,legendShrink=0.7,legendMar=2))}
  title( main =paste0("Corr. Yield ",crop[i]), line=-0.5, cex.main=2,font.main=2)
}
scorrbg<-list()
X11(width=4,height=9)
par(mfrow=c(3,1),omi=c(0.0,0.0,0.0,0.0),mai=c(0.65,0.5,0.65,0.1))
for (i in 1:3){ 
  scorrbg[[i]]<-corrbg[[i]]
  temp<-na.omit(scorrbg[[i]])
  temp$CORy.sig<-as.factor(ifelse(temp$COR.p.value>0.05,0,1))
  scorrbg[[i]]<-temp
  linesc[[i]]<-join(scorrbg[[i]][,c(4,5,2,6)],Prod_Ensemdf1[[i]][,c(2,3,6)],type = "inner")
  plot(linesc[[i]][,5]*1000,linesc[[i]][,3],xlab="", ylab="", main="",col=c("gray","black")[linesc[[i]]$CORy.sig])
  title( main=paste0(crop[i]),line=2,font=2, font.lab=2,cex.lab=2,cex.main=2,cex.axis=1.5)
  title( ylab="Corr. Coeff (r)", xlab=if(i==3){"Production (Thousand Tons)"},line=2.5,cex.lab=1.5,font=2, font.lab=2)
}
