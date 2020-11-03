library(raster)
library(ncdf4)
library(plyr)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(forecast)

crop<-c("Wheat","Maize","Rice")
yield_LPJR<-list()
ensemYi<-list()
for (i in 1:3){
  yield_LPJR[[i]]<-stack(paste0("C:/Users/hac809/Documents/Pughtam-cropozone/Global_evaluation_outputs/LPJYield_",crop[i],"_newcrops_1970-2010.nc")) 
  yield_LPJR[[i]]<-yield_LPJR[[i]][[11:41]]
  ensemYi[[i]]<-stack(paste0("C:/Users/hac809/Documents/Pughtam-cropozone/Best_Ensemble/Best_ensemble_outputs/",crop[i],"best_ensemb_1980-2010.nc"))
  }


Deraster<-function(x){
  df<-list()
  for(i in 1:nlayers(x)){
    df[[i]]<-as.data.frame(x[[i]],xy=TRUE)
    df[[i]]$year<-1979+i
  }
  df<-lapply(df, setNames, c("Lon","Lat","Yield","Years"))
  return(bind_rows(df, .id = "label"))
}

#memory.size()
#memory.limit(size=7000)
yield_LPJdf<-lapply(yield_LPJR,FUN=Deraster)
yield_BEdf<-lapply(ensemYi,FUN=Deraster)



cut<-function(x){
  cu<-x[c(-1,-2,-(length(x$year)),-(length(x$year)-1)),]
  return(data.frame(cu))
}
y_ext<-list()
yield_bothdf<-list()
yearsg_ma<-list()
for (i in 1:3){
  y_ext[[i]]<-cbind(yield_LPJdf[[i]][,-1],BE=yield_BEdf[[i]][,4])
  y_ext[[i]]$label<-rep(1:259200,31 )
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
  dtcropbg[[j]]<-cbind(year=yearsg_ma[[j]][,2],do.call("cbind",dtbg)[c(-3)])
  colnames(dtcropbg[[j]])<-c("year","label","LPJ","B_ens")
}
names(dtcropbg)<-crop

#CORRELATION
corrfun<-function(x){
  COR=cor.test(x$LPJ,x$B_ens)
  return(data.frame(COR$estimate,COR$p.value))
}
corrbg<-list()
corrbg1<-list()
Cropc<-list()
y_exts<-list()
for ( i in 1:3){
  corrbg[[i]]<-ddply(dtcropbg[[i]], .(label),corrfun)
  corrbg1[[i]]<-corrbg[[i]]
  corrbg1[[i]]$COR.estimate<-ifelse(corrbg1[[i]]$COR.p.value>0.05,0,corrbg1[[i]]$COR.estimate)
  y_exts[[i]]<-subset(y_ext[[i]],Years==1980)
  for(j in 1:nrow(corrbg1[[i]])){
    corrbg1[[i]]$Lon[j]<-y_exts[[i]][which(y_exts[[i]]$label==corrbg1[[i]]$label[j]),1]
    corrbg1[[i]]$Lat[j]<-y_exts[[i]][which(y_exts[[i]]$label==corrbg1[[i]]$label[j]),2]
  }
  
  for ( i in 1:3){
  Crop_corr<-ggplot(corrbg1[[i]],aes(xmin=Lon,ymin=Lat,xmax=Lon+0.5, ymax=Lat+0.5, fill=COR.estimate))
    Cropc[[i]]<-Crop_corr+geom_rect()+
    guides(fill=guide_colourbar(label=TRUE))+
    ggtitle(crop[i])+
    borders('world', colour='black',size = 0.00000001)+
    xlab("") +
    ylab("") +
    scale_fill_distiller(palette="RdBu", direction=1, limits=c(-1,1),values=NULL,na.value="gray50")+
    theme(plot.title = element_text(size = 25, face = "bold", hjust=0.5),
          legend.title = element_blank(),
          legend.key.width = unit(3,"cm"),
          legend.key.height = unit(1,"cm"),
          legend.position = "bottom",
          legend.direction = "horizontal",
          axis.text = element_text(size = 15),
          axis.title.x = element_text(size = 40, vjust = -0.5),
          axis.title.y = element_text(size = 40, vjust = 0.2),
          legend.text = element_text(size = 20,vjust = -0.3),panel.background = element_rect(fill = "#BFD5E3"))
}
Cropc[[3]]
#REGRESSION
regfun<-function(x){
  GRY=lm(Yield~yearn*Dummy,data=x)
  return(data.frame(Slope=GRY$coefficients[4],pvi=summary(GRY)$coefficients[4, 4]))
}
yield_reg<-list()
regGRY<-list()
regGRY1<-list()
Cropr<-list()
for (i in 1:3){
  yieldsdf_<-na.omit(y_ext[[i]])
  yieldsdf_<-yieldsdf_[yieldsdf_$label %in% names(which(table(yieldsdf_$label)>3)), ]
  yieldsdf_<-yieldsdf_[,c(6,4,1,2,3,5)]
  yieldsdf_$yearn<-yieldsdf_$Years-1969
  yieldsdf_L<-yieldsdf_[,c(1,7,3,4,5)]#LPJ data
  yieldsdf_B<-yieldsdf_[,c(1,7,3,4,6)]#BE data
  colnames(yieldsdf_B)<-names(yieldsdf_L)
  yieldsdf_B$Dummy<-1
  yieldsdf_L$Dummy<-0
  yield_reg[[i]]<-rbind(yieldsdf_B,yieldsdf_L)
  regGRY[[i]]<-ddply(yield_reg[[i]], .(label), regfun)
  regGRY1[[i]]<-regGRY[[i]]
  regGRY1[[i]]$Slope<-ifelse(regGRY1[[i]]$pvi>0.05,0,regGRY1[[i]]$Slope)
  for(j in 1:nrow(regGRY1[[i]])){
    regGRY1[[i]]$Lon[j]<-y_exts[[i]][which(y_exts[[i]]$label==regGRY1[[i]]$label[j]),1]
    regGRY1[[i]]$Lat[j]<-y_exts[[i]][which(y_exts[[i]]$label==regGRY1[[i]]$label[j]),2]
  }}
  for (i in 1:3){
    Crop_reg<-ggplot(regGRY1[[i]],aes(xmin=Lon,ymin=Lat,xmax=Lon+0.5, ymax=Lat+0.5, fill=Slope))
    Cropr[[i]]<-Crop_reg+geom_rect()+
      guides(fill=guide_colourbar(label=TRUE))+
      ggtitle(crop[i])+
      borders('world', colour='black',size = 0.00000001)+
      xlab("") +
      ylab("") +
      scale_fill_distiller(palette="RdBu", direction=1, limits=c(-0.25,0.25),values=NULL,na.value="gray50")+
      theme(plot.title = element_text(size = 25, face = "bold", hjust=0.5),
            legend.title = element_blank(),
            legend.key.width = unit(3,"cm"),
            legend.key.height = unit(1,"cm"),
            legend.position = "bottom",
            legend.direction = "horizontal",
            axis.text = element_text(size = 15),
            axis.title.x = element_text(size = 40, vjust = -0.5),
            axis.title.y = element_text(size = 40, vjust = 0.2),
            legend.text = element_text(size = 20,vjust = -0.3),panel.background = element_rect(fill = "#e3ddbf"))
  }
  Cropr[[3]]
  
  