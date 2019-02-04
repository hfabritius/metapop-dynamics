library(sp)
library(maptools) 
library(rgeos) 
library(fields)
library(graphics)
library(ncdf)
library(utils)
library(INLA)
library(spdep)
library(akima)  # interp
library(raster) # raster
library(rgdal)  # write .tif raster
library(boot)   # inv.logit
library(rgl)    # rgl
library(lattice)
library(grid) 

landscapechange<-function(data,LinearPrediction,Area_LinearPrediction,phi,spde,twi,time,RFmean,AreaRFmean,RFvar,AreaRFvar){
  # Landscape
  # L (patch centroids = 0)
  # New Random Field
  # Area
  # Area_RandomField
  
  if(phi<1){
  H<-data[[1]]; d<-data[[2]]; M<-data[[3]]; L<-data[[4]]; events<-data[[5]]; Landscape<-data[[6]]; RandomField<-data[[7]]; Area<-data[[8]]; Area_RandomField<-data[[9]]; D<-data[[10]]; newpatches<-data[[11]]

  # HABITAT SUITABILITY MODEL
  # ---------------------------------------------------------------------------
  theta1<-(-7.8991) # locations model
  theta2<-5.979
  theta1D<-theta1-log(1-phi^2)/2
  
  samplemap<-Epsilon<-Landscape
  Q <- inla.spde2.precision(spde, theta=c(theta1D,theta2))
  
  sample<-inla.qsample(n=1, Q=Q) # inla.qsample(n = 1L, Q=Q, b, mu, sample, constr, reordering = inla.reorderings(), seed = 0L, logdens = ifelse(missing(sample), FALSE, TRUE))
  samplemap <-interp(x=(mesh$loc[,1]*1000000),y=(mesh$loc[,2]*1000000),z=sample[,1] ,xo=seq(xmin(twi),xmax(twi),length=ncol(twi)), yo=seq(ymin(twi),ymax(twi),length=nrow(twi)), linear=FALSE, extrap=TRUE)
  samplemap <-flip(raster(t(samplemap$z)), direction='y') # samplemap is still a matrix
  Epsilon<-twi
  Epsilon[]<-samplemap[] # make it a raster
  Epsilon<-aggregate(Epsilon,fact=20)
  Epsilon<-crop(Epsilon,extent(RandomField))
  
  RandomField[]<-(phi*RandomField[] + Epsilon[])
  RandomField[]<-RandomField[]/sqrt(var(RandomField[!is.na(RandomField)])/RFvar) # adjust variance to that of the original random field
  Landscape[]<-inv.logit(LinearPrediction[] + RFmean +  RandomField[]) # shift the mean as much as that of the original random field
  Landscape[is.na(Landscape)] <- 0 # sea and lakes: set emergence rate to zero  
  L<-Landscape
  if(sum(H$PatchArea)>0) L[cellFromXY(L, coordinates(H))]<-0 # identify patch locations in the landscape

  # AREA MODEL
  # ---------------------------------------------------------------------------
  theta1<-(-7.595) # area model
  theta2<-6.914
  theta1D<-theta1-log(1-phi^2)/2
  
  samplemap<-Epsilon<-Landscape
  Q <- inla.spde2.precision(spde, theta=c(theta1D,theta2))
  
  sample<-inla.qsample(n=1, Q=Q) # inla.qsample(n = 1L, Q=Q, b, mu, sample, constr, reordering = inla.reorderings(), seed = 0L, logdens = ifelse(missing(sample), FALSE, TRUE))
  samplemap <-interp(x=(mesh$loc[,1]*1000000),y=(mesh$loc[,2]*1000000),z=sample[,1] ,xo=seq(xmin(twi),xmax(twi),length=ncol(twi)), yo=seq(ymin(twi),ymax(twi),length=nrow(twi)), linear=FALSE, extrap=TRUE)
  samplemap <-flip(raster(t(samplemap$z)), direction='y') # samplemap is still a matrix
  Epsilon<-twi
  Epsilon[]<-samplemap[] # make it a raster
  Epsilon<-aggregate(Epsilon,fact=20)
  Epsilon<-crop(Epsilon,extent(Area_RandomField))
  
  Area_RandomField[]<-(phi*Area_RandomField[] + Epsilon[])
  RandomField[]<-RandomField[]/sqrt(var(RandomField[!is.na(RandomField)])/RFvar) # adjust variance to that of the original random field
  Area[]<-exp(Area_LinearPrediction[] + AreaRFmean + Area_RandomField[]) # shift the mean as much as that of the original random field
  Area[is.na(Area)] <- 0 # sea and lakes: set emergence rate to zero  

  events<-rbind(events,c("Landscape change",mean(RandomField[!is.na(RandomField)]),mean(Landscape[!is.na(Landscape)]),0,time,NA,NA)) #store events: patch colonization at a given location
  data<-list(H,d,M,L,events,Landscape,RandomField,Area,Area_RandomField,D,newpatches)
  }
  
  return(data)
}
