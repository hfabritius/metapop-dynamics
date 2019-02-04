# Henna Fabritius, 2.12.2014

library(sp)
library(maptools) 
library(rgdal)
library(rgeos) 
library(fields)
library(raster)
library(graphics)
library(ncdf)
library(utils)
library(boot)

colonization<-function(data,x,alpha,im,em,ex,extRate,colRate,time){


  H<-data[[1]]; d<-data[[2]]; M<-data[[3]]; L<-data[[4]]; events<-data[[5]]; Landscape<-data[[6]]; RandomField<-data[[7]]; Area<-data[[8]]; Area_RandomField<-data[[9]]; D<-data[[10]]; newpatches<-data[[11]]
  y<-which(cumsum(H$ce[H$Occupied==0])>x[length(x)])[1] # index number of the habitat patch that gets colonized
  y<-which(cumsum(H$Occupied==0)==y)[1] # correction for the whole list of habitats
  
  H$Occupied[y]<-1 # update occupancy status to occupied
  H$ce[y]<-extRate/(H$PatchArea[y]^ex) # update the colonization rate to extinction rate
  if(class(M)=="matrix") M[,y]<-exp(-alpha*d[,y])*((H$PatchArea)^im)*(H$PatchArea[y]^em) # add the colonized patch to other patches' connectivity measures
  if(class(M)=="numeric") M[y]<-exp(-alpha*d[y])*((H$PatchArea)^im)*(H$PatchArea[y]^em) # add the colonized patch to other patches' connectivity measures
  if(class(M)=="matrix") M[y,y]<-0 # remove patch connectivity to self
  if(class(M)=="numeric") M[y]<-0 # remove patch connectivity to self

  if(sum(H$PatchArea[H$Occupied==0])>0){
    if(nrow(H)>1){ #if there is more than one patch,
      if(nrow(H[H$Occupied==0,])>1) H$ce[H$Occupied==0]<-colRate*rowSums(M[H$Occupied==0,]) # sum values to get connectivity measures per empty patch
      if(nrow(H[H$Occupied==0,])==1) H$ce[H$Occupied==0]<-colRate*rowSums(t(matrix(M[H$Occupied==0,]))) # sum values to get connectivity measures per empty patch
    }
   if(nrow(H)==1) H$ce<-0 #=> the only patch is empty and thus no patch can be colonized
  }
  
  H$ce[is.na(H$ce)]<-1
  
  nearest<-conn<-rdist(coordinates(H))/1000 # centroid-to-centroid distances between patches
  colnames(nearest)<-colnames(conn)<-rownames(conn)<-H$PatchNr
  diag(nearest)<-diag(conn)<-NA # ??
  conn<-exp(-alpha*conn)
  conn<-t(conn*sqrt(H$PatchArea)) # each ROW (after transpose) now contains remote items to sum for connectivity
  for(i in 1:nrow(H)) conn[i,which(H$Occupied==0)]<-0 # for calculating connectivity measures
  for(i in 1:nrow(H)) nearest[i,which(H$Occupied==0)]<-NA # for estimating distance to the nearest occupied patch
  connectivity<-sqrt(H$PatchArea)
  for(i in 1:length(connectivity)) connectivity[i]<-connectivity[i]*rowSums(conn,na.rm=TRUE)[i]
    
  events<-rbind(events,c("Patch colonization",floor(coordinates(H[y,])),as.character(H$PatchNr[y]),time,connectivity[y],min(nearest[y,!is.na(nearest[y,])]))) #store events: patch colonization at a given location
  data<-list(H,d,M,L,events,Landscape,RandomField,Area,Area_RandomField,D,newpatches)
  
  return(data)
}
