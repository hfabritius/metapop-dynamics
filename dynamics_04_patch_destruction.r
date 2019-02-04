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

destruction<-function(data,x,time){
  
  H<-data[[1]]; d<-data[[2]]; M<-data[[3]]; L<-data[[4]]; events<-data[[5]]; Landscape<-data[[6]]; RandomField<-data[[7]]; Area<-data[[8]]; Area_RandomField<-data[[9]]; D<-data[[10]]; newpatches<-data[[11]]
  y<-which(sum(H$ce)+cumsum(H$d)>x[length(x)])[1] # index number of the habitat patch that gets destructed
  
  L[cellFromXY(L, coordinates(H[y,]))]<-Landscape[cellFromXY(Landscape, coordinates(H[y,]))] # return the emergence rate of the patch
    
  nearest<-conn<-rdist(coordinates(H))/1000 # centroid-to-centroid distances between patches
  colnames(nearest)<-colnames(conn)<-rownames(conn)<-H$PatchNr
  diag(nearest)<-diag(conn)<-NA # ??
  conn<-exp(-alpha*conn)
  conn<-t(conn*sqrt(H$PatchArea)) # each ROW (after transpose) now contains remote items to sum for connectivity
  if(H$Occupied[y]==1) for(i in 1:nrow(H)) conn[i,which(H$Occupied==0)]<-0
  if(H$Occupied[y]==0) {
    H$Occupied[y]==1
    for(i in 1:nrow(H)) conn[i,which(H$Occupied==0)]<-0
    H$Occupied[y]==0
  }
  for(i in 1:nrow(H)) nearest[i,which(H$Occupied==0)]<-NA
  connectivity<-sqrt(H$PatchArea)
  for(i in 1:length(connectivity)) connectivity[i]<-connectivity[i]*rowSums(conn,na.rm=TRUE)[i]
    
  events<-rbind(events,c("Patch destruction",floor(coordinates(H[y,])),as.character(H$PatchNr[y]),time,connectivity[y],min(nearest[y,!is.na(nearest[y,])]))) # store events: patch destruction at a given location
  H$DestrYear[y]<-time
  D<-spRbind(D,H[y,]) # move the deleted patch to the list of deleted patches
  H<-H[-y,] # remove patch from the patch network
  if(class(M)=="matrix") M<-M[,-y] # remove patch from connectivity table & remove connectivity effect to other patches
  if(class(M)=="matrix") M<-M[-y,] # remove patch from connectivity table & remove connectivity effect to other patches
  if(class(M)=="numeric") M<-M[-y] # remove patch from connectivity table & remove connectivity effect to other patches
  if(class(d)=="matrix") d<-d[,-y] # remove patch from distance to other patches
  if(class(d)=="matrix") d<-d[-y,] # remove patch from distance to other patches
  if(class(d)=="numeric") d<-d[-y] # remove patch from distance to other patches
  if(sum(H$PatchArea[H$Occupied==0])>0){
    if(nrow(H)>1){ #if there is more than one patch,
      if(nrow(H[H$Occupied==0,])>1) H$ce[H$Occupied==0]<-rowSums(M[H$Occupied==0,]) # update connectivity values after patch removal
      if(nrow(H[H$Occupied==0,])==1) H$ce[H$Occupied==0]<-rowSums(t(matrix(M[H$Occupied==0,]))) # update connectivity values after patch removal
    }
    if(nrow(H)==1) H$ce<-0 #=> the only patch is empty and thus no patch can be colonized
  }
  
  H$ce[is.na(H$ce)]<-1
  
  data<-list(H,d,M,L,events,Landscape,RandomField,Area,Area_RandomField,D,newpatches)
  
  return(data)
}
