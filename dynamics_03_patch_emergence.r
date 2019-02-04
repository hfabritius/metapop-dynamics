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

polygonProjection = "+proj=utm +zone=35 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"

# 4. Patch emergence (length=probability equals destruction vector):
emergence<-function(data,x,emergRate,destrRate,alpha,im,em,AreaNuggetVariance,z,time){

  H<-data[[1]]; d<-data[[2]]; M<-data[[3]]; L<-data[[4]]; events<-data[[5]]; Landscape<-data[[6]]; RandomField<-data[[7]]; Area<-data[[8]]; Area_RandomField<-data[[9]]; D<-data[[10]]; newpatches<-data[[11]]
  created=0
  
  while(created==0){ # A loop to make sure that the emerging patch does not overlap with other patches

    # find the cell ID of the landscape raster that needs to be udpated
    y<-sample(x=1:ncell(L),size=1,prob=abs(emergRate*(destrRate*L[]/z)/(1-(L[]/z))))
    
    # Create a new spatial polygon for the new patch:
    newpoly<-Polygon(xyFromCell(L,y),hole=FALSE) # Polygon object located based on the raster cell ID (L[y])
    newpoly@area<-exp(log(Area[y])+rnorm(n=1,mean=0,sd=sqrt(AreaNuggetVariance)))/10^6
    newpoly<-Polygons(list(newpoly),paste("new",newpatches,sep="")) # Polygons object
    newpoly<-SpatialPolygons(Srl=list(newpoly),proj4string=CRS(polygonProjection)) # SpatialPolygons object
    newdata<-as.data.frame(H[1,]); newdata[1,]<-NA
    newdata$PatchArea<-sapply(slot(newpoly,"polygons"),slot,"area")
    newdata$PatchNr<-paste("new",newpatches,sep="")
    newdata$Occupied<-newdata$Protected<-newdata$Protectable<-0; newdata$d<-destrRate; newdata$New<-1
    newdata$Xcoord<-xyFromCell(L,y)[1]; newdata$Ycoord<-xyFromCell(L,y)[2]
    newdata$EmergYear<-time
    newdata$DestrYear<-0
    rownames(newdata)<-newdata$PatchNr
    newpoly<-SpatialPolygonsDataFrame(newpoly,newdata,match.ID=FALSE) # SpatialPolygonsDataFrame object
    rownames(newpoly@data)<-newdata$PatchNr

    # Evaluate whether the patch "fits" to the landscape
    distances<-rdist(coordinates(newpoly),coordinates(H))
    neighbours<-H[(which(rank(distances)<=5)),] # Find five closest patches
    
    overlap=0
    for(n in 1:nrow(neighbours)) if(sum(sqrt(neighbours$PatchArea[n]/pi),sqrt(newpoly$PatchArea/pi))>rdist(coordinates(newpoly),coordinates(neighbours[n,]))/1000) overlap==1 #patch radii overlap
    if(overlap==0) created=1 # The patch fits and the loop can finish
  }

  H<-spRbind(H,newpoly)
  newpatches<-newpatches+1

  # update the connectivity table, new patch as the last item (=> gets a colonization value)
  d<-rdist(coordinates(H))/1000 # update centroid-to-centroid distances between patches
  colnames(d)<-rownames(d)<-H$PatchNr
  diag(d)<-0 # ??
  M<-exp(-alpha*d)*(H$PatchArea)^im # each row has one focal patch
  diag(M)<-0 # remove patch colonization rate effect to self
  M[,H$Occupied==0]<-0 # remove patch colonization rate effect to empty patches
  for(i in 1:nrow(d)) M[,i]<-M[,i]*H$PatchArea[i]^em # columns, remote patch
  if(sum(H$PatchArea[H$Occupied==0])>0){
    if(nrow(H)>1){ #if there is more than one patch,
      if(nrow(H[H$Occupied==0,])>1) H$ce[H$Occupied==0]<-colRate*rowSums(M[H$Occupied==0,]) # sum values to get colonization rates per patch
      if(nrow(H[H$Occupied==0,])==1) H$ce[H$Occupied==0]<-colRate*rowSums(t(matrix(M[H$Occupied==0,]))) # sum values to get colonization rates per patch
    }
    if(nrow(H)==1) H$ce<-0 #=> the only patch is empty and thus no patch can be colonized
  }
  
  H$ce[is.na(H$ce)]<-1

  H$Occupied[nrow(H)]<-1
  nearest<-conn<-rdist(coordinates(H))/1000 # centroid-to-centroid distances between patches
  colnames(nearest)<-colnames(conn)<-rownames(conn)<-H$PatchNr
  diag(nearest)<-diag(conn)<-NA # ??
  conn<-exp(-alpha*conn)
  conn<-t(conn*sqrt(H$PatchArea)) # each ROW (after transpose) now contains remote items to sum for connectivity
  for(i in 1:nrow(H)) conn[i,which(H$Occupied==0)]<-0
  for(i in 1:nrow(H)) nearest[i,which(H$Occupied==0)]<-NA
  connectivity<-sqrt(H$PatchArea)
  for(i in 1:length(connectivity)) connectivity[i]<-connectivity[i]*rowSums(conn,na.rm=TRUE)[i]
  H$Occupied[nrow(H)]<-0
  
  L[y]<-0 # identify the patch location in the landscape & set emergence rate to zero
  events<-rbind(events,c("Patch emergence",xyFromCell(L,y),as.character(newdata$PatchNr),time,connectivity[length(connectivity)],min(nearest[length(connectivity),!is.na(nearest[length(connectivity),])]))) #store events: patch creation at a given location
  data<-list(H,d,M,L,events,Landscape,RandomField,Area,Area_RandomField,D,newpatches)
  
  return(data)
}
