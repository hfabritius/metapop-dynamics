# *******************************************************************************************************
# Metapopulation simulation with a changing landscape. Henna Fabritius, 12.3.2015
# *******************************************************************************************************

setwd(choose.dir())
polygonProjection = "+proj=utm +zone=35 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"

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
library(INLA)
library(spdep)
library(akima)  # interp
library(rgl)    # rgl
library(lattice)
library(grid) 

source("duplicate_patches.r")
source("patch_colonization_rates.r")
source("dynamics_01_patch_colonization.r")
source("dynamics_02_patch_extinction.r")
source("dynamics_03_patch_emergence.r")
source("dynamics_04_patch_destruction.r")
source("dynamics_05_landscape_change_annual.r")

# DEFINE PARAMETERS
# **********************************************************************************

# simulation parameters
repeatnumber<-100 # how many times each analysis is repeated
DestrRates  <-c(0.015, 0.03, 0.045, 0.06, 0.075, 0.09, 0.105, 0.12, 0.135, 0.15) # 0.015 steps, roughly equivalent to the range of 0.033 to 0.163 field rates
phis<-c(1,0.98,0.94,0.62) # tested values of phi (how fast the landscape changes), smallest phi means fastest change
extent<-extent(199995, 236015, 6874500, 6930005) # spatial extent of analysis
emergRate<-1 # patch emergence rate is not scaled in these simulations

# species-specifie parameters
alpha<-0.7 # dispersal parameter for M. diamina
ex<-0.17   # patch area vs. extinction scaling parameter for M.diamina based on M.cinxia
im<-0.30   # patch area vs. immigration scaling parameter for M.diamina based on M.cinxia
em<-0.07   # patch area vs. emigration scaling parameter for M.diamina based on M.cinxia
extRate<-0.1542996  # extinction rate scaling factor (Ovaskainen & Hanski 2003, estimated with an ABC algorithm) * FIXED MEAN 10.3.2016
colRate<-0.7613443 # colonization rate scaling factor, estimated with an ABC algorithm based on field data * FIXED MEAN 10.3.2016
z<-0.3866092 # mean patch area vs. landscape raster resolution

# LOAD & PREPARE BACKGROUND DATA
# **********************************************************************************

# 1. Prediction layer of patch locations (= occurrence probability)
LinearPrediction     <-raster("LinearPrediction_200x200.tif")
RandomField_original <-raster("RandomField_200x200.tif")
projection(LinearPrediction)    <-polygonProjection
projection(RandomField_original)<-polygonProjection
LinearPrediction    <-crop(LinearPrediction,extent)
RandomField_original<-crop(RandomField_original,extent)
LinearPrediction[]<-logit(LinearPrediction[])
RandomField_original[]<-logit(RandomField_original[])
RFmean<-mean(RandomField_original[!is.na(RandomField_original)])
RFvar<-var(RandomField_original[!is.na(RandomField_original)])

# 2. Prediction layer of patch areas
Area_LinearPrediction     <-raster("AREA_LinearPrediction_200x200.tif")
Area_RandomField_original <-raster("AREA_RandomField_200x200.tif")
projection(Area_LinearPrediction)    <-polygonProjection
projection(Area_RandomField_original)<-polygonProjection
Area_LinearPrediction    <-crop(Area_LinearPrediction,extent)
Area_RandomField_original<-crop(Area_RandomField_original,extent)
Area_LinearPrediction[]  <-log(Area_LinearPrediction[])
AreaNuggetVariance<-0.61821
AreaRFmean<-mean(Area_RandomField_original[!is.na(Area_RandomField_original)])
AreaRFvar<-var(Area_RandomField_original[!is.na(Area_RandomField_original)])

# 3. The habitat patch network
Habitat <- readShapeSpatial("habitat_network.shp", proj4string=CRS(polygonProjection), IDvar="PatchNr")
Habitatdata = read.csv2("Habitat.csv",header=TRUE)
Habitat@data<-Habitatdata
Habitat<-Habitat[((coordinates(Habitat)[,1]>xmin(extent))&(coordinates(Habitat)[,1]<xmax(extent))&(coordinates(Habitat)[,2]>ymin(extent))&(coordinates(Habitat)[,2]<ymax(extent))),] # Remove patches that are not covered by the raster
Habitat$PatchArea<-Habitat$PatchArea/10^6
rownames(Habitat@data)<-Habitat$PatchNr
Habitat$Protected<-0
Habitat$PatchNr<-as.character(Habitat$PatchNr)

# Split partially protectable sites into two, one of which is protectable (although in this analysis nothing will be protected, this is for comparability with the protection analysis)
partial<-(Habitat$Protectable>0 & Habitat$Protectable<1)
Duplicates<-Habitat[partial,]
Habitat<-Habitat[!partial,]
rownames(Duplicates@data)<-Duplicates$PatchNr
Habitat<-duplicates(Duplicates,Habitat)
Habitat$EmergYear<-Habitat$DestrYear<-0

# Create an spde object
twi = raster("twi11.tif") # TWI >> USE AS BASE RASTER (resolution,extent etc)
data<-read.csv2("randomdata_137000.csv")[,2:9]
locations <- as.matrix(data[,c("x","y")])
locations <- locations * 1e-6
mesh <- inla.mesh.2d(locations, max.edge=c(0.005,2), cutoff=0.001) # mesh <- inla.mesh.create.helper(points.domain=locations, max.edge=c(50,10000))  
spde <- inla.spde2.matern(mesh, alpha=2) # precision matrix # A <- inla.spde.make.A(mesh,locations) # observation matrix #image(spde$param.inla$M2), #image(spde$param.inla$M0) #image(A)

# REPEAT SIMULATION AS MANY TIMES AS POSSIBLE
# **********************************************************************************
for(repeats in 901:1000){
Results<-data.frame()
  
for(destr in 1:10){ # Three scenarios of analysis, each repeated for XX times
destrRate<-DestrRates[destr]

for(landscRate in 1:length(phis)){ # Three scenarios of analysis, each repeated for XX times
phi<-phis[landscRate]

# INITIALIZE DATA FOR SIMULATION
# **********************************************************************************

# Create working versions for analyses
Landscape<-L<-LinearPrediction
Landscape[]<-L[]<-inv.logit(LinearPrediction[] + RandomField_original[]) # create the original landscape
Landscape[is.na(Landscape)]<-0
L[is.na(L)]<-0
L[cellFromXY(L, coordinates(Habitat))]<-0
RandomField<-RandomField_original

Area<-Area_LinearPrediction
Area[]<-exp(Area_LinearPrediction[] + Area_RandomField_original[])
Area[is.na(Area)]<-0
Area_RandomField<-Area_RandomField_original

H<-Habitat # create an iterative working version
rownames(H@data)<-H$PatchNr

# Estimate extinction rates for occupied patches & patch destruction rates
H$ce<-NA
H$ce[(H$Occupied)==1]<-extRate/(H$PatchArea[H$Occupied==1]^ex) # Ovaskainen & Hanski 2003
H$d[H$Protected==0]<-destrRate # unconserved patches may disappear from the system
H$d[H$Protected>0]<-0 # conserved patches do not get destructed
D<-H[0,] # an empty placeholder for deleted fields

# Estimate colonization rates for unoccupied patches (Hanski & Ovaskainen 2000: The metapopulation capacity of a fragmented landscape)
source("patch_colonization_rates.r")
colonizationrates<-colonizationrates(H,alpha,im,em)
d<-colonizationrates[[1]]
M<-colonizationrates[[2]]
H$ce[H$Occupied==0]<-colonizationrates[[3]]

# set patch emergence probability to zero at sites that already have patches
if(sum(H$PatchArea)>0) L[cellFromXY(L, coordinates(H))]<-0 # identify patch locations in the landscape

# Set up a time counter t to run from 0 to 50 years
t<-0; x<-0; events<-t(matrix(c(0,0,0,0,0,0))) # keep track of total time & year & what happens in the network & where
colnames(events)<-c("Event","x coordinate","y coordinate","PatchName","Time","Connectivity")
newpatches<-1 # set up a counter for setting names/IDs to new patches
data<-list(H,d,M,L,events,Landscape,RandomField,Area,Area_RandomField,D,newpatches)

# RUN SIMULATION FOR 50 YEARS & STORE THE RESULTS
# **********************************************************************************
while(sum(t)<50){

# continuous-time simulation
data[[4]][data[[4]]>z]<-z
t<-c(t,rexp(1,     sum(data[[1]]$ce[data[[1]]$Occupied==0],data[[1]]$ce[data[[1]]$Occupied==1],data[[1]]$d,abs(emergRate*(destrRate*data[[4]][]/z)/(1-(emergRate*data[[4]][]/z)))))) # choose a time until something happens & store to the time vector
x<-c(x,runif(1,max=sum(data[[1]]$ce[data[[1]]$Occupied==0],data[[1]]$ce[data[[1]]$Occupied==1],data[[1]]$d,abs(emergRate*(destrRate*data[[4]][]/z)/(1-(emergRate*data[[4]][]/z)))))) # choose what happens & store info

if(floor(sum(t))>floor(sum(t[1:length(t)-1])))  data<-landscapechange(data,LinearPrediction,Area_LinearPrediction,phi,spde,twi,sum(t),RFmean,AreaRFmean,RFvar,AreaRFvar) # annual change to the landscape
# add hear if you want to record yearly development of the network

if(x[length(x)]<=sum(data[[1]]$ce[data[[1]]$Occupied==0]))                                   event="Patch colonization"
if(x[length(x)] >sum(data[[1]]$ce[data[[1]]$Occupied==0]) & x[length(x)]<=sum(data[[1]]$ce)) event="Patch extinction"
if(x[length(x)] >sum(data[[1]]$ce) & x[length(x)]<=sum(data[[1]]$ce,data[[1]]$d))            event="Patch destruction"
if(x[length(x)] >sum(data[[1]]$ce,data[[1]]$d))                                              event="Patch emergence"

if(event=="Patch colonization") data<-colonization(data,x,alpha,im,em,ex,extRate,colRate,sum(t))
if(event=="Patch extinction")   data<-extinction(data,x,colRate,sum(t))
if(event=="Patch destruction")  data<-destruction(data,x,sum(t))
if(event=="Patch emergence")    data<-emergence(data,x,emergRate,destrRate,alpha,im,em,AreaNuggetVariance,z,sum(t))

# interrupt the simulation if the patch network disappears or population goes extinct
if(sum(data[[1]]$PatchArea)==0) break     
if(sum(data[[1]]$PatchArea[data[[1]]$Occupied==1])==0) break
} # 50 years

H<-data[[1]]; d<-data[[2]]; M<-data[[3]]; L<-data[[4]]; events<-data[[5]]; Landscape<-data[[6]]; RandomField<-data[[7]]; Area<-data[[8]]; Area_RandomField<-data[[9]]; D<-data[[10]]; newpatches<-data[[11]]

d<-rdist(coordinates(H))/1000 # centroid-to-centroid distances between patches
colnames(d)<-rownames(d)<-H$PatchNr
diag(d)<-NA # ??
d<-exp(-alpha*d)
d<-t(d*sqrt(H$PatchArea)) # each ROW (after transpose) now contains remote items to sum for connectivity
H$Connectivity<-sqrt(H$PatchArea)
for(i in 1:nrow(H)) H$Connectivity[i]<-H$Connectivity[i]*rowSums(d,na.rm=TRUE)[i]

# Store the results of an individual run with DM_results_landscapechange_repeat (create a file before first repeat!!)
if(sum(H$PatchArea)>0){
if(sum(H$PatchArea[H$Occupied==1])>0)  result<-c(destr,phi,100*nrow(H[H$Occupied==1,])/nrow(H),nrow(H),sum(t),median(rdist(coordinates(H))/1000),mean(H$Connectivity),mean(H$PatchArea),mean(Landscape[!is.na(Landscape)])) # destruction rate, phi, % occupied, Nr. of patches, year, median distance bw/ patches, mean connectivity, mean patch area
if(sum(H$PatchArea[H$Occupied==1])==0) result<-c(destr,phi,0,nrow(H),sum(t),median(rdist(coordinates(H))/1000),mean(H$Connectivity),mean(H$PatchArea),mean(Landscape[!is.na(Landscape)]))
}
if(sum(H$PatchArea)==0) result<-c(destr,phi,0,0,0,sum(t),0,0,0,mean(Landscape[!is.na(Landscape)]))
Results<-rbind(Results,result)
names(result)<-colnames(Results)<-c("Destruction Rate","Phi","% occupied","Nr. of patches","Year", "Median interpatch distance","Connectivity","Mean patch area","Landscape suitability")

PatchNetwork<-rbind(H@data[,c("PatchArea","PatchNr","EmergYear","DestrYear","Occupied")],D@data[,c("PatchArea","PatchNr","EmergYear","DestrYear","Occupied")])
write.csv2(PatchNetwork,paste("Patch_Network_Rep",repeats,"_Phi",phi,"_Destr",destrRate,".csv",sep=""))
write.csv2(events,paste("Events_Rep",repeats,"_Phi",phi,"_Destr",destrRate,".csv",sep=""))
} # landscRate (phis)
} # destr
write.csv2(Results,paste("SimulationResults_Rep",repeats,".csv",sep=""))    
} # repeats
