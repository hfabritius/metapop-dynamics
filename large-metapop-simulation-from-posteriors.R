# ****************************************************************************************************************
# Code for running a large-scale metapopulation simulation from posterior distributions -Henna Fabritius
# Lots of code optimization & data table handling
# ****************************************************************************************************************
library(coda) # read.openbugs
library(fields) # rdist
library(VGAM) # cloglog
library(LaplacesDemon) # rbern, Bernoulli distribution
library(boot) # inv.logit
library(data.table) # data.table, fread (data tables faster to process than data frames)

# Study area data: a distance matrix and a data table of stand IDs and nr. of colonization sources per stand for aggregated connectivity estimation (to be filled during the simulation)
scenario<-"scenario_name"
treelist_name<-"trees_of_the_simulation"
subset_area<-"subarea_3000.txt" # for restricting the simulation size

stand_distances<-paste("stand_distances_RetPat",scenario,"_3000.txt",sep="") # read file name only (the file is scanned row by row later) # stand_distances<-file(description=stand_distance_file,open="r"), close(stand_distances)
stand_IDs<-fread(paste("stand_attributes_RetPat",scenario,"_3000.txt",sep=""))
colnames(stand_IDs)<-c("ID","x","y")
stand_IDs<-stand_IDs$ID # This includes all stand IDs in the shapefile, also the retention patches (with "_NC")

# Scenario folder that includes the treelist & possible t0 occurrences & will be used for storing related results
setwd(paste("E:/New_Indata/",scenario,sep="")) # used for storing the results
trees<-read.csv2(paste(treelist_name,".csv",sep=""),dec=",",header=T)
t0_occurrences<-list.files(pattern ="occurrences",recursive=TRUE) # length == 0 if t0 occurrences haven't been saved for use in repeated analyses
setDT(trees) # Set class into data table
trees<-trees[(Year %in% c(0,10,20,30,40,50,60,70,80,90,100)),]

# Include only trees that are part of the subset (in case we want to use a smaller subset)
subset<-read.table(subset_area,header=T)$x # alternative subset sizes
trees<-trees[(gsub("_NC","",trees[,StandID]) %in% subset),] # in case the subarea is smaller than the tree list

# Cut the stand attributes & distances into a subset if needed
incl_stands<-which(stand_IDs %in% trees[,StandID]) # For chopping the distance matrix rows to the right length
stand_IDs<-stand_IDs[stand_IDs %in% trees[,StandID]] # Reduce stands to only those that are part of the analysis (source trees cannot develop elsewhere). NOTE: This would require saving the stand ID incides for cutting also the distance matrix
NrSources<-rep(0,times=length(stand_IDs)) # Nr of stands with potential colonization sources (used for aggregated connectivity). Vectorized for faster computing

# Paths for reading in the study species model files (occupancy model for initializing the landscape)
model<-file.path("col_ext_model_")
occupancymodel<-file.path("occupancy_model_")

# Simulation parameters
TimeStepLength = 10 # Years between forest data time periods
NoRepl = 100 # How many times the stochastic analysis is repeated
Maturity<-10 # Age at which the modelled sp. reaches maturity and becomes a colonization source
MaturSize<-100 # cm2
ExactConnectivity<-2000 # Distance (m) within which connectivities are calculated exactly; for a larger radius, connectivities are calculated based on stand centroids & NrSources. 0 if not used.

# --------------------------------------------------------------------------------------------------------
# Read in study species model CODA file based on what was set above. NOTE: These may trigger a warning that is not fatal.
CODA_files<-read.openbugs(stem=occupancymodel,quiet=T)
OccModel<-rbind(as.matrix(CODA_files[1]),as.matrix(CODA_files[2]))
CODA_files<-read.openbugs(stem=model,quiet=T)
Model<-rbind(as.matrix(CODA_files[1]),as.matrix(CODA_files[2]))

# Estimate tree connectivity & add last year of the tree to all tree instances (to identify when trees are about to die)
trees[,c("connectivity","SpeciesOccupancy","deviation")]<-as.numeric() # Add columns here
trees[,"log_dbh"]<-(log(trees[,DBH])-3.248098)/0.383203 # predictive models have been built using normalized data 
trees[,"StandAgeC"]<-(trees[,StandAge]-92.29932)/70.29898 # predictive models have been built using normalized data
trees[StandAgeC>((200-92.29932)/70.29898),"StandAgeC"]<-(200-92.29932)/70.29898 # To prevent very old trees getting unreasonably high colonization rates
colonizable<-(30-92.29932)/70.29898

# --------------------------------------------------------------------------------------------------------
# METAPOPULATION SIMULATION
# --------------------------------------------------------------------------------------------------------
for(n in 50:NoRepl){ # Run multiple replicates of the analysis
i = sample(1:nrow(OccModel),1) # For each iteration, randomize the model parameter combinations used for the occupancy & col-ext models
ii = sample(1:nrow(Model),1) # For each iteration, randomize the model parameter combinations used for the occupancy & col-ext models
OccModel_i<-OccModel[i,] # OR: <-matrixStats::colMedians(OccModel). Take vectors of model parameters for faster computing
Model_ii<-Model[ii,] # Take vectors of model parameters for faster computing
trees[,c("connectivity","SpeciesOccupancy","deviation")]<-NA

# Define stand-specific deviations from the mean colonization rate
trees[,"deviation"][order(trees[,StandID])]<-rep(rnorm(n=length(unique(trees[,StandID])),mean=0,sd=(Model_ii["sigma.alpha"])),times=as.numeric(table(trees[,StandID]))[as.numeric(table(trees[,StandID]))!=0])

# 1. Model species occupancy at the first time point based on a non-spatial occupancy model
# --------------------------------------------------------------------------------------------------------
if(length(t0_occurrences)>0) trees[Year==min(Year),"SpeciesOccupancy"]<-fread(t0_occurrences)$SpeciesOccupancy # Read in t0 occurrences from a file
if(length(t0_occurrences)==0){ # If there are no t0 occurrences, define based on connectivity values

  # Vectorize data table columns for faster calculation
  nr<-nrow(trees[Year==min(Year)])
  SpeciesOccupancy<-rep(NA,times=nr)
  Species<-trees[Year==min(Year),Species]
  StandAgeC<-trees[Year==min(Year),StandAgeC]
  StandID<-trees[Year==min(Year),StandID]
  log_dbh<-trees[Year==min(Year),log_dbh]

  # Estimate occupancy based on occupancy model parameters
  if(("Ob0" %in% colnames(OccModel))==TRUE) SpeciesOccupancy<-rep(OccModel_i["Ob0"],times=nr)
  if(("Ob0[1]" %in% colnames(OccModel))==TRUE) SpeciesOccupancy<-OccModel_i["Ob0[1]"] * ifelse(Species=="PopTre",1,0) + OccModel_i["Ob0[2]"] * ifelse(Species=="SalCap",1,0)

  if(("b0" %in% colnames(OccModel))==TRUE) SpeciesOccupancy<-rep(OccModel_i["b0"],times=nr)
  if(("b1" %in% colnames(OccModel))==TRUE) SpeciesOccupancy<-SpeciesOccupancy + OccModel_i["b1"] * StandAgeC
  if(("b2" %in% colnames(OccModel))==TRUE) SpeciesOccupancy<-SpeciesOccupancy + OccModel_i["b2"] * StandAgeC^2
  if(("sigma.alpha" %in% colnames(OccModel))==TRUE) SpeciesOccupancy[order(StandID)]<-SpeciesOccupancy[order(StandID)]+rep(rnorm(n=length(unique(StandID)),mean=0,sd=(OccModel_i["sigma.alpha"])),times=as.numeric(table(StandID))[as.numeric(table(StandID))!=0])
  
  if(("beta_dbh" %in% colnames(OccModel))==TRUE) SpeciesOccupancy<-SpeciesOccupancy + OccModel_i["beta_dbh"] * log_dbh
  
  SpeciesOccupancy<-VGAM::cloglog(SpeciesOccupancy,inv=T)
  trees[Year==min(Year),"SpeciesOccupancy"]<-SpeciesOccupancy<-as.numeric(sapply(SpeciesOccupancy,function(x) rbern(1,x)))*Maturity # Randomize Lobaria occupancy based on the model prediction
  NrSources[stand_IDs %in% StandID][order(stand_IDs[stand_IDs %in% StandID])]<-aggregate(ifelse(SpeciesOccupancy>=Maturity,1,0),by=list(StandID),FUN=sum)$x[order(unique(StandID))] # Sum of colonization sources per included stand
  rm(SpeciesOccupancy,Species,StandAgeC,StandID,log_dbh); gc() # Clean memory
}

# 2. Simulate extinction-colonization based on the ext-col model
# --------------------------------------------------------------------------------------------------------
for(j in as.numeric(names(table(trees$Year)[2:length(table(trees$Year))]))){ if(nrow(trees[(Year==(j-TimeStepLength) & SpeciesOccupancy>0),])>0){ # For the remaining time points of the simulation, if the spp hasn't gone extinct,

  # The tree file could be read row by row to save memory: trees<-fread(stand_distances,skip=(j-1),nrows=1,header=T)
  # Add species occupancy from the past time step
  trees[Year==j & TreeID %in% trees[Year==(j-TimeStepLength) & SpeciesOccupancy>0,TreeID],"SpeciesOccupancy"][order(match(trees[Year==j & TreeID %in% trees[Year==(j-TimeStepLength) & SpeciesOccupancy>0,TreeID],SpeciesOccupancy],trees[Year==j & TreeID %in% trees[Year==(j-TimeStepLength) & SpeciesOccupancy>0,TreeID],TreeID]))]<-
  trees[Year==(j-TimeStepLength) & SpeciesOccupancy>0 & TreeID %in% trees[Year==j,TreeID],"SpeciesOccupancy"][order(match(trees[Year==(j-TimeStepLength) & SpeciesOccupancy>0 & TreeID %in% trees[trees$Year==j,TreeID],SpeciesOccupancy],trees[trees$Year==(j-TimeStepLength) & trees$SpeciesOccupancy>0 & TreeID %in% trees[trees$Year==j,TreeID],TreeID]))]+TimeStepLength # Occupancy of trees that were occupied at the previous time step & their IDs occur also this year
  trees[Year==j & is.na(SpeciesOccupancy),"SpeciesOccupancy"]<-0 # Set trees unoccupied if they were not occupied at the previous decade

  # Vectorize data table columns for faster calculation
  nr<-nrow(trees[Year==j])
  SpeciesOccupancy<-trees[Year==j,SpeciesOccupancy]
  StandAgeC<-trees[Year==j,StandAgeC]
  StandID<-trees[Year==j,StandID]
  Species<-trees[Year==j,Species]
  deviation<-trees[Year==j,deviation]
  log_dbh<-trees[Year==j,log_dbh]
  Keta<-Peta<-eta<-rep(NA,times=nr)
  tree_connectivity<-rep(0,times=nr)

  if(ExactConnectivity>0){ for(r in which(stand_IDs %in% StandID)){ if(trees[Year==j & StandID==stand_IDs[r],StandAgeC][1]>colonizable){# This gives the row numbers of the stand_distances table that need to be read
  distances<-as.numeric(gsub(",",".",fread(stand_distances,skip=(incl_stands[r]-1),nrows=1,header=T)))[incl_stands] # Distances of the stand r from other stands
  nearby_stands<-stand_IDs[which(distances<=ExactConnectivity & (stand_IDs %in% trees[,StandID]))] # IDs of the stands of the year for which exact connectivities are calculated
  
  # Calculate exact connectivity of focal trees to trees of the stands that are within ExactConnectivity
  if(sum(NrSources[stand_IDs %in% nearby_stands])>0){ # If there are any source trees within nearby stands,
  connectivity<-fields::rdist(rbind(trees[Year==j & StandID==stand_IDs[r],.(x,y)],trees[(Year==(j-TimeStepLength) & StandID %in% nearby_stands & trees$SpeciesOccupancy>=Maturity),.(x,y)]),compact=FALSE) # dist matrix of present & source trees
  if(ncol(connectivity)>nrow(trees[Year==j & StandID==stand_IDs[r]])) connectivity<-matrix(connectivity[1:nrow(trees[Year==j & StandID==stand_IDs[r]]),(nrow(trees[Year==j & StandID==stand_IDs[r]])+1):(ncol(connectivity))],nrow=nrow(trees[Year==j & StandID==stand_IDs[r]])) # remove distances to present trees, works also with just one source tree
  if(("al" %in% colnames(Model))==TRUE & ("c" %in% colnames(Model))==TRUE)  tree_connectivity[StandID==stand_IDs[r]]<-rowSums(exp(-(Model_ii["al"]*connectivity)^exp(Model_ii["c"])))
  if(("al" %in% colnames(Model))==TRUE & ("c" %in% colnames(Model))==FALSE) tree_connectivity[StandID==stand_IDs[r]]<-rowSums(exp(-Model_ii["al"]*connectivity))
  rm(connectivity); gc()
  }
  
  # Calculate aggregated connectivity to far-away stands
  if(("al" %in% colnames(Model))==TRUE & ("c" %in% colnames(Model))==TRUE)  tree_connectivity[StandID==stand_IDs[r]]<-tree_connectivity[StandID==stand_IDs[r]]+sum(exp(-(Model_ii["al"]*distances[which(distances>ExactConnectivity)])^exp(Model_ii["c"]))*NrSources[which(distances>ExactConnectivity)])
  if(("al" %in% colnames(Model))==TRUE & ("c" %in% colnames(Model))==FALSE) tree_connectivity[StandID==stand_IDs[r]]<-tree_connectivity[StandID==stand_IDs[r]]+sum(exp(-Model_ii["al"]*distances[which(distances>ExactConnectivity)])*NrSources[which(distances>ExactConnectivity)])
  rm(distances); gc()
  }
  if(trees[Year==j & StandID==stand_IDs[r],StandAgeC][1]<=colonizable) tree_connectivity[StandID==stand_IDs[r]]<-0
  }}
    
  # Estimate connectivity of present trees to previous time step's source trees IF there were mature trees
  if(ExactConnectivity==0){ if(nrow(trees[(Year==(j-TimeStepLength) & SpeciesOccupancy>=Maturity),])>0){
  connectivity<-rdist(rbind(trees[trees$Year==j,.(x,y)],trees[(trees$Year==(j-TimeStepLength) & trees$SpeciesOccupancy>=Maturity),.(x,y)]),compact=FALSE) # dist matrix of present & source trees
  connectivity<-matrix(connectivity[1:nrow(trees[trees$Year==j]),(nrow(trees[trees$Year==j])+1):(ncol(connectivity))],nrow=nrow(trees[trees$Year==j])) # remove distances to present trees, works also with just one source tree
  if(("al" %in% colnames(Model))==TRUE & ("c" %in% colnames(Model))==TRUE)  tree_connectivity<-rowSums(exp(-(Model_ii["al"]*connectivity)^exp(Model_ii["c"])))
  if(("al" %in% colnames(Model))==TRUE & ("c" %in% colnames(Model))==FALSE) tree_connectivity<-rowSums(exp(-Model_ii["al"]*connectivity))
  rm(connectivity); gc()
  }}

  # Estimate Keta depending on colonization model variables. NOTE: This assumes that Kb0[1] is for non-shaded trees and Kb0[2] is for non-shaded trees
  if(sum(NrSources)>0){ # Estimate colonization rate only if there were dispersal sources available
  if(("Kb0" %in% colnames(Model))==TRUE) Keta<-Model_ii["Kb0"] + ifelse(log(tree_connectivity)>0,log(tree_connectivity),0)
  if(("Kb0[1]" %in% colnames(Model))==TRUE) Keta<-Model_ii["Kb0[1]"] * ifelse(StandAgeC<(colonizable),1,0) + Model_ii["Kb0[2]"] * ifelse(StandAgeC>=(colonizable),1,0) + ifelse(log(tree_connectivity)>0,log(tree_connectivity),0) # open vs. shaded
  
  if(("b0" %in% colnames(Model))==TRUE) Keta<-rep(Model_ii["b0"],times=nr) + ifelse(log(tree_connectivity)>0,log(tree_connectivity),0)
  if(("b1" %in% colnames(Model))==TRUE) Keta<-Keta + Model_ii["b1"] * StandAgeC
  if(("b2" %in% colnames(Model))==TRUE) Keta<-Keta + Model_ii["b2"] * StandAgeC^2
  if(("sigma.alpha" %in% colnames(Model))==TRUE) Keta<-Keta + deviation
  
  if(("beta_dbh" %in% colnames(Model))==TRUE) Keta<-Keta + Model_ii["beta"] * log_dbh
  Keta<-VGAM::cloglog(Keta,inv=T)
  }
  if(sum(NrSources)==0) trees[Year==j,"Keta"]<-0 # If there were no dispersal sources, there cannot be colonizations
  Peta<-inv.logit(Model_ii["Pb0"]) * ifelse(SpeciesOccupancy>0,1,0) # If the tree was occupied at the previous time step, its occupancy has already been added
  Keta<-1-((1-Keta)^(TimeStepLength/10))
  Peta<-Peta^(TimeStepLength/10)
  eta<-Keta*(1-ifelse(SpeciesOccupancy>0,1,0))+Peta*ifelse(SpeciesOccupancy>0,1,0)
  eta<-sapply(eta,function(x) rbern(1,x)) # Turn prediction into integers
  SpeciesOccupancy[eta==0]<-0 # trees go extinct or stay unoccupied according to prediction
  SpeciesOccupancy[SpeciesOccupancy==0 & eta==1]<-1 # unoccupied trees get colonized according to prediction. NOTE: This way trees that were already occupied don't lose their age
  trees[Year==j,"SpeciesOccupancy"]<-SpeciesOccupancy
  trees[Year==j,"connectivity"]<-log(tree_connectivity)
  NrSources[stand_IDs %in% StandID][order(stand_IDs[stand_IDs %in% StandID])]<-aggregate(ifelse(SpeciesOccupancy>=Maturity,1,0),by=list(StandID),FUN=sum)$x[order(unique(StandID))] # Sum of colonization sources per included stand
  rm(SpeciesOccupancy,StandAgeC,StandID,Species,log_dbh,Keta,Peta,eta,tree_connectivity); gc()
  
  # Write the trees vector to a file: fwrite(trees[Year==j],paste("Heureka_Lobaria_sim_results_",n,"_i_",i,"_ii_",ii,".txt",sep=""),append=T)
}
if(nrow(trees[(Year==(j-TimeStepLength) & SpeciesOccupancy>0),])==0){ # If the metapopulation has gone exinct, fill the table with zeros for easier post-processing
  trees[Year==j,"SpeciesOccupancy"]<-0
  trees[Year==j,"connectivity"]<-(-Inf)
}}

# Save results - if Lobaria has gone extinct, fill empty slots with 0 occupancy
trees[is.na(SpeciesOccupancy),"SpeciesOccupancy"]<-0
write.table(trees,paste("LPUL_",scenario,"_sim_results_3000_",n,"_i_",i,"_ii_",ii,".txt",sep=""))
}
