library(lidR)
library(sf)
library(raster)
library(ForestTools)
library(readxl)
library(lwgeom)
library(ggplot2)

# Load the manual delineation and transform to CrownList
ManualCrownsFN <- "../data/manual_seg/manual_labels.shp"  
ManualCrowns <- st_read(ManualCrownsFN)
ManualCrowns$UniqueID <- 1:length(ManualCrowns$UniqueID)
ManualCentroids <- st_centroid(ManualCrowns)


# calculate the maximum radius, this takes a while to run (~3 mins)
ManualMaxRadius <- sapply(1:length(ManualCrowns$UniqueID),
                  function(i){
                    max((ManualCrowns[i,] %>% 
                      st_cast("MULTIPOINT") %>% 
                      st_cast("POINT") %>% 
                      st_distance(ManualCentroids[i,]))[[1]])
                    }
                  )

# load the CHM
CHM <- raster("../data/CHM/lfdp_chm.tif")


# Extract tree height from CHM
# this takes a while to run (10+ mins)
# one may also use the height of the centroids as tree heights, but that is slightly smaller 
ManualHeight <- sapply(1:length(ManualCrowns$UniqueID),
                       function(x){max(extract(CHM, ManualCrowns[x,])[[1]],na.rm = T)})



# Build a dataframe for the crowns with coordinate and other info
ManualCrownList = data.frame(ID = ManualCrowns$UniqueID,
                             X = st_coordinates(ManualCentroids)[,1],
                             Y = st_coordinates(ManualCentroids)[,2],
                             Height = ManualHeight,
                             CD = ManualMaxRadius * 2)
 

#### Build allometric models using AlloData dataset ####

# Load ground obs and data for allometric modelling
AlloData <- read_excel("./data/ground_obs/Crown allometry summary.xlsx")
GroundObs <- read.csv("./data/ground_obs/tree_census_2016_georeferenced_20cm_trans.csv")

# remove invalid rows
AlloData <- AlloData[(AlloData$DIAM9899 > 0) & (AlloData$MAXRAD > 0) & (AlloData$HEIGHT > 0),]

# establish species list of interest
species_list <-  c("CECSCH","CASARB","DACEXC","GUAGUI","MANBID","PREMON","SCHMOR","SLOBER","OTHERS")
AlloData$SPP[!AlloData$SPP %in% species_list] <- "OTHERS"

Allo_models <- list()
sigmas <- rep(0,length(species_list))
for (i in 1:length(species_list)){
  # subset dataset by species
  species <- species_list[i]
  AlloData_species <- AlloData[AlloData$SPP == species,]
  # calculate the log transformation of tree parameters
  logDBH <- log(AlloData_species$DIAM9899)
  logH <- log(AlloData_species$HEIGHT)
  logCD <- log(AlloData_species$MAXRAD*2)
  # build the model
  Allo_models[[i]] <- lm(logDBH ~ logH + logCD)
  sigmas[i] <- summary(Allo_models[[i]])$sigma
}


# allometric model at overall level without species info
logDBH <- log(AlloData$DIAM9899)
logH <- log(AlloData$HEIGHT)
logCD <- log(AlloData$MAXRAD*2)
Allo_model_overall <- lm(logDBH ~ logH + logCD)

# compare residuals
residual_species <- NULL
for (i in 1:9) {residual_species <- c(residual_species,Allo_models[[i]]$residuals)}

residual_overall <- Allo_model_overall$residuals
sigma_overall <- summary(Allo_model_overall)$sigma

#### Crown-stem matching #####

# Crop and get the ground observation in the same area
if (!is.null(bounds)){
  index <- (GroundObs$x_trans>bounds[1]) & (GroundObs$y_trans>bounds[3]) & (GroundObs$x_trans<bounds[2]) & (GroundObs$y_trans<bounds[4])
  Cropped_GroundObs <- GroundObs[index,]
} else {
  Cropped_GroundObs = GroundObs
}
Cropped_GroundObs <- Cropped_GroundObs[order(Cropped_GroundObs[,"DIAM"],decreasing = T),]

# 2.1 GET THE RELATIVE DISTANCE
Ncrowns <- dim(CrownList)[1]
Nstems <- dim(Cropped_GroundObs)[1]
if (Ncrowns < Nstems) {
  Cropped_GroundObs[1:Ncrowns,]
  Nstems <- Ncrowns
}
DMatrix <- matrix(nrow=Ncrowns, ncol=Nstems)
for (i in 1:Ncrowns){
  for (j in 1:Nstems){
    DMatrix[i,j] <- sqrt(sum((CrownList[i,c('X','Y')] - Cropped_GroundObs[j,c('x_trans','y_trans')])^2))
  }
}

alpha = 1
RankMatrix = apply(DMatrix,2,rank)
RDMatrix <- alpha * exp(-alpha * RankMatrix)


LMatrix <- matrix(nrow=Ncrowns, ncol=Nstems)
LMatrix_overall <- LMatrix
for (i in 1:Ncrowns){
  if (! i %% 100) {print(i)}
  for (j in 1:Nstems){
    templogDBH <- log(Cropped_GroundObs$DIAM[j])
    templogH <- log(CrownList[i,"Height"])
    templogCD <- log(CrownList[i,"CD"])
    if (Cropped_GroundObs$SPECIES[j] %in% species_list) {
      index <- which(species_list==Cropped_GroundObs$SPECIES[j])} else {
        index <- length(species_list)
      }
    HatlogDBH <- predict(Allo_models[[index]],newdata=list(logH=templogH,logCD=templogCD))
    HatlogDBH_overall <- predict(Allo_model_overall,newdata=list(logH=templogH,logCD=templogCD))
    
    LMatrix[i,j] <-  exp(-((HatlogDBH-templogDBH)/sigma_overall)^2/2)
    LMatrix_overall[i,j] <-  exp(-((HatlogDBH_overall-templogDBH)/sigma_overall)^2/2)
  }
}

# FINAL MATCHING
matching <- rep(0,Nstems)
ScoreMatrix <- LMatrix*RDMatrix
max_scores <- rep(0,Nstems)
for (j in 1:Nstems){
  if (max(ScoreMatrix[,j])>0){
    matching[j] <- which.max(ScoreMatrix[,j])[1]
    if (DMatrix[matching[j],j]<=20){
      ScoreMatrix[matching[j],] <- rep(0,Nstems) # zero the score of the assigned crowns, so that it will not be assigned again repeatly
    } else {matching[j] = 0}
  }}

index = which(matching > 0)
dist_rank = rep(0,Nstems) 
dist = rep(0,Nstems)
for (i in index) {dist_rank[i] =RankMatrix[matching[i],i]}
for (i in index) {dist[i] = DMatrix[matching[i],i]}
matching_overall <- rep(0,Nstems)

ScoreMatrix_overall <- LMatrix_overall*RDMatrix
for (j in 1:Nstems){
  matching_overall[j] <- which.max(ScoreMatrix_overall[,j])[1]
  ScoreMatrix_overall[matching_overall[j],] <- rep(0,Nstems)
}


##### PLOTTING #####

# distance plots
for (i in index) {dist_rank[i] =RankMatrix[matching[i],i]}
dist_rank = dist_rank[index]
for (i in index) {dist[i] =DMatrix[matching[i],i]}
dist = dist[index]

p1 = ggplot(as.data.frame(dist),aes(x=dist)) + geom_histogram() + xlab("Crown-stem distance")
p2 = ggplot(as.data.frame(dist_rank),aes(x=dist_rank)) + geom_histogram() + xlab("Distance ranking")


# allometric plots
M = matrix(rep(0,length(index)*6), ncol=6)
for (i in index) {
  M[count,] = c(i, 
                matching[i], 
                CrownList$Height[matching[i]], 
                CrownList$CD[matching[i]],
                Cropped_GroundObs$DIAM[i], 
                Cropped_GroundObs$SPECIES[i])
  count = count + 1
}
M = as.data.frame(M)

colnames(M) = c("stem_ID", "crown_ID", "Height","CD","DBH","Species")
M$Species[!M$Species %in% species_list] <- "OTHERS"
M$Height = as.numeric(M$Height)
M$DBH = as.numeric(M$DBH)
M$CD = as.numeric(M$CD)
p3 = ggplot(data = M, aes(x = log(DBH), y = log(CD))) + geom_point(size=0.3) + geom_smooth(method=lm) + facet_wrap(~Species)
p4 = ggplot(data = M, aes(x = log(DBH), y = log(Height))) + geom_point(size=0.3) + geom_smooth(method=lm) + facet_wrap(~Species)


