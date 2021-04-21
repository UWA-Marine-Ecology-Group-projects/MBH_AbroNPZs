###   ###   ###   MBH design clustered BRUVs    ###   ###   ###


# clear environment ----
rm(list = ls())

# libraries ----
#install.packages("MBHdesign")
library( MBHdesign)
library( parallel)
library( class)
library( fields)
#install.packages("pdist")
library( pdist)
library( raster)
library( rgdal)
library( sp)
library( rgeos)

# Directories ----
w.dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
#w.dir <- "~/MBH_AbroNPZs"
p.dir <- paste(w.dir, "plots", sep = '/')
d.dir <- paste(w.dir, "data", sep='/')
s.dir <- paste(w.dir, "shapefiles", sep='/')
r.dir <- paste(w.dir, "rasters", sep='/')
o.dir <- paste(w.dir, "outputs", sep='/')

####    NPZ 6   ####

# Read in the inclusion probs ----
inclProbs <- raster(paste(r.dir, "inclProbs_zone6.50_deployments.tif", sep='/'))
plot(inclProbs)
inclProbs <- setValues( inclProbs, values( inclProbs) / sum( values( inclProbs), na.rm=TRUE))
plot(inclProbs)

# check sun of incl probs --
cellStats(inclProbs, 'sum')

rootInclProbs <- inclProbs

rootInclProbs <- setValues( rootInclProbs, sqrt( values( rootInclProbs)))

cellStats(rootInclProbs, 'sum')
plot(rootInclProbs)


# Read data ----
zones <- readRDS(paste(d.dir, "Zones_Abro_NPZs.RDS", sep='/')) # this one in different folder
#Deans <- readRDS( "DeansPoints_forinNOutMP-d3.RDS")
rast <- readRDS(paste(d.dir, "abro_rasters_forInNOutNPZ.RDS", sep='/'))
#if( class( BRUVS) != "SpatialPointsDataFrame")
#Deans <- SpatialPointsDataFrame( coords=Deans[,c("Longitude","Latitude")], data=Deans, proj4string = CRS( proj4string( zones[[1]])))
#proj4string(Deans) <- proj4string(swrast$bathy)
straw.nums <- readRDS(paste(d.dir, "StrawmanNumbers_zones06.50_deployments.RDS", sep ='/'))


################################

###   ###   Part where legacy sites are added     ###     ###

####  choose reference sites
####  This is a one-step sample and hence
####  uses orignal inclProbs
################################

# working out density dependent probs of inclusion
# tmpD <- spTransform( Deans, CRS="+init=epsg:3577")
# tmpDist <- as.matrix( dist( coordinates( tmpD)))
# tmpDist1 <- apply( tmpDist, 1, function(x) sum( x<1000))
# Deans@data$sampleProbs <- 1 / sqrt( tmpDist1)  #to alter the importance a little bit
# Deans@data$sampleProbs <- Deans@data$sampleProbs / sum( Deans@data$sampleProbs, na.rm=TRUE)
# Deans@data$inclProbs <- raster::extract(inclProbs, Deans)
# Deans@data$sampleProbs <- Deans@data$sampleProbs * Deans@data$inclProbs #TPI inclusion probs are zero in most places

# numRef <- rep( NA, 2)
# names( numRef) <- c("InsideMP","OutsideMP")
# 
# for( zz in c( "InsideMP","OutsideMP")){
#   myZone <- zones[[zz]]
#   #if( zz == "MUS")
#   #myZone = zones$AMP - zones$IUCN2
#   tmpDrop <- as.vector( as.matrix( over( Deans, myZone)))
#   length1 <- length(tmpDrop)
#   wna <- which(is.na(tmpDrop))
#   lengthna <- length(wna)
#   length2 <- length1-lengthna
#   numRef[zz] <- min( floor( straw.nums[zz]/2), length2)
#   #numRef[zz] <- min( floor( straw.nums[zz]), length2)
#   #numRef[zz] <- min( floor( straw.nums[zz]/2), sum( length2, na.rm=TRUE))
#   prob <- Deans@data[!is.na( tmpDrop),"sampleProbs"]
#   prob <- na.exclude(prob)
#   prob <- as.vector(prob)
#   length3 <- length(prob)
#   
#   Deans@data[!is.na( tmpDrop), "sampleProbs"][sample.int(length3, numRef[zz], prob, replace=TRUE)] <- TRUE
#   #Deans@data[!is.na( tmpDrop), "f"][sample.int(length3, numRef[zz], prob, replace=FALSE)] <- TRUE
#   
#   #Deans@data[!is.na( tmpDrop), "f"][sample.int(length2, numRef[zz], prob=Deans@data[!is.na( tmpDrop),"sampleProbs"], replace=FALSE)] <- TRUE
#   #Deans@data[!is.na( tmpDrop), "f"][sample.int( sum( tmpDrop, na.rm=TRUE), numRef[zz], replace=FALSE)] <- TRUE
# }
# 
# Deans@data$sampleProbs
# 
# # load legacy sites ----
# #legacySites <- readOGR(dsn="C:/Users/00093391/Dropbox/UWA/Research Associate/MBHpackage/Ningaloo19_Data/legacySites_2019-12-23.shp")
# 
# # Deans Points are trated as legacy sites --
# Deans2 <- Deans
# table( Deans@data$f) / table( Deans2@data$f)
# 
# plot( inclProbs)
# plot( Deans, pch=20, col ='red', add=TRUE)
# #points( coordinates( Deans)[Deans@data$f,], col='red')
# 
# #legacySites <- Deans@data[Deans@data$f,]
# #legacySites <- SpatialPointsDataFrame( coords=legacySites[,c("Longitude","Latitude")], data=legacySites, proj4string=CRS(proj4string(inclProbs)))
# 
# legacySites <- Deans # 35 legacy sites




############################
####  Spatial sample of new sites ----
####  from altered incl. probs.
############################

### Here use quasiSamp to get random points ####
## these points will be the center of buffer for transects ###

####  Set the seed for reproducability
#set.seed( 777)
#### HAVE NOT BEEN ABLE TO MAKE THIS FUNCTION WORK ----
newSites <- list(npz6 = NULL, out6 = NULL)

for( zz in c("npz6", "out6")){
  print( zz)
  #the number of samples to take (specified minus the legacy number)
  #numby <- floor( (straw.nums[zz])/4)  # for clustered cluster - without legacy sites
  numby <- floor( (straw.nums[zz])) # for not clustered sites
  #numby <- floor( (straw.nums[zz] - numRef[zz])/2)
  #numby <- floor( (straw.nums[zz] - numRef[zz])) # with legacy sites 
  #set up spatial domain
  myZone <- zones[[zz]]
  #if( zz == "AMP"){
  # myZone = zones$AMP - zones$IUCN2
  #set.seed( 747)
  #}
  #tmpIP <- mask( rootInclProbs_agg_100m, myZone)
  tmpIP <- mask( inclProbs, myZone)
  tmpIP <- crop( tmpIP, myZone)
  #take the sample of clusters based on root incl probs
  newSites[[zz]] <- quasiSamp( n=numby, potential.sites=coordinates( tmpIP), inclusion.probs=values(tmpIP), nSampsToConsider=5000)
  
  #plotting (maybe remove at a later date?)
  tmpIPFull <- mask( inclProbs, myZone)
  tmpIPFull <- crop( tmpIPFull, myZone)
  plot( tmpIPFull)
  #plot( legacySites, add=TRUE, pch=1, col='red')
  points( newSites[[zz]][,c("x","y")], pch=20, col='black')
}
newSites <- do.call( "rbind", newSites)
head(newSites)

# Give id to sites and zones --
site.names <-row.names(newSites)
newSites$site <- as.factor(site.names)
#zone.names <- gsub('.{3}$', '', site.names) # remove last 3 characters
zone.names <- substr(site.names, 1, 3) # extract first three characters
newSites$zone <- as.factor(zone.names)
newSites$zone
newSites <- SpatialPointsDataFrame( coords=newSites[,c("x","y")], data=newSites, proj4string=CRS(proj4string(inclProbs)))
#some of the spatial balance is not great...  Presumably because the balance of the reference sites is also not great...

# Plot --
plot(inclProbs)
plot(rast$slope6, main = "Slope")
plot(zones$Both6, add=T)
plot(newSites, col=newSites$zone, pch = 20, add=T) # 41
newSites$zone



### Make sure the clusters centres are ~ 1 km apart ----


## Get CRS in utm ----
crs1 <- CRS("+init=epsg:32750") # WGS 84 / UTM zone 50S


## transform the points into UTM --
p1u <- spTransform(newSites, crs1)

## calculate if 2 points fall within 1500 m of eachother ----
# https://gis.stackexchange.com/questions/102796/remove-points-within-x-distance

dist1 <- gDistance(p1u, byid =T)
dist1
max(dist1)
min(dist1[dist1 > 0]) # minimum distance other than 0

## p1 ----
p1_matrix <- gWithinDistance(p1u, dist = 50, byid = TRUE)
diag(p1_matrix) <- NA
p1_matrix

# extract the upper triangular part of matrix and use the column sums as a criterion to remove the points:

p1_matrix[lower.tri(p1_matrix, diag=TRUE)] <- NA
p1_matrix

colSums(p1_matrix, na.rm=TRUE) == 0
v1 <- colSums(p1_matrix, na.rm=TRUE) == 0
p1u[v1, ] # 98 features left

remaining.sites <- p1u[v1, ]
remaining.sites <- spTransform(remaining.sites, proj4string(inclProbs))

# plot --
plot(inclProbs)
plot(rast$slope6, main = "Slope")
plot(zones$Both6, add=T)
plot(remaining.sites, col=remaining.sites$zone, pch = 20, add=T) # 41
remaining.sites$zone

## Save --
site <- "Abrolhos"
NPZ <- "npz6"
design <- "50Boss"

writeOGR(remaining.sites, o.dir, paste(site, NPZ, design, sep='-'), driver = "ESRI Shapefile")


