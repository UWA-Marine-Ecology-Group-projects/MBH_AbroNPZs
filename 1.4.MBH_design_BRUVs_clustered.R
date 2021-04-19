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


####    NPZ 6   ####

# Read in the inclusion probs ----
inclProbs <- raster(paste(r.dir, "inclProbs_zone6.tif", sep='/'))
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
straw.nums <- readRDS(paste(d.dir, "StrawmanNumbers_zones06.RDS", sep ='/'))


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
  numby <- floor( (straw.nums[zz])/4)  # 5 is the number of drops per cluster - without legacy sites
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
zone.names <- gsub('.{3}$', '', site.names)
newSites$zone <- as.factor(zone.names)
newSites <- SpatialPointsDataFrame( coords=newSites[,c("x","y")], data=newSites, proj4string=CRS(proj4string(inclProbs)))
#some of the spatial balance is not great...  Presumably because the balance of the reference sites is also not great...

# Plot --
plot(inclProbs)
plot(newSites, col=newSites$zone, pch = 20, add=T) # 41






### Make sure the clusters centres are ~ 1 km apart ----

# Join Deans and New sites data --
# newSites$class <- "MBH"
# Deans$class <- "Legacy"
# all <- union(newSites, Deans)
# all # 80 clusters

## Get CRS in utm ----
crs1 <- CRS("+init=epsg:32750") # WGS 84 / UTM zone 50S


## transform the points into UTM --
p1u <- spTransform(newSites, crs1)

## calculate if 2 points fall within 1500 m of eachother ----
# https://gis.stackexchange.com/questions/102796/remove-points-within-x-distance

## p1 ----
p1_matrix <- gWithinDistance(p1u, dist = 1500, byid = TRUE)
diag(p1_matrix) <- NA
p1_matrix

# extract the upper triangular part of matrix and use the column sums as a criterion to remove the points:

p1_matrix[lower.tri(p1_matrix, diag=TRUE)] <- NA
p1_matrix

colSums(p1_matrix, na.rm=TRUE) == 0
v1 <- colSums(p1_matrix, na.rm=TRUE) == 0
p1u[v1, ] # 98 features left

# plot --
# plot(inclProbs)
# plot(all[v1, ], pch=20, col='black', add=T)
# 
# # save points in utm and latlong ----
# ## Save points of design ---
# all2 <- all[v1, ]
# all2 # 70 clusters
# # transform again
# all2u <- spTransform(all2, gacrs)
# 
# plot(inclProbs)
# plot(all2[all2$class=="MBH",], pch=20, col='black', add=T)
# plot(all2[all2$class=="Legacy",], pch=20, col='red', add=T)
# 
# #writeOGR(all2, d.dir, "Allclusters-filtered-d7", driver="ESRI Shapefile")
# 
# ## separate Deans and MBH again
# 
# MBHsites <- all2[all2$class=="MBH",] #41
# Deansites <- all2[all2$class=="Legacy",] #29




###############################
####  Choose new points within clusters
####  Here I need to choose transects not points
##############################

getlocal <- function(ii){
  point <- newSites[ii,c("x","y")]
  r2 <- rasterize( point, inclProbs, field=1)
  pbuf <- buffer( r2, width=900) ## units are in metres
  buf <- mask( inclProbs, pbuf)
  buffer <- trim(buf, pad=0)
  return( buffer)
}


sampWithOver <- 25


fullSample <- list()
fullZones <- list()

## Get bruvs accoring to clusters ---
for( ii in 1:nrow(newSites)){
  tmp <- getlocal(ii)
  fullZones[[ii]] <- rownames( newSites@data)[ii]
  tmpm <- raster::as.matrix(tmp)
  tmpm <- t(tmpm)
  tmpdf <- as.data.frame (
    cbind (coordinates (tmp), as.numeric (tmpm)))
  colnames(tmpdf) <- c("x", "y", "inclProbs_design1")
  tmpdf <- tmpdf[ order(tmpdf$y, tmpdf$x),]  # order ascending first by northing and then by easting
  fullSample[[ii]] <- quasiSamp( n=sampWithOver, potential.sites=coordinates(tmp), inclusion.probs=values(tmp), nSampsToConsider=10000)
  fullSample[[ii]]$site <- newSites$site[ii]
  fullSample[[ii]]$zone <- newSites$zone[ii]
  #fullSample[[ii]] <- transectSamp( n=sampWithOver, potential.sites = tmpdf[,c(1,2)], 
  #potential.sites= tmpdf[,c("x","y")],
  #inclusion.probs= incprobdf[,3],
  #inclusion.probs= tmpdf[,3],
  #control=gb.control
  #constrainedSet=gb.constraints.bool
  
  plot( tmp)
  points( fullSample[[ii]]$points[,c("x","y")], pch=20, col='red')
  #plot( legacySites, add=TRUE, pch=4, col='blue')
}


fullSample <- do.call( "rbind", fullSample[1:length(fullSample)])
head(fullSample)
bruvs <- SpatialPointsDataFrame( coords=fullSample[,c("x","y")], data=fullSample, proj4string=CRS(proj4string(inclProbs)))

plot(inclProbs)
points(bruvs, col=bruvs$site, pch =20)
plot(newSites, col = 'black', add=T)


# Check that point are more than 400 m apart ----

## Get CRS in utm ----
crs1 <- CRS("+init=epsg:32750") # WGS 84 / UTM zone 50S


## transform the points into UTM --
p1u <- spTransform(bruvs, crs1)

## calculate if 2 points fall within 1400 m of eachother ----
# https://gis.stackexchange.com/questions/102796/remove-points-within-x-distance

## p1 ----
p1_matrix <- gWithinDistance(p1u, dist = 400, byid = TRUE)
diag(p1_matrix) <- NA
p1_matrix

# extract the upper triangular part of matrix and use the column sums as a criterion to remove the points:

p1_matrix[lower.tri(p1_matrix, diag=TRUE)] <- NA
p1_matrix

colSums(p1_matrix, na.rm=TRUE) == 0
v1 <- colSums(p1_matrix, na.rm=TRUE) == 0
x <- p1u[v1, ] # features left within right distance
x
# transform again from ploting
x2 <- spTransform(x, proj4string(inclProbs))

plot(inclProbs)
plot(x2, pch =20, col = x2$site, add=T)
plot(newSites, col = 'blue', add=T)








#### Try another way to filter points furhter than 400m apart ----
duplicated(p1u$ID) # check for duplicates in IDs
uniqueID <- paste0("bruv", 1:120) # create unique Id vector
p1u$uniqueID <- uniqueID

# split bruvs by site --
sites.list <- split(p1u, p1u$site, drop = F)


dist1 <- gDistance(sites.list[[1]], byid =T)
dist1
dist1[,1]
if(min(dist1[,1]) > 400){
  print("no")
} else {
  print("yes")
}

dist1[,1]


if(min(dist1[,1]) > 400){
  dist1[,1]
} else {
  dist1[,-1]
}


######################

p1_matrix <- gWithinDistance(sites.list[[2]], dist = 400, byid = TRUE)
diag(p1_matrix) <- NA
p1_matrix
siteID <- sites.list[[1]]$uniqueID
#colnames(p1_matrix) <- siteID
#rownames(p1_matrix) <- siteID
#p1_matrix <- cbind(replacement.col, p1_matrix)
p1_matrix
diag(p1_matrix) <- NA

# the test --

#replacement.col <- replicate(14, NA)

close.points <- vector()

if(any(p1_matrix[,2] == 'TRUE', na.rm = T) == TRUE){
  sites.list2 <- sites.list[[1]][-2,] # recalculate without the problem sites
  #current.site <- sites.list[[1]]$uniqueID[10]
  px_matrix <- gWithinDistance(sites.list2, dist = 400, byid = TRUE)
  diag(px_matrix) <- NA
  #replacement.col <- replicate(repl.no, "NA")
  #p1_matrix <- cbind(px_matrix[,1:(2-1)], replacement.col, p1_matrix[,2:14])
  index.no <- paste0(1)
  close.points[[1]]=index.no
  #p2_matrix <- p1_matrix
} else {
  p1_matrix <- p1_matrix
}
p1_matrix

# Test 2 --

if(any(p1_matrix[,2] == 'TRUE', na.rm = T) == TRUE){
  #sites.list2 <- sites.list[[1]][-2,] # recalculate without the problem sites
  #current.site <- sites.list[[1]]$uniqueID[10]
  remove.row <- which(p1_matrix[,2] == 'TRUE')
  p1_matrix <- p1_matrix[-remove.row,]
  #diag(p1_matrix) <- NA
  #replacement.col <- replicate(repl.no, "NA")
  #p1_matrix <- cbind(px_matrix[,1:(2-1)], replacement.col, p1_matrix[,2:14])
  # index.no <- paste0(1)
  # close.points[[1]]=index.no
  #p2_matrix <- p1_matrix
} else {
  p1_matrix <- p1_matrix
}


p1_matrix <- gWithinDistance(sites.list[[4]], dist = 400, byid = TRUE)
diag(p1_matrix) <- NA
p1_matrix



for( ii in 1:ncol(p1_matrix)){
if(any(p1_matrix[,ii] == 'TRUE', na.rm = T) == TRUE){
  remove.row <- which(p1_matrix[,ii] == 'TRUE')
  p1_matrix <- p1_matrix[-remove.row,]
} else {
  p1_matrix <- p1_matrix
}
}

p1_matrix <- gWithinDistance(sites.list[[4]], dist = 400, byid = TRUE)
diag(p1_matrix) <- NA
p1_matrix

p1_matrix
rownames(p1_matrix) <- paste0('row', rownames(p1_matrix))
rowname.ls=rownames(p1_matrix)

colnames(p1_matrix)

for( ii in colnames(p1_matrix)){
  if(any(p1_matrix[,ii] == 'TRUE', na.rm = T) == TRUE){
    p1_matrix <- p1_matrix[-(which(rownames(p1_matrix)==paste0('row',ii))),]
  } else {
    p1_matrix <- p1_matrix
  }
}

p1_matrix

### Test 3 ----

p1_matrix <- gWithinDistance(sites.list[[4]], dist = 400, byid = TRUE)
diag(p1_matrix) <- NA
p1_matrix

no.trues <- colSums(p1_matrix, na.rm=TRUE) # calculate no. of trues
col.order <- order(no.trues, decreasing = FALSE) # order cols by no. of trues
order.mat <- p1_matrix[,col.order] # reorder matrix according to no. of trues
no.trues2 <- colSums(order.mat, na.rm=TRUE) # check
order.mat

order.mat
rownames(order.mat) <- paste0('row', rownames(p1_matrix))
rowname.ls=rownames(order.mat)
order.mat

for( ii in 1:5){
  if(any(order.mat[,ii] == 'TRUE', na.rm = T) == TRUE){
    # remove row
    remove.row <- which(order.mat[,ii] == 'TRUE')
    order.mat <- order.mat[-remove.row,]
    # remove col
    remove.col.names <- names(remove.row)
    remove.col <- gsub('row', '', remove.col.names)
    #p1_matrix <- p1_matrix[-(which(rownames(p1_matrix)==paste0('row',ii))),]
  } else {
    order.mat <- order.mat
  }
}

order.mat

colSums(order.mat, na.rm =TRUE) == 0
v1 <- order.mat[1:4,]
v2 <- rownames(v1)
keep.rows <- gsub('row', '', v2)
x <- sites.list[[4]][keep.rows,]
x
plot(x)
test <- gWithinDistance(x, dist = 400, byid = TRUE)
test



v1 <- colSums(order.mat, na.rm=TRUE) == 0
v1.col.name <- names(v1)
keep.cols <- v1.col.name[1:4]
x <- sites.list[[4]][keep.cols,]
x <- sites.list[[4]][c(1:4), ] # choose the fi
plot(x)
test <- gWithinDistance(x, dist = 400, byid = TRUE)
test

remove.row <- which(order.mat[,2] == 'TRUE')
order.mat <- order.mat[-remove.row,]



rownames(p1_matrix) <- rownames(p1_matrix)
no.trues <- colSums(p1_matrix, na.rm=TRUE) 
col.order <- order(no.trues, decreasing = FALSE)
order.mat <- p1_matrix[,col.order] # reorder matrix according to no. of trues
no.trues2 <- colSums(order.mat, na.rm=TRUE) # check

max.col <- which(no.trues==(max(no.trues)))
many.col <- which(no.trues > 3 )
remove.row <- names(many.col)
p1_matrix <- p1_matrix[-many.col,]
p1_matrix <- p1_matrix[-(which(no.trues==(max(no.trues)))),]
p1_matrix <- p1_matrix[-(which(rownames(p1_matrix)==paste0('row',ii))),]


for( ii in colnames(p1_matrix)){
  if(any(p1_matrix[,ii] == 'TRUE', na.rm = T) == TRUE){
    p1_matrix <- p1_matrix[-(which(rownames(p1_matrix)==paste0('row',ii))),]
  } else {
    p1_matrix <- p1_matrix
  }
}

 # the loop 

for( ii in 1:ncol(p1_matrix)){
  if(ii < 15) {
  if(any(p1_matrix[,ii] == 'TRUE', na.rm = T) == TRUE){
  sites.list2 <- sites.list[[1]][-ii,]
  #current.site <- sites.list[[1]]$uniqueID[10]
  p1_matrix <- gWithinDistance(sites.list2, dist = 400, byid = TRUE)
  diag(p1_matrix) <- NA
  #replacement.col <- replicate(repl.no, "NA")
  p1_matrix <- cbind(p1_matrix[,1:(ii-1)], replacement.col, p1_matrix[,ii:14])
  p2_matrix <- p1_matrix
  } else {
    p1_matrix <- p1_matrix
    #p2_matrix <- p1_matrix
  }
  }
    else {
      if(any(p1_matrix[,ii] == 'TRUE', na.rm = T) == TRUE){
        sites.list2 <- sites.list[[1]][-ii,]
        #current.site <- sites.list[[1]]$uniqueID[10]
        p1_matrix <- gWithinDistance(sites.list2, dist = 400, byid = TRUE)
        diag(p1_matrix) <- NA
        #replacement.col <- replicate(repl.no, "NA")
        p1_matrix <- cbind(p2_matrix[,1:(ii-1)], replacement.col)
        p2_matrix <- p1_matrix
      } else {
        p1_matrix <- p1_matrix
      }
      
    }
}


p1_matrix
left1 <- colSums(p1_matrix, na.rm=TRUE) == 0

first4pts <- p1_matrix[,1:4] 


# extract the upper triangular part of matrix and use the column sums as a criterion to remove the points:

p1_matrix[lower.tri(p1_matrix, diag=TRUE)] <- NA
p1_matrix

colSums(p1_matrix, na.rm=TRUE) == 0
v1 <- colSums(p1_matrix, na.rm=TRUE) == 0
x <- p1u[v1, ] # features left within right distance




