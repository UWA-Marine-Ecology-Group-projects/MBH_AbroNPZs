# Get inclusion probabilities ---

## Create inclusion probabilities ####

library( rgdal)
library( sp)
library( raster)
library( rgeos)
library( sf)

# clear environment ----
rm(list = ls())

# Set working directories ----
w.dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
p.dir <- paste(w.dir, "plots", sep = '/')
s.dir <- paste(w.dir, "shapefiles", sep='/')
r.dir <- paste(w.dir, "rasters", sep='/')
o.dir <- paste(w.dir, "outputs", sep='/')
d.dir <- paste(w.dir, "data", sep='/')


###   ###     NPZ06        ####

# Read polygon of NPZ ----

npz6 <- readOGR(paste(s.dir, "Ab_NPZ06.shp", sep='/'))


# Read rasters ----
ders <- stack(paste(r.dir, "Ab_NPZ06_ders.tif", sep='/'))

ab_rasters <- list()

ab_rasters$depth <- ders$Ab_NPZ06_ders.1
ab_rasters$slope <- ders$Ab_NPZ06_ders.2
ab_rasters$tpi <- ders$Ab_NPZ06_ders.3
ab_rasters$aspect <- ders$Ab_NPZ06_ders.4


plot(ab_rasters$slope)


###################################################
#### converting polygons to a common raster.

r <- ab_rasters$slope
#plot( extent( r), add=TRUE)
#survArea first
npz_raster <- rasterize( x=npz6, y=r, field=npz6@data[,2], bkg.value=NA, fun="first")


###################################
#convert and combine
tmp1 <- as.data.frame( npz_raster, xy=TRUE)
tmp2 <- as.data.frame( ab_rasters$depth, xy=TRUE)
tmp3 <- as.data.frame( ab_rasters$slope, xy=TRUE)
tmp4 <- as.data.frame( ab_rasters$tpi, xy=TRUE)
tmp5 <- as.data.frame( ab_rasters$aspect, xy=TRUE)

ABDat <- cbind( tmp1, tmp2[,3])
ABDat <- cbind( ABDat, tmp3[,3])
ABDat <- cbind( ABDat, tmp4[,3])
ABDat <- cbind( ABDat, tmp5[,3])


colnames( ABDat) <- c("x", "y", "BATHY", "SLOPE", "TPI", "ASPECT")

# save data as RDS ----
saveRDS( ABDat, paste(d.dir, "ABData_forNPZ6.RDS", sep ='/'))
saveRDS( ab_rasters, paste(d.dir, "ABRasters_forNPZ6.RDS", sep='/'))
saveRDS( npz6,  paste(d.dir, "ABZones_forNPZ6.RDS", sep='/'))

# Read in data ----
ABDat <- readRDS(paste(d.dir, "ABData_forNPZ6.RDS", sep='/'))
ab_rasters <- readRDS(paste(d.dir,"ABRasters_forNPZ6.RDS", sep='/'))
zones <- readRDS(paste(d.dir, "ABZones_forNPZ6.RDS", sep='/'))

# Abrolhos NPZ6 polygon ----
plot(zones)

# Get stata based on slope ----

Bathy.quant <- c(0,0.5,0.85,0.95,1)
Bathy.cuts <- quantile( ab_rasters$slope, Bathy.quant)#c( -Inf,0.02,0.04,0.08,0.16,Inf)
Bathy.cuts
#trying to make it so there is no hand-picking (except for the hand-picked function)
tmp <- cumsum( Bathy.quant)

# Raster of strata ----
catB <- cut(ab_rasters$slope, breaks=Bathy.cuts, na.rm=TRUE)
plot(catB)
plot(zones, add=T)


## Proportions per strata ----
Bathy.targetProps <- c(0.2,0.2,0.3,0.3)


## Get inclusion probabilities ----


inclProbs <- catB

zoneID <- extract( x=catB, y=zones, cellnumbers=TRUE)
propsOfbathy <- table( catB@data@values[zoneID[[1]][,"cell"]])
propsOfbathy <- propsOfbathy / sum( propsOfbathy)
tmp <- Bathy.targetProps / propsOfbathy #the desired inclusion probs (unstandardised)
for( ii in 1:length( propsOfbathy)){
    inclProbs[zoneID[[1]][,"cell"]][zoneID[[1]][,"value"]==ii] <- tmp[ii]
  }
  inclProbs[zoneID[[1]][,"cell"]][is.na( inclProbs[zoneID[[1]][,"cell"]])] <- 0
  inclProbs[zoneID[[1]][,"cell"]] <- inclProbs[zoneID[[1]][,"cell"]] / sum( inclProbs[zoneID[[1]][,"cell"]])

inclProbs@data@values[inclProbs@data@values %in% c(0,1,2,3,4,5,6,7,8)] <- NA  #cheats way to crop
plot( inclProbs)


# save raster of inclusion probabilities ----
writeRaster(inclProbs,paste(d.dir, 'inclProbs_forNPZ6.tif', sep='/'), overwrite=TRUE)


# to check if sum of cells (probs) = 1
sumr <- cellStats(inclProbs, 'sum')
sumr


###   ###   ###     ###

##    MBH design   ####

numby <- 40 # how many samples

newSites <- quasiSamp( n=numby, potential.sites=coordinates(inclProbs), 
                       inclusion.probs=values(inclProbs), nSampsToConsider=5000)
newSites
# plot design
plot(inclProbs)
points( newSites[,c("x","y")], pch=20, col='black')

# Make points sp ----
newSites <- SpatialPointsDataFrame( coords=newSites[,c("x","y")], data=newSites, proj4string=CRS(proj4string(inclProbs)))
writeOGR(newSites, o.dir, "AB_NPZ6_BOSS", driver = "ESRI Shapefile")


### Make sure the drops are ~ 50 m apart ----

inclProbs <- raster(paste(d.dir, "inclProbs_forBOSS.tif", sep='/'))


#s <- readOGR(paste(o.dir, "GB_BOSS_d1.shp", sep='/'))
s # make sure is in UTM
s <- newSites
proj4string(s) <- proj4string(inclProbs)
coordinates(s) <- ~x+y
proj4string(s2)  
s2 <- spTransform(s, CRSobj = "+proj=utm +zone=50 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")

## calculate if 2 points fall within 25m of eachother ----
# https://gis.stackexchange.com/questions/102796/remove-points-within-x-distance


p1_matrix <- gWithinDistance(s2, dist = 50, byid = TRUE)
diag(p1_matrix) <- NA
p1_matrix

# extract the upper triangular part of matrix and use the column sums as a criterion to remove the points:
p1_matrix[lower.tri(p1_matrix, diag=TRUE)] <- NA
p1_matrix

colSums(p1_matrix, na.rm=TRUE) == 0
v1 <- colSums(p1_matrix, na.rm=TRUE) == 0
s2[v1, ] # 40 features left


## add unique id to points ----
IDnumber <- paste0(1:40)
s$IDnumber <- IDnumber
s


writeOGR(s, o.dir, "Ab_NPZ6_BOSS", driver = "ESRI Shapefile", overwrite = TRUE)
