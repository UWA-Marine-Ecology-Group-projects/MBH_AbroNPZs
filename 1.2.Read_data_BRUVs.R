###   ###   ###   Read data for MBH analysis    ###   ###   ###


# clear environment ----
rm(list = ls())


# libraries ----
library( rgdal)
library( sp)
library( raster)


# Directories ----
w.dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
#w.dir <- "~/MBH_AbroNPZs"
p.dir <- paste(w.dir, "plots", sep = '/')
d.dir <- paste(w.dir, "data", sep='/')
s.dir <- paste(w.dir, "shapefiles", sep='/')
r.dir <- paste(w.dir, "rasters", sep='/')



# Read polygons of NPZs and adjacent areas ----
npz6 <- readOGR(paste(s.dir, "Ab_NPZ06.shp", sep = '/'))
npz9 <- readOGR(paste(s.dir, "Ab_NPZ09.shp", sep = '/'))
out6 <- readOGR(paste(s.dir, "Area_next_npz6.shp", sep = '/'))
out9 <- readOGR(paste(s.dir, "Area_next_npz9.shp", sep = '/'))

# check crs --
proj4string(npz6)
proj4string(npz9)
proj4string(out6)
proj4string(out9)

# make same crs --
npz6 <- spTransform(npz6, proj4string(out6))
npz9 <- spTransform(npz9, proj4string(out6))

# Prepare list of polygons ----
zones <- list()
zones$npz6 <- npz6
zones$out6 <- out6
zones$npz9 <- npz9
zones$out9 <- out9
# join both polygons
zones$Both6 <- raster::union(zones$npz6, zones$out6)
zones$Both9 <- raster::union(zones$npz9, zones$out9)
zones$All <- raster::union(zones$Both6, zones$Both9)

#intial look to see area
plot( zones$All, border='black')
plot( zones$Both6, add=TRUE, col='orange')
plot( zones$Both9, add=TRUE, col='green')
plot( zones$npz6, add=TRUE, col='blue')
plot( zones$npz9, add=TRUE, col='red')

# Save zones rds ----

saveRDS(zones, file= paste(d.dir, "Zones_Abro_NPZs.RDS", sep='/'))


## Read raster data ----
ders06 <- stack(paste(r.dir, "Ab_NPZ06_ders.tif", sep ='/'))
plot(ders06)
names(ders06) <- c("depth", "slope", "tpi", "aspect")
b6 <- ders06$depth
s6 <- ders06$slope
ders09 <- stack(paste(r.dir, "Ab_NPZ09_ders.tif", sep ='/'))
plot(ders09)
b9 <- ders09$depth
s9 <- ders09$slope

# set names of derivatives --
der.names <- c("depth", "slope", "tpi", "aspect")

names(ders06) <- der.names
names(ders09) <- der.names

# Save rasters rds ----
abro_rasters <- list()
abro_rasters$depth6 <- b6
abro_rasters$slope6 <- s6
abro_rasters$depth9 <- b9
abro_rasters$slope9 <- s9
saveRDS(abro_rasters, file= paste(d.dir, "abro_rasters_forInNOutNPZ.RDS", sep='/'))



## Converting polygons to a common raster ----

###       ###       ### this takes a while for fine res data  ###      ###       ###
npz6_raster <- rasterize(x=zones$npz6, y=b6, field=zones$npz6@data[,1], bkg.value=NA, fun="first")
plot(npz6_raster)
npz9_raster <- rasterize(x=zones$npz9, y=b9, field=zones$npz9@data[,1], bkg.value=NA, fun="first")
plot(npz9_raster)
out6_raster <- rasterize(x=zones$out6, y=b6, field=zones$out6@data[,1], bkg.value=NA, fun="first")
plot(out6_raster)
out9_raster <- rasterize(x=zones$out9, y=b9, field=zones$out9@data[,1], bkg.value=NA, fun="first")
plot(out9_raster)


#convert and combine --
tmp1 <- as.data.frame( npz6_raster, xy=TRUE)
tmp2 <- as.data.frame( npz9_raster, xy=TRUE)
tmp3 <- as.data.frame( out6_raster, xy=TRUE)
tmp4 <- as.data.frame( out9_raster, xy=TRUE)
tmp5 <- as.data.frame( b6, xy=TRUE)
tmp6 <- as.data.frame( b9, xy=TRUE)
tmp7 <- as.data.frame( s6, xy=TRUE)
tmp8 <- as.data.frame( s9, xy=TRUE)

# Join data for NPZ6 and adjacent analysis --

# NPZ 6 --
npz6Dat <- cbind( tmp1, tmp3[,3])
npz6Dat <- cbind( npz6Dat, tmp5[,3])
npz6Dat <- cbind( npz6Dat, tmp7[,3])
head(npz6Dat)

# NPZ9 --
npz9Dat <- cbind( tmp2, tmp4[,3])
npz9Dat <- cbind( npz9Dat, tmp6[,3])
npz9Dat <- cbind( npz9Dat, tmp6[,3])
head(npz9Dat)

# Set column names --
df.names <- c("Eastern", "Northing", "npz6", "out6", "depth", "slope")

names(npz6Dat) <- df.names
names(npz9Dat) <- df.names


# Save raster dfs rds ----
saveRDS(npz6Dat, file= paste(d.dir, "npz6Dat_forInNOutNPZ.RDS", sep='/'))
saveRDS(npz9Dat, file=paste(d.dir, "npz9Dat_forInNOutNPZ.RDS", sep ='/'))






