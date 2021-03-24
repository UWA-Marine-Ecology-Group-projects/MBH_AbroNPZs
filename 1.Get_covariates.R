library(raster)
library(rgdal)
library(sp)
library(ggplot2)
library(pals)


w.dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
w.dir <- "~/MBH_AbroNPZs"
p.dir <- paste(w.dir, "plots", sep = '/')
dt.dir <- paste(w.dir, "Data/Tidy", sep='/')
s.dir <- paste(w.dir, "shapefiles", sep='/')
r.dir <- paste(w.dir, "rasters", sep='/')

b <- raster(paste(r.dir, "WA_500m_Bathy.tif", sep='/'))
plot(b)

ab <- readOGR(paste(s.dir, "Abrolhos_AMP_part1.shp", sep='/'))
plot(ab, add=T)
head(ab)


# separate polys of NPZ ----
swabrnpz06 <- ab[ab$POLYGONID=="swabrnpz06",]
#plot(swabrnpz06)
writeOGR(swabrnpz06, s.dir, "Ab_NPZ06", driver = "ESRI Shapefile")

swabrnpz09 <- ab[ab$POLYGONID=="swabrnpz09",]
#plot(swabrnpz09, add=T)
writeOGR(swabrnpz09, s.dir, "Ab_NPZ09", driver = "ESRI Shapefile")

npzs <- raster::union(swabrnpz06, swabrnpz09)
plot(npzs)

writeOGR(npzs, s.dir, "Ab_NPZs", driver = "ESRI Shapefile")

# Bathy 06 ----
e <- drawExtent()
b2 <- crop(b, e)
plot(b2)
plot(swabrnpz06, add=T) # south NPZ
plot(swabrnpz09, add=T) # North NPZ


# Crop bathy for each NPZ ----

# NPZ06 ----
b06 <- crop(b, swabrnpz06)
b06 <- b06*(-1)
b09 <- crop(b, swabrnpz09)

bcuts <- c(0, 20, 40, 80, 100, 120, 140, 160, 180, 200, 220)
pal <- colorRampPalette(c("red","blue"))

plot(b06,
     breaks = bcuts,
     col = pal(10))

r <- b06 #raster object
#preparing raster object to plot with geom_tile in ggplot2
r_points = rasterToPoints(r)
r_df = data.frame(r_points)
head(r_df) #breaks will be set to column "layer"
r_df$cuts=cut(r_df$WA_500m_Bathy, breaks=c(0, 20, 40, 80, 100, 120, 140, 160, 180, 200, 220)) #set breaks

ggplot(data=r_df) + 
  geom_tile(aes(x=x,y=y,fill=cuts)) + 
  scale_fill_brewer("depth",type = "seq", palette = "Blues") +
  coord_equal() +
  theme_bw() +
  theme(panel.grid.major = element_blank()) +
  xlab("Longitude") + ylab("Latitude")

# NPZ09 ----
b09 <- crop(b, swabrnpz09)
b09 <- b09*(-1)


bcuts <- c(0, 20, 40, 80, 100, 120, 140, 160, 180, 200, 220)
pal <- colorRampPalette(c("red","blue"))

plot(b09,
     breaks = bcuts,
     col = pal(10))

r <- b09 #raster object
#preparing raster object to plot with geom_tile in ggplot2
r_points = rasterToPoints(r)
r_df = data.frame(r_points)
head(r_df) #breaks will be set to column "layer"
r_df$cuts=cut(r_df$WA_500m_Bathy, breaks=c(0, 20, 40, 80, 100, 120, 140, 160, 180, 200, 220)) #set breaks

ggplot(data=r_df) + 
  geom_tile(aes(x=x,y=y,fill=cuts)) + 
  scale_fill_brewer("depth",type = "seq", palette = "Blues") +
  coord_equal() +
  theme_bw() +
  theme(panel.grid.major = element_blank()) +
  xlab("Longitude") + ylab("Latitude")

plot(b09)


### Calculate slope, aspect, tpi ----

# NPZ 06
slope06 <- terrain(b06, 'slope')
plot(slope06)

aspect06 <- terrain(b06, 'aspect', unit = 'degrees')
plot(aspect06, col=rainbow(100))

tpi06 <- terrain(b06, 'TPI')
plot(tp106)

depth06 <- b06

ders6 <- stack(depth06, slope06, tpi06, aspect06)
names(ders6) <- c("depth09", "slope09", "tpi09", "aspect09")

writeRaster(ders6, paste(r.dir, "Ab_NPZ06_ders.tif", sep ='/'), overwrite = TRUE)



# NPZ 09
slope09 <- terrain(b09, 'slope')
plot(s09)

aspect09 <- terrain(b09, 'aspect', unit = 'degrees')
plot(aspect09, col=rainbow(100))

tpi09 <- terrain(b09, 'TPI')
plot(tpi09)

depth09 <- b09

ders9 <- stack(depth09, slope09, tpi09, aspect09)
names(ders9) <- c("depth09", "slope09", "tpi09", "aspect09")

writeRaster(ders9, paste(r.dir, "Ab_NPZ09_ders.tif", sep ='/'), overwrite = TRUE)

plot(depth06)
plot(depth09, add=T)

origin(depth06)
origin(depth09)

both.ders <- raster::merge(ders6, ders9)
plot(both.ders)
names(both.ders) <- c("depth", "slope", "tpi", "aspect")
writeRaster(both.ders, paste(r.dir, "Ab_NPZs_ders.tif", sep ='/'), overwrite = TRUE)
