
###   ###   ###    Prepare spatial data for analysis    ###   ###   ###

# Libraries ----
library(raster)
library(rgdal)
library(sp)
library(ggplot2)
library(pals)
library(RColorBrewer)

# Directories ----
w.dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
#w.dir <- "~/MBH_AbroNPZs"
p.dir <- paste(w.dir, "plots", sep = '/')
d.dir <- paste(w.dir, "data", sep='/')
s.dir <- paste(w.dir, "shapefiles", sep='/')
r.dir <- paste(w.dir, "rasters", sep='/')


## Load bathy and park data ----

# load West coast bathy --
b <- raster(paste(r.dir, "WA_500m_Bathy.tif", sep='/'))
plot(b)

# load part of Abrolhos Marine Park --
ab <- readOGR(paste(s.dir, "Abrolhos_AMP_part1.shp", sep='/'))
plot(ab, add=T)
head(ab)

e <- extent(112.6301, 113.8376, -28.22494, -26.75744) #  (xmin, xmax, ymin, ymax)
b.cropped <- crop(b, e)
plot(b.cropped)
plot(ab, add=T)



# Separate polys of NPZ ----
# swabrnpz06 <- ab[ab$POLYGONID=="swabrnpz06",]
# #plot(swabrnpz06)
# writeOGR(swabrnpz06, s.dir, "Ab_NPZ06", driver = "ESRI Shapefile")
# 
# swabrnpz09 <- ab[ab$POLYGONID=="swabrnpz09",]
# #plot(swabrnpz09, add=T)
# writeOGR(swabrnpz09, s.dir, "Ab_NPZ09", driver = "ESRI Shapefile")
# 
# npzs <- raster::union(swabrnpz06, swabrnpz09)
# plot(npzs)



# Save --
#writeOGR(npzs, s.dir, "Ab_NPZs", driver = "ESRI Shapefile")

# Load NPZ polys ----
npz06 <- readOGR(paste(s.dir, "Ab_NPZ06.shp", sep='/'))
npz09 <- readOGR(paste(s.dir, "Ab_NPZ09.shp", sep='/'))
npz.both <- readOGR(paste(s.dir, "Ab_NPZs.shp", sep='/'))


# Bathy 06 ----
plot(b.cropped)
plot(npz.both, add=T)

e <- extent(113.1751, 113.6251, -28.21244, -27.92744)
b.06 <- crop(b.cropped, e)
plot(b.06)
s.06 <- terrain(b.06, "slope")
plot(s.06, main = "Slope")
plot(npz06, add=T)

bcuts <- c(-40, -60, -80, -100, -120, -140, -160, -180, -200, -220)
pal <- colorRampPalette(c("dark blue",  "light blue", "white"))

raster::plot(b.06,
             breaks = bcuts,
             col = pal(10),
             main = "Depth")
plot(npz06, add=T)

# Plottting with ggplot ----

# r <- b.06 #raster object
# #preparing raster object to plot with geom_tile in ggplot2
# r_points = rasterToPoints(r)
# r_df = data.frame(r_points)
# head(r_df) #breaks will be set to column "layer"
# r_df$cuts=cut(r_df$WA_500m_Bathy, breaks=c(-40, -60, -80, -100, -120, -140, -160, -180, -200, -220)) #set breaks
# 
# ggplot(data=r_df) + 
#   geom_tile(aes(x=x,y=y,fill=cuts)) + 
#   scale_fill_brewer("depth",type = "seq", palette = "Blues") +
#   coord_equal() +
#   theme_bw() +
#   theme(panel.grid.major = element_blank()) +
#   xlab("Longitude") + ylab("Latitude")

## Bathy 09 ----
plot(b.cropped)
plot(npz.both, add=T)

e <- extent(112.9476, 113.3876, -27.28244, -26.97244)
b.09 <- crop(b.cropped, e)
plot(b.09)
s.09 <- terrain(b.09, "slope")
plot(s.09, main = "Slope")
plot(npz09, add=T)

bcuts <- c(-40, -60, -80, -100, -120, -140, -160)
pal <- colorRampPalette(c("dark blue",  "light blue", "white"))

raster::plot(b.09,
     breaks = bcuts,
     col = pal(10),
     main = "Depth")
plot(npz09, add=T)

# Plottting with ggplot ----

# r <- b.09 #raster object
# #preparing raster object to plot with geom_tile in ggplot2
# r_points = rasterToPoints(r)
# r_df = data.frame(r_points)
# head(r_df) #breaks will be set to column "layer"
# r_df$cuts=cut(r_df$WA_500m_Bathy, breaks=c(-40, -60, -80, -100, -120, -140, -160)) #set breaks
# 
# ggplot(data=r_df) + 
#   geom_tile(aes(x=x,y=y,fill=cuts)) + 
#   scale_fill_brewer("depth",type = "seq", palette = "Blues") +
#   coord_equal() +
#   theme_bw() +
#   theme(panel.grid.major = element_blank()) +
#   xlab("Longitude") + ylab("Latitude")


# Load adjacent polys ----
p06 <- readOGR(paste(s.dir, "Area_next_npz6.shp", sep='/'))
p09 <- readOGR(paste(s.dir, "Area_next_npz9.shp", sep='/'))


# Join npz and adjacent area polys ----

# NPZ 6 --
plot(b.cropped)
plot(npz06, add= T)
plot(p06, add=T)


npz06 <- spTransform(npz06, proj4string(b))
zone6 <- union(npz06, p06)
plot(zone6)

# crop bathy for zone 6 --
b.zone6 <- crop(b.cropped, zone6)
plot(b.zone6)

# NPZ 9 --
plot(b.cropped)
plot(npz09, add= T)
plot(p09, add=T)


npz09 <- spTransform(npz09, proj4string(b))
zone9 <- union(npz09, p09)
plot(zone9)

# crop bathy for zone 9 --
b.zone9 <- crop(b.cropped, zone9)
plot(b.zone9)



### Calculate slope, aspect, tpi ----

# NPZ 06
slope06 <- terrain(b.zone6, 'slope')
plot(slope06)

aspect06 <- terrain(b.zone6, 'aspect', unit = 'degrees')
plot(aspect06, col=rainbow(100))

tpi06 <- terrain(b.zone6, 'TPI')
plot(tpi06)

depth06 <- b.zone6

ders6 <- stack(depth06, slope06, tpi06, aspect06)
names(ders6) <- c("depth09", "slope09", "tpi09", "aspect09")

writeRaster(ders6, paste(r.dir, "Ab_NPZ06_ders.tif", sep ='/'), overwrite = TRUE)



# NPZ 09
slope09 <- terrain(b.zone9, 'slope', neighbors = 8)
plot(slope09)

aspect09 <- terrain(b.zone9, 'aspect', unit = 'degrees')
plot(aspect09, col=rainbow(100))

tpi09 <- terrain(b.zone9, 'TPI')
plot(tpi09)

depth09 <- b.zone9

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
