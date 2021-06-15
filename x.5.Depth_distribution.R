### Calculate depth distribution of designs ###

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


## Set parameters of design ---
# site <- "Abrolhos"
# NPZ <- "npz6"
# design <- "40Bruvs"
# version <- "v2"


# Read bathy ----

b <- raster(paste(r.dir, "WA_500m_Bathy.tif", sep='/'))


# Read design shapefile ----

#d.shp <- readOGR(paste(o.dir, "Abrolhos-npz6-40Bruvs-v2.shp", sep='/'))
#d.shp <- readOGR(paste(o.dir, "Abrolhos-npz6-50Boss-v2.shp", sep='/'))
#d.shp <- readOGR(paste(o.dir, "Abrolhos-npz9-40Bruvs-v2.shp", sep='/'))
d.shp <- readOGR(paste(o.dir, "Abrolhos-npz9-50boss-v2.shp", sep='/'))


# Extract bathy for each point ----

point.depth <- raster::extract(b, d.shp, sp = T)
head(point.depth)

point.depth <- as.data.frame(point.depth)
head(point.depth)
max(point.depth$WA_500m_Bathy)
min(point.depth$WA_500m_Bathy)

# NPZ 6 ----

breaks <- c(-10,-20,-30,-40,-50,-60,-70,-80,-90,-100,-110,-120,-130,-140,-150,-160,-170,-180,-190,-200,-210,-220,-230,-250)
hist(point.depth$WA_500m_Bathy, breaks = breaks, freq = T, main = "10m intervals")

breaks <- c(-10,-30,-50,-70,-90,-110,-130,-150,-170,-190,-210)
hist(point.depth$WA_500m_Bathy, breaks = breaks)

breaks <- c(-10,-30,-70,-250)
hist(point.depth$WA_500m_Bathy, breaks = breaks, freq = T, main = "light classification intervals")


# NPZ 9 ----

breaks <- c(-80,-90,-100,-110,-120,-130,-140,-150)
hist(point.depth$WA_500m_Bathy, breaks = breaks, freq = T, main = "10m intervals")


