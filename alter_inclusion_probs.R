

#### Alter inclusion probabilities for legacy sites ----


# libraries ----
library( rgdal)
library( sp)
library( raster)
library( MBHdesign)

# clear environment ----
rm(list = ls())

# Directories ----
w.dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
#w.dir <- "~/MBH_AbroNPZs"
p.dir <- paste(w.dir, "plots", sep = '/')
d.dir <- paste(w.dir, "data", sep='/')
s.dir <- paste(w.dir, "shapefiles", sep='/')
r.dir <- paste(w.dir, "rasters", sep='/')


# Read legacy sites ----
legacy <- readOGR(paste(s.dir, "Legacy_sites.shp", sep='/'))
proj4string(legacy)
plot(legacy, pch=20, add=T)

legacys <- as.data.frame(legacy)
legacyss <- legacys[,c(6,7)]
# in ascending longitude
legacyss <- arrange(legacyss, coords.x1)
names(legacyss) <- c("x", "y")


# Read inclusion probabilities ----
inclProbs <- raster(paste(r.dir, "inclProbs_zone9.50_deployments.v2.tif", sep='/'))
plot(inclProbs)



# potential sites ----
pot.sites <- as.data.frame(inclProbs, xy = TRUE)
class(pot.sitess)

# as df
pot.sitess <- pot.sites[,c(1,2)]

# inclusion probs as numeric or df
ip <- pot.sites[,3]
class(ip)



# alter inclProbs test1----
altInclProbs <- alterInclProbs(legacy.sites = legacyss, 
                               potential.sites = pot.sitess,
                               #n = 40,
                               inclusion.probs = ip,
                               mc.cores = 6)




plot(altInclProbs)

#visualise
image( x=unique( pot.sitess[,1]), y=unique( pot.sitess[,2]),
       z=matrix( altInclProbs, nrow=sqrt(nrow(pot.sitess)), ncol=sqrt(nrow( pot.sitess))),
       main="Adjusted Inclusion Probabilities",
       ylab=colnames( pot.sitess)[2], xlab=colnames( pot.sitess)[1])

