

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
class(legacyss)

# test with just one
legacyss <- legacyss[c(3,8),]


# Read inclusion probabilities ----
inclProbs <- raster(paste(r.dir, "inclProbs_zone9.50_deployments.v2.tif", sep='/'))
plot(inclProbs)



# potential sites ----
pot.sites <- as.data.frame(inclProbs, xy = TRUE)
head(pot.sites)
pot.sites <- arrange(pot.sites, y)
head(pot.sites)
class(pot.sitess)

# as df
pot.sitess <- pot.sites[,c(1,2)]

# inclusion probs as numeric or df
ip <- as.data.frame(pot.sites[,3])
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
       z=matrix( ip, nrow=100, ncol=43),
       main="Inclusion Probabilities (Undadjusted)",
       ylab=colnames( pot.sitess)[2], xlab=colnames( pot.sitess)[1])

image( x=unique( X[,1]), y=unique( X[,2]),
       z=matrix( inclProbs, nrow=sqrt(nrow(X)), ncol=sqrt(nrow( X))),
       main="Adjusted Inclusion Probabilities",
       ylab=colnames( X)[2], xlab=colnames( X)[1])

