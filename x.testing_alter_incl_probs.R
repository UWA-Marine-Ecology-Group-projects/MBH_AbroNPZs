# test alter inclusion probs ###

#number of samples
n <- 10
#number of points to sample from
N <- 100^2
#the sampling grid (offset so that the edge locations have same area)
offsetX <- 1/(2*sqrt( N))
my.seq <- seq( from=offsetX, to=1-offsetX, length=sqrt(N))
X <- expand.grid( my.seq, my.seq)
#the legacy sites (three of them)
legacySites <- matrix( runif( 6), ncol=2, byrow=TRUE)
#names can be useful
colnames( X) <- colnames( legacySites) <- c("X1","X2")

#non-uniform inclusion probabilities
inclP <- 1-exp(-X[,1])
#scaling to enforce summation to n
inclP<- n * inclProbs / sum( inclP)
#uniform inclusion probabilities would be inclProbs <- rep( n/N, times=N)
#visualise
image( x=unique( X[,1]), y=unique( X[,2]),
       z=matrix( inclP, nrow=sqrt(nrow(X)), ncol=sqrt(nrow( X))),
       main="(Undadjusted) Inclusion Probabilities",
       ylab=colnames( X)[2], xlab=colnames( X)[1])
#The legacy locations
points( legacySites, pch=21, bg=grey(0.75), cex=1.5)

#alter inclusion probabilities
# so that new samples should be well-spaced from legacy
altInclProbs <- alterInclProbs( legacy.sites=legacySites,
                                potential.sites=X, inclusion.probs = inclP)
#visualise
image( x=unique( X[,1]), y=unique( X[,2]),
       z=matrix( altInclProbs, nrow=sqrt(nrow(X)), ncol=sqrt(nrow( X))),
       main="Adjusted Inclusion Probabilities",
       ylab=colnames( X)[2], xlab=colnames( X)[1])
#The legacy locations
points( legacySites, pch=21, bg=grey(0.75), cex=1.5)
