library(fields)
library(LatticeKrig)
# Load ozone data set
data(ozone2)  
x<-ozone2$lon.lat
y<- ozone2$y[16,]
# Find location that are not 'NA'. 
# LKrig is not set up to handle missing observations.   
good <-  !is.na(y)
x<- x[good,]
y<- y[good]
# Visualize the ozone data set
quilt.plot(x,y)
US(add=TRUE)

# Predict ozone concentrations over the domain using a single level lattice 
# with a maximum of 20 grid points in x and y direction.
obj<- LKrig(x,y,NC=30, lambda=.01,  a.wght=5)
# Plot fitted surface.
surface(obj) # see also predict.surface to just evaluate on a grid
US(add=TRUE)
# Check effective degrees of freedom in surface.
obj$trA.est

# Search over lambda and a.wght to find MLE.
# (these ranges have been selected to give a nice result)
NC<- 30
a.wght<- 4+ 1/(seq(1,5,,8))**2
lambda<- 10**seq( -6,-2,,8)*NC**2
alpha<- 1
out<- MLE.LKrig( x=x,y=y,NC=NC,alpha=alpha,a.wght=a.wght, lambda=lambda)
# Interpolate these values with Tps and find MLE
temp<-Tps(cbind( log10(out[,3]), 1/sqrt(out[,2]-4)) , out$lnProfileLike,lambda=0)
out.p<- predict.surface( temp)
# Plot ln profile likelihood  
image.plot( out.p, ylab="equivalent lattice range: 1/sqrt(a-4)",
            xlab="log10 lambda")
points( temp$x, pch="o", col="grey50")
title("ln Profike Likelihood")
# 'max.image' is a handy little function to find maximum in an image object
max.image<- function( obj){
  ind<- which.max( obj$z)
  ix<- row(obj$z)[ind]
  iy<- col(obj$z)[ind]
  list(x= obj$x[ix], y=obj$y[iy], z= obj$z[ind], ind= cbind( ix,iy))}
# Apply 'max.image' to find MLE over lambda and a.wght grid
max.out<- max.image(out.p)
points( max.out$x, max.out$y, pch=1, col="magenta")
contour( out.p, level= max.out$z - 2, add=TRUE, col="magenta", lwd=2)
# Refit surface using MLE of lambda and a.wght
objMLE<- LKrig( x,y,NC=NC, lambda=10**max.out$x, alpha=1, a.wght= 4 + 1/max.out$y^2 )
set.panel(2,1)
surface(obj,zlim=c(0,220)) 
title("Original Fit")
US( add=TRUE)
surface(objMLE,zlim=c(0,220)) 
title("MLE Fit")
US( add=TRUE)
# Note how the MLE fit is much smoother due to the larger value of lambda 

## Another example, this time with a covariate:
data(COmonthlyMet)
y<- CO.tmin.MAM.climate
good<- !is.na( y)
y<-y[good]
x<- as.matrix(CO.loc[good,])
Z<- CO.elev[good]
out<- LKrig( x,y,Z=Z,NC=30, lambda=.1, a.wght=5)
set.panel(2,1)
# Plot data
quilt.plot(x,y)
US(add=TRUE)
title("Minimum spring temperature in Colorado")
# Plot elevation
quilt.plot(x,Z)
US(add=TRUE)
title("Elevation at observation locations")
# Plot predictions with linear elevation term included
quilt.plot( x, predict(out),zlim=c(-12,12))
title("Predicted temperature with elevation term")
US(add=TRUE)
# Plot predictions without linear elevation term included
quilt.plot( x, predict(out, drop.Z=TRUE),zlim=c(-12,12))
title("Predicted temperature without elevation term")
US(add=TRUE)
# Note the larger gradients when elevation is included  
set.panel()

# A bigger problem: fitting takes about 30 seconds on fast laptop
data(CO2)
obj1<- LKrig( CO2$lon.lat,CO2$y,NC=100, lambda=5, a.wght=5)
# 4600 basis functions 100X46 lattice.
obj1$trA.est # about 1040 effective degrees of freedom 
#
glist<- list( x= seq( -180,180,,200),y=seq( -80,80,,100) )
xg<-  make.surface.grid(glist)
fhat<- predict( obj1, xg)
fhat <- matrix( fhat,200,100) # convert to image
#Plot data and gap-filled estimate
set.panel(2,1)
quilt.plot(CO2$lon.lat,CO2$y,zlim=c(373,381))
world(add=TRUE,col="magenta")
title("Simulated CO2 satellite observations")
image.plot( glist$x, glist$y, fhat,zlim=c(373,381))
world( add=TRUE, col="magenta")
title("Gap-filled global predictions")  

set.panel()

# Here is an illustration of using the fields function mKrig with this covariance
# to reproduce the computations of LKrig. The difference is that mKrig can
# not take advantage of any sparsity in the precision matrix 
# as the inverse covariance matrix may not be sparse.
# But this example reinforces the concept that LKrig finds the
# the standard geostatistical estiamte but just uses a
# particular covariance function
#
a.wght<- 5
lambda <-  1.5
obj1<- LKrig( x,y,NC=16, lambda=lambda, a.wght=5, NtrA=20,iseed=122)

# in both calls iseed is set the same so we can compare 
# Monte Carlo estimates of effective degrees of freedom
obj2<- mKrig( x,y, lambda=lambda, m=2, cov.function="LKrig.cov",
              cov.args=list( LKinfo=obj1$LKinfo), NtrA=20, iseed=122)
# The covariance "parameters" are all in the list LKinfo
# to create this outside of a call to LKrig use
LKinfo.test <- LKrig.setup( x, NC=16,  a.wght=5)

# compare the two results this is also an
# example of how tests are automated in fields
# set flag to have tests print results
test.for.zero.flag<- TRUE
test.for.zero( obj1$fitted.values, obj2$fitted.values,
               tag="comparing predicted values LKrig and mKrig")
test.for.zero( unlist(LKinfo.test), unlist(obj1$LKinfo),
               tag="comparing two ways of creating covariance list")