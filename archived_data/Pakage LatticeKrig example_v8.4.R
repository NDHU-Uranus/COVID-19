library(fields)
library(LatticeKrig)
# NOT RUN {
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
# thin plate spline-like model with the lambda parameter estimated by
# maximum likelihood. Default choices are made for a.wght, nlevel, NC
# and alpha.

obj<- LatticeKrig( x, y)
# }
# NOT RUN {
# summary of fit and a plot of fitted surface
print( obj)
surface( obj )
US(add=TRUE)
points(x)
# prediction standard errors
out.se<- predictSE( obj, xnew= x)
# predict at observations:
out.fhat<- predict( obj, xnew= x)
# conveniently predict on a 100X100 grid for plotting
out.surf<- predictSurface( obj, nx=100, ny=100)
# image.plot( out.surf) 
# }
# NOT RUN {
# running an example by first setting up the model object
# }
# NOT RUN {
# this is just a small model to run quickly
# compare the LKinfo object here  to one created implicitly:  obj$LKinfo
LKinfo1<- LKrigSetup( x, NC=5, nlevel=3, a.wght=4.1, nu=1.0)
obj1<- LatticeKrig( x,y, LKinfo= LKinfo1)
# }
# NOT RUN {
#
# In this example lon/lat are treated as just Euclidean coordinates 
# a quick adjustment for small regions is to account for the difference
# in physical distance in N-S verses E_W
# is to just scale the longitude degrees to be comparable to degrees in latitude
# at least in the middle of the domain. The assumption is that for small spatial
# domains this approximation will not be bad for the coordinates at the edges too.
# You accomplish this by adding a scaling, V matrix:
# Here the V argument is rolled into the LKinfo object created within the function
#
# }
# NOT RUN {
meanLat<- mean( x[,2])*pi/180
Vlonlat <- diag(  c( 1/cos(meanLat), 1) )
obj1<- LatticeKrig( x, y, V = Vlonlat )
# }
# NOT RUN {
# }
# NOT RUN {
# Refit using with just one level of  basis functions
# on a 20X20 grid within the spatial domain ( so about 400) 
# actually number is 720 ( see obj1b$LKinfo) due adding edge nodes
# Add an aspect ratio of spatial domain 
# and find the a.wght parameter along with nugget and process variances.
# this takes a while partly because LatticeKrig model is not optimized for small data sets!
obj1b<- LatticeKrig( x, y, nlevel=1, NC=20, findAwght=TRUE)
# rudimentary look at how likelihood was optimized
#log lambda and omega =  log(a.wght-4)/2 are useful parameterization ...
quilt.plot( obj1b$MLE$lnLike.eval[,c("logLambda","omega")],
            obj1b$MLE$lnLike.eval[,"lnProfileLike.FULL"], 
            xlab="loglamda", ylab="omega",
            zlim =c(-640,-612))
points( obj1b$MLE$lnLike.eval[,c("logLambda","omega")],cex=.25)

# }
# NOT RUN {
# fitting replicate spatial data sets
# here we use the common observations over days for the ozone
# data set. Whether these are true replicated fields is in question
# but the analysis is still useful

# }
# NOT RUN {
Y<-  na.omit( t( ozone2$y) ) 
ind<- attr( Y,"na.action")
X<- ozone2$lon.lat[-ind, ]

out1<- LatticeKrig( X, Y, nlevel=1, NC=20, findAwght=TRUE)
out2<- LatticeKrig( X, Y, nlevel=1, NC=20, findAwght=TRUE,
                    collapseFixedEffect=TRUE)
# compare the two models 
# Note second a.wght reflects more spatial correlation when individual 
# fixed effect is not removed ( 4.4 verses 4.07)
# nugget variance is nearly the same!
out1$MLE$summary[1:7]                        
out2$MLE$summary[1:7]



# }
# NOT RUN {
# Refit using the tensor product type of basis functions
# (default is "Radial"). An example how an additional argument that is 
# passed to the LKrigSetup function to create the LKinfo object.
obj2<- LatticeKrig( x, y, BasisType="Tensor")
# }
# NOT RUN {
#
# A 1-d example with 3 levels of basis functions
# See LKrig for an explanation if nlevel, NC,  alpha and a.wght 
# covariance parameters.


# }
# NOT RUN {
x<- matrix(rat.diet$t)
y<- rat.diet$trt
fitObj<- LatticeKrig( x, y)
# NOTE lots of defaults are set for the model! See print( fitObj)
plot( x,y)
xg<- matrix(seq( 0,105,,100))
lines( xg, predict(fitObj, xg) )
# }
# NOT RUN {
# }
# NOT RUN {
#  a 3D example
set.seed( 123)
N<- 1000
x<-  matrix( runif(3* N,-1,1), ncol=3, nrow=N)
y<-   10*exp( -rdist( x, rbind( c(.5,.5,.6) ) )/.5)

# NOTE setting of memory size for Cholesky. This avoids some warnings and
# extra computation by the spam package
LKinfo<- LKrigSetup( x,  nlevel=1,  a.wght= 6.01, NC=6, NC.buffer=2,
                     LKGeometry="LKBox", normalize=FALSE, mean.neighbor=200,
                     choleskyMemory=list(nnzR= 2E6) )                                      
out1<- LatticeKrig( x,y, LKinfo=LKinfo)

glist<- list( x1=seq( -1,1,,30), x2=seq( -1,1,,30), x3 = 0)
xgrid<- make.surface.grid( glist)

yhat<- predict( out1, xgrid)
# compare yhat to true function created above
image.plot( as.surface( glist, yhat))

# }
# NOT RUN {
#
###########################################################################
# Including a covariate (linear fixed part in spatial model)
########################################################################## 
# }
# NOT RUN {
data(COmonthlyMet)

obj  <- LatticeKrig(CO.loc,  CO.tmin.MAM.climate, Z=CO.elev)
obj2 <- LatticeKrig(CO.loc, CO.tmin.MAM.climate)

# compare with and without linear covariates
set.panel(1,2)
surface(obj)
US(add=TRUE)
title("With Elevation Covariate")

surface(obj2)
US(add=TRUE)
title("Without Elevation Covariate")

# }
# NOT RUN {
data(COmonthlyMet)
# Examining a few different "range" parameters
a.wghtGrid<-  4  +  c(.05, .1, .5, 1, 2, 4)^2

#NOTE smallest is "spline like" the largest is essentially independent
# coefficients at each level.  In this case the "independent" end is
# favored but the eff df. of the surface is very similar across models
# indicating about the same separate of the estimates into spatial
# signal and noise
#
for( k in 1:5 ){
  obj  <- LatticeKrig(CO.loc,  CO.tmin.MAM.climate, Z=CO.elev, 
                      a.wght=a.wghtGrid[k])
  cat( "a.wght:", a.wghtGrid[k], "ln Profile Like:",
       obj$lnProfileLike, "Eff df:", obj$trA.est, fill=TRUE)
}

# MLE
obj0  <- LatticeKrig(CO.loc,  CO.tmin.MAM.climate, Z=CO.elev, 
                     findAwght=TRUE)
print(obj0$MLE$summary)
# }
# NOT RUN {
#########################################################################
# Reproducing some of the analysis for the example in the
# JCGS LatticeKrig paper.
#########################################################################

#### Here is an example of dealing with approximate spherical geometry.
# }
# NOT RUN {
data(NorthAmericanRainfall)
library(mapproj)
x<- cbind(NorthAmericanRainfall$longitude, NorthAmericanRainfall$latitude)
y<- NorthAmericanRainfall$precip
log.y<- log(y)
elev<- NorthAmericanRainfall$elevation
# this is a simple projection as part of this and handled by the mapproj package
x.s<- mapproject( x[,1], x[,2], projection="stereographic")
x.s<- cbind( x.s$x, x.s$y)

# an alternative is to transform coordinates using another projection,
# e.g. a Lambert conformal projection
# with the project function from the rgdal package
# library( rgdal)
# x.s<- project(x,"+proj=lcc +lat_1=22 +lat_2=58 +lon_0=-93 +ellps=WGS84")
# this package has the advantage that the inverse projection is also 
# included ( inv=TRUE) so it is easy to evaluate the surface back on a Mercator grid.

obj0<- LatticeKrig(x.s, log.y, Z=elev )

fitSurface<- predictSurface( obj0, drop.Z=TRUE)
fitSurface$z<-  exp(fitSurface$z)/100
colorTable<- designer.colors( 256, c("red4", "orange", "yellow","green1", "green4", "blue"))
image.plot( fitSurface, col=colorTable)
map( "world", add=TRUE, col="grey30", lwd=3, proj="") 

# }
# NOT RUN {
# }