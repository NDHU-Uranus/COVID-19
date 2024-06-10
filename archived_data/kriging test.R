library(geoR)
library(ggplot2)
library(maptools)
library(rgdal)
library(dismo)
library(tmap)
califor_location <- read.csv("D:\\研究所Meeting\\COVID19\\ca location\\county location.csv")
califor_conf <- read.csv("D:\\研究所Meeting\\COVID19\\csse_covid_19_data\\California data\\statewide_cases(2021_2_27modify).csv")
califor_confirm <- califor_conf[c(-3,-4,-5)]
califor_confirm$confirm_rate <- califor_confirm[,2]/sum(califor_confirm[,2])

califor_location$county <- tolower(califor_location$county)
califor_confirm$county <- tolower(califor_confirm$county)
cali_confir_longlat <- merge(califor_location,califor_confirm,
                             by.x = 'county',
                             by.y = 'county')
head(cali_confir_longlat)

muldata <- cali_confir_longlat[5]
coor <-  cali_confir_longlat[c(2,3)]

# = = = = model fit (use (weighted)Least Square)= = = = #

# uvec : a vector with values used to define the variogram binning
vario.b <- variog(coords = coor, data = muldata,uvec=seq(0,12.25762, l=15)) #, max.dist=1)
# variogram cloud
vario.c <- variog(coords = coor, data = muldata, op="cloud")
#binned variogram and stores the cloud
vario.bc <- variog(coords = coor, data = muldata, bin.cloud=TRUE)

par(mfrow=c(2,2))
plot(vario.b, main="binned variogram") 
plot(vario.c, main="variogram cloud")
plot(vario.bc, bin.cloud=TRUE, main="clouds for binned variogram")  

ini.vals <- expand.grid(seq(0.005,0.1,l=20), seq(3,12,l=50)) # sigma^2 (partial sill) and phi (range parameter)
ols <- variofit(vario.b,fix.nug=TRUE, wei="equal", ini=ini.vals)
summary(ols)
lines(ols,lty=3,col=3,lwd=2)
# weight have "npairs" ,"cressie" ,"equal"
wls <- variofit(vario.b, ini=ini.vals, fix.nug=TRUE,wei="cressie")
summary(wls)
lines(wls, lty=2,col=4,lwd=2)

# = = = = = = = =simple kriging model = = = = = = = = = #
#convert this basic data frame into a spatial points data frame
coordinates(cali_confir_longlat) <- ~ long + lat
plot(cali_confir_longlat)

#establish the Coordinate Reference System (CRS).
proj4string(cali_confir_longlat) <- CRS('+proj=longlat +datum=NAD83')

#Then reproject to a more appropriate CRS, such as Teale Albers(亚尔勃斯投影). 
#Note the units=km, which is needed to fit the variogram.
TA <- CRS("+proj=aea +lat_1=34 +lat_2=40.5 +lat_0=0 +lon_0=-120 +x_0=0 +y_0=-4000000 +datum=NAD83 +units=km +ellps=GRS80")
aq <- spTransform(cali_confir_longlat, TA)

# create a template raster to interpolate to. 
# interpolate across California, so bring in the file counties.shp.
boroughs <- readOGR(dsn ="D:\\研究所Meeting\\COVID19\\R code\\CA_Counties( outline boundary)", layer = "CA_Counties_TIGER2016") 
boroughoutline <- fortify(boroughs, region = "COUNTYFP")

# Reproject to have the same CRS as TA.
ca <- spTransform(boroughs, TA)

# plot the points on CA 
plot(ca, border='gray')
points(aq, cex=.5, col='red')
#size of raster(點陣圖)
r <- raster(ca)
res(r) <- 10  # 10 km if your CRS's units are in km
#Go from raster to a SpatialGrid object
g <- as(r, 'SpatialGrid')
#variogram cloud
vcloud <- variogram(confirm_rate~1, locations=aq, width=20, cloud = TRUE)
plot(vcloud)
#sample variogram
sam_vari <- variogram(confirm_rate~1, locations=aq, width=20)
plot(sam_vari)
# fit model and give the original point
expfit_samviri <- fit.variogram(sam_vari, model = vgm(psill = 0.008, model = "Exp", range = 150, nugget = 0))
expfit_samviri 

plot(variogramLine(expfit_samviri, 400), type='l', ylim=c(0,0.005), col='blue', main = 'Exponential variogram model')
points(sam_vari[,2:3], pch=20, col='red')

ordikrig <- gstat(formula = confirm_rate~1, locations = aq, model=expfit_samviri)
# predict or interpolate for our grid g
pred_ordikrig <- predict(ordikrig, g)
# Convert kriged surface to a raster object for clipping(轉換為點陣圖)
rastpred_ordikrig <- raster(pred_ordikrig)
rastpred_ordikrig <- mask(rastpred_ordikrig, ca)

krigi_ca<-tm_shape(rastpred_ordikrig) + 
  tm_raster(n=10, palette="RdBu", auto.palette.mapping=FALSE, 
            title=" kriging in covid-19") +
  tm_legend(legend.outside=TRUE)
krigi_ca

#run 5-fold cross-validation to estimate the test prediction error
RMSE <- function(observed, predicted) {
  sqrt(mean((predicted - observed)^2, na.rm=TRUE))
}

set.seed(1234)
kf <- kfold(nrow(aq))
rmseok <- rep(NA, 5)
for (k in 1:5) {
  test <- aq[kf == k, ]
  train <- aq[kf != k, ]
  gscv <- gstat(formula = confirm_rate~1, locations = train, model=fve.o)
  p <- predict(gscv, newdata = test, debug.level=0)$var1.pred
  rmseok[k] <- RMSE(test$confirm_rate, p)
}
# 5-fold root mean squared error(RMSE)
mean(rmseok)

#------------another way to draw kriging ----------------#
## 1. Create a grid from the values in your points dataframe
## first get the range in data
x.range <- as.integer(range(cali_confir_longlat@coords[,1]))
y.range <- as.integer(range(cali_confir_longlat@coords[,2]))
##2. Create a grid with a slightly larger extent
plot(cali_confir_longlat)
#use the locator to click 4 points beyond the extent of the plot
#and use those to set your x and y extents
locator(4)
x.range <- as.integer(c(-123,-115))
y.range <- as.integer(c(33,41))
## now expand your range to a grid with spacing that you'd like to use in your interpolation
#here we will use 200m grid cells:
grd <- expand.grid(x=seq(from=x.range[1], to=x.range[2], by=200), y=seq(from=y.range[1], to=y.range[2], by=200))

## convert grid to SpatialPixel class
coordinates(grd) <- ~ x+y
gridded(grd) <- TRUE

## test it out - this is a good way of checking that your sample points are all well within your grid. If they are not, try some different values in you r x and y ranges:
plot(grd, cex=1.5)
points(cali_confir_longlat, pch=1, col='red', cex=1)
title("Interpolation Grid and Sample Points")


#sigma(sill), phi(range),tau(nugget)
model.variog<-vgm(psill=0.01, model="Exp", nugget=0, range=4.65)
krig<-krige(formula=confirm_rate ~ 1, locations=cali_confir_longlat, newdata=grd, model=model.variog)

krig.output=as.data.frame(krig)
names(krig.output)[1:3]<-c("long","lat","var1.pred")

plot<-ggplot(data=krig.output,aes(x=long,y=lat))#start with the base-plot and add the Kriged data to it
layer1<-c(geom_tile(data=krig.output,aes(fill=var1.pred)))#then create a tile layer and fill with predicted
layer2<-c(geom_path(data=boroughoutline,aes(long, lat, group=group),colour = "grey40", size=1))#then create an outline
plot+layer1+layer2+scale_fill_gradient(low="#FEEBE2", high="#7A0177")+coord_equal()

#----------------kriging by myself--------------------------#


#---------------Spatially Varying Coefficients(SVC) model---------------------------------#





