library(sp)
library(ggplot2)
library(rgdal)
CA_case <- read.csv("D:\\研究所Meeting\\COVID19\\CA government data\\2022_4_18 CA case.csv")[-59,3:4]
coordinates(CA_case) <- ~ long + lat
# First, make a rectangular grid over your `SpatialPolygonsDataFrame`
grd <- makegrid(CA_case, n = 500)
colnames(grd) <- c("x", "y")

# Next, convert the grid to `SpatialPoints` and subset these points by the polygon.
grd_pts <- SpatialPoints(
  coords      = grd, 
  proj4string = CRS(proj4string(CA_case))
)

# subset all points in `grd_pts` that fall within `spdf`
grd_pts_in <- grd_pts[CA_case, ]

ggplot(as.data.frame(coordinates(grd_pts_in))) +
  geom_point(aes(x, y))




boroughs <- readOGR(dsn ="D:\\研究所Meeting\\COVID19\\R code\\CA_Counties( outline boundary)", layer = "CA_Counties_TIGER2016") 
boroughoutline <- fortify(boroughs, region = "COUNTYFP")
TA <- CRS("+proj=aea +lat_1=34 +lat_2=40.5 +lat_0=0 +lon_0=-120 +x_0=0 +y_0=-4000000 +datum=NAD83 +units=km +ellps=GRS80")
# Reproject to have the same CRS as TA.
ca <- spTransform(boroughs, TA)
plot(ca, border='gray')
points(grd_pts, cex=.5, col='red')


