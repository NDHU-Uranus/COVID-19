# import data and libaries 
library(ggplot2)
library(ggmap)
library(data.table)
#Yunan<-read.csv("C:\\Program Files\\RStudio\\data\\pb_sp\\sample_data.csv", header = TRUE)
CAdata<-setDT(read.csv("D:\\研究所Meeting\\COVID19\\Kriging predict\\predkrig.csv", header = TRUE))
register_google("AIzaSyBc5bqZv4d85TVa1fsFwOXBHbcHTY-AUXc")
# call the map to see point distribution
#map type="terrain"
california_map<-get_map(location="california",zoom=6,maptype="toner-lite",scale=2)
ggmap(california_map)+geom_point(data=CAdata,aes(x=CAdata$long,y=CAdata$lat,fill="red",alpha=0.3,size=0.05,shape=21))+scale_shape_identity()


# 1.  generate bins for x, y coordinates (unit=decimal degree)
xbreaks <- seq(floor(min(CAdata$lat,na.rm=TRUE)), ceiling(max(CAdata$lat,na.rm=TRUE)), by = 0.3)
ybreaks <- seq(floor(min(CAdata$long,na.rm=TRUE)), ceiling(max(CAdata$long,na.rm=TRUE)), by = 0.3)

# 2.  allocate the data points into the bins
CAdata$longbin <- ybreaks[cut(CAdata$long, breaks = ybreaks, labels=F)]
CAdata$latbin <- xbreaks[cut(CAdata$lat, breaks = xbreaks, labels=F)]

# 3.  summarise the data for each bin (use the median)
datamat <- CAdata[, list(simu_y_hat= median(simu_y_hat)), 
                 by = c("longbin","latbin")]

# 4. Merge the summarised data with all possible x, y coordinate combinations to get 
# a value for every bin
datamat <- merge(setDT(expand.grid(latbin = xbreaks, longbin = ybreaks)), datamat, 
                 by = c("latbin", "longbin"), all.x = TRUE, all.y = FALSE)

# 5. Fill up the empty bins 0 to smooth the contour plot
datamat[is.na(simu_y_hat), ]$simu_y_hat <- 0

my_palette <- colorRampPalette(c("blue", "white", "red"))
# 6. Plot the contours
ggmap(california_map,extent ="device") +
  stat_contour(data = datamat, aes(x = longbin, y = latbin, z = simu_y_hat, 
                                   fill = ..level.., alpha = ..level..), geom = 'polygon', binwidth = 0.01) +
  scale_fill_gradientn(name = "simu_y_hat", colors = my_palette(100)) +
  guides(alpha = FALSE)

