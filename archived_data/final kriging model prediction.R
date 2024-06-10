library(modelr)
library(imbalance)
library(autoFRK)
library(geoR)
library(MLmetrics)
library(stats)
library(mgcv) 
library(rlmDataDriven)
library(plot3D)
library(ggplot2)
library(plotly)
library(fields)
library(splines)
library(ie2misc) #MAE 
library(proxy) #calculate distance matrix
library(BBmisc) #normalization
library(MASS) #covaariance matrix
library(dplyr)#full_join
library(ggmap)
library(data.table)
#Confirmed rate(單日確診率/per 100000 peoples) and normalization(or standardization)
CA_case <- read.csv("D:\\研究所Meeting\\COVID19\\CA government data\\2022_4 time series CA case.csv")
vari_unemp <-read.csv("D:\\研究所Meeting\\COVID19\\Unemployment\\California_U(modify).csv ")
vari_pover <-read.csv("D:\\研究所Meeting\\COVID19\\Poverty\\california_P(modify).csv ")
vari_educa <-read.csv("D:\\研究所Meeting\\COVID19\\Education\\california_E(modify).csv ")
vari_bed <-read.csv("D:\\研究所Meeting\\COVID19\\Hospital beds\\Hospitals(california).csv ")
vari_popula <- read.csv("D:\\研究所Meeting\\COVID19\\populations\\California counties population(2020).csv")
#vari_houseunits <-read.csv("D:\\研究所Meeting\\COVID19\\housing units\\housing_units_-_single_multi_mobile_-_by_county_2018(California).csv ")
CA_pm2.5_total_out <- read.csv("D:\\研究所Meeting\\COVID19\\PM2.5 in CA\\CA_pm2.5_total_out.csv")
#CA_pm2.5_total_out$several_days <- normalize(as.numeric(apply(CA_pm2.5_total_out[,4:19], 1, mean)),method = "range", range = c(0, 1))
CA_pm2.5_total_out$several_days <- scale(as.numeric(apply(CA_pm2.5_total_out[,4:19], 1, mean)))
CA_pm2.5_total_KNNout <- CA_pm2.5_total_out[,c(1:3,20)]
# "apply(df, 1, mean)" 1 means Calculate the average of each row
CA_case$Confirmed_rate <- scale(as.numeric(apply(CA_case[,23:39], 1, mean)))      
#CA_case$Confirmed_rate <- scale(as.numeric(apply(CA_case[,23:39], 1, mean)))
#hospital beds rate(病床持有率)
CA_case$Hospital_beds_rate <- scale(as.numeric(vari_bed$BEDS/sum(vari_bed$BEDS)))
CA_case$Unemployment <- scale(vari_unemp$Unemployment.2019.)
CA_case$Poverty <- scale(vari_pover$poverty.percent.)
CA_case$Education <- scale(vari_educa$education.2015.2019.)
CA_case$population_rate <- scale(as.numeric(vari_popula[,2]/sum(vari_popula[,2])))
CA_case <- cbind(CA_case[,c(1,4,5,40:45)],CA_pm2.5_total_KNNout[,4])
CA_casee <- CA_case[,-1]
#summary(lm(data=CA_case,Confirmed_rate~Hospital_beds+Unemployment+Poverty+Education+CA_pm2.5_total_KNNout[, 4]))
#----------generate grid points in CA------------------#
library(sf)
library(raster)
shp <- getData('GADM', country = 'USA', level = 2) %>%
  subset(NAME_1 == 'California') %>% # 選擇加州（California）
  st_as_sf()

grid <- st_make_grid(shp, n = c(45, 45),
                     what = "centers",
                     square = TRUE) %>% st_intersection(shp)

ggplot() +
  geom_sf(data = shp) +
  geom_sf(data = grid, col = "navy")

# 提取經緯度
grid_coords <- st_coordinates(grid)
# 創建經緯度資料框
coords_df <- data.frame(long = grid_coords[, "X"],
                        lat = grid_coords[, "Y"])
# 匯出經緯度資料框到CSV檔案
write.csv(coords_df, "D:\\研究所Meeting\\COVID19\\CA grid points\\grid_coordinates(864).csv", row.names = TRUE)
plot(coords_df)
#some unknow point to predict
simu_grid_poin <-read.csv("D:\\研究所Meeting\\COVID19\\CA grid points\\grid_coordinates(6168).csv ")
#---------------plot sth-----------------------------
# show some information in map
states <- map_data("state")
ca_df <- states %>%
  filter(region == "california")
# california county map
counties <- map_data("county")
ca_county <- counties %>%
  filter(region == "california")

ca_base <- ggplot(data = ca_df, mapping = aes(x = long, y = lat, group = group)) + 
  coord_quickmap() + 
  geom_polygon(color = "black", fill = "gray") 
ca_base + theme_void()

ca_basemap<-ca_base + theme_void() + 
  geom_polygon(data = ca_county, fill = NA, color = "white") +
  geom_polygon(color = "black", fill = NA)  # get the state border back on top
ca_basemap

# add confirmed and merge with county location 
CA_case$subregion <- tolower(CA_case$subregion)
ca_fact <- full_join(ca_county, CA_case[,c(-2,-3)], by = "subregion")[-2978,]
#base dot map
ca_basemap+
  geom_point(aes(x = long, y = lat, group = NULL),data = CA_case)

#-----------------------------population---------------------------------------#
califor_population_rate_map <- ca_base + 
  geom_polygon(data =  ca_fact, aes(fill = population_rate), color = "white") +
  geom_polygon(color = 'pink', fill = NA) +
  scale_fill_gradientn(colours = rocket(5))  +
  labs(title = 'califor population rate') + theme(plot.title = element_text(hjust = 0.5))
califor_population_rate_map 
#-----------------------------hospital beds rate------------------------#
califor_hospital_beds_rate_map <- ca_base + 
  geom_polygon(data =  ca_fact, aes(fill = Hospital_beds_rate), color = "white") +
  geom_polygon(color = 'pink', fill = NA) +
  scale_fill_gradientn(colours = mako(5))  +
  labs(title = 'califor hospital beds rate') + theme(plot.title = element_text(hjust = 0.5))
califor_hospital_beds_rate_map 
#-----------------------------PM2.5------------------------#
califor_PM2.5_map <- ca_base + 
  geom_polygon(data =  ca_fact, aes(fill = `CA_pm2.5_total_KNNout[, 4]`), color = "white") +
  geom_polygon(color = 'pink', fill = NA) +
  scale_fill_gradientn(colours = cividis(5))  +
  labs(title = 'califor PM2.5') + theme(plot.title = element_text(hjust = 0.5))
califor_PM2.5_map 


#---------------------------------[star here]---------------------------
#--------------------------------------------
#-------------------weight Huber m estimate package(WHM) --------------
kriging.predict.WHM. <- function(fold,basis,simu_grid_poin,vari.dist){
 
    nn <- nrow(CA_casee)
    ### number of basis = 5
    # use simulation location in mrts to generate basis
    #bbaa <- 10 # number of basis
    simubas <- as.matrix(cbind(mrts(knot = CA_casee[1:2], basis),CA_casee$Hospital_beds
                               ,CA_casee$Unemployment,CA_casee$Poverty,CA_casee$Education
                               ,CA_casee$population_rate,CA_casee$`CA_pm2.5_total_KNNout[, 4]`))
    
    CA_case_Design.matrix <- model.matrix(CA_casee$Confirmed_rate  ~ simubas)
    response.variable <- CA_casee$Confirmed_rate
    #oc:outcome
    #ite : Number of iterations for the estimation procedure.
    #var.function = "power" or "exponential"
    WHM_p <- whm(response.variable, CA_case_Design.matrix[,-1], var.function = "power", tuning.para = 1.345, ite = 200)
    WHM_p[[1]]
    irls.simresidu <-WHM_p[[2]]
    
    #covariavce data = (observation data) - mean
    covadata<- CA_casee[,3] - simubas%*%WHM_p[[1]]
    
    # uvec : a vector with values used to define the variogram binning
    #variogram()
    #max.dist : length of variogram
    jj=max(abs(CA_casee[,1:2]))
    simuvario.b <- variog(coords = CA_casee[1:2], data = covadata,uvec=seq(0,jj, l=vari.dist), max.dist=6)
    #plot(simuvario.b, main="binned variogram") 
    
    simuini.vals <- expand.grid(seq(0.05,simuvario.b$var.mark,l=jj), seq(0.05,jj,l=jj)) # sigma^2 (partial sill) and phi (range parameter)
    simuols <- variofit(simuvario.b,fix.nug = F, wei="equal",cov.model="exp", ini=simuini.vals)
    #lines(simuols,lty=3,col=3,lwd=2)
    
    #weight = "npairs"、"equal"、"cressie"
    simuini.valss <- expand.grid(seq(0.05,simuols$cov.pars[1],l=jj), seq(0.05,jj,l=jj)) # sigma^2 (partial sill) and phi (range parameter)
    simuwls <- variofit(simuvario.b, ini=simuini.valss, fix.nug = F,wei="cressie",cov.model="exp")
    #lines(simuwls, lty=2,col=4,lwd=2)
    
    ###-------- use kriging model  (weight=wls) and see the model fit (MSE) -------#
    #simuweight = wei.select
    simuweight = simuwls
    # mean function residu
    # τ^2 + σ^2, corresponds to the variance of the observation process Y
    simresidu<-as.matrix(irls.simresidu)
    TPS_simu_beta_hat <- WHM_p[[1]]
    # simu_b_grid:grid point(simu_grid_poin) into TPS basis
    simu_b_grid <- as.matrix(cbind(mrts(knot = simu_grid_poin[,2:3], length(CA_casee))))
    simu_phi_s0 <- t(t(TPS_simu_beta_hat) %*% t(simu_b_grid))
    simu_c_s0 <- simuweight$cov.pars[1]*exp(-abs(as.data.frame.matrix(dist(simu_grid_poin[,2:3], CA_casee[,1:2] , method="euclidean")))/(simuweight$cov.pars[2]+1e-20))
    simu_zigma_theta <- simuweight$cov.pars[1]*as.matrix(exp(-abs(dist(CA_casee[,1:2], method = "euclidean", diag = F, upper = FALSE))/(simuweight$cov.pars[2]+1e-20))) #+1e-20: prevent simu_zigma_theta = 0 or Inf
    simu_y_hat =  simu_phi_s0 + as.matrix(simu_c_s0)  %*% ginv(simu_zigma_theta + (simuweight$nugget+simuweight$cov.pars[1])*diag(nn)) %*% simresidu
   
    # grid point & y_hat
    krig.pred.data <- cbind(simu_grid_poin[,2:3],simu_y_hat)
    
 
  return(list(krig.pred.data))
}
#---------------------------------------------
krig.pre <- kriging.predict.WHM.(fold=5,basis=3,simu_grid_poin=simu_grid_poin,vari.dist=600) 

# grid point & y_hat
#simu_krig_data <- cbind(simu_grid_poin[1:3],simu_y_hat)

# predict dot map 
predot<-krig.pre[[1]]%>% as.data.frame %>% 
  ggplot(aes(long, lat)) + geom_point(aes(colour = simu_y_hat)) + 
  labs(title="pre dot") +
  scale_color_viridis(option = "H")
predot
#predict map
#first add "group"in simu_krig_data
# group<-as.data.frame(seq(from=1, to=length(simu_grid_poin$NAME))) 
# colnames(group)[1]  <- "group"
# simu_krig_data <- cbind(simu_krig_data,group)
# predict_map <- ca_base + 
#   geom_point(data =  simu_krig_data, aes(x=long, y=lat,colour = simu_y_hat)) +
#   geom_polygon(color = 'pink', fill = NA) +
#   scale_color_viridis(option = "H") +
#   labs(title = 'califor cumulative confirmed case') + theme(plot.title = element_text(hjust = 0.5))
# predict_map
outputKrig <- write.csv(krig.pre[[1]], "D:\\研究所Meeting\\COVID19\\Kriging predict\\predkrig.csv", row.names=FALSE)

#predict dot plot
# ca_basemap+
#     geom_point(aes(x = long, y = lat, group = NULL),data = simu_krig_data)

# 3D plot in estimate value---1
#scatter3D(krig.pre$long, krig.pre$lat, krig.pre$simu_y_hat,zlab="value",phi = 10, theta = 10,main="estimate value Plot")

##--------------------plot on google map----------------------------------#
CAdata<-setDT(read.csv("D:\\研究所Meeting\\COVID19\\Kriging predict\\predkrig.csv", header = TRUE))
register_google("AIzaSyBc5bqZv4d85TVa1fsFwOXBHbcHTY-AUXc")
# call the map to see point distribution
#map type="terrain"
california_map<-get_map(location="california",zoom=6,maptype="toner-lite",scale=2)
ggmap(california_map)+geom_point(data=CAdata,aes(x=CAdata$long,y=CAdata$lat,fill="red",alpha=0.3,size=0.05,shape=21))+scale_shape_identity()


# 1.  generate bins for x, y coordinates (unit=decimal degree)
xbreaks <- seq(floor(min(CAdata$lat,na.rm=TRUE)), ceiling(max(CAdata$lat,na.rm=TRUE)), by = 0.4)
ybreaks <- seq(floor(min(CAdata$long,na.rm=TRUE)), ceiling(max(CAdata$long,na.rm=TRUE)), by = 0.4)

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
                                   fill = ..level.., alpha = ..level..), geom = 'polygon', binwidth = 0.2) +
  scale_fill_gradientn(name = "simu_y_hat", colors = my_palette(100)) +
  guides(alpha = FALSE)


# -----------------------------------test-------------------------------------#
{
# library(ggplot2)
# 
# # 假設kriging的預測結果存儲在data中，包含座標（x、y）和預測值（pred）
# # data可以是一個data.frame，包含x、y和pred列
# 
# # 設定熱度圖顏色調色板
# my_palette <- colorRampPalette(c("blue", "white", "red"))
# 
# # 使用ggplot2繪製熱度圖
# ggplot(krig.pre[[1]], aes(long, lat, fill = simu_y_hat)) +
#   geom_tile(width = 1, height = 1) +
#   scale_color_viridis(option = "H")+
#   #scale_fill_gradientn(colors = my_palette(100)) +
#   theme_minimal()
# 
# library(ggplot2)
# 
# # 假設你有一個data.frame或資料集包含了座標和值的數據，命名為data
# # data中應該包含x座標、y座標和值（例如z）
# 
# # 使用geom_tile函數繪製熱度圖
# ggplot(krig.pre[[1]], aes(x = long, y = lat, fill = simu_y_hat)) +
#   geom_bin2d() +
#   scale_fill_gradient(low = "blue", high = "red") +
#   theme_gray()
# 
# aabba <- krig.pre[[1]]
# library(sp)
# coordinates(aabba) <- ~long+lat
# proj4string(aabba) <- CRS("+init=epsg:28992")
# 
# spplot(aabba, "simu_y_hat")+
#   layer(panel.points(long, lat, col="green", pch=19), data=aabba)
  }
