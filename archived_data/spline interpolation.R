library(mvtnorm)
CA_case <- read.csv("D:\\研究所Meeting\\COVID19\\CA government data\\2022_4_18 CA case.csv")[-59,]
#Confirmed rate(單日確診率)
CA_case$Confirmed_rate <- as.numeric(CA_case$reported_cases/CA_case$population)
# add label in data which used to general data with package "pdfos"
CA_case$Class <-rep(as.character(c('positive', 'negative')), times=29)
# select CA_case : long、lat、Confirmed_rate & label
CA_case <- CA_case[,c(3,4,16)]
# 原始數據
height <- c(170, 172, 168, 175, 180)
weight <- c(70, 75, 68, 80, 85)

# 計算平均值和標準差
height_mean <- mean(height)
height_sd <- sd(height)
weight_mean <- mean(weight)
weight_sd <- sd(weight)

# 設計擾動項
n <- length(height)
set.seed(123)
height_perturb <- rnorm(n, mean = 0, sd = height_sd * 0.1)
weight_perturb <- rnorm(n, mean = 0, sd = weight_sd * 0.1)

# 擾動後的數據
height_perturbed <- height + height_perturb
weight_perturbed <- weight + weight_perturb

# 印出擾動後的數據
print(height_perturbed)
print(weight_perturbed)


#---------------------
# 创建一个简单的数据集
x <- seq(0, 10, length.out = 50)
y <- sin(x)
data <- data.frame(x, y)

# 对 x 和 y 进行样条插值，生成 100 个新的观测点
new_x <- seq(0, 10, length.out = 100)
new_y <- spline(x, y, xout = new_x)$y
new_data <- data.frame(new_x, new_y)b


#-----------------------(主要研究)--------------
library(splines)

# 對經度和緯度進行樣條內插
new_longitude <- with(CA_case,spline(CA_case$long))
new_latitude <- with(CA_case, spline(CA_case$lat))

# 創建新的數據框，包含新的觀測點的經度和緯度
new_data <- data.frame(
  longitude = new_longitude$y, # 新的觀測點的經度
  latitude = new_latitude$y   # 新的觀測點的緯度
)
xy <- data.frame(longitude = CA_case$long, latitude = CA_case$lat)

# 使用smooth.spline对观测值进行平滑
fit <- with(CA_case, smooth.spline(x = xy, y = Confirmed_rate, spar = 0.5))

# 對新的觀測點進行預測，生成新的觀測值
new_data$value <- predict(fit,y= CA_case$Confirmed_rate, newdata = newdata)


# 對新的觀測點進行預測，生成新的觀測值
new_data$value <- predict(with(new_data, smooth.spline(longitude, latitude, value)), newdata = new_data[, c("longitude", "latitude")])

#---------------------------

# 使用predict函数预测新的观测值
fit <- with(CA_case, spline(x=cbind(long,lat),y=Confirmed_rate))
new_data$Confirmed_rate <- predict(fit, newdata = new_data)

# 输出新的数据框
print(new_data)

#--------------------------------
library(splines)

# 生成50個新經度和50個新緯度
new_longitude <- with(CA_case, seq(min(long), max(long), length.out = 50))
new_latitude <- with(CA_case, seq(min(lat), max(lat), length.out = 50))

# 將新經度和新緯度進行網格化組合
new_grid <- expand.grid(longitude = new_longitude, latitude = new_latitude)

# 對新的觀測點進行預測，生成新的觀測值
fit <- with(CA_case, smooth.spline(x = cbind(long, lat), y = Confirmed_rate))
new_grid$value <- predict(fit, newdata = new_grid[, c("longitude", "latitude")])$y

# 新數據集
new_data <- data.frame(longitude = new_grid$longitude, 
                       latitude = new_grid$latitude, 
                       value = new_grid$value)
#--------------------------------
library(splines)

# 将经度和纬度合并为一个矩阵
xy <- cbind(CA_case$long, CA_case$lat)

# 使用smooth.spline对观测值进行平滑
fit <- smooth.spline(x = xy[,1], y = CA_case$Confirmed_rate, spar = 0.5)

# 生成新的观测点
new_longitude <- seq(min(CA_case$long), max(CA_case$long), length.out = 50)
new_latitude <- seq(min(CA_case$lat), max(CA_case$lat), length.out = 50)
new_data <- expand.grid(longitude = new_longitude, latitude = new_latitude)

# 对新的观测点进行预测，生成新的观测值
new_data$Confirmed_rate <- predict(fit, newdata = new_data)

#---------------TPS generate-----------
# 載入fields套件
library(fields)
library(splines)
library(mvtnorm)
CA_case <- read.csv("D:\\研究所Meeting\\COVID19\\CA government data\\2022_4_18 CA case.csv")[-59,]
#Confirmed rate(單日確診率)
CA_case$Confirmed_rate <- as.numeric(CA_case$reported_cases/CA_case$population)*100000
# select CA_case : long、lat、Confirmed_rate & label
CA_case <- CA_case[,c(3,4,16)]
# 對經度和緯度進行樣條內插
new_longitude <- with(CA_case,spline(CA_case$long, n = 100))
new_latitude <- with(CA_case, spline(CA_case$lat, n = 100))

# 創建新的數據框，包含新的觀測點的經度和緯度
new_data <- data.frame(
  longitude = new_longitude$y, # 新的觀測點的經度
  latitude = new_latitude$y   # 新的觀測點的緯度
)
# 建立TPS模型
tps_model <- Tps(CA_case[,1:2], CA_case[,3])

# 預測新的50個經緯度和觀測值

new_data$obs <- predict(tps_model, new_data[,1:2])
