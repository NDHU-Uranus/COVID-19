library(MASS)
library(RobRSVD)
library(autoFRK)
library(rlmDataDriven)
#generate a t distrbution matrix
x<-matrix(rt(100, 1), nrow=10)

#generate the huber weight matrix with k=1.345
y=huberWeightLS(x, k=1.345)


library(yardstick)

# 創建一個 5x2 的矩陣
data <- as.data.frame(matrix(c(1,2,3,4,5,6,7,8,9,10), ncol = 2)) 

# 計算 Huber 損失
huber_loss(data[,1], data[,2], c = 1.5)
#--------------------------------------------
library(MASS)
CA_case <- read.csv("D:\\研究所Meeting\\COVID19\\CA government data\\2022_4_18 CA case.csv")[-59,]
#Confirmed rate(單日確診率)
CA_case$Confirmed_rate <- as.numeric(CA_case$reported_cases/CA_case$population)
# add label in data which used to general data with package "pdfos"
CA_case$Class <-rep(as.character(c('positive', 'negative')), times=29)
# select CA_case : long、lat、Confirmed_rate & label
CA_case <- CA_case[,c(3,4,16,17)]

simubas <- as.matrix(mrts(knot = CA_case[1:2], 10))

oc<- IRLS(X = simubas,Y = CA_case$Confirmed_rate,max.iter = 100,conv.eps = 1e-5)
oc[[1]]

#---------------package about huber m estimate 
vbnm <- model.matrix(CA_case$Confirmed_rate  ~ simubas)
mth <- CA_case$Confirmed_rate

#With power function as variance function
WHM_p <- whm(mth, vbnm[,-1], var.function = "power", tuning.para = 1.345)
#With exponential function as variance function
WHM_e <- whm(mth, vbnm[,-1], var.function = "exponential", tuning.para = 1.345)

#--------------------增加擾動項
library(mvtnorm)

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

