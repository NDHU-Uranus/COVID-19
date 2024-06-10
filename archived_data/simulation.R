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

#-------------- generate grid points in CA------------------------#
library(sf)
library(raster)
library(ggplot2)
shp <- getData('GADM', country = 'USA', level = 2) %>%
  subset(NAME_1 == 'California') %>% # 選擇加州（California）
  st_as_sf()

grid <- st_make_grid(shp, n = c(100, 100),
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
write.csv(coords_df, "D:\\研究所Meeting\\COVID19\\CA grid points\\grid_coordinates(5179).csv", row.names = TRUE)
plot(coords_df)
#--------------------------------------
#some unknow point to predict
simu_grid_poin <-read.csv("D:\\研究所Meeting\\COVID19\\CA grid points\\grid_coordinates(864).csv ")
#---------------------------------[star here]---------------------------
#-----------generate data with multivariate normal distribution------------------------------#
##--------------(((((111111))))))------------------*
library(mvtnorm)
library(Matrix)
mul <- function(n, jj, extreme_prob, extreme_range) {
  # location
  long <- runif(n, max = jj, min = 0)
  lat <- runif(n, max = jj, min = 0)
  coords <- data.frame(long, lat)
  # distance matrix
  dist_matrix <- as.matrix(dist(coords))
  # covariance matrix
  sigma <- exp(-dist_matrix)
  # multivariate normal
  random_data <- t(rmvnorm(1, mean = rep(0, n), sigma = sigma))
  # Calculate the number of extreme values based on extreme_prob
  num_extreme <- floor(n * extreme_prob)
  # Determine the positions to replace with extreme values
  extreme_positions <- sample(n, num_extreme)
  # Generate extreme values
  #extreme_vals <- rnorm(num_extreme, mean = 10, sd = 5)
  extreme_vals <- runif(num_extreme, min = extreme_range[1], max = extreme_range[2])
  # Replace the values at extreme positions with extreme values
  random_data[extreme_positions] <- extreme_vals
  # return
  return(list(location = coords, data = random_data))
}
# generate data with extreme values
# Add extreme values with a probability of 0.05
#result <- mul(n = 800, jj = 100, extreme_prob = 0.05) 
result <- mul(n = 800, jj = 100, extreme_prob = 0.00, extreme_range = c(-8,8))
CA_case <-cbind(as.data.frame(result$location),as.data.frame(result$data))
names(CA_case)[3] <- "value"
plot(CA_case$value)
#--------------------------------------------
#-------------------weight Huber m estimate package(WHM) --------------
k_folds_WHM.cv <- function(fold,bbaa,vari.dist){
  k_fold_mse <-c()
  k_fold_mae <-c()
  #pdfos: Probability density function estimation based oversampling
  #set.seed(6)
  #new_data <- pdfos(data = CA_case, numInstances = datanum,classAttr = "Class")
  cv  <- crossv_kfold(data = CA_case, k = fold ) # k-fold cross-Validation
  for (t in 1:length(cv$.id)) {
    #CA_case[,-4] is to elimate label
    case_train <- CA_case[cv$train[[t]]$idx,]
    case_test <- CA_case[cv$test[[t]]$idx,]
    nn <- length(case_train$long)
    ### number of basis = 5
    # use simulation location in mrts to generate basis
    #bbaa <- 10 # number of basis
    simubas <- as.matrix(mrts(knot = case_train[1:2], bbaa))
    
    CA_case_Design.matrix <- model.matrix(case_train$value  ~ simubas)
    response.variable <- case_train$value
    #oc:outcome
    #ite : Number of iterations for the estimation procedure.
    #var.function = "power" or "exponential"
    WHM_p <- whm(response.variable, CA_case_Design.matrix[,-1], var.function = "power", tuning.para = 1.345, ite = 500)
    WHM_p[[1]]
    irls.simresidu <-WHM_p[[2]]
    
    #covariavce data = (observation data) - mean
    covadata<- case_train[,3] - simubas%*%WHM_p[[1]]
    
    # uvec : a vector with values used to define the variogram binning
    #variogram()
    #max.dist : length of variogram
    jj=max(abs(case_train[,1:2]))
    simuvario.b <- variog(coords = case_train[1:2], data = covadata,uvec=seq(0,jj, l=vari.dist), max.dist=6)
    #plot(simuvario.b, main="binned variogram") 
    
    simuini.vals <- expand.grid(seq(0.01,simuvario.b$var.mark,l=jj), seq(0.1,jj,l=jj)) # sigma^2 (partial sill) and phi (range parameter)
    simuols <- variofit(simuvario.b,fix.nug=T, wei="equal",cov.model="exp", ini=simuini.vals)
    #lines(simuols,lty=3,col=3,lwd=2)
    
    simuini.valss <- expand.grid(seq(0.01,simuols$cov.pars[1],l=jj), seq(0.1,jj,l=jj)) # sigma^2 (partial sill) and phi (range parameter)
    simuwls <- variofit(simuvario.b, ini=simuini.valss, fix.nug=F,wei="cressie",cov.model="exp")
    #lines(simuwls, lty=2,col=4,lwd=2)
    
    ###-------- use kriging model  (weight=wls) and see the model fit (MSE) -------#
    #simuweight = wei.select
    simuweight = simuwls
    # mean function residu
    # τ^2 + σ^2, corresponds to the variance of the observation process Y
    simresidu<-as.matrix(irls.simresidu)
    TPS_simu_beta_hat <- WHM_p[[1]]
    # simu_b_grid:grid point(simu_grid_poin) into TPS basis
    simu_b_grid <- as.matrix(mrts(knot = case_test[1:2], bbaa))
    simu_phi_s0 <- t(t(TPS_simu_beta_hat) %*% t(simu_b_grid))
    simu_c_s0 <- simuweight$cov.pars[1]*exp(-abs(as.data.frame.matrix(dist(case_test, case_train , method="euclidean")))/(simuweight$cov.pars[2]))
    simu_zigma_theta <- simuweight$cov.pars[1]*as.matrix(exp(-abs(dist(case_train[-3], method = "euclidean", diag = F, upper = FALSE))/(simuweight$cov.pars[2]))) #+1e-20: prevent simu_zigma_theta = 0 or Inf
    simu_y_hat =  simu_phi_s0 + as.matrix(simu_c_s0)  %*% ginv(simu_zigma_theta + (simuweight$nugget+simuweight$cov.pars[1])*diag(nn)) %*% simresidu
    
    # simu_c_s <-c()
    # for (i in 1:length(case_train[,c("long")])){
    #   for (j in 1:length(case_test[,c("long")])){
    #     # if(i==j){next}
    #     #else{
    #     simu_disten<-sqrt((case_train[,c("long")][i]-case_test[,c("long")][j])^2+(case_train[,c("lat")][i]-case_test[,c("lat")][j])^2) 
    #     simu_c_s <- c(simu_c_s ,simu_disten)
    #     # }
    #   }
    #   if(i%%100 == 0) cat("running: ", i, "\n")
    # }
    # simu_c_s0 <- simuweight$cov.pars[1]*exp(-abs(matrix(simu_c_s,nrow = nrow(case_test))/(simuweight$cov.pars[2]+1e-20)))
    # simu_zigma_theta <- simuweight$cov.pars[1]*as.matrix(exp(-dist(case_train[-3], method = "euclidean", diag = F, upper = FALSE, p = 2))/(simuweight$cov.pars[2]+1e-20)) #+1e-20: prevent simu_zigma_theta = 0 or Inf
    # simu_y_hat =  simu_phi_s0 + simu_c_s0  %*% solve(simu_zigma_theta + (simuweight$nugget+simuweight$cov.pars[1])*diag(nn)) %*% simresidu
    
    # grid point & y_hat
    simu_krig_data <- cbind(case_test,simu_y_hat)
    
    # check MSE & resudual
    simu_MSE = MSE(case_test$value, simu_krig_data$simu_y_hat)
    k_fold_mse <- c(k_fold_mse ,simu_MSE)
    simu_MAE = mae(case_test$value, simu_krig_data$simu_y_hat)
    k_fold_mae <- c(k_fold_mae ,simu_MAE) 
  }
  k_fold_mMSE <- mean(k_fold_mse)
  k_fold_RMSE <- sqrt(k_fold_mse)
  k_fold_mRMSE <- mean(k_fold_RMSE)
  k_fold_mMAE <- mean(k_fold_mae)
  return(list(k_fold_mse,k_fold_mMSE,k_fold_RMSE,k_fold_mRMSE,k_fold_mae,k_fold_mMAE,simu_krig_data))
}

#-------------------------------------------------
#WHM With power function as variance function
n.times_tempWHM.kfoldcv <-function(repeat.time,fold,basis,vari.dist){
  n.times_MSE.WHM.kfoldcv <-c()
  n.times_RMSE.WHM.kfoldcv <-c()  
  n.times_MAE.WHM.kfoldcv <-c()    
  for (i in 1:repeat.time) {
    outWHM <- k_folds_WHM.cv(fold = fold,bbaa = basis,vari.dist=vari.dist)
    n.times_MSE.WHM.kfoldcv<- c(n.times_MSE.WHM.kfoldcv ,outWHM[[2]])
    n.times_RMSE.WHM.kfoldcv<- c(n.times_RMSE.WHM.kfoldcv ,outWHM[[4]])
    n.times_MAE.WHM.kfoldcv<- c(n.times_MAE.WHM.kfoldcv ,outWHM[[6]])
  }
  return(list(n.times_MSE.WHM.kfoldcv, n.times_RMSE.WHM.kfoldcv,n.times_MAE.WHM.kfoldcv))
}

# reapet v times 
system.time({n.times_WHM.kfoldcv <- n.times_tempWHM.kfoldcv(repeat.time=10,fold=5,basis=3,vari.dist=600)})

par(mfrow = c(1, 3))
# k fold mse
plot(n.times_WHM.kfoldcv[[1]], ylab="MSE",
     type = "b",col = "turquoise", pch = 20,main = "rep 10 times MSE")
# k fold RMSE
plot(n.times_WHM.kfoldcv[[2]], ylab="RMSE",
     type = "b",col = "turquoise", pch = 20,main = "rep 10 timesd RMSE")
# k fold MAE
plot(n.times_WHM.kfoldcv[[3]], ylab="MAE",
     type = "b",col = "turquoise", pch = 20,main = "rep 10 times MAE")
mean.whm.MSE=mean(n.times_WHM.kfoldcv[[1]])
mean.whm.RMSE=mean(n.times_WHM.kfoldcv[[2]])
mean.whm.MAE=mean(n.times_WHM.kfoldcv[[3]])

#------------------------------------------------------------------------------
#---------------  k-flods cross-validation(LM in mean not huber) (with original data)--------
k_folds_lm.ori.cv <- function(fold,bbaa,vari.dist){
  k_fold_mse <-c()
  k_fold_mae <-c()
  #set.seed(6)
  cv  <- crossv_kfold(data = CA_case, k = fold ) # k-fold cross-Validation
  for (t in 1:length(cv$.id)) {
    #CA_case[,-4] is to elimate label
    case_train <- CA_case[cv$train[[t]]$idx,]
    case_test <- CA_case[cv$test[[t]]$idx,]
    nn <- length(case_train$long)
    ### number of basis = 5
    # use simulation location in mrts to generate basis
    #bbaa <- 10 # number of basis
    simubas <- as.matrix(mrts(knot = case_train[1:2], bbaa))
    
    #oc:outcome
    lmmodel<- lm(case_train$value ~ simubas-1)
    lmmodel[[1]]
    irls.simresidu <-lmmodel[[2]]
    
    #covariavce data = (observation data) - mean
    covadata<- case_train[,3] - simubas%*%lmmodel[[1]]
    
    # uvec : a vector with values used to define the variogram binning
    #variogram()
    #max.dist : length of variogram
    jj=max(abs(case_train[,1:2]))
    simuvario.b <- variog(coords = case_train[1:2], data = covadata,uvec=seq(0,jj, l=vari.dist), max.dist=6)
    #plot(simuvario.b, main="binned variogram") 
    
    simuini.vals <- expand.grid(seq(0.01,simuvario.b$var.mark,l=jj), seq(0.1,jj,l=jj)) # sigma^2 (partial sill) and phi (range parameter)
    simuols <- variofit(simuvario.b,fix.nug=T, wei="equal",cov.model="exp", ini=simuini.vals)
    #lines(simuols,lty=3,col=3,lwd=2)
    
    simuini.valss <- expand.grid(seq(0.01,simuols$cov.pars[1],l=jj), seq(0.1,jj,l=jj)) # sigma^2 (partial sill) and phi (range parameter)
    simuwls <- variofit(simuvario.b, ini=simuini.valss, fix.nug=F,wei="cressie",cov.model="exp")
    #lines(simuwls, lty=2,col=4,lwd=2)
    
    ###-------- use kriging model  (weight=wls) and see the model fit (MSE) -------#
    #simuweight = wei.select
    simuweight = simuwls
    # mean function residu
    # τ^2 + σ^2, corresponds to the variance of the observation process Y
    simresidu<-as.matrix(irls.simresidu)
    TPS_simu_beta_hat <- lmmodel[[1]]
    # simu_b_grid:grid point(simu_grid_poin) into TPS basis
    simu_b_grid <- as.matrix(mrts(knot = case_test[1:2], bbaa))
    simu_phi_s0 <- t(t(TPS_simu_beta_hat) %*% t(simu_b_grid))
    simu_c_s0 <- simuweight$cov.pars[1]*exp(-abs(as.data.frame.matrix(dist(case_test, case_train , method="euclidean")))/(simuweight$cov.pars[2]+1e-20))
    simu_zigma_theta <- simuweight$cov.pars[1]*as.matrix(exp(-abs(dist(case_train[-3], method = "euclidean", diag = F, upper = FALSE))/(simuweight$cov.pars[2]+1e-20))) #+1e-20: prevent simu_zigma_theta = 0 or Inf
    simu_y_hat =  simu_phi_s0 + as.matrix(simu_c_s0)  %*% ginv(simu_zigma_theta + (simuweight$nugget+simuweight$cov.pars[1])*diag(nn)) %*% simresidu
    
    # grid point & y_hat
    simu_krig_data <- cbind(case_test,simu_y_hat)
    
    # check MSE & resudual
    simu_MSE = MSE(case_test$value, simu_krig_data$simu_y_hat)
    k_fold_mse <- c(k_fold_mse ,simu_MSE)
    simu_MAE = mae(case_test$value, simu_krig_data$simu_y_hat)
    k_fold_mae <- c(k_fold_mae ,simu_MAE)  
  }
  k_fold_mMSE <- mean(k_fold_mse)
  k_fold_RMSE <- sqrt(k_fold_mse)
  k_fold_mRMSE <- mean(k_fold_RMSE)
  k_fold_mMAE <- mean(k_fold_mae)
  return(list(k_fold_mse,k_fold_mMSE,k_fold_RMSE,k_fold_mRMSE,k_fold_mae,k_fold_mMAE,simu_krig_data))
} 
#--------------  
n.times_tempLM.kfoldcv <-function(repeat.time,fold,basis,vari.dist){
  n.times_MSE.LM.kfoldcv <-c()
  n.times_RMSE.LM.kfoldcv <-c()  
  n.times_MAE.LM.kfoldcv <-c()    
  for (i in 1:repeat.time) {
    outLM <- k_folds_lm.ori.cv(fold = fold,bbaa = basis,vari.dist=vari.dist)
    n.times_MSE.LM.kfoldcv<- c(n.times_MSE.LM.kfoldcv ,outLM[[2]])
    n.times_RMSE.LM.kfoldcv<- c(n.times_RMSE.LM.kfoldcv ,outLM[[4]])
    n.times_MAE.LM.kfoldcv<- c(n.times_MAE.LM.kfoldcv ,outLM[[6]])
  }
  return(list(n.times_MSE.LM.kfoldcv, n.times_RMSE.LM.kfoldcv,n.times_MAE.LM.kfoldcv))
}

# reapet v times 
system.time({n.times_LM.kfoldcv <- n.times_tempLM.kfoldcv(repeat.time=10,fold=5,basis=3,vari.dist=600)})

par(mfrow = c(1, 3))
# k fold mse
plot(n.times_LM.kfoldcv[[1]], ylab="MSE",
     type = "b",col = "orange", pch = 20,main = "rep 10 times lm MSE")
# k fold RMSE
plot(n.times_LM.kfoldcv[[2]], ylab="RMSE",
     type = "b",col = "orange", pch = 20,main = "rep 10 times lm RMSE")
# k fold RMSE
plot(n.times_LM.kfoldcv[[3]], ylab="MAE",
     type = "b",col = "orange", pch = 20,main = "rep 10 times lm RMSE")
mean.lm.MSE=mean(n.times_LM.kfoldcv[[1]])
mean.lm.RMSE=mean(n.times_LM.kfoldcv[[2]])
mean.lm.MAE=mean(n.times_LM.kfoldcv[[3]])

#--------------------------------------------------------------------------------  
#---------------  k-folds cross-validation(GLM in mean not huber) (with original data)--------    
k_folds_glm.ori.cv <- function(fold,bbaa,vari.dist){
  k_fold_mse <-c()
  k_fold_mae <-c()
  #set.seed(6)
  cv  <- crossv_kfold(data = CA_case, k = fold ) # k-fold cross-Validation
  for (t in 1:length(cv$.id)) {
    #CA_case[,-4] is to elimate label
    case_train <- CA_case[cv$train[[t]]$idx,]
    case_test <- CA_case[cv$test[[t]]$idx,]
    nn <- length(case_train$long)
    ### number of basis = 5
    # use simulation location in mrts to generate basis
    #bbaa <- 10 # number of basis
    simubas <- as.matrix(mrts(knot = case_train[1:2], bbaa))
    
    #oc:outcome
    glmmodel<- glm(case_train$value ~ simubas-1)
    glmmodel[[1]]
    irls.simresidu <-glmmodel[[2]]
    
    #covariavce data = (observation data) - mean
    covadata<- case_train[,3] - simubas%*%glmmodel[[1]]
    
    
    # uvec : a vector with values used to define the variogram binning
    #variogram()
    #max.dist : length of variogram
    jj=max(abs(case_train[,1:2]))
    simuvario.b <- variog(coords = case_train[1:2], data = covadata,uvec=seq(0,jj, l=vari.dist), max.dist=6)
    #plot(simuvario.b, main="binned variogram") 
    
    simuini.vals <- expand.grid(seq(0.01,simuvario.b$var.mark,l=jj), seq(0.1,jj,l=jj)) # sigma^2 (partial sill) and phi (range parameter)
    simuols <- variofit(simuvario.b,fix.nug=T, wei="equal",cov.model="exp", ini=simuini.vals)
    #lines(simuols,lty=3,col=3,lwd=2)
    
    simuini.valss <- expand.grid(seq(0.01,simuols$cov.pars[1],l=jj), seq(0.1,jj,l=jj)) # sigma^2 (partial sill) and phi (range parameter)
    simuwls <- variofit(simuvario.b, ini=simuini.valss, fix.nug=F,wei="cressie",cov.model="exp")
    #lines(simuwls, lty=2,col=4,lwd=2)
    
    ###-------- use kriging model  (weight=wls) and see the model fit (MSE) -------#
    #simuweight = wei.select
    simuweight = simuwls
    # mean function residu
    # τ^2 + σ^2, corresponds to the variance of the observation process Y
    simresidu<-as.matrix(irls.simresidu)
    TPS_simu_beta_hat <- glmmodel[[1]]
    # simu_b_grid:grid point(simu_grid_poin) into TPS basis
    simu_b_grid <- as.matrix(mrts(knot = case_test[1:2], bbaa))
    simu_phi_s0 <- t(t(TPS_simu_beta_hat) %*% t(simu_b_grid))
    simu_c_s0 <- simuweight$cov.pars[1]*exp(-abs(as.data.frame.matrix(dist(case_test, case_train , method="euclidean")))/(simuweight$cov.pars[2]+1e-20))
    simu_zigma_theta <- simuweight$cov.pars[1]*as.matrix(exp(-abs(dist(case_train[-3], method = "euclidean", diag = F, upper = FALSE))/(simuweight$cov.pars[2]+1e-20))) #+1e-20: prevent simu_zigma_theta = 0 or Inf
    simu_y_hat =  simu_phi_s0 + as.matrix(simu_c_s0)  %*% ginv(simu_zigma_theta + (simuweight$nugget+simuweight$cov.pars[1])*diag(nn)) %*% simresidu
    
    # grid point & y_hat
    simu_krig_data <- cbind(case_test,simu_y_hat)
    
    # check MSE & resudual
    simu_MSE = MSE(case_test$value, simu_krig_data$simu_y_hat)
    k_fold_mse <- c(k_fold_mse ,simu_MSE)
    simu_MAE = mae(case_test$value, simu_krig_data$simu_y_hat)
    k_fold_mae <- c(k_fold_mae ,simu_MAE)  
  }
  k_fold_mMSE <- mean(k_fold_mse)
  k_fold_RMSE <- sqrt(k_fold_mse)
  k_fold_mRMSE <- mean(k_fold_RMSE)
  k_fold_mMAE <- mean(k_fold_mae)
  return(list(k_fold_mse,k_fold_mMSE,k_fold_RMSE,k_fold_mRMSE,k_fold_mae,k_fold_mMAE,simu_krig_data))
} 
#--------------  
n.times_tempGLM.kfoldcv <-function(repeat.time,fold,basis,vari.dist){
  n.times_MSE.GLM.kfoldcv <-c()
  n.times_RMSE.GLM.kfoldcv <-c()  
  n.times_MAE.GLM.kfoldcv <-c()    
  for (i in 1:repeat.time) {
    outGLM <- k_folds_glm.ori.cv(fold = fold,bbaa = basis,vari.dist=vari.dist)
    n.times_MSE.GLM.kfoldcv<- c(n.times_MSE.GLM.kfoldcv ,outGLM[[2]])
    n.times_RMSE.GLM.kfoldcv<- c(n.times_RMSE.GLM.kfoldcv ,outGLM[[4]])
    n.times_MAE.GLM.kfoldcv<- c(n.times_MAE.GLM.kfoldcv ,outGLM[[6]])
  }
  return(list(n.times_MSE.GLM.kfoldcv, n.times_RMSE.GLM.kfoldcv,n.times_MAE.GLM.kfoldcv))
}

# reapet v times 
system.time({n.times_GLM.kfoldcv <- n.times_tempGLM.kfoldcv(repeat.time=10,fold=5,basis=3,vari.dist=600)})


par(mfrow = c(1, 3))
# k fold mse
plot(n.times_GLM.kfoldcv[[1]], ylab="MSE",
     type = "b",col = "purple", pch = 20,main = "rep 10 times glm MSE")
# k fold RMSE
plot(n.times_GLM.kfoldcv[[2]], ylab="RMSE",
     type = "b",col = "purple", pch = 20,main = "rep 10 times glm MSE")
# k fold RMSE
plot(n.times_GLM.kfoldcv[[3]], ylab="MAE",
     type = "b",col = "purple", pch = 20,main = "rep 10 times glm MSE")
mean.glm.MSE=mean(n.times_GLM.kfoldcv[[1]])
mean.glm.RMSE=mean(n.times_GLM.kfoldcv[[2]])
mean.glm.MAE=mean(n.times_GLM.kfoldcv[[3]])

#-----------merge all plot ----------------

par(mfrow = c(1, 3))
# k fold mse
plot(n.times_WHM.kfoldcv[[1]], ylim = c(0.7,1),ylab="MSE",
     type = "b",col = "turquoise", pch = 20,main = "rep 10 times MSE")
lines(n.times_LM.kfoldcv[[1]], col = "orange", pch = 20,type = "b")
lines(n.times_GLM.kfoldcv[[1]], col = "purple", pch = 20,type = "b")
#lines(n.times_GAM.kfoldcv[[1]], col = "dodgerblue", pch = 20,type = "b")
# "legend" : icon in plot 
legend("bottomright",pch=c(20,20,20),legend=c("Huber","LM","GLM"),
       col=c("turquoise","orange","purple"),bty="n",box.lty = 0,cex=1)

# k fold RMSE
plot(n.times_WHM.kfoldcv[[2]], ylim = c(0.85,1),ylab="RMSE",
     type = "b",col = "turquoise", pch = 20,main = "rep 10 times RMSE")
lines(n.times_LM.kfoldcv[[2]], col = "orange", pch = 20,type = "b")
lines(n.times_GLM.kfoldcv[[2]], col = "purple", pch = 20,type = "b")
#lines(n.times_GAM.kfoldcv[[2]], col = "dodgerblue", pch = 20,type = "b")
legend("bottomright",pch=c(20,20,20),legend=c("Huber","LM","GLM"),
       col=c("turquoise","orange","purple"),bty="n",box.lty = 0,cex=1)

# k fold MAE
plot(n.times_WHM.kfoldcv[[3]], ylim = c(0.65,0.8),ylab="MAE",
     type = "b",col = "turquoise", pch = 20,main = "rep 10 timesMAE")
lines(n.times_LM.kfoldcv[[3]], col = "orange", pch = 20,type = "b")
lines(n.times_GLM.kfoldcv[[3]], col = "purple", pch = 20,type = "b")
#lines(n.times_GAM.kfoldcv[[3]], col = "dodgerblue", pch = 20,type = "b")
legend("bottomright",pch=c(20,20,20),legend=c("Huber","LM","GLM"),
       col=c("turquoise","orange","purple"),bty="n",box.lty = 0,cex=1)

