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

CA_case <- read.csv("D:\\研究所Meeting\\COVID19\\CA government data\\2022_4_18 CA case.csv")[-59,]
#Confirmed rate(單日確診率/per 100000 peoples) and normalization(or standardization)
CA_case$Confirmed_rate <- normalize(as.numeric(CA_case$reported_cases/CA_case$population)*100000, method = "range", range = c(0, 1))
#CA_case$Confirmed_rate <- scale(as.numeric(CA_case$reported_cases/CA_case$population)*100000)
#CA_case$Confirmed_rate <- as.numeric(CA_case$reported_cases/CA_case$population)*100000
# add label in data which used to general data with package "pdfos"
#CA_case$Class <-rep(as.character(c('positive', 'negative')), times=29)
# select CA_case : long、lat、Confirmed_rate & label
CA_case <- CA_case[,c(3,4,16)]
# normal test
#shapiro.test(CA_case$Confirmed_rate)
# add some perturbed in data
# CA_case_Confirmed_sd <- sd(CA_case$Confirmed_rate)
# n <- length(CA_case$Confirmed_rate)
# height_perturb <- rnorm(n, mean = 0, sd = CA_case_Confirmed_sd  * 0.1)
#-------------------------

#Leave one out cross-validation (LOOCV)

#--------------------------------mean function----------------------

#Standardization
stan_e<-function(ei){
  sigma<-median(abs(ei-median(ei)))/0.6745
  e_stan<-ei/sigma   
  return(e_stan)
}

new_residu<-function(X,Y,beta){
  tt<- X %*% beta
  new_res <- as.matrix(Y-tt)
  return(new_res)
}

# huber weight function(second differential)
hu_weifun <-function(residual.wei, k ){
  temp_2<-c() 
  # j=1
  # k=1.345
  for (j in 1:length(residual.wei)) {
    if(abs(residual.wei[j]) <= k){
      temp_2[j] <- 1
    } 
    else {
      temp_2[j] <- k/abs(residual.wei[j])
    }
  }
  return(temp_2)
}

# weighted least squares
WLS<-function(X,Y,W){
  w_b = solve(t(X) %*% W %*% X) %*% t(X) %*% W %*% Y
  return(w_b)
}

#iterative reweighted least squares
IRLS<-function(X,Y,max.iter,conv.eps){
  beta_init <- solve(t(X) %*% X) %*% t(X) %*% Y #initialize beta
  # lmmodel<- lm(formula= Y ~ X,data=data)[[1]]
  # names(lmmodel)<- NULL
  # beta_init<-lmmodel[-1]
  beta_prev <- beta_init               #beta_{t-1} (for comparisons)
  mean_fun<- X %*% beta_init
  e_i.prev<- as.matrix(Y-mean_fun)
  
  for(ii in 1:max.iter){
    standar_ei <- stan_e(new_residu(X,Y,beta_prev)) 
    W <-  diag(hu_weifun(residual.wei=standar_ei, k = 1.345))
    new_beta <- WLS(X=X,Y=Y,W=W)
    if( sum((beta_prev-new_beta)^2) <= conv.eps){
      break
    } else {
      e_i.prev <- standar_ei
      beta_prev <- new_beta
    }
    if(ii%%5 == 0) cat("iterative: ", ii, "\n")
  }
  return(list(new_beta,e_i.prev)) 
}

#------------------- k-flods cross-validation (with TPS general Standardize data) -----------------#
#[step1]use k-flods CV to divid data in train & test
#[step2]use spline function to generate new data in train data 
# "bbaa" is the number of basis
# "datanum" Decide how much data to generate
k_flods_gen.cv <- function(flod,bbaa,general.datanum){
  k_fold_mse <-c()
  k_fold_mae <-c()
  cv  <- crossv_kfold(data = CA_case, k = flod ) # k-fold cross-Validation  
  #set.seed(6)
 
  for (t in 1:length(cv$.id)) {
     # Spline interpolation for longitude and latitude
    case_train_temp <- CA_case[cv$train[[t]]$idx,]
    case_test <- CA_case[cv$test[[t]]$idx,]
    #set.seed(6)
    new_longitude <- with(case_train_temp,spline(case_train_temp$long, n = general.datanum))
    new_latitude <- with(case_train_temp, spline(case_train_temp$lat, n = general.datanum))
    # Create a new data frame containing the longitude and latitude of the new observation points
    case_train <- data.frame(
      long = new_longitude$y, # Longitude of the new observation point
      lat = new_latitude$y)   # Latitude of the new observation point
    # Build a TPS model
    tps_model <- Tps(case_train_temp[,1:2], case_train_temp[,3])
    # Predict the new number of "general.datanum" latitudes and longitudes and observations
    case_train$Confirmed_rate <-as.numeric(predict(tps_model, case_train[,1:2]))
    case_train <- rbind(case_train,case_train_temp)
    nn <- length(case_train$long)
    ### number of basis = 5
    # use simulation location in mrts to generate basis
    #bbaa <- 10 # number of basis
    simubas <- as.matrix(mrts(knot = case_train[1:2], bbaa))
    
    #oc:outcome
    oc<- IRLS(X = simubas,Y = case_train[,3],max.iter = 100,conv.eps = 1e-5)
    oc[[1]]
    irls.simresidu <-oc[[2]]
    
    #covariavce data = (observation data) - mean
    covadata<- case_train[,3] - simubas%*%oc[[1]]

    # uvec : a vector with values used to define the variogram binning
    #variogram()
    #max.dist : length of variogram
    jj=max(abs(case_train[,1:2]))
    simuvario.b <- variog(coords = case_train[1:2], data = covadata,uvec=seq(0,jj, l=300), max.dist=6)
    #plot(simuvario.b, main="binned variogram") 
    
    simuini.vals <- expand.grid(seq(0.1,simuvario.b$var.mark,l=jj), seq(0.05,jj,l=jj)) # sigma^2 (partial sill) and phi (range parameter)
    simuols <- variofit(simuvario.b,fix.nug=T, wei="equal",cov.model="exp", ini=simuini.vals)
    #lines(simuols,lty=3,col=3,lwd=2)
    
    simuini.valss <- expand.grid(seq(0.1,simuols$cov.pars[1],l=jj), seq(0.05,jj,l=jj)) # sigma^2 (partial sill) and phi (range parameter)
    simuwls <- variofit(simuvario.b, ini=simuini.valss, fix.nug=T,wei="npairs",cov.model="exp")
    #lines(simuwls, lty=2,col=4,lwd=2)
    
    ###-------- use kriging model  (weight=wls) and see the model fit (MSE) -------#
    #simuweight = wei.select
    simuweight = simuwls
    # mean function residu
    # τ^2 + σ^2, corresponds to the variance of the observation process Y
    simresidu<-as.matrix(irls.simresidu)
    TPS_simu_beta_hat <- oc[[1]]
    # simu_b_grid:grid point(simu_grid_poin) into TPS basis
    simu_b_grid <- as.matrix(mrts(knot = case_test[1:2], bbaa))
    simu_phi_s0 <- t(t(TPS_simu_beta_hat) %*% t(simu_b_grid))
    simu_c_s0 <- simuweight$cov.pars[1]*exp(-abs(as.data.frame.matrix(dist(case_test, case_train , method="euclidean")))/(simuweight$cov.pars[2]+1e-20))
    simu_zigma_theta <- simuweight$cov.pars[1]*as.matrix(exp(-abs(dist(case_train[-3], method = "euclidean", diag = F, upper = FALSE))/(simuweight$cov.pars[2]+1e-20))) #+1e-20: prevent simu_zigma_theta = 0 or Inf
    simu_y_hat =  simu_phi_s0 + as.matrix(simu_c_s0)  %*% ginv(simu_zigma_theta + (simuweight$nugget+simuweight$cov.pars[1])*diag(nn)) %*% simresidu
    
    # grid point & y_hat
    simu_krig_data <- cbind(case_test,simu_y_hat)
    
    # check MSE & resudual
    simu_MSE = MSE(case_test$Confirmed_rate, simu_krig_data$simu_y_hat)
    k_fold_mse <- c(k_fold_mse ,simu_MSE)
    simu_MAE = mae(case_test$Confirmed_rate, simu_krig_data$simu_y_hat)
    k_fold_mae <- c(k_fold_mae ,simu_MAE)   
  }
  k_fold_mMSE <- mean(k_fold_mse)
  k_fold_RMSE <- sqrt(k_fold_mse)
  k_fold_mRMSE <- mean(k_fold_RMSE)
  k_fold_mMAE <- mean(k_fold_mae)
  return(list(k_fold_mse,k_fold_mMSE,k_fold_RMSE,k_fold_mRMSE,k_fold_mae,k_fold_mMAE,simu_krig_data))
}
#----------------------------------general data outcome----------------- 
# "system.time({  })" Calculate running time 
#outcome_k_flods_cv <- k_flods_cv(flod = 10,bbaa = 10,wei.select = simuwls)
system.time({outcome_k_flods_gen.cv <- k_flods_gen.cv(flod = 5,bbaa = 10,general.datanum = 100)})



par(mfrow = c(1, 3))
# k fold mse
plot(outcome_k_flods_gen.cv[[1]], ylab="MSE",
     type = "b",col = "red", pch = 20,main = "5 fold MSE")
# k fold RMSE
plot(outcome_k_flods_gen.cv[[3]], ylab="RMSE",
     type = "b",col = "red", pch = 20,main = "5 fold RMSE")
# k fold MAE
plot(outcome_k_flods_gen.cv[[5]], ylab="MAE",
     type = "b",col = "red", pch = 20,main = "5 fold MAE")


#mean MSE
outcome_k_flods_gen.cv[[2]]
#mean RMSE
outcome_k_flods_gen.cv[[4]]
#mean MAE
outcome_k_flods_gen.cv[[6]]

#-------------------weight Huber m estimate package(WHM) (with original data)--------------
k_flods_WHM.cv <- function(flod,bbaa){
  k_fold_mse <-c()
  k_fold_mae <-c()
  #pdfos: Probability density function estimation based oversampling
  #set.seed(6)
  #new_data <- pdfos(data = CA_case, numInstances = datanum,classAttr = "Class")
  cv  <- crossv_kfold(data = CA_case[,1:3], k = flod ) # k-fold cross-Validation
  for (t in 1:length(cv$.id)) {
    #CA_case[,-4] is to elimate label
    case_train <- CA_case[cv$train[[t]]$idx,]
    case_test <- CA_case[cv$test[[t]]$idx,]
    nn <- length(case_train$long)
    ### number of basis = 5
    # use simulation location in mrts to generate basis
    #bbaa <- 10 # number of basis
    simubas <- as.matrix(mrts(knot = case_train[1:2], bbaa))
    
    CA_case_Design.matrix <- model.matrix(case_train$Confirmed_rate  ~ simubas)
    response.variable <- case_train$Confirmed_rate
    #oc:outcome
    WHM_p <- whm(response.variable, CA_case_Design.matrix[,-1], var.function = "power", tuning.para = 1.345)
    WHM_p[[1]]
    irls.simresidu <-WHM_p[[2]]
    
    #covariavce data = (observation data) - mean
    covadata<- case_train[,3] - simubas%*%WHM_p[[1]]

    # uvec : a vector with values used to define the variogram binning
    #variogram()
    #max.dist : length of variogram
    jj=max(abs(case_train[,1:2]))
    simuvario.b <- variog(coords = case_train[1:2], data = covadata,uvec=seq(0,jj, l=300), max.dist=6)
    #plot(simuvario.b, main="binned variogram") 
    
    simuini.vals <- expand.grid(seq(0.1,simuvario.b$var.mark,l=jj), seq(0.05,jj,l=jj)) # sigma^2 (partial sill) and phi (range parameter)
    simuols <- variofit(simuvario.b,fix.nug=T, wei="equal",cov.model="exp", ini=simuini.vals)
    #lines(simuols,lty=3,col=3,lwd=2)
    
    simuini.valss <- expand.grid(seq(0.1,simuols$cov.pars[1],l=jj), seq(0.05,jj,l=jj)) # sigma^2 (partial sill) and phi (range parameter)
    simuwls <- variofit(simuvario.b, ini=simuini.valss, fix.nug=T,wei="npairs",cov.model="exp")
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
    simu_c_s0 <- simuweight$cov.pars[1]*exp(-abs(as.data.frame.matrix(dist(case_test, case_train , method="euclidean")))/(simuweight$cov.pars[2]+1e-20))
    simu_zigma_theta <- simuweight$cov.pars[1]*as.matrix(exp(-abs(dist(case_train[-3], method = "euclidean", diag = F, upper = FALSE))/(simuweight$cov.pars[2]+1e-20))) #+1e-20: prevent simu_zigma_theta = 0 or Inf
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
    simu_MSE = MSE(case_test$Confirmed_rate, simu_krig_data$simu_y_hat)
    k_fold_mse <- c(k_fold_mse ,simu_MSE)
    simu_MAE = mae(case_test$Confirmed_rate, simu_krig_data$simu_y_hat)
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

system.time({outcome_k_flods_WHM.cv <- k_flods_WHM.cv(flod = 5,bbaa = 10)})

par(mfrow = c(1, 3))
# k fold mse
plot(outcome_k_flods_WHM.cv[[1]], ylab="MSE",
     type = "b",col = "turquoise", pch = 20,main = "5 fold MSE")
# k fold RMSE
plot(outcome_k_flods_WHM.cv[[3]], ylab="RMSE",
     type = "b",col = "turquoise", pch = 20,main = "5 fold RMSE")
# k fold MAE
plot(outcome_k_flods_WHM.cv[[3]], ylab="MAE",
     type = "b",col = "turquoise", pch = 20,main = "5 fold MAE")

#mean MSE
outcome_k_flods_WHM.cv[[2]]
#mean RMSE
outcome_k_flods_WHM.cv[[4]]
#mean MAE
outcome_k_flods_WHM.cv[[6]]

#------------------- k-flods cross-validation (with original data) -----------------#  
k_flods_ori.cv <- function(flod,bbaa){
  k_fold_mse <-c()
  k_fold_mae <-c()
  #set.seed(6)
  cv  <- crossv_kfold(data = CA_case[,1:3], k = flod ) # k-fold cross-Validation
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
    oc<- IRLS(X = simubas,Y = case_train[,3],max.iter = 100,conv.eps = 1e-5)
    oc[[1]]
    irls.simresidu <-oc[[2]]
    
    #covariavce data = (observation data) - mean
    covadata<- case_train[,3] - simubas%*%oc[[1]]

    
    # uvec : a vector with values used to define the variogram binning
    #variogram()
    #max.dist : length of variogram
    jj=max(abs(case_train[,1:2]))
    simuvario.b <- variog(coords = case_train[1:2], data = covadata,uvec=seq(0,jj, l=300), max.dist=6)
    #plot(simuvario.b, main="binned variogram") 
    
    simuini.vals <- expand.grid(seq(0.1,simuvario.b$var.mark,l=jj), seq(0.05,jj,l=jj)) # sigma^2 (partial sill) and phi (range parameter)
    simuols <- variofit(simuvario.b,fix.nug=T, wei="equal",cov.model="exp", ini=simuini.vals)
    #lines(simuols,lty=3,col=3,lwd=2)
    
    simuini.valss <- expand.grid(seq(0.1,simuols$cov.pars[1],l=jj), seq(0.05,jj,l=jj)) # sigma^2 (partial sill) and phi (range parameter)
    simuwls <- variofit(simuvario.b, ini=simuini.valss, fix.nug=T,wei="npairs",cov.model="exp")
    #lines(simuwls, lty=2,col=4,lwd=2)
    
    ###-------- use kriging model  (weight=wls) and see the model fit (MSE) -------#
    #simuweight = wei.select
    simuweight = simuwls
    # mean function residu
    # τ^2 + σ^2, corresponds to the variance of the observation process Y
    simresidu<-as.matrix(irls.simresidu)
    TPS_simu_beta_hat <- oc[[1]]
    # simu_b_grid:grid point(simu_grid_poin) into TPS basis
    simu_b_grid <- as.matrix(mrts(knot = case_test[1:2], bbaa))
    simu_phi_s0 <- t(t(TPS_simu_beta_hat) %*% t(simu_b_grid))
    simu_c_s0 <- simuweight$cov.pars[1]*exp(-abs(as.data.frame.matrix(dist(case_test, case_train , method="euclidean")))/(simuweight$cov.pars[2]+1e-20))
    simu_zigma_theta <- simuweight$cov.pars[1]*as.matrix(exp(-abs(dist(case_train[-3], method = "euclidean", diag = F, upper = FALSE))/(simuweight$cov.pars[2]+1e-20))) #+1e-20: prevent simu_zigma_theta = 0 or Inf
    simu_y_hat =  simu_phi_s0 + as.matrix(simu_c_s0)  %*% ginv(simu_zigma_theta + (simuweight$nugget+simuweight$cov.pars[1])*diag(nn)) %*% simresidu
    
    # grid point & y_hat
    simu_krig_data <- cbind(case_test,simu_y_hat)
    
    # check MSE & resudual
    simu_MSE = MSE(case_test$Confirmed_rate, simu_krig_data$simu_y_hat)
    k_fold_mse <- c(k_fold_mse ,simu_MSE)
    simu_MAE = mae(case_test$Confirmed_rate, simu_krig_data$simu_y_hat)
    k_fold_mae <- c(k_fold_mae ,simu_MAE) 
  }
  k_fold_mMSE <- mean(k_fold_mse)
  k_fold_RMSE <- sqrt(k_fold_mse)
  k_fold_mRMSE <- mean(k_fold_RMSE)
  k_fold_mMAE <- mean(k_fold_mae)
  return(list(k_fold_mse,k_fold_mMSE,k_fold_RMSE,k_fold_mRMSE,k_fold_mae,k_fold_mMAE,simu_krig_data))
} 
#-----------------------------------original data outcome--------------------
system.time({outcome_k_flods_ori.cv <- k_flods_ori.cv(flod = 5,bbaa = 10)})

ori.realdot<-outcome_k_flods_ori.cv[[7]]%>% as.data.frame %>% 
  ggplot(aes(long, lat)) + geom_point(aes(colour = Confirmed_rate)) + 
  labs(title="real dot") + scale_color_viridis(option = "H")
ori.predot<-outcome_k_flods_ori.cv[[7]]%>% as.data.frame %>% 
  ggplot(aes(long, lat)) + geom_point(aes(colour = simu_y_hat)) + 
  labs(title="ori pre dot") + scale_color_viridis(option = "H")

ori.dotvs <-subplot(ori.predot,ori.realdot)%>% 
  layout(title = "ori predict   &    real")
ori.dotvs


par(mfrow = c(1, 3))
# k fold mse
plot(outcome_k_flods_ori.cv[[1]], ylab="MSE",
     type = "b",col = "magenta", pch = 20,main = "5 fold ori MSE")
# k fold RMSE
plot(outcome_k_flods_ori.cv[[3]], ylab="RMSE",
     type = "b",col = "magenta", pch = 20,main = "5 fold ori RMSE")
# k fold MAE
plot(outcome_k_flods_ori.cv[[5]], ylab="MAE",
     type = "b",col = "magenta", pch = 20,main = "5 fold ori RMSE")

#mean MSE
outcome_k_flods_ori.cv[[2]]
#mean RMSE
outcome_k_flods_ori.cv[[4]]
#mean MAE
outcome_k_flods_ori.cv[[6]]


#---------------  k-flods cross-validation(LM in mean not huber) (with original data)--------
k_flods_lm.ori.cv <- function(flod,bbaa){
  k_fold_mse <-c()
  k_fold_mae <-c()
  #set.seed(6)
  cv  <- crossv_kfold(data = CA_case[,1:3], k = flod ) # k-fold cross-Validation
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
    lmmodel<- lm(case_train$Confirmed_rate ~ simubas-1)
    lmmodel[[1]]
    irls.simresidu <-lmmodel[[2]]
    
    #covariavce data = (observation data) - mean
    covadata<- case_train[,3] - simubas%*%lmmodel[[1]]

    # uvec : a vector with values used to define the variogram binning
    #variogram()
    #max.dist : length of variogram
    jj=max(abs(case_train[,1:2]))
    simuvario.b <- variog(coords = case_train[1:2], data = covadata,uvec=seq(0,jj, l=300), max.dist=6)
    #plot(simuvario.b, main="binned variogram") 
    
    simuini.vals <- expand.grid(seq(0.1,simuvario.b$var.mark,l=jj), seq(0.05,jj,l=jj)) # sigma^2 (partial sill) and phi (range parameter)
    simuols <- variofit(simuvario.b,fix.nug=T, wei="equal",cov.model="exp", ini=simuini.vals)
    #lines(simuols,lty=3,col=3,lwd=2)
    
    simuini.valss <- expand.grid(seq(0.1,simuols$cov.pars[1],l=jj), seq(0.05,jj,l=jj)) # sigma^2 (partial sill) and phi (range parameter)
    simuwls <- variofit(simuvario.b, ini=simuini.valss, fix.nug=T,wei="npairs",cov.model="exp")
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
    simu_MSE = MSE(case_test$Confirmed_rate, simu_krig_data$simu_y_hat)
    k_fold_mse <- c(k_fold_mse ,simu_MSE)
    simu_MAE = mae(case_test$Confirmed_rate, simu_krig_data$simu_y_hat)
    k_fold_mae <- c(k_fold_mae ,simu_MAE)  
  }
  k_fold_mMSE <- mean(k_fold_mse)
  k_fold_RMSE <- sqrt(k_fold_mse)
  k_fold_mRMSE <- mean(k_fold_RMSE)
  k_fold_mMAE <- mean(k_fold_mae)
  return(list(k_fold_mse,k_fold_mMSE,k_fold_RMSE,k_fold_mRMSE,k_fold_mae,k_fold_mMAE,simu_krig_data))
} 
#--------------  
system.time({outcome_k_flods_lm.ori.cv <- k_flods_lm.ori.cv(flod = 5,bbaa = 10)})

lm.ori.realdot<-outcome_k_flods_lm.ori.cv[[5]]%>% as.data.frame %>% 
  ggplot(aes(long, lat)) + geom_point(aes(colour = Confirmed_rate)) + 
  labs(title="real dot") + scale_color_viridis(option = "H")
lm.ori.predot<-outcome_k_flods_lm.ori.cv[[5]]%>% as.data.frame %>% 
  ggplot(aes(long, lat)) + geom_point(aes(colour = simu_y_hat)) + 
  labs(title="lm pre dot") + scale_color_viridis(option = "H")

lm.ori.dotvs <-subplot(lm.ori.predot,lm.ori.realdot)%>% 
  layout(title = "lm predict   &    real")
lm.ori.dotvs


par(mfrow = c(1, 3))
# k fold mse
plot(outcome_k_flods_lm.ori.cv[[1]], ylab="MSE",
     type = "b",col = "orange", pch = 20,main = "5 fold lm MSE")
# k fold RMSE
plot(outcome_k_flods_lm.ori.cv[[3]], ylab="RMSE",
     type = "b",col = "orange", pch = 20,main = "5 fold lm RMSE")
# k fold RMSE
plot(outcome_k_flods_lm.ori.cv[[5]], ylab="MAE",
     type = "b",col = "orange", pch = 20,main = "5 fold lm RMSE")

#mean MSE
outcome_k_flods_lm.ori.cv[[2]]
#mean RMSE
outcome_k_flods_lm.ori.cv[[4]]
#mean MAE
outcome_k_flods_lm.ori.cv[[6]]
#--------------------------------------------------------------------------------  

#---------------  k-flods cross-validation(GLM in mean not huber) (with original data)--------    
k_flods_glm.ori.cv <- function(flod,bbaa){
  k_fold_mse <-c()
  k_fold_mae <-c()
  #set.seed(6)
  cv  <- crossv_kfold(data = CA_case[,1:3], k = flod ) # k-fold cross-Validation
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
    glmmodel<- glm(case_train$Confirmed_rate ~ simubas-1)
    glmmodel[[1]]
    irls.simresidu <-glmmodel[[2]]
    
    #covariavce data = (observation data) - mean
    covadata<- case_train[,3] - simubas%*%glmmodel[[1]]

    
    # uvec : a vector with values used to define the variogram binning
    #variogram()
    #max.dist : length of variogram
    jj=max(abs(case_train[,1:2]))
    simuvario.b <- variog(coords = case_train[1:2], data = covadata,uvec=seq(0,jj, l=300), max.dist=6)
    #plot(simuvario.b, main="binned variogram") 
    
    simuini.vals <- expand.grid(seq(0.1,simuvario.b$var.mark,l=jj), seq(0.05,jj,l=jj)) # sigma^2 (partial sill) and phi (range parameter)
    simuols <- variofit(simuvario.b,fix.nug=T, wei="equal",cov.model="exp", ini=simuini.vals)
    #lines(simuols,lty=3,col=3,lwd=2)
    
    simuini.valss <- expand.grid(seq(0.1,simuols$cov.pars[1],l=jj), seq(0.05,jj,l=jj)) # sigma^2 (partial sill) and phi (range parameter)
    simuwls <- variofit(simuvario.b, ini=simuini.valss, fix.nug=T,wei="npairs",cov.model="exp")
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
    simu_MSE = MSE(case_test$Confirmed_rate, simu_krig_data$simu_y_hat)
    k_fold_mse <- c(k_fold_mse ,simu_MSE)
    simu_MAE = mae(case_test$Confirmed_rate, simu_krig_data$simu_y_hat)
    k_fold_mae <- c(k_fold_mae ,simu_MAE)  
  }
  k_fold_mMSE <- mean(k_fold_mse)
  k_fold_RMSE <- sqrt(k_fold_mse)
  k_fold_mRMSE <- mean(k_fold_RMSE)
  k_fold_mMAE <- mean(k_fold_mae)
  return(list(k_fold_mse,k_fold_mMSE,k_fold_RMSE,k_fold_mRMSE,k_fold_mae,k_fold_mMAE,simu_krig_data))
} 
#--------------  
system.time({outcome_k_flods_glm.ori.cv <- k_flods_glm.ori.cv(flod = 5,bbaa = 10)})

glm.ori.realdot<-outcome_k_flods_glm.ori.cv[[7]]%>% as.data.frame %>% 
  ggplot(aes(long, lat)) + geom_point(aes(colour = Confirmed_rate)) + 
  labs(title="real dot") + scale_color_viridis(option = "H")
glm.ori.predot<-outcome_k_flods_glm.ori.cv[[7]]%>% as.data.frame %>% 
  ggplot(aes(long, lat)) + geom_point(aes(colour = simu_y_hat)) + 
  labs(title="glm pre dot") + scale_color_viridis(option = "H")

glm.ori.dotvs <-subplot(glm.ori.predot,glm.ori.realdot)%>% 
  layout(title = "glm predict   &    real")
glm.ori.dotvs


par(mfrow = c(1, 3))
# k fold mse
plot(outcome_k_flods_glm.ori.cv[[1]], ylab="MSE",
     type = "b",col = "purple", pch = 20,main = "5 fold glm MSE")
# k fold RMSE
plot(outcome_k_flods_glm.ori.cv[[3]], ylab="RMSE",
     type = "b",col = "purple", pch = 20,main = "5 fold glm MSE")
# k fold RMSE
plot(outcome_k_flods_glm.ori.cv[[5]], ylab="MAE",
     type = "b",col = "purple", pch = 20,main = "5 fold glm MSE")


#mean MSE
outcome_k_flods_glm.ori.cv[[2]]
#mean RMSE
outcome_k_flods_glm.ori.cv[[4]]
#mean MAE
outcome_k_flods_glm.ori.cv[[6]]

#------------------------------
#---------------  k-flods cross-validation(GAM in mean not huber) (with original data)--------    
k_flods_gam.ori.cv <- function(flod,bbaa){
  k_fold_mse <-c()
  k_fold_mae <-c()
  #set.seed(7)
  cv  <- crossv_kfold(data = CA_case[,1:3], k = flod ) # k-fold cross-Validation
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
    gammodel<- gam(case_train$Confirmed_rate ~ simubas-1)
    gammodel[[1]]
    irls.simresidu <-gammodel[[2]]
    
    #covariavce data = (observation data) - mean
    covadata<- case_train[,3] - simubas%*%gammodel[[1]]

    # uvec : a vector with values used to define the variogram binning
    #variogram()
    #max.dist : length of variogram
    jj=max(abs(case_train[,1:2]))
    simuvario.b <- variog(coords = case_train[1:2], data = covadata,uvec=seq(0,jj, l=300), max.dist=6)
    #plot(simuvario.b, main="binned variogram") 
    
    simuini.vals <- expand.grid(seq(0.1,simuvario.b$var.mark,l=jj), seq(0.05,jj,l=jj)) # sigma^2 (partial sill) and phi (range parameter)
    simuols <- variofit(simuvario.b,fix.nug=T, wei="equal",cov.model="exp", ini=simuini.vals)
    #lines(simuols,lty=3,col=3,lwd=2)
    
    simuini.valss <- expand.grid(seq(0.1,simuols$cov.pars[1],l=jj), seq(0.05,jj,l=jj)) # sigma^2 (partial sill) and phi (range parameter)
    simuwls <- variofit(simuvario.b, ini=simuini.valss, fix.nug=T,wei="npairs",cov.model="exp")
    #lines(simuwls, lty=2,col=4,lwd=2)
    
    ###-------- use kriging model  (weight=wls) and see the model fit (MSE) -------#
    #simuweight = wei.select
    simuweight = simuwls
    # mean function residu
    # τ^2 + σ^2, corresponds to the variance of the observation process Y
    simresidu<-as.matrix(irls.simresidu)
    TPS_simu_beta_hat <- gammodel[[1]]
    # simu_b_grid:grid point(simu_grid_poin) into TPS basis
    simu_b_grid <- as.matrix(mrts(knot = case_test[1:2], bbaa))
    simu_phi_s0 <- t(t(TPS_simu_beta_hat) %*% t(simu_b_grid))
    simu_c_s0 <- simuweight$cov.pars[1]*exp(-abs(as.data.frame.matrix(dist(case_test, case_train , method="euclidean")))/(simuweight$cov.pars[2]+1e-20))
    simu_zigma_theta <- simuweight$cov.pars[1]*as.matrix(exp(-abs(dist(case_train[-3], method = "euclidean", diag = F, upper = FALSE))/(simuweight$cov.pars[2]+1e-20))) #+1e-20: prevent simu_zigma_theta = 0 or Inf
    simu_y_hat =  simu_phi_s0 + as.matrix(simu_c_s0)  %*% ginv(simu_zigma_theta + (simuweight$nugget+simuweight$cov.pars[1])*diag(nn)) %*% simresidu
    
    # grid point & y_hat
    simu_krig_data <- cbind(case_test,simu_y_hat)
    
    # check MSE & resudual
    simu_MSE = MSE(case_test$Confirmed_rate, simu_krig_data$simu_y_hat)
    k_fold_mse <- c(k_fold_mse ,simu_MSE)
    simu_MAE = mae(case_test$Confirmed_rate, simu_krig_data$simu_y_hat)
    k_fold_mae <- c(k_fold_mae ,simu_MAE) 
  }
  k_fold_mMSE <- mean(k_fold_mse)
  k_fold_RMSE <- sqrt(k_fold_mse)
  k_fold_mRMSE <- mean(k_fold_RMSE)
  k_fold_mMAE <- mean(k_fold_mae)
  return(list(k_fold_mse,k_fold_mMSE,k_fold_RMSE,k_fold_mRMSE,k_fold_mae,k_fold_mMAE,simu_krig_data))
} 
#--------------   
system.time({outcome_k_flods_gam.ori.cv <- k_flods_gam.ori.cv(flod = 5,bbaa = 10)})

gam.ori.realdot<-outcome_k_flods_gam.ori.cv[[7]]%>% as.data.frame %>% 
  ggplot(aes(long, lat)) + geom_point(aes(colour = Confirmed_rate)) + 
  labs(title="real dot") + scale_color_viridis(option = "H")
gam.ori.predot<-outcome_k_flods_gam.ori.cv[[7]]%>% as.data.frame %>% 
  ggplot(aes(long, lat)) + geom_point(aes(colour = simu_y_hat)) + 
  labs(title="gam pre dot") + scale_color_viridis(option = "H")

gam.ori.dotvs <-subplot(gam.ori.predot,gam.ori.realdot)%>% 
  layout(title = "gam predict   &    real")
gam.ori.dotvs

# par(mfrow = c(1, 2))
# # 3D plot in estimate value---1
# scatter3D(outcome_k_flods_ori.cv[[5]]$long, outcome_k_flods_ori.cv[[5]]$lat, 
#           outcome_k_flods_ori.cv[[5]]$Confirmed_rate,
#           zlab="value",phi = 10, theta = 60,main="real value Plot")
# scatter3D(outcome_k_flods_ori.cv[[5]]$long, outcome_k_flods_ori.cv[[5]]$lat, 
#           outcome_k_flods_ori.cv[[5]]$simu_y_hat,
#           zlab="value",phi = 10, theta = 60,main="estimate value Plot")

par(mfrow = c(1, 3))
# k fold mse
plot(outcome_k_flods_gam.ori.cv[[1]], ylab="MSE",
     type = "b",col = "dodgerblue", pch = 20,main = "5 fold gam MSE")
# k fold RMSE
plot(outcome_k_flods_gam.ori.cv[[3]], ylab="RMSE",
     type = "b",col = "dodgerblue", pch = 20,main = "5 fold gam MSE")
# k fold RMSE
plot(outcome_k_flods_gam.ori.cv[[5]], ylab="MAE",
     type = "b",col = "dodgerblue", pch = 20,main = "5 fold gam MSE")


#mean MSE
outcome_k_flods_gam.ori.cv[[2]]
#mean RMSE
outcome_k_flods_gam.ori.cv[[4]]
#mean RMSE
outcome_k_flods_gam.ori.cv[[6]]



#-----------merge all plot ----------------
par(mfrow = c(1, 3))
# k fold mse
plot(outcome_k_flods_WHM.cv[[1]], ylim = c(0,0.4),ylab="MSE",
     type = "b",col = "turquoise", pch = 20,main = "5 fold MSE")
lines(outcome_k_flods_lm.ori.cv[[1]], col = "orange", pch = 20,type = "b")
lines(outcome_k_flods_glm.ori.cv[[1]], col = "purple", pch = 20,type = "b")
lines(outcome_k_flods_gam.ori.cv[[1]], col = "dodgerblue", pch = 20,type = "b")
# "legend" : icon in plot 
legend("bottomright",pch=c(20,20,20),legend=c("Huber","LM","GLM","GAM"),
       col=c("turquoise","orange","purple","dodgerblue"),bty="n",box.lty = 0,cex=1)
#legend("topleft",pch=c(20,20,20),legend=c("my huber","LM","GLM","GAM","huber in package"),
#       col=c("magenta","gray","purple","dodgerblue","turquoise"),bty="n",box.lty = 0,cex=1)


# k fold RMSE
plot(outcome_k_flods_WHM.cv[[3]], ylim = c(0,0.4),ylab="RMSE",
     type = "b",col = "turquoise", pch = 20,main = "5 fold RMSE")
lines(outcome_k_flods_lm.ori.cv[[3]], col = "orange", pch = 20,type = "b")
lines(outcome_k_flods_glm.ori.cv[[3]], col = "purple", pch = 20,type = "b")
lines(outcome_k_flods_gam.ori.cv[[3]], col = "dodgerblue", pch = 20,type = "b")
#lines(outcome_k_flods_gen.cv[[3]], col = "red", pch = 20,type = "b")
#lines(outcome_k_flods_WHM.cv[[3]], col = "turquoise", pch = 20,type = "b")
#lines(outcome_k_flods_ori.cv[[3]], col = "magenta", pch = 20,type = "b")
legend("bottomright",pch=c(20,20,20),legend=c("Huber","LM","GLM","GAM"),
       col=c("turquoise","orange","purple","dodgerblue"),bty="n",box.lty = 0,cex=1)
#legend("topleft",pch=c(20,20,20),legend=c("my huber","LM","GLM","GAM","huber in package"),
#       col=c("magenta","gray","purple","dodgerblue","turquoise"),bty="n",box.lty = 0,cex=0.55)

# k fold MAE
plot(outcome_k_flods_WHM.cv[[5]], ylim = c(0,0.4),ylab="MAE",
     type = "b",col = "turquoise", pch = 20,main = "5 fold MAE")
lines(outcome_k_flods_lm.ori.cv[[5]], col = "orange", pch = 20,type = "b")
lines(outcome_k_flods_glm.ori.cv[[5]], col = "purple", pch = 20,type = "b")
lines(outcome_k_flods_gam.ori.cv[[5]], col = "dodgerblue", pch = 20,type = "b")
legend("bottomright",pch=c(20,20,20),legend=c("Huber","LM","GLM","GAM"),
       col=c("turquoise","orange","purple","dodgerblue"),bty="n",box.lty = 0,cex=1)
