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

#Confirmed rate(單日確診率/per 100000 peoples) and normalization(or standardization)
CA_case <- read.csv("D:\\研究所Meeting\\COVID19\\CA government data\\2022_4 time series CA case.csv")
CA_pm2.5_total_out <- read.csv("D:\\研究所Meeting\\COVID19\\PM2.5 in CA\\CA_pm2.5_total_out.csv")
CA_pm2.5_total_out$several_days <- scale(as.numeric(apply(CA_pm2.5_total_out[,4:19], 1, mean)))
CA_pm2.5_total_KNNout <- CA_pm2.5_total_out[,c(1:3,20)]
# "apply(df, 1, mean)" 1 means Calculate the average of each row
#CA_case$Confirmed_rate <- normalize(as.numeric(apply(CA_case[,23:39], 1, mean)),method = "range", range = c(0, 1))      
CA_case$Confirmed_rate <- scale(as.numeric(apply(CA_case[,23:39], 1, mean)))      
CA_case <- cbind(CA_case[,c(4,5,40)],CA_pm2.5_total_KNNout[,4])

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

#---------------------------------[star here]---------------------------
#--------------------------------------------
#-------------------weight Huber m estimate package(WHM) --------------
k_folds_WHM.cv <- function(fold,bbaa){
  k_fold_mse <-c()
  k_fold_mae <-c()
  #pdfos: Probability density function estimation based oversampling
  #set.seed(6)
  #new_data <- pdfos(data = CA_case, numInstances = datanum,classAttr = "Class")
  cv  <- crossv_kfold(data = CA_case[,1:3], k = fold ) # k-fold cross-Validation
  for (t in 1:length(cv$.id)) {
    #CA_case[,-4] is to elimate label
    case_train <- CA_case[cv$train[[t]]$idx,]
    case_test <- CA_case[cv$test[[t]]$idx,]
    nn <- length(case_train$long)
    ### number of basis = 5
    # use simulation location in mrts to generate basis
    #bbaa <- 10 # number of basis
    simubas <- as.matrix(cbind(mrts(knot = case_train[1:2], bbaa),case_train$`CA_pm2.5_total_KNNout[, 4]`))
    
    CA_case_Design.matrix <- model.matrix(case_train$Confirmed_rate  ~ simubas)
    response.variable <- case_train$Confirmed_rate
    #oc:outcome
    #ite : Number of iterations for the estimation procedure.
    #var.function = "power" or "exponential"
    WHM_p <- whm(response.variable, CA_case_Design.matrix[,-1], var.function = "power", tuning.para = 1.345, ite = 200)
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
    simu_b_grid <- as.matrix(cbind(mrts(knot = case_test[1:2], bbaa),case_test$`CA_pm2.5_total_KNNout[, 4]`))
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

n.times_tempWHM.kfoldcv <-function(repeat.time,fold,basis){
  n.times_MSE.WHM.kfoldcv <-c()
  n.times_RMSE.WHM.kfoldcv <-c()  
  n.times_MAE.WHM.kfoldcv <-c()    
  for (i in 1:repeat.time) {
    outWHM <- k_folds_WHM.cv(fold = fold,bbaa = basis)
    n.times_MSE.WHM.kfoldcv<- c(n.times_MSE.WHM.kfoldcv ,outWHM[[2]])
    n.times_RMSE.WHM.kfoldcv<- c(n.times_RMSE.WHM.kfoldcv ,outWHM[[4]])
    n.times_MAE.WHM.kfoldcv<- c(n.times_MAE.WHM.kfoldcv ,outWHM[[6]])
  }
  return(list(n.times_MSE.WHM.kfoldcv, n.times_RMSE.WHM.kfoldcv,n.times_MAE.WHM.kfoldcv))
}

# reapet v times 
system.time({n.times_WHM.kfoldcv <- n.times_tempWHM.kfoldcv(repeat.time=10,fold=5,basis=10)})

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
 

#---------------  k-flods cross-validation(LM in mean not huber) (with original data)--------
k_folds_lm.ori.cv <- function(repeat.time,fold,bbaa){
  k_fold_mse <-c()
  k_fold_mae <-c()
  #set.seed(6)
  cv  <- crossv_kfold(data = CA_case[,1:3], k = fold ) # k-fold cross-Validation
  for (t in 1:length(cv$.id)) {
    #CA_case[,-4] is to elimate label
    case_train <- CA_case[cv$train[[t]]$idx,]
    case_test <- CA_case[cv$test[[t]]$idx,]
    nn <- length(case_train$long)
    ### number of basis = 5
    # use simulation location in mrts to generate basis
    #bbaa <- 10 # number of basis
    simubas <- as.matrix(cbind(mrts(knot = case_train[1:2], bbaa),case_train$`CA_pm2.5_total_KNNout[, 4]`))
    
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
    simu_b_grid <- as.matrix(cbind(mrts(knot = case_test[1:2], bbaa),case_test$`CA_pm2.5_total_KNNout[, 4]`))
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
n.times_tempLM.kfoldcv <-function(repeat.time,fold,basis){
  n.times_MSE.LM.kfoldcv <-c()
  n.times_RMSE.LM.kfoldcv <-c()  
  n.times_MAE.LM.kfoldcv <-c()    
  for (i in 1:repeat.time) {
    outLM <- k_folds_lm.ori.cv(fold = fold,bbaa = basis)
    n.times_MSE.LM.kfoldcv<- c(n.times_MSE.LM.kfoldcv ,outLM[[2]])
    n.times_RMSE.LM.kfoldcv<- c(n.times_RMSE.LM.kfoldcv ,outLM[[4]])
    n.times_MAE.LM.kfoldcv<- c(n.times_MAE.LM.kfoldcv ,outLM[[6]])
  }
  return(list(n.times_MSE.LM.kfoldcv, n.times_RMSE.LM.kfoldcv,n.times_MAE.LM.kfoldcv))
}

# reapet v times 
system.time({n.times_LM.kfoldcv <- n.times_tempLM.kfoldcv(repeat.time=10,fold=5,basis=10)})

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


#--------------------------------------------------------------------------------  

#---------------  k-folds cross-validation(GLM in mean not huber) (with original data)--------    
k_folds_glm.ori.cv <- function(fold,bbaa){
  k_fold_mse <-c()
  k_fold_mae <-c()
  #set.seed(6)
  cv  <- crossv_kfold(data = CA_case[,1:3], k = fold ) # k-fold cross-Validation
  for (t in 1:length(cv$.id)) {
    #CA_case[,-4] is to elimate label
    case_train <- CA_case[cv$train[[t]]$idx,]
    case_test <- CA_case[cv$test[[t]]$idx,]
    nn <- length(case_train$long)
    ### number of basis = 5
    # use simulation location in mrts to generate basis
    #bbaa <- 10 # number of basis
    simubas <- as.matrix(cbind(mrts(knot = case_train[1:2], bbaa),case_train$`CA_pm2.5_total_KNNout[, 4]`))
    
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
    simu_b_grid <- as.matrix(cbind(mrts(knot = case_test[1:2], bbaa),case_test$`CA_pm2.5_total_KNNout[, 4]`))
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
n.times_tempGLM.kfoldcv <-function(repeat.time,fold,basis){
  n.times_MSE.GLM.kfoldcv <-c()
  n.times_RMSE.GLM.kfoldcv <-c()  
  n.times_MAE.GLM.kfoldcv <-c()    
  for (i in 1:repeat.time) {
    outGLM <- k_folds_glm.ori.cv(fold = fold,bbaa = basis)
    n.times_MSE.GLM.kfoldcv<- c(n.times_MSE.GLM.kfoldcv ,outGLM[[2]])
    n.times_RMSE.GLM.kfoldcv<- c(n.times_RMSE.GLM.kfoldcv ,outGLM[[4]])
    n.times_MAE.GLM.kfoldcv<- c(n.times_MAE.GLM.kfoldcv ,outGLM[[6]])
  }
  return(list(n.times_MSE.GLM.kfoldcv, n.times_RMSE.GLM.kfoldcv,n.times_MAE.GLM.kfoldcv))
}

# reapet v times 
system.time({n.times_GLM.kfoldcv <- n.times_tempGLM.kfoldcv(repeat.time=10,fold=5,basis=10)})


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



#------------------------------
#---------------  k-folds cross-validation(GAM in mean not huber) (with original data)--------    
k_folds_gam.ori.cv <- function(fold,bbaa){
  k_fold_mse <-c()
  k_fold_mae <-c()
  #set.seed(6)
  cv  <- crossv_kfold(data = CA_case[,1:3], k = fold ) # k-fold cross-Validation
  for (t in 1:length(cv$.id)) {
    #CA_case[,-4] is to elimate label
    case_train <- CA_case[cv$train[[t]]$idx,]
    case_test <- CA_case[cv$test[[t]]$idx,]
    nn <- length(case_train$long)
    ### number of basis = 5
    # use simulation location in mrts to generate basis
    #bbaa <- 10 # number of basis
    simubas <- as.matrix(cbind(mrts(knot = case_train[1:2], bbaa),case_train$`CA_pm2.5_total_KNNout[, 4]`))
    
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
    simu_b_grid <- as.matrix(cbind(mrts(knot = case_test[1:2], bbaa),case_test$`CA_pm2.5_total_KNNout[, 4]`))
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
n.times_tempGAM.kfoldcv <-function(repeat.time,fold,basis){
  n.times_MSE.GAM.kfoldcv <-c()
  n.times_RMSE.GAM.kfoldcv <-c()  
  n.times_MAE.GAM.kfoldcv <-c()    
  for (i in 1:repeat.time) {
    outGAM <- k_folds_gam.ori.cv(fold = fold,bbaa = basis)
    n.times_MSE.GAM.kfoldcv<- c(n.times_MSE.GAM.kfoldcv ,outGAM[[2]])
    n.times_RMSE.GAM.kfoldcv<- c(n.times_RMSE.GAM.kfoldcv ,outGAM[[4]])
    n.times_MAE.GAM.kfoldcv<- c(n.times_MAE.GAM.kfoldcv ,outGAM[[6]])
  }
  return(list(n.times_MSE.GAM.kfoldcv, n.times_RMSE.GAM.kfoldcv,n.times_MAE.GAM.kfoldcv))
}

# reapet v times 
system.time({n.times_GAM.kfoldcv <- n.times_tempGAM.kfoldcv(repeat.time=10,fold=5,basis=10)})

par(mfrow = c(1, 3))
# k fold mse
plot(n.times_GAM.kfoldcv[[1]], ylab="MSE",
     type = "b",col = "dodgerblue", pch = 20,main = "rep 10 times gam MSE")
# k fold RMSE
plot(n.times_GAM.kfoldcv[[2]], ylab="RMSE",
     type = "b",col = "dodgerblue", pch = 20,main = "rep 10 times gam MSE")
# k fold RMSE
plot(n.times_GAM.kfoldcv[[3]], ylab="MAE",
     type = "b",col = "dodgerblue", pch = 20,main = "rep 10 times gam MSE")


#-----------merge all plot ----------------
par(mfrow = c(1, 3))
# k fold mse
plot(n.times_WHM.kfoldcv[[1]], ylim = c(0.5,1.2),ylab="MSE",
     type = "b",col = "turquoise", pch = 20,main = "rep 10 times MSE")
lines(n.times_LM.kfoldcv[[1]], col = "orange", pch = 20,type = "b")
lines(n.times_GLM.kfoldcv[[1]], col = "purple", pch = 20,type = "b")
lines(n.times_GAM.kfoldcv[[1]], col = "dodgerblue", pch = 20,type = "b")
# "legend" : icon in plot 
legend("bottomright",pch=c(20,20,20),legend=c("Huber","LM","GLM","GAM"),
       col=c("turquoise","orange","purple","dodgerblue"),bty="n",box.lty = 0,cex=1)

# k fold RMSE
plot(n.times_WHM.kfoldcv[[2]], ylim = c(0.5,1.2),ylab="RMSE",
     type = "b",col = "turquoise", pch = 20,main = "rep 10 times RMSE")
lines(n.times_LM.kfoldcv[[2]], col = "orange", pch = 20,type = "b")
lines(n.times_GLM.kfoldcv[[2]], col = "purple", pch = 20,type = "b")
lines(n.times_GAM.kfoldcv[[2]], col = "dodgerblue", pch = 20,type = "b")
legend("bottomright",pch=c(20,20,20),legend=c("Huber","LM","GLM","GAM"),
       col=c("turquoise","orange","purple","dodgerblue"),bty="n",box.lty = 0,cex=1)

# k fold MAE
plot(n.times_WHM.kfoldcv[[3]], ylim = c(0.5,1.2),ylab="MAE",
     type = "b",col = "turquoise", pch = 20,main = "rep 10 timesMAE")
lines(n.times_LM.kfoldcv[[3]], col = "orange", pch = 20,type = "b")
lines(n.times_GLM.kfoldcv[[3]], col = "purple", pch = 20,type = "b")
lines(n.times_GAM.kfoldcv[[3]], col = "dodgerblue", pch = 20,type = "b")
legend("bottomright",pch=c(20,20,20),legend=c("Huber","LM","GLM","GAM"),
       col=c("turquoise","orange","purple","dodgerblue"),bty="n",box.lty = 0,cex=1)

