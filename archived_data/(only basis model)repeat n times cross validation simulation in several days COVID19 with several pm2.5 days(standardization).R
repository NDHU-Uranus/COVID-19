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
system.time({n.times_WHM.kfoldcv1 <- n.times_tempWHM.kfoldcv(repeat.time=10,fold=5,basis=10)})

system.time({n.times_WHM.kfoldcv2 <- n.times_tempWHM.kfoldcv(repeat.time=10,fold=5,basis=5)})

system.time({n.times_WHM.kfoldcv3 <- n.times_tempWHM.kfoldcv(repeat.time=10,fold=5,basis=3)})

system.time({n.times_WHM.kfoldcv4 <- n.times_tempWHM.kfoldcv(repeat.time=10,fold=4,basis=10)})

system.time({n.times_WHM.kfoldcv5 <- n.times_tempWHM.kfoldcv(repeat.time=10,fold=4,basis=5)})

system.time({n.times_WHM.kfoldcv6 <- n.times_tempWHM.kfoldcv(repeat.time=10,fold=4,basis=3)})

system.time({n.times_WHM.kfoldcv7 <- n.times_tempWHM.kfoldcv(repeat.time=10,fold=3,basis=10)})

system.time({n.times_WHM.kfoldcv8 <- n.times_tempWHM.kfoldcv(repeat.time=10,fold=3,basis=5)})

system.time({n.times_WHM.kfoldcv9 <- n.times_tempWHM.kfoldcv(repeat.time=10,fold=3,basis=3)})

#-----------merge all plot ----------------
par(mfrow = c(1, 3))
# k fold mse
plot(n.times_WHM.kfoldcv1[[1]],ylim = c(0.3,1.5),ylab="MSE",
     type = "b",col = "turquoise", pch = 20,main = "rep 10 times MSE")
lines(n.times_WHM.kfoldcv2[[1]], col = "orange", pch = 20,type = "b")
lines(n.times_WHM.kfoldcv3[[1]], col = "purple", pch = 20,type = "b")
lines(n.times_WHM.kfoldcv4[[1]], col = "green", pch = 20,type = "b")
lines(n.times_WHM.kfoldcv5[[1]], col = "red", pch = 20,type = "b")
lines(n.times_WHM.kfoldcv6[[1]], col = "grey", pch = 20,type = "b")
lines(n.times_WHM.kfoldcv7[[1]], col = "black", pch = 20,type = "b")
lines(n.times_WHM.kfoldcv8[[1]], col = "pink", pch = 20,type = "b")
lines(n.times_WHM.kfoldcv9[[1]], col = "dodgerblue", pch = 20,type = "b")
# "legend" : icon in plot 
#legend("topright",pch=c(20,20,20),legend=c("5f10b","5f5b","5f3b","4f10b","4f5b","4f3b","3f10b","3f5b","3f3b"),
#       col=c("turquoise","orange","purple","green","red","grey","black","pink","dodgerblue"),bty="n",box.lty = 0,cex=1)

# k fold RMSE
plot(n.times_WHM.kfoldcv1[[2]],ylim = c(0.3,1.5),ylab="RMSE",
     type = "b",col = "turquoise", pch = 20,main = "rep 10 times RMSE")
lines(n.times_WHM.kfoldcv2[[2]], col = "orange", pch = 20,type = "b")
lines(n.times_WHM.kfoldcv3[[2]], col = "purple", pch = 20,type = "b")
lines(n.times_WHM.kfoldcv4[[2]], col = "green", pch = 20,type = "b")
lines(n.times_WHM.kfoldcv5[[2]], col = "red", pch = 20,type = "b")
lines(n.times_WHM.kfoldcv6[[2]], col = "grey", pch = 20,type = "b")
lines(n.times_WHM.kfoldcv7[[2]], col = "black", pch = 20,type = "b")
lines(n.times_WHM.kfoldcv8[[2]], col = "pink", pch = 20,type = "b")
lines(n.times_WHM.kfoldcv9[[2]], col = "dodgerblue", pch = 20,type = "b")
#legend("topright",pch=c(20,20,20),legend=c("5f10b","5f5b","5f3b","4f10b","4f5b","4f3b","3f10b","3f5b","3f3b"),
#       col=c("turquoise","orange","purple","green","red","grey","black","pink","dodgerblue"),bty="n",box.lty = 0,cex=1)

# k fold MAE
plot(n.times_WHM.kfoldcv1[[3]],ylim = c(0.3,1.5),ylab="MAE",
     type = "b",col = "turquoise", pch = 20,main = "rep 10 timesMAE")
lines(n.times_WHM.kfoldcv2[[3]], col = "orange", pch = 20,type = "b")
lines(n.times_WHM.kfoldcv3[[3]], col = "purple", pch = 20,type = "b")
lines(n.times_WHM.kfoldcv4[[3]], col = "green", pch = 20,type = "b")
lines(n.times_WHM.kfoldcv5[[3]], col = "red", pch = 20,type = "b")
lines(n.times_WHM.kfoldcv6[[3]], col = "grey", pch = 20,type = "b")
lines(n.times_WHM.kfoldcv7[[3]], col = "black", pch = 20,type = "b")
lines(n.times_WHM.kfoldcv8[[3]], col = "pink", pch = 20,type = "b")
lines(n.times_WHM.kfoldcv9[[3]], col = "dodgerblue", pch = 20,type = "b")
legend("topright",pch=c(20,20,20),legend=c("5f10b","5f5b","5f3b","4f10b","4f5b","4f3b","3f10b","3f5b","3f3b"),
       col=c("turquoise","orange","purple","green","red","grey","black","pink","dodgerblue"),bty="n",box.lty = 0,cex=1.25)

