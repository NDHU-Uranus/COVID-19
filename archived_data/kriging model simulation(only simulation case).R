library(mvtnorm)
library(MASS)
library(fields)
library(plotly)
library(magrittr)
library(geoR)
library(autoFRK)
library(graphics)
library(Metrics)
library(ggplot2)
library(viridis)
library(plot3D)
library(plot3Drgl)
library(dplyr)
library(rgdal)
library(sp)
library(graphics)
library(raster)
library(gstat)
library(rgdal)
library(geoR)
library(tictoc)

#Simulation Case : random coordinary random valus(mean=linear)
#data1000 train00 test100 sigma^2=1 noise:N(0,1)

set.seed(788965)
mul<-function(n,mm,jj) {
  #random generate coords
  #mm,is sigma^2 ->mm*(exp(...))
  x <- runif(n,max = jj,min = 0)
  y <- runif(n,max = jj,min = 0)
  xxx<-as.data.frame(cbind(x,y)) 
  #covariance is stationary
  #ele_cov = stationary.cov(x1 = xxx,x2 =NULL, theta = 1, Distance = "rdist", Covariance="Exponential")
  #just covariance no any situation
  ele_cov = cov(xxx)
  asd =as.vector(t(rmvnorm(n/2, sigma = mm*ele_cov)))  
  return(list(x=xxx,data=asd,ele_cov,jj)) 
}
# to plot multivariate normal distribution in ggplot
# ggplot(m, aes(x=V1, y=V2))+
#   geom_point(alpha = .2) +
#   geom_density_2d()+
#   theme_bw()

multi_data<-mul(1000,1,10)

jj<-multi_data[[4]]
covar<- multi_data[[2]]
multicoor <- multi_data[[1]] 
grid_po <-multicoor

##---------observation = linear function (mean tern) + covariance tern
# (linear mean tern)
simean<-function(n,x,y){
  # z1=3*x+4*y
  err<- rnorm(n, mean = 0, sd = 1)
  #1.  z2 = 10 + 2*x + 0.7*y^2  +  err
  # select function
  #2.z3 = 0.75/pi*0.3*0.4*exp(-(x-0.2)^2/(0.3)^2)-(y-0.3)^2/0.4^2 + err
  #3.z4 = 1.9*(1.45+exp(x)*sin(13*(x-0.6)^2))*exp(-y)*sin(7*y) + err 
  z5 = exp((-(x-0.25)^2-(y-0.25)^2)/0.1)+0.5*exp((-(x-0.7)^2-(y-0.7)^2)/0.07)
  # sim_mean<-cbind(z1,z2,z3)
  
  #  return(sim_mean)
  return(z5)
}

simumean <- simean(length(multicoor[,1]),multicoor[1],multicoor[2])

muldata <- simumean + covar

# random value
#multidata <- as.data.frame(sample(muldata)) 
multidata <- as.data.frame(muldata) 

summary(multidata)

#merge coordinates & value
gridpo<-cbind(grid_po,multidata)
colnames(gridpo)[3] <- "value"
# color plot
colop<-ggplot(gridpo, aes(x, y))+
  geom_point(aes(color = value)) +
  scale_color_viridis(option = "H")
colop

###-------diving train and test
#random draw the grid point(with sample)
# This random sample size will equal to line 1248(identity matrix:diag)
nn <- 900 
rangridpo<-sample(1:nrow(gridpo), nn)
#simudata:simulation data(from random draw)
simudata<-gridpo[rangridpo,]
plot(grid_po,cex = 2 ,lwd = 2,main = 'Regular points with Random selected')
# random pick the points and marked in red 
points(simudata[1:2],cex = 2,lwd = 2, col = "red")
# other points are not be select will marked in blue 
simu_grid_poin<-gridpo[-rangridpo,]
points(simu_grid_poin[1:2],cex = 2,lwd = 2, col = "blue")

realdot<-simu_grid_poin%>% as.data.frame %>% 
  ggplot(aes(x, y)) + geom_point(aes(colour = value)) + 
  labs(title="real dot") + 
  scale_color_viridis(option = "H")
realdot

###----------number of basis = 5
# use simulation location in mrts to generate basis
bbaa <- 5 # number of basis
simubas <- as.matrix(mrts(knot = simudata[1:2], bbaa))

###------------Huber M estimate to estimate β
#------Huber----------#

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

#oc:outcome
oc<- IRLS(X=simubas,Y=simudata[,3],max.iter=100,conv.eps=1e-10)
oc[[1]]
irls.simresidu <-oc[[2]]

#covariavce data = (observation data) - mean
covadata<- simudata[,3] - simubas%*%oc[[1]]
plot(covadata, main = "Mean residual plot")
hist(covadata, main = "Mean residual histogram")
#residual normal test if p-value > 0.05,residual consistent with normal distribution
shapiro.test(covadata)
# ASD plot(A Standard Deviation)(Cumulative Distribution):
# the greater the area under the curve, the smaller the residual!
#As you go from 0 to the right on the X-axis, the curve will pull up very quickly,
#which means that most of the residuals are within a very small standard deviation
std <- sd(simudata$value) # Standard Deviation with real value
plot(ecdf(abs(covadata/std)),main="ASD mean reseidual chart",xlab="Residual/STD",ylab="Cumulative distrubtion",col="pink")
summary(covadata)
#=================once variogram -> to eatimate nugget effect==============#
# uvec : a vector with values used to define the variogram binning
#variogram()
#max.dist : length of variogram
#fix.nug=F ->means have nugget effect
simuvario.b <- variog(coords = simudata[1:2], data = covadata,uvec=seq(0,jj, l=35), max.dist=6)
plot(simuvario.b, main="binned variogram") 
 # sigma^2 (partial sill) and phi (range parameter)
simuini.vals <- expand.grid(seq(6,simuvario.b$var.mark,l=jj), seq(0.1,jj,l=jj))
simuols <- variofit(simuvario.b,fix.nug=F, wei="equal",cov.model="exp", ini=simuini.vals)
lines(simuols,lty=3,col=3,lwd=2)

#weighgt = 'cressie' or 'npairs'
simuini.valss <- expand.grid(seq(6,simuols$cov.pars[1],l=jj), seq(0.1,jj,l=jj)) # sigma^2 (partial sill) and phi (range parameter)
simuwls <- variofit(simuvario.b, ini=simuini.valss, fix.nug=F,wei="npairs",cov.model="exp")
lines(simuwls, lty=2,col=4,lwd=2)

Name=c("OLS","WLS")
legend("bottomright", Name, ncol = 1, cex = 1, col=c(3,4),
       lty = c(3,2), lwd = c(2,2), bg = 'gray95')
 #=================twice variogram -> deduct eatimate nugget effect==============#
### residual data deduct(去除) nugget effect ,and estimate again
# NOnugg_covadata <- covadata - simuols$nugget
# simuvario.b <- variog(coords = simudata[1:2], data = NOnugg_covadata,uvec=seq(0,jj, l=35), max.dist=6)
# plot(simuvario.b, main="binned variogram") 
# 
# simuini.vals <- expand.grid(seq(0.1,simuvario.b$var.mark,l=jj), seq(0.1,jj,l=jj))
# simuools <- variofit(simuvario.b,fix.nug=T, wei="equal",cov.model="exp", ini=simuini.vals)
# lines(simuools,lty=3,col=3,lwd=2)
# 
# simuini.valss <- expand.grid(seq(0.1,simuools$cov.pars[1],l=jj), seq(0.1,jj,l=jj)) # sigma^2 (partial sill) and phi (range parameter)
# simuwwls <- variofit(simuvario.b, ini=simuini.valss, fix.nug=T,wei="npairs",cov.model="exp")
# lines(simuwwls, lty=2,col=4,lwd=2)


###--------use kriging model  (weight=wls) and see the model fit (MSE)
simuweight <- simuwls
# mean function residu
# τ^2 + σ^2, corresponds to the variance of the observation process Y
simresidu<-as.matrix(irls.simresidu)
TPS_simu_beta_hat <- oc[[1]]
# simu_b_grid:grid point(simu_grid_poin) into TPS basis
simu_b_grid <- as.matrix(mrts(knot = simu_grid_poin[1:2], bbaa))
simu_phi_s0 <- t(t(TPS_simu_beta_hat) %*% t(simu_b_grid))

simu_c_s <-c()
# report running time :system.time
running_time <- system.time({
for (i in 1:length(simudata[,c("x")])){
  for (j in 1:length(simu_grid_poin[,c("x")])){
    # if(i==j){next}
    #else{
    simu_disten<-sqrt((simudata[,c("x")][i]-simu_grid_poin[,c("x")][j])^2+(simudata[,c("y")][i]-simu_grid_poin[,c("y")][j])^2) 
    simu_c_s <- c(simu_c_s ,simu_disten)
    # }
  }
  if(i%%50 == 0) cat("running: ", i, "\n")
}
                        })
running_time

simu_c_s0 <- simuweight$cov.pars[1]*exp(-abs(matrix(simu_c_s,nrow = nrow(simu_grid_poin))/(simuweight$cov.pars[2]+1e-10)))

simu_zigma_theta <- simuweight$cov.pars[1]*as.matrix(exp(-dist(simudata[-3], method = "euclidean", diag = F, upper = FALSE, p = 2))/(simuweight$cov.pars[2]+1e-10)) #+1e-10: prevent simu_zigma_theta = 0 or Inf
#use y to deduct simu_c_s0  %*% solve(simu_zigma_theta +.....,and we can get the
#estimation about mean part
#simu_mean <- simudata[,3] -  

simu_y_hat =  simu_phi_s0 + simu_c_s0  %*% solve(simu_zigma_theta + (simuweight$nugget+simuweight$cov.pars[1])*diag(nn)) %*% simresidu
# grid point & y_hat
simu_krig_data <- cbind(simu_grid_poin[1:2],simu_y_hat)

# check MSE & resudual
simu_MSE = mse(simu_grid_poin$value, simu_krig_data$simu_y_hat)
cat("simulation MSE:",simu_MSE) 
simu_residual <- simu_grid_poin$value - simu_krig_data$simu_y_hat
summary(simu_residual)

# Compare y,y hat,residual
simu_comptable <-cbind(simu_grid_poin$value,simu_krig_data$simu_y_hat,simu_residual)
head(simu_comptable,10)

# predict dot map vs real dot map
predot<-simu_krig_data%>% as.data.frame %>% 
  ggplot(aes(x, y)) + geom_point(aes(colour = simu_y_hat)) + 
  labs(title="pre dot") +
  scale_color_viridis(option = "H")

dotvs <-subplot(predot,realdot)%>% 
  layout(title = 'predict   &    real')
dotvs

# 3D in real value---1
scatter3D(simu_grid_poin$x, simu_grid_poin$y, simu_grid_poin$value,
          zlab="value",phi = 10, theta = 45,main="real value Plot")
# 3D plot in estimate value---1
scatter3D(simu_krig_data$x, simu_krig_data$y, simu_krig_data$simu_y_hat,
          zlab="value",phi = 10, theta = 45,main="estimate value Plot")
#contour-plot
fig3<- plot_ly(z = ~volcano, type = "contour")

