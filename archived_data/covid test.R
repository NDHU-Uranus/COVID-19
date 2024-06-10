#-------------------covid-19 US data  use in glm.nb & glm & lm-----------------------#
require(ggplot2)
require(sandwich)
require(msm)
library(MASS)
library(magrittr)
library(ggmap)
require(car)
library(tidyverse)

covid1231_to <- read.csv("D:\\研究所Meeting\\COVID19\\csse_covid_19_data\\csse_covid_19_daily_reports_us\\12-31-2020(修改版).csv")
popula_2019 <- read.csv("D:\\研究所Meeting\\COVID19\\populations\\nst-est2019-01(修改版).csv") %>% as.data.frame()
house_income <- read.csv("D:\\研究所Meeting\\COVID19\\house income\\Median Household Income by State(已修改).csv") %>% as.data.frame()
covid1231_to[is.na(covid1231_to)] <- 0
covid1231_to_sub <- covid1231_to[c(-2,-3,-13,-14,-38),]

#request data from large to small
covid1231_to_sub1 <- covid1231_to[c(-2,-3,-13,-14,-38),c(-2,-3,-5,-6,-7)] %>%
  arrange(desc(Confirmed))

head(covid1231_to_sub1)
summary(covid1231_to_sub1)

# merge data
covid1231_to_sub_p <- merge(covid1231_to_sub,popula_2019,
               by.x = 'Province_State',
               by.y = 'state')

covid1231_to_sub_p_i <- merge(covid1231_to_sub_p,house_income,
                              by.x = 'Province_State',
                              by.y = 'State')

summary(covid1231_to_sub_p_i)

# creat model
m1 <- glm.nb(Confirmed ~ Bed + population + income, data = covid1231_to_sub_p_i)
summary(m1)
with(m1, cbind(res.deviance = deviance, df = df.residual,
               covid1231_to_sub_glm = pchisq(deviance, df.residual, lower.tail=FALSE)))

m2 <- glm(Confirmed ~ Bed + population + income, data = covid1231_to_sub_p_i)
summary(m2)
with(m2, cbind(res.deviance = deviance, df = df.residual,
               covid1231_to_sub_glm = pchisq(deviance, df.residual, lower.tail=FALSE)))

m3 <- lm(Confirmed~Bed + population + income, data=covid1231_to_sub_p_i)
summary(m3)


#Normal test, nihility hypothesis H0: the residual is subject to normal distribution
#because the p-value > 0.05 means that H0 will not be rejected
shapiro.test(m1$residual)
qqplot(m1$residual)
qqnorm(m1$residual)
qqline(m1$residual)
#Independent residual test
durbinWatsonTest(m1) 

# simular with other model(m2 & m3)
shapiro.test(m2$residual)
qqnorm(m2$residual)
qqline(m2$residual)
durbinWatsonTest(m2) 

shapiro.test(m3$residual)
qqnorm(m3$residual)
qqline(m3$residual)
durbinWatsonTest(m3) 


#draw the US map
s <- map_data('state')
ggplot(s,aes(x = long, y = lat, group = group)) + 
  geom_polygon(color = 'pink') +
  coord_map('polyconic') 
#changing to lower case
covid1231_to_sub1$Province_State <- tolower(covid1231_to_sub1$Province_State)
#merge data
dataa <- merge(s,covid1231_to_sub1,
               by.x = 'region',
               by.y = 'Province_State')
ggplot(dataa,aes(x = long, y = lat,
                 group = group,
                 fill = Confirmed)) + 
  geom_polygon(color = 'gray') +
  coord_map('polyconic') +
  scale_fill_gradient2()      
        
                                                                                                                                                                                                                                                                          
#### kriging
library(fields)
library(LatticeKrig)
us_lon_lat <- as.matrix(covid1128_to_sub[,2:3])
us_confir <- covid1128_to_sub[,4]
obj_us<- LatticeKrig(x = us_lon_lat, y = us_confir)
print(obj_us)
# and a plot of fitted surface
surface( obj_us )
US(add=TRUE)
points(us_lon_lat)


# ------------------------------- use pm2.5 model ------------------------------- #

# = = = = = = = =  Huber's M-estimator = = = = = = = = #
X <- covid1231_to_sub_p_i[,c(7,8,9)] %>% as.matrix()
Y <- covid1231_to_sub_p_i[,4] %>% as.matrix()

#Standardization
stan_e<-function(ei){
  sigma<-median(abs(ei-median(ei)))/0.6745
  e_stan<-ei/sigma   
  return(e_stan)
}

# new error from iterative
new_residu<-function(beta){
  tt<- X %*% beta
  new_res <- as.matrix(Y-tt)
  
  return(new_res)
}

# = = = = weight function(second differential) = = = = #
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

#initial value (least square)
b0 = solve(t(X) %*% X) %*% t(X) %*% Y
b0
lsmean_fun<- X %*% b0
ls.residu=as.matrix(Y-lsmean_fun)
#check transpose
#dd<-as.data.frame(t(X))

# = = = weighted least squares = = = = #
WLS<-function(X,Y,W){
  w_b = solve(t(X) %*% W %*% X) %*% t(X) %*% W %*% Y
  return(w_b)
}

#i=1
w <-0L
for (i in 1:length(ls.residu)) {
  w[i] <- 1/new_residu(b0)[i]
}
ww<-diag(w)
b1=WLS(X=X,Y=Y,W=ww)
b1
wlsmean_fun<- X %*% b1
wls.residu=as.matrix(Y-wlsmean_fun)


# = = = iterative reweighted least squares = = = = #

IRLS<-function(X,Y,max.iter,conv.eps){
  beta_init <- solve(t(X) %*% X) %*% t(X) %*% Y #initialize beta
  beta_prev <- beta_init               #beta_{t-1} (for comparisons)
  mean_fun<- X %*% beta_init
  e_i.prev<- as.matrix(Y-mean_fun)
  
  for(ii in 1:max.iter){
    standar_ei <- stan_e(new_residu(beta_prev)) 
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

outcome <- IRLS(X=X,Y=Y,max.iter=300,conv.eps=1e-10)
outcome[[1]]
irls.residu <-outcome[[2]]
qqnorm(irls.residu)
qqline(irls.residu)
#durbinWatsonTest(outcome) 

#par(new=TRUE)
par(mfrow = c(2,2))
#,ylim = c(min(ls.error),max(ls.error))
plot(ls.residu,type = "l" ,col="blue",lwd=1)
plot(wls.residu,type = "l" ,col="red",lwd=1)
plot(irls.residu,type = "l" ,col="orange",lwd=1)
# legend("topright",                         
#        pch = 1,cex=1.5, lwd=6,                              
#        col = c("blue", "red", "green"),          
#        legend = c("LS", "WLS", "IRLS"),           
#        lty = c(1,1),x.intersp=0.05,y.intersp=0.15,bty="n"
# )
#real=c(-26,23,-6.99,3.254,8)
#out_table <-as.data.frame(cbind(b0,b1,outcome[[1]],real))
#names(out_table) <- c("LS", "WLS", "IRLS","real")
out_table <-as.data.frame(cbind(b0,b1,outcome[[1]]))
names(out_table) <- c("LS", "WLS", "IRLS")
out_table
