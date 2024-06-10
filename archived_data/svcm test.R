
library(varycoef)
library(sp)
help("varycoef")

set.seed(123)

# attach sp and load data
library(sp)
data("meuse")

# documentation
help("meuse")

# overview
summary(meuse)

par(mfrow = c(1, 1))

meuse$l_cad <- log(meuse$cadmium)

# construct spatial object
sp.meuse <- meuse
coordinates(sp.meuse) <- ~x+y
proj4string(sp.meuse) <- CRS("+init=epsg:28992")

# using package tmap
library(tmap)
# producing an interactive map
tmap_leaflet(tm_shape(sp.meuse) + tm_dots("l_cad", style = "cont"))

# linear model (LM)
lm.fit <- lm(l_cad ~ 1+dist+lime+elev, data = meuse)
summary(lm.fit)

####svc model-----------------
# response variable
y <- meuse$l_cad
# covariates for fixed effects
X <- model.matrix(~1+dist+lime+elev, data = meuse)
# locations
locs <- as.matrix(meuse[, 1:2])/1000

# covariates for SVC (random effects)
W <- model.matrix(~1+dist+lime, data = meuse)

# construct initial value (recall transformation from meters to kilometers)
(init <- c(
  # 3 times for 3 SVC
  rep(c(
    # range
    fV$range[2]/1000,
    # variance
    fV$psill[2]),     
    3), 
  # nugget
  fV$psill[1]
))

# control settings vor MLE
control <- SVC_mle_control(
  # profile likelihood optimization
  profileLik = TRUE,
  # initial values
  init = init
)

# MLE
VC.fit <- SVC_mle(y = y, X = X, W = W, locs = locs,
                  control = control)
# outcome
summary(VC.fit)

# residuals
oldpar <- par(mfrow = c(1, 2))
plot(VC.fit, which = 1:2)

par(mfrow = c(1, 1))
plot(VC.fit, which = 3)

par(oldpar)

newlocs <- coordinates(meuse.grid)/1000
# prediciton
VC.pred <- predict(VC.fit, newlocs = newlocs)
# outcome
head(VC.pred)

# transformation to a spatial points data frame
colnames(VC.pred)[1:3] <- c("Intercept", "dist", "lime")
sp.VC.pred <- VC.pred
coordinates(sp.VC.pred) <- coordinates(meuse.grid)
proj4string(sp.VC.pred) <- proj4string(sp.meuse)

# points of interest (POI), we come back to them later
POI1.id <- 235
POI2.id <- 2016
POI1 <- meuse.grid[POI1.id, ]
POI2 <- meuse.grid[POI2.id, ]

tm_POIs <- 
  tm_shape(POI1) +
  # square (pch = 15)
  tm_symbols(shape = 15, col = "black") +
  tm_shape(POI2) +
  # triangle (pch = 17)
  tm_symbols(shape = 17, col = "black")

# tm.GS: GS fit
tm.GS <- tm_shape(GS.fit) + tm_dots("var1.pred", style = "cont")

# tm1: Intercept
tm1 <- tm_shape(sp.VC.pred) + 
  tm_dots("Intercept", style = "cont") +
  tm_POIs
# tm2: dist
tm2 <- tm_shape(sp.VC.pred) + 
  tm_dots("dist", style = "cont") +
  tm_POIs
# tm1: Intercept
tm3 <- tm_shape(sp.VC.pred) + 
  tm_dots("lime", style = "cont") +
  tm_POIs

tmap_arrange(list(tm.GS, tm1, tm2, tm3), ncol = 2, nrow = 2)

coef(VC.fit)
