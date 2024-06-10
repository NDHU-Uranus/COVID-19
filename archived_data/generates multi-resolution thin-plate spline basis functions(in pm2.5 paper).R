library(autoFRK)
originalPar <- par(no.readonly = TRUE)
knot <- seq(0, 1, l = 25)
b <- mrts(knot, 25)
x0 <- seq(0, 1, l = 200)
bx <- predict(b, x0)
par(mfrow = c(5, 5), mar = c(0, 0, 0, 0))
for (i in 1:25) {
  plot(bx[, i], type = "l", axes = FALSE)
  box()
}
par(originalPar)


# ----------------mean function estimate----------------#
#Standardization
stan_e<-function(ei){
  sigma<-median(abs(ei-median(ei)))/0.6745
  e_stan<-ei/sigma   
  return(e_stan)
}

new_residu<-function(beta){
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

#par(mfrow=c(2,2))
# outcome(original data)
outcome1 <- IRLS(X=X,Y=Y,max.iter=300,conv.eps=1e-10)
outcome1[[1]]
irls.residu1 <-outcome1[[2]]
qqnorm(irls.residu1)
qqline(irls.residu1, col = "steelblue", lwd = 2)

fitt_valu1 <- X %*% outcome1[[1]]
plot(x = fitt_valu1 , y = irls.residu1,main="Residual vs fitted(original)")




#----------mrts struction--------------------#
mrts<-function (knot, k, x = NULL, maxknot = 5000) 
{
  # to make sure your R is in 64-bit
  is64bit <- length(grep("64-bit", sessionInfo()$platform)) > 
    0
  if ((!is64bit) & (max(NROW(x), NROW(knot)) > 20000)) {
    stop("Use 64-bit version of R for such a volume of data!")
  }
  # number of columns 
  if (NCOL(knot) == 1) {
    xobs <- as.matrix(as.double(as.matrix(knot)))
  }
  else {
    xobs <- apply(knot, 2, as.double)
  }
  #uniquecombs:find the unique rows in a matrix
  Xu <- uniquecombs(cbind(xobs))
  if (is.null(x) & length(Xu) != length(xobs)) {
    x <- xobs
  }
  colnames(Xu) <- NULL
  n <- n.Xu <- NROW(Xu)
  ndims <- NCOL(Xu)
  if (k < (ndims + 1)) {
    stop("k-1 can not be smaller than the number of dimensions!")
  }
  if (maxknot < n) {
    bmax <- maxknot
    Xu <- subKnot(Xu, bmax)
    Xu <- as.matrix(Xu)
    if (is.null(x)) {
      x <- knot
    }
    n <- NROW(Xu)
    n.Xu <- n
  }
  xobs_diag <- diag(sqrt(n/(n - 1))/apply(xobs, 2, sd), ndims)
  if (!is.null(x)) {
    if (NCOL(x) == 1) {
      x <- as.matrix(as.double(as.matrix(x)))
    }
    else {
      x <- as.matrix(array(as.double(as.matrix(x)), dim(x)))
    }
    if (k - ndims - 1 > 0) {
      result <- mrtsrcpp_predict0(Xu, xobs_diag, x, k - 
                                    ndims - 1)
    }
    else {
      X2 <- scale(Xu, scale = FALSE)
      shift <- colMeans(Xu)
      nconst <- sqrt(diag(t(X2) %*% X2))
      X2 <- cbind(1, t((t(x) - shift)/nconst) * sqrt(n))
      result <- list(X = X2[, 1:k])
      x <- NULL
    }
  }
  else {
    if (k - ndims - 1 > 0) {
      result <- mrtsrcpp(Xu, xobs_diag, k - ndims - 1)
    }
    else {
      X2 <- scale(Xu, scale = FALSE)
      shift <- colMeans(Xu)
      nconst <- sqrt(diag(t(X2) %*% X2))
      X2 <- cbind(1, t((t(Xu) - shift)/nconst) * sqrt(n))
      result <- list(X = X2[, 1:k])
    }
  }
  obj <- result$X
  if (is.null(result$nconst)) {
    X2 <- scale(Xu, scale = FALSE)
    result$nconst <- sqrt(diag(t(X2) %*% X2))
  }
  attr(obj, "UZ") <- result$UZ
  attr(obj, "Xu") <- Xu
  attr(obj, "nconst") <- result$nconst
  attr(obj, "BBBH") <- result$BBBH
  attr(obj, "class") <- c("matrix", "mrts")
  class(obj) <- "mrts"
  if (is.null(x)) {
    return(obj)
  }
  else {
    shift <- colMeans(attr(obj, "Xu"))
    X2 <- sweep(cbind(x), 2, shift, "-")
    X2 <- cbind(1, sweep(X2, 2, attr(obj, "nconst"), 
                         "/"))
    if (k - ndims - 1 > 0) {
      obj0 <- as.matrix(cbind(X2, result$X1))
    }
    else {
      obj0 <- as.matrix(X2)
    }
    dimnames(obj) <- NULL
    aname <- names(attributes(obj))
    attributes(obj0) <- c(attributes(obj0), attributes(obj)[setdiff(aname, 
                                                                    c("dim", "dimnames"))])
    return(obj0)
  }
}
