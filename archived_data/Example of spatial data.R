library(MASS)
cov.matern=function(x, nu = 2, alpha = 1, vars=1)
{
  if(nu == 0.5)
    return(vars*exp( - x * alpha))
  ismatrix <- is.matrix(x)
  if(ismatrix){nr=nrow(x); nl=ncol(x)}
  x <- c(alpha * x)
  output <- rep(1, length(x))
  n <- sum(x > 0)
  if(n > 0) {
    x1 <- x[x > 0]
    output[x > 0] <-
      (1/((2^(nu - 1)) * gamma(nu))) * (x1^nu) * besselK(x1, nu)
  }
  if(ismatrix){
    output <- matrix(output, nr, nl)
  }
  vars*output
}
locs=cbind(rep(0:20, 21)/20, rep(0:20, each=21)/20)
V=cov.matern(as.matrix(dist(locs)),nu=1/2, alpha=7*sqrt(3))
set.seed(20)
z=mvrnorm(mu=rep(0, 21^2), Sigma=V)
z=matrix(z, ncol=21)
persp(x=(0:20)/21, y=(0:20)/21, z, theta=45, phi=35, r=5, expand=0.6, axes=T,
      ticktype="detailed", xlab="x", ylab="y", zlab="z")
filled.contour(x=0:20, y=0:20, z, color.palette=gray.colors)
V=cov.matern(as.matrix(dist(locs)),nu=2, alpha=7*sqrt(3))
set.seed(20)
z=mvrnorm(mu=rep(0, 21^2), Sigma=V)
z=matrix(z, ncol=21)
