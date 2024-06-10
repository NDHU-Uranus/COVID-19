#--------------example-----------#
require(ggplot2)
require(sandwich)
require(msm)
library(MASS)
#The number of awards earned by students at one high school. 
#Predictors of the number of awards earned include the type of program in which the student was enrolled 
#(e.g., vocational, general or academic) and the score on their final exam in math.
require(ggplot2)
require(sandwich)
require(msm)
p <- read.csv("https://stats.idre.ucla.edu/stat/data/poisson_sim.csv")
p <- within(p, {
  prog <- factor(prog, levels=1:3, labels=c("General", "Academic", 
                                            "Vocational"))
  id <- factor(id)
})
summary(p)

with(p, tapply(num_awards, prog, function(x) {
  sprintf("M (SD) = %1.2f (%1.2f)", mean(x), sd(x))
}))

ggplot(p, aes(num_awards, fill = prog)) +
  geom_histogram(binwidth=.5, position="dodge")

m1 <- glm.nb(num_awards ~ prog, data = p)
m2 <- glm(num_awards ~ prog, family=poisson(link = log), data=p)
summary(m2)
summary(m1)
with(m1, cbind(res.deviance = deviance, df = df.residual,
               p = pchisq(deviance, df.residual, lower.tail=FALSE)))

with(m2, cbind(res.deviance = deviance, df = df.residual,
               p = pchisq(deviance, df.residual, lower.tail=FALSE)))

# Visualize the data set
x <- as.matrix(covid1128_to_sub2[,2:3])
y <- as.matrix(covid1128_to_sub2[4])

quilt.plot(x,y)
US(add=TRUE)
# thin plate spline-like model with the lambda parameter estimated by
# maximum likelihood. Default choices are made for a.wght, nlevel, NC
# and alpha.

obj<- LatticeKrig( x, y)
# }
# NOT RUN {
# summary of fit and a plot of fitted surface
print( obj)
surface( obj )
US(add=TRUE)
points(x)
