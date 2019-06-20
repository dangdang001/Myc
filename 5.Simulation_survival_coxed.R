# simulate survival data with the sim.survdata function

# Typical methods for generating survival often employ known distributions - such as the exponential, Weibull, or Gompertz 
# that imply specific shapes for the baseline hazard function. 
# This approach is potentially problematic because it contradicts a key advantage of the Cox model the ability to leave the distribution of the baseline hazard unspecified. 


# a unique baseline hazard by fitting a cubic spline to randomly-drawn points

#install.packages("coxed")
library(coxed)
n=1000
x1=rnorm(n,0,1)
x2=rnorm(n,0,1)
x3=runif(n,0,1)
x4=rbinom(n,1,0.5)
x5=rbinom(n,1,0.2)

df=data.frame(x1,x2,x3,x4,x5)
beta=c(1.2,0.5,0.1,0.5,0.6)


simdata <- sim.survdata(N=n, T=100, num.data.frames=1,beta=beta,X=df,xvars = 5,censor=0.2)

simdata <- sim.survdata(N=1000, T=100, num.data.frames=1)

head(simdata$data, 10)

hist(simdata$data$y)

survsim.plot(simdata, df=1, type="baseline")

