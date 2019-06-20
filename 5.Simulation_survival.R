# Simulating survival time with time-invariant covarites 

n=1000

x1=rnorm(n,0,1)
x2=rnorm(n,0,1)
x3=runif(n,0,1)
x4=rbinom(n,1,0.5)
x5=rbinom(n,1,0.2)

df=data.frame(x1,x2,x3,x4,x5)

# 1. Exponential

simu_exp<-function(x,lambda){
    u=runif(1,0,1)
    t=-log(u)/(lambda*exp(1.2*x[1]+0.5*x[2]+0.1*x[3]+0.5*x[4]+0.6*x[5]))
    return(t)
}

surv_exp=apply(df,1,simu_exp,lambda=2)

hist(surv_exp)



# 2. Weibull

simu_weibull<-function(x,lambda,v){
  u=runif(1,0,1)
  t=-(log(u)/(lambda*exp(1.2*x[1]+0.5*x[2]+0.1*x[3]+0.5*x[4]+0.6*x[5])))^(1/v)
  return(t)
}

surv_weibull=apply(df,1,simu_weibull,lambda=1,v=1)

hist(surv_weibull)


# 3. Gompertz distribution

simu_gompertz<-function(x,lambda,alpha){
  u=runif(1,0,1)
  t=1/alpha*log(1-(alpha*log(u))/(lambda*exp(1.2*x[1]+0.5*x[2]+0.1*x[3]+0.5*x[4]+0.6*x[5])))
  return(t)
}

surv_gompertz=apply(df,1,simu_gompertz,lambda=1,alpha=2)

hist(surv_gompertz)


df$surv_exp=surv_exp
df$surv_weibull=surv_weibull
df$surv_gompertz=surv_gompertz

df$status=1

plot(surv_gompertz)




# Fit a Cox proportional hazards model

library(survival)

fit_exp <- coxph(Surv(time = df$surv_exp, event = df$status) ~ x1+x2+x3+x4+x5, data = df)
fit_weibull <- coxph(Surv(time = df$surv_weibull, event = df$status) ~ x1+x2+x3+x4+x5, data = df)
fit_gompertz <- coxph(Surv(time = df$surv_gompertz, event = df$status) ~ x1+x2+x3+x4+x5, data = df)

fit_exp
fit_weibull
fit_gompertz 
