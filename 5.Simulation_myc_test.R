# K: number of myc-target genes
# n: number of patients


simu_myc<-function(n,K,mean,myc,size){
  df<-NULL
  for(k in 1:K){
    scale=runif(1,0,4)
    effect=myc*scale
    df=cbind(df,rnbinom(n, mu = mean*effect, size = size))
  }
  return(df)
}

df=simu_myc(n=1000,K=100,mean=5,myc=1,size=1)


# 1. Exponential


n=1000

x1=rnorm(n,0,1)
x2=rnorm(n,0,1)
x3=runif(n,0,1)
x4=rbinom(n,1,0.5)
x5=rbinom(n,1,0.2)
myc=rnorm(n,0,1)

df=data.frame(x1,x2,x3,x4,x5,myc)
coef=c(1,0.5,0.2,0.8,1.5,5)

simu_exp<-function(x,lambda){
  u=runif(1,0,1)
  t=-log(u)/(lambda*exp(coef %*% x))
  return(t)
}

surv_exp=apply(df,1,simu_exp,lambda=2)

hist(surv_exp)



# 2. Weibull

simu_weibull<-function(x,lambda,v){
  u=runif(1,0,1)
  t=-(log(u)/(lambda*exp(coef %*% x)))^(1/v)
  return(t)
}

surv_weibull=apply(df,1,simu_weibull,lambda=1,v=1)

hist(surv_weibull)


# 3. Gompertz distribution

simu_gompertz<-function(x,lambda,alpha){
  u=runif(1,0,1)
  t=1/alpha*log(1-(alpha*log(u))/(lambda*exp(coef %*% x)))
  return(t)
}

surv_gompertz=apply(df1,1,simu_gompertz,lambda=1,alpha=1)

hist(surv_gompertz)


df1$surv_gompertz=surv_gompertz

df1$status=1


simu_myc<-function(n,K,mean,myc,size){
  df<-NULL
  for(k in 1:K){
    scale=runif(1,0,4)
    effect=exp(myc)*scale
    temp=c()
    for(i in 1:n){
      cur=rnbinom(1, mu = mean*effect[i], size = size)
      temp=c(temp,cur)
    }
    df=cbind(df,temp)
    colnames(df)[k]<-paste0('g',k)
  }
  df=as.data.frame(df)
  return(df)
}

df2=simu_myc(n=1000,K=10,mean=5,myc=myc,size=1)

dim(df1)

library(survival)


fit_gompertz <- coxph(Surv(time = df$surv_gompertz, event = df$status) ~ x1+x2+x3+x4+x5+g1+g2+g3+g4+g5+g6+g7+g8+g9+g10, data = df)

fit_gompertz 


library(rms)

f <- cph(Surv(time = df$surv, event = df$status) ~ x1+x2+x3+x4+x5+myc, x=TRUE, y=TRUE,surv = TRUE)

f$stats['Dxy']/2+0.5

validate(f, B=100)

cal <- calibrate(f, u=1.5, cmethod='KM', cuts=seq(0,1,0.1), B=20)  # usually B=200 or 300
plot(cal, add=pa)