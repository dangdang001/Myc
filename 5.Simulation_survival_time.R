# Name: Donglei Yin
# Date: 2019-06-22
# Purpose: Simulation of survival time and censoring based on Cox PH model


library(ggplot2)
library(gridExtra)

#n:number of patients/rows
#K:number of myc-target genes
#lambda:parameter of exponential distribution
#coef_myc: relative coeffient magnitude of myc effect in the true model
#theta: censoring distribution parameter, C follows U(0,theta)

#theta=100000: censoring rate=0
#theta=200: censoring rate=0.1
#theta=100: censoring rate=0.2
#theta=60: censoring rate=0.3
#theta=30: censoring rate=0.5
#theta=20: censoring rate=0.6
#theta=10: censoring rate=0.7
#theta=1: censoring rate=0.9

n=1000
K=10
coef_myc=1
  
# 1. Simulate survival time using Cox PH model based on exponential distribution as basedline hazard distribution

# baseline covariates: x1,x2,x3,x4,x5,myc

x1=rnorm(n,0,1)
x2=rnorm(n,0,1)
x3=rnorm(n,0,1)
x4=rbinom(n,1,0.2)
x5=rbinom(n,1,0.5)
myc=rnorm(n,0,1)

df=data.frame(x1,x2,x3,x4,x5,myc)

coef=c(1,1,1,1,1,coef_myc) # coeffient of baseline covariates: x1,x2,x3,x4,x5,myc

# use inverese cdf formula to generate survival time

simu_exp<-function(x,lambda){
   u=runif(1,0,1)
   t=-log(u)/(lambda*exp(coef %*% x))
   return(t)
}

simu_weibull<-function(x,lambda,v){
  u=runif(1,0,1)
  t=-(log(u)/(lambda*exp(coef %*% x)))^(1/v)
  return(t)
}

simu_gompertz<-function(x,lambda,alpha){
  u=runif(1,0,1)
  t=1/alpha*log(1-(alpha*log(u))/(lambda*exp(coef %*% x)))
  return(t)
}

#true survival time

lambda=2.5
surv_exp=apply(df,1,simu_exp,lambda)

lambda=4;v=1
surv_weibull=apply(df,1,simu_weibull,lambda,v)

lambda=0.01;alpha=0.1
surv_gompertz=apply(df,1,simu_gompertz,lambda,alpha)

df$surv_exp=surv_exp
df$surv_weibull=surv_weibull
df$surv_gompertz=surv_gompertz



h1<-ggplot(df, aes(surv_exp)) + 
  geom_histogram(bins = 20,colour="white",fill = "lightblue") +
  theme_bw()+
  labs(x = "Exponential",y="Frequency")

h2<-ggplot(df, aes(surv_weibull)) + 
  geom_histogram(bins = 20,colour="white",fill = "lightblue") +
  theme_bw()+
  labs(x = "Weibull",y="Frequency")

h3<-ggplot(df, aes(surv_gompertz)) + 
  geom_histogram(bins = 20,colour="white",fill = "lightblue") +
  theme_bw()+
  labs(x = "Gompertz",y="Frequency")

h1
h2
h3


h=grid.arrange(h1,h2,h3, ncol=3)

ggsave("C:/Users/Donglei/Desktop/Simulation_result//histogram_no_censoring.png", h, height=6, width=15, units='in', dpi=600)



# 2. Simulating observed time

# simulate censoring time based on uniform distribution

#theta=100000: censoring rate=0
#theta=200: censoring rate=0.1
#theta=60: censoring rate=0.3
#theta=30: censoring rate=0.5
#theta=10: censoring rate=0.7
#theta=1: censoring rate=0.9

theta=200
C=runif(n,0,theta)
time=pmin(df$surv_gompertz,C)
status=ifelse(df$surv_gompertz<=C,1,0)
data=data.frame(status,time)

g1<-ggplot(data, aes(time)) + 
  geom_histogram(bins = 20,colour="white",fill = "lightblue") +
  theme_bw()+
  xlim(c(0,110))+
  labs(x = "Censoring Rate 0.1",y="Frequency")


theta=60
C=runif(n,0,theta)
time=pmin(df$surv_gompertz,C)
status=ifelse(df$surv_gompertz<=C,1,0)
data=data.frame(status,time)

g2<-ggplot(data, aes(time)) + 
  geom_histogram(bins = 20,colour="white",fill = "lightblue") +
  theme_bw()+
  xlim(c(0,110))+
  labs(x = "Censoring Rate 0.3",y="Frequency")

theta=30
C=runif(n,0,theta)
time=pmin(df$surv_gompertz,C)
status=ifelse(df$surv_gompertz<=C,1,0)
data=data.frame(status,time)

g3<-ggplot(data, aes(time)) + 
  geom_histogram(bins = 20,colour="white",fill = "lightblue") +
  theme_bw()+
  xlim(c(0,110))+
  labs(x = "Censoring Rate 0.5",y="Frequency")

theta=10
C=runif(n,0,theta)
time=pmin(df$surv_gompertz,C)
status=ifelse(df$surv_gompertz<=C,1,0)
data=data.frame(status,time)

g4<-ggplot(data, aes(time)) + 
  geom_histogram(bins = 20,colour="white",fill = "lightblue") +
  theme_bw()+
  xlim(c(0,110))+
  labs(x = "Censoring Rate 0.7",y="Frequency")


g=grid.arrange(g1,g2,g3,g4, ncol=2)

ggsave("C:/Users/Donglei/Desktop/Simulation_result//histogram_censoring.png", g, height=12, width=12, units='in', dpi=600)

  
