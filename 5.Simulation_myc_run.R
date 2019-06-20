library(survival)


n=1000 # number of patients/rows
K=10 # number of myc-target genes
lambda=1 # parameter of exponential distribution
nb_mean=5 # parameter of NB distribution
nb_size=1 # parameter of NB distribution
coef=c(1,1,1,1,1,5) # coeffient of baseline covariates: x1,x2,x3,x4,x5,myc


# 1.Simulate survival time using Cox PH model based on exponential distribution as basedline hazard distribution

# baseline covariates: x1,x2,x3,x4,x5,myc

x1=rnorm(n,0,1)
x2=rnorm(n,0,1)
x3=runif(n,0,1)
x4=rbinom(n,1,0.5)
x5=rbinom(n,1,0.2)
myc=rnorm(n,0,1)

df1=data.frame(x1,x2,x3,x4,x5,myc)

# use inverese cdf formula to generate survival time

simu_exp<-function(x,lambda){
  u=runif(1,0,1)
  t=-log(u)/(lambda*exp(coef %*% x))
  return(t)
}

surv_exp=apply(df1,1,simu_exp,lambda=lambda)

#hist(surv_exp)

df1$surv=surv_exp
df1$status=1



# 2. Simulate myc's effect on target gene expression

simu_myc<-function(n,K,mean,myc,size){
  df<-NULL
  weight=runif(K,0,2) # weight assigned to each of the K gene, indicating myc's effect size on the gene
  for(k in 1:K){
    effect=exp(myc)*weight[k]
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

df2=simu_myc(n,K,mean=nb_mean,myc=myc,size=nb_size)

df<-cbind(df1,df2)




fit_myc <- coxph(Surv(time = df$surv, event = df$status) ~ x1+x2+x3+x4+x5+myc, data = df)
summary(fit_myc) 


fit_baseline <- coxph(Surv(time = df$surv, event = df$status) ~ x1+x2+x3+x4+x5, data = df)
summary(fit_baseline) 


fit_myc_gene <- coxph(Surv(time = df$surv, event = df$status) ~.-myc-surv-status, data = df)
summary(fit_myc_gene) 




library(rms)

f <- cph(Surv(time = df$surv, event = df$status) ~ x1+x2+x3+x4+x5+myc, x=TRUE, y=TRUE,surv = TRUE)

f$stats['Dxy']/2+0.5

validate(f, B=100)

cal <- calibrate(f, u=1.5, cmethod='KM', cuts=seq(0,1,0.1), B=20)  # usually B=200 or 300
plot(cal, add=pa)