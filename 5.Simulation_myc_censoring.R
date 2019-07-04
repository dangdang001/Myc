library(survival)
library(ggplot2)


#n:number of patients/rows
#K:number of myc-target genes
#lambda:parameter of exponential distribution
#nb_mean:parameter of NB distribution
#nb_size:parameter of NB distribution
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

n=10000
K=10
lambda=0.01
alpha=0.1
nb_mean=5
nb_size=1
coef_myc=2
theta=100

simu_myc_effect<-function(n,K,lambda,alpha,nb_mean,nb_size,coef_myc,theta){
  
  # 1.Simulate survival time using Cox PH model based on exponential distribution as basedline hazard distribution
  
  # baseline covariates: x1,x2,x3,x4,x5,myc
  
  x1=rnorm(n,0,1)
  x2=rnorm(n,0,1)
  x3=rnorm(n,0,1)
  x4=rbinom(n,1,0.2)
  x5=rbinom(n,1,0.5)
  myc=rnorm(n,0,1)
  
  df1=data.frame(x1,x2,x3,x4,x5,myc)
  
  coef=c(1,1,1,1,1,coef_myc) # coeffient of baseline covariates: x1,x2,x3,x4,x5,myc
  
  # use inverese cdf formula to generate survival time
  
  # simu_exp<-function(x,lambda){
  #   u=runif(1,0,1)
  #   t=-log(u)/(lambda*exp(coef %*% x))
  #   return(t)
  # }
  simu_gompertz<-function(x,lambda,alpha){
    u=runif(1,0,1)
    t=1/alpha*log(1-(alpha*log(u))/(lambda*exp(coef %*% x)))
    return(t)
  }

  #true survival time
  #surv=apply(df1,1,simu_exp,lambda)
  surv=apply(df1,1,simu_gompertz,lambda,alpha)
  
  #censoring time
  C=runif(n,0,theta)
  
  #hist(surv)

  df1$time=pmin(surv,C)
  df1$status=ifelse(surv<=C,1,0)
  
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
    return(list(df,weight))
  }
  
  
  
  res=simu_myc(n,K,mean=nb_mean,myc=myc,size=nb_size)
  
  df2=res[[1]]
  weight=res[[2]]
  
  # myc-scores:
  
  # 1. simple weighted sum
  myc_score1=as.matrix(df2) %*% weight
  
  # 2. signed weighted sum
  sign=ifelse(weight>1,1,-1)
  myc_score2=as.matrix(df2) %*% (weight*sign)
  
  # 3. down and up seperated sum
  ind=weight>1
  myc_score3_up=as.matrix(df2[ind]) %*% weight[ind]
  myc_score3_down=as.matrix(df2[-ind]) %*% weight[-ind]
  
  # 4. normalized weighted sum
  
  #myc_score4=(as.matrix(df2) %*% weight)/sum(weight^2)
  
  
  df<-cbind(df1,df2,myc_score1,myc_score2,myc_score3_up,myc_score3_down)
  
  
  fit_myc <- coxph(Surv(time = df$time, event = df$status) ~ x1+x2+x3+x4+x5+myc, data = df)
  #summary(fit_myc) 
  
  
  fit_baseline <- coxph(Surv(time = df$time, event = df$status) ~ x1+x2+x3+x4+x5, data = df)
  #summary(fit_baseline) 
  
  fit_myc_gene <- coxph(Surv(time = df$time, event = df$status) ~.-myc-time-status-myc_score1-myc_score2-myc_score3_up-myc_score3_down, data = df)
  #summary(fit_myc_gene) 
  
  fit_myc_score1 <- coxph(Surv(time = df$time, event = df$status) ~x1+x2+x3+x4+x5+myc_score1, data = df)
  #summary(fit_myc_score1) 
  
  fit_myc_score2 <- coxph(Surv(time = df$time, event = df$status) ~x1+x2+x3+x4+x5+myc_score2, data = df)
  #summary(fit_myc_score2) 
  
  fit_myc_score3 <- coxph(Surv(time = df$time, event = df$status) ~x1+x2+x3+x4+x5+myc_score3_up+myc_score3_down, data = df)
  #summary(fit_myc_score3) 
  
  
  out<-data.frame(rbind(summary(fit_myc)$concordance,
                        summary(fit_baseline)$concordance,
                        summary(fit_myc_gene)$concordance,
                        summary(fit_myc_score1)$concordance,
                        summary(fit_myc_score2)$concordance,
                        summary(fit_myc_score3)$concordance))
  
  colnames(out)<-c('C_index','SD')
  out$Model<-c('True Model','Baseline','Myc_gene','Myc_score1','Myc_score2','Myc_score3')
  out$Myc_effect<-coef_myc
  return(out)
}

#theta=100000: censoring rate=0
#theta=200: censoring rate=0.1
#theta=100: censoring rate=0.2
#theta=60: censoring rate=0.3
#theta=30: censoring rate=0.5
#theta=20: censoring rate=0.6
#theta=10: censoring rate=0.7

para=1

res1=simu_myc_effect(n=1000,K=10,lambda=0.01,alpha=0.1,nb_mean=5,nb_size=1,coef_myc=1,theta=para)
res2=simu_myc_effect(n=1000,K=10,lambda=0.01,alpha=0.1,nb_mean=5,nb_size=1,coef_myc=2,theta=para)
res3=simu_myc_effect(n=1000,K=10,lambda=0.01,alpha=0.1,nb_mean=5,nb_size=1,coef_myc=3,theta=para)
res4=simu_myc_effect(n=1000,K=10,lambda=0.01,alpha=0.1,nb_mean=5,nb_size=1,coef_myc=4,theta=para)
res5=simu_myc_effect(n=1000,K=10,lambda=0.01,alpha=0.1,nb_mean=5,nb_size=1,coef_myc=5,theta=para)
res6=simu_myc_effect(n=1000,K=10,lambda=0.01,alpha=0.1,nb_mean=5,nb_size=1,coef_myc=6,theta=para)
res7=simu_myc_effect(n=1000,K=10,lambda=0.01,alpha=0.1,nb_mean=5,nb_size=1,coef_myc=7,theta=para)
res8=simu_myc_effect(n=1000,K=10,lambda=0.01,alpha=0.1,nb_mean=5,nb_size=1,coef_myc=8,theta=para)
res9=simu_myc_effect(n=1000,K=10,lambda=0.01,alpha=0.1,nb_mean=5,nb_size=1,coef_myc=9,theta=para)
res10=simu_myc_effect(n=1000,K=10,lambda=0.01,alpha=0.1,nb_mean=5,nb_size=1,coef_myc=10,theta=para)


final=rbind(res1,res2,res3,res4,res5,res6,res7,res8,res9,res10)

outdir=paste0('C:/Users/Donglei/Desktop/Simulation_result/censoring_0.9','.png')
png(outdir, height = 7, width = 9, units = 'in', res = 600)

p <- ggplot(final, aes(x=Myc_effect, y=C_index, group = Model, color=Model))+ 
  geom_errorbar(aes(ymin=C_index-SD, ymax=C_index+SD), width=.1, 
                position=position_dodge(0.05)) +
  geom_line(aes(linetype=Model)) + 
  geom_point(aes(shape=Model))+
  ylim(c(0.5,1))+
  scale_x_continuous(breaks = c(1,2,3,4,5,6,7,8,9,10))+
  labs(x="Myc-effect multiplier", y = "C-index")+
  theme_bw()
p
dev.off()
