# Updates on 2019-07-04:
# Add censoring 


# Updates on 2019-06-19:
# Add noisy Myc-target genes which might not contribute to survival, control the proportion to see difference

# Name: Donglei Yin
# Date: 2019-06-19
# Purpose: Simulation of Myc effect on gene regulation and survival


library(survival)
library(ggplot2)

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
m=10
K=10
lambda=0.01
alpha=0.1
nb_mean=5
nb_size=1
coef_myc=2
theta=100

simu_myc_effect<-function(n,m,K,lambda,alpha,nb_mean,nb_size,coef_myc,theta){
  
  
  # Part 1. Simulate myc's effect on target gene expression
  
  simu_myc<-function(n,m,mean,myc,size){
    df<-NULL
    weight=runif(m,0,2) # weight assigned to each of the K gene, indicating myc's effect size on the gene
    for(k in 1:m){
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
  
  
  myc=rnorm(n,0,1)
  res=simu_myc(n,m,mean=nb_mean,myc=myc,size=nb_size)
  
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
  
  df2$myc_score1=myc_score1
  df2$myc_score2=myc_score2
  df2$myc_score3_up=myc_score3_up
  df2$myc_score3_down=myc_score3_down
  
  # 4. normalized weighted sum
  
  #myc_score4=(as.matrix(df2) %*% weight)/sum(weight^2)
  
  
  # Part 2.Simulate survival time using Cox PH model based on exponential distribution as basedline hazard distribution
  
  # True model covariates: x1,x2,x3,x4,x5,myc-target gene(1,2,...K)
  
  x1=rnorm(n,0,1)
  x2=rnorm(n,0,1)
  x3=runif(n,0,1)
  x4=rbinom(n,1,0.5)
  x5=rbinom(n,1,0.2)
  
  
  df1=data.frame(x1,x2,x3,x4,x5,g1=df2$g1,g2=df2$g2,g3=df2$g3,g4=df2$g4,g5=df2$g5,g6=df2$g6,g7=df2$g7,g8=df2$g8,g9=df2$g9,g10=df2$g10)
  
  coef_base=c(1,1,1,1,1)
  coef=c(coef_base,rep(coef_myc,K)) # coeffient of baseline covariates: x1,x2,x3,x4,x5,myc-target gene(1,2,...K)
  
  # use inverese cdf formula to generate survival time
  
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
  
  
  #hist(surv)
  
  df1$time=surv
  df1$status=1
  
  
  df<-merge(df1,df2)
  
  
  fit_myc <- coxph(Surv(time = df$time, event = df$status) ~ x1+x2+x3+x4+x5+g1+g2+g3+g4+g5+g6+g7+g8+g9+g10, data = df)
  #summary(fit_myc) 
  
  
  fit_baseline <- coxph(Surv(time = df$time, event = df$status) ~ x1+x2+x3+x4+x5, data = df)
  #summary(fit_baseline) 
  
  fit_myc_gene <- coxph(Surv(time = df$time, event = df$status) ~.-time-status-myc_score1-myc_score2-myc_score3_up-myc_score3_down, data = df)
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
  out$Proportion_of_noise<-m/K
  return(out)
}

#theta=100000: censoring rate=0
#theta=200: censoring rate=0.1
#theta=100: censoring rate=0.2
#theta=60: censoring rate=0.3
#theta=30: censoring rate=0.5
#theta=20: censoring rate=0.6
#theta=10: censoring rate=0.7
#theta=1: censoring rate=0.9

para=100

res1=simu_myc_effect(n=1000,m=10,K=10,lambda=0.01,alpha=0.1,nb_mean=5,nb_size=1,coef_myc=1,theta=para)
res2=simu_myc_effect(n=1000,m=100,K=10,lambda=0.01,alpha=0.1,nb_mean=5,nb_size=1,coef_myc=1,theta=para)
res3=simu_myc_effect(n=1000,m=200,K=10,lambda=0.01,alpha=0.1,nb_mean=5,nb_size=1,coef_myc=1,theta=para)
res4=simu_myc_effect(n=1000,m=300,K=10,lambda=0.01,alpha=0.1,nb_mean=5,nb_size=1,coef_myc=1,theta=para)
res5=simu_myc_effect(n=1000,m=400,K=10,lambda=0.01,alpha=0.1,nb_mean=5,nb_size=1,coef_myc=1,theta=para)
res6=simu_myc_effect(n=1000,m=500,K=10,lambda=0.01,alpha=0.1,nb_mean=5,nb_size=1,coef_myc=1,theta=para)


final=rbind(res1,res2,res3,res4,res5,res6,res7)


outdir=paste0('C:/Users/Donglei/Desktop/Simulation_result/noisy_censoring_0.9','.png')
png(outdir, height = 7, width = 9, units = 'in', res = 600)

p <- ggplot(final, aes(x=Proportion_of_noise, y=C_index, group = Model, color=Model))+ 
  geom_errorbar(aes(ymin=C_index-SD, ymax=C_index+SD), width=.1, 
                position=position_dodge(0.05)) +
  geom_line(aes(linetype=Model)) + 
  geom_point(aes(shape=Model))+
  scale_x_continuous(breaks = c(1,5,10,20,30,40,50))+
  labs(x="Nosiy Myc-gene Multiplier", y = "C-index")+
  theme_classic()+
  theme(panel.border = element_rect(colour = "black", fill=NA))
p

dev.off()
