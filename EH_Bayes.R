EH_Bayes = function(data,nsamp, conf=c(0.025,0.975), k=100){
require(LearnBayes)
require(lhs)


##log transform data
IH_data<-log(data)

 


 S = sum((IH_data - mean(IH_data))^2) 
 n = length(IH_data)
 r_s <- randomLHS(nsamp/k,k)
 sigma2 = as.vector(S/qchisq(r_s, n - 1))

 r_m <- randomLHS(nsamp/k, k)
 mu = as.vector(qnorm(r_m, mean = mean(IH_data), sd = sqrt(sigma2)/sqrt(n)))

y_bar<-exp(mu+0.5*sigma2)

LCL95<-quantile(y_bar,conf[[1]])
UCL95<-quantile(y_bar,conf[[2]])
CI95<-c(LCL95,UCL95)
names(CI95)<-c("LCL","UCL")

return(CI95)
}




EH_Cen_Bayes = function(obs, cen, nsamp=2000, conf=c(0.025,0.975), k=100, restriction = 5){
  require(LearnBayes)
  require(lhs)
  require(NADA)
  
  rosData <- as.data.frame(ros(obs,cen))
  data <- rosData$modeled
  
  ##log transform data
  IH_data<-log(data)
  restriction<-log(restriction)
  
  
  
  
  S = sum((IH_data - mean(IH_data))^2) 
  n = length(IH_data)
  r_s <- randomLHS(nsamp/k,k)
  sigma2 = as.vector(S/qchisq(r_s, n - 1))
  
  r_m <- randomLHS(nsamp/k, k)
  mu = as.vector(qnorm(r_m, mean = mean(IH_data), sd = sqrt(sigma2)/sqrt(n)))
  
  
  eval_set<-cbind(mu,sigma2)
  eval_set<-data.frame(eval_set)
  
  eval_set<-eval_set[eval_set$sigma2<restriction,]
  mu1<-eval_set$mu
  sigma21<-eval_set$sigma2
  
  
  y_bar<-exp(mu1+0.5*sigma21)
  
  LCL95<-quantile(y_bar,conf[[1]])
  UCL95<-quantile(y_bar,conf[[2]])
  CI95<-c(LCL95,UCL95)
  names(CI95)<-c("LCL","UCL")
  
  return(CI95)
}











EH_Bayes_restricted=function(data,nsamp,restriction){
  library(LearnBayes)
  
  
  ##log transform data
  IH_data<-log(data)
  restriction<-log(restriction)
  
  
  
  S = sum((IH_data - mean(IH_data))^2) 
  n = length(IH_data)
  sigma2 = S/rchisq(nsamp, n - 1)
  
  
  mu = rnorm(nsamp, mean = mean(IH_data), sd = sqrt(sigma2)/sqrt(n))
  
  eval_set<-cbind(mu,sigma2)
  eval_set<-data.frame(eval_set)
  
  eval_set<-eval_set[eval_set$sigma2<restriction,]
  mu1<-eval_set$mu
  sigma21<-eval_set$sigma2
  
  
  y_bar<-exp(mu1+0.5*sigma21)
  
  LCL95<-quantile(y_bar,0.025)
  UCL95<-quantile(y_bar,0.975)
  CI95<-c(LCL95,UCL95)
  names(CI95)<-c("LCL","UCL")
  
  return(CI95)
}


# Examples

EH_Bayes(c(0.1,0.2,0.05),1000)

# restrict GSD <5

EH_Bayes_restricted(c(0.1,0.2,0.05),1000,5)

