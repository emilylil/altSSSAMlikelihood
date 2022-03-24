#catchfitfunctions.R

#Multivariate lognormal catch-at-age
##################################################################################
fitMVLN<-function(sd,rho)
{
  # varvec = rep(0.2,ages)
  sdvec = rep(sd,ages)
  sddiag = matrix(0,ages,ages)
  diag(sddiag)=sdvec
  
  # rho=0.3
  rho = rho
  cormat = matrix(0,ages,ages)
  for(i in 1:ages)
  {
    for(j in 1:ages)
    {
      cormat[i,j]=rho^abs(i-j)
    }
  }
  varcovmat = sddiag%*%cormat%*%sddiag
  
  #Either one produces Observed Catch on regular scale
  NLLfunction<-function(param,ObsCatage,varcovmat)
  {
    EstCatage<-param
    NLLsum <- 0
    for(i in 1:nrow(ObsCatage))
    {
      NLLsum = NLLsum - dmvnorm(log(ObsCatage[i,]+1),log(EstCatage+1),varcovmat,log=T)
    }
    return(NLLsum)
  }
  #Set initial parameters:
  param = N*0.05
  EstC<-nlm(NLLfunction,param,ObsCatage=ObsCatage,varcovmat=varcovmat,hessian=F,print.level=1,iterlim=10000)
  return(EstC$estimate)
}

#Log normal total catch
##################################################################################
fitLN<- function(LNSD)
{
  #Either one produces Observed Catch on regular scale
  NLLfunction<-function(param,ObsC,LNSD)
  {
    TotC<-param
    return(sum(-dnorm(log(ObsC),log(TotC),sd=LNSD,log=T)))
  }
  #Set initial parameters:
  param = 100000
  EstC<-nlm(NLLfunction,param,ObsC=ObsC,LNSD=LNSD,hessian=T,print.level=1,iterlim=10000)
  return(EstC$estimate)
}

#Multinomial proportions at age
##################################################################################
fitMN<- function(MNESS)
{
  ObsPropCadj = ObsPropC
  for(i in 1:years)
  {
    ObsPropCadj[i,] = (ObsPropC[i,]+0.001)/sum(ObsPropC[i,]+0.001)
  }
  #Either one produces Observed Catch on regular scale
  NLLfunction<-function(param,ObsPropCadj,MNESS)
  {
    probs<-exp(param)
    probs[ages]=exp(mean(param))
    probs <- probs/sum(probs)
    
    MNESS <- MNESS
    NLLsum <- 0
    for(i in 1:nrow(ObsPropCadj))
    {
      NLLsum = NLLsum - dmultinom((ObsPropCadj[i,]*MNESS),prob=probs,log=T)
    }
    return(NLLsum)

  }
  #Set initial parameters:
  # param = PropC[1:11]
  param = rep(0.1,(ages-1))
  
  EstC<-nlm(NLLfunction,param,ObsPropCadj=ObsPropCadj,MNESS=MNESS,hessian=T,print.level=1,iterlim=1000)
  
  estimates = exp(EstC$estimate)
  estimates[ages] = exp(mean(param))
  estimates <- estimates/sum(estimates)
  return(estimates)
}

#Dirichlet multinomial proportions at age
##################################################################################
fitDMN<- function(ISS)
{
  NLLfunction<-function(param,ObsPropC,ISS)
  {
    
    probs<-exp(param[1:ages])
    probs[ages]=exp(mean(param[1:ages]))
    probs <- probs/sum(probs)
    
    theta <- exp(param[ages])/(1+exp(param[ages]))
    
    ISS <- ISS
    NLLsum <- 0
    for(i in 1:nrow(ObsPropC))
    {
      NLLsum = NLLsum - lgamma(ISS*theta)
      NLLsum = NLLsum + lgamma(ISS+ISS*theta)
      
      for(j in 1:ncol(ObsPropC))
      {
        # if(ObsPropC>0)
        # {
          NLLsum = NLLsum - lgamma(ISS*(ObsPropC[i,j]+0.001)+theta*ISS*probs[j])
          NLLsum = NLLsum + lgamma(theta*ISS*probs[j])
        # }
      }
    }
    return(NLLsum)
    
  }
  #Set initial parameters:
  # param = PropC[1:11]
  param = rep(0.1,(ages-1))
  param[ages] = 0
  
  EstC<-nlm(NLLfunction,param,ObsPropC=ObsPropC,ISS=ISS,hessian=T,print.level=1,iterlim=1000)
  
  estimates = exp(EstC$estimate[1:(ages-1)])
  estimates[ages] = exp(mean(EstC$estimate[1:(ages-1)]))
  estimates <- estimates/sum(estimates)
  
  esttheta <- exp(EstC$estimate[ages])/(1+exp(EstC$estimate[ages]))
  
  estimatesandtheta = c(estimates,esttheta)
  return(estimatesandtheta)
  
}

#Logistic normal proportions at age
##################################################################################
fitLogisN<- function(LNSD)
{}