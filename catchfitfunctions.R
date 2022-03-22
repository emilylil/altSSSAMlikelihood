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
  #Either one produces Observed Catch on regular scale
  NLLfunction<-function(param,ObsPropC,MNESS)
  {
    probs<-exp(param+1)
    probs[12]=exp(1)
    probs <- probs/sum(probs)
    
    MNESS <- MNESS
    NLLsum <- 0
    for(i in 1:nrow(ObsPropC))
    {
      NLLsum = NLLsum - dmultinom((ObsPropC[i,]*MNESS),prob=probs,log=T)
    }
    return(NLLsum)

  }
  #Set initial parameters:
  # param = PropC[1:11]
  param = rep((1/12),11)
  EstC<-nlm(NLLfunction,param,ObsPropC=ObsPropC,MNESS=MNESS,hessian=T,print.level=1,iterlim=10000)
  
  estimates = exp(EstC$estimate+1)
  estimates[12] = exp(1)
  estimates <- estimates/sum(estimates)
  return(estimates)
}

#Dirichlet multinomial proportions at age
##################################################################################
fitDMN<- function(LNSD)
{}

#Logistic normal proportions at age
##################################################################################
fitLogisN<- function(LNSD)
{}