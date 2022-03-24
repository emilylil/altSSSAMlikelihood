#catchsimfunctions.R

#Multivariate lognormal catch-at-age
##################################################################################
simMVLN<- function(sd,rho)
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
  ObsC = rlmvnorm(years,LogC,varcovmat)
  # ObsC = exp(rmvnorm(years,LogC,varcovmat))
  return(ObsC)
}

#Log normal total catch
##################################################################################
# Standard deviation of total catch and ESS of multinomial
simLN<- function(LNSD)
{
  ObsTotC = rlnorm(years,TotLogC,LNSD)
  return(ObsTotC)
}

#Multinomial proportions at age
##################################################################################
# Standard deviation of total catch and ESS of multinomial
simMN<- function(MNESS)
{
  ObsPropC = t(rmultinom(years,MNESS,prob=PropC)/MNESS)
  return(ObsPropC)
}

#Dirichlet multinomial proportions at age
##################################################################################
simDMN<- function(ISS,theta)
{
  ISS = ISS
  theta = theta
  
  beta = ISS*theta
  alpha = beta*PropC
  
  p = rdirichlet(n=years, alpha=alpha)
  x=matrix(nrow=years,ncol=ages)
  for(i in 1:years)
  {
    x[i,] = t(rmultinom(1,ISS,p[i,])/ISS)
  }
  return(x)
}

#Logistic normal proportions at age
##################################################################################
simLogisN<- function()
{
  ObsPropC = 0
  return(ObsPropC)
}