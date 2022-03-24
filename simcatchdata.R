# Simulate/Estimate catch data with different observation error
# Emily Morgan Liljestrand
# Created: 3/14/22
# Updated: 3/24/22

# Simulate 32 years of "real" catch data with 12 ages
# Simulate observed catch data using several obs error structure options
# Fit observed catch data assuming several obs error structure options
# Compare the fit to the real
rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(MASS)
library(dmutate)
library(emdbook)

##################################### REAL DATA #####################################
# Variables for real data
years = 32
ages = 12
R = 1500000 # Recruitment
M = 0.2 # Natural Mortality
Sel = c(0.1, 0.4, 1, 1, 0.9, 0.75, 0.7, 0.6, 0.4, 0.2, 0.1, 0.05) # Selectivity schedule (domed)
FFull = 0.5 # Fishing mortality at full vulnerability
Z = M + Sel*FFull # Total Mortality

N = rep(0,ages) # Abundance
C = rep(0,ages) # Catch

#Calculate true abundance and catch
for(i in 1:ages)
{
  N[i] = R*((exp(-Z[i])))^(i-1)
  C[i] = N[i]*((FFull*Sel[i])/Z[i])*(1-exp(-Z[i]))
}

# The true total catch, on regular or log scale, and the true proportions
LogC=log(C)
TotC = sum(C)
TotLogC =  log(TotC)
PropC = C/TotC

##################################### SIMULATED DATA #####################################

# Load all the functions to simulate catch data
source("catchsimfunctions.R")
#####################################
#Multivariate lognormal catch-at-age (sd for log-normal, rho)
ObsCatage = simMVLN(0.06,0.3)

#Deconstruct as necessary:
ObsC=rowSums(ObsCatage)
ObsPropC=ObsCatage/ObsC

#####################################
#Log normal catch (sd for total catch)
ObsC = simLN(0.06676)

#Multinomial proportion-at-age (effective sample size)
ObsPropC = simMN(500)
#Dirichlet multinomial proportion-at-age (effective sample size, "theta")
ObsPropC = simDMN(300,0.5)
#Logistic Normal proportion-at-age (effective sample size)
ObsPropC = simLogisN()

#Reconstruct as necessary:
ObsCatage = ObsC*ObsPropC


##################################### ESTIMATED DATA #####################################

# Load all the functions to fit catch data
source("catchfitfunctions.R")

#Multivariate lognormal catch-at-age (sd for log-normal, rho)
EstCatAge = fitMVLN(0.06,0.5)
#Deconstruct as necessary:
EstC=rowSums(EstCatage)
EstPropC=EstCatage/EstC

#####################################
#Log normal catch (sd for total catch)
EstC = fitLN(0.06676)

#Multinomial proportion-at-age (effective sample size)
EstPropC = fitMN(500)

#Dirichlet multinomial proportion-at-age (input sample size)
EstPropC = fitDMN(500)
EstTheta = EstPropC[13]
EstPropC = EstPropC[1:12]
EstESS = 1/(1+EstTheta)+500*(EstTheta/(1+EstTheta))

#Logistic Normal proportion-at-age (effective sample size)
EstPropC = fitLogisN()

#Reconstruct as necessary:
EstCatage = EstC*EstPropC

##################################### COMPARISON #####################################
source("catchsimfunctions.R")
source("catchfitfunctions.R")

nsims = 10
AbsoluteErrorTot=rep(0,nsims)
PercentErrorTot=rep(0,nsims)
AbsoluteErrorProp=rep(0,nsims)

#1. Simulate MVLN, fit MVLN
for(i in 1:nsims)
{
  #Multivariate lognormal catch-at-age (sd for log-normal, rho)
  ObsCatage = simMVLN(0.06,0.3)
  EstCatAge = fitMVLN(0.06,0.5)
  
  AbsoluteErrorTot[i] = sum(EstCatAge)-sum(C)
  PercentErrorTot[i] = (sum(EstCatAge)-sum(C))/sum(C)
}
AbsoluteErrorTot1 = AbsoluteErrorTot
PercentErrorTot1 = PercentErrorTot
hist(AbsoluteErrorTot1)
hist(PercentErrorTot1)

#2. Simulate MVLN, fit LN and MN
for(i in 1:nsims)
{
  #Multivariate lognormal catch-at-age (sd for log-normal, rho)
  ObsCatage = simMVLN(0.06,0.3)
  #Deconstruct:
  ObsC=rowSums(ObsCatage)
  ObsPropC=ObsCatage/ObsC
  
  EstC = fitLN(0.06)
  EstPropC = fitMN(300)
  
  AbsoluteErrorTot[i] = EstC-TotC
  PercentErrorTot[i] = (EstC-TotC)/TotC
  
  AbsoluteErrorProp[i] = EstPropC[2]-PropC[2]
}
AbsoluteErrorTot2 = AbsoluteErrorTot
PercentErrorTot2 = PercentErrorTot
AbsoluteErrorProp2 =AbsoluteErrorProp
hist(AbsoluteErrorTot2)
hist(PercentErrorTot2)
hist(AbsoluteErrorProp2)

#3. Simulate MVLN, fit LN and DMN
EstThetavec=rep(0,nsims)
for(i in 1:nsims)
{
  #Multivariate lognormal catch-at-age (sd for log-normal, rho)
  ObsCatage = simMVLN(0.06,0.3)
  
  #Deconstruct:
  ObsC=rowSums(ObsCatage)
  ObsPropC=ObsCatage/ObsC
  
  EstC = fitLN(0.06)
  
  EstPropC = fitDMN(500)
  EstTheta = EstPropC[13]
  EstPropC = EstPropC[1:12]
  EstESS = 1/(1+EstTheta)+500*(EstTheta/(1+EstTheta))
  EstThetavec[i] = EstTheta
  
  AbsoluteErrorTot[i] = EstC-TotC
  PercentErrorTot[i] = (EstC-TotC)/TotC
  
  AbsoluteErrorProp[i] = EstPropC[2]-PropC[2]
}
AbsoluteErrorTot3 = AbsoluteErrorTot
PercentErrorTot3 = PercentErrorTot
AbsoluteErrorProp3 =AbsoluteErrorProp
hist(AbsoluteErrorTot3)
hist(PercentErrorTot3)
hist(AbsoluteErrorProp3)

#4. Simulate LN and MN, fit MVLN
for(i in 1:nsims)
{
  ObsC = simLN(0.06676)
  ObsPropC = simMN(500)
  
  # Reconstruct:
  ObsCatage = ObsC*ObsPropC
  
  EstCatAge = fitMVLN(0.06,0.5)
  AbsoluteErrorTot[i] = sum(EstCatAge)-sum(C)
  PercentErrorTot[i] = (sum(EstCatAge)-sum(C))/sum(C)
}
AbsoluteErrorTot4 = AbsoluteErrorTot
PercentErrorTot4 = PercentErrorTot
hist(AbsoluteErrorTot4)
hist(PercentErrorTot4)

#5. Simulate LN and MN, fit LN and MN
for(i in 1:nsims)
{
  ObsC = simLN(0.06676)
  ObsPropC = simMN(500)
  
  EstC = fitLN(0.06676)
  EstPropC = fitMN(500)
  
  AbsoluteErrorTot[i] = EstC-TotC
  PercentErrorTot[i] = (EstC-TotC)/TotC
  
  AbsoluteErrorProp[i] = EstPropC[2]-PropC[2]
}
AbsoluteErrorTot5 = AbsoluteErrorTot
PercentErrorTot5 = PercentErrorTot
AbsoluteErrorProp5 =AbsoluteErrorProp
hist(AbsoluteErrorTot5)
hist(PercentErrorTot5)
hist(AbsoluteErrorProp5)

#6. Simulate LN and MN, fit LN and DMN
for(i in 1:nsims)
{
  ObsC = simLN(0.06676)
  ObsPropC = simMN(300)
  
  EstC = fitLN(0.06676)
  EstPropC = fitDMN(500)
  EstTheta = EstPropC[13]
  EstPropC = EstPropC[1:12]
  EstESS = 1/(1+EstTheta)+500*(EstTheta/(1+EstTheta))
  EstThetavec[i] = EstTheta
  
  AbsoluteErrorTot[i] = EstC-TotC
  PercentErrorTot[i] = (EstC-TotC)/TotC
  
  AbsoluteErrorProp[i] = EstPropC[2]-PropC[2]
}
AbsoluteErrorTot6 = AbsoluteErrorTot
PercentErrorTot6 = PercentErrorTot
AbsoluteErrorProp6 =AbsoluteErrorProp
hist(AbsoluteErrorTot6)
hist(PercentErrorTot6)
hist(AbsoluteErrorProp6)

#7. Simulate LN and DMN, fit MVLN
for(i in 1:nsims)
{
  ObsC = simLN(0.06676)
  ObsPropC = simDMN(500)
  
  # Reconstruct:
  ObsCatage = ObsC*ObsPropC
  
  EstCatAge = fitMVLN(0.06,0.5)
  AbsoluteErrorTot[i] = sum(EstCatAge)-sum(C)
  PercentErrorTot[i] = (sum(EstCatAge)-sum(C))/sum(C)
}
AbsoluteErrorTot7 = AbsoluteErrorTot
PercentErrorTot7 = PercentErrorTot
hist(AbsoluteErrorTot7)
hist(PercentErrorTot7)

#8. Simulate LN and DMN, fit LN and MN
for(i in 1:nsims)
{
  ObsC = simLN(0.06676)
  ObsPropC = simDMN(500)
  
  EstC = fitLN(0.06676)
  EstPropC = fitMN(300)
  
  AbsoluteErrorTot[i] = EstC-TotC
  PercentErrorTot[i] = (EstC-TotC)/TotC
  
  AbsoluteErrorProp[i] = EstPropC[2]-PropC[2]
}
AbsoluteErrorTot8 = AbsoluteErrorTot
PercentErrorTot8 = PercentErrorTot
AbsoluteErrorProp8 =AbsoluteErrorProp
hist(AbsoluteErrorTot8)
hist(PercentErrorTot8)
hist(AbsoluteErrorProp8)

#9. Simulate LN and DMN, fit LN and DMN
for(i in 1:nsims)
{
  ObsC = simLN(0.06676)
  ObsPropC = simDMN(300)
  
  EstC = fitLN(0.06676)
  EstPropC = fitDMN(500)
  EstTheta = EstPropC[13]
  EstPropC = EstPropC[1:12]
  EstESS = 1/(1+EstTheta)+500*(EstTheta/(1+EstTheta))
  EstThetavec[i] = EstTheta
  
  AbsoluteErrorTot[i] = EstC-TotC
  PercentErrorTot[i] = (EstC-TotC)/TotC
  
  AbsoluteErrorProp[i] = EstPropC[2]-PropC[2]
}
AbsoluteErrorTot9 = AbsoluteErrorTot
PercentErrorTot9 = PercentErrorTot
AbsoluteErrorProp9 =AbsoluteErrorProp
hist(AbsoluteErrorTot9)
hist(PercentErrorTot9)
hist(AbsoluteErrorProp9)