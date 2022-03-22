# Simulate/Estimate catch data with different observation error
# Emily Morgan Liljestrand
# Created: 3/14/22
# Updated: 3/14/22

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
#Dirichlet multinomial proportion-at-age (effective sample size)
ObsPropC = simDMN()
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
#Dirichlet multinomial proportion-at-age (effective sample size)
EstPropC = fitDMN()
#Logistic Normal proportion-at-age (effective sample size)
EstPropC = fitLogisN()

#Reconstruct as necessary:
EstCatage = EstC*EstPropC

##################################### COMPARISON #####################################

mean(EstCatAge-C)/C
hist(EstCatAge-C)

mean(EstPropC-PropC)
hist(EstPropC-PropC)

mean(EstC-TotC)
hist(EstC-TotC)

