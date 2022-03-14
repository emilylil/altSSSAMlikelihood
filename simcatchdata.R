# Simulate/Estimate catch data with different observation error
# Emily Morgan Liljestrand
# Created: 3/14/22
# Updated: 3/14/22

# Simulate 32 years of catch data with 12 ages
# Option 1: Multivariate lognormal catch-at-age
# Option 2: Log normal catch, multinomial proportions
# Option 3: Log normal catch, dirichlet multinomial proportions
# Option 4: Log normal catch, logistic normal proportions

# Variables
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

#Calculate observed catch
ObsC=matrix(0,nrow=years,ncol=ages)

# Option 1: Multivariate lognormal catch-at-age
var = rep(0,ages)
rho = 0.3

for(i in 1:years)
{
  ObsC[i,] = 0
}

# Option 2: Log normal catch, multinomial proportions


# Option 3: Log normal catch, dirichlet multinomial proportions


# Option 4: Log normal catch, logistic normal proportions
