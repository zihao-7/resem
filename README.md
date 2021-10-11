# Rerandomized Survey Experiment using Mahalanobis distances  
 
 This package provides functions for conducting ReSEM(Rerandomized Survey Experiment using Mahalanobis distances). You can use functions in this package to generate sampling and treatment assignment indicators. You can also compute difference-in-means estimators and regression adjusted difference-in-means estimators for the average treatment effect, as well as the corresponding confidence intervals.

## Installation
library(devtools)

install_github("zihao-7/resem")

## Library the package
library(ReSEM)

## Explanation for the main functions
?assign.ReSEM.sampling

?assign.ReSEM.treatment

?dif.means

?dif.reg

?CI.ReSEM

## A simple example

### total size, sample size, and treated size 
N = 1000

n = 100

n1 = 50

### asymptotic acceptance probabilities
p.S = 0.01

p.T = 0.01

### generate covariate
W = matrix(NA, nrow = N, ncol = 2)

X = matrix(NA, nrow = N, ncol = 4)

C = matrix(NA, nrow = N, ncol = 6)

C[,1] = as.numeric( runif(N) <= 0.5 )

C[,2] = rnorm(N)

C[,3] = as.numeric( runif(N) <= 0.5 )

C[,4] = rnorm(N)

C[,5] = as.numeric( runif(N) <= 0.5 )

C[,6] = rnorm(N)

W = C[,1:2]

X = C[,1:4]

E = C[,1:2]

### generate potential outcomes
Y0 = -0.5 * rowSums(C) + 0.1 * rnorm(N)

tau = 0.6 * rowSums(C)

Y1 = Y0 + tau

### generate constrained Gaussian variables to compute quantile range
num.simu = 10^5

L.J.aS = generate.constrained.Gaussian(num.simu, ncol(W), qchisq(p.S, df = ncol(W)))

L.K.aT = generate.constrained.Gaussian(num.simu, ncol(X), qchisq(p.T, df = ncol(X)))

epsilon = rnorm(num.simu, mean = 0, sd = 1)

### get 100 treatment assignment from ReSEM
iter.max = 100

ind.ReSEM.all = assign.ReSEM.sampling(num.assign = iter.max, N, n, W, prob.S=p.T)$ind.assign

for (k in 1:nrow(ind.ReSEM.all)) {

  sampled = (ind.ReSEM.all[k,] == 0)
  
  treatment.ind = assign.ReSEM.treatment(n, n1, X = X[sampled,], prob.T=p.S)$ind.assign
  
  ind.ReSEM.all[k, sampled] = treatment.ind
  
}

### difference-in-means estimator, and regression adjusted difference-in-means estimators
tau.diff.ReSEM = rep(NA, iter.max)

tau.reg.ReSEM = rep(NA, iter.max)

for(iter in 1:iter.max){

  tau.diff.ReSEM[iter] = dif.means( ind.ReSEM.all[iter,], Y = obs.outcome(ind.ReSEM.all[iter,], Y1, Y0) )
  
  tau.reg.ReSEM[iter] = dif.reg( ind.ReSEM.all[iter,], Y = obs.outcome(ind.ReSEM.all[iter,], Y1, Y0), C, E )
  
}

### CIs for the two estimators
alpha = 0.05

CI.diff.all = matrix(NA, nrow=iter.max, ncol=2)

CI.reg.all = matrix(NA, nrow=iter.max, ncol=2)

for (k in 1:iter.max) {

  assignment = ind.ReSEM.all[k,]
  
  Y = obs.outcome(assignment, Y1, Y0)
  
  CI.diff.all[k,] = CI.ReSEM(assignment, Y, X, W, C, E, alpha, design ="ST") 
  
  CI.reg.all[k,] = CI.ReSEM(assignment, Y, X, W, C, E, alpha, design ="ST.adjusted") 
  
}

coverage.CI.diff = mean( mean(tau) > CI.diff.all[,1] & mean(tau) < CI.diff.all[,2] )

coverage.CI.reg = mean( mean(tau) > CI.reg.all[,1] & mean(tau) < CI.reg.all[,2] )


