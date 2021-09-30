#' Rejective Sampling
#'
#' Generates sampling indicators under rejective sampling using the Mahalanobis distance.
#'
#' @param num.assign Number of assignments to be generated.
#' @param N Size of the finite population.
#' @param n Number of the sampled units.
#' @param W An \code{N} by J matrix which contains J dimensional covariate vector for each unit, available at the sampling stage. 
#' @param prob.S Asymptotic acceptance probability for the sampling stage.
#'
#' @return A list with the elements
#' \item{ind.assign}{A \code{num.assign} by \code{N} sampling matrix. Each row of the sampling matrix is an \code{N}-dimensional vector that represents one sampling result, where we use 0 for sampled units, and -1 for units that will not enter the experiment. }
#' \item{M.S.all}{An vector of length \code{num.assign}, which stores the Mahalanobis distances corresponding to each sampling. }
#' @export
assign.ReSEM.sampling <- function(num.assign, N, n, W, prob.S) {
  ind.assign = matrix(-1, nrow = num.assign, ncol = N)
  M.S.all = rep(NA, num.assign)
  
  a.S = qchisq(prob.S, df = ncol(W))
  S2.W.inv = solve(var(W))
  W.bar = colMeans(W)
  
  k = 1
  while(k <= num.assign){
    selected = sample(c(1:N), n)
    tau.hat.W = colMeans( W[selected, ] ) - W.bar
    M.S = as.numeric( 1/(1/n - 1/N) * t(tau.hat.W) %*% S2.W.inv %*% tau.hat.W )
    if(M.S <= a.S){
      ind.assign[k, selected] = 0
      M.S.all[k] = M.S
      k = k+1
    }
  }
  return(list(ind.assign=ind.assign, M.S.all=M.S.all))
}

#' Rerandomized treatment assignment
#'
#' Generates treatment-control assignment under rerandomization using the Mahalanobis distance.
#'
#' @param n Number of units to assign.
#' @param n1  Number of the treated units.
#' @param X An \code{n} by K matrix which contains K dimensional covariate vector for each unit, available at the treatment assignment stage.
#' @param prob.T Asymptotic acceptance probability for the treatment assignment stage.
#'
#' @return A list with the elements
#' \item{ind.assign}{An assignment vector of length \code{n}, where we use 1 for treated units, and 0 for control units. }
#' \item{M.T}{The corresponding Mahalanobis distance. }
#' @export
assign.ReSEM.treatment <- function(n, n1, X, prob.T) {
  ind.assign = rep(0, n)
  a.T = qchisq(prob.T, df = ncol(X))
  S2.X.inv = solve(var(X))
  repeat {
    treated = sample(n, n1)
    tau.hat.X = colMeans(X[treated, ]) - colMeans(X[-treated, ])
    M.T = as.numeric( 1/(n/(n1*(n-n1))) * t(tau.hat.X) %*% S2.X.inv %*% tau.hat.X )
    if (M.T <= a.T) {
      break
    }
  }
  ind.assign[treated] = 1
  return(list(ind.assign=ind.assign, M.T=M.T))
}

#' Observed outcomes
#' 
#' Computes observed outcomes using potential outcomes and treatment assignment.
#'
#' @param assignment Treatment assignment vector, which contains -1, 0, 1 values. We use 1 for treated units, 0 for control units, and -1 for units that will not enter the experiment.  
#' @param Y1 Potential outcomes under treatment.
#' @param Y0 Potential outcomes under control.
#' 
#' @return A vector of the same length as \code{assignment}, \code{Y1} and \code{Y0}, which contains observed outcomes for all units. If a unit will not enter the experiment, its observed outcome will be NA.
#' @export
obs.outcome <- function(assignment, Y1, Y0){
  N = length(Y1)
  Y = rep(NA, N)
  Y[assignment == 1] = Y1[assignment == 1]
  Y[assignment == 0] = Y0[assignment == 0]
  return(Y)
}

#' Difference in means
#'
#' Calculates the difference-in-means estimator for the average treatment effect.
#' 
#' @param assignment Treatment assignment vector, which contains -1, 0, 1 values. We use 1 for treated units, 0 for control units, and -1 for units that will not enter the experiment. 
#' @param Y Observed outcomes for all units. Will be NA for units that will not enter the experiment.
#' 
#' @return Difference-in-means estimator for the average treatment effect.
#' @export 
dif.means <- function(assignment, Y){
  return( mean(Y[assignment==1]) - mean( Y[assignment==0] ) )
}

#' Regression-adjusted difference-in-means
#'
#' Calculates the regression-adjusted difference-in-means estimator for the average treatment effect.
#'
#' @param assignment Treatment assignment vector, which contains -1, 0, 1 values. We use 1 for treated units, 0 for control units, and -1 for units that will not enter the experiment. 
#' @param Y Observed outcomes for all units. Will be NA for units that will not enter the experiment.
#' @param C Available covariate vector for sampled units, at the analysis stage.
#' @param E Available covariate vector for all units, at the analysis stage.
#
#' @return Regression-adjusted difference-in-means estimator for the average treatment effect, with estimated adjustment coefficients.
#' @export 
dif.reg <- function(assignment, Y, C, E){
  delta.E = colMeans( E[assignment >= 0, ] ) - colMeans(E)
  tau.C = colMeans( C[assignment == 1, ] ) - colMeans( C[assignment == 0, ] )
  beta.1 = as.numeric( lm(Y[assignment==1] ~ C[assignment==1, ])$coefficients )[-1]
  beta.0 = as.numeric( lm(Y[assignment==0] ~ C[assignment==0, ])$coefficients )[-1]
  gamma.1 = as.numeric( lm(Y[assignment==1] ~ E[assignment==1, ])$coefficients )[-1]
  gamma.0 = as.numeric( lm(Y[assignment==0] ~ E[assignment==0, ])$coefficients )[-1]
  r1 = sum(assignment==1)/sum(assignment>=0)
  r0 = 1 - r1
  beta = r0 * beta.1 + r1 * beta.0
  gamma = gamma.1 -  gamma.0
  tau.hat = mean(Y[assignment==1]) - mean( Y[assignment==0] ) - sum(beta * tau.C) - sum(gamma * delta.E)
  return(tau.hat)
}

#' Confidence interval under ReSEM
#'
#' Computes confidence interval for the average treatment effect under ReSEM.
#'
#' @param assignment A vector of length N, which contains -1, 0, 1 values. We use 1 for treated units, 0 for control units, and -1 for units that will not enter the experiment. 
#' @param Y A vector of length N, which contains observed outcomes for all units. Will be NA for units that will not enter the experiment.
#' @param X An n by K matrix which contains K dimensional covariate vector for sampled units, available at the treatment assignment stage.
#' @param W An N by J matrix which contains J dimensional covariate vector for each unit, available at the sampling stage.
#' @param C An n by K' matrix which contains K' dimensional covariate vector for sampled units, available at the analysis stage. Assume X is a subset of C.
#' @param E An N by J' matrix which contains J' dimensional covariate vector for each unit, available at the analysis stage. Assume W is a subset of E.
#' @param alpha Significance level. The computed CI will have 1 - \code{alpha} confidence level .
#' @param design Rerandomization design used to compute the CI. Can be "ST", "S", "T", "CRSE", "ST.adjusted".
#'
#' @return 1 - \code{alpha} confidence interval for the average treatment effect.
#' @export 
CI.ReSEM <- function(assignment, Y, X, W, C, E, alpha=0.05, design){ 

  r1 = sum(assignment==1)/sum(assignment>=0)
  r0 = 1 - r1
  f = sum(assignment>=0)/length(assignment)
  
  s2.1 = var(Y[assignment==1])
  s2.0 = var(Y[assignment==0])
  s.1X = cov(Y[assignment==1], X[assignment==1,])
  s.0X = cov(Y[assignment==0], X[assignment==0,])
  s.1W = cov(Y[assignment==1], W[assignment==1,])
  s.0W = cov(Y[assignment==0], W[assignment==0,])
  s.1C = cov(Y[assignment==1], C[assignment==1,])
  s.0C = cov(Y[assignment==0], C[assignment==0,])
  s.1E = cov(Y[assignment==1], E[assignment==1,])
  s.0E = cov(Y[assignment==0], E[assignment==0,])
  
  s2.1midX = var( lm(Y[assignment==1] ~ X[assignment==1,])$fitted.value )
  s2.0midX = var( lm(Y[assignment==0] ~ X[assignment==0,])$fitted.value )
  s2.1midC = var( lm(Y[assignment==1] ~ C[assignment==1,])$fitted.value )
  s2.0midC = var( lm(Y[assignment==0] ~ C[assignment==0,])$fitted.value )

  s2.X.1 = var(X[assignment == 1, ])
  s2.X.0 = var(X[assignment == 0, ])
  s2.X.1.neg.half.power = expm::expm(-0.5*expm::logm(s2.X.1))
  s2.X.0.neg.half.power = expm::expm(-0.5*expm::logm(s2.X.0))
  s2.taumidX = t( s2.X.1.neg.half.power %*% t(s.1X) - s2.X.0.neg.half.power %*% t(s.0X) ) %*% 
    ( s2.X.1.neg.half.power %*% t(s.1X) - s2.X.0.neg.half.power %*% t(s.0X) )
  
  s2.W.1 = var(W[assignment == 1, ])
  s2.W.0 = var(W[assignment == 0, ])
  s2.W.1.neg.half.power = expm::expm(-0.5*expm::logm(s2.W.1))
  s2.W.0.neg.half.power = expm::expm(-0.5*expm::logm(s2.W.0))
  s2.taumidW = t( s2.W.1.neg.half.power %*% t(s.1W) - s2.W.0.neg.half.power %*% t(s.0W) ) %*% 
    ( s2.W.1.neg.half.power %*% t(s.1W) - s2.W.0.neg.half.power %*% t(s.0W) )
  
  s2.C.1 = var(C[assignment == 1, ])
  s2.C.0 = var(C[assignment == 0, ])
  s2.C.1.neg.half.power = expm::expm(-0.5*expm::logm(s2.C.1))
  s2.C.0.neg.half.power = expm::expm(-0.5*expm::logm(s2.C.0))
  s2.taumidC = t( s2.C.1.neg.half.power %*% t(s.1C) - s2.C.0.neg.half.power %*% t(s.0C) ) %*% 
    ( s2.C.1.neg.half.power %*% t(s.1C) - s2.C.0.neg.half.power %*% t(s.0C) )
  
  s2.E.1 = var(E[assignment == 1, ])
  s2.E.0 = var(E[assignment == 0, ])
  s2.E.1.neg.half.power = expm::expm(-0.5*expm::logm(s2.E.1))
  s2.E.0.neg.half.power = expm::expm(-0.5*expm::logm(s2.E.0))
  s2.taumidE = t( s2.E.1.neg.half.power %*% t(s.1E) - s2.E.0.neg.half.power %*% t(s.0E) ) %*% 
    ( s2.E.1.neg.half.power %*% t(s.1E) - s2.E.0.neg.half.power %*% t(s.0E) )
  
  V.tt = as.numeric( s2.1/r1 + s2.0/r0 - f * s2.taumidC )
  R2.S = as.numeric( (1-f) / V.tt * s2.taumidW )
  R2.T = 1/V.tt * as.numeric( s2.1midX/r1 + s2.0midX/r0 - s2.taumidX )
  R2.E = as.numeric( (1-f) / V.tt * s2.taumidE )
  R2.C = 1/V.tt * as.numeric( s2.1midC/r1 + s2.0midC/r0 - s2.taumidC )
  
  if (design == "CRSE") {
    QR = c(qnorm(alpha/2), qnorm(1-alpha/2))
    CI = dif.means(assignment, Y) + sqrt(V.tt) * unname(QR) / sqrt(n)
  }
  if (design == "S") {
    QR = quantile( ifelse(R2.S<=1, sqrt(1-R2.S), 0) * epsilon + sqrt(R2.S) * L.J.aS,
                   probs = c(alpha/2, 1-alpha/2) ) 
    CI = dif.means(assignment, Y) + sqrt(V.tt) * unname(QR) / sqrt(n)
  }
  if (design == "T") {
    QR = quantile( ifelse(R2.T<=1, sqrt(1-R2.T), 0) * epsilon + sqrt(R2.T) * L.K.aT,
                   probs = c(alpha/2, 1-alpha/2) )
    CI = dif.means(assignment, Y) + sqrt(V.tt) * unname(QR) / sqrt(n)
  }
  if (design == "ST") {
    QR = quantile( ifelse(R2.S+R2.T<=1, sqrt(1-R2.S-R2.T), 0) * epsilon 
                   + sqrt(R2.S) * L.J.aS + sqrt(R2.T) * L.K.aT, 
                   probs = c(alpha/2, 1-alpha/2) ) 
    CI = dif.means(assignment, Y) + sqrt(V.tt) * unname(QR) / sqrt(n)
  }
  if (design == "ST.adjusted") {
    QR = ifelse(R2.E+R2.C<=1, sqrt(1-R2.E-R2.C), 0) * c(qnorm(alpha/2), qnorm(1-alpha/2))
    CI = dif.reg(assignment, Y, C, E) + sqrt(V.tt) * unname(QR) / sqrt(n)
  }
  return(CI)
}

#' Generate constrained Gaussian
#' 
#' Generates constrained Gaussian random variables.
#'
#' @param num Number of constrained Gaussian random variables to generate.
#' @param K Dimension of the Gaussian vector used to generate constrained Gaussian random variable.
#' @param a Threshold used to generate constrained Gaussian random variable.
#
#' @return A vector of length \code{num}, which contains iid constrained Gaussian random variables.
#' @export 
generate.constrained.Gaussian <- function(num, K, a){
  if (! "Runuran" %in% installed.packages()) { 
    install.packages("Runuran")
  }
  library(Runuran)
  
  chi_aK = sqrt(urchisq(num, df=K, lb=0, ub=a))
  S = 2 * rbinom(num, 1, prob=0.5) - 1
  if (K >= 2) {
    beta_K = rbeta(num, shape1=1/2, shape2=(K-1)/2)
  } else {
    beta_K = 1
  }
  L = chi_aK * S * sqrt(beta_K)
  return(L)
}



