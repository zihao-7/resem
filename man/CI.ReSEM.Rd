% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ReSEM.R
\name{CI.ReSEM}
\alias{CI.ReSEM}
\title{Confidence interval under ReSEM}
\usage{
CI.ReSEM(assignment, Y, X, W, C, E, alpha = 0.05, design)
}
\arguments{
\item{assignment}{A vector of length N, which contains -1, 0, 1 values. We use 1 for treated units, 0 for control units, and -1 for units that will not enter the experiment.}

\item{Y}{A vector of length N, which contains observed outcomes for all units. Will be NA for units that will not enter the experiment.}

\item{X}{An n by K matrix which contains K dimensional covariate vector for sampled units, available at the treatment assignment stage.}

\item{W}{An N by J matrix which contains J dimensional covariate vector for each unit, available at the sampling stage.}

\item{C}{An n by K' matrix which contains K' dimensional covariate vector for sampled units, available at the analysis stage. Assume X is a subset of C.}

\item{E}{An N by J' matrix which contains J' dimensional covariate vector for each unit, available at the analysis stage. Assume W is a subset of E.}

\item{alpha}{Significance level. The computed CI will have 1 - \code{alpha} confidence level .}

\item{design}{Rerandomization design used to compute the CI. Can be "ST", "S", "T", "CRSE", "ST.adjusted".}
}
\value{
1 - \code{alpha} confidence interval for the average treatment effect.
}
\description{
Computes confidence interval for the average treatment effect under ReSEM.
}
