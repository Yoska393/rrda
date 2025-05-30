% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cv.R
\name{MSE_lambda_rank}
\alias{MSE_lambda_rank}
\title{Compute MSE for different ranks of the coefficient Bhat and lambda.}
\usage{
MSE_lambda_rank(Bhat_comp_1L, X, Y, nrank)
}
\arguments{
\item{Bhat_comp_1L}{A list containing components of Bhat.}

\item{X}{A numeric matrix of explanatory variables.}

\item{Y}{A numeric matrix of response variables.}

\item{nrank}{A numeric vector indicating the rank(s) of Bhat.}
}
\value{
A numeric vector of MSE values for each rank.
}
\description{
For cross-validation, compute MSE for different ranks of the coefficient Bhat and lambda.
}
