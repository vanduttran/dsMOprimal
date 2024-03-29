% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/server_func.R
\name{.estimateR}
\alias{.estimateR}
\title{RCCA tuning}
\usage{
.estimateR(
  loginFD,
  logins,
  funcPreProc,
  querytables,
  mc.cores = 1,
  nfold = 5,
  grid1 = seq(0.001, 1, length = 5),
  grid2 = seq(0.001, 1, length = 5)
)
}
\arguments{
\item{loginFD}{Login information of the FD server (one of the servers containing cohort data).}

\item{logins}{Login information of servers containing cohort data.}

\item{mc.cores}{Number of cores for parallel computing. Default: 1}

\item{grid1}{Checking values for \code{lambda1}.}

\item{grid2}{Checking values for \code{lambda2}.}

\item{func}{Encoded definition of a function for preparation of raw data matrices. 
Two arguments are required: conns (list of DSConnection-classes), 
symbol (names of the two R symbols) (see datashield.assign).}

\item{symbol}{Encoded vector of names of the two R symbols to assign in the Datashield R session on each server in \code{logins}.
The two assigned R variables will be used as the input raw data to compute covariance matrices for CCA.
Other assigned R variables in \code{func} are ignored.}

\item{lambda1}{Regularized parameter value for first data set. Default, 0.}

\item{lambda2}{Regularized parameter value for second data set. Default, 0.}

\item{tune}{Logical value indicating whether the tuning for lambda values will be performed. Default, FALSE, no tuning.}

\item{tune_param}{Tuning parameters. \code{nfold} n-fold cross-validation.}
}
\value{
Optimal values of \code{lambda1} and \code{lambda2}.
}
\description{
Estimate optimized parameters of regulation lambda1 and lambda2
}
\keyword{internal}
