% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/server_func.R
\name{dscPush}
\alias{dscPush}
\title{Description of a pushed value}
\usage{
dscPush(conns, expr, async = T)
}
\arguments{
\item{conns}{A list of Opal connections.}

\item{expr}{An encoded expression to evaluate.}

\item{async}{See DSI::datashield.aggregate options. Default, TRUE.}
}
\value{
Returned value of given expression on opal
}
\description{
Description of a pushed value
}
