% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/server_func.R
\name{crossLogin}
\alias{crossLogin}
\title{Cross login}
\usage{
crossLogin(logins)
}
\arguments{
\item{logins}{An encoded dataframe with server, url, user, password, driver,
and options fields.}
}
\value{
Object(s) of class DSConnection.
}
\description{
Call datashield.login on remote servers
}
