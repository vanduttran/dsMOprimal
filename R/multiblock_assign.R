#' @title Random matrix
#'
#' Create a random data matrix
#' @param nrow Number of rows
#' @param ncol Number of colums
#' @return A matrix of dimension (nrow, ncol)
#' @importFrom stats runif
#' @export
randomData <- function(nrow, ncol) {
    return (matrix(runif(nrow*ncol, min=0, max=1), nrow, ncol))
}

#' @title Trace
#'
#' Trace of a matrix
#' @param x A symmetric numeric matrix
#' @return tr(x)
tr <- function(x) {
    if (nrow(x) != ncol(x)) {
    	stop("x should be symmetric.")
    }
    return (sum(diag(x)))
}

#' @title Matrix scaling
#'
#' Scale a matrix to trace of 1
#' @param x A numeric matrix
#' @return A scaled matrix with tr(xx') = 1
#' @export
trscale <- function(x) {
    return (x/sqrt(tr(crossprod(x))))
}
