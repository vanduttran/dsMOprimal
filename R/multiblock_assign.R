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

#' @title Cross login
#'
#' Cross login # rather on client side
#' @param server ..
#' @param url
#' 
#' @export
crossLogin <- function(server, url) {
    url <- dsCDISC:::.decode.arg(url)
    myDf <- data.frame(server=server,
                       url=url,
                       user='test',
                       password='test123',
                       table='test.TEST')
    opal::datashield.login(myDf)
}


#' @title Cross aggregate
#'
#' Cross aggregate # rather on client side
#' @param server
#' @param func
#' @export
crossAggregate <- function(servers, func, wait = F, async = T) {
    func <- dsCDISC:::.decode.arg(func)
    opal::datashield.aggregate(servers, as.symbol(func), wait, async)
}
