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
#' @param x A symmetric numeric matrix or data frame.
#' @return tr(x)
tr <- function(x) {
    if (nrow(x) != ncol(x)) {
    	stop("x should be symmetric.")
    }
    y <- apply(x, c(1,2), as.numeric)
    return (sum(diag(y)))
}

#' @title Matrix scaling
#'
#' Scale a matrix to trace of 1
#' @param x A numeric matrix or data frame.
#' @return A scaled matrix with tr(xx') = 1
#' @export
trscale <- function(x) {
    y <- apply(x, c(1,2), as.numeric)
    return (y/sqrt(tr(crossprod(y))))
}

#' @title Cross login
#'
#' Cross login # rather on client side
#' @param logins An encoded dataframe with server, url, user, password, and table fields.
#' 
#' @export
crossLogin <- function(logins) {
    loginfo <- dsCDISC:::.decode.arg(logins)
    myDf <- data.frame(server=loginfo$server,
                       url=loginfo$url,
                       user=loginfo$user,
                       password=loginfo$password,
                       table=loginfo$table)
    opal::datashield.login(myDf)
}


#' @title Cross aggregate
#'
#' Cross aggregate # rather on client side
#' @param servers
#' @param func
#' @export
crossAggregate <- function(servers, func, wait = F, async = T) {
    func <- dsCDISC:::.decode.arg(func)
    opal::datashield.aggregate(servers, as.symbol(func), wait, async)
}
