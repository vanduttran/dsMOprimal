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
#' @param opal An opal object or list of opal objects.
#' @param expr An encoded expression to evaluate.
#' @param wait See opal::datashield.aggreate options. Default: FALSE.
#' @param async See opal::datashield.aggreate options. Default: TRUE.
#' @export
crossAggregate <- function(opal, expr, wait = F, async = T) {
    expr <- dsCDISC:::.decode.arg(expr)
    opal::datashield.aggregate(opal, as.symbol(expr), wait, async)
}


#' @title Cross assign
#'
#' Cross assign # rather on client side
#' @param opal An opal object or list of opal objects.
#' @param symbol Name of an R symbol.
#' @param value An encoded expression with allowed assign function calls.
#' @param wait See opal::datashield.aggreate options. Default: FALSE.
#' @param async See opal::datashield.aggreate options. Default: TRUE.
#' @export
crossAssign <- function(opal, symbol, value, wait = F, async = T) {
    value <- dsCDISC:::.decode.arg(value)
    opal::datashield.assign(opal, symbol, as.symbol(value), wait, async)
}
