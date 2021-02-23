#' @title Row means
#'
#' Row means of linear transformations of a matrix 
#' @param x A numeric matrix
#' @param y A list of symmetric matrices of dimension (ncol(x), ncol(x))
#' @return Row means of x if y = NULL, or row means x %*% yy for each matrix yy in y otherwise
#' @export
rowmeans <- function(x, y = NULL) {
    ## if (is.null(dim(x)) || min(dim(x)) < 10) {
    ##     stop("x should be a matrix with two dimensions higher than 10.")
    ## }
    if (is.null(y)) {
        return (matrix(rowMeans(x), ncol=1, dimnames=list(rownames(x), "mean")))
    } else if (!all(sapply(y, isSymmetric))) {
        stop("y is not all symmetric.")
    } else {
        return (lapply(y, function(yy) matrix(rowMeans(tcrossprod(x, yy)), ncol=1, dimnames=list(rownames(x), "mean"))))
    }
}


#' @title Column means, deprecated
#'
#' Column means of a matrix 
#' @param x A numeric matrix
#' @return Column means of x
#' @export
colmeans <- function(x) {
    return (matrix(colMeans(x), nrow=1, dimnames=list("mean", colnames(x))))
}


#' @title Matrix centering
#'
#' Center matrix columns to 0
#' @param x A numeric matrix or data frame.
#' @return A centered matrix with column mean = 0
#' @export
center <- function(x, na.rm = TRUE) {
    y <- apply(x, c(1,2), as.numeric)
    if (na.rm) {
        y <- y[!is.na(rowSums(y)), , drop=F]
    } else {
        y[is.na(y)] <- 0
    }
    y <- y[order(rownames(y)), ]
    return (scale(y, center=TRUE, scale=FALSE))
}


#' @title Product of x' and first column of xx'
#'
#' Product of x' and first column of xx'
#' @param x A numeric matrix
#' @return t(x) %*% (x %*% t(x))[,1]
#' @export
singularProd <- function(x) {
    return (crossprod(tcrossprod(x)[, 1, drop=F], x))
}


#' @title Matrix cross product
#' 
#' Calculates the cross product t(x) \%*\% x
#' @param x A numeric matrix
#' @return t(x) \%*\% x
#' @export
#crossProd <- function(x) {
    ## if (is.null(dim(x)) || min(dim(x)) < 10) {
    ##     stop("x should be a matrix with two dimensions higher than 10.")
    ## }
#    return (crossprod(x))
#}
crossProd <- function(x, y = NULL) {
    ## if (is.null(dim(x)) || min(dim(x)) < 10) {
    ##     stop("x should be a matrix with two dimensions higher than 10.")
    ## }
    if (is.null(y)) {print(x); return (crossprod(x))}
    yd <- dsSwissKnife:::.decode.arg(y)
    if (is.list(yd)) yd <- do.call(rbind, yd)
    cat("x: ", dim(x), "\n")
    cat("y: ", dim(yd), "\n")
    return (crossprod(x, yd)) #(lapply(y, function(yy) matrix(crossprod(x, yy))))
}


#' @title Matrix cross product
#' 
#' Calculates the cross product x \%*\% t(x)
#' @param x A numeric matrix
#' @return x \%*\% t(x)
#' @export
#tcrossProd <- function(x) {
    ## if (is.null(dim(x)) || min(dim(x)) < 10) {
    ##     stop("x should be a matrix with two dimensions higher than 10.")
    ## }
#    return (tcrossprod(x))
#}
tcrossProd <- function(x, y = NULL) {
    ## if (is.null(dim(x)) || min(dim(x)) < 10) {
    ##     stop("x should be a matrix with two dimensions higher than 10.")
    ## }
    #yd <- dsSwissKnife:::.decode.arg(y)
    if (is.null(y)) return (tcrossprod(x))
    return (lapply(y, function(yy) matrix(tcrossprod(x, yy))))
}


#' @title Matrix triple product
#'
#' Calculate the triple product
#' @param x A numeric matrix
#' @param y A list of symmatric numeric matrices of dimension (ncol(x), ncol(x))
#' @return List of x %*% y %*% t(x)
#' @export
tripleProd <- function(x, y) {
    print(head(y[[1]]))
    print(class(y[[1]]))
    if (!all(sapply(y, isSymmetric))) {
        stop("y is not all symmetric.")
    } else {
        return (lapply(y, function(yy) tcrossprod(x, tcrossprod(x, yy))))
    }
}


#' @title Cross login
#'
#' Cross login # rather on client side
#' @param logins An encoded dataframe with server, url, user, password, and driver fields.
#' @export
crossLogin <- function(logins) {
    require(DSOpal)
    loginfo <- dsSwissKnife:::.decode.arg(logins)
    myDf <- data.frame(server=loginfo$server,
                       url=loginfo$url,
                       user=loginfo$user,
                       password=loginfo$password,
                       driver=loginfo$driver)
    DSI::datashield.login(myDf)
    #x <- tryCatch(DSI::datashield.login(myDf), error=function(e) return (sessionInfo()))
    #save(x, file = '/srv_local/session.Rdata')
}


#' @title Cross aggregate
#'
#' Cross aggregate
#' @param opal A list of opal objects.
#' @param expr An encoded expression to evaluate.
#' @param wait See DSI::datashield.aggreate options. Default: FALSE.
#' @param async See DSI::datashield.aggreate options. Default: TRUE.
#' @import DSI
#' @export
#crossAggregate <- function(opal, expr, wait = F, async = T) {
#    expr <- dsSwissKnife:::.decode.arg(expr)
#    ## only allow: crossProd, singularProd
#    stopifnot(grepl("^crossProd\\(|^singularProd\\(", dsSwissKnife:::.decode.arg(expr)))
#    DSI::datashield.aggregate(conns=opal, expr=as.symbol(expr), async=async)
#}
crossAggregate <- function(opal, expr, wait = F, async = T) {
    expr <- dsSwissKnife:::.decode.arg(expr)
    if (grepl("^as.call", expr)) {
        expr <- eval(str2expression(expr))
        stopifnot(!is.call(expr))
        DSI::datashield.aggregate(conns=opal, expr=expr, async=async)
    } else {
        ## only allow: crossProd, singularProd
        stopifnot(grepl("^crossProd\\(|^singularProd\\(", dsSwissKnife:::.decode.arg(expr)))
        DSI::datashield.aggregate(conns=opal, expr=as.symbol(expr), async=async)
    }
}


#' @title Cross assign
#'
#' Cross assign 
#' @param opal A list of opal objects.
#' @param symbol Name of an R symbol.
#' @param value A variable name or an R epxression with allowed assign function calls.
#' @param value.call A logical value, TRUE if value is function call, FALSE if value is a variable name.
#' @param wait See DSI::datashield.aggreate options. Default: FALSE.
#' @param async See DSI::datashield.aggreate options. Default: TRUE.
#' @import DSI
#' @export
crossAssign <- function(opal, symbol, value, value.call, variables = NULL, wait = F, async = T) {
    value <- dsSwissKnife:::.decode.arg(value)
    variables <- dsSwissKnife:::.decode.arg(variables)
    DSI::datashield.assign(conns=opal, symbol=symbol, value=ifelse(value.call, as.symbol(value), value), variables=variables, async=async)
}


#' @title Cross push
#' @export
pushValue <- function(value) {
    print(value)
    print(dsSwissKnife:::.decode.arg(value))
    pid <- Sys.getpid()
    save(value, file=paste0("/tmp/RtmpBCgKyI/",pid))
    return (pid)
}


#' @title Push
#'
#' Push the output of a function call to other servers 
#' @param opal A list of opal objects.
#' @export
push <- function(opal, symbol, value, value.call, variables = NULL, wait = F, async = T) {
    value <- dsSwissKnife:::.decode.arg(value)
    if (value.call && is.list(value)) {
	value <- do.call(rbind, value)
    	print(value)
    }
    variables <- dsSwissKnife:::.decode.arg(variables)
    DSI::datashield.assign(conns=opal, symbol=symbol, value=ifelse(value.call, as.symbol(value), value), variables=variables, async=async)
    datashield.errors()

    return (NULL)
}

