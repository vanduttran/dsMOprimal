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
    #y <- head(y, 171)
    return (scale(y, center=TRUE, scale=FALSE))
}


#' @title Matrix partition
#' @description Partition a symmetric matrix into square blocks
#' @param x A symmetric matrix
#' @param sep A numeric vectors indicating sizes of square blocks
#' @return List of blocks
#' @keywords internal
partitionMatrix <- function(x, sep) {
    cssep <- cumsum(sep)
    ind <- mclapply(1:length(sep), mc.cores=min(length(sep), detectCores()),  function(i) {
        return (c(ifelse(i==1, 0, cssep[i-1])+1, cssep[i]))
    })
    parMat <- mclapply(1:length(ind), mc.cores=min(length(sep), detectCores()), function(i) {
        lapply(i:length(ind), function(j) {
            return (x[ind[[i]][1]:ind[[i]][2], ind[[j]][1]:ind[[j]][2]])
        })
    })
    return (parMat)
}


#' @title Product of x' and first column of xx'
#' @description Product of x' and first column of xx'
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
    if (is.null(y)) {return (crossprod(x))}
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
tcrossProd <- function(x, y = NULL, chunk=500) {
    ## if (is.null(dim(x)) || min(dim(x)) < 10) {
    ##     stop("x should be a matrix with two dimensions higher than 10.")
    ## }
    #yd <- dsSwissKnife:::.decode.arg(y)
    #if (is.null(y)) return (tcrossprod(x))
    if (is.null(y)) {
        nblocks <- ceiling(nrow(x)/chunk)
        sepblocks <- rep(ceiling(nrow(x)/nblocks), nblocks-1)
        sepblocks <- c(sepblocks, nrow(x) - sum(sepblocks))
        tcpblocks <- partitionMatrix(tcrossprod(x), sep=sepblocks)
        return (lapply(tcpblocks, function(tcpb) {
            return (lapply(tcpb, function(tcp) {
                .encode.arg(tcp)
            }))
        }))
    }
    return (lapply(y, function(yy) .encode.arg(matrix(tcrossprod(x, yy)))))
}

    
#' @title 
pushTCrossProd <- function(x, y = NULL) {
    if (is.null(y)) return (tcrossprod(x))
    return (lapply(y, function(yy) matrix(tcrossprod(x, yy))))
}


#' @title Matrix triple product
#'
#' Calculate the triple product
#' @param x A numeric matrix
#' @param y A list of symmatric numeric matrices of dimension (ncol(x), ncol(x))
#' @return List of x %*% y %*% t(x)
#' @import bigmemory
#' @export
#tripleProdrm <- function(x, y) {
#    print(head(y[[1]]))
#    print(class(y[[1]]))
#    if (!all(sapply(y, isSymmetric))) {
#        stop("y is not all symmetric.")
#    } else {
#        return (lapply(y, function(yy) tcrossprod(x, tcrossprod(x, yy))))
#    }
#}
tripleProd <- function(x, pids) {
    pids <- dsSwissKnife:::.decode.arg(pids)
    tp <- lapply(pids, function(pid) {
        print(pid)
        if (file.exists(paste0("/tmp/",pid))) {
            load(paste0("/tmp/",pid))
            y <- as.matrix(attach.big.matrix(dscbigmatrix))
            #yd <- dsSwissKnife:::.decode.arg(value)
            #if (is.list(yd)) y <- do.call(rbind, yd)
            stopifnot(isSymmetric(y))
            return (tcrossprod(x, tcrossprod(x, y)))
        } else {
            return (NULL)
        }
    })
    names(tp) <- pids
    return (tp)
}
# tripleProd <- function(x, pids) {
#     pids <- dsSwissKnife:::.decode.arg(pids)
#     tp <- lapply(pids, function(pid) {
#         print(pid)
#         if (file.exists(paste0("/tmp/",pid))) {
#             load(paste0("/tmp/",pid))
#             yd <- dsSwissKnife:::.decode.arg(value)
#             if (is.list(yd)) y <- do.call(rbind, yd)
#             stopifnot(isSymmetric(y))
#             return (tcrossprod(x, tcrossprod(x, y)))
#         } else {
#             return (NULL)
#         }
#     })
#     names(tp) <- pids
#     return (tp)
# }

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
                       driver=loginfo$driver,
                       options=loginfo$options)
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
        print(expr)
        print(is.call(expr))
        stopifnot(is.call(expr))
        DSI::datashield.aggregate(conns=opal, expr=expr, async=async)
    } else {
        ## only allow: crossProd, singularProd
        print(expr)
        stopifnot(grepl("^crossProd\\(|^singularProd\\(", expr)) #dsSwissKnife:::.decode.arg(expr)))
        DSI::datashield.aggregate(conns=opal, expr=as.symbol(expr), async=async)
    }
}


#' @title Description of a pushed value
#' @description Description of a pushed value
#' @param opal A list of opal objects.
#' @param expr An encoded expression to evaluate.
#' @param async See DSI::datashield.aggreate options. Default: TRUE.
#' @import DSI
#' @export
dscPush <- function(opal, expr, async = T) {
  expr <- dsSwissKnife:::.decode.arg(expr)
  stopifnot(grepl("^as.call", expr))
  expr <- eval(str2expression(expr))
  return (DSI::datashield.aggregate(conns=opal, expr=expr, async=async))
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
#' @import bigmemory
#' @export
# pushValuesave <- function(value, name) {
#     print(value)
#     pid <- Sys.getpid()
#     save(value, file=paste0("/tmp/", dsSwissKnife:::.decode.arg(name)))
#     return (pid)
# }
pushValue <- function(value, name) {
    valued <- dsSwissKnife:::.decode.arg(value)
    if (is.list(valued)) valued <- do.call(rbind, valued)
    stopifnot(isSymmetric(valued))
    dscbigmatrix <- describe(as.big.matrix(valued))
    save(dscbigmatrix, file=paste0("/tmp/", dsSwissKnife:::.decode.arg(name)))
    return (dscbigmatrix)
}


#' @title Encode function  arguments
#' @description Serialize to JSON, then encode base64,
#'  then replace '+', '/' and '=' in the result in order to play nicely with the opal sentry.
#'  Used to encode non-scalar function arguments prior to sending to the opal server.
#'  There's a corresponding function in the server package calle .decode_args
#' @param some.object the object to be encoded
#' @return encoded text with offending characters replaced by strings
#' @keywords internal
#'
.encode.arg <- function(some.object){
  encoded <- RCurl::base64Encode(jsonlite::toJSON(some.object, null = 'null'));
  # go fishing for '+', '/' and '=', opal rejects them :
  my.dictionary <- c('\\/' = '-slash-', '\\+' = '-plus-', '\\=' = '-equals-')
  sapply(names(my.dictionary), function(x){
    encoded[1] <<- gsub(x, my.dictionary[x], encoded[1])
  })
  return(paste0(encoded[1],'base64'))

}


