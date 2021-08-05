
#' @title Ranking
#'
#' Ranking of features in each sample 
#' @param x A matrix or data frame, samples in rows and features in columns
#' @return Ranking of features
#' @export
dsRank <- function(x) {
    return (t(apply(x, 1, rank)))
}


#' @title Matrix dimension
#' @description Dimension of a data frame or matrix
#' @param x A matrix or data frame
#' @return Dimension of x
#' @export
dsDim <- function(x) {
    return (dim(x))
}


#' @title Row means, deprecated
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
#' @param na.rm A logical value indicating NA values should be removed. Default, FALSE, NA set to 0.
#' @return A centered matrix with column mean = 0
#' @export
center <- function(x, na.rm = FALSE) {
    #y <- apply(head(x[order(rownames(x)), ], 101), c(1,2), as.numeric)
    y <- apply(x, c(1,2), as.numeric)
    if (na.rm) {
        y <- y[!is.na(rowSums(y)), , drop=F]
    } else {
        y[is.na(y)] <- 0
    }
    y <- y[order(rownames(y)), ]
    y <- head(y, 101)
    return (scale(y, center=TRUE, scale=FALSE))
}


#' @title Matrix partition
#' @description Partition a symmetric matrix into square blocks
#' @param x A symmetric matrix
#' @param sep A numeric vectors indicating sizes of square blocks
#' @return List of blocks
#' @import parallel
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


#' @title Loadings
#' 
#' Loadings of features in a new feature basis
#' @param x A numeric matrix
#' @param y A numeric matrix for the new basis of the same nrow to x.
#' @param operator An operation to compute the loadings (\code{crossprod}, \code{cor}). Default, \code{crossprod}
#' @return operator(x, y)
#' @export
loadings <- function(x, y, operator = 'crossprod') {
    stopifnot(operator %in% c('crossprod', 'cor'))
    yd <- dsSwissKnife:::.decode.arg(y)
    if (is.list(yd)) yd <- do.call(rbind, yd)
    if (operator=='cor') return(cor(x, yd))
    return (crossprod(x, yd))
}


#' @title Matrix cross product
#' 
#' Calculates the cross product t(x) \%*\% x
#' @param x A numeric matrix
#' @return t(x) \%*\% x
#' @export
crossProd <- function(x, y = NULL) {
    ## if (is.null(dim(x)) || min(dim(x)) < 10) {
    ##     stop("x should be a matrix with two dimensions higher than 10.")
    ## }
    if (is.null(y)) {return (crossprod(x))}
    yd <- dsSwissKnife:::.decode.arg(y)
    if (is.list(yd)) yd <- do.call(rbind, yd)
    cat("x: ", dim(x), "\n")
    cat("y: ", dim(yd), "\n")
    return (crossprod(x, yd))
}


#' @title Matrix cross product
#' @description Calculates the cross product t(x) \%*\% y
#' @param x A numeric matrix
#' @param y A list of numeric matrices. Default, y = x.
#' @return \code{t(x) \%*\% y}
#' @export
crossProdnew <- function(x, y = NULL, chunk=500) {
    if (is.null(y)) {
        nblocks <- ceiling(ncol(x)/chunk)
        sepblocks <- rep(ceiling(ncol(x)/nblocks), nblocks-1)
        sepblocks <- c(sepblocks, ncol(x) - sum(sepblocks))
        tcpblocks <- partitionMatrix(crossprod(x), sep=sepblocks)
        return (lapply(tcpblocks, function(tcpb) {
            return (lapply(tcpb, function(tcp) {
                .encode.arg(tcp)
            }))
        }))
    }
    return (lapply(y, function(yy) .encode.arg(matrix(crossprod(x, yy)))))
}


#' @title Matrix cross product
#' @description Calculates the cross product x \%*\% t(y)
#' @param x A numeric matrix
#' @param y A list of numeric matrices. Default, y = x.
#' @return \code{x \%*\% t(y)}
#' @export
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


#' @title Rebuild matrix from blocks
#' @description Rebuild a matrix from its blocks
#' @param chunks List of list of encoded matrix blocks, obtained from crossProd or tcrossProd
#' @keywords internal
rebuildMatrix <- function(blocks) {
    ## decode matrix blocks
    matblocks <- mclapply(blocks, mc.cores=length(blocks), function(y) {
        mclapply(y, mc.cores=length(y), function(x) {
            return (do.call(rbind, dsSwissKnife:::.decode.arg(x)))
        })
    })
    uptcp <- lapply(matblocks, function(bl) do.call(cbind, bl))
    ## combine the blocks into one matrix
    if (length(uptcp)>1) {
        ## without the first layer of blocks
        no1tcp <- lapply(2:length(uptcp), function(i) {
            cbind(do.call(cbind, lapply(1:(i-1), function(j) {
                t(matblocks[[j]][[i-j+1]])
            })), uptcp[[i]])
        })
        ## with the first layer of blocks
        tcp <- rbind(uptcp[[1]], do.call(rbind, no1tcp))
        rm(list=c("no1tcp"))
    } else {
        tcp <- uptcp[[1]]
    }
    rm(list=c("matblocks", "uptcp"))
    return (tcp)
}


#' @title Push a symmetric matrix
#' @description Push symmetric matrix data into the federated server
#' @param value An encoded value to be pushed
#' @import bigmemory parallel
#' @return Description of the pushed value
#' @export
pushSymmMatrix <- function(value) {
    print("symmetric")
    valued <- dsSwissKnife:::.decode.arg(value)
    print("decoded")
    stopifnot(is.list(valued) && length(valued)>0)
    if (FALSE) {#is.list(valued[[1]])) {
        dscbigmatrix <- mclapply(valued, mc.cores=min(length(valued), detectCores()), function(x) {
            x.mat <- do.call(rbind, x)
            stopifnot(ncol(x.mat)==1)
            return (describe(as.big.matrix(x.mat)))
        })
    } else {
        # dscbigmatrix <- mclapply(valued, mc.cores=length(valued), function(y) {
        #     ## N.B. mclapply with length(y) cores allows allocating memory for all blocks. 
        #     ##      or only last mc.cores blocks are allocated.
        #     ##      lapply allocates memory only for the last block in the list.
        #     return (mclapply(y, mc.cores=length(y), function(x) {
        #         x.mat <- do.call(rbind, dsSwissKnife:::.decode.arg(x))
        #         return (describe(as.big.matrix(x.mat)))
        #     }))
        # })
        ## Possible solution: Rebuild the whole matrix here, and return its only allocation
        # matblocks <- mclapply(valued, mc.cores=length(valued), function(y) {
        #     mclapply(y, mc.cores=length(y), function(x) {
        #         return (do.call(rbind, dsSwissKnife:::.decode.arg(x)))
        #     })
        # })
        # rm(list=c("valued"))
        # uptcp <- lapply(matblocks, function(bl) do.call(cbind, bl))
        # ## combine the blocks into one matrix
        # if (length(uptcp)>1) {
        #     ## without the first layer of blocks
        #     no1tcp <- lapply(2:length(uptcp), function(i) {
        #         cbind(do.call(cbind, lapply(1:(i-1), function(j) {
        #             t(matblocks[[j]][[i-j+1]])
        #         })), uptcp[[i]])
        #     })
        #     ## with the first layer of blocks
        #     tcp <- rbind(uptcp[[1]], do.call(rbind, no1tcp))
        #     rm(list=c("no1tcp"))
        # } else {
        #     tcp <- uptcp[[1]]
        # }
        tcp <- rebuildMatrix(valued)
        stopifnot(isSymmetric(tcp))
        dscbigmatrix <- describe(as.big.matrix(tcp))
        rm(list=c("tcp"))
    }
    return (dscbigmatrix)
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
#' @return Returned value of given expression on opal
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


#' @title Sum matrices
#' @description Compute the sum of a matrix and those stored in bigmemory
#' @param x A symmetric matrix
#' @param dsc A list of big memory descriptions
#' @return Sum of x and those stored in dsc
#' @import bigmemory
#' @export
sumMatrices <- function(x, dsc = NULL) {
    stopifnot(isSymmetric(x))
    dsclist <- dsc #dsSwissKnife:::.decode.arg(dsc)
    print(length(dsclist))
    print(names(dsclist))
    dscmat <- lapply(dsclist, function(dscblocks) {
        print(class(dscblocks))
        y <- as.matrix(attach.big.matrix(dscblocks))
        stopifnot(isSymmetric(y))
        return (y)
    })
    return (Reduce("+", c(x, dscmat)))
}


#' @title Encode function  arguments
#' @description Serialize to JSON, then encode base64,
#'  then replace '+', '/' and '=' in the result in order to play nicely with the opal sentry.
#'  Used to encode non-scalar function arguments prior to sending to the opal server.
#'  There's a corresponding function in the server package calle .decode_args
#' @param some.object the object to be encoded
#' @return encoded text with offending characters replaced by strings
#' @keywords internal
.encode.arg <- function(some.object){
  encoded <- RCurl::base64Encode(jsonlite::toJSON(some.object, null = 'null'));
  # go fishing for '+', '/' and '=', opal rejects them :
  my.dictionary <- c('\\/' = '-slash-', '\\+' = '-plus-', '\\=' = '-equals-')
  sapply(names(my.dictionary), function(x){
    encoded[1] <<- gsub(x, my.dictionary[x], encoded[1])
  })
  return(paste0(encoded[1],'base64'))

}


#' @title Federated covariance matrix
#' @description Compute the covariance matrix for the virtual cohort
#' @param logins Login information of the servers containing cohort data
#' @param querytab Encoded name of a table reference in data repositories
#' @param queryvar Encoded variables from the table reference
#' @param nameFD Name of the server to federate, among those in logins. Default, the first one in logins.
#' @import DSI parallel bigmemory
#' @export
federateCov <- function(x, loginFD, logins, querytab, queryvar) {
    require(DSOpal)
    loginFDdata    <- dsSwissKnife:::.decode.arg(loginFD)
    logindata      <- dsSwissKnife:::.decode.arg(logins)
    querytable     <- dsSwissKnife:::.decode.arg(querytab)
    queryvariables <- dsSwissKnife:::.decode.arg(queryvar)

    ## assign Cov matrix on each individual server
    opals <- DSI::datashield.login(logins=logindata)
    nNode <- length(opals)
    DSI::datashield.assign(opals, "rawData", querytable, variables=queryvariables, async=T)
    DSI::datashield.assign(opals, "centeredData", as.symbol('center(rawData)'), async=T)
    DSI::datashield.assign(opals, "crossProdSelf", as.symbol('crossProdnew(centeredData, chunk=50)'), async=T)
    
    ## push data from non-FD servers to FD-assigned server: user and password for login between servers are required
    loginFDdata$user     <- loginFDdata$userserver
    loginFDdata$password <- loginFDdata$passwordserver
    DSI::datashield.assign(opals, 'FD', as.symbol(paste0("crossLogin('", .encode.arg(loginFDdata), "')")), async=T)
    command <- paste0("dscPush(FD, '", 
                      .encode.arg(paste0("as.call(list(as.symbol('pushSymmMatrix'), dsSSCP:::.encode.arg(crossProdSelf)", "))")), 
                      "', async=T)")
    cat("Command: ", command, "\n")
    crossProdSelfDSC <- DSI::datashield.aggregate(opals, as.symbol(command), async=T)
    crossProdSelfDSC <- mclapply(crossProdSelfDSC, mc.cores=min(length(crossProdSelfDSC), detectCores()), function(dscblocks) {
        return (dscblocks[[1]])
    })
    print("OK")
    return (sumMatrices(rebuildMatrix(x), crossProdSelfDSC))
    
    crossProdSelf <- mclapply(crossProdSelfDSC, mc.cores=min(length(crossProdSelfDSC), detectCores()), function(dscblocks) {
        print(dscblocks)
        return (as.matrix(attach.big.matrix(dscblocks[[1]]))) ## it is on server3, cannot read memory on server1!!!
        ## retrieve the blocks as matrices: on FD
        matblocks <- lapply(dscblocks[[1]], function(dscblock) {
            lapply(dscblock, function(dsc) {
                as.matrix(attach.big.matrix(dsc))
            })
        })
        uptcp <- lapply(matblocks, function(bl) do.call(cbind, bl))
        ## combine the blocks into one matrix
        if (length(uptcp)>1) {
            ## without the first layer of blocks
            no1tcp <- lapply(2:length(uptcp), function(i) {
                cbind(do.call(cbind, lapply(1:(i-1), function(j) {
                    t(matblocks[[j]][[i-j+1]])
                })), uptcp[[i]])
            })
            ## with the first layer of blocks
            tcp <- rbind(uptcp[[1]], do.call(rbind, no1tcp))
        } else {
            tcp <- uptcp[[1]]
        }
        stopifnot(isSymmetric(tcp))
        return (tcp)
    })
    gc(reset=F)
    return (crossProdSelf)
    # logindata.FD <- logindata[logindata$server == nameFD, , drop=F]
    # logindata.FD$user <- logindata.FD$userserver
    # logindata.FD$password <- logindata.FD$passwordserver
    # ##- received by FD from other nodes ----
    # invisible(mclapply(setdiff(names(opals), nameFD), mc.cores=1, function(opn) {
    #     opals.loc <- paste0("crossLogin('", .encode.arg(logindata.FD), "')")
    #     datashield.assign(opals[opn], 'mates', as.symbol(opals.loc), async=F)
    #     
    #     command.opn <- paste0("crossAggregate(mates, '", 
    #                           .encode.arg(paste0("as.call(list(as.symbol('pushValue'), dsSSCP:::.encode.arg(crossProdSelf), dsSSCP:::.encode.arg('", opn, "')))")), 
    #                           "', async=F)")
    #     cat("Command: ", command.opn, "\n")
    #     print(datashield.assign(opals[opn], "pidMate", as.symbol(command.opn), async=F))
    # }))
    # datashield.symbols(opals)
    
    #-----
    
    
}


