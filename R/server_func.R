
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
partitionMatrix <- function(x, seprow, sepcol=seprow) {
    print("HERE")
    print(class(x))
    print(summary(x))
    print(x)
    stopifnot(sum(seprow)==nrow(x) && sum(sepcol)==ncol(x))
    csseprow <- cumsum(seprow)
    indrow <- mclapply(1:length(seprow), mc.cores=min(length(seprow), detectCores()),  function(i) {
        return (c(ifelse(i==1, 0, csseprow[i-1])+1, csseprow[i]))
    })
    cssepcol <- cumsum(sepcol)
    indcol <- mclapply(1:length(sepcol), mc.cores=min(length(sepcol), detectCores()),  function(i) {
        return (c(ifelse(i==1, 0, cssepcol[i-1])+1, cssepcol[i]))
    })
    parMat <- mclapply(1:length(indrow), mc.cores=min(length(sepcol), detectCores()), function(i) {
        lapply(ifelse(isSymmetric(x), i, 1):length(indcol), function(j) {
            return (x[indrow[[i]][1]:indrow[[i]][2], indcol[[j]][1]:indcol[[j]][2]])
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
crossProdnew <- function(x, y = NULL, chunk = 500) {
    nblocksrow <- ceiling(ncol(x)/chunk)
    sepblocksrow <- rep(ceiling(ncol(x)/nblocksrow), nblocksrow-1)
    sepblocksrow <- c(sepblocksrow, ncol(x) - sum(sepblocksrow))
    if (is.null(y)) {
        tcpblocks <- partitionMatrix(crossprod(x), seprow=sepblocksrow)
        return (lapply(tcpblocks, function(tcpb) {
            return (lapply(tcpb, function(tcp) {
                .encode.arg(tcp)
            }))
        }))
    } else {
        nblockscol <- ceiling(ncol(y)/chunk)
        sepblockscol <- rep(ceiling(ncol(y)/nblockscol), nblockscol-1)
        sepblockscol <- c(sepblockscol, ncol(y) - sum(sepblockscol))
        tcpblocks <- partitionMatrix(crossprod(x, y), seprow=sepblocksrow, sepcol=sepblockscol)
        return (lapply(tcpblocks, function(tcpb) {
            return (lapply(tcpb, function(tcp) {
                .encode.arg(tcp)
            }))
        }))
    }
    # else {
    #     return (lapply(y, function(yy) {
    #         nblockscol <- ceiling(ncol(yy)/chunk)
    #         sepblockscol <- rep(ceiling(ncol(yy)/nblockscol), nblockscol-1)
    #         sepblockscol <- c(sepblockscol, ncol(yy) - sum(sepblockscol))
    #         tcpblocks <- partitionMatrix(crossprod(x, yy), seprow=sepblocksrow, sepcol=sepblockscol)
    #         return (lapply(tcpblocks, function(tcpb) {
    #             return (lapply(tcpb, function(tcp) {
    #                 .encode.arg(tcp)
    #             }))
    #         }))
    #     }))
    # }
}


#' @title Matrix cross product
#' @description Calculates the cross product x \%*\% t(y)
#' @param x A numeric matrix
#' @param y A list of numeric matrices. Default, y = x.
#' @return \code{x \%*\% t(y)}
#' @export
tcrossProd <- function(x, y = NULL, chunk = 500) {
    if (is.null(y)) {
        nblocks <- ceiling(nrow(x)/chunk)
        sepblocks <- rep(ceiling(nrow(x)/nblocks), nblocks-1)
        sepblocks <- c(sepblocks, nrow(x) - sum(sepblocks))
        tcpblocks <- partitionMatrix(tcrossprod(x), seprow=sepblocks)
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
        if (length(unique(sapply(uptcp, ncol)))==1) {
            tcp <- do.call(rbind, uptcp)
        } else {
            ## without the first layer of blocks
            no1tcp <- lapply(2:length(uptcp), function(i) {
                cbind(do.call(cbind, lapply(1:(i-1), function(j) {
                    t(matblocks[[j]][[i-j+1]])
                })), uptcp[[i]])
            })
            ## with the first layer of blocks
            tcp <- rbind(uptcp[[1]], do.call(rbind, no1tcp))
            rm(list=c("no1tcp"))
        }
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
    valued <- dsSwissKnife:::.decode.arg(value)
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
        #stopifnot(isSymmetric(tcp))
        dscbigmatrix <- describe(as.big.matrix(tcp))
        rm(list=c("tcp"))
    }
    return (dscbigmatrix)
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
#' @param wait See DSI::datashield.aggregate options. Default: FALSE.
#' @param async See DSI::datashield.aggregate options. Default: TRUE.
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
#' @param async See DSI::datashield.aggregate options. Default: TRUE.
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
#' @param value A variable name or an R expression with allowed assign function calls.
#' @param value.call A logical value, TRUE if value is function call, FALSE if value is a variable name.
#' @param wait See DSI::datashield.aggregate options. Default: FALSE.
#' @param async See DSI::datashield.aggregate options. Default: TRUE.
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
#' @param dsc A list of big memory descriptions
#' @return Sum of x and those stored in dsc
#' @import bigmemory
#' @keywords internal
sumMatrices <- function(dsc = NULL) {
    dscmat <- lapply(dsc, function(dscblocks) {
        y <- as.matrix(attach.big.matrix(dscblocks))
        #stopifnot(isSymmetric(y))
        return (y)
    })
    return (Reduce("+", dscmat))
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
#' @param loginFD Login information of the FD server (one of the servers containing cohort data)
#' @param logins Login information of other servers containing cohort data
#' @param querytable Name of table references in data repositories
#' @param queryvariables List of variable sets from the table references.
#' @param nameFD Name of the server to federate, among those in logins. Default, the first one in logins.
#' @return Covariance matrix of the virtual cohort
#' @import DSI parallel bigmemory
#' @export
federateCov <- function(loginFD, logins, querytable, queryvariables) {
    require(DSOpal)
    stopifnot((length(queryvariables) %in% c(1,2)) && (length(querytable) %in% c(1,2)))
    if (length(querytable)==1) querytable <- rep(querytable, length(queryvariables))
    
    loginFDdata    <- dsSwissKnife:::.decode.arg(loginFD)
    logindata      <- dsSwissKnife:::.decode.arg(logins)
    
    ## assign Cov matrix on each individual server
    opals <- DSI::datashield.login(logins=logindata)
    nNode <- length(opals)
    if (length(queryvariables)==1) {
        DSI::datashield.assign(opals, "rawData", querytable[[1]], variables=queryvariables[[1]], async=T)
        DSI::datashield.assign(opals, "centeredData", as.symbol('center(rawData)'), async=T)
        DSI::datashield.assign(opals, "crossProdSelf", as.symbol('crossProdnew(centeredData, chunk=50)'), async=T)
    } else {
        DSI::datashield.assign(opals, "rawData1", querytable[[1]], variables=queryvariables[[1]], async=T)
        DSI::datashield.assign(opals, "rawData2", querytable[[2]], variables=queryvariables[[2]], async=T)
        DSI::datashield.assign(opals, "centeredData1", as.symbol('center(rawData1)'), async=T)
        DSI::datashield.assign(opals, "centeredData2", as.symbol('center(rawData2)'), async=T)
        DSI::datashield.assign(opals, "crossProdSelf", as.symbol('crossProdnew(centeredData1, centeredData2, chunk=50)'), async=T)
    }
    
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
    rescov <- sumMatrices(crossProdSelfDSC)
    gc(reset=F)
    return (rescov)
}

#' @title Federated PCA
#' @description Perform the principal component analysis for the virtual cohort
#' @param loginFD Login information of the FD server (one of the servers containing cohort data)
#' @param logins Login information of other servers containing cohort data
#' @param querytab Encoded name of a table reference in data repositories.
#' @param queryvar Encoded value of a list of a variable set from the table references.
#' @param nameFD Name of the server to federate, among those in logins. Default, the first one in logins.
#' @return PCA object
#' @import DSI parallel bigmemory
#' @export
federatePCA <- function(loginFD, logins, querytab, queryvar) {
    querytable     <- dsSwissKnife:::.decode.arg(querytab)
    queryvariables <- dsSwissKnife:::.decode.arg(queryvar)
    covmat <- federateCov(loginFD, logins, querytable, queryvariables)
    return (princomp(covmat=covmat))
}


#' @title Federated RCCA
#' @description Perform the regularized canonical correlation analysis for the virtual cohort
#' @param loginFD Login information of the FD server (one of the servers containing cohort data).
#' @param logins Login information of servers containing cohort data.
#' @param querytab Encoded value of a vector containing names of one or two table references in data repositories. If one, 
#' the variables will be taken from the common table. Otherwise, the two sets of variables 
#' will be taken from corresponding tables.
#' @param queryvar Encoded value of a list of two variable sets from the table references.
#' @return RCCA object
#' @import DSI parallel bigmemory CCA
#' @export
federateRCCA <- function(loginFD, logins, querytab, queryvar) {
    querytable     <- dsSwissKnife:::.decode.arg(querytab)
    queryvariables <- dsSwissKnife:::.decode.arg(queryvar)
    stopifnot(length(queryvariables)==2 && (length(querytable) %in% c(1,2)))
    ## if only one table is given, it is duplicated
    if (length(querytable)==1) querytable <- rep(querytable, 2)
    Cxx <- federateCov(loginFD, logins, querytable[1], queryvariables[1])
    Cyy <- federateCov(loginFD, logins, querytable[2], queryvariables[2])
    Cxy <- federateCov(loginFD, logins, querytable, queryvariables)
    return(list(Cxx=Cxx, Cyy=Cyy, Cxy=Cxy))
    res <- geigen(Cxy, Cxx, Cyy)
    names(res) <- c("cor", "xcoef", "ycoef")
    return (res)
}
