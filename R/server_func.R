#' @title Ranking, deprecated
#' @description Ranking of features in each sample 
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


#' @keywords internal
.printTime <- function(message = "") {
    cat(message, "---", as.character(Sys.time()), "\n")
}


#' @title Row means, deprecated
#' @description Row means of linear transformations of a matrix 
#' @param x A numeric matrix
#' @param y A list of symmetric matrices of dimension (ncol(x), ncol(x))
#' @return Row means of x if y = NULL, or row means of x %*% yy for each matrix yy in y otherwise
#' @keywords internal
.rowmeans <- function(x, y = NULL) {
    if (is.null(y)) {
        return (matrix(rowMeans(x), ncol=1, dimnames=list(rownames(x), "mean")))
    } else if (!all(sapply(y, isSymmetric))) {
        stop("y is not all symmetric.")
    } else {
        return (lapply(y, function(yy) matrix(rowMeans(tcrossprod(x, yy)), ncol=1, dimnames=list(rownames(x), "mean"))))
    }
}


#' @title Column means, deprecated
#' @description Column means of a matrix 
#' @param x A numeric matrix
#' @return Column means of x
#' @keywords internal
.colmeans <- function(x) {
    return (matrix(colMeans(x), nrow=1, dimnames=list("mean", colnames(x))))
}


#' @title Assign rownames
#' @description Assign row names to a matrix 
#' @param x A numeric matrix
#' @param row.names An encoded vector of names with the length of nrow(x), or a name of variable in \code{envir}.
#' @param envir A data frame or environment where row.names will be queried if it is a variable name. Default, \code{.GlobalEnv}
#' @return Matrix with rownames
#' @export
setRowNames <- function(x, row.names, envir = .GlobalEnv) {
    if (length(row.names)==1) {
        if (grepl("base64$", row.names)) rn <- .decode.arg(row.names)
        else {
            if (typeof(envir)!='environment') envir <- as.environment(envir)
            rn <- get(row.names, envir=envir)
        }        
        if (nrow(x) == length(rn)) rownames(x) <- rn
        else stop("Cannot assign row.names of length different from nrow(x)")
    } else stop("Wrong row.names provided!")
    return (x)
}


#' @title Rownames
#' @description Row names of a matrix 
#' @param x A numeric matrix
#' @return rownames of x
#' @export
rowNames <- function(x) {
    return (rownames(x))
}


#' @title Colnames
#' @description Col names of a matrix 
#' @param x A numeric matrix
#' @return colnames of x
#' @export
colNames <- function(x) {
    return (colnames(x))
}


#' @title Matrix centering
#' @description Center matrix columns or rows to 0
#' @param x A numeric matrix or data frame. A list of matrices or data frames can also be provided,
#' where they will be concatenated to perform the centering.
#' @param subset Encoded value of an index vector indicating the subset of individuals to consider. 
#' Default, NULL, all individuals are considered.
#' @param byColumn A logical value indicating whether the input data is centered by column or row.
#' Default, TRUE, centering by column. Constant variables across samples are removed. 
#' If FALSE, centering and scaling by row. Constant samples across variables are removed.
#' @param scale  A logical value indicating whether the variables should be scaled to have unit variance. Default, FALSE.
#' @param na.rm A logical value indicating whether NA values should be removed. Default, FALSE, NA set to 0.
#' @return The centered matrix.
#' @export
center <- function(x, subset = NULL, byColumn = TRUE, scale = FALSE, na.rm = FALSE) {
    ## horizontally combine x, if x is a list of data frames
    if (is.list(x) && !is.data.frame(x)) {
        y <- apply(do.call(cbind, lapply(x, function(xx) xx[order(rownames(xx)), ])), c(1,2), as.numeric)
    } else {
        y <- apply(x, c(1,2), as.numeric)
    }
    ## missing values
    if (na.rm) {
        y <- y[!is.na(rowSums(y)), , drop=F]
    } else {
        y[is.na(y)] <- 0
    }
    ## constant value across samples or variables
    if (isTRUE(byColumn)) {
        y <- y[, apply(y, 2, function(yc) length(unique(yc))) != 1]
    } else {
        y <- y[apply(y, 1, function(yc) length(unique(yc))) != 1, ]
    }
    ## ordering
    y <- y[order(rownames(y)), ]

    ## subseting
    subset <- .decode.arg(subset)
    if (!is.null(subset)) y <- y[subset, , drop=F]
    
    if (isTRUE(byColumn)) return (scale(y, center=TRUE, scale=scale))
    else return (t(scale(t(y), center=TRUE, scale=TRUE)))
}


#' @title Matrix partition
#' @description Partition a symmetric matrix into square blocks
#' @param x A symmetric matrix
#' @param sep A numeric vectors indicating sizes of square blocks
#' @return List of blocks
#' @import parallel
#' @keywords internal
.partitionMatrix <- function(x, seprow, sepcol = seprow) {
    stopifnot(sum(seprow)==nrow(x) && sum(sepcol)==ncol(x))
    csseprow <- cumsum(seprow)
    indrow <- lapply(1:length(seprow), function(i) {
        return (c(ifelse(i==1, 0, csseprow[i-1])+1, csseprow[i]))
    })
    cssepcol <- cumsum(sepcol)
    indcol <- lapply(1:length(sepcol), function(i) {
        return (c(ifelse(i==1, 0, cssepcol[i-1])+1, cssepcol[i]))
    })
    parMat <- lapply(1:length(indrow), function(i) {
        lapply(ifelse(isSymmetric(x), i, 1):length(indcol), function(j) {
            return (x[indrow[[i]][1]:indrow[[i]][2], indcol[[j]][1]:indcol[[j]][2]])
        })
    })
    
    # indrow <- mclapply(1:length(seprow), mc.cores=max(2, min(length(seprow), detectCores())),  function(i) {
    #     return (c(ifelse(i==1, 0, csseprow[i-1])+1, csseprow[i]))
    # })
    # cssepcol <- cumsum(sepcol)
    # indcol <- mclapply(1:length(sepcol), mc.cores=max(2, min(length(sepcol), detectCores())),  function(i) {
    #     return (c(ifelse(i==1, 0, cssepcol[i-1])+1, cssepcol[i]))
    # })
    # parMat <- mclapply(1:length(indrow), mc.cores=max(2,min(length(sepcol), detectCores())), function(i) {
    #     lapply(ifelse(isSymmetric(x), i, 1):length(indcol), function(j) {
    #         return (x[indrow[[i]][1]:indrow[[i]][2], indcol[[j]][1]:indcol[[j]][2]])
    #     })
    # })
    
    return (parMat)
}


#' @title Singular product
#' @description Product of t(x) and first column of x \%*\% t(x)
#' @param x A numeric matrix
#' @return t(x) \%*\% (x \%*\% t(x))[,1]
#' @export
singularProd <- function(x) {
    return (crossprod(tcrossprod(x)[, 1, drop=F], x))
}


## TOCHECK: check y, check number of queries for security issue
#' @title Loadings
#' @description Loadings of features in a new feature basis
#' @param x A numeric matrix
#' @param y A numeric matrix for the new basis of the same nrow to x.
#' @param operator An operation to compute the loadings (\code{crossprod}, \code{cor}). Default, \code{crossprod}
#' @return operator(x, y)
#' @export
loadings <- function(x, y, operator = 'crossprod') {
    operator <- match.arg(operator, choices=c('crossprod', 'cor', "prod"))
    yd <- .decode.arg(y)
    if (is.list(yd)) yd <- do.call(rbind, yd)
    
    return (switch(operator,
                   cor=cor(x, yd),
                   prod=crossprod(t(x), yd),
                   crossprod=crossprod(x, yd)))
}


#' @title Matrix cross product
#' @description Calculates the cross product t(x) \%*\% y
#' @param x A numeric matrix
#' @param y An encoded matrix. Default, NULL, y = x.
#' @return t(x) \%*\% y
#' @export
crossProdrm <- function(x, y = NULL) {
    if (is.null(y)) return (crossprod(x))
    yd <- .decode.arg(y)
    if (is.list(yd)) yd <- do.call(rbind, yd)

    return (crossprod(x, yd))
}


#' @title Matrix cross product
#' @description Calculates the cross product t(x) \%*\% y
#' @param x A numeric matrix
#' @param y A list of numeric matrices. Default, NULL, y = x.
#' @param chunk Size of chunks into what the resulting matrix is partitioned. Default: 500.
#' @return \code{t(x) \%*\% y}
#' @export
crossProd <- function(x, y = NULL, chunk = 500) {
    nblocksrow <- ceiling(ncol(x)/chunk)
    sepblocksrow <- rep(ceiling(ncol(x)/nblocksrow), nblocksrow-1)
    sepblocksrow <- c(sepblocksrow, ncol(x) - sum(sepblocksrow))

    if (is.null(y)) {
        tcpblocks <- .partitionMatrix(crossprod(x), seprow=sepblocksrow)
        return (lapply(tcpblocks, function(tcpb) {
            return (lapply(tcpb, function(tcp) {
                .encode.arg(tcp)
            }))
        }))
    } else {
        nblockscol <- ceiling(ncol(y)/chunk)
        sepblockscol <- rep(ceiling(ncol(y)/nblockscol), nblockscol-1)
        sepblockscol <- c(sepblockscol, ncol(y) - sum(sepblockscol))
        tcpblocks <- .partitionMatrix(crossprod(x, y), seprow=sepblocksrow, sepcol=sepblockscol)
        return (lapply(tcpblocks, function(tcpb) {
            return (lapply(tcpb, function(tcp) {
                .encode.arg(tcp)
            }))
        }))
    }
}


#' @title Matrix cross product
#' @description Calculates the cross product x \%*\% t(y)
#' @param x A numeric matrix
#' @param y A list of numeric matrices. Default, NULL, y = x.
#' @param chunk Size of chunks into what the resulting matrix is partitioned. Default: 500.
#' @return \code{x \%*\% t(y)}
#' @export
tcrossProd <- function(x, y = NULL, chunk = 500) {
    if (is.null(y)) {
        nblocks <- ceiling(nrow(x)/chunk)
        sepblocks <- rep(ceiling(nrow(x)/nblocks), nblocks-1)
        sepblocks <- c(sepblocks, nrow(x) - sum(sepblocks))
        tcpblocks <- .partitionMatrix(tcrossprod(x), seprow=sepblocks)
        etcpblocks <- lapply(tcpblocks, function(tcpb) {
            return (lapply(tcpb, function(tcp) {
                .encode.arg(tcp)
            }))
        })
        return (etcpblocks)
    }
    return (lapply(y, function(yy) .encode.arg(matrix(tcrossprod(x, yy)))))
}


#' @title Bigmemory description of a matrix
#' @description Bigmemory description of a matrix
#' @param value Encoded value of a matrix
#' @import bigmemory
#' @return Bigmemory description of the given matrix
#' @export
matrix2Dsc <- function(value) {
    valued <- .decode.arg(value)
    tcp <- do.call(rbind, .decode.arg(valued))
    dscbigmatrix <- describe(as.big.matrix(tcp, backingfile = ""))
    rm(list=c("valued", "tcp"))
    return (dscbigmatrix)
}


#' @title Symmetric matrix reconstruction
#' @description Rebuild a matrix from its partition
#' @param matblocks List of lists of matrix blocks, obtained from .partitionMatrix
#' @param mc.cores Number of cores for parallel computing. Default: 1
#' @return The complete symmetric matrix
#' @keywords internal
.rebuildMatrix <- function(matblocks, mc.cores = 1) {
    uptcp <- lapply(matblocks, function(bl) do.call(cbind, bl))
    ## combine the blocks into one matrix
    if (length(uptcp)>1) {
        if (length(unique(sapply(uptcp, ncol)))==1) {
            tcp <- do.call(rbind, uptcp)
        } else {
            ## without the first layer of blocks
            no1tcp <- mclapply(2:length(uptcp), mc.cores=mc.cores, function(i) {
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
    print(class(tcp))
    print(dim(tcp))
    print(tcp[1:min(2, nrow(tcp)),1:min(2, nrow(tcp))])
    stopifnot(isSymmetric(tcp))
    rm(list=c("uptcp"))
    return (tcp)
}


#' @title Symmetric matrix reconstruction
#' @description Rebuild a matrix from its encoded partition
#' @param blocks List of list of encoded matrix blocks
#' @param mc.cores Number of cores for parallel computing. Default: 1
#' @return The complete symmetric matrix
#' @keywords internal
.rebuildMatrixEnc <- function(blocks, mc.cores = 1) {
    ## decode matrix blocks
    matblocks <- mclapply(blocks, mc.cores=mc.cores, function(y) {
        lapply(y, function(x) {
            return (do.call(rbind, .decode.arg(x)))
        })
    })
    tcp <- .rebuildMatrix(matblocks, mc.cores=mc.cores)
    return (tcp)
}


#' @title Symmetric matrix reconstruction
#' @description Rebuild a symmetric matrix from its partition through variables in the environment
#' @param symbol Generic variable name
#' @param len Length of list of lists of matrix blocks. Variables were generated as symbol11, symbol12, etc.
#' @param mc.cores Number of cores for parallel computing. Default: 1
#' @return The complete symmetric matrix
#' @export
rebuildMatrixVar <- function(symbol, len, mc.cores = 1) {
    ## access to matrix blocks
    matblocks <- mclapply(1:len, mc.cores=mc.cores, function(i) {
        lapply(1:(len-i+1), function(j) {
            dscblock <- get(paste(c(symbol, i, j), collapse="__"), envir = .GlobalEnv) #parent.frame())
            return (as.matrix(attach.big.matrix(dscblock)))
        })
    })
    tcp <- .rebuildMatrix(matblocks, mc.cores=mc.cores)
    return (tcp)
}


#' @title Push a symmetric matrix
#' @description Push symmetric matrix data into the federated server
#' @param value An encoded value to be pushed
#' @import bigmemory
#' @return Description of the pushed value
#' @export
pushSymmMatrixServer <- function(value) {
    valued <- .decode.arg(value)
    stopifnot(is.list(valued) && length(valued)>0)
    
    tcp <- .rebuildMatrixEnc(valued)
    dscbigmatrix <- describe(as.big.matrix(tcp))
    rm(list=c("valued", "tcp"))
    
    return (dscbigmatrix)
}


#' @title Matrix triple product
#' @description Calculate the triple product x \%*\% y \%*\% t(x)
#' @param x A numeric matrix
#' @param mate Encoded value of a vector of mate servers' names
#' @param chunk Size of chunks into what the resulting matrix is partitioned. Default: 500.
#' @param mc.cores Number of cores for parallel computing. Default: 1
#' @return \code{x \%*\% t(y)}
#' @export
tripleProdChunk <- function(x, mate, chunk = 500, mc.cores = 1) {
    pids <- .decode.arg(mate)
    
    nblocks <- ceiling(nrow(x)/chunk)
    sepblocks <- rep(ceiling(nrow(x)/nblocks), nblocks-1)
    sepblocks <- c(sepblocks, nrow(x) - sum(sepblocks))
    
    tps <- mclapply(pids, mc.cores=mc.cores, function(pid) {
        y <- get(paste("crossProdSelf", pid, sep='_'), envir = .GlobalEnv) #parent.frame())

        stopifnot(isSymmetric(y))
        # NB. this computation of tcpblocks could be done more efficiently with y being a chunked matrix in bigmemory
        tp <- tcrossprod(x, tcrossprod(x, y))
        tcpblocks <- .partitionMatrix(tp, seprow=sepblocks)
        etcpblocks <- lapply(tcpblocks, function(tcpb) {
            return (lapply(tcpb, function(tcp) {
                .encode.arg(tcp)
            }))
        })
        return (etcpblocks)
    })
    names(tps) <- pids
    return (tps)
}


#' @title Cross login
#' @description Call datashield.login on remote servers
#' @param logins An encoded dataframe with server, url, user, password, driver, and options fields.
#' @export
crossLogin <- function(logins) {
    require(DSOpal)
    loginfo <- .decode.arg(logins)
    myDf <- data.frame(server=loginfo$server,
                       url=loginfo$url,
                       user=loginfo$user,
                       password=loginfo$password,
                       driver=loginfo$driver,
                       options=loginfo$options)
    datashield.login(myDf)
}


#' @title Cross logout
#' @description Call datashield.logout on remote servers.
#' @param opals A list of opal objects
#' @export
crossLogout <- function(opals) {
    require(DSOpal)
    DSI::datashield.logout(opals)
}


#' @title Cross aggregate
#' @description Call datashield.aggregate on remote servers.
#' @param conns A list of DSConnection-class.
#' @param expr An encoded expression to evaluate.
#' @param async See DSI::datashield.aggregate options. Default: TRUE.
#' @import DSI
#' @export
crossAggregatePrimal <- function(conns, expr, async = T) {
    expr <- .decode.arg(expr)
    if (grepl("^singularProd\\(", expr)) {
        DSI::datashield.aggregate(conns=conns, expr=as.symbol(expr), async=async)
    } else {
        stop(paste0("Failed to execute: ", expr))
    }
}


#' @title Cross aggregate through two-layer connection
#' @description Call datashield.aggregate on remote servers through two-layer connection.
#' @param conns A list of DSConnection-class.
#' @param expr An encoded expression to evaluate. This is restricted to pushToDscDual.
#' @param async See DSI::datashield.aggregate options. Default: TRUE.
#' @import DSI
#' @export
crossAggregateDual <- function(conns, expr, async = T) {
    expr <- .decode.arg(expr)
    if (grepl("^pushToDscDual\\(", expr)) {
        DSI::datashield.aggregate(conns=conns, expr=as.symbol(expr), async=async)
    } else {
        stop(paste0("Failed to execute: ", expr))
    }
}


#' @title Description of a pushed value
#' @description Description of a pushed value
#' @param conns A list of DSConnection-class.
#' @param expr An encoded expression to evaluate.
#' @param async See DSI::datashield.aggregate options. Default: TRUE.
#' @return Returned value of given expression on opal
#' @import DSI
#' @export
dscPush <- function(conns, expr, async = T) {
  expr <- .decode.arg(expr)
  stopifnot(grepl("^as.call", expr))
  expr <- eval(str2expression(expr))
  return (DSI::datashield.aggregate(conns=conns, expr=expr, async=async))
}


#' @title Bigmemory description of a pushed object
#' @description Bigmemory description of a pushed object
#' @param conns A one-element list of DSConnection-class.
#' @param symbol Name of an object to be pushed.
#' @param async See DSI::datashield.aggregate options. Default: TRUE.
#' @return Bigmemory description of the pushed object on conns
#' @import DSI
#' @export
pushToDscDual <- function(conns, symbol, async = T) {
    ## TODO: check for allowed conns
    stopifnot(is.list(conns) && length(conns)==1 && class(conns[[1]])=="OpalConnection")
    
    chunkList <- get(symbol, envir = .GlobalEnv) #parent.frame())
    stopifnot(is.list(chunkList) && is.list(chunkList[[1]]))
    if (is.list(chunkList[[1]][[1]])) {
        dsc <- lapply(chunkList, function(z) {
            return (lapply(z, function(x) {
                return (lapply(x, function(y) {
                    expr <- list(as.symbol("matrix2Dsc"), y)
                    y.dsc <- DSI::datashield.aggregate(conns=conns, expr=as.call(expr), async=async)
                    return (y.dsc[[1]])
                }))
            }))
        })
    } else {
        dsc <- lapply(chunkList, function(x) {
            return (lapply(x, function(y) {
                expr <- list(as.symbol("matrix2Dsc"), y)
                y.dsc <- DSI::datashield.aggregate(conns=conns, expr=as.call(expr), async=async)
                return (y.dsc[[1]])
            }))
        })
    }
    return (dsc)
}


#' @title Bigmemory description of a pushed object
#' @description Bigmemory description of a pushed object
#' @param conns A list of DSConnection-classes.
#' @param symbol Name of an object to be pushed.
#' @param sourcename Name of the pushed object source.
#' @param async See DSI::datashield.aggregate options. Default: TRUE.
#' @return Bigmemory description of the pushed object on conns
#' @import DSI
#' @export
pushToDscPrimal <- function(conns, symbol, sourcename, async = T) {
    ## TODO: check for allowed conns
    stopifnot(is.list(conns) && length(setdiff(unique(sapply(conns, class)), "OpalConnection"))==0)
    
    chunkList <- get(symbol, envir = .GlobalEnv) #parent.frame())
    invisible(lapply(1:length(chunkList), function(i) {
        lapply(1:length(chunkList[[i]]), function(j) {
            DSI::datashield.assign(conns, paste(c(sourcename, i, j), collapse="__"), 
                                   as.call(list(as.symbol("matrix2Dsc"), 
                                                chunkList[[i]][[j]])), 
                                   async=async)
        })
    }))
    DSI::datashield.assign(conns, paste(symbol, sourcename, sep="_"), 
                           as.call(list(as.symbol("rebuildMatrixVar"),
                                        sourcename,
                                        length(chunkList))),
                           async=async)
}


#' @title Cross assign
#' @description Call datashield.assign on remote servers.
#' @param conns A list of DSConnection-class.
#' @param symbol Name of an R symbol.
#' @param value A variable name or an R expression with allowed assign function calls.
#' @param value.call A logical value, TRUE if value is function call, FALSE if value is a variable name.
#' @param async See DSI::datashield.assign options. Default: TRUE.
#' @import DSI
#' @export
crossAssign <- function(conns, symbol, value, value.call, variables = NULL, async = T) {
    valued <- .decode.arg(value)
    variablesd <- .decode.arg(variables)
    DSI::datashield.assign(conns=conns, symbol=symbol, 
                           value=ifelse(value.call, as.symbol(valued), valued), 
                           variables=variablesd, async=async)
}


#' @title Cross assign a preprocessing function
#' @description Call preprocessing function on remote servers.
#' @param conns A list of DSConnection-class.
#' @param func Encoded definition of a function for preparation of raw data matrices. 
#' Two arguments are required: conns (list of DSConnection-classes), 
#' symbol (name of the R symbol) (see datashield.assign).
#' @param symbol Encoded vector of names of R symbols to assign in the Datashield R session on each server in \code{logins}.
#' The assigned R variables will be used as the input raw data to compute covariance matrix.
#' Other assigned R variables in \code{func} are ignored.
#' @export
crossAssignFunc <- function(conns, func, symbol) {
    funcPreProc <- .decode.arg(func)
    querytables <- .decode.arg(symbol)
    tryCatch({
        ## take a snapshot of the current session
        safe.objs <- .ls.all()
        safe.objs[['.GlobalEnv']] <- setdiff(safe.objs[['.GlobalEnv']], '.Random.seed')  # leave alone .Random.seed for sample()
        ## lock everything so no objects can be changed
        .lock.unlock(safe.objs, lockBinding)
        
        ## apply funcPreProc for preparation of querytables on conns
        ## TODO: control hacking!
        ## TODO: control identical colnames!
        funcPreProc(conns=conns, symbol=querytables)
        
        ## unlock back everything
        .lock.unlock(safe.objs, unlockBinding)
        ## get rid of any sneaky objects that might have been created in the filters as side effects
        .cleanup(safe.objs)
    }, error=function(e) {
        print(paste0("DATA MAKING CROSS PROCESS: ", e))
        DSI::datashield.logout(conns)
    })
    return (NULL)
}


#' @title Sum matrices
#' @description Compute the sum of a matrix and those stored in bigmemory
#' @param dsc A list of big memory descriptions
#' @return Sum of x and those stored in dsc
#' @import bigmemory
#' @keywords internal
.sumMatrices.rm <- function(dsc = NULL) {
    dscmat <- lapply(dsc, function(dscblocks) {
        y <- as.matrix(attach.big.matrix(dscblocks))
        return (y)
    })
    return (Reduce("+", dscmat))
}


#' @title Sum matrices
#' @description Compute the sum of crossProd matrices from mates
#' @param mate Encoded value of a vector of mate servers' names
#' @param mc.cores Number of cores for parallel computing. Default: 1
#' @return A matrix
#' @import parallel
#' @keywords internal
.sumMatrices <- function(mate, mc.cores=1) {
    pids <- .decode.arg(mate)
    covmat <- mclapply(pids, mc.cores=mc.cores, function(pid) {
        y <- get(paste("crossProdSelf", pid, sep='_'), envir = .GlobalEnv)
        return (y)
    })
    return (Reduce("+", covmat))
}


#' @title Federated covariance matrix
#' @description Compute the covariance matrix for the virtual cohort
#' @param loginFD Login information of the FD server (one of the servers containing cohort data)
#' @param logins Login information of other servers containing cohort data
#' @param funcPreProc Definition of a function for preparation of raw data matrices. 
#' Two arguments are required: conns (list of DSConnection-classes), 
#' symbol (name of the R symbol) (see datashield.assign).
#' @param querytables Name (or a vector of two names) of the R symbol(s) to assign in the Datashield R session on each server in \code{logins}.
#' The assigned R variable(s) will be used as the input raw data to compute covariance matrix.
#' Other assigned R variables in \code{func} are ignored.
#' @param querysubset A list of index vectors indicating the subsets of individuals to consider. 
#' Default, NULL, all individuals are considered.
#' @param covSpace The space of variables where covariance matrix is computed. If \code{length(querytables)=1},
#' \code{covSpace} is always \code{"X"}. If \code{length(querytables)=2}, it can be \code{"X"} for the first querytable,
#' \code{"Y"} for the second querytable, and \code{"XY"} for covariance between the two querytables. Default, \code{"X"}.
#' @param chunk Size of chunks into what the resulting matrix is partitioned. Default: 500.
#' @return Covariance matrix of the virtual cohort
#' @import DSI parallel bigmemory
#' @keywords internal
.federateCov <- function(loginFD, logins, funcPreProc, querytables, querysubset = NULL, covSpace = "X", chunk = 500, mc.cores = 1) {
    require(DSOpal)
    ## covariance of only one matrix or between two matrices
    stopifnot(length(querytables) %in% c(1,2))
    covSpace <- match.arg(covSpace, choices=c('X', 'Y', "XY"))
    loginFDdata <- .decode.arg(loginFD)
    logindata   <- .decode.arg(logins)
    
    ## assign crossprod matrix on each individual server
    opals <- datashield.login(logins=logindata)
    
    tryCatch({
        ## take a snapshot of the current session
        safe.objs <- .ls.all()
        safe.objs[['.GlobalEnv']] <- setdiff(safe.objs[['.GlobalEnv']], '.Random.seed')  # leave alone .Random.seed for sample()
        ## lock everything so no objects can be changed
        .lock.unlock(safe.objs, lockBinding)
        
        ## apply funcPreProc for preparation of querytables on opals
        ## TODO: control hacking!
        ## TODO: control identical colnames!
        funcPreProc(conns=opals, symbol=querytables)
        ## unlock back everything
        .lock.unlock(safe.objs, unlockBinding)
        ## get rid of any sneaky objects that might have been created in the filters as side effects
        .cleanup(safe.objs)
    }, error=function(e) {
        print(paste0("DATA MAKING PROCESS: ", e))
        return (paste0("DATA MAKING PROCESS: ", e, ' --- ', datashield.symbols(opals), ' --- ', datashield.errors(), ' --- ', datashield.logout(opals)))
    })
    tryCatch({
        if (is.null(querysubset)) {
            datashield.assign(opals, "centeredData", as.symbol(paste0('center(', querytables[1], ')')), async=F)
        } else {
            stopifnot(all(names(opals)==names(querysubset)))
            lapply(names(opals), function(opn) {
                datashield.assign(opals[opn], "centeredData", as.symbol(paste0("center(", querytables[1], ", subset='", .encode.arg(querysubset[[opn]]), "')")), async=F)
            })
        }
        size <- sapply(datashield.aggregate(opals, as.symbol('dsDim(centeredData)'), async=F), function(x) x[1])
        if (length(querytables)==1 || covSpace=="X") {
            datashield.assign(opals, "crossProdSelf", as.symbol(paste0('crossProd(x=centeredData, y=NULL, chunk=', chunk, ')')), async=T)
        } else {
            if (is.null(querysubset)) {
                datashield.assign(opals, "centeredData2", as.symbol(paste0('center(', querytables[2], ')')), async=T)
            } else {
                stopifnot(all(names(opals)==names(querysubset)))
                lapply(names(opals), function(opn) {
                    datashield.assign(opals[opn], "centeredData2", as.symbol(paste0("center(", querytables[2], ", subset='", .encode.arg(querysubset[[opn]]), "')")), async=T)
                })
            }
            size2 <- sapply(datashield.aggregate(opals, as.symbol('dsDim(centeredData2)'), async=T), function(x) x[1])
            stopifnot(all(size==size2))
            if (covSpace=="Y") {
                datashield.assign(opals, "crossProdSelf", as.symbol(paste0('crossProd(x=centeredData2, y=NULL, chunk=', chunk, ')')), async=T)
            } else {
                datashield.assign(opals, "crossProdSelf", as.symbol(paste0('crossProd(x=centeredData, y=centeredData2, chunk=', chunk, ')')), async=T)
            }
        }
        
        ## push data from non-FD servers to FD-assigned server: user and password for login between servers are required
        loginFDdata$user     <- loginFDdata$userserver
        loginFDdata$password <- loginFDdata$passwordserver
        datashield.assign(opals, 'FD', as.symbol(paste0("crossLogin('", .encode.arg(loginFDdata), "')")), async=T)
        ## send X'X from opn to FD
        tryCatch({
            lapply(names(opals), function(opn) {
                command.opn <- list(as.symbol("pushToDscPrimal"),
                                    as.symbol("FD"),
                                    'crossProdSelf',
                                    opn,
                                    async=T)
                cat("Command: pushToDscPrimal(FD, crossProdSelf...", "\n")
                invisible(datashield.aggregate(opals[opn], as.call(command.opn), async=F))
            })
            .printTime(paste0(".federateSSCP pairwise X'X communicated: ", opn))
            
            # command <- paste0("dscPush(FD, '", 
            #                   .encode.arg(paste0("as.call(list(as.symbol('pushSymmMatrixServer'), dsMOprimal:::.encode.arg(crossProdSelf)", "))")), 
            #                   "', async=T)")
            # cat("Command: ", command, "\n")
            
            # crossProdSelfDSC <- DSI::datashield.aggregate(opals, as.symbol(command), async=T)
            # print(covSpace)
            # crossProdSelfDSC <- lapply(crossProdSelfDSC, function(dscblocks) {
            #     print(dscblocks[[1]])
            #     return (dscblocks[[1]])
            # })
            # rescov <- .sumMatrices(crossProdSelfDSC)/(sum(size)-1)

            rescov <- .sumMatrices(.encode.arg(names(opals)), mc.cores=mc.cores)/(sum(size)-1)
            ## set dimnames to covariance matrix
            if (covSpace=="X" || covSpace=="XY") {
                rownames(rescov) <- DSI::datashield.aggregate(opals[1], as.symbol('colNames(centeredData)'), async=T)[[1]]
            } else {
                rownames(rescov) <- DSI::datashield.aggregate(opals[1], as.symbol('colNames(centeredData2)'), async=T)[[1]]
            }
            if (covSpace=="Y" || covSpace=="XY") {
                colnames(rescov) <- DSI::datashield.aggregate(opals[1], as.symbol('colNames(centeredData2)'), async=T)[[1]]
            } else {
                colnames(rescov) <- DSI::datashield.aggregate(opals[1], as.symbol('colNames(centeredData)'), async=T)[[1]]
            }
        }, 
        error=function(e) {
            print(paste0("COVARIATES PUSH PROCESS: ", e));
            return(paste0("COVARIATES PUSH PROCESS: ", e))
        }, 
        finally=datashield.assign(opals, 'crossEnd', as.symbol("crossLogout(FD)"), async=T))
        
    }, 
    error=function(e) {
        print(paste0("COVARIATES PROCESS: ", e))
        return(paste0("COVARIATES PROCESS: ", e, ' --- ', datashield.symbols(opals), ' --- ', datashield.errors()))
    }, 
    finally=datashield.logout(opals))
    gc(reset=F)
    return (rescov)
}


#' @title Federated PCA
#' @description Perform the principal component analysis for the virtual cohort
#' @param loginFD Login information of the FD server (one of the servers containing cohort data)
#' @param logins Login information of other servers containing cohort data
#' @param func Encoded definition of a function for preparation of raw data matrices. 
#' Two arguments are required: conns (list of DSConnection-classes), 
#' symbol (name of the R symbol) (see datashield.assign).
#' @param symbol Encoded name of the R symbol to assign in the Datashield R session on each server in \code{logins}.
#' The assigned R variable will be used as the input raw data to compute covariance matrix for PCA.
#' Other assigned R variables in \code{func} are ignored.
#' @param ncomp Number of components. Default: 2.
#' @param chunk Size of chunks into what the resulting matrix is partitioned. Default: 500.
#' @return PCA object
#' @import DSI parallel bigmemory
#' @examples
#' dataProc <- function(conns, symbol) {
#'     DSI::datashield.assign(conns, symbol, 'test.CNSIM', variables=c('LAB_TSC', 'LAB_TRIG', 'LAB_HDL', 'LAB_GLUC_ADJUSTED', 'PM_BMI_CONTINUOUS'), async=T)
#' }
#' dataProc(conns=opals, symbol="rawData")
#' federatePCA(.encode.arg(loginFD), .encode.arg(logins), .encode.arg(dataProc, serialize.it = T), .encode.arg("rawData"))
#' @export
federatePCA <- function(loginFD, logins, func, symbol, ncomp = 2, chunk = 500) {
    funcPreProc <- .decode.arg(func)
    querytables <- .decode.arg(symbol)
    if (length(querytables) != 1) {
        stop("One data matrix is required!")
    }
    covmat <- .federateCov(loginFD, logins, funcPreProc, querytables, chunk=chunk)
    pcaObj <- princomp(covmat=covmat)

    if (ncomp > 2 && length(which(pcaObj$sdev > 1e-6)) <= ncomp) {
        print(paste0("Security issue: maximum ", length(which(pcaObj$sdev > 1e-6)) - 1, " components could be inquired. ncomp will be set to 2."))
        ncomp <- 2
    }
    if (ncomp < 1) {
        print("ncomp should be at least 1. ncomp will be set to 2.")
        ncomp <- 2
    }
    pcaObj$loadings <- pcaObj$loadings[, 1:ncomp, drop=F]
    
    ## compute loadings
    logindata <- .decode.arg(logins)
    opals <- datashield.login(logins=logindata)
    
    out <- tryCatch({
        ## take a snapshot of the current session
        safe.objs <- .ls.all()
        safe.objs[['.GlobalEnv']] <- setdiff(safe.objs[['.GlobalEnv']], '.Random.seed')  # leave alone .Random.seed for sample()
        ## lock everything so no objects can be changed
        .lock.unlock(safe.objs, lockBinding)
        
        ## apply funcPreProc for preparation of querytables on opals
        ## TODO: control hacking!
        ## TODO: control identical colnames!
        funcPreProc(conns=opals, symbol=querytables)
        ## unlock back everything
        .lock.unlock(safe.objs, unlockBinding)
        ## get rid of any sneaky objects that might have been created in the filters as side effects
        .cleanup(safe.objs)
        
        datashield.assign(opals, "centeredData", as.symbol(paste0('center(', querytables[1], ')')), async=T)
        expr <- list(as.symbol("loadings"),
                     as.symbol("centeredData"),
                     .encode.arg(pcaObj$loadings), #TODO: this can be bugged due to large loadings matrix
                     "prod")
        pcaObj$scores <- do.call(rbind, datashield.aggregate(opals,
                                                             as.call(expr),
                                                             async=T))
    },
    error=function(e) {
        print(paste0("LOADINGS MAKING PROCESS: ", e))
        return (paste0("LOADINGS MAKING PROCESS: ", e, ' --- ', datashield.symbols(opals), ' --- ', datashield.errors(), ' --- ', datashield.logout(opals)))
    },
    finally=datashield.logout(opals))
    
    return (pcaObj)
}


#' @title RCCA tuning
#' @description Estimate optimized parameters of regulation lambda1 and lambda2
#' @param loginFD Login information of the FD server (one of the servers containing cohort data).
#' @param logins Login information of servers containing cohort data.
#' @param func Encoded definition of a function for preparation of raw data matrices. 
#' Two arguments are required: conns (list of DSConnection-classes), 
#' symbol (names of the two R symbols) (see datashield.assign).
#' @param symbol Encoded vector of names of the two R symbols to assign in the Datashield R session on each server in \code{logins}.
#' The two assigned R variables will be used as the input raw data to compute covariance matrices for CCA.
#' Other assigned R variables in \code{func} are ignored.
#' @param lambda1 Regularized parameter value for first data set. Default, 0.
#' @param lambda2 Regularized parameter value for second data set. Default, 0.
#' @param tune Logical value indicating whether the tuning for lambda values will be performed. Default, FALSE, no tuning.
#' @param tune_param Tuning parameters. \code{nfold} n-fold cross-validation. 
#' @param grid1 Checking values for \code{lambda1}.
#' @param grid2 Checking values for \code{lambda2}.
#' @return Optimal values of \code{lambda1} and \code{lambda2}.
#' @keywords internal
.estimateR <- function(loginFD, logins, funcPreProc, querytables,
                       nfold = 5, grid1 = seq(0.001, 1, length = 5), grid2 = seq(0.001, 1, length = 5)) {
    stopifnot(length(querytables) == 2)
    loginFDdata <- .decode.arg(loginFD)
    logindata   <- .decode.arg(logins)
    
    opals <- datashield.login(logins=.decode.arg(logins))
    nNode <- length(opals)
    
    tryCatch({
        ## take a snapshot of the current session
        safe.objs <- .ls.all()
        safe.objs[['.GlobalEnv']] <- setdiff(safe.objs[['.GlobalEnv']], '.Random.seed')  # leave alone .Random.seed for sample()
        ## lock everything so no objects can be changed
        .lock.unlock(safe.objs, lockBinding)
        
        ## apply funcPreProc for preparation of querytables on opals
        ## TODO: control hacking!
        ## TODO: control identical colnames!
        funcPreProc(conns=opals, symbol=querytables)
        
        ## unlock back everything
        .lock.unlock(safe.objs, unlockBinding)
        ## get rid of any sneaky objects that might have been created in the filters as side effects
        .cleanup(safe.objs)
    }, error=function(e) {
        print(paste0("DATA MAKING PROCESS: ", e))
        datashield.logout(opals)
    })
    
    tryCatch({
        DSI::datashield.assign(opals, "centeredDatax", as.symbol(paste0('center(', querytables[1], ')')), async=T)
        DSI::datashield.assign(opals, "centeredDatay", as.symbol(paste0('center(', querytables[2], ')')), async=T)
        sizex <- sapply(DSI::datashield.aggregate(opals, as.symbol('dsDim(centeredDatax)'), async=T), function(x) x[1])
        sizey <- sapply(DSI::datashield.aggregate(opals, as.symbol('dsDim(centeredDatay)'), async=T), function(x) x[1])
        stopifnot(all(sizex==sizey))
        
        ## random Mfold partitions to leave out
        names(sizex) <- names(opals)
        foldspar <- split(unlist(lapply(names(opals), function(opn) {
            paste(opn, 1:sizex[opn], sep="_")
        }))[sample(1:sum(sizex))], rep(1:nfold, length = sum(sizex)))
        foldslef <- lapply(foldspar, function(fl) {
            setNames(mclapply(names(opals), mc.cores=min(nNode, detectCores()), function(opn) {
                sort(as.numeric(sub(fl[grep(opn, fl)], pattern=paste0(opn,"_"), replacement='')))
            }), names(opals))
        })
        ## remaining individuals on each cohort
        foldsrem <- lapply(foldslef, function(fl) {
            setNames(mclapply(names(opals), mc.cores=min(nNode, detectCores()), function(opn) {
                setdiff(1:sizex[opn], fl[[opn]])
            }), names(opals))
        })
        grid <- expand.grid(grid1, grid2)
        cv.score <- apply(grid, 1, function(lambda) {
            xscore <- NULL
            yscore <- NULL
            for (m in 1:nfold) {
                ## covariance matrices for the virtual cohort
                Cxx <- .federateCov(loginFD, logins, funcPreProc, querytables, querysubset=foldsrem[[m]], covSpace="X")
                Cyy <- .federateCov(loginFD, logins, funcPreProc, querytables, querysubset=foldsrem[[m]], covSpace="Y")
                Cxy <- .federateCov(loginFD, logins, funcPreProc, querytables, querysubset=foldsrem[[m]], covSpace="XY")
                
                ## add parameters of regularization
                Cxx <- Cxx + diag(lambda[1], ncol(Cxx))
                Cyy <- Cyy + diag(lambda[2], ncol(Cyy))
                ## CCA core call
                res <- fda::geigen(Cxy, Cxx, Cyy)
                names(res) <- c("cor", "xcoef", "ycoef")
                rownames(res$xcoef) <- rownames(Cxx)
                rownames(res$ycoef) <- rownames(Cyy)
                ## tuning scores
                lapply(names(opals), function(opn) {
                    DSI::datashield.assign(opals[opn], "centeredDataxm", as.symbol(paste0("center(", querytables[1], ", subset='", .encode.arg(foldslef[[m]][[opn]]), "')")), async=T)
                    DSI::datashield.assign(opals[opn], "centeredDataym", as.symbol(paste0("center(", querytables[2], ", subset='", .encode.arg(foldslef[[m]][[opn]]), "')")), async=T)
                })
                cvx <- do.call(rbind, DSI::datashield.aggregate(opals, as.call(list(as.symbol("loadings"),
                                                                                    as.symbol("centeredDataxm"),
                                                                                    .encode.arg(res$xcoef[,1,drop=F]),
                                                                                    "prod")), async=T))
                cvy <- do.call(rbind, DSI::datashield.aggregate(opals, as.call(list(as.symbol("loadings"),
                                                                                    as.symbol("centeredDataym"),
                                                                                    .encode.arg(res$ycoef[,1,drop=F]),
                                                                                    "prod")), async=T))
                xscore <- c(xscore, cvx)
                yscore <- c(yscore, cvy)
            }
            return (cor(xscore, yscore, use = "pairwise"))
        })
    }, finally=DSI::datashield.logout(opals))
    cv.score.grid <- cbind(grid, cv.score)
    mat <- matrix(cv.score, nrow=length(grid1), ncol=length(grid2))
    plot <- FALSE
    if (isTRUE(plot)) image(list(grid1 = grid1, grid2 = grid2, mat = mat))
    opt <- cv.score.grid[which.max(cv.score.grid[,3]), ]
    out <- list(opt.lambda1 = opt[[1]], 
                opt.lambda2 = opt[[2]], 
                opt.score   = opt[[3]],
                grid1       = grid1, 
                grid2       = grid2, 
                mat         = mat)
    out$call <- match.call()
    class(out) <- "estimateR"
    return (out)
}


#' @title Federated RCCA
#' @description Perform the regularized canonical correlation analysis for the virtual cohort
#' @param loginFD Login information of the FD server (one of the servers containing cohort data).
#' @param logins Login information of servers containing cohort data.
#' @param func Encoded definition of a function for preparation of raw data matrices. 
#' Two arguments are required: conns (list of DSConnection-classes), 
#' symbol (names of the two R symbols) (see datashield.assign).
#' @param symbol Encoded vector of names of the two R symbols to assign in the Datashield R session on each server in \code{logins}.
#' The two assigned R variables will be used as the input raw data to compute covariance matrices for CCA.
#' Other assigned R variables in \code{func} are ignored.
#' @param lambda1 Non-negative regularized parameter value for first data set. Default, 0. If there are more variables than samples, it should be > 0.
#' @param lambda2 Non-negative regularized parameter value for second data set. Default, 0. If there are more variables than samples, it should be > 0.
#' @param chunk Size of chunks into what the SSCP matrix is partitioned. Default: 500.
#' @param tune Logical value indicating whether the tuning for lambda values will be performed. Default, FALSE, no tuning.
#' @param tune_param Tuning parameters. \code{nfold} n-fold cross-validation. \code{grid1} checking values for \code{lambda1}.
#' \code{grid2} checking values for \code{lambda2}.
#' @return RCCA object
#' @import DSI parallel bigmemory
#' @importFrom fda geigen
#' @examples
#' dataProc <- function(conns, symbol) {
#'     DSI::datashield.assign(conns, symbol[1], 'test.CNSIM', variables=c('LAB_TSC', 'LAB_TRIG', 'LAB_HDL'), async=T)
#'     DSI::datashield.assign(conns, symbol[2], 'test.CNSIM', variables=c('LAB_GLUC_ADJUSTED', 'PM_BMI_CONTINUOUS'), async=T)
#' }
#' dataProc(conns=opals, symbol=c("rawDataX", "rawDataY"))
#' federateRCCA(.encode.arg(loginFD), .encode.arg(logins), .encode.arg(dataProc, serialize.it = T), .encode.arg(c("rawDataX", "rawDataY")))
#' @export
federateRCCA <- function(loginFD, logins, func, symbol, lambda1 = 0, lambda2 = 0, chunk = 500,
                         tune = FALSE, tune_param = .encode.arg(list(nfold = 5, grid1 = seq(0.001, 1, length = 5), grid2 = seq(0.001, 1, length = 5)))) {
    require(DSOpal)
    funcPreProc <- .decode.arg(func)
    querytables <- .decode.arg(symbol)
    if (length(querytables) != 2) {
        stop("Two data matrices are required!")
    }

    ## estimating the parameters of regularization
    if (isTRUE(tune)) {
        tune_param <- .decode.arg(tune_param)
        tuneres <- .estimateR(loginFD, logins, funcPreProc, querytables,
                              nfold=tune_param$nfold, grid1=tune_param$grid1, grid2=tune_param$grid2)
        lambda1 <- tuneres$opt.lambda1
        lambda2 <- tuneres$opt.lambda2
    }
    ## covariance matrices for the virtual cohort
    Cxx <- .federateCov(loginFD, logins, funcPreProc, querytables, covSpace="X", chunk = chunk)
    Cyy <- .federateCov(loginFD, logins, funcPreProc, querytables, covSpace="Y", chunk = chunk)
    Cxy <- .federateCov(loginFD, logins, funcPreProc, querytables, covSpace="XY", chunk = chunk)

    ## add parameters of regularization
    Cxx <- Cxx + diag(lambda1, ncol(Cxx))
    Cyy <- Cyy + diag(lambda2, ncol(Cyy))
    
    ## CCA core call
    res <- fda::geigen(Cxy, Cxx, Cyy)
    names(res) <- c("cor", "xcoef", "ycoef")
    rownames(res$xcoef) <- rownames(Cxx)
    rownames(res$ycoef) <- rownames(Cyy)
    res$lambda <- list(lambda1=lambda1, lambda2=lambda2)

    ## assign centered data on each individual server
    ## NB: this block only works with some call a priori, e.g. .federateCov, or with require(DSOpal) !!!
    opals <- datashield.login(logins=.decode.arg(logins))
    nNode <- length(opals)
    
    tryCatch({
        ## take a snapshot of the current session
        safe.objs <- .ls.all()
        safe.objs[['.GlobalEnv']] <- setdiff(safe.objs[['.GlobalEnv']], '.Random.seed')  # leave alone .Random.seed for sample()
        ## lock everything so no objects can be changed
        .lock.unlock(safe.objs, lockBinding)
        
        ## apply funcPreProc for preparation of querytables on opals
        ## TODO: control hacking!
        ## TODO: control identical colnames!
        funcPreProc(conns=opals, symbol=querytables)
        
        ## unlock back everything
        .lock.unlock(safe.objs, unlockBinding)
        ## get rid of any sneaky objects that might have been created in the filters as side effects
        .cleanup(safe.objs)
    }, error=function(e) {
        print(paste0("DATA MAKING PROCESS: ", e))
        datashield.logout(opals)
    })
    
    tryCatch({
        DSI::datashield.assign(opals, "centeredDatax", as.symbol(paste0('center(', querytables[1], ')')), async=T)
        DSI::datashield.assign(opals, "centeredDatay", as.symbol(paste0('center(', querytables[2], ')')), async=T)
        sizex <- sapply(DSI::datashield.aggregate(opals, as.symbol('dsDim(centeredDatax)'), async=T), function(x) x[1])
        sizey <- sapply(DSI::datashield.aggregate(opals, as.symbol('dsDim(centeredDatay)'), async=T), function(x) x[1])
        stopifnot(all(sizex==sizey))
        sampleNames <- DSI::datashield.aggregate(opals, as.symbol('rowNames(centeredDatax)'), async=T)
        sampleNames <- unlist(lapply(names(sampleNames), function(x) paste0(x, "_", sampleNames[[x]])), use.names=F)
        res$names <- list(Xnames=rownames(Cxx),
                          Ynames=rownames(Cyy),
                          ind.names=sampleNames)
        ## canonical covariates
        cvx <- do.call(rbind, DSI::datashield.aggregate(opals, as.call(list(as.symbol("loadings"),
                                                                       as.symbol("centeredDatax"),
                                                                       .encode.arg(res$xcoef),
                                                                       "prod")), async=T))
        cvy <- do.call(rbind, DSI::datashield.aggregate(opals, as.call(list(as.symbol("loadings"),
                                                                       as.symbol("centeredDatay"),
                                                                       .encode.arg(res$ycoef),
                                                                       "prod")), async=T))
    }, error=function(e) print(paste0("COVARIATES PROCESS: ", e)), finally=DSI::datashield.logout(opals))
    
    ## loadings: correlation between raw data and canonical covariates
    ## formula: cor(a,b) = diag(1/sqrt(diag(cov(a)))) %*% cov(a,b) %*% diag(1/sqrt(diag(cov(b))))
    invdiagcovx <- diag(1/sqrt(diag(Cxx)))
    invdiagcovy <- diag(1/sqrt(diag(Cyy)))
    invdiagcovcvx <- diag(1/sqrt(diag(cov(cvx))))
    invdiagcovcvy <- diag(1/sqrt(diag(cov(cvy))))
    # xxscores <- Reduce('+', lapply(names(opals), function(opn) {
    #     xx <- datashield.aggregate(opals[opn], as.call(list(as.symbol("loadings"),
    #                                                         as.symbol("centeredDatax"),
    #                                                         .encode.arg(cvx[[opn]]), ## cvx obtained without rbind
    #                                                         "crossprod")), 
    #                                async=T)
    #     return (xx[[1]])
    # }))
    xxscores <- invdiagcovx %*% Cxx %*% res$xcoef %*% invdiagcovcvx
    # yxscores <- Reduce('+', lapply(names(opals), function(opn) {
    #     yx <- datashield.aggregate(opals[opn], as.call(list(as.symbol("loadings"),
    #                                                         as.symbol("centeredDatay"),
    #                                                         .encode.arg(cvx[[opn]]),
    #                                                         "crossprod")),
    #                                async=T)
    #     return (yx[[1]])
    # }))
    yxscores <- invdiagcovy %*% t(Cxy) %*% res$xcoef %*% invdiagcovcvx
    # xyscores <- Reduce('+', lapply(names(opals), function(opn) {
    #     xy <- datashield.aggregate(opals[opn], as.call(list(as.symbol("loadings"),
    #                                                         as.symbol("centeredDatax"),
    #                                                         .encode.arg(cvy[[opn]]), ## cvy obtained without rbind
    #                                                         "crossprod")),
    #                                async=T)
    #     return (xy[[1]])
    # }))
    xyscores <- invdiagcovx %*% Cxy %*% res$ycoef %*% invdiagcovcvy
    # yyscores <- Reduce('+', lapply(names(opals), function(opn) {
    #     yy <- datashield.aggregate(opals[opn], as.call(list(as.symbol("loadings"),
    #                                                         as.symbol("centeredDatay"),
    #                                                         .encode.arg(cvy[[opn]]),
    #                                                         "crossprod")),
    #                                async=T)
    #     return (yy[[1]])
    # }))
    yyscores <- invdiagcovy %*% Cyy %*% res$ycoef %*% invdiagcovcvy
    
    res$scores <- list(xscores=cvx,
                       yscores=cvy,
                       corr.X.xscores=xxscores,
                       corr.Y.xscores=yxscores,
                       corr.X.yscores=xyscores,
                       corr.Y.yscores=yyscores)
    
    return (res)
}
