#' @title Matrix dimension
#' @description Dimension of a data frame or matrix
#' @param x A matrix, a data frame, or a list of matrices or data frames
#' @returns Dimension of x
#' @export
dsDim <- function(x) {
    if (is.list(x) && !is.data.frame(x)) {
        return (lapply(x, dim))
    }
    return (dim(x))
}


#' @title Range of a variable in a data frame
#' @description Return random approximated range of a numeric variable in a
#' data frame
#' @param x A variable of a data frame in form of 'data.frame$variable'.
#' @returns Approximate range of x
#' @importFrom stats runif
#' @export
dsRange <- function(x) {
    if (!grepl("$", deparse(substitute(x))))
        stop("x should be a variable of a data frame in form of
             'data.frame$variable'.")
    rx <- range(x)
    drx <- diff(rx)
    range.min <- rx[1] - drx*runif(1, 1e-2, 1)/100
    range.max <- rx[2] + drx*runif(1, 1e-2, 1)/100
    return (c(range.min, range.max))
}


#' @title Factor coercion
#' @description Factor coercion of a categorical variable in a data frame
#' @param x A variable of a data frame in form of 'data.frame$variable'.
#' @returns An object of class \code{factor}
#' @export
dsFactor <- function(x) {
    if (!grepl("$", deparse(substitute(x))))
        stop("x should be a variable of a data frame in form
             of 'data.frame$variable'.")
    return (as.factor(x))
}


#' @title Levels attributes
#' @description Return levels of a categorical variable in a data frame
#' @param x A variable of a data frame in form of 'data.frame$variable'
#' @returns Levels of x
#' @export
dsLevels <- function(x) {
    if (!is.factor(x)) stop("A factor is required.")
    return (levels(x))
}


#' @title Row names
#' @description Assign row names to a matrix 
#' @param x A matrix, a data frame, or a list of matrices or data frames.
#' @param row.names An encoded vector of names with the length of nrow(x)
#' @returns Matrix with row names
#' @export
setRowNames <- function(x, row.names) {
    rn <- .decode.arg(row.names)
    if (nrow(x) != length(rn))
        stop("Cannot assign row.names of length different from nrow(x).")
    if (any(duplicated(rn)))
        stop("Repeated rownames.")
    if (is.list(x) && !is.data.frame(x)) {
        x <- lapply(x, function(xx) {
            rownames(xx) <- rn
            return (xx)
        })
    } else {
        rownames(x) <- rn
    }
    return (x)
}


#' @title Row names
#' @description Get row names of a matrix 
#' @param x A matrix, a data frame, or a list of matrices or data frames.
#' @returns Row names of x
#' @export
rowNames <- function(x) {
    if (is.list(x) && !is.data.frame(x)) {
        return (lapply(x, rownames))
    }
    return (rownames(x))
}


#' @title Colnames
#' @description Get column names of a matrix 
#' @param x A matrix, a data frame, or a list of matrices or data frames.
#' @returns Column names of x
#' @export
colNames <- function(x) {
    if (is.list(x) && !is.data.frame(x)) {
        return (lapply(x, colnames))
    }
    return (colnames(x))
}


#' @title Matrix centering
#' @description Center matrix columns or rows to 0
#' @param x A numeric matrix or data frame. A list of matrices or data frames
#' can also be provided, where they will be cbind-ed to perform the centering.
#' @param subset Encoded value of an index vector indicating the subset of
#' individuals to consider. Default, NULL, all individuals are considered.
#' @param byColumn A logical value indicating whether the input data is
#' centered by column or row. Default, TRUE, centering by column. Constant
#' variables across samples are removed. If FALSE, centering and scaling by
#' row. Constant samples across variables are removed.
#' @param scale A logical value indicating whether the variables should be
#' scaled to have unit variance. Default, FALSE.
#' @returns The centered matrix.
#' @export
center <- function(x, subset = NULL, byColumn = TRUE, scale = FALSE) {
    ## convert x to numeric matrix
    if (is.list(x) && !is.data.frame(x)) {
        y <- lapply(x, function(xx) apply(xx, c(1,2), as.numeric))
    } else {
        y <- list(apply(x, c(1,2), as.numeric))
    }
    rn <- lapply(y, rownames)
    cn <- lapply(y, colnames)
    if (min(lengths(rn))==0 || min(lengths(cn))==0)
        stop("Input data should have column (variable)
             names and row (sample) names.")
    if (min(lengths(rn)) < 3 || min(lengths(cn)) < 3)
        stop("Input data should have at least 3 columns (variables) and
             3 rows (samples).")
    if (min(lengths(rn))!=max(lengths(rn)))
        stop("Input data should have the same number of rows (samples).")
    if (max(apply(do.call(cbind, rn),
                  1,
                  function(rni) length(unique(rni)))) > 1)
        stop("Input data blocks should have the same order or rows (samples).")
    if (length(cn) > 1 && length(Reduce(intersect, cn)) > 0)
        stop("Input data blocks should have different column names.")
    ## check for constant values across samples or variables
    if (isTRUE(byColumn)) {
        invisible(lapply(y, function(yy) {
            lenuni <- apply(yy, 2, function(yc) length(unique(yc)))
            if (min(lenuni)==1) stop(paste0("Constant variables: ",
                                            colnames(yy)[which(lenuni==1)]))
        }))
    } else {
        invisible(lapply(y, function(yy) {
            lenuni <- apply(yy, 1, function(yc) length(unique(yc)))
            if (min(lenuni)==1) stop(paste0("Constant samples: ",
                                            rownames(yy)[which(lenuni==1)]))
        }))
    }
    ## check for missing values
    invisible(lapply(names(y), function(yy) {
        if (any(is.na(y[[yy]]))) stop(paste0("Missing values in ", yy))
    }))
    ## subsetting
    subset <- .decode.arg(subset)
    if (!is.null(subset)) {
        y <- lapply(y, function(yy) {
            yy[subset, , drop=F]
        })
    }
    ## ordering
    y <- lapply(y, function(yy) {
        yy[order(rownames(yy)), , drop=F]
    })
    ## centering
    y <- lapply(y, function(yy) {
        if (isTRUE(byColumn)) {
            scale(yy, center=TRUE, scale=scale)
        } else {
            t(scale(t(yy), center=TRUE, scale=TRUE))
        }
    })

    return (y)
}


#' @title Matrix partition
#' @description Partition a matrix into blocks
#' @param x A matrix
#' @param seprow A numeric vectors indicating sizes of blocks in rows
#' @param sepcol A numeric vectors indicating sizes of blocks in columns
#' @returns List of blocks
#' @importFrom arrow arrow_table
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
            xij <- x[indrow[[i]][1]:indrow[[i]][2],
                     indcol[[j]][1]:indcol[[j]][2], drop=F]
            return (arrow_table(as.data.frame(xij)))
        })
    })
    return (parMat)
}


#' @title Singular product
#' @description Product of t(x) and first column of x \%*\% t(x)
#' @param x A list of numeric matrices
#' @returns t(x) \%*\% (x \%*\% t(x))[,1]
#' @export
singularProd <- function(x) {
    sp <- lapply(x, function(xx) {
        return (crossprod(tcrossprod(xx)[, 1, drop=F], xx))
    })
    return (sp)
}


#' @title Matrix cross product
#' @description Calculates the cross product of given matrices
#' @param x A list of numeric matrices.
#' @param y A list of numeric matrices. Default, NULL, \code{y} = \code{x}.
#' @param pair A logical value indicating pairwise cross products are computed.
#' Default, FALSE. Ignored if \code{y} is provided.
#' @param chunk Size of chunks into what the resulting matrix is partitioned.
#' Default, 500L.
#' @returns A list of cross product matrices.
#' @importFrom utils combn
#' @importFrom arrow write_to_raw
#' @importFrom parallel mclapply detectCores
#' @export
crossProd <- function(x, y = NULL, pair = FALSE, chunk = 500L) {
    if (!is.list(x) || is.data.frame(x)) stop('x should be a list of matrices')
    if (is.null(y)) {
        xblocks <- mclapply(
            1:length(x),
            mc.cores=min(length(x), detectCores()),
            function(i) {
                nblockscol <- ceiling(ncol(x[[i]])/chunk)
                sepblockscol <- rep(ceiling(ncol(x[[i]])/nblockscol),
                                    nblockscol-1)
                sepblockscol <- c(sepblockscol,
                                  ncol(x[[i]]) - sum(sepblockscol))
                
                tcpblocks <- .partitionMatrix(crossprod(x[[i]]),
                                              seprow=sepblockscol)
                
                return (lapply(tcpblocks, function(tcpb) {
                    return (lapply(tcpb, function(tcp) {
                        return (.encode.arg(write_to_raw(tcp)))
                    }))
                }))
            })
        if (pair) {
            xpairs <- combn(1:length(x), 2)
            xblocks.cross <- mclapply(
                1:ncol(xpairs),
                mc.cores=min(ncol(xpairs), detectCores()),
                function(xc) {
                    x1 <- x[[xpairs[1,xc]]]
                    x2 <- x[[xpairs[2,xc]]]
                    
                    nblocksrow <- ceiling(ncol(x1)/chunk)
                    sepblocksrow <- rep(ceiling(ncol(x1)/nblocksrow),
                                        nblocksrow-1)
                    sepblocksrow <- c(sepblocksrow,
                                      ncol(x1) - sum(sepblocksrow))
                    
                    nblockscol <- ceiling(ncol(x2)/chunk)
                    sepblockscol <- rep(ceiling(ncol(x2)/nblockscol),
                                        nblockscol-1)
                    sepblockscol <- c(sepblockscol,
                                      ncol(x2) - sum(sepblockscol))
                    
                    tcpblocks <- .partitionMatrix(crossprod(x1, x2), 
                                                  seprow=sepblocksrow, 
                                                  sepcol=sepblockscol)
                    
                    return (lapply(tcpblocks, function(tcpb) {
                        return (lapply(tcpb, function(tcp) {
                            return (.encode.arg(write_to_raw(tcp)))
                        }))
                    }))
                })
            names(xblocks.cross) <- sapply(1:ncol(xpairs), function(xc) {
                paste(names(x)[xpairs[,xc]], collapse='__')
            })
            xblocks <- c(xblocks, xblocks.cross)
        }
    } else {
        if (!is.list(y) || is.data.frame(y))
            stop('y should be a list of matrices.')
        xpairs <- t(expand.grid(1:length(x), 1:length(y)))
        xblocks <- mclapply(
            1:ncol(xpairs),
            mc.cores=min(ncol(xpairs), detectCores()),
            function(xc) {
                x1 <- x[[xpairs[1,xc]]]
                x2 <- y[[xpairs[2,xc]]]
                
                nblocksrow <- ceiling(ncol(x1)/chunk)
                sepblocksrow <- rep(ceiling(ncol(x1)/nblocksrow), nblocksrow-1)
                sepblocksrow <- c(sepblocksrow, ncol(x1) - sum(sepblocksrow))
                
                nblockscol <- ceiling(ncol(x2)/chunk)
                sepblockscol <- rep(ceiling(ncol(x2)/nblockscol), nblockscol-1)
                sepblockscol <- c(sepblockscol, ncol(x2) - sum(sepblockscol))
                
                tcpblocks <- .partitionMatrix(crossprod(x1, x2), 
                                              seprow=sepblocksrow, 
                                              sepcol=sepblockscol)
                
                return (lapply(tcpblocks, function(tcpb) {
                    return (lapply(tcpb, function(tcp) {
                        return (.encode.arg(write_to_raw(tcp)))
                    }))
                }))
            })
        names(xblocks) <- sapply(1:ncol(xpairs), function(xc) {
            paste(c(names(x)[xpairs[1,xc]], names(y)[xpairs[2,xc]]),
                  collapse='__')
        })
    }
    return (xblocks)
}


#' @title Matrix cross product
#' @description Calculates the cross product x \%*\% t(y).
#' @param x A list of numeric matrices.
#' @param y A list of numeric matrices. Default, NULL, y = x.
#' @param chunk Size of chunks into what the resulting matrix is partitioned.
#' Default, 500L.
#' @returns \code{x \%*\% t(y)}
#' @importFrom arrow write_to_raw
#' @importFrom parallel mclapply detectCores
#' @export
tcrossProd <- function(x, y = NULL, chunk = 500L) {
    if (!is.list(x) || is.data.frame(x))
        stop('x should be a list of matrices.')
    if (is.null(y)) {
        etcpblocks <- mclapply(
            1:length(x),
            mc.cores=min(length(x), detectCores()),
            function(i) {
                nblocksrow <- ceiling(nrow(x[[i]])/chunk)
                sepblocksrow <- rep(ceiling(nrow(x[[i]])/nblocksrow),
                                    nblocksrow-1)
                sepblocksrow <- c(sepblocksrow,
                                  nrow(x[[i]]) - sum(sepblocksrow))
                
                tcpblocks <- .partitionMatrix(tcrossprod(x[[i]]),
                                              seprow=sepblocksrow,
                                              sepcol=sepblocksrow)
                
                return (lapply(tcpblocks, function(tcpb) {
                    return (lapply(tcpb, function(tcp) {
                        return (.encode.arg(write_to_raw(tcp)))
                    }))
                }))
        })
        names(etcpblocks) <- names(x)
    } else {
        if (!is.list(y) || is.data.frame(y))
            stop('y should be a list of matrices.')
        if (is.list(y[[1]])) {
            etcpblocks <- mclapply(
                1:length(y),
                mc.cores=min(length(y), detectCores()),
                function(k) {
                    etcpblockk <- mclapply(
                        1:length(x),
                        mc.cores=min(length(x), detectCores()),
                        function(i) {
                            nblocksrow <- ceiling(nrow(x[[i]])/chunk)
                            sepblocksrow <- rep(
                                ceiling(nrow(x[[i]])/nblocksrow),
                                nblocksrow-1)
                            sepblocksrow <- c(
                                sepblocksrow,
                                nrow(x[[i]]) - sum(sepblocksrow))
                            
                            nblockscol <- ceiling(nrow(y[[k]][[i]])/chunk)
                            sepblockscol <- rep(
                                ceiling(nrow(y[[k]][[i]])/nblockscol),
                                nblockscol-1)
                            sepblockscol <- c(
                                sepblockscol,
                                nrow(y[[k]][[i]]) - sum(sepblockscol))
                            
                            tcpblocks <- .partitionMatrix(
                                tcrossprod(x[[i]], y[[k]][[i]]),
                                seprow=sepblocksrow,
                                sepcol=sepblockscol)
                            
                            return (lapply(tcpblocks, function(tcpb) {
                                return (lapply(tcpb, function(tcp) {
                                    return (.encode.arg(write_to_raw(tcp)))
                                }))
                            }))
                        })
                    names(etcpblockk) <- names(x)
                    return (etcpblockk)
                })
            names(etcpblocks) <- names(y)
        } else if (all(names(x)==names(y))) {
            etcpblocks <- mclapply(
                1:length(x),
                mc.cores=min(length(x), detectCores()),
                function(i) {
                    nblocksrow <- ceiling(nrow(x[[i]])/chunk)
                    sepblocksrow <- rep(ceiling(nrow(x[[i]])/nblocksrow),
                                        nblocksrow-1)
                    sepblocksrow <- c(sepblocksrow,
                                      nrow(x[[i]]) - sum(sepblocksrow))
                    
                    nblockscol <- ceiling(nrow(y[[i]])/chunk)
                    sepblockscol <- rep(ceiling(nrow(y[[i]])/nblockscol),
                                        nblockscol-1)
                    sepblockscol <- c(sepblockscol,
                                      nrow(y[[i]]) - sum(sepblockscol))
                    
                    tcpblocks <- .partitionMatrix(tcrossprod(x[[i]], y[[i]]),
                                                  seprow=sepblocksrow,
                                                  sepcol=sepblockscol)
                    
                    return (lapply(tcpblocks, function(tcpb) {
                        return (lapply(tcpb, function(tcp) {
                            return (.encode.arg(write_to_raw(tcp)))
                        }))
                    }))
                })
            names(etcpblocks) <- names(x)
        } else {
            stop("Wrong format of x and y.")
        }
    }
    return (etcpblocks)
}


#' @title Bigmemory description of a matrix
#' @description Bigmemory description of a matrix
#' @param value Encoded value of a matrix
#' @returns Bigmemory description of the given matrix
#' @importFrom arrow read_ipc_stream
#' @importFrom bigmemory as.big.matrix describe
#' @export
matrix2DscFD <- function(value) {
    tcp <- as.matrix(read_ipc_stream(.decode.arg(value)))
    dscbigmatrix <- describe(as.big.matrix(tcp, backingfile = ""))
    rm(list=c("tcp"))
    return (dscbigmatrix)
}


#' @title Bigmemory description of a matrix
#' @description Bigmemory description of a matrix
#' @param value Encoded value of a matrix
#' @returns Bigmemory description of the given matrix
#' @importFrom arrow read_ipc_stream
#' @importFrom bigmemory as.big.matrix describe
#' @export
matrix2DscMate <- function(value) {
    tcp <- as.matrix(read_ipc_stream(.decode.arg(value)))
    dscbigmatrix <- describe(as.big.matrix(tcp, backingfile = ""))
    rm(list=c("tcp"))
    return (dscbigmatrix)
}


#' @title Matrix reconstruction
#' @description Rebuild a matrix from its partition.
#' @param matblocks List of lists of matrix blocks, obtained from
#' .partitionMatrix.
#' @param mc.cores Number of cores for parallel computing. Default, 1.
#' @returns The reconstructed matrix.
#' @importFrom parallel mclapply
#' @keywords internal
.rebuildMatrix <- function(matblocks, mc.cores = 1) {
    uptcp <- lapply(matblocks, function(bl) do.call(cbind, bl))
    ## combine the blocks into one matrix
    if (length(uptcp)>1) {
        if (length(unique(sapply(uptcp, ncol)))==1) {
            tcp <- do.call(rbind, uptcp)
        } else {
            ## without the first layer of blocks
            no1tcp <- mclapply(
                2:length(uptcp),
                mc.cores=mc.cores,
                function(i) {
                    cbind(do.call(cbind,
                                  lapply(1:(i-1),
                                         function(j) {
                                             t(matblocks[[j]][[i-j+1]])
                                         })), uptcp[[i]])
                })
            ## with the first layer of blocks
            tcp <- rbind(uptcp[[1]], do.call(rbind, no1tcp))
            rm(list=c("no1tcp"))
            rownames(tcp) <- colnames(tcp)
            stopifnot(isSymmetric(tcp))
        }
    } else {
        ## for either asymmetric matrix or column-wise blocks
        tcp <- uptcp[[1]]
    }
    rm(list=c("uptcp"))
    
    return (tcp)
}


#' @title Matrix reconstruction
#' @description Rebuild a matrix from its partition bigmemory objects.
#' @param dscblocks List of lists of bigmemory objects pointed to matrix
#' blocks.
#' @param mc.cores Number of cores for parallel computing. Default, 1.
#' @returns The complete symmetric matrix
#' @importFrom parallel mclapply
#' @importFrom bigmemory attach.big.matrix
#' @keywords internal
.rebuildMatrixDsc <- function(dscblocks, mc.cores = 1) {
    ## access to matrix blocks 
    matblocks <- mclapply(dscblocks, mc.cores=mc.cores, function(y) {
        lapply(y, function(x) {
            bmx <- (attach.big.matrix(x))[,,drop=F]
            rm(x)
            gc()
            return (bmx)
        })
    })
    tcp <- .rebuildMatrix(matblocks, mc.cores=mc.cores)
    return (tcp)
}


#' @title Matrix reconstruction
#' @description Rebuild a matrix from its partition through variables
#' in the environment.
#' @param symbol Generic variable name
#' @param len1 First-order length of lists of matrix blocks. Variables were
#' generated as symbol__1__1__1, symbol__1__1__2, etc.
#' @param len2 Second-order length of lists of matrix blocks.
#' @param len3 Third-order length of lists of matrix blocks.
#' @param querytables Names to be assigned to the result. 
#' @param mc.cores Number of cores for parallel computing. Default, 1.
#' @returns The complete symmetric matrix
#' @importFrom bigmemory attach.big.matrix
#' @importFrom parallel mclapply
#' @export
rebuildMatrixVar <- function(symbol, len1, len2, len3,
                             querytables = NULL,
                             mc.cores = 1) {
    ## access to matrix blocks
    matblocks <- mclapply(1:len1, mc.cores=mc.cores, function(i) {
        lapply(1:len2[i], function(j) {
            lapply(1:len3[[i]][j], function(k) {
                varname <- paste(c(symbol, i, j, k), collapse="__")
                dscblock <- get(varname, pos=1)
                bmk <- (attach.big.matrix(dscblock))[,,drop=F]
                rm(varname)
                rm(dscblock)
                gc()
                return (bmk)
            })
        })
    })
    tcp <- mclapply(1:len1, mc.cores=mc.cores, function(i) {
        .rebuildMatrix(matblocks[[i]], mc.cores=mc.cores)
    })
    names(tcp) <- querytables
    
    return (tcp)
}


#' @title Matrix triple product
#' @description Calculate the triple product x \%*\% y \%*\% t(x).
#' @param x A list of numeric matrices.
#' @param mate Mate server name, from which y was pushed.
#' @param chunk Size of chunks into what the resulting matrix is partitioned.
#' Default, 500L.
#' @param mc.cores Number of cores for parallel computing. Default, 1.
#' @returns \code{x \%*\% t(y)}
#' @importFrom arrow write_to_raw
#' @importFrom parallel mclapply detectCores
#' @export
tripleProdChunk <- function(x, mate, chunk = 500L, mc.cores = 1) {
    nblocks <- ceiling(nrow(x[[1]])/chunk)
    sepblocks <- rep(ceiling(nrow(x[[1]])/nblocks), nblocks-1)
    sepblocks <- c(sepblocks, nrow(x[[1]]) - sum(sepblocks))
    mc.cores <- min(mc.cores, detectCores())
    
    tpcs <- mclapply(
        1:length(x),
        mc.cores=min(length(x), mc.cores),
        function(i) {
            y <- get(paste("pushed", mate, sep='_'), pos=1)
            stopifnot(isSymmetric(y[[i]], check.attributes=F))
            ## NB. this computation of tcpblocks could be done more efficiently
            ## with y being a chunked matrix in bigmemory
            tp <- tcrossprod(x[[i]], tcrossprod(x[[i]], y[[i]]))
            tcpblocks <- .partitionMatrix(tp,
                                          seprow=sepblocks,
                                          sepcol=sepblocks)
            etcpblocks <- lapply(tcpblocks, function(tcpb) {
                return (lapply(tcpb, function(tcp) {
                    return (.encode.arg(write_to_raw(tcp)))
                }))
            })
            return (etcpblocks)
        })
    names(tpcs) <- names(x)
    
    return (tpcs)
}


#' @title Cross login
#' @description Call datashield.login on remote servers
#' @param logins An encoded dataframe with server, url, user, password, driver,
#' and options fields.
#' @returns Object(s) of class OpalConnection.
#' @export
crossLogin <- function(logins) {
    loginfo <- .decode.arg(logins)
    myDf <- data.frame(server=loginfo$server,
                       url=loginfo$url,
                       user=loginfo$user,
                       password=loginfo$password,
                       driver=loginfo$driver,
                       options=loginfo$options)
    return (.login(myDf)) #datashield.login(myDf)
}


#' @title Cross logout
#' @description Call datashield.logout on remote servers.
#' @param conns A list of Opal connections.
#' @importFrom DSI datashield.logout
#' @export
crossLogout <- function(conns) {
    datashield.logout(conns)
}


#' @title Cross aggregate
#' @description Call datashield.aggregate on remote servers.
#' @param conns A list of Opal connections.
#' @param expr An encoded expression to evaluate.
#' @param async See DSI::datashield.aggregate options. Default, TRUE.
#' @importFrom DSI datashield.aggregate
#' @export
crossAggregatePrimal <- function(conns, expr, async = T) {
    expr <- .decode.arg(expr)
    if (grepl("^singularProd\\(", expr)) {
        return (datashield.aggregate(conns=conns,
                                     expr=as.symbol(expr),
                                     async=async))
    } else {
        stop(paste0("Failed to execute: ", expr))
    }
}


#' @title Cross aggregate through two-layer connection
#' @description Call datashield.aggregate on remote servers through two-layer
#' connection.
#' @param conns A list of Opal connections.
#' @param expr An encoded expression to evaluate. This is restricted to
#' \code{pushToDscFD}.
#' @param async See DSI::datashield.aggregate options. Default, TRUE.
#' @importFrom DSI datashield.aggregate
#' @export
crossAggregateDual <- function(conns, expr, async = T) {
    expr <- .decode.arg(expr)
    if (grepl("^pushToDscFD\\(", expr)) {
        return (datashield.aggregate(conns=conns,
                                     expr=as.symbol(expr),
                                     async=async))
    } else {
        stop(paste0("Failed to execute: ", expr))
    }
}


#' @title Bigmemory description of a pushed object
#' @description Bigmemory description of a pushed object
#' @param conns A list of Opal connections. Only one connection is accepted.
#' @param object The object to be pushed.
#' @param async See DSI::datashield.aggregate options. Default, TRUE.
#' @returns Bigmemory description of the pushed object on \code{conns}.
#' @importFrom DSI datashield.aggregate
#' @export
pushToDscFD <- function(conns, object, async = T) {
    ## TODO: check for allowed conns
    stopifnot(is.list(conns) && length(conns)==1 &&
                  class(conns[[1]])=="OpalConnection")
    ## check object format
    if (!is.list(object))
        stop("object is not a list.")
    if (length(setdiff(unlist(lapply(object, class)), "list")) > 0)
        stop("object is not a list of lists")
    if (length(setdiff(unlist(lapply(object, function(x)
        lapply(x, class))), "list")) > 0)
        stop("object is not a list of lists of lists")
    
    chunkList <- object
    ## TODO: mclapply
    if (!is.list(object[[1]][[1]][[1]])) {
        dsc <- lapply(chunkList, function(clomics) {
            return (lapply(clomics, function(clrow) {
                return (lapply(clrow, function(clcol) {
                    expr <- list(as.symbol("matrix2DscFD"), clcol)
                    clcoldsc <- datashield.aggregate(conns=conns,
                                                     expr=as.call(expr),
                                                     async=async)
                    .printTime(paste0("pushToDscFD ", datashield.errors()))
                    return (clcoldsc[[1]])
                }))
            }))
        })
    } else {
        dsc <- lapply(chunkList, function(clnodes) {
            return (lapply(clnodes, function(clomics) {
                return (lapply(clomics, function(clrow) {
                    return (lapply(clrow, function(clcol) {
                        expr <- list(as.symbol("matrix2DscFD"), clcol)
                        clcoldsc <- datashield.aggregate(conns=conns,
                                                         expr=as.call(expr),
                                                         async=async)
                        .printTime(paste0("pushToDscFD ", datashield.errors()))
                        return (clcoldsc[[1]])
                    }))
                }))
            }))
        })
    }
    names(dsc) <- names(chunkList)

    return (dsc)
}


#' @title Bigmemory description of a pushed object
#' @description Bigmemory description of a pushed object
#' @param conns A list of Opal connections.
#' @param object The object to be pushed.
#' @param sourcename Name of the pushed object source.
#' @param async See DSI::datashield.assign options. Default, TRUE.
#' @returns Bigmemory description of the pushed object on \code{conns}.
#' @importFrom DSI datashield.assign
#' @export
pushToDscMate <- function(conns, object, sourcename, async = T) {
    ## TODO: check for allowed conns
    stopifnot(is.list(conns) &&
                  length(setdiff(unique(sapply(conns, class)),
                                 "OpalConnection"))==0)

    ## check object format
    if (!is.list(object)) 
        stop("object is not a list.")
    if (length(setdiff(unlist(lapply(object, class)), "list")) > 0)
        stop("object is not a list of lists")
    if (length(setdiff(unlist(lapply(object, function(x) 
        lapply(x, class))), "list")) > 0)
        stop("object is not a list of lists of lists")
    
    chunkList <- object
    
    invisible(lapply(1:length(chunkList), function(i) {
        lapply(1:length(chunkList[[i]]), function(j) {
            lapply(1:length(chunkList[[i]][[j]]), function(k) {
                datashield.assign(conns, paste(c(sourcename, i, j, k),
                                               collapse="__"),
                                  as.call(list(as.symbol("matrix2DscMate"),
                                               chunkList[[i]][[j]][[k]])),
                                  async=async)
            })
        })
    }))
    datashield.assign(conns, paste("pushed", sourcename, sep="_"),
                      as.call(list(as.symbol("rebuildMatrixVar"),
                                   symbol=sourcename,
                                   len1=length(chunkList),
                                   len2=lengths(chunkList),
                                   len3=lapply(chunkList, lengths),
                                   querytables=names(chunkList))),
                      async=async)
}


#' @title Cross assign
#' @description Call datashield.assign on remote servers.
#' @param conns A list of Opal connections.
#' @param symbol Name of an R symbol.
#' @param value A variable name or an R expression with allowed assign function
#' calls.
#' @param value.call A logical value, TRUE if value is function call, FALSE
#' if value is a variable name.
#' @param async See DSI::datashield.assign options. Default, TRUE.
#' @importFrom DSI datashield.assign
#' @export
crossAssign <- function(conns, symbol, value, value.call, async = T) {
    if (is.call(value)) {
        datashield.assign(conns=conns,
                          symbol=symbol,
                          value=value,
                          async=async)
    } else {
        valued <- .decode.arg(value)
        datashield.assign(conns=conns, symbol=symbol,
                          value=ifelse(value.call, as.symbol(valued), valued),
                          async=async)
    }
}


#' @title Cross assign a preprocessing function
#' @description Call preprocessing function on remote servers.
#' @param conns A list of Opal connections.
#' @param func Encoded definition of a function for preparation of raw data
#' matrices. Two arguments are required: conns (list of Opal connections), 
#' symbol (name of the R symbol) (see datashield.assign).
#' @param symbol Encoded vector of names of R symbols to assign in the
#' DataSHIELD R session on each server in \code{logins}.
#' The assigned R variables will be used as the input raw data to compute
#' covariance matrix. Other assigned R variables in \code{func} are ignored.
#' @importFrom DSI datashield.logout
#' @export
crossAssignFunc <- function(conns, func, symbol) {
    funcPreProc <- .decode.arg(func)
    querytables <- .decode.arg(symbol)
    
    ## create data X with funcPreProc
    tryCatch({
        ## take a snapshot of the current session
        safe.objs <- .ls.all()
        ## leave alone .Random.seed for sample()
        safe.objs[['.GlobalEnv']] <- setdiff(safe.objs[['.GlobalEnv']],
                                             '.Random.seed')
        ## lock everything so no objects can be changed
        .lock.unlock(safe.objs, lockBinding)
        
        ## apply funcPreProc for preparation of querytables on conns
        ## TODO: control hacking!
        ## TODO: control identical colnames!
        funcPreProc(conns=conns, symbol=querytables)
        
        ## unlock back everything
        .lock.unlock(safe.objs, unlockBinding)
        ## get rid of any sneaky objects that might have been created in the
        ## filters as side effects
        .cleanup(safe.objs)
    }, error=function(e) {
        .printTime(paste0("DATA MAKING CROSS PROCESS: ", e))
        datashield.logout(conns)
    })
    return (NULL)
}


#' @title Federated covariance matrix
#' @description Compute the covariance matrix for the virtual cohort
#' @param loginFD Login information of the FD server (one of the servers
#' containing cohort data)
#' @param logins Login information of other servers containing cohort data
#' @param funcPreProc Definition of a function for preparation of raw data
#' matrices. Two arguments are required: conns (list of Opal connections),
#' symbol (name of the R symbol) (see datashield.assign).
#' @param querytables Name (or a vector of two names) of the R symbol(s) to
#' assign in the Datashield R session on each server in \code{logins}.
#' The assigned R variable(s) will be used as the input raw data to compute
#' covariance matrix. Other assigned R variables in \code{func} are ignored.
#' @param querysubset A list of index vectors indicating the subsets of
#' individuals to consider. Default, NULL, all individuals are considered.
#' @param pair A logical value indicating pairwise cross products are computed.
#' Default, FALSE.
#' @param chunk Size of chunks into what the resulting matrix is partitioned.
#' Default, 500L.
#' @param mc.cores Number of cores for parallel computing. Default, 1.
#' @param connRes A logical value indicating if the connection to \code{logins}
#' is returned. Default, no. 
#' @returns Covariance matrix of the virtual cohort
#' @importFrom DSI datashield.aggregate datashield.assign
#' datashield.logout datashield.symbols datashield.errors
#' @importFrom parallel mclapply
#' @importFrom utils combn
#' @keywords internal
.federateCov <- function(loginFD, logins, funcPreProc, querytables,
                         querysubset = NULL, pair = FALSE,
                         chunk = 500L, mc.cores = 1, connRes = FALSE) {
    loginFDdata <- .decode.arg(loginFD)
    logindata   <- .decode.arg(logins)
    opals <- .login(logins=logindata)
    nnode <- length(opals)
    mc.cores <- min(mc.cores, detectCores())
    mc.nodes <- min(mc.cores, nnode)
    .printTime(".federateCov Login-ed")
    
    ## create data X with funcPreProc
    tryCatch({
        ## take a snapshot of the current session
        safe.objs <- .ls.all()
        ## leave alone .Random.seed for sample()
        safe.objs[['.GlobalEnv']] <- setdiff(safe.objs[['.GlobalEnv']],
                                             '.Random.seed')
        ## lock everything so no objects can be changed
        .lock.unlock(safe.objs, lockBinding)
        ## apply funcPreProc for preparation of querytables on opals
        ## TODO: control hacking!
        ## TODO: control identical colnames!
        funcPreProc(conns=opals, symbol=querytables)
        ## unlock back everything
        .lock.unlock(safe.objs, unlockBinding)
        ## get rid of any sneaky objects that might have been created in the
        ## filters as side effects
        .cleanup(safe.objs)
    }, error=function(e) {
        .printTime(paste0("DATA MAKING PROCESS: ", e))
        return (paste0("DATA MAKING PROCESS: ", e,
                       ' --- ', datashield.symbols(opals),
                       ' --- ', datashield.errors(),
                       ' --- ', datashield.logout(opals)))
    })
    .printTime(".federateCov Input data processed")
    
    ## compute X'X on opals
    tryCatch({
        ## center data
        if (is.null(querysubset)) {
            datashield.assign(opals,
                              "centeredData",
                              as.call(c(as.symbol("center"),
                                        x=as.call(c(as.symbol("list"),
                                                    setNames(
                                                        lapply(querytables,
                                                               as.symbol),
                                                        querytables))))),
                              async=T)
        } else {
            stopifnot(all(names(opals)==names(querysubset)))
            invisible(mclapply(names(opals), mc.cores=mc.nodes, function(opn) {
                datashield.assign(
                    opals[opn],
                    "centeredData",
                    as.call(c(as.symbol("center"),
                              x=as.call(c(as.symbol("list"),
                                          setNames(
                                              lapply(querytables,
                                                     as.symbol),
                                              querytables)
                              )),
                              subset=.encode.arg(querysubset[[opn]]))),
                    async=T)
            }))
        }
        
        ## number of samples
        nsamples <- sapply(
            datashield.aggregate(opals,
                                 as.symbol('dsDim(centeredData)'),
                                 async=T),
            function(x) x[[1]][1])
        
        ## variables
        variables <- datashield.aggregate(opals[1],
                                          as.symbol('colNames(centeredData)'),
                                          async=T)[[1]]

        .printTime(".federateCov Dimension retrieved")
        
        ## compute X'X
        datashield.assign(opals, "crossProdSelf", 
                          as.call(list(as.symbol("crossProd"),
                                       x=as.symbol("centeredData"),
                                       pair=pair,
                                       chunk=chunk)),
                          async=T)
        .printTime(".federateCov X'X computed")
        
        ## connection from non-FD servers to FD-assigned server:
        ## user and password for login between servers are required
        loginFDdata$user     <- loginFDdata$userserver
        loginFDdata$password <- loginFDdata$passwordserver
        datashield.assign(opals, 'FD',
                          as.symbol(paste0("crossLogin('",
                                           .encode.arg(loginFDdata), "')")),
                          async=T)
    }, error=function(e) {
        .printTime(paste0("COV MAKING PROCESS: ", e))
        return (paste0("COV MAKING PROCESS: ", e,
                       ' --- ', datashield.symbols(opals),
                       ' --- ', datashield.errors(),
                       ' --- ', datashield.logout(opals)))
    })
    
    ## send X'X from opals to FD
    tryCatch({
        ## send decomposed X'X to FD memory
        command <- list(as.symbol("pushToDscFD"),
                        as.symbol("FD"),
                        as.symbol("crossProdSelf"),
                        async=T)
        .printTime("Command: pushToDscFD(FD, 'crossProdSelf')")
        crossProdSelfDSC <- datashield.aggregate(opals,
                                                 as.call(command),
                                                 async=T)
        .printTime("X'X communicated to FD: ")
        ## names of X'X
        if (pair) {
            xpairs <- combn(1:length(querytables), 2)
            crossNames <- sapply(1:ncol(xpairs), function(xc) {
                paste(querytables[xpairs[,xc]], collapse='__')
            })
        } else {
            crossNames <- c()
        }
        crossProdNames <- c(querytables, crossNames)
        ## rebuild X'X
        crossProdSelf <- mclapply(
            names(opals),
            mc.cores=mc.nodes,
            function(opn) {
                cps <- lapply(crossProdSelfDSC[[opn]], function(dscblocki) {
                    return (.rebuildMatrixDsc(
                        dscblocki,
                        mc.cores=max(1, floor(mc.cores/mc.nodes))))
                })
                if (is.null(names(cps))) names(cps) <- crossProdNames
                return (cps)
            })
        ## compute X'X on virtual cohort
        rescov <- mclapply(
            crossProdNames,
            mc.cores=mc.cores, function(qtabi) {
                Reduce("+",
                       lapply(crossProdSelf,
                              function(cps)
                                  cps[[qtabi]]))/(sum(nsamples)-1)
            })
        names(rescov) <- crossProdNames
        ## set rownames and colnames of X'X
        for (crn in querytables) {
            if (is.null(rownames(rescov[[crn]])))
                rownames(rescov[[crn]]) <- variables[[crn]]
            if (is.null(colnames(rescov[[crn]])))
                colnames(rescov[[crn]]) <- variables[[crn]]
        }
        for (crn in crossNames) {
            crntype <- strsplit(crn, split="__")[[1]]
            if (is.null(rownames(rescov[[crn]])))
                rownames(rescov[[crn]]) <- 
                    rownames(rescov[[crntype[1]]])
            if (is.null(colnames(rescov[[crn]])))
                colnames(rescov[[crn]]) <- 
                    colnames(rescov[[crntype[2]]])
            
        }
    }, error=function(e) {
        .printTime(paste0("COV PUSH PROCESS: ", e,
                          ' --- ', datashield.errors()))
        datashield.assign(opals, 'crossEnd',
                          as.symbol("crossLogout(FD)"), async=T)
        return (paste0("COV PUSH PROCESS: ", e,
                       ' --- ', datashield.symbols(opals),
                       ' --- ', datashield.errors(),
                       ' --- ', datashield.logout(opals)))
    })
    
    if (!connRes) {
        datashield.assign(opals, 'crossEnd',
                          as.symbol("crossLogout(FD)"), async=T)
        datashield.logout(opals)
    }

    return (list(cov=rescov,
                 conns=opals))
}


#' @title Federated PCA
#' @description Perform the principal component analysis for the virtual cohort.
#' @usage federatePCA(loginFD,
#'                    logins,
#'                    func,
#'                    symbol,
#'                    ncomp = 2,
#'                    chunk = 500L,
#'                    mc.cores = 1)
#' @param loginFD Login information of the FD server (one of the servers
#' containing cohort data).
#' @param logins Login information of other servers containing cohort data.
#' @param func Encoded definition of a function for preparation of raw data
#' matrices. Two arguments are required: conns (list of Opal connections), 
#' symbol (name of the R symbol) (see datashield.assign).
#' @param symbol Encoded vector of names of the R symbols to assign in the
#' DataSHIELD R session on each server in \code{logins}.
#' The assigned R variable will be used as the input raw data to compute
#' covariance matrix for PCA. Other assigned R variables in \code{func} are
#' ignored.
#' @param ncomp Number of components. Default, 2.
#' @param chunk Size of chunks into what the resulting matrix is partitioned.
#' Default, 500L.
#' @param mc.cores Number of cores for parallel computing. Default, 1.
#' @returns PCA object
#' @importFrom DSI datashield.aggregate datashield.assign datashield.logout
#' @importFrom stats princomp
#' @importFrom arrow write_to_raw
#' @importFrom parallel mclapply
#' @examples
#' \dontrun{
#' dataProc <- function(conns, symbol) {
#'     DSI::datashield.assign(conns, symbol, 'test.CNSIM',
#'         variables=c('LAB_TRIG', 'LAB_GLUC_ADJUSTED', 'PM_BMI_CONTINUOUS'))
#' }
#' federatePCA(.encode.arg(loginFD), .encode.arg(logins),
#'             .encode.arg(dataProc, serialize.it = T),
#'             .encode.arg("rawData"))
#' }
#' @export
federatePCA <- function(loginFD, logins, func, symbol, ncomp = 2,
                        chunk = 500L, mc.cores = 1) {
    .printTime("federatePCA started")
    if (ncomp < 2) {
        print("ncomp should be at least 2. ncomp will be set to 2.")
        ncomp <- 2
    }
    funcPreProc <- .decode.arg(func)
    querytables <- .decode.arg(symbol)
    ntab <- length(querytables)
    mc.cores <- min(mc.cores, detectCores())
    
    ## compute covariance matrices for the virtual cohort
    fedCov <- .federateCov(loginFD, logins, funcPreProc, querytables,
                           chunk=chunk, mc.cores=mc.cores, connRes=T)
    
    ## pca object
    pcaObjs <- mclapply(
        querytables,
        mc.cores=min(ntab, mc.cores),
        function(tab) {
            pcaObj <- princomp(covmat=fedCov$cov[[tab]])
            
            nsdevpos <- length(which(pcaObj$sdev > 1e-6))
            if (nsdevpos <= ncomp) {
                print(paste0("Security issue: maximum ", nsdevpos - 1, 
                             " components could be inquired."))
                ncomp <- nsdevpos - 1
                if (ncomp < 2) stop("Less than 2 components were found.")
            }
            
            pcaObj$loadings <- pcaObj$loadings[, 1:ncomp, drop=F]
            return (pcaObj)
        })
    names(pcaObjs) <- querytables
    
    ## pca loadings
    loadings <- mclapply(
        querytables,
        mc.cores=min(ntab, mc.cores),
        function(tab) {
            xx <- t(pcaObjs[[tab]]$loadings)
            nblockscol <- ceiling(ncol(xx)/chunk)
            sepblockscol <- rep(ceiling(ncol(xx)/nblockscol), nblockscol-1)
            sepblockscol <- c(sepblockscol, ncol(xx) - sum(sepblockscol))
            tcpblocks <- .partitionMatrix(xx, seprow=ncomp, sepcol=sepblockscol)
            return (lapply(tcpblocks, function(tcpb) {
                return (lapply(tcpb, function(tcp) {
                    return (.encode.arg(write_to_raw(tcp)))
                }))
            }))
        })
    names(loadings) <- querytables
    
    ## send loadings back to non-FD servers
    tryCatch({
        opals <- fedCov$conns
        nnode <- length(opals)
        pushToDscMate(conns=opals, object=loadings, sourcename='FD', async=T)
        
        ## compute X*loadings'
        datashield.assign(opals, "scores", 
                          as.call(list(as.symbol("tcrossProd"),
                                       x=as.symbol("centeredData"),
                                       y=as.symbol("pushed_FD"),
                                       chunk=chunk)),
                          async=T)
        
        command <- list(as.symbol("pushToDscFD"),
                        as.symbol("FD"),
                        as.symbol("scores"),
                        async=T)
        cat("Command: pushToDscFD(FD, scores)", "\n")
        scoresDSC <- datashield.aggregate(opals, as.call(command), async=T)
        .printTime(paste0("scores communicated to FD: "))
        
        mc.nodes <- min(nnode, mc.cores)
        scoresLoc <- mclapply(
            names(opals),
            mc.cores=mc.nodes,
            function(opn) {
                mc.tabs <- max(1, min(ntab, floor(mc.cores/mc.nodes)))
                cps <- mclapply(
                    querytables,
                    mc.cores=mc.tabs,
                    function(tab) {
                        cpsi <- .rebuildMatrixDsc(
                            scoresDSC[[opn]][[tab]],
                            mc.cores=max(1,
                                         floor(mc.cores/(mc.nodes*mc.tabs))))
                        colnames(cpsi) <- paste0("Comp.", 1:ncomp)
                        return (cpsi)
                    })
                names(cps) <- querytables
                return (cps)
            })
        names(scoresLoc) <- names(opals)
        for (qtabi in querytables) {
            pcaObjs[[qtabi]]$scores <- do.call(
                rbind,
                lapply(scoresLoc, function(sl) sl[[qtabi]])
            )
        }
    }, error = function(e) {
        .printTime(paste0("LOADINGS MAKING PROCESS: ", e))
        return (paste0("LOADINGS MAKING PROCESS: ", e))
    }, finally = {
        datashield.assign(opals, 'crossEnd',
                          as.symbol("crossLogout(FD)"), async=T)
        datashield.logout(opals)
    })
    
    class(pcaObjs) <- "federatePCA"
    return (pcaObjs)
}


#' @title RCCA tuning
#' @description Estimate optimized regulation parameters lambda1 and lambda2.
#' @param loginFD Login information of the FD server (one of the servers
#' containing cohort data).
#' @param logins Login information of servers containing cohort data.
#' @param funcPreProc Definition of a function for preparation of raw data
#' matrices. Two arguments are required: conns (list of Opal connections),
#' symbol (name of the R symbol) (see datashield.assign).
#' @param querytables Name (or a vector of two names) of the R symbol(s) to
#' assign in the Datashield R session on each server in \code{logins}.
#' The assigned R variable(s) will be used as the input raw data to compute
#' covariance matrix. Other assigned R variables in \code{func} are ignored.
#' @param chunk Size of chunks into what the resulting matrix is partitioned.
#' Default, 500L.
#' @param mc.cores Number of cores for parallel computing. Default, 1.
#' @param nfold n-fold cross-validation. Default, 5.
#' @param grid1 Tuning values for \code{lambda1}.
#' @param grid2 Tuning values for \code{lambda2}.
#' @returns Optimal values of \code{lambda1} and \code{lambda2}.
#' @importFrom parallel mclapply
#' @importFrom graphics image
#' @importFrom stats cor setNames
#' @importFrom DSI datashield.logout datashield.assign datashield.aggregate
#' @importFrom arrow write_to_raw
#' @keywords internal
.estimateR <- function(loginFD, logins, funcPreProc, querytables,
                       chunk = 500L, mc.cores = 1,
                       nfold = 5,
                       grid1 = seq(0.001, 1, length = 5),
                       grid2 = seq(0.001, 1, length = 5)) {
    stopifnot(length(querytables) == 2)
    loginFDdata <- .decode.arg(loginFD)
    logindata   <- .decode.arg(logins)
    
    opals <- .login(logins=logindata)
    nnode <- length(opals)
    
    ## create data X with funcPreProc
    tryCatch({
        ## take a snapshot of the current session
        safe.objs <- .ls.all()
        ## leave alone .Random.seed for sample()
        safe.objs[['.GlobalEnv']] <- setdiff(safe.objs[['.GlobalEnv']],
                                             '.Random.seed')
        ## lock everything so no objects can be changed
        .lock.unlock(safe.objs, lockBinding)
        
        ## apply funcPreProc for preparation of querytables on opals
        ## TODO: control hacking!
        ## TODO: control identical colnames!
        funcPreProc(conns=opals, symbol=querytables)
        
        ## unlock back everything
        .lock.unlock(safe.objs, unlockBinding)
        ## get rid of any sneaky objects that might have been created in the
        ## filters as side effects
        .cleanup(safe.objs)
    }, error=function(e) {
        .printTime(paste0("DATA MAKING PROCESS: ", e))
        datashield.logout(opals)
    })
    
    ## tuning with nfold
    tryCatch({
        ## center data
        datashield.assign(opals,
                          "centeredData",
                          as.call(c(as.symbol("center"),
                                    x=as.call(c(as.symbol("list"),
                                                setNames(
                                                    lapply(querytables,
                                                           as.symbol),
                                                    querytables)))
                          )),
                          async=T)
        
        ## number of samples
        nsamples <- sapply(
            datashield.aggregate(opals,
                                 as.symbol('dsDim(centeredData)'),
                                 async=T),
            function(x) x[[1]][1])
        if (is.null(names(nsamples))) names(nsamples) <- names(opals)
        
        ## variables
        variables <- datashield.aggregate(opals[1],
                                          as.symbol('colNames(centeredData)'),
                                          async=T)[[1]]
        
        ## random Mfold partitions to leave out
        foldspar <- split(
            unlist(mclapply(names(opals),
                            mc.cores=min(nnode, mc.cores),
                            function(opn) {     
                                paste(opn, 1:nsamples[opn], sep="_")
                            }))[sample(1:sum(nsamples))], 
            rep(1:nfold, length = sum(nsamples)))
        foldslef <- lapply(foldspar, function(fl) {
            setNames(mclapply(names(opals),
                              mc.cores=min(nnode, mc.cores),
                              function(opn) {
                                  sort(as.numeric(sub(fl[grep(opn, fl)],
                                                      pattern=paste0(opn,"_"),
                                                      replacement='')))
            }), names(opals))
        })
        ## remaining individuals on each cohort
        foldsrem <- lapply(foldslef, function(fl) {
            setNames(mclapply(names(opals),
                              mc.cores=min(nnode, mc.cores),
                              function(opn) {
                                  setdiff(1:nsamples[opn], fl[[opn]])
                              }), names(opals))
        })
        grid <- expand.grid(grid1, grid2)
        cv.score <- apply(grid, 1, function(lambda) {
            xscore <- NULL
            yscore <- NULL
            for (m in 1:nfold) {
                ## compute covariance matrices for the virtual cohort
                fedCov <- .federateCov(loginFD, logins, funcPreProc,
                                       querytables, querysubset=foldsrem[[m]],
                                       pair=T, chunk=chunk, mc.cores=mc.cores,
                                       connRes=T)
                Cxx <- fedCov$cov[[1]]
                Cyy <- fedCov$cov[[2]]
                Cxy <- fedCov$cov[[3]]
                
                ## add parameters of regularization
                Cxx <- Cxx + diag(lambda[1], nrow=nrow(Cxx), ncol=ncol(Cxx))
                Cyy <- Cyy + diag(lambda[2], nrow=nrow(Cyy), ncol=ncol(Cyy))
                
                ## CCA core call
                res <- geigen(Cxy, Cxx, Cyy)
                names(res) <- c("cor", "xcoef", "ycoef")
                ncomp <- 1
                res$xcoef <- res$xcoef[, 1:ncomp, drop=F]
                res$ycoef <- res$ycoef[, 1:ncomp, drop=F]
                rownames(res$xcoef) <- rownames(Cxx)
                rownames(res$ycoef) <- rownames(Cyy)
                colnames(res$xcoef) <- colnames(res$ycoef) <- 
                    paste0("Comp.", 1:ncomp)
                
                ## rcca coefs
                loadings <- mclapply(
                    res[c("xcoef", "ycoef")],
                    mc.cores=mc.cores,
                    function(ccacoef) {
                        xx <- t(ccacoef)
                        nblockscol <- ceiling(ncol(xx)/chunk)
                        sepblockscol <- rep(ceiling(ncol(xx)/nblockscol),
                                            nblockscol-1)
                        sepblockscol <- c(sepblockscol,
                                          ncol(xx) - sum(sepblockscol))
                        tcpblocks <- .partitionMatrix(xx,
                                                      seprow=ncomp,
                                                      sepcol=sepblockscol)
                        return (lapply(tcpblocks, function(tcpb) {
                            return (lapply(tcpb, function(tcp) {
                                return (.encode.arg(write_to_raw(tcp)))
                            }))
                        }))
                    })
                names(loadings) <- querytables
                
                ## compute scores on foldslef[[m]]
                tryCatch({
                    mopals <- fedCov$conns
                    ## send coefs back to non-FD servers
                    pushToDscMate(conns=mopals, object=loadings,
                                   sourcename='FD', async=T)
                    
                    ## centeredData for foldslef
                    invisible(mclapply(
                        names(mopals), mc.cores=length(mopals), function(opn) {
                            datashield.assign(
                                mopals[opn],
                                "centeredDatam",
                                as.call(c(as.symbol("center"),
                                          x=as.call(c(as.symbol("list"),
                                                      setNames(
                                                          lapply(querytables,
                                                                 as.symbol),
                                                          querytables)
                                          )),
                                          subset=.encode.arg(
                                              foldslef[[m]][[opn]])
                                )),
                                async=T)
                        }))
                    
                    ## compute X*coefs'
                    datashield.assign(
                        mopals, "scores", 
                        as.call(list(as.symbol("tcrossProd"),
                                     x=as.symbol("centeredDatam"),
                                     y=as.symbol("pushed_FD"),
                                     chunk=chunk)),
                        async=T)
                    
                    command <- list(as.symbol("pushToDscFD"),
                                    as.symbol("FD"),
                                    as.symbol("scores"),
                                    async=T)
                    cat("Command: pushToDscFD(FD, scores)", "\n")
                    scoresDSC <- datashield.aggregate(mopals, as.call(command),
                                                      async=T)
                    .printTime(paste0("scores communicated to FD: "))
                    
                    scoresLoc <- lapply(scoresDSC, function(dscblocks) {
                        cps <- lapply(dscblocks, function(dscblocki) {
                            cpsi <- .rebuildMatrixDsc(dscblocki,
                                                      mc.cores=mc.cores)
                            colnames(cpsi) <- paste0("Comp.", 1:ncomp)
                            return (cpsi)
                        })
                        names(cps) <- c("xcoef", "ycoef")
                        return (cps)
                    })
                    cvx <- do.call(rbind, lapply(scoresLoc, function(sl)
                        sl$xcoef))
                    cvy <- do.call(rbind, lapply(scoresLoc, function(sl)
                        sl$ycoef))
                }, error = function(e) {
                    .printTime(paste0("ESTIMATE COEF PUSH PROCESS: ", m, e))
                    return (paste0("ESTIMATE COEF PUSH PROCESS: ", m, e))
                }, finally = {
                    datashield.assign(mopals, 'crossEnd',
                                      as.symbol("crossLogout(FD)"), async=T)
                    datashield.logout(mopals)
                })
                
                xscore <- c(xscore, cvx)
                yscore <- c(yscore, cvy)
            }
            return (cor(xscore, yscore, use = "pairwise"))
        })
    }, finally = datashield.logout(opals))
    
    cv.score.grid <- cbind(grid, cv.score)
    mat <- matrix(cv.score, nrow=length(grid1), ncol=length(grid2))
    #plot <- FALSE
    #if (isTRUE(plot))
    #    image(list(grid1 = grid1, grid2 = grid2, mat = mat))
    opt <- cv.score.grid[which.max(cv.score.grid[,3]), ]
    out <- list(opt.lambda1 = opt[[1]],
                opt.lambda2 = opt[[2]],
                opt.score   = opt[[3]],
                grid1       = grid1,
                grid2       = grid2,
                mat         = mat)
    class(out) <- ".estimateR"
    
    return (out)
}


#' @title Federated RCCA
#' @description Perform the regularized canonical correlation analysis for the
#' virtual cohort.
#' @usage federateRCCA(loginFD,
#'                     logins,
#'                     func,
#'                     symbol,
#'                     ncomp = 2,
#'                     lambda1 = 0,
#'                     lambda2 = 0,
#'                     chunk = 500L,
#'                     mc.cores = 1,
#'                     tune = FALSE,
#'                     tune_param = .encode.arg(
#'                         list(nfold = 5,
#'                              grid1 = seq(0.001, 1, length = 5),
#'                              grid2 = seq(0.001, 1, length = 5))))
#' @param loginFD Login information of the FD server (one of the servers
#' containing cohort data).
#' @param logins Login information of servers containing cohort data.
#' @param func Encoded definition of a function for preparation of raw data
#' matrices. Two arguments are required: conns (list of Opal connections), 
#' symbol (names of the two R symbols) (see datashield.assign).
#' @param symbol Encoded vector of names of the two R symbols to assign in the
#' DataSHIELD R session on each server in \code{logins}.
#' The two assigned R variables will be used as the input raw data to compute
#' covariance matrices for CCA. Other assigned R variables in \code{func} are
#' ignored.
#' @param ncomp Number of components (covariates). Default, 2.
#' @param lambda1 Non-negative regularized parameter value for first data set.
#' Default, 0. If there are more variables than samples, it should be > 0.
#' @param lambda2 Non-negative regularized parameter value for second data set.
#' Default, 0. If there are more variables than samples, it should be > 0.
#' @param chunk Size of chunks into what the SSCP matrix is partitioned.
#' Default, 500L.
#' @param mc.cores Number of cores for parallel computing. Default, 1.
#' @param tune Logical value indicating whether the tuning for lambda values
#' will be performed. Default, FALSE, no tuning.
#' @param tune_param Tuning parameters.
#' \code{nfold} n-fold cross-validation.
#' \code{grid1} tuning values for \code{lambda1}.
#' \code{grid2} tuning values for \code{lambda2}.
#' @returns RCCA object
#' @importFrom fda geigen
#' @importFrom stats cov
#' @importFrom DSI datashield.assign datashield.aggregate
#' @importFrom parallel mclapply
#' @importFrom arrow write_to_raw
#' @examples
#' \dontrun{
#' dataProc <- function(conns, symbol) {
#'     datashield.assign(conns, symbol[1], 'test.CNSIM', variables=c('LAB_TSC', 'LAB_TRIG', 'LAB_HDL'))
#'     datashield.assign(conns, symbol[2], 'test.CNSIM', variables=c('LAB_GLUC_ADJUSTED', 'PM_BMI_CONTINUOUS'))
#' }
#' federateRCCA(.encode.arg(loginFD), .encode.arg(logins), 
#'              .encode.arg(dataProc, serialize.it = T), 
#'              .encode.arg(c("rawDataX", "rawDataY")))
#' }
#' @export
federateRCCA <- function(loginFD, logins, func, symbol, ncomp = 2,
                         lambda1 = 0, lambda2 = 0,
                         chunk = 500L, mc.cores = 1,
                         tune = FALSE,
                         tune_param = .encode.arg(list(
                             nfold = 5,
                             grid1 = seq(0.001, 1, length = 5), 
                             grid2 = seq(0.001, 1, length = 5)))) {
    .printTime("federateRCCA started")
    if (ncomp < 2) {
        print("ncomp should be at least 2. ncomp will be set to 2.")
        ncomp <- 2
    }
    funcPreProc <- .decode.arg(func)
    querytables <- .decode.arg(symbol)
    if (length(querytables) != 2) {
        stop("Two data matrices are required!")
    }

    ## estimate the parameters of regularization
    if (isTRUE(tune)) {
        tune_param <- .decode.arg(tune_param)
        tuneres <- .estimateR(loginFD, logins, funcPreProc, querytables,
                              chunk=chunk, mc.cores=mc.cores,
                              nfold=tune_param$nfold,
                              grid1=tune_param$grid1,
                              grid2=tune_param$grid2)
        lambda1 <- tuneres$opt.lambda1
        lambda2 <- tuneres$opt.lambda2
    }
    
    ## compute covariance matrices for the virtual cohort
    fedCov <- .federateCov(loginFD, logins, funcPreProc, querytables, pair=T,
                           chunk=chunk, mc.cores=mc.cores, connRes=T)
    
    ## compute rcca
    tryCatch({
        Cxx <- fedCov$cov[[1]]
        Cyy <- fedCov$cov[[2]]
        Cxy <- fedCov$cov[[3]]
        
        ## add parameters of regularization
        Cxx <- Cxx + diag(lambda1, nrow=nrow(Cxx), ncol=ncol(Cxx))
        Cyy <- Cyy + diag(lambda2, nrow=nrow(Cyy), ncol=ncol(Cyy))
        
        ## CCA core call
        rccaObj <- geigen(Cxy, Cxx, Cyy)
        names(rccaObj) <- c("cor", "xcoef", "ycoef")
        rccaObj$xcoef <- rccaObj$xcoef[, 1:ncomp, drop=F]
        rccaObj$ycoef <- rccaObj$ycoef[, 1:ncomp, drop=F]
        rownames(rccaObj$xcoef) <- rownames(Cxx)
        rownames(rccaObj$ycoef) <- rownames(Cyy)
        colnames(rccaObj$xcoef) <- colnames(rccaObj$ycoef) <- 
            paste0("Comp.", 1:ncomp)
        rccaObj$lambda <- list(lambda1=lambda1, lambda2=lambda2)
        
        ## rcca coefs
        loadings <- mclapply(
            rccaObj[c("xcoef", "ycoef")],
            mc.cores=mc.cores,
            function(ccacoef) {
                xx <- t(ccacoef)
                nblockscol <- ceiling(ncol(xx)/chunk)
                sepblockscol <- rep(ceiling(ncol(xx)/nblockscol), nblockscol-1)
                sepblockscol <- c(sepblockscol, ncol(xx) - sum(sepblockscol))
                tcpblocks <- .partitionMatrix(xx,
                                              seprow=ncomp,
                                              sepcol=sepblockscol)
                return (lapply(tcpblocks, function(tcpb) {
                    return (lapply(tcpb, function(tcp) {
                        return (.encode.arg(write_to_raw(tcp)))
                    }))
                }))
            })
        names(loadings) <- querytables
        
        opals <- fedCov$conns
        ## send coefs back to non-FD servers
        pushToDscMate(conns=opals, object=loadings, sourcename='FD', async=T)
        
        ## compute X*coefs'
        datashield.assign(opals, "scores", 
                          as.call(list(as.symbol("tcrossProd"),
                                       x=as.symbol("centeredData"),
                                       y=as.symbol("pushed_FD"),
                                       chunk=chunk)),
                          async=T)
        
        command <- list(as.symbol("pushToDscFD"),
                        as.symbol("FD"),
                        as.symbol("scores"),
                        async=T)
        cat("Command: pushToDscFD(FD, scores)", "\n")
        scoresDSC <- datashield.aggregate(opals, as.call(command), async=T)
        .printTime(paste0("scores communicated to FD: "))
        
        scoresLoc <- lapply(scoresDSC, function(dscblocks) {
            cps <- lapply(dscblocks, function(dscblocki) {
                cpsi <- .rebuildMatrixDsc(dscblocki, mc.cores=mc.cores)
                colnames(cpsi) <- paste0("Comp.", 1:ncomp)
                return (cpsi)
            })
            names(cps) <- c("xcoef", "ycoef")
            return (cps)
        })
        cvx <- do.call(rbind, lapply(scoresLoc, function(sl) sl$xcoef))
        cvy <- do.call(rbind, lapply(scoresLoc, function(sl) sl$ycoef))
        ## TODO: would rownames be necessary?
    }, error = function(e) {
        .printTime(paste0("COEF PUSH PROCESS: ", e))
        return (paste0("COEF PUSH PROCESS: ", e))
    }, finally = {
        datashield.assign(opals, 'crossEnd',
                          as.symbol("crossLogout(FD)"), async=T)
        datashield.logout(opals)
    })
    
    ## correlation between raw data and canonical covariates
    ## formula: cor(a,b) = diag(1/sqrt(diag(cov(a)))) %*% cov(a,b) %*% diag(1/sqrt(diag(cov(b))))
    invdiagcovx <- diag(1/sqrt(diag(Cxx)), nrow=nrow(Cxx), ncol=ncol(Cxx))
    invdiagcovy <- diag(1/sqrt(diag(Cyy)), nrow=nrow(Cyy), ncol=ncol(Cyy))
    invdiagcovcvx <- diag(1/sqrt(diag(cov(cvx))), nrow=ncomp, ncol=ncomp)
    invdiagcovcvy <- diag(1/sqrt(diag(cov(cvy))), nrow=ncomp, ncol=ncomp)
    
    # xxscores <- Reduce('+', lapply(names(opals), function(opn) {
    #     xx <- datashield.aggregate(opals[opn], as.call(list(as.symbol("loadings"),
    #                                                         as.symbol("centeredDatax"),
    #                                                         .encode.arg(cvx[[opn]]), ## cvx obtained without rbind
    #                                                         "crossprod")), 
    #                                async=T)
    #     return (xx[[1]])
    # }))
    #xxscores <- invdiagcovx %*% Cxx %*% rccaObj$xcoef %*% invdiagcovcvx
    xxscores <- crossprod(t(crossprod(invdiagcovx, Cxx)),
                          tcrossprod(rccaObj$xcoef, invdiagcovcvx))
    
    # yxscores <- Reduce('+', lapply(names(opals), function(opn) {
    #     yx <- datashield.aggregate(opals[opn], as.call(list(as.symbol("loadings"),
    #                                                         as.symbol("centeredDatay"),
    #                                                         .encode.arg(cvx[[opn]]),
    #                                                         "crossprod")),
    #                                async=T)
    #     return (yx[[1]])
    # }))
    #yxscores <- invdiagcovy %*% t(Cxy) %*% rccaObj$xcoef %*% invdiagcovcvx
    yxscores <- crossprod(t(tcrossprod(invdiagcovy, Cxy)),
                          tcrossprod(rccaObj$xcoef, invdiagcovcvx))
    
    # xyscores <- Reduce('+', lapply(names(opals), function(opn) {
    #     xy <- datashield.aggregate(opals[opn], as.call(list(as.symbol("loadings"),
    #                                                         as.symbol("centeredDatax"),
    #                                                         .encode.arg(cvy[[opn]]), ## cvy obtained without rbind
    #                                                         "crossprod")),
    #                                async=T)
    #     return (xy[[1]])
    # }))
    #xyscores <- invdiagcovx %*% Cxy %*% rccaObj$ycoef %*% invdiagcovcvy
    xyscores <- crossprod(t(crossprod(invdiagcovx, Cxy)),
                          tcrossprod(rccaObj$ycoef, invdiagcovcvy))
    
    # yyscores <- Reduce('+', lapply(names(opals), function(opn) {
    #     yy <- datashield.aggregate(opals[opn], as.call(list(as.symbol("loadings"),
    #                                                         as.symbol("centeredDatay"),
    #                                                         .encode.arg(cvy[[opn]]),
    #                                                         "crossprod")),
    #                                async=T)
    #     return (yy[[1]])
    # }))
    #yyscores <- invdiagcovy %*% Cyy %*% rccaObj$ycoef %*% invdiagcovcvy
    yyscores <- crossprod(t(crossprod(invdiagcovy, Cyy)),
                          tcrossprod(rccaObj$ycoef, invdiagcovcvy))
    
    rccaObj$scores <- list(xscores=cvx[, 1:ncomp],
                           yscores=cvy[, 1:ncomp],
                           corr.X.xscores=xxscores[, 1:ncomp],
                           corr.Y.xscores=yxscores[, 1:ncomp],
                           corr.X.yscores=xyscores[, 1:ncomp],
                           corr.Y.yscores=yyscores[, 1:ncomp])
    
    class(rccaObj) <- "federateRCCA"
    return (rccaObj)
}


#' @title Map data values to color codes
#' @description Produce color codes for plotting.
#' @usage mapColor(x,
#'                 range.min = NA,
#'                 range.max = NA,
#'                 levels = NA,
#'                 nbreaks = 10,
#'                 colors = .encode.arg(c('orange', 'blue')),
#'                 ...)
#' @param x A factor or a numeric vector.
#' @param range.min Global minimum of x, including other nodes.
#' Default, \code{min(x)}.
#' @param range.max Global maximum of x, including other nodes.
#' Default, \code{max(x)}.
#' @param levels Encoded value of a character vector indicating the global
#' levels of x, when x is a factor. Default, \code{.encode.arg(levels(x))}.
#' @param nbreaks An integer indicating the number of intervals into which x is
#' to be cut, less than \code{length(x)/10}, when x is numeric.
#' @param colors Encoded value of a vector of colors to interpolate, must be a
#' valid argument to col2rgb(). Default, \code{.encode.arg(c('orange', 'blue'))}.
#' @param ... arguments to pass to \code{colorRampPalette}
#' @returns Color codes
#' @importFrom grDevices colorRampPalette
#' @export
mapColor <- function(x, range.min = NA, range.max = NA, levels = NA,
                     nbreaks = 10,
                     colors = .encode.arg(c('orange', 'blue')), ...) {
    rbPal <- colorRampPalette(.decode.arg(colors), ...)
    if (is.factor(x)) {
        glevels <- union(.decode.arg(levels), levels(x))
        levels(x) <- glevels
        if (length(x) < 10*nlevels(x)) {
            stop("x should be longer than 10 times nlevels(x)")
        }
        colbreaks <- rbPal(nlevels(x))[x]
    } else if (is.numeric(x)) {
        if (length(x) < 10*nbreaks) {
            stop("x should be longer than 10 times nbreaks")
        }
        if (is.na(range.min)) range.min <- range(x)[1]
        if (is.na(range.max)) range.max <- range(x)[2]
        colbreaks <- as.vector(cut(c(range.min, x, range.max), nbreaks,
                                   labels=rbPal(nbreaks)))
        colbreaks <- colbreaks[-c(1, length(colbreaks))]
    } else {
        stop("A numeric vector or a factor is required.")
    }
    
    return (colbreaks)
}
