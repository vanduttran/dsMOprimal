#' @title Ranking
#'
#' Ranking of features in each sample 
#' @param x A matrix or data frame, samples in rows and features in columns
#' @return Ranking of features
#' @export
dsRank <- function(x) {
    out = t(apply(x, 1, rank))
    rownames(out) = rownames(x)
    colnames(out) = colnames(x)
    
    return(out)
}


#' @title computeWeights
#'
#' Weight most and least regulated features in original dataMatrix accordingly to ranking 
#' @param expressionData A matrix or data frame, samples in rows and features in columns
#' @param indexMatrix rank matrix of same dimension of expressionData matrix
#' @param top which most regulated genes we consider
#' @param bottom which least regulated genes we consider
#' @param topWeight weight associated to top elements
#' @param bottomWeight weight associated to bottom elements
#' @return outWeiight matrix of the same dimension of expressionData
#' @export
computeWeights <- function(expressionData, indexMatrix, top = 2, bottom= 2, topWeight = 10, bottomWeight= 0.1){
     
    indexMatrix_weight= matrix(1, nrow(indexMatrix), ncol(indexMatrix)) 

    indexMatrix_weight[apply(indexMatrix, 1, as.numeric)<=as.numeric(bottom)] <- bottomWeight 
    indexMatrix_weight[apply(indexMatrix, 1, as.numeric)>=as.numeric(top)] <- topWeight
    
    # I'm transposing original WeightexpressionData to give it the same dimension of expressionData
    
    WeightexpressionData = t(apply(expressionData, 1, as.numeric) * apply(indexMatrix_weight, 1, as.numeric))   
    
    rownames(WeightexpressionData) = as.vector(rownames(expressionData))
    colnames(WeightexpressionData) = as.vector(colnames(expressionData))
    
    print("......")
    print(dim(expressionData)); print(dim(WeightexpressionData));
    return(t(apply(WeightexpressionData,1, as.numeric)))
    
 }



#' @title aggRownames
#'
#' Weight most and least regulated features in original dataMatrix accordingly to ranking 
#' @param expressionData A matrix or data frame, samples in rows and features in columns
#' @return rownames of expressonData
#' @export
aggRownames <- function(expressionData){
    
    expD <- dsSwissKnife:::.decode.arg(expressionData)
    return(rownames(expD))

}


#' @title Federate SSCP on weighted data
#' @description Function for computing the federated SSCP matrix
#' @param loginFD Login information of the FD server
#' @param logins Login information of data repositories
#' @param querytab Encoded name of a table reference in data repositories
#' @param queryvar Encoded variables from the table reference
#' @param TOL Tolerance of 0
#' @import DSOpal parallel bigmemory
#' @keywords internal
federateSSCPweight <- function(loginFD, logins, querytab, queryvar, TOL = 1e-10) {
    require(DSOpal)

    loginFDdata    <- dsSwissKnife:::.decode.arg(loginFD)
    logindata      <- dsSwissKnife:::.decode.arg(logins)
    querytable     <- dsSwissKnife:::.decode.arg(querytab)
    queryvariables <- dsSwissKnife:::.decode.arg(queryvar)

    opals <- DSI::datashield.login(logins=logindata)
    nNode <- length(opals)

    datashield.assign(opals, "rawData", querytable, variables=queryvariables, async=T)
    datashield.assign(opals, "indexMatrix", as.symbol('dsRank(cent)'), async=T)
    datashield.assign(opals, "weightMatrix",as.symbol("computeWeights(rawData, indexMatrix)"), async = T)
    datashield.assign(opals, "centeredData", as.symbol('center(computeWeights)'), async=T)
    datashield.assign(opals, "crossProdSelf", as.symbol('crossProd(centeredData)'), async=T)
    datashield.assign(opals, "tcrossProdSelf", as.symbol('tcrossProd(centeredData, chunk=50)'), async=T)

    samples <- datashield.aggregate(opalborder, as.symbol('aggRownames(rawData)'), async=T)
    

       ##- received by each from other nodes ----
    invisible(mclapply(names(opals), mc.cores=1, function(opn) {
        logindata.opn <- logindata[logindata$server != opn, , drop=F]
        logindata.opn$user <- logindata.opn$userserver
        logindata.opn$password <- logindata.opn$passwordserver
        opals.loc <- paste0("crossLogin('", .encode.arg(logindata.opn), "')")
        datashield.assign(opals[opn], 'mates', as.symbol(opals.loc), async=F)

        command.opn <- list(paste0("crossAssign(mates, symbol='rawDataMate', value='",
                                   querytab,
                                   "', value.call=F, variables='",
                                   queryvar,
                                   "', async=F)"),
                            paste0("crossAssign(mates, symbol='centeredDataMate', value='",
                                   .encode.arg("center(rawDataMate)"),
                                   "', value.call=T, async=F)")
        )
        for (command in command.opn) {
            cat("Command: ", command, "\n")
            print(datashield.aggregate(opals[opn], as.symbol(command), async=F))
        }

        command.opn <- paste0("crossAggregate(mates, '", .encode.arg('singularProd(centeredDataMate)'), "', async=F)")
        cat("Command: ", command.opn, "\n")
        print(datashield.assign(opals[opn], "singularProdMate", as.symbol(command.opn), async=F))

        command.opn <- paste0("crossAggregate(mates, '",
                              .encode.arg(paste0("as.call(list(as.symbol('pushValue'), dsSSCP:::.encode.arg(crossProdSelf), dsSSCP:::.encode.arg('", opn, "')))")),
                              "', async=F)")
        cat("Command: ", command.opn, "\n")
        print(datashield.assign(opals[opn], "pidMate", as.symbol(command.opn), async=F))
    }))
    datashield.symbols(opals)
      


     #-----

    ## (X_i) * (X_i)': push this symmetric matrix from each node to FD
    #crossProdSelf     <- datashield.aggregate(opals, as.symbol('tcrossProd(centeredData)'), async=T)
    datashield.assign(opals, 'FD', as.symbol(paste0("crossLogin('", loginFD, "')")), async=T)

    # command <- paste0("crossAggregate(FD, '", 
    #                   .encode.arg(paste0("as.call(list(as.symbol('garbageCollect')", "))")), 
    #                   "', async=T)")
    # cat("Command: ", command, "\n")
    # datashield.assign(opals, "GC", as.symbol(command), async=T)

    command <- paste0("dscPush(FD, '",
                      .encode.arg(paste0("as.call(list(as.symbol('pushSymmMatrix'), dsSSCP:::.encode.arg(tcrossProdSelf)", "))")),
                      "', async=T)")
    cat("Command: ", command, "\n")
    crossProdSelfDSC <- datashield.aggregate(opals, as.symbol(command), async=T)



    crossProdSelf <- mclapply(crossProdSelfDSC, mc.cores=min(length(opals), detectCores()), function(dscblocks) {
        return (as.matrix(attach.big.matrix(dscblocks[[1]])))
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

    ## (X_i) * (X_j)' * ((X_j) * (X_j)')[,1]: push this single-column matrix from each node to FD
    #singularProdCross <- datashield.aggregate(opals, as.symbol('tcrossProd(centeredData, singularProdMate)'), async=T)
    datashield.assign(opals, "singularProdCross", as.symbol('tcrossProd(centeredData, singularProdMate)'), async=T)



    command <- paste0("dscPush(FD, '",
                      .encode.arg(paste0("as.call(list(as.symbol('pushSingMatrix'), dsSSCP:::.encode.arg(singularProdCross)", "))")),
                      "', async=T)")
    cat("Command: ", command, "\n")
    singularProdCrossDSC <- datashield.aggregate(opals, as.symbol(command), async=T)

    singularProdCross <- mclapply(singularProdCrossDSC, mc.cores=length(singularProdCrossDSC), function(dscbigmatrix) {
        dscMatList <- lapply(dscbigmatrix[[1]], function(dsc) {
            dscMat <- matrix(as.matrix(attach.big.matrix(dsc)), ncol=1) #TOCHECK: with more than 2 servers
            stopifnot(ncol(dscMat)==1)
            return (dscMat)
        })
        return (dscMatList)
    })
    gc(reset=F)
    #return(singularProdCross)
    ##  (X_i) * (X_j)' * (X_j) * (X_i)'
    #prodDataCross     <- datashield.aggregate(opals, as.symbol('tripleProd(centeredData, crossProdMate)'), async=F)
    ## N.B. save-load increase numeric imprecision!!!
    prodDataCross     <- datashield.aggregate(opals, as.call(list(as.symbol("tripleProd"),
                                                                  as.symbol("centeredData"),
                                                                  .encode.arg(names(opals)))), async=T)

        ## deduced from received info by federation
    crossProductPair <- lapply(1:(nNode-1), function(opi) {
        crossi <- lapply((opi+1):(nNode), function(opj) {
            opni <- names(opals)[opi]
            opnj <- names(opals)[opj]

            a1 <- solveSSCP(XXt=prodDataCross[[opni]][[opnj]],
                            XtX=prodDataCross[[opnj]][[opni]],
                            r=crossProdSelf[[opnj]][, 1, drop=F],
                            Xr=singularProdCross[[opni]][[opnj]],
                            TOL=TOL)
            a2 <- solveSSCP(XXt=prodDataCross[[opnj]][[opni]],
                            XtX=prodDataCross[[opni]][[opnj]],
                            r=crossProdSelf[[opni]][, 1, drop=F],
                            Xr=singularProdCross[[opnj]][[opni]],
                            TOL=TOL)
            cat("Precision on a1 = t(a2):", max(abs(a1 - t(a2))), "\n")
            return (a1)
        })
        names(crossi) <- names(opals)[(opi+1):(nNode)]
        return (crossi)
    })
    names(crossProductPair) <- names(opals)[1:(nNode-1)]
    

     ## SSCP whole matrix
    XXt <- do.call(rbind, lapply(1:nNode, function(opi) {
        upper.opi <- do.call(cbind, as.list(crossProductPair[[names(opals)[opi]]]))
        lower.opi <- do.call(cbind, lapply(setdiff(1:opi, opi), function(opj) {
            t(crossProductPair[[names(opals)[opj]]][[names(opals)[opi]]])
        }))
        return (cbind(lower.opi, crossProdSelf[[opi]], upper.opi))
    }))
    datashield.logout(opals)
    
    rownames(XXt) = samples
    colnames(XXt) = samples

    return (XXt)
}

