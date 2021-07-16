#' @title Ranking
#'
#' Ranking of features in each sample 
#' @param x A matrix or data frame, samples in rows and features in columns
#' @return Ranking of features
#' @export
dsRank <- function(x) {
    return (t(apply(x, 1, rank)))
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
computeWeights <- function(expressionData, indexMatrix, top, bottom, topWeight, bottomWeight){
    
    outWeight = lapply(names(opals), function(x){
    #x = names(opals)[1]
    
    indexMatrix_weight= matrix(1, nrow(indexMatrix[[x]]), ncol(indexMatrix[[x]]))
    indexMatrix_weight[indexMatrix[[x]]<=bottom] <- bottomWeight 
    indexMatrix_weight[indexMatrix[[x]]>=top] <- topWeight
    
    WeightexpressionData = expressionData[[x]] * indexMatrix_weight
    
    })
    
    return(outWeight)
    
 }




