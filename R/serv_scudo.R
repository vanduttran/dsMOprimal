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
computeWeights <- function(expressionData, indexMatrix, top = 2, bottom= 2, topWeight = 10, bottomWeight= 0.1){
    
    indexMatrix_weight= matrix(1, nrow(indexMatrix), ncol(indexMatrix))
    indexMatrix_weight[indexMatrix<=bottom] <- bottomWeight 
    indexMatrix_weight[indexMatrix>=top] <- topWeight

    WeightexpressionData = expressionData * indexMatrix_weight 
    return(WeightexpressionData)
    
 }




