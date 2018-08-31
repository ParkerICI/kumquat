
#' Select clusters for further analysis
#' 
#' Selects clusters from a set of clusters identified by a clustering method investigate for stratifying signal.
#' 
#' @param citrus.clustering A \code{citrus.clustering} object.
#' @param method Method for determining which clusters from a clustering should be analyzed for stratifying signal.
#' @param ... Other parameters passed to specific cluster selection methods.
#' 
#' @return A vector of cluster IDs. 
#' 
#' @author Robert Bruggner
#' @export
citrus.selectClusters <- function(citrus.clustering, method = "minimumClusterSize", ...) {
    do.call(paste0("citrus.selectClusters.", method), args = list(citrus.clustering = citrus.clustering, ...))
}

#' Selects clusters for endpoint analysis
#' 
#' Selects clusters for endpoint analysis by cluster size. Selected clusters must have a minimum number of 
#' cells in them as a proportion of the total number of clustered events. If \eqn{n} total events are clustered,
#' clusters contatining at least \eqn{n * minimumClusterSizePercent} events are selected.
#' 
#' @param citrus.clustering A \code{citrus.clustering} object.
#' @param minimumClusterSizePercent The percentage \eqn{0 < x < 1} of the total number of clustered events a cluster 
#' must contain in order to be selected. 
#' @param ... Other arguments (ignored).
#' 
#' @author Robert Bruggner
#' @export
citrus.selectClusters.minimumClusterSize <- function(citrus.clustering, minimumClusterSizePercent = 0.05, ...) {
    clusterSizes <- sapply(citrus.clustering$clusterMembership, length)
    minimumClusterSize <- (length(citrus.clustering$clusterMembership) + 1) * minimumClusterSizePercent
    return(which(clusterSizes >= minimumClusterSize))
}


citrus.convertConditionMatrix <- function(conditionMatrix) {
    conditions <- list()
    for (i in 1:nrow(conditionMatrix)) {
        for (j in 1:ncol(conditionMatrix)) {
            if (conditionMatrix[i, j]) {
                conditions <- append(conditions, list(unique(c(rownames(conditionMatrix)[i], colnames(conditionMatrix)[j]))))
            }
        }
    }
    return(conditions)
}

#' List possible model types
#' 
#' Returns valid model types that Citrus can compute.
#' 
#' @return Vector of model types that Citrus can compute.
#' 
#' @author Robert Bruggner
#' @export
citrus.modelTypes <- function() {
    return(c("pamr", "glmnet", "sam"))
}


