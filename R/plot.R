
#' Plotting stratifying clusters
#' 
#' @inheritParams plot_error_rate
#' @param clusters.data A \code{data.frame} containing the clusters data. Each row is a cell
#'   and column represents the expression of different markers (non numeric columns are dropped)
#' @param col.names Which columns of \code{clustersData} to include in the plot. The default is to include all
#'   columns in the plot
#' @param by.cluster If this is \code{TRUE} this function will generate a separate plot 
#'   for each cluster, otherwise or a single plot with all the clusters. 
#'   The latter option may take quite long to plot if there are several clusters   
#' @export    
plot_stratifying_clusters <- function(citrus.model, clusters.data, output.dir, col.names = names(clusters.data), 
                          by.cluster = FALSE) {
    differential.features <- citrus.model$differentialFeatures

    for (cv.point in names(differential.features)) {
        cluster.ids <- as.numeric(differential.features[[cv.point]][["clusters"]])
        out.dir <- file.path(output.dir, cv.point)
        dir.create(out.dir, recursive = TRUE, showWarnings = FALSE)
        
        output.file <- NULL
        if(!by.cluster)
            output.file <- file.path(out.dir, "clusters.pdf")
        citrus.plotClusters(cluster.ids, clusters.data, 
                            byCluster = by.cluster, outputDir = out.dir, outputFile = output.file)
    }
}

#' Plotting model error rate
#' 
#' @param citrus.model An object of class \code{citrus.regressionResult}, as returned from
#'   \code{\link{get_model}}
#' @param output.dir The output file path
#' @return Returns \code{NULL}
#' @export
plot_error_rate <- function(citrus.model, output.file) {
    citrus.plotTypeErrorRate(modelType = citrus.model$modelType, outputFile = output.file, 
                                 regularizationThresholds = citrus.model$regularizationThresholds, 
                                 thresholdCVRates = citrus.model$thresholdCVRates, finalModel = citrus.model$finalModel$model, 
                                 cvMinima = citrus.model$cvMinima, family = citrus.model$family)

    return(invisible(NULL))
}

#' @export
plot_stratifying_features <- function(citrus.model, output.dir, by.cluster = FALSE, all.features = FALSE) {
    do.call(paste("citrus.plotModelDifferentialFeatures", citrus.model$family, 
                  sep = "."), args = list(differentialFeatures = citrus.model$differentialFeatures, 
                                          features = citrus.model$allFeatures, modelOutputDirectory = output.dir, 
                                          labels = citrus.model$labels, byCluster = by.cluster, allFeatures = all.features))
}


