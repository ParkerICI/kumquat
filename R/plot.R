
    
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


#' @export
plot_error_rate <- function(citrus.model, output.dir) {
  citrus.plotTypeErrorRate(modelType = citrus.model$modelType, modelOutputDirectory = output.dir, 
                                 regularizationThresholds = citrus.model$regularizationThresholds, 
                                 thresholdCVRates = citrus.model$thresholdCVRates, finalModel = citrus.model$finalModel$model, 
                                 cvMinima = citrus.model$cvMinima, family = citrus.model$family)

}

#' @export
plot_stratifying_features <- function(citrus.model, output.dir, by.cluster = FALSE, all.features = FALSE) {
    do.call(paste("citrus.plotModelDifferentialFeatures", citrus.model$family, 
                  sep = "."), args = list(differentialFeatures = citrus.model$differentialFeatures, 
                                          features = citrus.model$allFeatures, modelOutputDirectory = output.dir, 
                                          labels = citrus.model$labels, byCluster = by.cluster, allFeatures = all.features))
}


