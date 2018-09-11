
    
#' @export    
plot_clusters <- function(citrus.model, clusters.data, output.dir, col.names = names(clusters.data), 
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
    
