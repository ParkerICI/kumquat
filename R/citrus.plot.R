
citrus.plotTypeErrorRate <- function(modelType, outputFile, regularizationThresholds, 
    thresholdCVRates, finalModel, cvMinima, family) {
    if (modelType == "sam") {
        message("Error rate plots are not available for sam models")
        return(invisible(NULL))
    }
    
    pdf(outputFile, width = 6, height = 6)
    thresholds <- regularizationThresholds
    errorRates <- thresholdCVRates[, "cvm"]
    ylim <- c(0, 1)
    if (family == "continuous") {
        ylim <- c(min(thresholdCVRates[, "cvm"] - thresholdCVRates[, "cvsd"]) * 0.9, 
            max(thresholdCVRates[, "cvm"] + thresholdCVRates[, "cvsd"]) * 1.1)
    }
    ylab <- "Model Cross Validation Error Rate"
    if (modelType == "glmnet") {
        thresholds <- log(thresholds)
        xlab <- "log(Regularization Threshold)"
        nonzeroCounts <- finalModel$df
        if (family == "survival") {
            ylim <- range(errorRates, na.rm = T)
            ylab <- "Model Cross Validation Partial Likelihood"
        }
    } else if (modelType == "pamr") {
        xlab <- "Regularization Threshold"
        nonzeroCounts <- finalModel$nonzero
    }
    plot(errorRates, type = "o", pch = 20, col = "red", main = "Number of model features\n", 
        axes = F, xlab = xlab, ylim = ylim, ylab = ylab)
    # Plot SEM
    for (i in 1:length(thresholds)) {
        lines(c(i, i), c(errorRates[i] + thresholdCVRates[, "cvsd"][i], errorRates[i] - 
            thresholdCVRates[, "cvsd"][i]), col = "red", lty = 3)
    }
    grid()
    axis(1, at = 1:length(errorRates), labels = sapply(thresholds, sprintf, fmt = "%1.2f"))
    if (family == "survival") {
        axis(2)
    } else {
        axis(2, at = c(0, 0.25, 0.5, 0.75, 1), labels = c(0, 25, 50, 75, 100))
    }
    
    axis(3, labels = nonzeroCounts, at = 1:length(errorRates))
    legendLabels <- c("Cross Validation Error Rate")
    legendColors <- c("red")
    legendPchs <- c(20)
    legendLty <- c(1)
    legendPtCex <- c(1)
    if (!is.null(thresholdCVRates$fdr)) {
        lines(thresholdCVRates$fdr, type = "o", pch = 1, cex = 1.5, col = "blue")
        legendLabels <- c(legendLabels, "Feature False Discovery Rate")
        legendColors <- c(legendColors, "blue")
        legendPchs <- c(legendPchs, 1)
        legendLty <- c(legendLty, 1)
        legendPtCex <- c(legendPtCex, 1)
    }
    
    cv.min <- cvMinima$cv.min.index
    cv.1se <- cvMinima$cv.1se.index
    if (!is.null(cv.min)) {
        points(c(cv.min, cv.min), y = c(errorRates[cv.min], errorRates[cv.min]), 
            col = "green", pch = 20, cex = 2)
    }
    if (!is.null(cv.1se)) {
        points(c(cv.1se, cv.1se), y = c(errorRates[cv.1se], errorRates[cv.1se]), 
            col = "orange", pch = 9, cex = 2)
    }
    
    legendLabels <- c(legendLabels, "cv.min", "cv.1se")
    legendColors <- c(legendColors, "green", "orange")
    legendPchs <- c(legendPchs, 20, 9)
    legendLty <- c(legendLty, 0, 0)
    legendPtCex <- c(legendPtCex, 2, 1.5)
    
    if ("cv.fdr.constrained" %in% names(cvMinima)) {
        
        cv.fdr.constrained <- cvMinima$cv.fdr.constrained.index
        points(c(cv.fdr.constrained, cv.fdr.constrained), y = c(errorRates[cv.fdr.constrained], 
            errorRates[cv.fdr.constrained]), col = "yellow", pch = 17, cex = 1.5)
        points(c(cv.fdr.constrained, cv.fdr.constrained), y = c(errorRates[cv.fdr.constrained], 
            errorRates[cv.fdr.constrained]), col = "black", pch = 2, cex = 1.5)
        legendLabels <- c(legendLabels, "cv.fdr.constrained")
        legendColors <- c(legendColors, "yellow")
        legendPchs <- c(legendPchs, 17)
        legendLty <- c(legendLty, 0)
        legendPtCex <- c(legendPtCex, 1.5)
        
    }
    legend(x = "topleft", legendLabels, col = legendColors, pch = legendPchs, lty = legendLty, 
        pt.cex = legendPtCex, cex = 0.8, bg = "white")
    dev.off()
}


getFeaturesPlot.classification <- function(melted) {
    p <- (ggplot2::ggplot(melted, ggplot2::aes(x = factor(labels), y = value)) 
          + ggplot2::facet_wrap(~variable, ncol = 1) 
          + ggplot2::geom_boxplot(outlier.colour = rgb(0, 0, 0, 0), colour = rgb(0, 0, 0, 0.3)) 
          + ggplot2::geom_point(ggplot2::aes(color = factor(labels)), alpha = I(0.25), shape = 19, size = I(2)) 
          + ggplot2::coord_flip() 
          + ggplot2::theme_bw() 
          + ggplot2::ylab("") 
          + ggplot2::xlab("") 
          + ggplot2::theme(legend.position = "none")
    )
    #if (any(grepl(pattern = "abundance", nonzeroFeatureNames))) {
    #    p <- p + ggplot2::scale_y_log10() + ggplot2::ylab("Log10 scale")
    #}
    return(p)
}

citrus.plotModelDifferentialFeatures.classification <- function(differentialFeatures, 
    features, modelOutputDirectory, labels, byCluster = FALSE, allFeatures = FALSE) {
    
    cvPoints <- names(differentialFeatures)
    if(allFeatures)
        cvPoints <- c(cvPoints, "allFeatures")
 
    for (cvPoint in cvPoints) {
        if(cvPoint == "allFeatures")
            nonzeroFeatureNames <- colnames(features)
        else
            nonzeroFeatureNames <- differentialFeatures[[cvPoint]][["features"]]
        outDir <- file.path(modelOutputDirectory, cvPoint, "features")
        dir.create(outDir, recursive = TRUE, showWarnings = FALSE)

        # Write features to file for easy parsing
        write.csv(features[, nonzeroFeatureNames], file = file.path(outDir, 
            paste("features_", cvPoint, ".csv", sep = "")), quote = F)
     
        melted <- reshape2::melt(data.frame(features[, nonzeroFeatureNames, drop = F], 
            labels = labels, check.names = F), id.vars = "labels")
        
        melted$cluster <- sapply(strsplit(as.character(melted$variable), "_"), "[", 2)
        
        p <- getFeaturesPlot.classification(melted)
        
        ggplot2::ggsave(file.path(outDir, paste("features_", cvPoint, ".pdf", sep = "")), 
                        plot = p, width = 4, height = length(nonzeroFeatureNames) * 1.5, limitsize = FALSE)
        
        if(byCluster) {
            outDir <- file.path(modelOutputDirectory, cvPoint, "clusters")
            dir.create(outDir, recursive = TRUE, showWarnings = FALSE)
            plyr::d_ply(melted, ~cluster, function(x) {
                p <- getFeaturesPlot.classification(x)
                ggplot2::ggsave(file.path(outDir, paste("c", x$cluster[1], ".pdf", sep = "")), plot = p)
            })
        }
        

    }
}


# This one still needs refactoring
citrus.plotModelDifferentialFeatures.continuous <- function(differentialFeatures, 
    features, modelOutputDirectory, labels, ...) {
    for (cvPoint in names(differentialFeatures)) {
        nonzeroFeatureNames <- differentialFeatures[[cvPoint]][["features"]]
        
        # Write features to file for easy parsing
        write.csv(features[, nonzeroFeatureNames], file = file.path(modelOutputDirectory, 
            paste("features_", cvPoint, ".csv", sep = "")), quote = F)
        
        melted <- reshape2::melt(data.frame(features[, nonzeroFeatureNames, drop = F], 
            labels = labels, check.names = F), id.vars = "labels")
        
        pdf(file.path(modelOutputDirectory, paste("features_", cvPoint, ".pdf", sep = "")), 
            width = 4, height = length(nonzeroFeatureNames) * 1.5)
        p <- (ggplot2::ggplot(melted, ggplot2::aes(x = value, y = labels)) 
                + ggplot2::facet_wrap(~variable, ncol = 1) 
                + ggplot2::geom_point(size = I(2)) 
                + ggplot2::theme_bw() 
                + ggplot2::ylab("") 
                + ggplot2::xlab("") 
                + ggplot2::theme(legend.position = "none")
        )
        if (any(grepl(pattern = "abundance", nonzeroFeatureNames))) {
            p <- p + ggplot2::scale_x_log10() + ggplot2::xlab("Log10 scale")
        }
        print(p)
        dev.off()
    }
    
}

citrus.densityPlotGrid <- function(clustersData, backgroundData) {
    clustersData <- reshape2::melt(clustersData, id.vars = "cellType")
    clustersData <- data.frame(clustersData, src = "Cluster", check.names = FALSE, stringsAsFactors = FALSE)
    
    backgroundData <- reshape2::melt(backgroundData)
    backgroundData <- data.frame(backgroundData, src = "Background", check.names = FALSE, stringsAsFactors = FALSE)
    
    p <- (ggplot2::ggplot(data = clustersData, ggplot2::aes(x = value, y = ..scaled.., fill = src)) 
            + ggplot2::geom_density() 
            + ggplot2::facet_grid(cellType ~ variable, scales = "free") 
            + ggplot2::geom_density(data = backgroundData) 
            + ggplot2::theme_bw() 
            + ggplot2::scale_fill_manual(values = c(Background = rgb(0.3, 0.3, 1, 0.2), Cluster = rgb(1, 0.3, 0.3, 0.5))) 
            + ggplot2::theme(legend.position = "bottom", axis.text.y = ggplot2::element_blank(), 
                axis.ticks.y = ggplot2::element_blank(), axis.title = ggplot2::element_blank()) 
            + ggplot2::labs(fill = "Distribution:")
        )

    return(p)
}

citrus.clusterDensityPlots <- function(clustersData, backgroundData) {
    clustersData$type <- "Cluster"
    backgroundData$type <- "Background"

    tab <- rbind(clustersData, backgroundData)
    tab$cellType <- NULL
    m <- reshape2::melt(tab, id.vars = "type")
    
    p <- (ggplot2::ggplot(data = m, ggplot2::aes(x = value, fill = type, y = ..scaled..)) 
          + ggplot2::geom_density() 
          + ggplot2::facet_wrap(~ variable, scales = "free") 
          + ggplot2::theme_bw() 
          + ggplot2::scale_fill_manual(values = c(Background = rgb(0.3, 0.3, 1, 0.2), Cluster = rgb(1, 0.3, 0.3, 0.5))) 
    )
    return(p)
    
}

#' Plot cluster histograms
#' 
#' Plot expression of markers in cluster cells relative to all cells
#' 
#' @param clusterIds Vector of cluster IDs to plot
#' @param clustersData \code{data.frame} containing the clusters data. Each row is a cell
#'   and column represents the expression of different markers (non numeric columns are dropped). 
#'   The \code{data.frame} must contain a column called \code{cellType} which indicates the cluster assignments. 
#'   The values in \code{clusterIds} must match the values in \code{cellType}.  
#' @param colNames Which columns of \code{clustersData} to include in the plot. The default is to include all
#'   columns in the plot
#' @param byCluster If this is \code{TRUE} this function will generate a separate plot 
#'   for each cluster, otherwise or a single plot with all the clusters. 
#'   The latter option may take quite long to plot if there are several clusters
#' @param outputDir The directory in which the plots will be written (only used if \code{byCluster == TRUE}).
#' @param outputFile The output file the plot will be written to (only used if \code{byCluster == FALSE}). Note that
#'   this needs to be the full path to the file, as \code{outputDir} is ignored when \code{byCluster == FALSE}. If
#'   this parameter is \code{NULL} and \code{byCluster == FALSE} no file will be written 
#'   (and the plot will be returned invisibly)
#'   
#' @return Either \code{NULL} if \code{byCluster == TRUE} or a \code{ggplot2} plot object otherwise
#' 
citrus.plotClusters <- function(clusterIds, clustersData, colNames = names(clustersData), 
                                byCluster = FALSE, outputDir = NULL, outputFile = NULL) {
    data <- clustersData[, colNames]
    data <- data[, sapply(data, is.numeric)]
    data$cellType <- clustersData$cellType
    
    clusterDataList <- list()
    bgData <- data
    bgData$cellType <- NULL
    
    if (nrow(bgData) > 5000) 
        bgData <- bgData[sample(1:nrow(bgData), 5000), ]
    
    
    for (clusterId in sort(clusterIds)) {
        temp <- data[data$cellType == clusterId, ]

        if (nrow(temp) > 2500) 
            temp <- temp[sample(1:nrow(temp), 2500), ]
      
        if(byCluster) {
            p <- citrus.clusterDensityPlots(temp, bgData)
            plotDim <- (ncol(data) / 4) + 1
            
            ggplot2::ggsave(file.path(outputDir, sprintf("c%d.pdf", clusterId)), plot = p, 
                            width = plotDim, height = plotDim, limitsize = FALSE)
        } else {
            clusterDataList[[as.character(clusterId)]] <- temp
        }
    }

    if(!byCluster) {
        clustersData <- do.call(rbind, clusterDataList)
        p <- citrus.densityPlotGrid(clustersData, bgData)
    
        if (!is.null(outputFile)) 
            ggplot2::ggsave(outputFile, p, width = (2.2 * ncol(data) + 2), height = (2* length(clusterIds)), limitsize = FALSE)
        return(invisible(p))
    
    } else {
        return(invisible(NULL))
    }
    
}



#' Plot a citrus.full.result
#' @name plot.citrus.full.result
#' @param citrus.full.result A citrus.full.result object
#' @param outputDirectory Full path to directory in which to place plot output. 
#' 
#' 
#' 
#' @author Robert Bruggner
plot.citrus.full.result <- function(citrus.full.result, outputDirectory) {
    
    if (!file.exists(outputDirectory)) {
        stop(paste0("Output directory '", outputDirectory, "' does not exist."))
    }
    # Should
    for (conditionName in names(results$conditions)) {
        cat(paste0("\nPlotting Results for ", conditionName, "\n"))
        conditionOutputDir <- file.path(outputDirectory, conditionName)
        dir.create(conditionOutputDir, showWarnings = F)
        parallel::mclapply(citrus.full.result$conditionRegressionResults[[conditionName]], 
            plot, outputDirectory = conditionOutputDir, citrus.foldClustering = citrus.full.result$citrus.foldClustering, 
            citrus.foldFeatureSet = citrus.full.result$conditionFoldFeatures[[conditionName]], 
            citrus.combinedFCSSet = citrus.full.result$citrus.combinedFCSSet, family = citrus.full.result$family, 
            labels = citrus.full.result$labels, conditions = citrus.full.result$conditions[[conditionName]])
        cat("\n")
    }
}



