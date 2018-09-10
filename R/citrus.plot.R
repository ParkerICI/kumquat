######################## Helper Functions
.graphColorPalette <- function(x, alpha = 1) {
    # rainbow(x,alpha=.8,start=.65,end=.15)
    topo.colors(x, alpha = alpha)
}

.formatDecimal <- function(x) {
    sprintf("%1.2f", x)
}

.getDisplayNames <- function(citrus.combinedFCSSet, clusteringColumns) {
    
    # Get channel names colLabels = citrus.combinedFCSSet$fileChannelNames[[1]][[1]]
    
    # Get reagent names reagentNames =
    # citrus.combinedFCSSet$fileReagentNames[[1]][[1]]
    
    # Set display names to channel names displayNames = colLabels
    
    # update display names to equal reagent names where reagent names > 2
    # displayNames[nchar(reagentNames)>2] = reagentNames[nchar(reagentNames)>2]
    
    displayNames <- citrus.combinedFCSSet$fileReagentNames[[1]][[1]]
    
    if (all(is.numeric(clusteringColumns))) {
        # If clustering columns is numeric, pass back the indices
        return(displayNames[clusteringColumns])
    } else {
        # Otherwise, set the names of display names to be the channel names
        # names(displayNames) = colLabels
        names(displayNames) <- citrus.combinedFCSSet$fileChannelNames[[1]][[1]]
        # and return the display names from the corresponding clustering columns
        return(as.vector(displayNames[clusteringColumns]))
    }
    
}

.getClusterMedians <- function(clusterId, clusterAssignments, clusterCols, data) {
    apply(data[clusterAssignments[[clusterId]], clusterCols], 2, median)
}

.decimalFormat <- function(x) {
    sprintf("%.2f", x)
}

.scaleToOne <- function(x) {
    x <- x - min(x)
    x <- x/max(x)
    return(x)
}

.getClusterFeatureMatrix <- function(featureVector) {
    df <- do.call("rbind", strsplit(gsub(pattern = "(cluster [0-9]+) ", replacement = "\\1~", 
        featureVector), "~"))
    return(cbind(cluster = do.call("rbind", strsplit(df[, 1], " "))[, 2], feature = df[, 
        2]))
}

citrus.plotTypeErrorRate <- function(modelType, modelOutputDirectory, regularizationThresholds, 
    thresholdCVRates, finalModel, cvMinima, family) {
    if (modelType == "sam") {
        return()
    }
    
    pdf(file.path(modelOutputDirectory, "ModelErrorRate.pdf"), width = 6, height = 6)
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
    axis(1, at = 1:length(errorRates), labels = sapply(thresholds, .formatDecimal))
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
        
        # Write features to file for easy parsing
        write.csv(features[, nonzeroFeatureNames], file = file.path(modelOutputDirectory, 
            paste("features_", cvPoint, ".csv", sep = "")), quote = F)
        
        melted <- reshape2::melt(data.frame(features[, nonzeroFeatureNames, drop = F], 
            labels = labels, check.names = F), id.vars = "labels")
        
        melted$cluster <- sapply(strsplit(as.character(melted$variable), "_"), "[", 2)
        
        cvPoint <- sub(pattern = "\\.", replacement = "_", x = cvPoint)
        
        
        if(byCluster) {
            outDir <- file.path(modelOutputDirectory, cvPoint)
            dir.create(outDir, recursive = TRUE, showWarnings = FALSE)
            plyr::d_ply(melted, ~cluster, function(x) {
                p <- getFeaturesPlot.classification(x)
                ggplot2::ggsave(file.path(outDir, paste("c", x$cluster[1], ".pdf", sep = "")), plot = p)
            })
        }
        
        p <- getFeaturesPlot.classification(melted)
        
        ggplot2::ggsave(file.path(modelOutputDirectory, paste("features-", cvPoint, ".pdf", sep = "")), 
                        plot = p, width = 4, height = length(nonzeroFeatureNames) * 1.5, limitsize = FALSE)
    }
}

citrus.plotModelDifferentialFeatures.continuous <- function(differentialFeatures, 
    features, modelOutputDirectory, labels, ...) {
    for (cvPoint in names(differentialFeatures)) {
        nonzeroFeatureNames <- differentialFeatures[[cvPoint]][["features"]]
        
        # Write features to file for easy parsing
        write.csv(features[, nonzeroFeatureNames], file = file.path(modelOutputDirectory, 
            paste("features_", cvPoint, ".csv", sep = "")), quote = F)
        
        melted <- reshape2::melt(data.frame(features[, nonzeroFeatureNames, drop = F], 
            labels = labels, check.names = F), id.vars = "labels")
        
        pdf(file.path(modelOutputDirectory, paste("features-", sub(pattern = "\\.", 
            replacement = "_", x = cvPoint), ".pdf", sep = "")), width = 4, height = length(nonzeroFeatureNames) * 
            1.5)
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

citrus.overlapDensityPlot <- function(clusterDataList, backgroundData) {
    combined <- data.frame(check.names = F, check.rows = F)
    for (clusterName in names(clusterDataList)) {
        combined <- rbind(combined, cbind(melt(clusterDataList[[clusterName]], varnames = c("row", 
            "marker")), clusterId = clusterName, src = "Cluster"))
    }
    p <- (ggplot2::ggplot(data = combined, ggplot2::aes(x = value, y = ..scaled.., fill = src)) 
            + ggplot2::geom_density() 
            + ggplot2::facet_grid(clusterId ~ marker, scales = "free") 
            + ggplot2::geom_density(data = cbind(reshape2::melt(backgroundData, varnames = c("row", "marker")), src = "Background")) 
            + ggplot2::theme_bw() 
            + ggplot2::scale_fill_manual(values = c(Background = rgb(0.3, 0.3, 1, 0.2), Cluster = rgb(1, 0.3, 0.3, 0.5))) 
            + ggplot2::theme(legend.position = "bottom", axis.text.y = ggplot2::element_blank(), 
                axis.ticks.y = ggplot2::element_blank(), axis.title = ggplot2::element_blank()) 
            + ggplot2::labs(fill = "Distribution:")
        )
    print(p)
}

citrus.plotModelClusters <- function(differentialFeatures, modelOutputDirectory, 
    clusterAssignments, citrus.combinedFCSSet, clusteringColumns, ...) {
    for (cvPoint in names(differentialFeatures)) {
        clusterIds <- as.numeric(differentialFeatures[[cvPoint]][["clusters"]])
        outputFile <- file.path(modelOutputDirectory, paste("clusters-", sub(pattern = "\\.", 
            replacement = "_", x = cvPoint), ".pdf", sep = ""))
        citrus.plotClusters(clusterIds, clusterAssignments = clusterAssignments, 
            citrus.combinedFCSSet, clusteringColumns, outputFile = outputFile, ...)
    }
}

#' Plot cluster histograms
#' 
#' Plot expression of markers in cluster cells relative to all cells
#' 
#' @param clusterIds Vector of cluster IDs to plot
#' @param clusterAssignments List containing indicies of cells assigned to each cluster.
#' @param citrus.combinedFCSSet Combined FCS data that was clustered.
#' @param clusteringColumns Columns for which to plot distributions
#' @param conditions Vector of conditions clustering was performed on.
#' @param outputFile If not \code{NULL}, plot is written to \code{outputFile}.
#' @param ... Other parameters (ignored).
#' 
#' @author Robert Bruggner
#' @export
citrus.plotClusters <- function(clusterIds, clusterAssignments, citrus.combinedFCSSet, 
    clusteringColumns, conditions = NULL, outputFile = NULL, ...) {
    
    data <- citrus.combinedFCSSet$data
    
    if (!is.null(outputFile)) {
        pdf(file = outputFile, width = (2.2 * length(clusteringColumns) + 2), height = (2 * 
            length(clusterIds)))
    }
    clusterDataList <- list()
    for (clusterId in sort(clusterIds)) {
        if (length(clusterAssignments[[clusterId]]) > 2500) {
            clusterDataList[[as.character(clusterId)]] <- data[clusterAssignments[[clusterId]], 
                clusteringColumns][sample(1:length(clusterAssignments[[clusterId]]), 
                2500), ]
        } else {
            clusterDataList[[as.character(clusterId)]] <- data[clusterAssignments[[clusterId]], 
                clusteringColumns]
        }
        
        colnames(clusterDataList[[as.character(clusterId)]]) <- .getDisplayNames(citrus.combinedFCSSet, 
            clusteringColumns)
        
    }
    if (nrow(data) > 2500) {
        bgData <- data[sample(1:nrow(data), 2500), clusteringColumns]
    } else {
        bgData <- data[, clusteringColumns]
    }
    colnames(bgData) <- colnames(clusterDataList[[1]])
    citrus.overlapDensityPlot(clusterDataList = clusterDataList, backgroundData = bgData)
    if (!is.null(outputFile)) {
        dev.off()
    }
}

#' Plot results of a Citrus regression analysis
#' 
#' Makes many plots showing results of a Citrus analysis
#' 
#' @name plot.citrus.regressionResult
#' @export 
#' 
#' @param citrus.regressionResult A \code{citrus.regressionResult} object.
#' @param outputDirectory Full path to output directory for plots.
#' @param citrus.foldClustering A \code{citrus.foldClustering} object.
#' @param citrus.foldFeatureSet A \code{citrus.foldFeatureSet} object. 
#' @param citrus.combinedFCSSet A \code{citrus.combinedFCSSet} object.
#' @param plotTypes Vector of plots types to make. Valid options are \code{errorRate} (Cross-validated error rates for predictive models),
#' \code{stratifyingFeatures} (plots of non-zero model features),\code{stratifyingClusters} (plots of clustering marker distributions in stratifying clusters)
#' @author Robert Bruggner
#' 
plot.citrus.regressionResult <- function(citrus.regressionResult, outputDirectory, 
    citrus.foldClustering, citrus.foldFeatureSet, citrus.combinedFCSSet, plotTypes = c("errorRate", 
        "stratifyingFeatures", "stratifyingClusters"), byCluster = FALSE, allFeatures = FALSE, ...) {
    addtlArgs <- list(...)
    
    theme <- "black"
    if ("theme" %in% names(addtlArgs)) {
        theme <- addtlArgs[["theme"]]
    }
    
    modelType <- citrus.regressionResult$modelType
    
    cat(paste("Plotting results for", modelType, "\n"))
    
    # Make model output directoy
    modelOutputDirectory <- file.path(outputDirectory, paste0(modelType, "_results"))
    dir.create(modelOutputDirectory, showWarnings = F, recursive = T)
    
    if ("errorRate" %in% plotTypes) {
        cat("Plotting Error Rate\n")
        citrus.plotTypeErrorRate(modelType = modelType, modelOutputDirectory = modelOutputDirectory, 
            regularizationThresholds = citrus.regressionResult$regularizationThresholds, 
            thresholdCVRates = citrus.regressionResult$thresholdCVRates, finalModel = citrus.regressionResult$finalModel$model, 
            cvMinima = citrus.regressionResult$cvMinima, family = citrus.regressionResult$family)
    }
    
    if ("stratifyingFeatures" %in% plotTypes) {
        cat("Plotting Stratifying Features\n")
        do.call(paste("citrus.plotModelDifferentialFeatures", citrus.regressionResult$family, 
            sep = "."), args = list(differentialFeatures = citrus.regressionResult$differentialFeatures, 
            features = citrus.foldFeatureSet$allFeatures, modelOutputDirectory = modelOutputDirectory, 
            labels = citrus.regressionResult$labels, byCluster = byCluster, allFeatures = allFeatures))
    }
    
    if ("stratifyingClusters" %in% plotTypes) {
        cat("Plotting Stratifying Clusters\n")
        citrus.plotModelClusters(differentialFeatures = citrus.regressionResult$differentialFeatures, 
            modelOutputDirectory = modelOutputDirectory, clusterAssignments = citrus.foldClustering$allClustering$clusterMembership, 
            citrus.combinedFCSSet = citrus.combinedFCSSet, clusteringColumns = citrus.foldClustering$allClustering$clusteringColumns, 
            ...)
    }
    
}

#' Plot a citrus.full.result
#' @name plot.citrus.full.result
#' @param citrus.full.result A citrus.full.result object
#' @param outputDirectory Full path to directory in which to place plot output. 
#' 
#' @export
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



