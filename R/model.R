#' Calculate a prediction model
#' 
#' This function is the main wrapper around the Citrus model building functionality
#' 
#' @param endpoint A vector of length equal to \code{nrow(features)} containing 
#'   the values to be used as prediction endpoints. If this vertor is numeric
#'   a model of family \code{"continuous"} is generated, otherwise if it is
#'   a character or factor vector the model will be of family 
#'   \code{"classification"}
#' @param model.type The type of model, either \code{"pamr"}, \code{"sam"} or
#'   \code{"glmnet"}
#' @param ... additional arguments
#' @inherit citrus.endpointRegress return
#' 
#' @export 
get_model <- function(features, endpoint, model.type, ...) {
    if(is.matrix(features))
        citrus.features <- list(allFeatures = features, nFolds = 1)
    else
        citrus.features <- features
    
    family <- "classification"
    
    if (is.character(endpoint) || is.factor(endpoint)) {
        family <- "classification"
        endpoint <- as.factor(endpoint)

    } else family <- "continuous"
    
    citrus.res <- citrus.endpointRegress(model.type, citrus.foldFeatureSet = citrus.features, 
                                                 labels = endpoint, family = family, ...)
    
    citrus.res$allFeatures <- citrus.features$allFeatures
    return(citrus.res)
    
}

#' High level wrapper
#' 
#' This function is a high-level wrapper around the whole model building pipeline. It will calculate a model
#' and write a series of standardized plots. If you want more control over the plotting output, please
#' run \code{\link{get_model}} to calculate the model, and the individual plotting functions
#' (\code{\link{plot_error_rate}}, \code{\link{plot_stratifying_features}}, \code{\link{plot_stratifying_clusters}})
#' 
#' @export
run_analysis <- function(features, endpoint, output.directory, model.type, clusters.data = NULL, ...) {
    citrus.features <- list(allFeatures = features, nFolds = 1)
    
    out.dir <- file.path(output.directory, sprintf("%s_results", model.type))
    dir.create(out.dir, recursive = TRUE, showWarnings = FALSE)
    message("Calculating model")
    citrus.res <- get_model(citrus.features, endpoint, model.type, ...)
    
    message("Plotting model error rate")
    plot_error_rate(citrus.res, file.path(out.dir, "model_error_rate.pdf"))
    
    message("Plotting stratifying features")
    plot_stratifying_features(citrus.res, out.dir, by.cluster = FALSE, all.features = FALSE)
    
    if(!is.null(clusters.data)) {
        message("Plotting clusters")
        plot_stratifying_clusters(citrus.res, clusters.data = clusters.data, output.dir = out.dir, 
                                  by.cluster = TRUE)
    }
    saveRDS(citrus.res, file = file.path(out.dir, "model.rds"))
    
    
    return(invisible(citrus.res))
    
}


