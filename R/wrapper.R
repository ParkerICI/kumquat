#' @export
convert_to_citrus_featureset <- function(tab, endpoint.grouping) {
    rnames <- do.call(paste, list(tab[, endpoint.grouping], sep = "_"))
    cnames <- setdiff(colnames(tab), endpoint.grouping)
    ret <- as.matrix(tab[, cnames])
    row.names(ret) <- rnames
    colnames(ret) <- cnames
    
    return(list(allFeatures = ret, nFolds = 1))
    
    return(ret)
}

#' @export
run_citrus_analysis <- function(citrus.features, endpoint, working.directory, model.type) {
    family <- NULL
    
    if (is.character(endpoint) || is.factor(endpoint)) {
        family <- "classification"
        endpoint <- as.factor(endpoint)
        
    } else family <- "continuous"
    
    citrus.res <- citrus::citrus.endpointRegress(model.type, citrus.foldFeatureSet = citrus.features, labels = endpoint, 
        family = family)
    
    plot(citrus.res, working.directory, citrus.foldClustering = NULL, citrus.foldFeatureSet = citrus.features, citrus.combinedFCSSet = NULL, 
        "stratifyingFeatures")
    return(citrus.res)
    
}
