#' @export
convert_to_citrus_featureset <- function(tab) {
    endpoint.grouping <- grep("cluster_", names(tab), invert = TRUE, value = TRUE)
    
    rnames <- do.call(paste, list(tab[, endpoint.grouping], sep = "_"))
    cnames <- setdiff(colnames(tab), endpoint.grouping)
    ret <- as.matrix(tab[, cnames])
    row.names(ret) <- rnames
    colnames(ret) <- cnames
    
    return(list(allFeatures = ret, nFolds = 1))
    
    return(ret)
}


#' @export 
get_model <- function(citrus.features, endpoint, model.type) {
    family <- NULL
    
    if (is.character(endpoint) || is.factor(endpoint)) {
        family <- "classification"
        endpoint <- as.factor(endpoint)
        
    } else family <- "continuous"
    
    citrus.res <- citrus.endpointRegress(model.type, citrus.foldFeatureSet = citrus.features, 
                                                 labels = endpoint, family = family)
    
    citrus.res$allFeatures <- citrus.features$allFeatures
    return(citrus.res)
    
}

#' High level wrapper
#' @export
run_analysis <- function(citrus.features, endpoint, output.directory, model.type, plot.by.cluster = FALSE,
                                plot.all.features = FALSE, clusters.data = NULL) {

    out.dir <- file.path(output.directory, sprintf("%s_results", model.type))
    dir.create(out.dir, recursive = TRUE, showWarnings = FALSE)
    message("Calculating model")
    citrus.res <- get_model(citrus.features, endpoint, model.type)
    
    message("Plotting model error rate")
    plot_error_rate(citrus.res, out.dir)
    
    message("Plotting stratifying features")
    plot_stratifying_features(citrus.res, out.dir)
    
    return(invisible(citrus.res))
    
}

#' @export
get_signifcant_features_matrix <- function(sam.model) {
    tabs.list <- sam.model$siggenes.table[c("genes.up", "genes.lo")]
    
    ret <- lapply(tabs.list, function(x) {
        if(is.vector(x))
            df <- data.frame(t(x), check.names = FALSE, stringsAsFactors = FALSE)
        else
            df <- data.frame(x, check.names = FALSE, stringsAsFactors = FALSE)
        
        df$"Gene Name" <- NULL
        
        temp <- sapply(setdiff(names(df), "Gene ID"), function(x) {as.numeric(df[, x])})
        df <- data.frame(feature = df$"Gene ID", temp, check.names = FALSE, stringsAsFactors = FALSE)
        return(df)
    })
    
    ret <- do.call(rbind, ret)
    row.names(ret) <- NULL
    
    names(ret) <- gsub("q-value\\(%\\)", "q-value", names(ret))
    ret[, "q-value"] <- ret[, "q-value"] / 100
    
    return(ret)    
}

#' @export
select_score_by_cluster <- function(tab, col.name, f, use.abs = TRUE) {
    cl <- sapply(strsplit(tab$feature, "_"), "[", 2)
    tab$cluster <- cl
    
    tab <- plyr::ddply(tab, ~cluster, function(x) {
        v <- x[, col.name]
        if(use.abs)
            v <- abs(v)
        w <- which(v == f(v))
        return(x[w, ])
        
    })
    return(tab)
}
