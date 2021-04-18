

#' @export
get_signifcant_features_matrix <- function(sam.model) {
    browser()
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

