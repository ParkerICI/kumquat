#' kumquat: Identification of stratifying subpopulations in Flow Cytometry Data (based on the Citrus package)
#' 
#' Automated identification of clusters whose behavior differs between sample groups.
#' 
#' @docType package
#' @name kumquat
NULL


#' Normalize numeric data in multiple steps
#'
#' This function normalizes numeric data in a table in multiple steps, based on groups defined by different combinations of categorical variables. The table
#' is assumed to be in molten format (e.g. see \code{reshape2::melt}), with variable and value columns identifying observations. 
#' This function will then proceed to normalize the data based on the value identified by a categorical variable, 
#' then normalize the normalized data again by another value etc. (see below for Details)
#' At each step of the normalization the table is grouped using the \code{variable.var}, \code{subject.col} and all the columns
#' in \code{names(norm.template)}. After this grouping, for every group, there can be only one row for the value of the current grouping
#' variable that has been selected as a basis for normalization. In other words the function will not allow you to normalize a vector of values
#' by another vector of values, it will only allow normalization of a vector by an individual number. This is done to prevent the result to depend
#' on the ordering of the table.
#' 
#' An example should help clarify the working of this function. Assume you have a dataset where different variables have been measured for multiple subjects,
#' under different stimulation conditions, and at different timepoints. For each variable you want the data at each timepoint to be normalized by the value in
#' the "unstim" condition. Then you want this data to be further normalized by the value at the "baseline" timepoint. Assume \code{tab} is in molten
#' format and has the following columns
#'   \itemize{
#'     \item{\code{variable}}: identifies the variable
#'     \item{\code{value}}: the corresponding value of the variable
#'     \item{\code{timepoint}}: categorical variable that identifies the timepoint
#'     \item{\code{condition}}: categorical variable that identified the condition
#'     \item{\code{subject}}: categorical variable that identifies data belonging to the same subject (all the normalization is done within subject)
#'   }
#' To achieve the result described above, you would invoke this function as \code{multistep_normalize(tab, list(condition = "unstim", timepoint = "baseline"), "subject")}.
#' Note that the function would fail if you only specify a single variable (either \code{condition} or \code{timepoint}), because a single variable is not enough
#' to identify a single value for normalization, since you have multiple conditions for each timepoint and viceversa.
#'
#' @param tab The input \code{data.frame} See the details for assumption about its structure
#' @param norm.template A named list identyfying which categorical variables should be used to group data for normalization. The values in the list
#'   represent the value of the corresponding variable that identify the rows that are used as reference for normalization at each step. The data will be normalized
#'   in the same order specified by this list (i.e. data will be normalized according to the first variable, then again according to the second etc.)
#' @param subject.var The name of the column that identifies different subjects in \code{tab}. All normalization operations are done within the subgroups
#'   identified by this variable (i.e. data will never be normalized across subsets identified by different values of subject.var)
#'
#'
#'
#' @export
multistep_normalize <- function(tab, norm.template, subject.col, variable.var = "variable", value.var = "value") {
    var.names <- names(norm.template)
    var.values <- unlist(norm.template, use.names = F)
    
    ret <- tab
    
    for(i in 1:length(var.names)) {
        variable.names <- c(variable.var, subject.col, var.names[var.names != var.names[i]])
        mutate.s <- sprintf("%s / %s[%s == '%s']", value.var, value.var, var.names[i], var.values[i])
        filter.s <- sprintf("%s != '%s'", var.names[i], var.values[i])
        
        # Check for uniqueness of normalization values
        dplyr::group_by_(ret, .dots = variable.names) %>%
            dplyr::do({
                x <- .[[var.names[i]]]
                if(length(x[x == var.values[i]]) != 1)
                    stop("This combination of variables does not identify a single reference value for normalization")
                .
            })
        
        
        ret <- dplyr::group_by_(ret, .dots = variable.names) %>%
                   dplyr::mutate_(.dots = setNames(mutate.s, value.var)) %>%
                   dplyr::filter_(.dots = filter.s)
        
    }
    
    return(ret)
}


#' Calculate cluster features for model building
#'
#' This function takes input data and (optionally) a table of sample metadata, and rearranges data into a matrix of features to be used for model building
#'
#' The input table needs to be in molten format (i.e. see \code{reshape2::melt}) with \code{variable.var} and
#' \code{value.var} columns identifying variables and their values (for instance cell population abundances). The
#' \code{metadata.tab}, if provided, must contain a column (identified by the \code{sample.col} function argument), which matches the names of the samples in
#' \code{tab} (i.e. the part after the \code{@}, "sample1" in the above example). The rest of the columns in \code{metadata.tab} represent file-level
#' metadata, which is used to identify the data corresponding to a given combination of predictors (see below)
#' An example will help clarify the working of this function. Suppose you have collected data from multiple patients at multiple timepoints and under multiple
#' stimulation conditions.
#' In this case the \code{metadata.tab} would look like this
#' \itemize{
#'   \item{\code{sample.id}}{ This is used to merge sample metadata with the input data (see above)}
#'   \item{\code{timepoint}}{ The timepoint information}
#'   \item{\code{condition}}{ The stimulation condition}
#'   \item{\code{subject}}{ The subjet each file was derived from}
#' }
#' Let's assume a few different scenarios.
#' \enumerate{
#'   \item You have subject level information (e.g. "responder" vs "non-responder") and you want to predict whether any combination of the \code{timepoint} and
#'         \code{condition} information predicts this outcome. In this case you would call the function with \code{predictors = c("condition", "timepoint")} and
#'         \code{endpoint.grouping = "sample"}. The features in the resulting output would look like \code{cluster_1_feature1_condition_timepoint}
#'   \item You have subject and timepoint level information, and you want to see if any of the stimulation conditions predicts it. In this case you would call
#'         the function with \code{predictors = c("condition")} and \code{endpoint.grouping = c("sample", "timepoint")}. The features in the resulting output
#'         would look like \code{cluster_1_feature1_condition}
#' }
#' Internally this function uses \code{reshape2::dcast} to structure the data in the appropriate format with the following formula (see the \code{reshape2::dcast}
#' documentation for details on how the formula is interpreted): 
#' \code{endpoint.grouping1 + endpoint.grouping2 + ... ~ variable.var + predictors1 + predictors2 + ...}
#'
#' @param tab A \code{data.frame} of data in "molten" format (see Details)
#' @param variable.var The column in \code{tab} that identifies the variable
#' @param value.var The column in \code{tab} that identfies the value
#' @param metadata.tab Optional. A \code{data.frame} containing sample-level metadata to be merged with \code{tab} (see Details)
#' @param predictors Columns in \code{tab} (after optionally merging with \code{metadata.tab}) that identify predictors. 
#'   The data will be processed using \code{reshape2::dcast}
#'   according to the formula \code{endpoint.grouping1 + endpoint.grouping2 + ... ~ variable.var + predictors1 + predictors2 + ...}
#' @param endpoint.grouping  Columns in \code{tab} (after optionally merging with \code{metadata.tab}) that identify the grouping of the endpoint
#'   (see Details). The data will be processed using \code{reshape2::dcast}
#'   according to the formula \code{endpoint.grouping1 + endpoint.grouping2 + ... ~ variable.var + predictors1 + predictors2 + ...}
#'   The combination of \code{predictors} and \code{endpoint.grouping} must uniquely identify every value in \code{tab}
#'   The function will throw an error if this is not the case.
#' @param sample.col Optional, only used if \code{metadata.tab} is provided. The name of the column that will be used to 
#'   merge \code{tab} with \code{metadata.tab} 
#' @return Returns a matrix where each row corresponds to a combination of the levels of the variables specified in \code{endpoint.grouping}, and the columns are
#'       numeric features corresponding to combinations of the levels of the \code{predictors} 
#' @export

get_cluster_features <- function(tab, predictors = NULL, metadata.tab = NULL, variable.var = "variable", value.var = "value", 
                                 endpoint.grouping = NULL, sample.col = "sample.id") {
    df <- tab
    if(!is.null(metadata.tab))
        df <- merge(df, metadata.tab, by = sample.col)
    
    ret <- NULL
    
    formula.exp <- as.formula(sprintf("%s ~ %s", paste(endpoint.grouping, collapse = "+"),
                                      paste(c(variable.var, predictors), collapse = "+")))
    
    if(nrow(df) != nrow(unique(df[, c(variable.var, predictors, endpoint.grouping)])))
        stop("The combination of predictors and endpoint.grouping does not uniquely identify every row in metadata.tab")
    ret <- reshape2::dcast(df, formula.exp)
    
    m <- ret[, !(names(ret) %in% endpoint.grouping)]
    row.names(m) <- paste(ret[, endpoint.grouping], sep = "_")
    return(as.matrix(m))
}




