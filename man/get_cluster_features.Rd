% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kumquat.R
\name{get_cluster_features}
\alias{get_cluster_features}
\title{Calculate cluster features for model building}
\usage{
get_cluster_features(
  tab,
  predictors = NULL,
  metadata.tab = NULL,
  variable.var = "variable",
  value.var = "value",
  endpoint.grouping = NULL,
  sample.col = "sample.id"
)
}
\arguments{
\item{tab}{A \code{data.frame} of data in "molten" format (see Details)}

\item{predictors}{Columns in \code{tab} (after optionally merging with \code{metadata.tab}) that identify predictors. 
The data will be processed using \code{reshape2::dcast}
according to the formula \code{endpoint.grouping1 + endpoint.grouping2 + ... ~ variable.var + predictors1 + predictors2 + ...}}

\item{metadata.tab}{Optional. A \code{data.frame} containing sample-level metadata to be merged with \code{tab} (see Details)}

\item{variable.var}{The column in \code{tab} that identifies the variable}

\item{value.var}{The column in \code{tab} that identfies the value}

\item{endpoint.grouping}{Columns in \code{tab} (after optionally merging with \code{metadata.tab}) that identify the grouping of the endpoint
(see Details). The data will be processed using \code{reshape2::dcast}
according to the formula \code{endpoint.grouping1 + endpoint.grouping2 + ... ~ variable.var + predictors1 + predictors2 + ...}
The combination of \code{predictors} and \code{endpoint.grouping} must uniquely identify every value in \code{tab}
The function will throw an error if this is not the case.}

\item{sample.col}{Optional, only used if \code{metadata.tab} is provided. The name of the column that will be used to 
merge \code{tab} with \code{metadata.tab}}
}
\value{
Returns a matrix where each row corresponds to a combination of the levels of the variables specified in \code{endpoint.grouping}, and the columns are
      numeric features corresponding to combinations of the levels of the \code{predictors}
}
\description{
This function takes input data and (optionally) a table of sample metadata, and rearranges data into a matrix of features to be used for model building
}
\details{
The input table needs to be in molten format (i.e. see \code{reshape2::melt}) with \code{variable.var} and
\code{value.var} columns identifying variables and their values (for instance cell population abundances). The
\code{metadata.tab}, if provided, must contain a column (identified by the \code{sample.col} function argument), which matches the names of the samples in
\code{tab} (i.e. the part after the \code{@}, "sample1" in the above example). The rest of the columns in \code{metadata.tab} represent file-level
metadata, which is used to identify the data corresponding to a given combination of predictors (see below)
An example will help clarify the working of this function. Suppose you have collected data from multiple patients at multiple timepoints and under multiple
stimulation conditions.
In this case the \code{metadata.tab} would look like this
\itemize{
  \item{\code{sample.id}}{ This is used to merge sample metadata with the input data (see above)}
  \item{\code{timepoint}}{ The timepoint information}
  \item{\code{condition}}{ The stimulation condition}
  \item{\code{subject}}{ The subjet each file was derived from}
}
Let's assume a few different scenarios.
\enumerate{
  \item You have subject level information (e.g. "responder" vs "non-responder") and you want to predict whether any combination of the \code{timepoint} and
        \code{condition} information predicts this outcome. In this case you would call the function with \code{predictors = c("condition", "timepoint")} and
        \code{endpoint.grouping = "sample"}. The features in the resulting output would look like \code{cluster_1_feature1_condition_timepoint}
  \item You have subject and timepoint level information, and you want to see if any of the stimulation conditions predicts it. In this case you would call
        the function with \code{predictors = c("condition")} and \code{endpoint.grouping = c("sample", "timepoint")}. The features in the resulting output
        would look like \code{cluster_1_feature1_condition}
}
Internally this function uses \code{reshape2::dcast} to structure the data in the appropriate format with the following formula (see the \code{reshape2::dcast}
documentation for details on how the formula is interpreted): 
\code{endpoint.grouping1 + endpoint.grouping2 + ... ~ variable.var + predictors1 + predictors2 + ...}
}
