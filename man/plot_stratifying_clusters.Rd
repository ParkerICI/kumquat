% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.R
\name{plot_stratifying_clusters}
\alias{plot_stratifying_clusters}
\title{Plotting stratifying clusters}
\usage{
plot_stratifying_clusters(
  citrus.model,
  clusters.data,
  output.dir,
  col.names = names(clusters.data),
  by.cluster = FALSE
)
}
\arguments{
\item{citrus.model}{An object of class \code{citrus.regressionResult}, as returned from
\code{\link{get_model}}}

\item{clusters.data}{A \code{data.frame} containing the clusters data. Each row is a cell
and column represents the expression of different markers (non numeric columns are dropped)}

\item{output.dir}{The directory that will contain the output. Additional directories and 
sub-directories will be created in \code{output.dir}}

\item{col.names}{Which columns of \code{clustersData} to include in the plot. The default is to include all
columns in the plot}

\item{by.cluster}{If this is \code{TRUE} this function will generate a separate plot 
for each cluster, otherwise or a single plot with all the clusters. 
The latter option may take quite long to plot if there are several clusters}
}
\description{
Plotting stratifying clusters
}
