% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/citrus.classModel.R,
%   R/citrus.continuousModel.R, R/citrus.model.R
\name{citrus.generateRegularizationThresholds}
\alias{citrus.generateRegularizationThresholds}
\alias{citrus.generateRegularizationThresholds.classification}
\alias{citrus.generateRegularizationThresholds.continuous}
\title{Generate model regularization thresholds}
\usage{
citrus.generateRegularizationThresholds.classification(
  features,
  labels,
  modelType,
  n = 100,
  ...
)

citrus.generateRegularizationThresholds.continuous(
  features,
  labels,
  modelType,
  n = 100,
  ...
)

citrus.generateRegularizationThresholds(
  features,
  labels,
  modelType,
  family,
  n = 100,
  ...
)
}
\arguments{
\item{features}{Features used to construct model}

\item{labels}{Endpoint lables for samples and features}

\item{modelType}{Method used to construct endpoint model. Valid options are: \code{pamr} and \code{glmnet}.}

\item{n}{Number of regularization thresholds to generate}

\item{...}{Other arguments passed to model-fitting methods}

\item{family}{Model family. Valid options are: \code{classification} and \code{continuous}.}
}
\value{
A vector of regularization threshold values.
}
\description{
Generate a range of regularization thresholds for model construction
}
\author{
Robert Bruggner
}
