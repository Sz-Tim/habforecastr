% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predictions.R
\name{summarise_predictions}
\alias{summarise_predictions}
\title{This function summarises and aggregates predictions for a set of models, including both PCA and non-PCA models.}
\usage{
summarise_predictions(d.y, dPCA.y, resp, fit.dir, y_i.i, suffix = NULL)
}
\arguments{
\item{d.y}{A list of data frames containing the data for each response variable.}

\item{dPCA.y}{A list of data frames containing the PCA-transformed data for each response variable.}

\item{resp}{A character string specifying the response variable.}

\item{fit.dir}{A character string specifying the directory where the model fits are stored.}

\item{y_i.i}{A data frame with information about the variables of interest.}

\item{suffix}{A character string specifying a suffix for the model ID. Default is NULL.}
}
\value{
A data frame with the summarised and aggregated predictions.
}
\description{
This function summarises and aggregates predictions for a set of models, including both PCA and non-PCA models.
}
\examples{
\dontrun{
library(tidyverse)
library(glue)
library(tidymodels)
d.y <- list(response = data.frame(obsid = 1:10, y = rnorm(10), date = Sys.Date() - 1:10, siteid = 1:10, response = rnorm(10)))
dPCA.y <- list(response = data.frame(obsid = 1:10, y = rnorm(10), date = Sys.Date() - 1:10, siteid = 1:10, response = rnorm(10)))
resp <- "response"
fit.dir <- "path/to/fits"
y_i.i <- data.frame(abbr = "response")
predictions <- summarise_predictions(d.y, dPCA.y, resp, fit.dir, y_i.i)
}
}
