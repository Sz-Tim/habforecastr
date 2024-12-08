% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/metrics.R
\name{calc_R2}
\alias{calc_R2}
\title{Calculate pseudo-R2s}
\usage{
calc_R2(dat.df, type = "mf", ...)
}
\arguments{
\item{dat.df}{A dataframe containing the data to be evaluated.}

\item{type}{A character string specifying the type of pseudo-R2 to calculate. Options are "mf" for McFadden and "vz" for Veall-Zimmermann. Default is "mf".}
}
\value{
A dataframe with the calculated pseudo-R2 values.
}
\description{
This function calculates pseudo-R2 values for logistic regression models. It currently supports McFadden and the Veall-Zimmermann correction of the Aldrich-Nelson pseudo-R2.
}
\examples{
\dontrun{
library(tidyverse)
dat.df <- tibble(alert = factor(c("A1", "A0", "A1", "A0")), prA1 = runif(4), model = factor(c("model1", "model1", "model2", "model2")), PCA = FALSE, covSet = "set1")
r2_values <- calc_R2(dat.df)
}
}
