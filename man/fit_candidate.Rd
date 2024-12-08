% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fit_candidate.R
\name{fit_candidate}
\alias{fit_candidate}
\title{Wrapper to fit a model}
\usage{
fit_candidate(mod, resp, form.ls, d.ls, opts, tunes, out.dir, y, suffix = NULL)
}
\arguments{
\item{mod}{A character string specifying the model to fit (e.g., `"Ridge"`, `"ENet"`, `"RF"`, `"NN"`, `"MARS"`, `"Boost"`, `"lgbm"`, `"HB"`).}

\item{resp}{A character string specifying the response variable.}

\item{form.ls}{A list of formulas for the models.}

\item{d.ls}{A list of data frames for the models.}

\item{opts}{A list of options for model fitting, including resampling and control parameters.}

\item{tunes}{A list of tuning parameters for the models.}

\item{out.dir}{A character string specifying the output directory for saving model results.}

\item{y}{A character string specifying the target variable.}

\item{suffix}{A character string specifying a suffix for the model ID. Default is `NULL`.}
}
\value{
None. Fitted object is stored.
}
\description{
This function fits a specified model to the data, including machine learning and hierarchical Bayesian models.
}
\examples{
\dontrun{
library(tidymodels)
form.ls <- list(alert = list(ML = y ~ ., ML_PCA = y ~ .))
d.ls <- list(alert = data.frame(y = factor(c("A1", "A2")), x1 = rnorm(2), x2 = rnorm(2)))
opts <- vfold_cv(d.ls$alert, v = 5)
tunes <- list(Ridge = 10, ENet = 10, RF = 10, NN = 10, MARS = 10, Boost = 10, lgbm = 10)
out.dir <- "model_output"
fit_candidate("Ridge", "alert", form.ls, d.ls, opts, tunes, out.dir, "y")
}
}
