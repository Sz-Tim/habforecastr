
#' Create recipe and prepare using training data
#'
#' This function creates a recipe for data preprocessing and prepares it using the training data.
#'
#' @param train.df A data frame containing the training data.
#' @param response A character string specifying the response variable.
#' @param covsExclude A character vector of covariates to exclude from the recipe. Default is NULL.
#' @param dimReduce A logical value indicating whether to apply dimensionality reduction (PCA). Default is FALSE.
#'
#' @return A prepared recipe object.
#' @export
#'
#' @examples
#' \dontrun{
#' library(tidymodels)
#' train.df <- data.frame(obsid = 1:10, y = rnorm(10), date = Sys.Date() - 1:10, siteid = 1:10, year = 2023, lon = runif(10), lat = runif(10), yday = 1:10, response = rnorm(10))
#' response <- "response"
#' recipe <- prep_recipe(train.df, response)
#' }
prep_recipe <- function(train.df, response, covsExclude=NULL, dimReduce=FALSE) {
  respExclude <- grep(response, c("lnN", "tl", "alert"), value=T, invert=T)
  pred_vars <- names(train.df |> select(-all_of(response)))
  rec <- recipe(alert ~ ., train.df) |>
    update_role(all_of(pred_vars), new_role="predictor") |>
    update_role(obsid, y, date, siteid, year, new_role="ID") |>
    update_role(lon, lat, yday, new_role="RE") |>
    step_select(-any_of(respExclude)) |>
    step_dummy(all_factor_predictors()) |>
    step_logit(starts_with("prAlert"), offset=0.01) |>
    step_logit(ends_with("A1"), offset=0.01) |>
    step_interact(terms=~lon:lat, sep="X") |>
    step_interact(terms=~lnNWt1:all_predictors(), sep="X") |>
    # Wind x fetch x non-wind/non-interactions
    # step_interact(terms=~UWk:fetch:matches("^(?![UV])(?=.*Dir)(?!.*X).*", perl=T), sep="X") |>
    # step_interact(terms=~VWk:fetch:matches("^(?![UV])(?=.*Dir)(?!.*X).*", perl=T), sep="X") |>

    step_interact(terms=~UWk:fetch:matches("^[Precip|Shortwave|sst].*Dir(?!.*X).*", perl=T), sep="X") |>
    step_interact(terms=~VWk:fetch:matches("^[Precip|Shortwave|sst].*Dir(?!.*X).*", perl=T), sep="X") |>
    step_YeoJohnson(all_predictors()) |>
    step_normalize(all_predictors()) |>
    step_harmonic(yday, frequency=1, cycle_size=365, keep_original_cols=T) |>
    step_rename(ydayCos=yday_cos_1, ydaySin=yday_sin_1) |>
    step_interact(terms=~ydaySin:ydayCos, sep="X") |>
    step_bs(lon, deg_free=8, keep_original_cols=T) |>
    step_bs(lat, deg_free=8, keep_original_cols=T) |>
    step_bs(lonXlat, deg_free=8) |>
    step_rename_at(contains("_"), fn=~str_remove_all(.x, "_")) |>
    step_select(-all_of(covsExclude))
  if(dimReduce) {
    rec <- rec |>
      step_pca(all_predictors(), threshold=0.8)
  }
  rec |>
    prep(training=train.df)
}








#' Filter list of covariates following recipe thinning
#'
#' This function filters the list of covariates following recipe thinning.
#'
#' @param all_covs A list of all possible covariates.
#' @param data.y A list of data frames containing the data for each response variable.
#' @param covsExclude A character string specifying the covariates to exclude. Default is "NA".
#'
#' @return A list of filtered covariates.
#' @export
#'
#' @examples
#' \dontrun{
#' all_covs <- list(main = c("cov1", "cov2"), interact = c("cov3", "cov4"), nonHB = c("cov5", "cov6"))
#' data.y <- list(response1 = data.frame(cov1 = rnorm(10), cov2 = rnorm(10)), response2 = data.frame(cov3 = rnorm(10), cov4 = rnorm(10)))
#' filtered_covs <- filter_corr_covs(all_covs, data.y)
#' }
filter_corr_covs <- function(all_covs, data.y, covsExclude="NA") {
  uncorr_covs <- unique(unlist(map(data.y, ~unlist(map(.x, names)))))
  if(any(grepl("^PC", uncorr_covs))) {
    covs_incl <- list(main=grep("^PC", uncorr_covs, value=T),
                      interact=NULL,
                      nonHB=NULL)
  } else {
    uncorr_covs <- unique(c(uncorr_covs, "yday", "ydayCos", "ydaySin", "lon", "lat"))
    covs_incl <- map(all_covs, ~.x[.x %in% uncorr_covs])
  }

  covs_incl <- map(covs_incl, ~grep(covsExclude, .x, value=T, invert=T))
  return(covs_incl)
}







#' Get regex for excluded covariates given a covariate set
#'
#' This function generates a regular expression for excluded covariates based on a specified covariate set.
#'
#' @param covSet A list containing the covariate set information.
#'
#' @return A character string containing the regular expression for excluded covariates.
#' @export
#'
#' @examples
#' \dontrun{
#' covSet <- list(Avg = 0, Xf = 1, XN = 0, Del = 1)
#' regex <- get_excluded_cov_regex(covSet)
#' }
get_excluded_cov_regex <- function(covSet) {
  with(covSet, glue("Dt|",
                    "{ifelse(Avg==0, 'Avg|', 'NA|')}",
                    "{ifelse(Xf==0, 'Xfetch|', 'NA|')}",
                    "{ifelse(XN==0, 'lnNWt1X|', 'NA|')}",
                    "{ifelse(Del==0, 'Delta', 'NA')}"))
}
