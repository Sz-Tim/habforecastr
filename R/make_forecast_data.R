

#' Make forecast data
#'
#' This function loads compiled datasets and applies the appropriate recipe to
#' prep a new dataset for a specific covSet.
#'
#' @param y_i
#' @param covSet
#' @param data.dir
#' @param responses
#'
#' @return
#' @export
#'
#' @examples
make_forecast_data <- function(y_i, covSet, data.dir, responses=c(alert="alert")) {

  library(tidyverse); library(tidymodels)
  dir.create(glue("{data.dir}/compiled/"), recursive=T, showWarnings=F)

  # covariate set / response info
  id <- covSet$id
  y.i <- covSet$y
  y_i.i <- y_i |> filter(abbr==y.i)

  # covariate type columns
  col_metadata <- c("obsid", "y", "date", "year", "yday", "siteid", "lon", "lat")
  col_resp <- c("lnN", "tl", "alert")
  col_cmems <- readRDS("data/cmems_vars.rds")
  col_wrf <- readRDS("data/wrf_vars.rds")

  # All possible covariates
  all_covs <- make_all_covs(col_cmems, col_wrf, y_i)

  # Load full forecast dataset
  obs.ls <- load_dataset_y(data.dir, y.i, col_metadata, col_resp, all_covs)
  if(!all(col_resp %in% names(obs.ls))) {
    missing_cols <- which(! col_resp %in% names(obs.ls))
    obs.ls[col_resp[missing_cols]] <- NA
  }

  # Load recipes
  prep.ls <- readRDS(glue("data/0_init/compiled/{y.i}_{id}_dy_recipePrepped.rds"))
  prepPCA.ls <- readRDS(glue("data/0_init/compiled/{y.i}_{id}_dPCAy_recipePrepped.rds"))

  d.y <- list(test=map(prep.ls, ~bake(.x, obs.ls)))
  dPCA.y <- list(test=map(prepPCA.ls, ~bake(.x, obs.ls)))
  saveRDS(d.y, glue("{data.dir}/compiled/{y.i}_{id}_dy_forecast-{format(max(d.y$test$alert$date), '%F')}.rds"))
  saveRDS(dPCA.y, glue("{data.dir}/compiled/{y.i}_{id}_dPCAy_forecast-{format(max(dPCA.y$test$alert$date), '%F')}.rds"))

}
