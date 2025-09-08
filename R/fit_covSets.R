


#' Fit a subset of available covariates
#'
#' This function subsets the available covariates, fitting a specified model with that subset.
#'
#' @param y_i A data frame with information about the variables of interest, including abbreviations.
#' @param run_type A character string specifying the run type. Default is "0_init".
#' @param covSet A list containing the covariate set information.
#' @param mod A character string specifying the model to fit (e.g., "Ridge", "ENet", "RF", "NN", "MARS", "Boost", "lgbm", "HB").
#' @param train_prop A numeric value specifying the proportion of data to use for training. Default is 0.75.
#' @param nTuneVal A numeric value specifying the number of tuning values. Default is 2.
#' @param prior_strength A numeric value specifying the strength of the priors. Default is 1.
#' @param ncores A numeric value specifying the number of cores to use for parallel processing. Default is 4.
#' @param responses A character vector specifying the response variables. Default is c(alert = "alert").
#' @param rebalance_thresh Imbalance threshold determining whether to rebalance the training data using Synthetic Minority Oversampling TEchnique (SMOTE). If the proportion of alerts is less than the threshold, the training data will be rebalanced to the threshold. Default is 0, resulting in no rebalancing.
#' @param HB_CV Logical: Perform cross-validation with Bayesian models?
#'
#' @return None. Data, fitted objects, and logs are stored.
#' @export
#'
#' @examples
#' \dontrun{
#' library(tidymodels)
#' y_i <- data.frame(abbr = c("HAB", "Toxin"))
#' covSet <- list(id = "example", y = "HAB", prop_covs = 0.5, seed = 123)
#' mod <- "Ridge"
#' fit_covSet(y_i, run_type = "0_init", covSet, mod, train_prop = 0.75, nTuneVal = 2, prior_strength = 1, ncores = 4, responses = c(alert = "alert"))
#' }
fit_covSet <- function(y_i, run_type="0_init", covSet, mod, train_prop=0.75,
                       nTuneVal=2, prior_strength=1, ncores=4,
                       responses=c(alert="alert"), rebalance_thresh=0, HB_CV=FALSE) {

  # covariate set / response info
  id <- covSet$id
  y.i <- covSet$y
  y_i.i <- y_i |> filter(abbr==y.i)

  # directories
  data.dir <- glue("data/{run_type}/")
  fit.dir <- glue("out/{run_type}/model_fits/{id}/")
  cv.dir <- glue("{fit.dir}/cv/")
  ens.dir <- glue("out/{run_type}/ensembles/")
  out.dir <- glue("out/{run_type}/compiled/{id}/")
  log.dir <- glue("out/logs/{run_type}/")
  dir.create(glue("{data.dir}/compiled/"), recursive=T, showWarnings=F)
  dir.create(ens.dir, recursive=T, showWarnings=F)
  dir.create(cv.dir, recursive=T, showWarnings=F)
  dir.create(out.dir, recursive=T, showWarnings=F)
  dir.create(glue("{fit.dir}/vi/"), recursive=T, showWarnings=F)
  dir.create(log.dir, recursive=T, showWarnings=F)

  # covariate type columns
  col_metadata <- c("obsid", "y", "date", "year", "yday", "siteid", "lon", "lat")
  col_resp <- c("lnN", "tl", "alert")
  col_cmems <- readRDS("data/cmems_vars.rds")
  col_wrf <- readRDS("data/wrf_vars.rds")

  # All possible covariates
  all_covs <- make_all_covs(col_cmems, col_wrf, y_i)

  # Randomly select covariate subset
  set.seed(covSet$seed)
  if(grepl("fish", y_i.i$type)) {
    all_covs <- map(all_covs, ~grep("PrevYr", .x, value=TRUE, invert=TRUE))
  }
  if(grepl("hab", y_i.i$type)) {
    covs_exclude <- unlist(all_covs[2:3], use.names=F) %>%
      sample(x=., size=floor(length(.)*(1-covSet$prop_covs)), replace=F)
  } else {
    covs_exclude <- unlist(all_covs[2:4], use.names=F) %>%
      sample(x=., size=floor(length(.)*(1-covSet$prop_covs)), replace=F)
  }

  # testing/training splits
  obs.ls <- load_dataset_y(data.dir, y.i, col_metadata, col_resp, all_covs)

  set.seed(1003)
  if(train_prop < 1) {
    if(grepl("fish", y_i.i$type)) {
      obs.ls <- obs.ls |> arrange(date, siteid)
      obs.split <- initial_time_split(obs.ls, prop=train_prop)
    } else {
      obs.split <- group_initial_split(obs.ls, group=year, prop=train_prop)
    }
    saveRDS(obs.split, glue("{data.dir}/compiled/{y.i}_{id}_dataSplit.rds"))
    obs.train <- training(obs.split)
    obs.test <- testing(obs.split)
  } else {
    obs.train <- obs.ls
  }

  if(n_distinct(obs.train$alert)==1) {
    cat("No alerts in training data for", y.i, "!\n",
        "  Data range:", format(min(obs.ls$date), "%F"), "to", format(max(obs.ls$date), "%F"), "\n",
        "  Removing dataSplit and moving to next row of covSet.df\n")
    file.remove(glue("{data.dir}/compiled/{y.i}_{id}_dataSplit.rds"))
    return()
  }

  # rebalance training data using SMOTE to rebalance_thresh
  prop_alert <- mean(obs.train$alert == "A1")
  set.seed(1003)
  if(prop_alert < rebalance_thresh) {
    obs.train <- rebalance_smote(obs.train, dup_size=rebalance_thresh/prop_alert)
  }

  # prepare dataset
  prep.ls <- map(responses, ~prep_recipe(obs.train, .x, covs_exclude))
  prepPCA.ls <- map(responses, ~prep_recipe(obs.train, .x, covs_exclude, TRUE))
  d.y <- list(train=map(prep.ls, ~bake(.x, obs.train)))
  dPCA.y <- list(train=map(prepPCA.ls, ~bake(.x, obs.train)))
  if(train_prop < 1) {
    d.y$test <- map(prep.ls, ~bake(.x, obs.test))
    dPCA.y$test <- map(prepPCA.ls, ~bake(.x, obs.test))
  }
  saveRDS(prep.ls, glue("{data.dir}/compiled/{y.i}_{id}_dy_recipePrepped.rds"))
  saveRDS(prepPCA.ls, glue("{data.dir}/compiled/{y.i}_{id}_dPCAy_recipePrepped.rds"))
  saveRDS(d.y, glue("{data.dir}/compiled/{y.i}_{id}_dy_testPct-{train_prop}.rds"))
  saveRDS(dPCA.y, glue("{data.dir}/compiled/{y.i}_{id}_dPCAy_testPct-{train_prop}.rds"))
  covs <- filter_corr_covs(all_covs, d.y) |> map(~.x[! .x %in% covs_exclude])
  covsPCA <- names(dPCA.y$train[[1]] |> select(starts_with("PC")))

  # formulas
  form.ls <- map(
    responses,
    ~list(ML=formula(glue("{.x} ~ .")),
          ML_PCA=formula(glue("{.x} ~ .")),
          HB=make_HB_formula(.x, c(covs$main, covs$interact)),
          HB_PCA=make_HB_formula(.x, covsPCA),
          HB_vars=formula(glue("{.x} ~ yday + lon + lat + siteid +",
                               "{paste(unlist(covs), collapse='+')}")),
          HB_vars_PCA=formula(glue("{.x} ~ yday + lon + lat + siteid +",
                                   "{paste(covsPCA, collapse='+')}"))
    )
  )

  cat("Starting", id, "for", y.i, ":", as.character(Sys.time()), "\n",
      file=glue("{log.dir}/{id}_{y.i}_{mod}.log"))

  for(r in responses) {

    set.seed(1003)
    folds_og_HB <- d.y$train[[r]] |>
      vfold_cv(strata=r)
    set.seed(1003)
    folds_PCA_HB <- dPCA.y$train[[r]] |>
      vfold_cv(strata=r)
    set.seed(1003)
    folds_og_ML <- d.y$train[[r]] |>
      select(-obsid, -y, -date, -year, -yday, -siteid, -lon, -lat) |>
      vfold_cv(strata=r)
    set.seed(1003)
    folds_PCA_ML <- dPCA.y$train[[r]] |>
      select(-obsid, -y, -date, -year, -yday, -siteid, -lon, -lat) |>
      vfold_cv(strata=r)

    if(mod == "HB") {
      # HB models --------------------------------------------------------------
      priStr <- switch(
        prior_strength,
        "1"=list(r1=0.2, r2=2, hs1=0.5, hs2=0.6, b=0.75, de=0.3, i=1),
        "2"=list(r1=0.1, r2=2, hs1=3, hs2=0.2, b=0.2, de=0.1, i=2),
        "3"=list(r1=0.1, r2=1, hs1=5, hs2=0.5, b=0.5, de=0.05, i=3)
      )
      opts <- list(
        iter=1500,
        warmup=1000,
        chains=3,
        cores=3,
        ctrl=list(adapt_delta=0.8, max_treedepth=10),
        prior_i=priStr$i
      )
      priors <- map(responses,
                    ~list(HB=make_HB_priors(priStr, "HB", .x, covs),
                          HBL_PCA=make_HB_priors(priStr, "HB", .x, covsPCA, PCA=T))
      )

      # fit models
      fit_candidate(mod, r, form.ls, d.y$train, opts, priors, fit.dir, y.i)
      fit_candidate(mod, r, form.ls, dPCA.y$train, opts, priors, fit.dir, y.i, "_PCA")

      if(HB_CV) {
        # run CV
        HB_run_CV("HB", folds_og_HB, cv.dir, y.i, y_i.i, r, form.ls, opts, priors, PCA=F)
        HB_run_CV("HB", folds_PCA_HB, cv.dir, y.i, y_i.i, r, form.ls, opts, priors, PCA=T)
      }
    } else {

      # ML models --------------------------------------------------------------
      if(.Platform$OS.type=="unix") {
        plan(multicore, workers=ncores)
      } else {
        plan(multisession, workers=ncores)
      }
      tunes <- list(nTuneVal) |> set_names(mod)
      fit_candidate(mod, r, form.ls, d.y$train, folds_og_ML, tunes, fit.dir, y.i)
      fit_candidate(mod, r, form.ls, dPCA.y$train, folds_PCA_ML, tunes, fit.dir, y.i, "_PCA")
      plan(sequential)
    }
  }

  cat("Finished", id, "for", y.i, ":", as.character(Sys.time()), "\n")
}
