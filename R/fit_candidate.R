#' Wrapper to fit a model
#'
#' This function fits a specified model to the data, including machine learning and hierarchical Bayesian models.
#'
#' @param mod A character string specifying the model to fit (e.g., `"Ridge"`, `"ENet"`, `"RF"`, `"NN"`, `"MARS"`, `"Boost"`, `"lgbm"`, `"HB"`).
#' @param resp A character string specifying the response variable.
#' @param form.ls A list of formulas for the models.
#' @param d.ls A list of data frames for the models.
#' @param opts A list of options for model fitting, including resampling and control parameters.
#' @param tunes A list of tuning parameters for the models.
#' @param out.dir A character string specifying the output directory for saving model results.
#' @param y A character string specifying the target variable.
#' @param suffix A character string specifying a suffix for the model ID. Default is `NULL`.
#'
#' @return None. Fitted object is stored.
#' @export
#'
#' @examples
#' \dontrun{
#' library(tidymodels)
#' form.ls <- list(alert = list(ML = y ~ ., ML_PCA = y ~ .))
#' d.ls <- list(alert = data.frame(y = factor(c("A1", "A2")), x1 = rnorm(2), x2 = rnorm(2)))
#' opts <- vfold_cv(d.ls$alert, v = 5)
#' tunes <- list(Ridge = 10, ENet = 10, RF = 10, NN = 10, MARS = 10, Boost = 10, lgbm = 10)
#' out.dir <- "model_output"
#' fit_candidate("Ridge", "alert", form.ls, d.ls, opts, tunes, out.dir, "y")
#' }
fit_candidate <- function(mod, resp, form.ls, d.ls, opts, tunes, out.dir, y, suffix=NULL) {
  library(glue); library(tidymodels)
  dir.create(glue("{out.dir}/meta/"), showWarnings=F, recursive=T)
  dir.create(glue("{out.dir}/vi/"), showWarnings=F, recursive=T)
  PCA_run <- all(!is.null(suffix), grepl("PCA", suffix))
  # Fit ML models
  if(mod %in% c("Ridge", "ENet", "RF", "NN", "MARS", "Boost", "lgbm")) {
    fit_ID <- glue("{y}_{resp}_{mod}{ifelse(is.null(suffix),'',suffix)}")
    if(file.exists(glue("{out.dir}/{fit_ID}.rds"))) {
      cat("File already exists:", glue("{out.dir}/{fit_ID}.rds"), "\n")
      return()
    }
    mod.prefix <- ifelse(PCA_run, "PCA.", "")
    ML_form <- ifelse(PCA_run, "ML_PCA", "ML")
    ML_spec <- switch(mod,
                      Ridge=logistic_reg(penalty=tune(),
                                         mixture=0) |>
                        set_engine("glmnet") |> set_mode("classification"),
                      ENet=logistic_reg(penalty=tune(),
                                        mixture=tune()) |>
                        set_engine("glmnet") |> set_mode("classification"),
                      RF=rand_forest(trees=tune(),
                                     min_n=tune()) |>
                        set_engine("randomForest") |> set_mode("classification"),
                      NN=mlp(hidden_units=tune(),
                             penalty=tune(),
                             epochs=tune()) |>
                        set_engine("nnet", maxNWts=1e4) |> set_mode("classification"),
                      MARS=mars(num_terms=tune(),
                                prod_degree=tune()) |>
                        set_engine("earth") |> set_mode("classification"),
                      Boost=boost_tree(trees=tune(),
                                       tree_depth=tune(),
                                       min_n=tune(),
                                       learn_rate=tune(),
                                       loss_reduction=tune()) |>
                        set_engine("xgboost") |> set_mode("classification"),
                      lgbm=boost_tree(tree_depth=tune(),
                                      trees=tune(),
                                      learn_rate=tune(),
                                      mtry=tune(),
                                      min_n=tune(),
                                      loss_reduction=tune()) |>
                        set_engine("lightgbm") |> set_mode("classification")
    )
    avg_prec2 <- metric_tweak("avg_prec2", average_precision, event_level="second")
    pr_auc2 <- metric_tweak("pr_auc2", pr_auc, event_level="second")
    roc_auc2 <- metric_tweak("roc_auc2", roc_auc, event_level="second")
    wf <- workflow() |>
      add_model(ML_spec) |>
      add_formula(form.ls[[resp]][[ML_form]])
    set.seed(1003)
    if(mod=="lgbm") {
      out_tune <- wf |>
        tune_grid(resamples=opts,
                  grid=grid_latin_hypercube(
                    extract_parameter_set_dials(ML_spec) |>
                      update(mtry=finalize(mtry(), opts)),
                    size=tunes[[mod]]),
                  metrics=metric_set(roc_auc2, pr_auc2, avg_prec2),
                  control=control_grid(save_pred=T,
                                       parallel_over="everything"))
    } else {
      out_tune <- wf |>
        tune_grid(resamples=opts,
                  grid=grid_latin_hypercube(extract_parameter_set_dials(ML_spec),
                                            size=tunes[[mod]]),
                  metrics=metric_set(roc_auc2, pr_auc2, avg_prec2),
                  control=control_grid(save_pred=T,
                                       parallel_over="everything"))
      # saveRDS(out_tune |> butcher, glue("{out.dir}/meta/{fit_ID}_tune.rds"))
    }
    best <- select_best(out_tune, metric="avg_prec2")
    out_tune |>
      collect_predictions() |>
      filter(.config==best$.config) |>
      arrange(.row) |>
      mutate(obsid=d.ls[[resp]]$obsid[.row],
             y=y) |>
      select(y, obsid, .pred_A1) |>
      rename_with(~glue("{mod.prefix}{mod}_{resp}_A1"), .cols=".pred_A1") |>
      saveRDS(glue("{out.dir}/cv/{fit_ID}_CV.rds"))
    out <- wf |>
      finalize_workflow(best) |>
      fit(data=d.ls[[resp]] |>
            select(-obsid, -y, -date, -year, -yday, -siteid, -lon, -lat))
    out |>
      extract_fit_engine() |>
      vip::vi(scale=T) |>
      saveRDS(glue("{out.dir}/vi/{fit_ID}_vi.rds"))
    out <- out |>
      butcher()
  }

  # Fit Hierarchical Bayesian models
  if(mod == "HB") {
    library(brms)
    fit_ID <- glue("{y}_{resp}_{mod}{opts$prior_i}{ifelse(is.null(suffix),'',suffix)}")
    if(file.exists(glue("{out.dir}/{fit_ID}.rds"))) {
      cat("File already exists:", glue("{out.dir}/{fit_ID}.rds"), "\n")
      return()
    }
    HB.family <- switch(resp,
                        lnN=hurdle_lognormal,
                        tl=cumulative,
                        alert=bernoulli)
    HB_form <- ifelse(PCA_run, paste0(mod, "_PCA"), paste0(mod))
    HB_form_dummy <- ifelse(PCA_run, "HB_vars_PCA", "HB_vars")
    wf <- workflow() |>
      add_model(bayesian(mode="classification", engine="brms",
                         formula.override=bayesian_formula(form.ls[[resp]][[HB_form]]),
                         family=HB.family,
                         prior=tunes[[resp]][[HB_form]],
                         init=0,
                         iter=opts$iter,
                         warmup=opts$warmup,
                         control=opts$ctrl,
                         chains=opts$chains,
                         cores=opts$cores,
                         save_model=glue("{out.dir}/meta/{fit_ID}.stan")),
                formula=form.ls[[resp]]$HB_vars) |>
      add_recipe(recipe(d.ls[[resp]], formula=form.ls[[resp]][[HB_form_dummy]]))
    out <- wf |>
      fit(data=d.ls[[resp]])
    out |>
      extract_fit_engine() |>
      fixef() |>
      as_tibble(rownames="Variable") |>
      saveRDS(glue("{out.dir}/vi/{fit_ID}_vi.rds"))
    out <- out |>
      axe_env_bayesian() |>
      axe_env_bayesian()
  }
  saveRDS(out, glue("{out.dir}/{fit_ID}.rds"))
  cat("Saved ", y, "_", resp, "_", mod, " as ", out.dir, "*", suffix, "\n", sep="")
}
