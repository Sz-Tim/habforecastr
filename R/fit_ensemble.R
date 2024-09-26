#' Calculate ensemble model predictions
#'
#' @param out.ls
#' @param wt.ls
#' @param resp
#' @param y_i.i
#'
#' @return
#' @export
fit_ensemble <- function(out.ls, wt.ls, resp, y_i.i, method="wtmean", out.path=NULL, opt=NULL) {
  library(tidyverse); library(tidymodels)
  if(grepl("RF|GLM|HB", method) & resp != "alert") {
    stop("Models only implemented for alert")
  }
  if(resp=="alert") {
    if(method=="wtmean") {
      out <- left_join(
        out.ls[[resp]] |>
          pivot_longer(ends_with("_A1"), names_to="model", values_to="pr"),
        wt.ls[[resp]]
      ) |>
        mutate(pr_logit=brms::logit_scaled(pr, lb=-0.01, ub=1.01)) |>
        group_by(obsid) |>
        summarise(ens_alert_A1=sum(pr*wt, na.rm=T),
                  ensLogitMn_alert_A1=brms::inv_logit_scaled(
                    sum(pr_logit*wt, na.rm=T))) |>
        ungroup()
    } else if(grepl("[GLM|RF|HB]_fit", method)) {
      avg_prec2 <- metric_tweak("avg_prec2", average_precision, event_level="second")
      folds <- vfold_cv(wt.ls[[resp]], strata="alert")
      # folds <- group_vfold_cv(wt.ls[[resp]], group=year)
      ens_rec <- recipe(alert~., data=wt.ls[[resp]]) |>
        update_role(y, obsid, siteid, date, year, new_role="ID") |>
        step_logit(ends_with("_A1"), offset=0.01) |>
        step_normalize(all_predictors())
      ens_rec2 <- recipe(alert~., data=wt.ls[[resp]]) |>
        update_role(y, obsid, siteid, date, year, new_role="ID")

      if(method=="GLM_fit") {
        size <- ifelse(is.null(opt), 1e3, opt)
        GLM_spec <- logistic_reg(penalty=tune(), mixture=0) |>
          set_engine("glmnet", lower.limits=0) |> set_mode("classification")

        GLM_wf <- workflow() |>
          add_model(GLM_spec) |>
          add_recipe(ens_rec)
        GLM_tune <- GLM_wf |>
          tune_grid(resamples=folds,
                    grid=grid_latin_hypercube(extract_parameter_set_dials(GLM_spec),
                                              size=size),
                    metrics=metric_set(avg_prec2))
        GLM_out <- GLM_wf |>
          finalize_workflow(select_best(GLM_tune, metric="avg_prec2")) |>
          fit(wt.ls[[resp]]) |>
          butcher()
        saveRDS(GLM_out, glue("{out.path}/{y_i.i$abbr}_EnsGLM.rds"))

        GLM_wf2 <- workflow() |>
          add_model(GLM_spec) |>
          add_recipe(ens_rec2)
        GLM_tune2 <- GLM_wf2 |>
          tune_grid(resamples=folds,
                    grid=grid_latin_hypercube(extract_parameter_set_dials(GLM_spec),
                                              size=size),
                    metrics=metric_set(avg_prec2))
        GLM_out2 <- GLM_wf2 |>
          finalize_workflow(select_best(GLM_tune2, metric="avg_prec2")) |>
          fit(wt.ls[[resp]]) |>
          butcher()
        saveRDS(GLM_out2, glue("{out.path}/{y_i.i$abbr}_EnsGLM2.rds"))
      }

      if(method=="RF_fit") {
        size <- ifelse(is.null(opt), 1e2, opt)
        RF_spec <- rand_forest(trees=tune(),
                               min_n=tune()) |>
          set_engine("randomForest") |> set_mode("classification")

        RF_wf <- workflow() |>
          add_model(RF_spec) |>
          add_recipe(ens_rec)
        RF_tune <- RF_wf |>
          tune_grid(resamples=folds,
                    grid=grid_latin_hypercube(extract_parameter_set_dials(RF_spec),
                                              size=size),
                    metrics=metric_set(avg_prec2))
        RF_out <- RF_wf |>
          finalize_workflow(select_best(RF_tune, metric="avg_prec2")) |>
          fit(wt.ls[[resp]]) |>
          butcher()
        saveRDS(RF_out, glue("{out.path}/{y_i.i$abbr}_EnsRF.rds"))

        RF_wf2 <- workflow() |>
          add_model(RF_spec) |>
          add_recipe(ens_rec2)
        RF_tune2 <- RF_wf2 |>
          tune_grid(resamples=folds,
                    grid=grid_latin_hypercube(extract_parameter_set_dials(RF_spec),
                                              size=size),
                    metrics=metric_set(avg_prec2))
        RF_out2 <- RF_wf2 |>
          finalize_workflow(select_best(RF_tune2, metric="avg_prec2")) |>
          fit(wt.ls[[resp]]) |>
          butcher()
        saveRDS(RF_out2, glue("{out.path}/{y_i.i$abbr}_EnsRF2.rds"))
      }

      if(method=="HB_fit") {
        mn <- ifelse(is.null(opt), 0.5, opt)
        library(brms)
        wf <- workflow() |>
          add_model(bayesian(mode="classification", engine="brms",
                             family=bernoulli,
                             prior=prior_string(glue("exponential({mn})"), class="b", lb=0),
                             init=0,
                             iter=1500,
                             warmup=1000,
                             chains=3,
                             cores=3,
                             save_model=glue("{out.path}/{y_i.i$abbr}_EnsHB1.stan")),
                    formula=alert~.) |>
          add_recipe(ens_rec2)
        HB_out <- wf |>
          fit(data=wt.ls[[resp]])
        saveRDS(HB_out |> axe_env_bayesian() |> axe_env_bayesian(),
                glue("{out.path}/{y_i.i$abbr}_EnsHB.rds"))
      }
    }
    if(grepl("GLM", method)) {
      GLM_out <- readRDS(glue("{out.path}/{y_i.i$abbr}_EnsGLM.rds"))
      GLM_out2 <- readRDS(glue("{out.path}/{y_i.i$abbr}_EnsGLM2.rds"))
      out <- out.ls[[resp]] %>%
        mutate(ensGLM_alert_A1=predict(GLM_out, new_data=., type="prob")[[2]],
               ensGLM2_alert_A1=predict(GLM_out2, new_data=., type="prob")[[2]]) |>
        select(obsid, starts_with("ensGLM"))
    }
    if(grepl("RF", method)) {
      RF_out <- readRDS(glue("{out.path}/{y_i.i$abbr}_EnsRF.rds"))
      RF_out2 <- readRDS(glue("{out.path}/{y_i.i$abbr}_EnsRF2.rds"))
      out <- out.ls[[resp]] %>%
        mutate(ensRF_alert_A1=predict(RF_out, new_data=., type="prob")[[2]],
               ensRF2_alert_A1=predict(RF_out2, new_data=., type="prob")[[2]]) |>
        select(obsid, starts_with("ensRF"))
    }
    if(grepl("HB", method)) {
      HB_out <- readRDS(glue("{out.path}/{y_i.i$abbr}_EnsHB.rds"))
      pred <- parsnip::extract_fit_engine(HB_out) |>
        posterior_epred(out.ls[[resp]], allow_new_levels=T) |>
        summarise_post_preds(resp, y_i.i)
      out <- out.ls[[resp]] %>%
        mutate(ensHB_alert_A1=pred[,1]) |>
        select(obsid, starts_with("ensHB"))
    }
  }
  if(resp=="tl") {
    thresh <- as.numeric(str_sub(y_i.i$tl_thresh, -1, -1))
    alert_cols <- paste0("ens_tl_TL", thresh:3)
    out <- left_join(
      out.ls[[resp]] |>
        select(-ends_with("_A1")) |>
        pivot_longer(contains("_tl_TL"), names_to="ID", values_to="pr") |>
        mutate(model=str_split_fixed(ID, "_TL", 2)[,1],
               predCat=str_split_fixed(ID, "_tl_", 2)[,2]) |>
        select(-ID),
      wt.ls[[resp]]) |>
      group_by(obsid, predCat) |>
      summarise(ens_tl=sum(pr*wt, na.rm=T)) |>
      ungroup() |>
      pivot_wider(names_from=predCat, values_from=ens_tl, names_prefix="ens_tl_") |>
      mutate(ens_tl_A1=rowSums(pick(all_of(alert_cols))))
  }
  if(resp=="lnN") {
    thresh <- log1p(y_i.i$N_thresh)
    out <- left_join(
      out.ls[[resp]] |>
        select(-ends_with("_A1")) |>
        pivot_longer(ends_with("_lnN_lnN"), names_to="model", values_to="pr"),
      wt.ls[[resp]]) |>
      group_by(obsid) |>
      summarise(ens_lnN=sum(pr*wt, na.rm=T),
                ens_lnN_A1=sum((pr >= thresh)*wt, na.rm=T)) |>
      ungroup()
  }
  return(left_join(out.ls[[resp]], out))
}
