

#' Calculate J, F-1, precision, and recall for a dataframe
#'
#' This is a function to pass to future_map to avoid capturing large objects
#' in the environment which makes it extraordinarily slow
#'
#' @param dat
#'
#' @return
#' @export
get_metrics <- function(dat) {
  tibble(F1=f_meas(dat, alert, pred, event_level="second")$.estimate,
         precision=precision(dat, alert, pred, event_level="second")$.estimate,
         recall=recall(dat, alert, pred, event_level="second")$.estimate,
         mcc=mcc(dat, alert, pred)$.estimate)
}









#' Compute classification metrics across probability thresholds
#'
#' @param L.df
#' @param prSteps
#'
#' @return
#' @export
compute_thresholds <- function(L.df, prMin=0, prMax=1, prSteps=0.1, byPrevAlert=F, cores=1) {
  library(tidyverse); library(yardstick); library(furrr)
  pred.df <- map_dfr(seq(prMin, prMax, by=prSteps),
                     ~L.df |> mutate(thresh=.x)) |>
    mutate(pred=if_else(prA1 < thresh, "A0", "A1") |> factor(levels=c("A0", "A1")))
  col_to_drop <- c("obsid", "siteid", "date", "prA1")
  if(!byPrevAlert) {
    col_to_drop <- c(col_to_drop, "prevAlert")
  }
  if(cores > 1) {
    plan(multisession, workers=cores)
  } else {
    plan(sequential)
  }
  metric.df <- pred.df |>
    select(-all_of(col_to_drop)) |>
    nest(dat=c(alert, pred)) |>
    ungroup() |>
    mutate(metrics=future_map(dat, get_metrics)) |>
    select(-dat) |>
    unnest(metrics) |>
    mutate(mcc=if_else(is.na(mcc), 0, mcc))
  return(metric.df)
}









#' Find minimum possible PR-AUC
#'
#' Following Boyd et al 2012
#'
#' @param df Dataframe
#' @param ... Unquoted grouping variables
#'
#' @return
#' @export
find_AUCPR_min <- function(df, ...) {
  df_aucpr_min <- df |>
    filter(model==levels(model)[1]) |>
    group_by(...) |>
    summarise(prop=mean(alert=="A1")) |>
    mutate(AUCPR_min=1+((1-prop)*log(1-prop))/prop) |>
    ungroup() |>
    select(-prop)
  return(left_join(df, df_aucpr_min))
}









#' Calculate pseudo-R2s
#'
#' Currently supports McFadden and the Veall-Zimmermann correction of the Aldrich-Nelson pseudo-R2
#'
#' @param dat.df
#' @param type
#'
#' @return
#' @export
calc_R2 <- function(dat.df, type="mf", ...) {
  # See Smith & McKenna. 2013. A comparison of logistic regression pseudo-R2
  #   indices. Multiple Linear Regression Viewpoints, 39(2)
  if(type=="mf") {
    return(
      dat.df |>
        group_by(..., model, PCA, covSet) |>
        summarise(LL=sum(dbinom(alert, 1, prA1, log=T))) |>
        group_by(...) |>
        arrange(y, model) |>
        mutate(R2=1 - LL/first(LL)) |>
        select(-LL)
    )
  }
  if(type=="vz") {
    return(
      dat.df |>
        group_by(..., model, PCA, covSet) |>
        summarise(LL=sum(dbinom(alert, 1, prA1, log=T)),
                  N=n(),
                  p=mean(alert)) |>
        group_by(...) |>
        arrange(y, model) |>
        mutate(R2=(-2*(first(LL)-LL))/(-2*(first(LL)-LL)+N) /
                 ( (-2*(p*log(p)+(1-p)*log(1-p)))/(1-2*(p*log(p)+(1-p)*log(1-p))) )) |>
        select(-N, -p, -LL)
    )
  }
}





#' Title
#'
#' @param df
#' @param y
#' @param type
#'
#' @return
#' @export
get_intervals <- function(df, y, type="hdci") {
  library(tidybayes)

  ci_fun <- switch(type,
                   "hdci"=hdci,
                   "qi"=qi)
  df |>
    summarise(mn=mean({{y}}),
              med=median({{y}}),
              L025=ci_fun({{y}}, 0.95)[,1],
              L05=ci_fun({{y}}, 0.9)[,1],
              L10=ci_fun({{y}}, 0.8)[,1],
              L15=ci_fun({{y}}, 0.7)[,1],
              L20=ci_fun({{y}}, 0.6)[,1],
              L25=ci_fun({{y}}, 0.5)[,1],
              L30=ci_fun({{y}}, 0.4)[,1],
              L35=ci_fun({{y}}, 0.3)[,1],
              L40=ci_fun({{y}}, 0.2)[,1],
              L45=ci_fun({{y}}, 0.1)[,1],
              L55=ci_fun({{y}}, 0.1)[,2],
              L60=ci_fun({{y}}, 0.2)[,2],
              L65=ci_fun({{y}}, 0.3)[,2],
              L70=ci_fun({{y}}, 0.4)[,2],
              L75=ci_fun({{y}}, 0.5)[,2],
              L80=ci_fun({{y}}, 0.6)[,2],
              L85=ci_fun({{y}}, 0.7)[,2],
              L90=ci_fun({{y}}, 0.8)[,2],
              L95=ci_fun({{y}}, 0.9)[,2],
              L975=ci_fun({{y}}, 0.95)[,2])
}
