#' Summarise and aggregate predictions for a set of models
#'
#' @param d.y
#' @param set
#' @param resp
#' @param fit.dir
#' @param y_i.i
#' @param suffix
#'
#' @return
#' @export
summarise_predictions <- function(d.y, dPCA.y, resp, fit.dir, y_i.i, suffix=NULL) {
  library(tidyverse); library(glue); library(tidymodels)
  fits.f <- dirf(fit.dir, glue("{y_i.i$abbr[1]}_{resp}.*{ifelse(is.null(suffix),'',suffix)}"))
  names(fits.f) <- str_split_fixed(str_split_fixed(fits.f, glue("{resp}_"), 2)[,2], "_|\\.", 2)[,1]
  fits <- map(fits.f, readRDS)

  fits.dPCA <- grep("PCA", fits.f)
  fits.d <- grep("PCA", fits.f, invert=T)

  if(length(fits.dPCA) > 0) {
    preds.dPCA <- imap_dfc(fits[fits.dPCA], ~get_predictions(.x, .y, resp, dPCA.y, y_i.i)) |>
      rename_with(~glue("PCA.{.x}"))
  } else {
    preds.dPCA <- NULL
  }

  if(length(fits.d) > 0) {
    preds.d <- imap_dfc(fits[fits.d], ~get_predictions(.x, .y, resp, d.y, y_i.i))
  } else {
    preds.d <- NULL
  }

  out.df <- d.y[[resp]] |>
    select(y, obsid, siteid, date, {{resp}})
  if(!is.null(preds.d)) {
    out.df <- out.df |> bind_cols(preds.d)
  }
  if(!is.null(preds.dPCA)) {
    out.df <- out.df |> bind_cols(preds.dPCA)
  }
  return(out.df)
}





#' Generate predictions for a model
#'
#' @param fit
#' @param mod
#' @param resp
#' @param d.df
#' @param y_i.i
#'
#' @return
#' @export
get_predictions <- function(fit, mod, resp, d.df, y_i.i) {
  library(tidyverse); library(glue); library(tidymodels); library(brms)

  if(grepl("HB", mod)) {
    pred <- parsnip::extract_fit_engine(fit) |>
      posterior_epred(d.df[[resp]], allow_new_levels=T) |>
      summarise_post_preds(resp, y_i.i)
  } else {
    pred_type <- ifelse(resp=="lnN", "raw", "prob")
    preds <- predict(fit, d.df[[resp]], pred_type)
    pred <- summarise_ML_preds(preds, resp, y_i.i)
  }
  pred.df <- as_tibble(pred) |>
    rename_with(~glue("{mod}_{resp}_{.x}"))
  return(pred.df)
}




#' Calculate mean predictions for ML models
#'
#' @param preds
#' @param resp
#' @param y_i.i
#'
#' @return
#' @export
summarise_ML_preds <- function(preds, resp, y_i.i) {
  if(resp=="alert") {
    pred <- cbind(A1=preds$.pred_A1)
  }
  if(resp=="tl") {
    thresh <- as.numeric(str_sub(y_i.i$tl_thresh, -1, -1)) + 1
    pred <- cbind(preds,
                  rowSums(preds[,thresh:4]))
    colnames(pred) <- c(paste0("TL", 0:3), "A1")
  }
  if(resp=="lnN") {
    thresh <- log1p(y_i.i$N_thresh)
    pred <- cbind(lnN=preds,
                  A1=preds>=thresh)
  }
  return(pred)
}




#' Calculate mean predictions for HB models
#'
#' @param post
#' @param resp
#' @param y_i.i
#'
#' @return
#' @export
summarise_post_preds <- function(post, resp, y_i.i) {
  if(resp=="alert") {
    pred <- cbind(A1=colMeans(post))
  }

  if(resp=="tl") {
    thresh <- as.numeric(str_sub(y_i.i$tl_thresh, -1, -1)) + 1
    pred <- cbind(apply(post, 2:3, mean),
                  colMeans(apply(post[,,(thresh):4], 1:2, sum)))
    colnames(pred) <- c(paste0("TL", 0:3), "A1")
  }

  if(resp=="lnN") {
    thresh <- log1p(y_i.i$N_thresh)
    pred <- cbind(lnN=colMeans(post),
                  A1=colMeans(post>=thresh))
  }
  return(pred)
}









#' Merge summarised predictions across all models into single dataframes
#'
#' @param files vector of full file path
#'
#' @return
#' @export
merge_pred_dfs <- function(files, CV=NULL) {
  f.df <- tibble(f=files,
                 covSet=str_split(files, "/") |>
                   map_chr(~grep("^d[0-9][0-9]?-", .x, value=T) |>
                             str_split_fixed("-", 2) |>
                             magrittr::extract(,1)))
  if(is.null(CV)) {
    map(1:nrow(f.df),
        ~readRDS(f.df$f[.x]) |>
          lapply(function(x) {x |> mutate(covSet=paste0(f.df$covSet[.x], "."))})) |>
      list_transpose() |>
      map_depth(2, ~.x |>
                  pivot_longer(ends_with("_A1"), names_to="model", values_to="prA1") |>
                  mutate(model=paste0(covSet, model)) |>
                  select(-covSet) |>
                  pivot_wider(names_from="model", values_from="prA1")) |>
      map(~reduce(.x, full_join))
  } else if(CV=="HB") {
    map_dfr(1:nrow(f.df),
            ~readRDS(f.df$f[.x]) |> mutate(covSet=paste0(f.df$covSet[.x], "."))) |>
      pivot_longer(ends_with("A1"), names_to="model", values_to="prA1") |>
      na.omit() |>
      mutate(model=paste0(covSet, model)) |>
      select(-covSet) |>
      pivot_wider(names_from="model", values_from="prA1")
  } else if(CV=="ML") {
    map(1:nrow(f.df),
        ~readRDS(f.df$f[.x]) |>
          mutate(covSet=paste0(f.df$covSet[.x], ".")) |>
          pivot_longer(ends_with("A1"), names_to="model", values_to="prA1") |>
          mutate(model=paste0(covSet, model)) |>
          select(-covSet) |>
          pivot_wider(names_from="model", values_from="prA1")) |>
      reduce(full_join)
  } else {
    stop("CV must be 'HB', 'ML', or NULL")
  }
}












#' Calculate model weights based on mean log loss
#'
#' @param cv.df
#' @param resp
#' @param wt.penalty
#'
#' @return
#' @export
calc_LL_wts <- function(cv.df, resp, wt.penalty=2) {
  library(yardstick)
  if(resp=="alert") {
    wt.df <- cv.df |>
      pivot_longer(ends_with("_A1"), names_to="model", values_to="pr") |>
      group_by(model) |>
      average_precision(pr, truth=alert, event_level="second") |>
      mutate(.estimate=log(.estimate))
  }
  if(resp=="tl") {
    wt.df <- cv.df |>
      pivot_longer(contains("_tl_TL"), names_to="ID", values_to="pr") |>
      mutate(model=str_split_fixed(ID, "_TL", 2)[,1],
             predCat=str_split_fixed(ID, "_tl_", 2)[,2]) |>
      select(-ID) |>
      pivot_wider(names_from=predCat, values_from=pr) |>
      group_by(model) |>
      mn_log_loss(TL0, TL1, TL2, TL3, truth=tl)
  }
  if(resp=="lnN") {
    wt.df <- cv.df |>
      pivot_longer(ends_with("_lnN_lnN"), names_to="model", values_to="pr") |>
      group_by(model) |>
      rmse(truth=lnN, estimate=pr)
  }
  return(wt.df |>
           ungroup() |>
           mutate(wt=(1/.estimate^wt.penalty)/sum(1/.estimate^wt.penalty)))
}







