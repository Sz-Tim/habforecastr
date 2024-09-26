#' Run cross-validation for Bayesian models
#'
#' @param mod
#' @param folds
#' @param cv.dir
#' @param y
#' @param y_i.i
#' @param r
#' @param form.ls
#' @param HB.i
#' @param priors
#'
#' @return
#' @export
HB_run_CV <- function(mod, folds, cv.dir, y, y_i.i, r, form.ls, HB.i, priors, PCA=F, reverse=F) {
  for(f in 1:nrow(folds)) {
    if(reverse) {
      f <- (nrow(folds):1)[f]
    }
    f_ <- paste0("_f", str_pad(f, 2, side="left", pad="0"))
    f_ <- ifelse(PCA, paste0("_PCA", f_), f_)
    d.cv <- list(train=list(alert=training(folds$splits[[f]])),
                 test=list(alert=testing(folds$splits[[f]])))
    if(!file.exists(glue("{cv.dir}/{y}_{r}_{mod}_CV{f_}.rds"))) {
      fit_candidate(mod, r, form.ls, d.cv$train, HB.i, priors, cv.dir, y, f_)
    }
    if(file.exists(glue("{cv.dir}/{y}_{r}_{mod}1{f_}.rds"))) {
      summarise_predictions(d.cv$test, d.cv$test, r, cv.dir, y_i.i, glue("{mod}1{f_}")) |>
        saveRDS(glue("{cv.dir}/{y}_{r}_{mod}_CV{f_}.rds"))
      file.remove(glue("{cv.dir}/{y}_{r}_{mod}1{f_}.rds"))
    }
  }
}
