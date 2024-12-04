#' Run cross-validation for Bayesian models
#'
#' This function performs cross-validation for Bayesian models, fitting the model to training data and evaluating it on test data for each fold.
#'
#' @param mod A character string specifying the model to fit (e.g., "HB").
#' @param folds A resampling object containing the cross-validation folds.
#' @param cv.dir A character string specifying the directory to save cross-validation results.
#' @param y A character string specifying the target variable.
#' @param y_i.i A data frame with information about the variables of interest.
#' @param r A character string specifying the response variable.
#' @param form.ls A list of formulas for the models.
#' @param HB.i A list of options for the hierarchical Bayesian model.
#' @param priors A list of priors for the Bayesian model.
#' @param PCA A logical value indicating whether to use principal component analysis (PCA). Default is FALSE.
#' @param reverse A logical value indicating whether to reverse the order of the folds. Default is FALSE.
#'
#' @return None. The function saves the cross-validation results to the specified directory.
#' @export
#'
#' @examples
#' \dontrun{
#' library(tidymodels)
#' folds <- vfold_cv(mtcars, v = 5)
#' cv.dir <- "path/to/cv_results"
#' y <- "mpg"
#' y_i.i <- data.frame(variable = "mpg")
#' r <- "response"
#' form.ls <- list(response = mpg ~ wt + cyl)
#' HB.i <- list(iter = 1000, warmup = 500, chains = 4)
#' priors <- list(prior(normal(0, 1), class = "b"))
#' HB_run_CV("HB", folds, cv.dir, y, y_i.i, r, form.ls, HB.i, priors)
#' }
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
