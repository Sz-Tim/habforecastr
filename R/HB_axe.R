#' Adapted from butcher::axe_env for bayesian(engine='brms')
#'
#' This function adapts the `butcher::axe_env()` function for Bayesian models fitted with the `brms` engine, removing unnecessary environments to reduce the size of the model object.
#'
#' @param x A fitted Bayesian model object.
#' @param verbose A logical value indicating whether to print detailed information during the process. Default is FALSE.
#' @param ... Additional arguments (currently not used).
#'
#' @return The modified Bayesian model object with reduced environments.
#' @export
#'
#' @examples
#' \dontrun{
#' library(brms)
#' fit <- brm(mpg ~ wt + cyl, data = mtcars)
#' fit_reduced <- axe_env_bayesian(fit)
#' }
axe_env_bayesian <- function(x, verbose = FALSE, ...) {
  old <- x
  for(i in seq_along(x$fit$actions$model$spec$args)) {
    attr(x$fit$actions$model$spec$args[[i]], ".Environment") <- rlang::base_env()
  }
  for(i in seq_along(x$fit$fit$spec$args)) {
    attr(x$fit$fit$spec$args[[i]], ".Environment") <- rlang::base_env()
  }
  for(i in seq_along(x$fit$fit$spec$method$fit$args)) {
    attr(x$fit$fit$spec$method$fit$args[[i]], ".Environment") <- rlang::base_env()
  }
  attr(x$fit$actions$model$formula, ".Environment") <- rlang::base_env()

  return(x)
}
