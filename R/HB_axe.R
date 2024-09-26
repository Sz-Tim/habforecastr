#' Adapted from butcher::axe_env for bayesian(engine='brms')
#'
#' @param x
#' @param verbose
#' @param ...
#'
#' @return
#' @export
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
