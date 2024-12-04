

#' Create formulas for Hierarchical Bayesian models
#' This function creates formulas for Hierarchical Bayesian models, incorporating specified covariates and spline interactions.
#'
#' @param resp A character string specifying the response variable.
#' @param covs A character vector of covariates to include in the model.
#' @param splinesInt A character string specifying the type of spline interactions to include. Options are "time", "space", or "both". Default is "both".
#'
#' @return A `brmsformula` object for the specified Hierarchical Bayesian model.
#' @export
#'
#' @examples
#' \dontrun{
#' library(brms)
#' resp <- "lnN"
#' covs <- c("cov1", "cov2", "cov3")
#' formula <- make_HB_formula(resp, covs, splinesInt = "both")
#' }
make_HB_formula <- function(resp, covs, splinesInt="both") {
  library(tidyverse); library(brms); library(glue)

  splines_int <- switch(splinesInt,
                        "time"="s(yday, bs=('cc'))",
                        "space"="s(lon, lat)",
                        "both"="t2(yday, lon, lat, bs=c('cc','ts'), d=c(1,2))")

  if(resp=="lnN") {
    return(bf(glue("{resp} ~ 1 + {paste(covs, collapse='+')}",
                   "+ {splines_int}",
                   "+ (1 + {paste(covs, collapse='+')} | siteid)"),
              glue("hu ~ 1 + {paste(covs, collapse='+')}",
                   "+ {splines_int}",
                   "+ (1 + {paste(covs, collapse='+')} | siteid)")))
  } else {
    return(bf(glue("{resp} ~ 1 + {paste(covs, collapse='+')}",
                   "+ {splines_int}",
                   "+ (1 + {paste(covs, collapse='+')} | siteid)")))
  }
}








#' Create priors for each Hierarchical Bayesian model
#'
#' This function creates priors for each Hierarchical Bayesian model based on the specified parameters.
#'
#' @param prior_i A list containing the prior parameters.
#' @param mod A character string specifying the model type.
#' @param resp A character string specifying the response variable.
#' @param covs A character vector of covariates to include in the model.
#' @param PCA A logical value indicating whether to use principal component analysis (PCA). Default is FALSE.
#'
#' @return A list of priors for the specified Hierarchical Bayesian model.
#' @export
#'
#' @examples
#' \dontrun{
#' library(brms)
#' prior_i <- list(hs1 = 0.5, hs2 = 0.6, r1 = 0.2, r2 = 2)
#' mod <- "HB"
#' resp <- "lnN"
#' covs <- c("cov1", "cov2", "cov3")
#' priors <- make_HB_priors(prior_i, mod, resp, covs)
#' }
make_HB_priors <- function(prior_i, mod, resp, covs, PCA=F) {
  library(tidyverse); library(brms)
  p <- c(prior(normal(0, 1), class="Intercept"),
         prior(normal(0, 0.1), class="sd"))
  if(resp=="tl") {
    p <- c(p,
           prior_string(glue("horseshoe({prior_i$hs1}, par_ratio={prior_i$hs2})"), class="b"))
  } else {
    p <- c(p,
           prior_string(glue("R2D2({prior_i$r1}, {prior_i$r2})"), class="b"))
  }
  if(resp=="lnN") {
    p <- c(p,
           prior_string(glue("R2D2({prior_i$r1}, {prior_i$r2})"), class="b", dpar="hu"))
  }
  return(p)
}
