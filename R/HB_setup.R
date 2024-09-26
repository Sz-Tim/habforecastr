

#' Create formulas for Hierarchical Bayesian models
#'
#' @param resp
#' @param covs
#' @param splinesInt
#'
#' @return
#' @export
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
#' @param prior_i
#' @param mod
#' @param resp
#' @param covs
#'
#' @return
#' @export
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
