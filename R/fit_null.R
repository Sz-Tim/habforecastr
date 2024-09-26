#' Identify rows in y within thresh days of x
#'
#' @param x
#' @param y
#' @param thresh
#'
#' @return
#' @export
get_rows_by_yday <- function(x, y, thresh) {
  xy_diff <- x - y
  comp.mx <- cbind(abs(xy_diff),
                   abs(xy_diff - 365),
                   abs(xy_diff + 365))
  ids <- which(apply(comp.mx, 1, function(j) any(j <= thresh)))
  return(ids)
}






#' Calculate null model predictions
#'
#' @param obs.ls
#' @param resp
#'
#' @return
#' @export
calc_null <- function(obs.ls, resp) {
  library(tidyverse)
  obs.df <- obs.ls[[resp]] |>
    mutate(yday=yday(date)) |>
    select(yday, {{resp}})
  if(resp=="alert") {
    temp.df <- tibble(yday=1:365,
                      null4wk_alert_A1=NA_real_)
    for(i in 1:nrow(temp.df)) {
      obs.ids <- get_rows_by_yday(i, obs.df$yday, 14)
      temp.df$null4wk_alert_A1[i] <- mean(obs.df$alert[obs.ids] == "A1")
    }
    out.df <- obs.ls[[resp]] |>
      mutate(yday=yday(date),
             nullGrand_alert_A1=mean(alert=="A1"))
  }
  if(resp=="tl") {
    temp.df <- tibble(yday=1:365,
                      null4wk_tl_TL0=NA_real_,
                      null4wk_tl_TL1=NA_real_,
                      null4wk_tl_TL2=NA_real_,
                      null4wk_tl_TL3=NA_real_)
    for(i in 1:nrow(temp.df)) {
      obs.ids <- get_rows_by_yday(i, obs.df$yday, 14)
      temp.df$null4wk_tl_TL0[i] <- mean(obs.df$tl[obs.ids] == "TL0")
      temp.df$null4wk_tl_TL1[i] <- mean(obs.df$tl[obs.ids] == "TL1")
      temp.df$null4wk_tl_TL2[i] <- mean(obs.df$tl[obs.ids] == "TL2")
      temp.df$null4wk_tl_TL3[i] <- mean(obs.df$tl[obs.ids] == "TL3")
    }
    out.df <- obs.ls[[resp]] |>
      mutate(yday=yday(date),
             nullGrand_tl_TL0=mean(tl=="TL0"),
             nullGrand_tl_TL1=mean(tl=="TL1"),
             nullGrand_tl_TL2=mean(tl=="TL2"),
             nullGrand_tl_TL3=mean(tl=="TL3"))
  }
  if(resp=="lnN") {
    temp.df <- tibble(yday=1:365,
                      null4wk_lnN_lnN=NA_real_)
    for(i in 1:nrow(temp.df)) {
      obs.ids <- get_rows_by_yday(i, obs.df$yday, 14)
      temp.df$null4wk_lnN_lnN[i] <- mean(obs.df$lnN[obs.ids])
    }
    out.df <- obs.ls[[resp]] |>
      mutate(yday=yday(date),
             nullGrand_lnN_lnN=mean(lnN))
  }
  return(list(yday.df=temp.df,
              obs.df=left_join(out.df, temp.df) |> select(-yday)))
}

