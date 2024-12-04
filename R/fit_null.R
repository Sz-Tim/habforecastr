#' Identify rows in y within thresh days of x
#'
#' This function identifies rows where the absolute difference between the day of the year values in `x` and `y` is within a specified threshold.
#'
#' @param x A numeric vector representing day of the year values.
#' @param y A numeric vector representing day of the year values to compare against.
#' @param thresh A numeric value specifying the threshold for the absolute difference.
#'
#' @return A numeric vector of row indices where the absolute difference is within the threshold.
#' @export
#'
#' @examples
#' \dontrun{
#' x <- c(10, 50, 100, 200, 300)
#' y <- c(15, 55, 105, 205, 305)
#' thresh <- 5
#' result <- get_rows_by_yday(x, y, thresh)
#' }
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
#' This function calculates null model predictions for different response variables by averaging observations within a specified time window.
#'
#' @param obs.ls A list of data frames containing observations for different response variables.
#' @param resp A character string specifying the response variable (e.g., `"alert"`, `"tl"`, `"lnN"`).
#'
#' @return A list containing two data frames: `yday.df` with null model predictions for each day of the year, and `obs.df` with the original observations and corresponding null model predictions.
#' @export
#'
#' @examples
#' \dontrun{
#' library(tidyverse)
#' obs.ls <- list(alert = data.frame(date = Sys.Date() - 1:10, alert = sample(c("A1", "A2"), 10, replace = TRUE)))
#' result <- calc_null(obs.ls, "alert")
#' }
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

