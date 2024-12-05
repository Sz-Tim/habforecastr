#' Split circular buffer into NSEW quadrants
#'
#' This function splits a circular buffer into four quadrants: North, South, East, and West.
#'
#' @param sf An `sf` object representing the circular buffer.
#' @param radius Numeric value specifying the radius of the buffer in meters. Default is 110,000 meters.
#'
#' @return An `sf` object with the circular buffer split into NSEW quadrants.
#' @export
#' #'
#' @examples
#' \dontrun{
#' library(sf)
#' # Create a sample sf object
#' sf <- st_as_sf(data.frame(siteid = 1, geometry = st_sfc(st_point(c(0, 0)))), crs = 27700)
#' # Split the buffer into quadrants
#' quadrants <- split_to_NSEW(sf, radius = 110e3)
#' }
split_to_NSEW <- function(sf, radius=110e3) {
  library(tidyverse); library(sf); library(lwgeom)
  hub.df <- sf |>
    st_centroid() |>
    select(siteid, geometry)
  spoke.df <- hub.df |>
    st_buffer(dist=radius, nQuadSegs=2) |>
    st_cast("POINT") |>
    group_by(siteid) |>
    mutate(spoke.id=row_number()) |>
    filter(spoke.id %% 2 == 0) |>
    mutate(side=c("start", "start", "end", "end")) |>
    ungroup()
  coords <- cbind(st_coordinates(filter(spoke.df, side=="start")),
                  st_coordinates(filter(spoke.df, side=="end")))
  spoke.lines <- map(1:nrow(coords),
                     ~st_linestring(matrix(coords[.x,],ncol=2, byrow=T))) |>
    st_sfc() |> st_as_sf() |> st_set_crs(27700) |>
    rename(geometry=x) |>
    mutate(siteid=filter(spoke.df, side=="start")$siteid)
  sf.quad <- sf$siteid |>
    map_dfr(~st_split(filter(sf, siteid==.x), filter(spoke.lines, siteid==.x)) |>
              st_collection_extract() |>
              mutate(quadrant=c("E", "S", "W", "N")))

  return(sf.quad)
}





#' Lag multiple variables at once
#'
#' This function creates multiple lagged versions of specified variables in a data frame.
#'
#' @param data A data frame.
#' @param ... Unquoted variable names to lag.
#' @param n An integer specifying the number of lags (default is 2).
#'
#' @return A data frame with the original variables and their lagged versions.
#' @export
#'
#' @examples
#' # Example usage:
#' # df <- data.frame(time = 1:10, value = rnorm(10))
#' # get_lags(df, value, n = 3)
get_lags <- function(data, ..., n=2){
  library(tidyverse); library(rlang)
  variable <- enquos(...)

  indices <- seq_len(n)
  combos <- tidyr::crossing(indices, var=as.list(variable))

  quosures <- map2(combos$indices, combos$var,
                   ~quo(lag(!!.y, !!.x)) ) |>
    set_names(paste0(map_chr(combos$var, quo_text), combos$indices))
  mutate(data, !!!quosures )
}





#' Identify traffic light and alert status based on density/concentration
#'
#' This function assigns traffic light and alert statuses to each row in a data frame based on specified density or concentration thresholds.
#'
#' @param y.df A data frame containing the data to be evaluated.
#' @param N A numeric vector representing the density or concentration values.
#' @param tl_i A data frame containing the traffic light (`tl`), alert (`A`), abbreviation (`abbr`), and minimum threshold (`min_ge`) values.
#'
#' @return A data frame with additional columns for traffic light (`tl`) and alert status (`alert`).
#' @export
#' #'
#' @examples
#' \dontrun{
#' y.df <- data.frame(y = c("A", "B", "C"), N = c(10, 20, 30))
#' tl_i <- data.frame(tl = c("Green", "Yellow", "Red"), A = c("Low", "Medium", "High"), abbr = c("A", "B", "C"), min_ge = c(5, 15, 25))
#' result <- get_trafficLights(y.df, y.df$N, tl_i)
#' }
get_trafficLights <- function(y.df, N, tl_i) {
  library(tidyverse)
  y.df |>
    rowwise() |>
    mutate(tl = tl_i$tl[max(which(y == tl_i$abbr & {{N}} >= tl_i$min_ge))],
           alert = tl_i$A[max(which(y == tl_i$abbr & {{N}} >= tl_i$min_ge))]) |>
    ungroup()
}





#' Calculate local and regional autoregressive terms for HAB and toxin observations
#'
#' This function calculates local and regional autoregressive terms for Harmful Algal Bloom (HAB) and toxin observations.
#'
#' @param yRaw.df A data frame containing raw observations of HAB and toxins.
#' @param y_i A data frame with information about the variables of interest, including abbreviations.
#' @param tl_i A data frame containing traffic light (`tl`), alert (`A`), abbreviation (`abbr`), and minimum threshold (`min_ge`) values.
#' @param site.100km A data frame with site information within a 100 km radius.
#' @param dist.df A data frame containing distance information between sites. Default is `NULL`.
#' @param forecastStart A character string specifying the start date for forecasting. Default is `"3000-01-01"`.
#'
#' @return A data frame with calculated autoregressive terms and additional features for each observation.
#' @export
#' #'
#' @examples
#' \dontrun{
#' yRaw.df <- data.frame(siteid = rep(1:3, each = 10), date = rep(seq.Date(Sys.Date(), by = "day", length.out = 10), 3), N = rnorm(30))
#' y_i <- data.frame(abbr = c("HAB", "Toxin"))
#' tl_i <- data.frame(tl = c("Green", "Yellow", "Red"), A = c("Low", "Medium", "High"), abbr = c("HAB", "Toxin"), min_ge = c(0, 1, 2))
#' dist.df <- data.frame(origins = rep(1:3, each = 3), dest_c = list(1:3, 1:3, 1:3))
#' result <- calc_y_features(yRaw.df, y_i, tl_i, dist.df)
#' }
calc_y_features <- function(yRaw.df, y_i, tl_i, dist.df=NULL, forecastStart="3000-01-01") {
  y.ls <- yRaw.df |>
    group_by(siteid, date) |>
    slice_head(n=1) |>
    ungroup() |>
    pivot_longer(any_of(y_i$abbr), names_to="y", values_to="N") |>
    filter((!is.na(N)) | (date >= forecastStart)) |>
    mutate(lnN=log1p(N)) |>
    get_trafficLights(N, tl_i) |>
    arrange(y, siteid, date) |>
    group_by(y, siteid) |>
    get_lags(lnN, alert, date, n=2) |>
    ungroup() |>
    mutate(year=year(date),
           lnNWt1=lnN1/log1p(as.numeric(date-date1)),
           lnNWt2=lnN2/log1p(as.numeric(date-date2)),
           lnNAvg1=0, lnNAvg2=0, prAlertAvg1=0, prAlertAvg2=0,
           lnNAvgPrevYr=0, prAlertAvgPrevYr=0) |>
    group_by(y) |>
    group_split()
  y_new <- vector("list", length(y.ls))

  for(i in seq_along(y.ls)) {
    y.df_i <- y.ls[[i]] |> select(siteid, date, lnN1, lnN2, alert1, alert2)
    yYr_i <- y.ls[[i]] |>
      group_by(siteid, year) |>
      summarise(lnN_yr=mean(lnN, na.rm=T),
                prAlert_yr=mean(alert=="A1", na.rm=T)) |>
      ungroup() |>
      mutate(year_matchObs=year+1)
    y_new[[i]] <- y.ls[[i]] |>
      left_join(yYr_i |> select(year_matchObs, siteid, lnN_yr, prAlert_yr),
                by=c("year"="year_matchObs", "siteid"="siteid")) |>
      rename(lnNPrevYr=lnN_yr, prAlertPrevYr=prAlert_yr)

    for(j in 1:nrow(y.ls[[i]])) {
      site_j <- y.df_i$siteid[j]
      date_j <- y.df_i$date[j]
      yr_j <- year(date_j)
      wk.df <- y.df_i |>
        filter(siteid %in% dist.df$dest_c[dist.df$origins==site_j][[1]] &
                 date <= date_j &
                 date > date_j-7)
      yr.df <- yYr_i |>
        filter(siteid %in% dist.df$dest_c[dist.df$origins==site_j][[1]] &
                 year == yr_j - 1)
      yr.site_j <- yYr_i |> filter(siteid==site_j & year==yr_j-1)
      y_new[[i]]$lnNAvg1[j] <- mean(wk.df$lnN1, na.rm=T)
      y_new[[i]]$lnNAvg2[j] <- mean(wk.df$lnN2, na.rm=T)
      y_new[[i]]$prAlertAvg1[j] <- mean(wk.df$alert1 == "A1", na.rm=T)
      y_new[[i]]$prAlertAvg2[j] <- mean(wk.df$alert2 == "A1", na.rm=T)
      y_new[[i]]$lnNAvgPrevYr[j] <- mean(yr.df$lnN_yr, na.rm=T)
      y_new[[i]]$prAlertAvgPrevYr[j] <- mean(yr.df$prAlert_yr, na.rm=T)
      if(j %% 1000 == 0) {cat(i, ":", j, "of", nrow(y_new[[i]]), "\n")}
    }
  }
  y.df <- do.call('rbind', y_new)
  return(y.df)
}




#' Summarise HAB states for each toxin observation
#'
#' This function summarises Harmful Algal Bloom (HAB) states for each toxin observation by calculating the mean values of HAB-related variables within a specified time window.
#'
#' @param site_tox.sf An `sf` object representing the toxin observation sites.
#' @param site_hab.sf An `sf` object representing the HAB observation sites.
#' @param tox.obs A data frame containing toxin observations.
#' @param hab.df A data frame containing HAB observations.
#'
#' @return A data frame with toxin observations and summarised HAB states.
#' @export
#' #'
#' @examples
#' \dontrun{
#' library(sf)
#' # Create sample sf objects for toxin and HAB sites
#' site_tox.sf <- st_as_sf(data.frame(siteid = 1:3, geom = st_sfc(st_point(c(0, 0)), st_point(c(1, 1)), st_point(c(2, 2)))), crs = 4326)
#' site_hab.sf <- st_as_sf(data.frame(siteid = 4:6, geom = st_sfc(st_point(c(0.5, 0.5)), st_point(c(1.5, 1.5)), st_point(c(2.5, 2.5)))), crs = 4326)
#' # Create sample data frames for toxin and HAB observations
#' tox.obs <- data.frame(siteid = 1:3, date = Sys.Date() - 1:3)
#' hab.df <- data.frame(siteid = 4:6, date = Sys.Date() - 1:6, y = rep("HAB", 6), lnN = rnorm(6), alert = rep("A1", 6))
#' # Summarise HAB states for each toxin observation
#' result <- summarise_hab_states(site_tox.sf, site_hab.sf, tox.obs, hab.df)
#' }
summarise_hab_states <- function(site_tox.sf, site_hab.sf, tox.obs, hab.df) {
  library(tidyverse); library(sf)

  # identify hab sites within each toxin site buffer
  hab_ids <- site_tox.sf |>
    select(siteid, geom) %>%
    mutate(hab_id=st_intersects(., site_hab.sf)) |>
    st_drop_geometry()

  # match hab observations to toxin observations
  hab.df <- hab.df |>
    mutate(prA=as.numeric(alert=="A1")) |>
    select(siteid, date, y, lnN, prA) |>
    rename(hab_id=siteid, lnNAvg=lnN) |>
    pivot_wider(names_from=y, values_from=c(lnNAvg, prA), names_glue="{y}{.value}")
  hab_y_names <- names(hab.df)[-(1:2)]
  habSums <- vector("list", nrow(tox.obs))
  for(i in 1:nrow(tox.obs)) {
    hab_sites <- filter(hab_ids, siteid==tox.obs$siteid[i])$hab_id[[1]]
    habSums[[i]] <- hab.df |>
      filter(hab_id %in% hab_sites &
               date <= tox.obs$date[i] - 7*1 &
               date >= tox.obs$date[i] - 7*5) |>
      summarise(across(all_of(hab_y_names), ~mean(.x, na.rm=T)))
  }
  out.df <- tox.obs |> bind_cols(do.call('rbind', habSums))
  return(out.df)
}





#' Detrend observations using a loess smoother
#'
#' This function calculates residuals from a loess smoother to detrend observations.
#'
#' @param x A numeric vector representing the predictor variable.
#' @param y A numeric vector representing the response variable.
#' @param span A numeric value controlling the degree of smoothing. Default is 0.75.
#' @param robust A logical value indicating whether to use robust fitting. Default is TRUE.
#'
#' @return A numeric vector of residuals from the loess smoother.
#' @export
#'
#' @examples
#' \dontrun{
#' x <- 1:100
#' y <- sin(x / 10) + rnorm(100)
#' detrended_y <- detrend_loess(x, y, span = 0.5, robust = FALSE)
#' plot(x, y)
#' lines(x, detrended_y, col = "red")
#' }
detrend_loess <- function (x, y, span=0.75, robust=TRUE) {
  # modified from astsa::trend
  if(sum(!is.na(y)) < 10) {
    return(y)
  }
  if(length(y) < 10) {
    return(y)
  }
  fam <- ifelse(robust, "symmetric", "gaussian")
  lo <- stats::predict(stats::loess(y ~ x, span=span, family=fam), se = F)
  return(c(y - lo))
}
