#' Split circular buffer into NSEW quadrants
#'
#' @param sf
#' @param radius
#'
#' @return
#' @export
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
#' From https://stackoverflow.com/questions/55814028/multiple-lags-with-dplyr
#'
#' @param data Dataframe
#' @param ... Unquoted variable names to lag
#' @param n Number of lags
#'
#' @return
#' @export
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
#' @param y.df
#' @param N
#' @param tl_i
#'
#' @return
#' @export
get_trafficLights <- function(y.df, N, tl_i) {
  library(tidyverse)
  y.df |>
    rowwise() |>
    mutate(tl=tl_i$tl[max(which(y==tl_i$abbr & {{N}} >= tl_i$min_ge))],
           alert=tl_i$A[max(which(y==tl_i$abbr & {{N}} >= tl_i$min_ge))]) |>
    ungroup()
}





#' Calculate local and regional autoregressive terms for HAB and toxin observations
#'
#' @param yRaw.df
#' @param y_i
#' @param tl_i
#' @param site.100km
#'
#' @return
#' @export
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
#' @param site_tox.sf
#' @param site_hab.sf
#' @param tox.obs
#' @param hab.df
#'
#' @return
#' @export
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
#' @param x
#' @param y
#' @param span
#' @param robust
#'
#' @return
#' @export
detrend_loess <- function (x, y, span=0.75, robust=TRUE) {
  # modified from astsa::trend
  if(sum(!is.na(y)) < 10) {
    return(y)
  }
  if(length(y) < 10) {
    return(y)
  }
  fam = ifelse(robust, "symmetric", "gaussian")
  lo = stats::predict(stats::loess(y ~ x, span=span, family=fam), se = F)
  return(c(y - lo))
}
