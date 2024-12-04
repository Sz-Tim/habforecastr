#' Find id of environmental layer points nearest to each site point location
#'
#' This function finds the ID of the environmental layer points that are nearest to each site point location.
#'
#' @param site.df A data frame with site locations, containing columns `lon` and `lat`.
#' @param env.sf An `sf` object with environmental features.
#' @param id_env A character string specifying the ID column in `env.sf`.
#'
#' @return A data frame with the nearest environmental feature IDs added.
#' @export
find_nearest_feature_id <- function(site.df, env.sf, id_env) {
  site.df |>
    st_as_sf(coords=c("lon", "lat"), crs=27700, remove=F) |>
    st_transform(4326) %>%
    mutate(new_id=st_nearest_feature(., env.sf)) |>
    rename_with(~id_env, new_id) |>
    st_drop_geometry()
}




#' Find ids of environmental layer points within each site buffer
#'
#' This function finds the IDs of environmental layer points that fall within each site buffer.
#'
#' @param site.sf An `sf` object with site quadrants.
#' @param env.sf An `sf` object with environmental features.
#' @param id_env A character string specifying the ID column in `env.sf`.
#'
#' @return A data frame with the intersecting environmental feature IDs added.
#' @export
find_buffer_intersect_ids <- function(site.sf, env.sf, id_env) {
  site.sf |>
    select(siteid, quadrant, geom) |>
    st_transform(4326) |>
    st_make_valid() %>%
    mutate(new_id=st_intersects(., env.sf)) |>
    rename_with(~id_env, new_id) |>
    st_drop_geometry()
}


#' Extract CMEMS or WRF data to site point locations
#'
#' This function extracts CMEMS or WRF data to site point locations.
#'
#' @param site.df A data frame with site locations, containing columns `lon` and `lat`.
#' @param env_vars A character vector specifying the column names of environmental variables.
#' @param env.df A data frame of environmental variables.
#' @param id_env A character string specifying the ID column in `env.df`.
#' @param site.v A character string specifying the version for site alignment.
#'
#' @return A data frame with extracted environmental data for each site point location.
#' @export
extract_env_pts <- function(site.df, env_vars, env.df, id_env, site.v) {
  library(tidyverse); library(zoo)

  env.site <- env.df |>
    group_by(version) |>
    group_split() |>
    map2_dfr(.x=_, .y=site.v, ~.x |> filter({{id_env}} %in% site.df[[.y]])) |>
    arrange({{id_env}}, date) |>
    group_by({{id_env}}) |>
    mutate(across(any_of(env_vars),
                  ~rollmeanr(.x, k=7, na.pad=T, align="right"),
                  .names="{.col}Wk")) |>
    mutate(across(any_of(paste0(env_vars, "Wk")),
                  ~.x - lag(.x),
                  .names="{.col}Delta")) |>
    # ungroup() |>
    # mutate(yday=yday(date)) |>
    # group_by({{id_env}}) |>
    # mutate(across(any_of(env_vars),
    #               ~detrend_loess(yday, .x, span=0.3),
    #               .names="{.col}Dt")) |>
    ungroup()
  env.site <- env.site |>
    select({{id_env}}, version, date,
           any_of(paste0(env_vars, "Wk")),
           any_of(paste0(env_vars, "WkDelta")),
           any_of(paste0(env_vars, "Dt")))
  return(env.site)
}







#' Extract CMEMS or WRF data within site buffer quadrants
#'
#' This function extracts environmental data from CMEMS or WRF datasets within specified site buffer quadrants.
#'
#' @param site.buffer sf object with site quadrants.
#' @param vars List containing column names for environmental variables. Should include `all` and optionally `sea`.
#' @param env.df Data frame containing environmental data with a `date` column.
#' @param id_env Character vector of column names in `site.buffer` used to match with `env.df`.
#'
#' @return A data frame with environmental data aggregated within site buffer quadrants.
#' @export
#' #'
#' @examples
#' \dontrun{
#' site.buffer <- st_read("path/to/site_buffer.shp")
#' vars <- list(all = c("temp", "salinity"), sea = c("sst"))
#' env.df <- read.csv("path/to/env_data.csv")
#' id_env <- c("siteid", "quadrant")
#' result <- extract_env_buffers(site.buffer, vars, env.df, id_env)
#' }
extract_env_buffers <- function(site.buffer, vars, env.df, id_env) {

  library(tidyverse); library(zoo)

  env.df <- env.df |> arrange(date, pick(any_of(id_env)))
  env.buffer <- expand_grid(siteid=unique(site.buffer$siteid),
                            quadrant=unique(site.buffer$quadrant),
                            date=unique(env.df$date)) |>
    mutate(v=1 + ((date >= "2019-04-01")*(length(id_env)>1))) |>
    bind_cols(as_tibble(setNames(map(vars$all, ~NA_real_), vars$all)))

  env_id.ls <- map(id_env, ~map(site.buffer[[.x]], ~.x))
  env.df$date_id <- match(env.df$date, unique(env.df$date))
  env_dates.ls <- split(1:nrow(env.df), env.df$date_id)

  ij <- 1
  startTime <- Sys.time()
  if(is.null(vars$sea)) {
    for(i in 1:nrow(site.buffer)) {
      for(j in 1:length(env_dates.ls)) {
        if(length(env_id.ls[[env.buffer$v[ij]]][[i]]) > 0) {
          for(k in vars$all) {
            env.buffer[ij,k] <- mean(env.df[env_dates.ls[[j]],][env_id.ls[[env.buffer$v[ij]]][[i]],][[k]])
          }
        }
        ij <- ij+1
      }
      if(i %% 50 == 0) {
        cat(i, "of", nrow(site.buffer), "--",
            round(as.numeric(difftime(Sys.time(), startTime, units="min")), 2),
            "minutes \n")
      }
    }
  } else {
    for(i in 1:nrow(site.buffer)) {
      for(j in 1:length(env_dates.ls)) {
        if(length(env_id.ls[[env.buffer$v[ij]]][[i]]) > 0) {
          for(k in vars$sea) {
            env.buffer[ij,k] <- mean(env.df[env_dates.ls[[j]],][env_id.ls[[env.buffer$v[ij]]][[i]],][[k]])
          }
          sst_ij <- env.df[env_dates.ls[[j]],][env_id.ls[[env.buffer$v[ij]]][[i]],][["sst"]]
          elev_ij <- env.df[env_dates.ls[[j]],][env_id.ls[[env.buffer$v[ij]]][[i]],][["elev"]]
          env.buffer[ij,"sst"] <- mean(sst_ij[elev_ij==0], na.rm=T)
        }
        ij <- ij+1
      }
      if(i %% 50 == 0) {
        cat(i, "of", nrow(site.buffer), "--",
            round(as.numeric(difftime(Sys.time(), startTime, units="min")), 2),
            "minutes \n")
      }
    }
  }

  env.buffer <- env.buffer |>
    group_by(siteid, quadrant) |>
    mutate(across(any_of(vars$all),
                  ~rollmean(.x, k=7, na.pad=T, align="right"),
                  .names="{.col}AvgWk")) |>
    mutate(across(any_of(paste0(vars$all, "AvgWk")),
                  ~.x - lag(.x),
                  .names="{.col}Delta")) |>
    # ungroup() |>
    # mutate(yday=yday(date)) |>
    # group_by(siteid, quadrant) |>
    # mutate(across(any_of(vars$all),
    #               ~detrend_loess(yday, .x, span=0.3),
    #               .names="{.col}AvgDt")) |>
    ungroup() |>
    select(siteid, quadrant, date,
           any_of(paste0(vars$all, "AvgWk")),
           any_of(paste0(vars$all, "AvgWkDelta")),
           any_of(paste0(vars$all, "AvgDt"))) |>
    group_by(siteid, date) |>
    mutate(across(where(is.numeric), na.aggregate)) |>
    ungroup()
  return(env.buffer)
}
