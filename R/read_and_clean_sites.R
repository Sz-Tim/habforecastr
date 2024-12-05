# Read and clean sites

#'
#' This function reads and cleans site data from the HABReports database.
#'
#' @param url_sites A character string specifying the URL of the site table in the HABReports database.
#' @param dateStart A Date object specifying the start date; sites operating only before this date will be ignored.
#'
#' @return A data frame with columns: sin, east, north, and date.
#' @export
#'
#' @examples
#' # Example usage:
#' # url_sites <- "http://example.com/sites"
#' # dateStart <- as.Date("2023-01-01")
#' # read_and_clean_sites(url_sites, dateStart)
read_and_clean_sites <- function(url_sites, dateStart) {
  library(tidyverse)
  paste0(url_sites, "?fromdate=gte.", dateStart) |>
    url() |>
    readLines(warn=F) |>
    fromJSON() |> as_tibble() |>
    filter(east < 7e5 &
           north < 125e4 &
           !(east==0 & north==0) &
           sin != "-99" &
           !is.na(fromdate) & !is.na(todate) &
           fromdate != todate) |>
    mutate(fromdate=lubridate::date(fromdate), todate=lubridate::date(todate)) |>
    rowwise() |>
    mutate(date=list(seq(fromdate, todate, by=1))) |>
    ungroup() |>
    arrange(sin, fromdate) |>
    select(sin, east, north, date) |>
    unnest(date) |>
    arrange(sin, date) |>
    group_by(sin, date) |>
    slice_head(n=1) |>
    ungroup()
}



#' Read and clean FSA data
#'
#' This function reads and cleans FSA data from the specified URL.
#'
#' @param url_fsa A character string specifying the URL of the FSA data.
#' @param hab_i A data frame containing information about HAB indices.
#' @param sites A data frame containing site information.
#' @param dateStart A character string specifying the start date (default is "2016-01-01").
#'
#' @return A cleaned data frame with columns: obsid, lon, lat, sin, date, and HAB indices.
#' @export
read_and_clean_fsa <- function(url_fsa, hab_i, sites, dateStart="2016-01-01") {
  library(tidyverse)
  paste0(url_fsa, "?date_collected=gte.", dateStart) |>
    url() |>
    readLines(warn=F) |>
    fromJSON() |> as_tibble() |>
    filter(!is.na(date_collected)) |>
    mutate(datetime_collected=as_datetime(date_collected),
           date=lubridate::date(datetime_collected)) |>
    mutate(across(any_of(hab_i$full), ~na_if(.x, -99))) |>
    group_by(sin) |> mutate(N=n()) |> ungroup() |> filter(N > 2) |>
    select(oid, sin, date, easting, northing, all_of(hab_i$full)) |>
    left_join(sites, by=c("sin", "date")) |>
    mutate(east=if_else(is.na(east), easting, east),
           north=if_else(is.na(north), northing, north)) |>
    rename(obsid=oid) |>
    group_by(sin) |> mutate(lon=median(east), lat=median(north)) |> ungroup() |>
    rename(all_of(setNames(hab_i$full, hab_i$abbr))) |>
    select(obsid, lon, lat, sin, date, all_of(hab_i$abbr)) |>
    arrange(sin, date)
}


#' Read and clean CEFAS data
#'
#' This function reads and cleans CEFAS data from the specified URL.
#'
#' @param url_cefas A character string specifying the URL of the CEFAS data.
#' @param tox_i A data frame containing information about toxin indices.
#' @param sites A data frame containing site information.
#' @param dateStart A character string specifying the start date (default is "2016-01-01").
#'
#' @return A cleaned data frame with columns: obsid, lon, lat, sin, date, and toxin indices.
#' @export
read_and_clean_cefas <- function(url_cefas, tox_i, sites, dateStart="2016-01-01") {
  library(tidyverse)
  paste0(url_cefas, "?date_collected=gte.", dateStart) |>
    url() |>
    readLines(warn=F) |>
    fromJSON() |> as_tibble() |>
    filter(!is.na(date_collected) & sin != "-99") |>
    mutate(datetime_collected=as_datetime(date_collected),
           date=lubridate::date(datetime_collected)) |>
    mutate(across(any_of(tox_i$full), ~if_else(.x == -99, NA_real_, .x)),
           across(any_of(tox_i$full), ~if_else(.x < 0, 0, .x))) |>
    group_by(sin, date) |> slice_head(n=1) |> ungroup() |>
    group_by(sin) |> mutate(N=n()) |> ungroup() |> filter(N > 2) |>
    select(oid, sin, date, easting, northing, all_of(tox_i$full)) |>
    left_join(sites, by=c("sin", "date")) |>
    mutate(east=if_else(is.na(east), easting, east),
           north=if_else(is.na(north), northing, north)) |>
    rename(obsid=oid) |>
    group_by(sin) |> mutate(lon=median(east), lat=median(north)) |> ungroup() |>
    filter(lat > 500000) |>
    rename(all_of(setNames(tox_i$full, tox_i$abbr))) |>
    select(obsid, lon, lat, sin, date, all_of(tox_i$abbr)) |>
    arrange(sin, date)
}


#' Read and clean fish data
#'
#' This function reads and cleans fish data from the specified URLs.
#'
#' @param url_mowi A character string specifying the URL of the MOWI data.
#' @param url_ssf A character string specifying the URL of the SSF data.
#' @param fish_i A data frame containing information about fish indices.
#' @param sites A data frame containing site information.
#' @param dateStart A character string specifying the start date (default is "2016-01-01").
#'
#' @return A cleaned data frame with columns: obsid, lon, lat, sin, date, and fish indices.
#' @export
read_and_clean_fish <- function(url_mowi, url_ssf, fish_i, sites, dateStart="2016-01-01") {
  library(tidyverse)
  bind_rows(url(paste0(url_mowi, "?date_collected=gte.", dateStart)) |>
              readLines(warn=F) |>
              fromJSON() |> as_tibble(),
            url(paste0(url_ssf, "?date_collected=gte.", dateStart)) |>
              readLines(warn=F) |>
              fromJSON() |> as_tibble()) |>
    filter(!is.na(date_collected)) |>
    mutate(datetime_collected=as_datetime(date_collected),
           date=lubridate::date(datetime_collected)) |>
    mutate(across(any_of(fish_i$full), ~na_if(.x, -99))) |>
    group_by(sin) |> mutate(N=n()) |> ungroup() |> filter(N > 2) |>
    select(oid, sin, date, easting, northing, any_of(fish_i$full)) |>
    mutate(easting=if_else(easting==0 & northing==0, NA_real_, easting),
           northing=if_else(easting==0 & northing==0, NA_real_, northing)) |>
    left_join(sites, by=c("sin", "date")) |>
    mutate(east=if_else(is.na(east), easting, east),
           north=if_else(is.na(north), northing, north)) |>
    rename(obsid=oid) |>
    group_by(sin) |> mutate(lon=median(east, na.rm=T), lat=median(north, na.rm=T)) |> ungroup() |>
    filter(!is.na(lon)) |>
    rename(any_of(setNames(fish_i$full, fish_i$abbr))) |>
    select(obsid, lon, lat, sin, date, any_of(fish_i$abbr)) |>
    arrange(sin, date) |>
    filter(sin!=0) |>
    group_by(date, sin) |>
    summarise(across(where(is.numeric), ~mean(.x, na.rm=T)))
}
