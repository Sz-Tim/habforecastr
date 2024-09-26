# Read and clean sites

#' Read and clean sites
#'
#' @param url_sites URL of site table in HABReports database
#' @param dateStart Start date; ignores sites only operating before this
#'
#' @return data.frame with sin, east, north, and date
#' @export
read_and_clean_sites <- function(url_sites, dateStart) {
  paste0(url_sites, "?fromdate=gte.", dateStart) |>
    url() |>
    readLines(warn=F) |>
    fromJSON() |> as_tibble() |>
    filter(east < 7e5,
           north < 125e4,
           !(east==0 & north==0),
           sin != "-99",
           !is.na(fromdate) & !is.na(todate),
           fromdate != todate) |>
    mutate(fromdate=date(fromdate), todate=date(todate)) |>
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



#' Title
#'
#' @param url_fsa
#' @param hab_i
#' @param sites
#' @param dateStart
#'
#' @return
#' @export
read_and_clean_fsa <- function(url_fsa, hab_i, sites, dateStart="2016-01-01") {
  paste0(url_fsa, "?date_collected=gte.", dateStart) |>
    url() |>
    readLines(warn=F) |>
    fromJSON() |> as_tibble() |>
    filter(!is.na(date_collected)) |>
    mutate(datetime_collected=as_datetime(date_collected),
           date=date(datetime_collected)) |>
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


#' Title
#'
#' @param url_cefas
#' @param tox_i
#' @param sites
#' @param dateStart
#'
#' @return
#' @export
read_and_clean_cefas <- function(url_cefas, tox_i, sites, dateStart="2016-01-01") {
  paste0(url_cefas, "?date_collected=gte.", dateStart) |>
    url() |>
    readLines(warn=F) |>
    fromJSON() |> as_tibble() |>
    filter(!is.na(date_collected) & sin != "-99") |>
    mutate(datetime_collected=as_datetime(date_collected),
           date=date(datetime_collected)) |>
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


#' Title
#'
#' @param url_mowi
#' @param url_ssf
#' @param fish_i
#' @param sites
#' @param dateStart
#'
#' @return
#' @export
read_and_clean_fish <- function(url_mowi, url_ssf, fish_i, sites, dateStart="2016-01-01") {
  bind_rows(url(paste0(url_mowi, "?date_collected=gte.", dateStart)) |>
              readLines(warn=F) |>
              fromJSON() |> as_tibble(),
            url(paste0(url_ssf, "?date_collected=gte.", dateStart)) |>
              readLines(warn=F) |>
              fromJSON() |> as_tibble()) |>
    filter(!is.na(date_collected)) |>
    mutate(datetime_collected=as_datetime(date_collected),
           date=date(datetime_collected)) |>
    mutate(across(any_of(fish_i$full), ~na_if(.x, -99))) |>
    group_by(sin) |> mutate(N=n()) |> ungroup() |> filter(N > 2) |>
    select(oid, sin, date, easting, northing, any_of(fish_i$full)) |>
    left_join(sites, by=c("sin", "date")) |>
    mutate(east=if_else(is.na(east), easting, east),
           north=if_else(is.na(north), northing, north)) |>
    rename(obsid=oid) |>
    group_by(sin) |> mutate(lon=median(east), lat=median(north)) |> ungroup() |>
    rename(any_of(setNames(fish_i$full, fish_i$abbr))) |>
    select(obsid, lon, lat, sin, date, any_of(fish_i$abbr)) |>
    arrange(sin, date) |>
    filter(sin!=0) |>
    group_by(date, sin) |>
    summarise(across(where(is.numeric), ~mean(.x, na.rm=T)))
}
