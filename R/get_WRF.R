#' Download WRF data
#'
#' This function downloads Weather Research and Forecasting (WRF) model data based on specified parameters.
#'
#' @param wrf.dir A character string specifying the directory or URL where WRF files are located.
#' @param nDays_buffer An integer specifying the number of days to buffer around the date range.
#' @param dateRng A vector of two dates specifying the date range for the data.
#' @param out.dir A character string specifying the output directory for saving the downloaded data.
#' @param forecast A logical value indicating whether to download forecast data. Default is FALSE.
#'
#' @return A message indicating the completion of the download process.
#' @export
#'
#' @examples
#' \dontrun{
#' library(tidyverse)
#' library(ncdf4)
#' library(lubridate)
#' library(glue)
#' wrf.dir <- "https://thredds.sams.ac.uk/thredds/catalog/scoats-wrf/Archive/catalog.html"
#' nDays_buffer <- 5
#' dateRng <- as.Date(c("2023-01-01", "2023-01-14"))
#' out.dir <- "path/to/output"
#' get_WRF(wrf.dir, nDays_buffer, dateRng, out.dir)
#' }
get_WRF <- function(wrf.dir, nDays_buffer, dateRng, out.dir, forecast=F) {
  library(tidyverse); library(ncdf4); library(lubridate); library(glue);
  library(xml2); library(rvest)
  dir.create(out.dir, showWarnings=F)

  # metadata for all WRF files within timespan
  if(grepl("https", wrf.dir)) {
    thredd_base <- "https://thredds.sams.ac.uk/thredds/"
    Archive <- ifelse(forecast, "Archive_forecast", "Archive")
    catalog <- ifelse(forecast, "F/catalog.html", "/catalog.html")

    thredds_head <- glue("{thredd_base}catalog/scoats-wrf/{Archive}/netcdf_")
    wrf_links <- map(seq(year(dateRng[1]), year(dateRng[2]), by=1),
                     ~glue("{thredds_head}{.x}{catalog}") |>
                       read_html() |> html_nodes("a") |> html_attr("href") |>
                       grep("netcdf_20.*/wrf_", x=_, value=T) |>
                       str_split_fixed(glue("{Archive}/"), 2) |>
                       magrittr::extract(,2)) |>
      do.call('c', args=_)
    wrf_i <- tibble(fname=wrf_links)
    wrf_base <- glue("{thredd_base}dodsC/scoats-wrf/{Archive}")
  } else {
    wrf_i <- tibble(fname=dir(wrf.dir, ".nc$", recursive=T))
    wrf_base <- wrf.dir
  }
  wrf_i <- wrf_i |>
    mutate(res=str_sub(fname, -6-forecast, -4-forecast),
           day_0=str_sub(fname, 23+forecast, 24+forecast),
           month_0=str_sub(fname, 21+forecast, 22+forecast),
           year_0=str_sub(fname, 17+forecast, 20+forecast),
           day_1=str_sub(fname, 28+forecast, 29+forecast),
           month_1=str_sub(fname, 26+forecast, 27+forecast),
           year_1=if_else(month_0=="12" & month_1=="01",
                          as.character(as.numeric(year_0) + 1),
                          year_0),
           date_0=ymd(paste0(year_0, month_0, day_0)),
           date_1=ymd(paste0(year_1, month_1, day_1))) |>
    filter(date_0 >= dateRng[1]-nDays_buffer,
           date_1 <= dateRng[2]+nDays_buffer)
  wrf_dates <- unique(wrf_i$date_0)

  for(i in 1:length(wrf_dates)) {
    wrf_i.i <- filter(wrf_i, date_0==wrf_dates[i])
    nc_f.i <- map(wrf_i.i$fname, ~glue("{wrf_base}/{.x}")) |>
      set_names(wrf_i.i$res)
    nc.ls <- map(nc_f.i, nc_open)
    time.ls <- map(nc.ls,
                   ~tibble(Times=ncvar_get(.x, "Times")) |>
                     mutate(Time.dt=as_datetime(str_replace(Times, "_", " "))))
    var.ls <- map2_dfr(nc.ls, time.ls,
                       ~expand_grid(time=.y$Time.dt,
                                    lat_i=1:(.x$dim$south_north$len),
                                    lon_i=1:(.x$dim$west_east$len)) |>
                         mutate(date=date(time),
                                U=c(ncvar_get(.x, "U10")),
                                V=c(ncvar_get(.x, "V10")),
                                Shortwave=c(ncvar_get(.x, "Shortwave")),
                                Precipitation=c(ncvar_get(.x, "Precipitation")),
                                sst=c(ncvar_get(.x, "sst"))) |>
                         group_by(time) |>
                         mutate(i=row_number()) |>
                         ungroup(),
                       .id="res") |>
      group_by(res, date) |>
      group_split()
    for(j in seq_along(var.ls)) {
      j.fname <- glue("{out.dir}/wrf_{var.ls[[j]]$date[1]}_{var.ls[[j]]$res[1]}.rds")
      j.fnameF <- glue("{out.dir}/wrfF_{var.ls[[j]]$date[1]}_{var.ls[[j]]$res[1]}.rds")
      if(!file.exists(j.fname)) {
        var.ls[[j]] |>
          mutate(sst=if_else(sst > -100, sst, NA_real_)) |>
          group_by(i, lat_i, lon_i, date) |>
          summarise(U_mn=mean(U),
                    V_mn=mean(V),
                    UV_mn=mean(sqrt(U^2 + V^2)),
                    Shortwave=sum(Shortwave),
                    Precip=sum(Precipitation),
                    sst=mean(sst, na.rm=T)) |>
          ungroup() |>
          rename(U=U_mn, V=V_mn, UV=UV_mn) |>
          saveRDS(ifelse(forecast, j.fnameF, j.fname))
      }
    }
    # save domain extents IF all three domains are represented
    if(length(nc_f.i)==3 & all(c("d01", "d02", "d03") %in% names(nc_f.i))) {
      nest_WRF_domains(nc.ls) |>
        walk(~saveRDS(.x, glue("{out.dir}/wrfDomains_{wrf_dates[i]}_{.x$res[1]}.rds")))
    }
    walk(nc.ls, nc_close)
  }
}







#' Identify preserved points within nested WRF domains
#'
#' This function identifies preserved points within nested Weather Research and Forecasting (WRF) model domains.
#'
#' @param nc.ls A list of NetCDF file connections for the WRF domains.
#'
#' @return A list of data frames containing the preserved points for each WRF domain.
#' @export
#'
#' @examples
#' \dontrun{
#' library(ncdf4)
#' nc.ls <- list(d01 = nc_open("path/to/d01.nc"), d02 = nc_open("path/to/d02.nc"), d03 = nc_open("path/to/d03.nc"))
#' result <- nest_WRF_domains(nc.ls)
#' }
nest_WRF_domains <- function(nc.ls) {
  library(tidyverse); library(ncdf4); library(sf)
  nDomain <- length(nc.ls)
  lon <- map(nc.ls, ~ncvar_get(.x, "XLONG"))
  lat <- map(nc.ls, ~ncvar_get(.x, "XLAT"))
  elev <- map(nc.ls, ~ncvar_get(.x, "HGT"))
  coord.rows <- map(lon, ~matrix(1:nrow(.x), nrow=nrow(.x), ncol=ncol(.x)))
  coord.cols <- map(lon, ~matrix(1:ncol(.x), nrow=nrow(.x), ncol=ncol(.x), byrow=T))
  coord.wrf <- map(1:nDomain,
                   ~tibble(lon=c(lon[[.x]]),
                           lat=c(lat[[.x]]),
                           elev=c(elev[[.x]]),
                           row=c(coord.rows[[.x]]),
                           col=c(coord.cols[[.x]])) |>
                     mutate(i=row_number(),
                            res=names(nc.ls)[.x]) |>
                     st_as_sf(coords=c("lon", "lat"), crs=4326))
  hull.wrf <- map(1:nDomain,
                  ~sf::st_convex_hull(st_union(coord.wrf[[.x]])) |>
                    st_as_sf() |> mutate(res=names(nc.ls)[.x]))
  merge.wrf <- map(coord.wrf, ~.x |> mutate(in_d02=0, in_d03=0))
  if(nDomain==3) {
    merge.wrf <- map(merge.wrf,
                     ~.x %>%
                       mutate(in_d02=as.numeric(st_within(., hull.wrf[[2]], sparse=F)[,1]),
                              in_d03=as.numeric(st_within(., hull.wrf[[3]], sparse=F)[,1])))
  }
  if(nDomain==2) {
    merge.wrf <- map(merge.wrf,
                     ~.x %>%
                       mutate(in_d02=as.numeric(st_within(., hull.wrf[[2]], sparse=F)[,1])))
  }
  merge.wrf <- map(merge.wrf,
                   ~.x %>%
                     mutate(lon=st_coordinates(.)[,1],
                            lat=st_coordinates(.)[,2]) |>
                     st_drop_geometry() |>
                     filter(res=="d03" | (res=="d02" & !in_d03) | (res=="d01" & !in_d02)) |>
                     select(-in_d02, -in_d03))
  return(merge.wrf)
}





#' Filter WRF data to include only preserved points in each domain
#'
#' This function filters Weather Research and Forecasting (WRF) model data to include only preserved points in each domain.
#'
#' @param domain A character string specifying the WRF domain (e.g., "d01", "d02", "d03").
#' @param wrf.out A character string specifying the output directory where WRF data is stored.
#' @param v2_start A character string specifying the start date for version 2 of the data. Default is NULL.
#' @param refreshStart A character string specifying the start date for refreshing the data. Default is NULL.
#'
#' @return A data frame containing the filtered WRF data with preserved points.
#' @export
#'
#' @examples
#' \dontrun{
#' library(tidyverse)
#' domain <- "d01"
#' wrf.out <- "path/to/wrf_output"
#' result <- subset_WRF(domain, wrf.out)
#' }
subset_WRF <- function(domain, wrf.out, v2_start=NULL, refreshStart=NULL) {
  f.domain <- dirf(wrf.out, glue("wrfDomains_.*{domain}.rds"))
  f.wrf <- dirf(wrf.out, glue("wrfF?_.*{domain}.rds"))
  if(is.null(v2_start)) {
    domain.ls <- map(f.domain, readRDS)
    i_chg <- c(1, which(map_lgl(1:(length(domain.ls)-1),
                                ~!(identical(domain.ls[[.x]]$lon, domain.ls[[.x+1]]$lon) &
                                     identical(domain.ls[[.x]]$lat, domain.ls[[.x+1]]$lat))))+1)
    domain.ls <- domain.ls[i_chg]
  } else {
    i_chg <- c(1, grep(v2_start, f.domain))
    domain.ls <- map(f.domain[i_chg], readRDS)
  }
  v_dateRng <- tibble(v=seq_along(i_chg),
                      start=str_split_fixed(str_split_fixed(f.domain[i_chg], "Domains_", 2)[,2], "_d0", 2)[,1]) |>
    mutate(start=ymd(start),
           end=lead(start, default=ymd("3000-01-01")))
  v_dateRng$start[1] <- "2013-01-01"

  iwalk(domain.ls,
        ~.x |> mutate(version=.y) |>
          select(-row, -col) |>
          saveRDS(glue("{wrf.out}/domain_{domain}_{.y}.rds")))

  wrf.ls <- vector("list", length(f.wrf))
  for(i in seq_along(f.wrf)) {
    date_i <- str_split_fixed(str_split_fixed(f.wrf[i], "/wrfF?_", 2)[,2], "_d0", 2)[,1]
    if(date_i == "NA") next
    if(is.null(refreshStart) || date_i >= refreshStart) {
      v_i <- which(date_i >= v_dateRng$start & date_i < v_dateRng$end)
      wrf.ls[[i]] <- readRDS(f.wrf[i]) |>
        select(-lon_i, -lat_i) |>
        right_join(domain.ls[[v_i]] |> select(-row, -col), by="i") |>
        mutate(version=ifelse(is.null(v2_start), v_i, 1 + (date_i >= v2_start)))
    }
  }
  wrf.ls <- do.call('rbind', wrf.ls)
  return(wrf.ls)
}





#' Aggregate WRF Data
#'
#' This function aggregates Weather Research and Forecasting (WRF) model data from multiple domains, filters out invalid data, and applies transformations to certain variables.
#'
#' @param wrf.out A character string specifying the output directory where WRF data is stored.
#' @param v2_start A Date object specifying the start date for version 2 of the data. Default is `ymd("2019-04-01")`.
#' @param refreshStart A Date object specifying the start date for refreshing the data. Default is NULL.
#'
#' @return A data frame containing the aggregated WRF data with transformations applied.
#' @export
#'
#' @examples
#' \dontrun{
#' library(lubridate)
#' wrf.out <- "path/to/wrf_output"
#' result <- aggregate_WRF(wrf.out)
#' }
aggregate_WRF <- function(wrf.out, v2_start=ymd("2019-04-01"), refreshStart=NULL) {
  d01 <- subset_WRF("d01", wrf.out, v2_start=v2_start, refreshStart)
  d02 <- subset_WRF("d02", wrf.out, v2_start=v2_start, refreshStart)
  d03 <- subset_WRF("d03", wrf.out, v2_start=v2_start, refreshStart)
  wrf.df <-  bind_rows(d01, d02, d03) |>
    filter(!is.na(date)) |>
    arrange(date, res, i) |>
    group_by(date) |>
    mutate(wrf_id=row_number()) |>
    ungroup() |>
    mutate(across(where(is.numeric), ~if_else(.x > 1e30, NA, .x))) |>
    mutate(yday=yday(date)) |>
    group_by(wrf_id, version, yday) |>
    mutate(across(where(is.numeric), zoo::na.aggregate)) |>
    ungroup() |>
    mutate(Shortwave=log1p(Shortwave),
           Precip=log1p(pmax(Precip, 0)*3600*24*1000), # m/s to mm/day
           UV=log1p(UV))
  return(wrf.df)
}
