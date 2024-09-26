#' Download CMEMS layers
#'
#' @param userid
#' @param pw
#' @param i.df
#' @param bbox
#' @param nDays_buffer
#' @param dateRng
#' @param out.dir
#'
#' @return
#' @export
get_CMEMS <- function(userid, pw, i.df, bbox, nDays_buffer, dateRng, out.dir,
                      toolbox=TRUE) {
  if(toolbox) {
    save(list=ls(all.names=TRUE), file="temp/get_CMEMS.RData")
    system2("bash", paste0(getwd(), "/code/00_getCMEMS.sh"))
    file.remove("temp/get_CMEMS.RData")
    return("Finished running /code/getCMEMS.sh")
  }

  library(tidyverse); library(ncdf4); library(lubridate); library(glue)

  for(i in 1:nrow(i.df)) {
    # load .nc file
    url <- paste0("https://", userid, ":", pw, "@",
                  i.df$server[i], "/thredds/dodsC/", i.df$ID[i])
    nc <- nc_open(url)

    # get metadata and identify lon/lat/time indexes to extract
    nc.lon <- ncvar_get(nc, "longitude")
    nc.lat <- ncvar_get(nc, "latitude")
    nc.time <- ncvar_get(nc, "time")
    time_units <- ncatt_get(nc, "time", "units")
    nc.origin <- as_datetime(str_split_fixed(time_units$value, " ", 3)[3])
    nc.dt <- switch(str_split_fixed(time_units$value, " ", 3)[1],
                    seconds=dseconds,
                    minutes=dminutes,
                    hours=dhours)
    nc.date <- date(nc.origin + nc.dt(nc.time))
    lon_i <- which(between(nc.lon, bbox$xmin, bbox$xmax))
    lat_i <- which(between(nc.lat, bbox$ymin, bbox$ymax))
    time_i <- which(between(nc.date, dateRng[1]-nDays_buffer, dateRng[2]+nDays_buffer))
    var.start <- c(lon_i[1], lat_i[1], 1, time_i[1])
    var.count <- c(length(lon_i), length(lat_i), 1, length(time_i))

    # extract variable
    if(length(names(nc$var)) > 1) { cat("Too many vars in nc:", names(nc$var))}
    nc.var <- ncvar_get(nc, names(nc$var), start=var.start, count=var.count)
    nc_close(nc)

    # reshape and save
    expand_grid(date=nc.date[time_i],
                lat=c(nc.lat[lat_i]),
                lon=c(nc.lon[lon_i])) |>
      mutate(source=i.df$source[i],
             vals=c(nc.var)) |>
      rename_with(~gsub("vals", i.df$var[i], .x)) |>
      saveRDS(glue("{out.dir}/cmems_{i.df$var[i]}_{i.df$source[i]}.rds"))
  }
}
