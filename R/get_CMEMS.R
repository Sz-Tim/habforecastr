#' Download CMEMS layers
#'
#' TODO: This is deprecated following the updated access through the `copernicusmarine` toolbox.
#' This function downloads Copernicus Marine Environment Monitoring Service (CMEMS) layers based on specified parameters.
#'
#' @param userid A character string specifying the user ID for CMEMS access.
#' @param pw A character string specifying the password for CMEMS access.
#' @param i.df A data frame containing information about the CMEMS layers to download, including server and variable IDs.
#' @param bbox A bounding box specifying the geographical area of interest.
#' @param nDays_buffer An integer specifying the number of days to buffer around the date range.
#' @param dateRng A vector of two dates specifying the date range for the data.
#' @param out.dir A character string specifying the output directory for saving the downloaded data.
#' @param toolbox A logical value indicating whether to use a toolbox script for downloading. Default is TRUE.
#'
#' @return A message indicating the completion of the download process.
#' @export
#'
#' @examples
#' \dontrun{
#' userid <- "your_userid"
#' pw <- "your_password"
#' i.df <- data.frame(server = "my_server", ID = "my_id", source = "my_source", var = "my_var")
#' bbox <- list(xmin = -10, xmax = 10, ymin = -10, ymax = 10)
#' nDays_buffer <- 5
#' dateRng <- as.Date(c("2023-01-01", "2023-12-31"))
#' out.dir <- "path/to/output"
#' get_CMEMS(userid, pw, i.df, bbox, nDays_buffer, dateRng, out.dir)
#' }
get_CMEMS <- function(userid, pw, i.df, bbox, nDays_buffer, dateRng, out.dir,
                      toolbox=TRUE) {
  if(toolbox) {
    file.remove("temp/")
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
    nc.date <- lubridate::date(nc.origin + nc.dt(nc.time))
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
