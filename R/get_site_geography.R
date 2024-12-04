#' Find pairwise shortest in-ocean paths between sites
#'
#' This function calculates the pairwise shortest in-ocean paths between sites using a specified ocean raster.
#'
#' @param ocean.path A character string specifying the path to the ocean raster file.
#' @param site.df A data frame containing site information, including longitude and latitude.
#' @param transMx.path A character string specifying the path to save the transition matrix. Default is NULL.
#' @param recalc_transMx A logical value indicating whether to recalculate the transition matrix. Default is FALSE.
#' @param site_savePath A character string specifying the path to save the updated site data frame. Default is NULL.
#'
#' @return A list containing the updated site data frame and a data frame of pairwise shortest paths.
#' @export
#'
#' @examples
#' \dontrun{
#' ocean.path <- "path/to/ocean_raster.tif"
#' site.df <- data.frame(siteid = 1:5, lon = runif(5, -10, 10), lat = runif(5, -10, 10))
#' result <- get_shortestPaths(ocean.path, site.df)
#' }
get_shortestPaths <- function(ocean.path, site.df, site_savePath=NULL) {
  library(raster); library(gdistance);
  library(terra); library(spaths)
  library(tidyverse); library(sf); library(glue);

  # adapted from:
  # https://agrdatasci.github.io/gdistance/reference/index.html

  # load ocean raster and calculate transition matrix
  mesh.r <- raster(ocean.path)
  crs(mesh.r) <- CRS("+init=epsg:27700")

  # locate sites in mesh
  set_ll_warn(TRUE)
  site.spdf <- SpatialPointsDataFrame(site.df[,c("lon", "lat")],
                                      data=site.df[,"siteid"],
                                      proj4string=CRS("+init=epsg:27700")) |>
    points2nearestcell(mesh.r) |>
    as.data.frame()
  site.df_new <- site.df |> dplyr::select(siteid, sin) |> left_join(site.spdf)
  if(!is.null(site_savePath)) {
    saveRDS(site.df_new, site_savePath)
  }

  # Pairwise each i to all others within site.df
  mesh.r <- rast(ocean.path)
  crs(mesh.r) <- crs("+init=epsg:27700")
  site.sf <- site.df_new |>
    st_as_sf(coords=c("lon", "lat"), crs=27700) |>
    st_crop(mesh.r)
  out_paths <- spaths::shortest_paths(mesh.r,
                                      site.sf,
                                      output="distances")
  # dist.df <- map_dfr(1:nrow(site.spdf),
  #                    ~shortestPath(mesh.tmx,
  #                                  as.matrix(site.spdf[.x, c("lon", "lat")]),
  #                                  as.matrix(site.spdf[-.x, c("lon", "lat")]),
  #                                  output="SpatialLines") |>
  #                      st_as_sf() |>
  #                      mutate(origin=site.spdf$siteid[.x],
  #                             destination=site.spdf$siteid[-.x],
  #                             distance=st_length(.)) |>
  #                      st_drop_geometry())

  return(list(site.df=site.df_new, dist.df=out_paths))
}





#' Identify whether point locations fall within valid raster cells
#'
#' This function checks whether point locations fall within valid (non-NA) cells of a raster layer.
#'
#' @param locs A spatial object (e.g., `SpatialPointsDataFrame`) containing point locations.
#' @param ras A raster object.
#' @param layer An integer specifying the layer of the raster to check. Default is 1.
#'
#' @return A logical vector indicating whether each point location falls within a valid raster cell (FALSE) or not (TRUE).
#' @export
#'
#' @examples
#' \dontrun{
#' library(raster)
#' locs <- SpatialPointsDataFrame(coords = matrix(runif(10), ncol = 2), data = data.frame(id = 1:5))
#' ras <- raster(matrix(runif(100), 10, 10))
#' result <- point_in_cell(locs, ras)
#' }
point_in_cell <- function (locs, ras, layer=1) {
  # copied from rSDM since not available for newer R versions
  if (!isTRUE(raster::compareCRS(locs, ras))) {
    stop("Coordinate data and raster object must have the same projection. Check their CRS or proj4string")
  }
  else {
    if (nlayers(ras) > 1)
      ras <- raster(ras, layer)
    rasvals <- raster::extract(ras, locs)
    missing <- is.na(rasvals)
    missing
  }
}



#' Shift out-of-bounds points to nearest cell
#'
#' This function shifts point locations that fall outside valid raster cells to the nearest valid cell.
#'
#' @param locs A spatial object (e.g., `SpatialPointsDataFrame`) containing point locations.
#' @param ras A raster object.
#' @param layer An integer specifying the layer of the raster to check. Default is 1.
#' @param move A logical value indicating whether to move the points to the nearest valid cell. Default is TRUE.
#' @param distance A numeric value specifying the maximum distance to move points. Points farther than this distance will not be moved. Default is NULL.
#' @param showchanges A logical value indicating whether to print the changes made. Default is TRUE.
#'
#' @return The modified spatial object with points shifted to the nearest valid raster cell.
#' @export
#'
#' @examples
#' \dontrun{
#' library(raster)
#' locs <- SpatialPointsDataFrame(coords = matrix(runif(10), ncol = 2), data = data.frame(id = 1:5))
#' ras <- raster(matrix(runif(100), 10, 10))
#' result <- points2nearestcell(locs, ras)
#' }
points2nearestcell <- function (locs=NULL, ras=NULL, layer=1,
                                move=T, distance=NULL, showchanges=T) {
  # copied from rSDM since not available for newer R versions
  miss <- point_in_cell(locs, ras, layer)
  if (sum(miss) > 0) {
    coord.miss <- sp::coordinates(locs[miss, ])
    if (nlayers(ras) > 1)
      ras <- raster::raster(ras, layer)
    cells.notNA <- raster::rasterToPoints(ras, spatial = TRUE)
    coord.ras <- sp::coordinates(cells.notNA)
    cell.id <- factor(seq_len(nrow(coord.ras)))
    nearest.cell <- class::knn1(coord.ras, coord.miss, cell.id)
    new.coords <- matrix(coord.ras[nearest.cell, ], ncol = 2)
    colnames(new.coords) <- c("longitude_new", "latitude_new")
    if (!is.null(distance)) {
      distances <- raster::pointDistance(coord.miss, new.coords,
                                         lonlat = raster::isLonLat(locs))
      x <- ifelse(distances < distance, new.coords[, 1],
                  coord.miss[, 1])
      y <- ifelse(distances < distance, new.coords[, 2],
                  coord.miss[, 2])
      new.coords <- cbind(longitude_new = x, latitude_new = y)
    }
    if (isTRUE(move)) {
      locs@coords[miss, ] <- new.coords
    }
    if (isTRUE(showchanges)) {
      coords <- data.frame(coord.miss, new.coords)
      distances <- round(raster::pointDistance(coord.miss,
                                               new.coords, lonlat = raster::isLonLat(locs)))
      moved <- apply(coords, 1, function(x) {
        !isTRUE(identical(x[1], x[3]) & identical(x[2],
                                                  x[4]))
      })
      coords <- cbind(coords, distances, moved)
      print(coords)
      message(sum(moved), " out of ", nrow(coords),
              " points have been moved.")
    }
  }
  else message("All points fall within a raster cell")
  return(locs)
}





#' Extract fetch for each site point location
#'
#' This function extracts the fetch (distance over water) for each site point location from a specified raster file.
#'
#' @param site.df A data frame containing site information, including longitude and latitude.
#' @param fetch.path A character string specifying the path to the fetch raster file.
#'
#' @return A data frame with the original site information and an additional column for fetch.
#' @export
#'
#' @examples
#' \dontrun{
#' library(tidyverse)
#' library(sf)
#' site.df <- data.frame(siteid = 1:5, lon = runif(5, -10, 10), lat = runif(5, -10, 10))
#' fetch.path <- "path/to/fetch_raster.tif"
#' result <- get_fetch(site.df, fetch.path)
#' }
get_fetch <- function(site.df, fetch.path) {
  library(tidyverse); library(sf)
  site.df |>
    st_as_sf(coords=c("lon", "lat"), crs=27700, remove=F) %>%
    mutate(fetch=raster::extract(raster(fetch.path), .,
                                 small=T, buffer=1e2, fun=mean)) %>%
    mutate(fetch=if_else(is.na(fetch),
                         raster::extract(raster(fetch.path), .,
                                         small=T, buffer=5e2, fun=mean),
                         fetch)) %>%
    mutate(fetch=if_else(is.na(fetch),
                         raster::extract(raster(fetch.path), .,
                                         small=T, buffer=1e3, fun=mean),
                         fetch)) %>%
    mutate(fetch=if_else(is.na(fetch),
                         raster::extract(raster(fetch.path), .,
                                         small=T, buffer=5e3, fun=mean),
                         fetch)) %>%
    mutate(fetch=if_else(is.na(fetch),
                         raster::extract(raster(fetch.path), .,
                                         small=T, buffer=10e3, fun=mean),
                         fetch)) %>%
    mutate(fetch=if_else(is.na(fetch),
                         raster::extract(raster(fetch.path), .,
                                         small=T, buffer=50e3, fun=mean),
                         fetch)) %>%
    st_drop_geometry()
}





#' Identify most-open bearing for each site point location
#'
#' This function identifies the most-open bearing (direction with the least obstruction) for each site point location based on a specified coastal path.
#'
#' @param site.df A data frame containing site information, including longitude and latitude.
#' @param coast.path A character string specifying the path to the coastal shapefile.
#' @param buffer A numeric value specifying the buffer distance around each site point. Default is 10,000 meters.
#' @param nDir An integer specifying the number of directions to consider. Default is 120.
#'
#' @return A data frame with the original site information and an additional column for the most-open bearing.
#' @export
#'
#' @examples
#' \dontrun{
#' library(tidyverse)
#' library(sf)
#' site.df <- data.frame(siteid = 1:5, lon = runif(5, -10, 10), lat = runif(5, -10, 10))
#' coast.path <- "path/to/coast_shapefile.shp"
#' result <- get_openBearing(site.df, coast.path)
#' }
get_openBearing <- function(site.df, coast.path, buffer=10e3, nDir=120) {
  library(tidyverse); library(sf)

  hub.df <- site.df |>
    select(siteid, lon, lat) |>
    st_as_sf(coords=c("lon", "lat"), crs=27700)
  spoke.df <- hub.df |>
    st_buffer(dist=buffer, nQuadSegs=nDir/4) |>
    st_cast("POINT") |>
    group_by(siteid) |>
    mutate(spoke.id=row_number()) |>
    filter(spoke.id <= nDir)
  hubRep.df <- full_join(hub.df, spoke.df |> st_drop_geometry())
  coords <- cbind(st_coordinates(hubRep.df), st_coordinates(spoke.df))
  spoke.lines <- lapply(1:nrow(coords),
                        function(i){
                          st_linestring(matrix(coords[i,], ncol=2, byrow=TRUE))
                        }) |>
    st_sfc() |> st_as_sf() |> st_set_crs(27700) |>
    rename(geometry=x) |>
    mutate(siteid=spoke.df$siteid,
           spoke.id=spoke.df$spoke.id)
  spoke.mesh <- st_intersection(spoke.lines, st_read(coast.path)) |>
    st_cast("LINESTRING") %>%
    mutate(len=round(as.numeric(st_length(.)))) |>
    group_by(siteid) |>
    filter(len==max(len)) |>
    ungroup() |>
    st_cast("POINT")
  bearings <- as.numeric(lwgeom::st_geod_azimuth(st_transform(spoke.mesh, 4326)))
  bearing.df <- spoke.mesh[1:nrow(spoke.mesh) %% 2 == 0,] |>
    st_drop_geometry() |>
    mutate(bearing=bearings[1:length(bearings) %% 2 == 1]) |>
    group_by(siteid) |>
    summarise(openBearing=median(bearing)) |>
    ungroup()
  return(full_join(site.df, bearing.df, by="siteid"))
}
