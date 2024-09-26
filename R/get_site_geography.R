#' Find pairwise shortest in-ocean paths between sites
#'
#' @param ocean.path
#' @param site.df
#' @param transMx.path
#' @param recalc_transMx
#'
#' @return
#' @export
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
#' @param locs
#' @param ras
#' @param layer
#'
#' @return
#' @export
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



#' Shift out-of-bounds pointsto nearest cell
#'
#' @param locs
#' @param ras
#' @param layer
#' @param move
#' @param distance
#' @param showchanges
#'
#' @return
#' @export
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
#' @param site.df
#' @param fetch.path
#'
#' @return
#' @export
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
#' @param site.df
#' @param coast.path
#' @param buffer
#' @param nDir
#'
#' @return
#' @export
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
