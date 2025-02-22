% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_site_geography.R
\name{get_shortestPaths}
\alias{get_shortestPaths}
\title{Find pairwise shortest in-ocean paths between sites}
\usage{
get_shortestPaths(ocean.path, site.df, site_savePath = NULL)
}
\arguments{
\item{ocean.path}{A character string specifying the path to the ocean raster file.}

\item{site.df}{A data frame containing site information, including longitude and latitude.}

\item{site_savePath}{A character string specifying the path to save the updated site data frame. Default is NULL.}

\item{transMx.path}{A character string specifying the path to save the transition matrix. Default is NULL.}

\item{recalc_transMx}{A logical value indicating whether to recalculate the transition matrix. Default is FALSE.}
}
\value{
A list containing the updated site data frame and a data frame of pairwise shortest paths.
}
\description{
This function calculates the pairwise shortest in-ocean paths between sites using a specified ocean raster.
}
\examples{
\dontrun{
ocean.path <- "path/to/ocean_raster.tif"
site.df <- data.frame(siteid = 1:5, lon = runif(5, -10, 10), lat = runif(5, -10, 10))
result <- get_shortestPaths(ocean.path, site.df)
}
}
