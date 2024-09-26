#' Shortcut for dir(..., full.names=T)
#'
#' @param ...
#'
#' @return
#' @export
dirf <- function(...) {
  dir(..., full.names=T)
}



#' Shortcut for dir(..., full.names=T, recursive=T)
#'
#' @param ...
#'
#' @return
#' @export
dirrf <- function(...) {
  dir(..., full.names=T, recursive=T)
}



#' Shortcut for dir(..., recursive=T)
#'
#' @param ...
#'
#' @return
#' @export
dirr <- function(...) {
  dir(..., recursive=T)
}
