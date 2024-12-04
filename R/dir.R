#' Shortcut for dir(..., full.names=T)
#'
#' @param ... Arguments to `dir()`. Usually path and possibly pattern
#'
#' @return A character vector containing the names of the files in the specified
#' directories (empty if there were no files). If a path does not exist or is
#' not a directory or is unreadable it is skipped.
#' @export
dirf <- function(...) {
  dir(..., full.names=T)
}



#' Shortcut for dir(..., full.names=T, recursive=T)
#'
#' @param ... Arguments to `dir()`. Usually path and possibly pattern
#'
#' @return A character vector containing the names of the files in the specified
#' directories (empty if there were no files). If a path does not exist or is
#' not a directory or is unreadable it is skipped.
#' @export
dirrf <- function(...) {
  dir(..., full.names=T, recursive=T)
}



#' Shortcut for dir(..., recursive=T)
#'
#' @param ... Arguments to `dir()`. Usually path and possibly pattern
#'
#' @return A character vector containing the names of the files in the specified
#' directories (empty if there were no files). If a path does not exist or is
#' not a directory or is unreadable it is skipped.
#' @export
dirr <- function(...) {
  dir(..., recursive=T)
}
