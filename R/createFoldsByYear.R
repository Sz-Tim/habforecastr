#' Identify indexes for ML training
#'
#' @param data.df Dataframe
#'
#' @return List of 2 lists where `i.in` gives the row indexes to include in
#' each fold, and `i.out` gives the indexes to exclude in each fold.
#' @export
createFoldsByYear <- function(data.df) {
  folds_out <- data.df |> mutate(rowNum=row_number()) |>
    group_by(year) |> group_split() |> map(~.x$rowNum)
  folds_out <- folds_out[-1]
  folds_in <- map(folds_out, ~(1:nrow(data.df))[-.x])
  return(list(i.in=folds_in, i.out=folds_out))
}


