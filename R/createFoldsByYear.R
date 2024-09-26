#' Identify indexes for ML training
#'
#' @param data.df
#'
#' @return
#' @export
createFoldsByYear <- function(data.df) {
  folds_out <- data.df |> mutate(rowNum=row_number()) |>
    group_by(year) |> group_split() |> map(~.x$rowNum)
  folds_out <- folds_out[-1]
  folds_in <- map(folds_out, ~(1:nrow(data.df))[-.x])
  return(list(i.in=folds_in, i.out=folds_out))
}


