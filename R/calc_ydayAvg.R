
#' Calculate average for each day of year
#'
#' @param data_df Dataframe
#' @param out_f Output file
#' @param ... Columns to group by
#'
#' @return Written object to out_f only
#' @export
calc_ydayAvg <- function(data_df, out_f, ...) {
  data_df |>
    mutate(yday=yday(date)) |>
    group_by(...) |>
    summarise(across(where(is.numeric), ~mean(.x, na.rm=T))) |>
    ungroup() |>
    saveRDS(out_f)
}