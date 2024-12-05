
#' Get Traffic light threshold level information
#'
#' This function retrieves and processes threshold level information for a given set of target parameters.
#'
#' @param target_i A tibble containing the target parameters to be joined with the threshold level information.
#'
#' @return A tibble with the combined and processed threshold level information.
#' @export
#'
#' @examples
#' target_i <- tibble(hab_parameter = c("param1", "param2"), abbr = c("P1", "P2"))
#' get_tl_info(target_i)
get_tl_info <- function(target_i) {
  library(jsonlite)
  library(tidyverse)
  params <- url("https://www.habreports.org/dbdatastuff/hab_parameters") |>
    readLines(warn=F) |>
    fromJSON() |> as_tibble()
  tl <- url("https://www.habreports.org/dbdatastuff/hab_symbology") |>
    readLines(warn=F) |>
    fromJSON() |> as_tibble()
  tl_info <- inner_join(
    params |> select(id, hab_parameter, param_name, units, conversion),
    tl |> select(parameter_id, min_ge, max_lt, label, alert) |> rename(id=parameter_id)
  )
  inner_join(tl_info, target_i, by=join_by(hab_parameter)) |>
    filter(min_ge != -99) |>
    group_by(abbr) |>
    mutate(alert=case_when(is.na(alert)~0,
                           alert=="warn"~1,
                           alert=="alert"~2),
           A=paste0("A", as.numeric(alert>0)),
           tl=factor(alert, levels=0:2, labels=c("TL0", "TL2", "TL3"))) |>
    ungroup()

}
