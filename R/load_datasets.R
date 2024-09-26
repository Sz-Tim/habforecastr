#' Load datasets for compilation
#'
#' @param sub.dir
#' @param target
#'
#' @return
#' @export
load_datasets <- function(sub.dir, target, ydayAvg.dir) {
  d.ls <- list(
    site=readRDS(glue("data/site_{target}_df.rds")),
    obs=readRDS(glue("data/{sub.dir}/{target}_obs.rds")),
    cmems.pt=readRDS(glue("data/{sub.dir}/cmems_sitePt_{target}.rds")),
    cmems.buf=readRDS(glue("data/{sub.dir}/cmems_siteBufferNSEW_{target}.rds")),
    wrf.pt=readRDS(glue("data/{sub.dir}/wrf_sitePt_{target}.rds")),
    wrf.buf=readRDS(glue("data/{sub.dir}/wrf_siteBufferNSEW_{target}.rds")),
    fsa=readRDS(glue("data/{sub.dir}/fsa_df.rds")),
    cefas=readRDS(glue("data/{sub.dir}/cefas_df.rds"))
  ) |>
    map(~.x |> mutate(across(matches("date"), ~as_date(.x))))
  if(target=="tox") {
    d.ls$habAvg <- readRDS(glue("data/{sub.dir}/tox_habAvg.rds"))
  }
  yday_env <- list(
    cmems.pt=readRDS(glue("data/{ydayAvg.dir}/ydayAvg_cmems_sitePt_{target}.rds")) |>
      filter(version==max(version)) |> select(-version),
    cmems.buf=readRDS(glue("data/{ydayAvg.dir}/ydayAvg_cmems_siteBufferNSEW_{target}.rds")) |>
      pivot_wider(names_from="quadrant", values_from=-(1:3), names_sep="Dir"),
    wrf.pt=readRDS(glue("data/{ydayAvg.dir}/ydayAvg_wrf_sitePt_{target}.rds")) |>
      filter(version==max(version)) |> select(-version),
    wrf.buf=readRDS(glue("data/{ydayAvg.dir}/ydayAvg_wrf_siteBufferNSEW_{target}.rds")) |>
      pivot_wider(names_from="quadrant", values_from=-(1:3), names_sep="Dir")
  )
  d.ls$compiled <- d.ls$site |> select(-sin) |>
    right_join(d.ls$obs, by="siteid", multiple="all") |>
    left_join(d.ls$cmems.pt |> select(-ends_with("Dt")), by=c("cmems_id", "date")) |>
    left_join(d.ls$cmems.buf |> select(-ends_with("Dt")) |>
                pivot_wider(names_from="quadrant", values_from=-(1:3), names_sep="Dir"),
              by=c("siteid", "date")) |>
    mutate(wrf_id=if_else(date < "2019-04-01", wrf_id.1, wrf_id.2)) |>
    select(-wrf_id.1, -wrf_id.2, -version) |>
    left_join(d.ls$wrf.pt |> select(-version, -ends_with("Dt")), by=c("wrf_id", "date")) |>
    left_join(d.ls$wrf.buf |> select(-ends_with("Dt")) |>
                pivot_wider(names_from="quadrant", values_from=-(1:3), names_sep="Dir"),
              by=c("siteid", "date")) |>
    mutate(year=year(date),
           yday=yday(date))
  for(i in seq_along(yday_env)) {
    env_id_col <- switch(i,
                         "1"=c("cmems_id", "yday"),
                         "2"=c("siteid", "yday"),
                         "3"=c("wrf_id", "yday"),
                         "4"=c("siteid", "yday"))
    varNames_i <- names(yday_env[[i]] |> select(-all_of(env_id_col)))
    for(j in varNames_i) {
      na_rows <- which(is.na(d.ls$compiled[[j]]))
      if(length(na_rows > 0)) {
        d_meta <- d.ls$compiled[na_rows, c("siteid", "cmems_id", "wrf_id", "yday")]
        d.ls$compiled[[j]][na_rows] <- left_join(d_meta |>
                                                   select(all_of(env_id_col)),
                                                 yday_env[[i]] |>
                                                   select(all_of(c(env_id_col, j))))[[j]]
      }
    }
  }
  if(target=="tox") {
    d.ls$compiled <- d.ls$compiled |>
      left_join(d.ls$habAvg |> select(-date, -siteid), by="obsid")
  }

  return(d.ls)
}

