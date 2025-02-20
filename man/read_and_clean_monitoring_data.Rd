% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read_and_clean_sites.R
\name{read_and_clean_monitoring_data}
\alias{read_and_clean_monitoring_data}
\title{Read and clean monitoring data}
\usage{
read_and_clean_monitoring_data(i, urls, targ_i, sites, dateStart)
}
\arguments{
\item{i}{Monitoring target type: "hab", "tox", or "habfish"}

\item{urls}{named list of urls}

\item{targ_i}{named list of target info}

\item{sites}{monitoring sites from read_and_clean_sites()}

\item{dateStart}{ymd start date}
}
\value{
Dataframe with corresponding monitoring data
}
\description{
Wrapper for other read_and_clean functions
}
