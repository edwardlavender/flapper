% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dats.R
\docType{data}
\name{dat_acoustics}
\alias{dat_acoustics}
\title{Example passive acoustic telemetry dataset}
\format{A dataframe with 59,420 observations and 8 variables:
\describe{
  \item{timestamp}{A POSIXct object which defines the time of each observation.}
  \item{individual_id}{A unique identifier of the individual that was detected.}
  \item{transmitter_id}{The acoustic transmitter that was detected (see \code{\link{dat_ids}}).}
  \item{receiver_id}{A unique identifier of the receiver at which the individual was detected (see \code{\link{dat_moorings}}).}
  \item{receiver}{The receiver at which the individual was detected (see \code{\link{dat_moorings}}).}
  \item{receiver_long}{The longitude of the receiver at which the individual was detected (see \code{\link{dat_moorings}}).}
  \item{receiver_lat}{The latitude of the receiver at which the individual was detected (see \code{\link{dat_moorings}}).}
  \item{receiver_depth}{The depth of the receiver at which the individual was detected (see \code{\link{dat_moorings}}).}
}}
\source{
Data were collected by Marine Scotland Science and Scottish Natural Heritage. Data were processed by Edward Lavender.
}
\usage{
dat_acoustics
}
\description{
A dataset containing a sample of processed flapper skate (\emph{Dipturus intermedius}) detection timeseries. Data are arranged by \code{individual_id}, \code{timestamp} and then \code{receiver_id} (see below).
}
\keyword{datasets}