% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/assembles.R
\name{assemble_sentinel_counts}
\alias{assemble_sentinel_counts}
\title{Assemble counts of transmissions/detections from sentinel tags}
\usage{
assemble_sentinel_counts(
  sentinel,
  moorings,
  breaks = "hours",
  as_POSIXct = function(x, tz = "UTC", ...) fasttime::fastPOSIXct(x, tz = tz, ...),
  detection_range = 2000,
  dist_btw_receivers = NULL
)
}
\arguments{
\item{sentinel}{A dataframe that includes all sentinel tag transmissions and detections. This should include the following columns: `timestamp', the time of each observation; `type', the type of each observation ("transmission" or "detection"); `sink_id', a unique identifier of the receivers which received the transmissions; and `source_id', a unique identifier of the receivers from which the transmission was released (by the built-in sentinel tag). (When the `source_id' is the same as the `sink_id', the observation is of type "transmission"; otherwise, the observation must be of type "detection"). See \code{\link[flapper]{dat_sentinel}} for an example.}

\item{moorings}{A dataframe that includes receiver identifiers, deployment times and locations. This information is necessary so that, for each time window, only receivers that could have detected transmissions (i.e., were nearby to the transmitter and active at the time of transmission) are included in the processed dataframe. This should include the following columns: `receiver_id', a unique identifier for each receiver; `receiver_start_date', the deployment date for each receiver; `receiver_end_date', the deployment end date for each receiver; `receiver_long', the receiver longitude (decimal degrees); and `receiver_lat', the receiver latitude (decimal degrees). See \code{\link[flapper]{dat_moorings}} for an example.}

\item{breaks}{The time interval over which transmissions/detections are counted (e.g. \code{"hours"}). This is passed to the \code{breaks} argument of \code{\link[base]{cut.POSIXt}}.}

\item{as_POSIXct}{A function that defines how to convert character time bins (returned by \code{\link[base]{cut.POSIXt}}) to POSIXct format. Usually, \code{\link[base]{as.POSIXct}} is suitable but, for large time series, this is slow. In this case, \code{\link[fasttime]{fastPOSIXct}} can be used to improve algorithm speed, if suitable. The default function is \code{function(x, tz = "UTC",...) fasttime::fastPOSIXct(x, tz = tz,...)}.}

\item{detection_range}{A number that defines the maximum distance (m) between a transmitting receiver and other receivers within with detections may plausibly be received. For each transmission, all detections at active receivers within this distance are counted; other receivers are not considered.}

\item{dist_btw_receivers}{(optional) A dataframe that specifies the distances between all combinations of receivers. If not provided, this is computed internally by \code{\link[flapper]{dist_btw_receivers}}, which assumes Euclidean distances are appropriate. If provided, the dataframe should contain columns: `r1', `r2', `dist' (see \code{\link[flapper]{dist_btw_receivers}}). Note that, if provided, the `dist' column should be in m, not km, as returned by default by \code{\link[flapper]{dist_btw_receivers}}, to match the units of \code{detection_range}.}
}
\value{
The function returns a dataframe that, for each source receiver (i.e., sentinel tag), specifies the number of transmissions from that receiver and the number of detections at each nearby receiver over the specified time window. The dataframe has the following columns: `timestamp_bin', the time window; `source_id', the identifier of the source of transmission; `sink_id', the identifier of each potential recipient of each transmission; `n_trms', the number of transmissions from the source over each time window; `n_dets', the number of detections of each transmission at each potential recipient receiver; `dist_btw_receivers', the distance between the source and the sink receiver. Rows are ordered by `source_id', then `timestamp' and finally `sink_id'.
}
\description{
This function assembles a dataframe of counts of transmissions/detections from sentinel tags over user-defined time windows. For each sentinel tag, in each time window, the function counts the number of transmissions of that tag (i.e., the expected number of detections at nearby receivers). At the same time, the function determines the number of detections at all nearby receivers (defined as those within a user-specified distance) that were active at the same time as the sentinel tag. The result is a dataframe which comprises, for each tag and each time window, the number of transmissions (i.e., the expected number of detections) and the observed number of detections at each nearby receiver, arranged in long format.
}
\examples{
#### Example (1): Use default options and example dataframes, which contain
# ... all required columns:
sentinel_counts <- assemble_sentinel_counts(
  sentinel = dat_sentinel,
  moorings = dat_moorings
)
head(sentinel_counts)
tail(sentinel_counts)
#### Example (2): Adjust time window
sentinel_counts <- assemble_sentinel_counts(
  sentinel = dat_sentinel,
  moorings = dat_moorings,
  breaks = "days"
)
head(sentinel_counts)
tail(sentinel_counts)
#### Example (2): Adjust maximum detection distance and supply dist_btw_receivers
dist_btw_receivers_m <- dist_btw_receivers(
  dat_moorings[, c(
    "receiver_id",
    "receiver_long",
    "receiver_lat"
  )],
  f = function(x) {
    return(x * 1000)
  }
)
sentinel_counts <- assemble_sentinel_counts(
  sentinel = dat_sentinel,
  moorings = dat_moorings,
  detection_range = 1500,
  dist_btw_receivers = dist_btw_receivers_m
)
head(sentinel_counts)
tail(sentinel_counts)

}
\author{
Edward Lavender
}
