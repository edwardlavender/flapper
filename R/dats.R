#####################################
#####################################
#### dat_ids

#' @title Example tagged individuals dataset
#' @description  A dataset containing the characteristics of a sample of flapper skate (\emph{Dipturus intermedius}) that were tagged as part of the Movement Ecology of Flapper Skate (MEFS) project. Data are arranged by \code{individual_id} (see below).
#'
#' @format A dataframe with 3 observations and 6 variables:
#' \describe{
#'   \item{individual_id}{A number which distinguishes each unique individual.}
#'   \item{transmitter_id}{A character which defines the acoustic transmitter code.}
#'   \item{archival_id}{A character which defines the archival tag code.}
#'   \item{tag_start_date}{A date which defines the date of tag deployment.}
#'   \item{tag_long}{A number which defines the longitude (decimal degrees) of tag deployment.}
#'   \item{tag_lat}{A number which defines the latitude (decimal degrees) of tag deployment.}
#' }
#'
#'
#' @source Data were collected by Marine Scotland Science and Scottish Natural Heritage. Data were processed by Edward Lavender.
"dat_ids"



#####################################
#####################################
#### dat_moorings

#' @title Example passive acoustic telemetry receiver moorings dataset
#' @description A dataset containing a sample of passive acoustic telemetry receiver locations and associated information. Data are arranged by \code{receiver_id} (see below).
#'
#' @format A dataframe with 39 observations and 8 variables:
#' \describe{
#'   \item{receiver_id}{A number which distinguishes each unique receiver deployment.}
#'   \item{receiver}{A character which defines each receiver.}
#'   \item{sentinel_id}{A character which defines each receiver's built-in acoustic transmitter ('sentinel tag'). Some receivers did not have sentinel tags; their \code{sentinel_id} is \code{NA}.}
#'   \item{receiver_start_date}{A date which defines the start date of the receiver's deployment.}
#'   \item{receiver_end_date}{A date which defines the end date of the receiver's deployment.}
#'   \item{receiver_long}{A number which defines the latitude (decimal degrees) of the receiver.}
#'   \item{receiver_lat}{A number which defines the latitude (decimal degrees) of the receiver.}
#'   \item{receiver_depth}{A number which defines the approximate depth of each receiver below the surface (m).}
#' }
#'
#' @source Data were collected by Marine Scotland Science and Scottish Natural Heritage. Data were processed by Edward Lavender.
"dat_moorings"


#####################################
#####################################
#### dat_acoustics

#' @title Example passive acoustic telemetry dataset
#' @description A dataset containing a sample of processed flapper skate (\emph{Dipturus intermedius}) detection timeseries. Data are arranged by \code{individual_id}, \code{timestamp} and then \code{receiver_id} (see below).
#'
#' @format A dataframe with 59,420 observations and 8 variables:
#' \describe{
#'   \item{timestamp}{A POSIXct object which defines the time of each observation.}
#'   \item{individual_id}{A unique identifier of the individual that was detected.}
#'   \item{transmitter_id}{The acoustic transmitter that was detected (see \code{\link{dat_ids}}).}
#'   \item{receiver_id}{A unique identifier of the receiver at which the individual was detected (see \code{\link{dat_moorings}}).}
#'   \item{receiver}{The receiver at which the individual was detected (see \code{\link{dat_moorings}}).}
#'   \item{receiver_long}{The longitude of the receiver at which the individual was detected (see \code{\link{dat_moorings}}).}
#'   \item{receiver_lat}{The latitude of the receiver at which the individual was detected (see \code{\link{dat_moorings}}).}
#'   \item{receiver_depth}{The depth of the receiver at which the individual was detected (see \code{\link{dat_moorings}}).}
#' }
#'
#' @source Data were collected by Marine Scotland Science and Scottish Natural Heritage. Data were processed by Edward Lavender.
"dat_acoustics"


#####################################
#####################################
#### dat_sentinel

#' @title Example sentinel tag range testing dataset
#' @description A dataset containing a sample of transmissions and detections assembled from sentinel tags. Sentinel tags are built-in acoustic transmitters which release acoustic signals at programmed intervals. The receiver unit of which each sentinel tag is a part records the timing of these transmissions. At the same time, any nearby receivers will also record transmissions from these acoustic tags. After some processing, the result is a dataframe comprising transmissions and detections, which contains information on detection probability. Data are arranged by \code{source_id}, \code{timestamp} and then \code{sink_id}.
#'
#' @format A dataframe with 106,754 observations and 7 variables:
#' \describe{
#'   \item{timestamp}{A POSIXct object which defines the time of each observation.}
#'   \item{type}{A character which defines the type of observation: "transmission" or "detection". "transmission" corresponds }
#'   \item{sink_id}{A unique identifier of the receiver which received the transmission (see \code{\link{dat_moorings}}).}
#'   \item{sink}{The receiver which received the transmission (see \code{\link{dat_moorings}}).}
#'   \item{source_id}{A unique identifier of the receiver which released the transmission (see \code{\link{dat_moorings}}).}
#'   \item{source}{The receiver which released the transmission (see \code{\link{dat_moorings}}).}
#'   \item{sentinel_id}{A unique identifier of the sentinel tag which released the transmission (see \code{\link{dat_moorings}}).}
#' }
#'
#' @source Data were collected by Marine Scotland Science and Scottish Natural Heritage. Data were processed by Edward Lavender.
"dat_sentinel"


#####################################
#####################################
#### dat_archival

#' @title Example archival dataset
#' @description A dataset containing a sample of flapper skate (\emph{Dipturus intermedius}) depth (m) timeseries. Observations were sampled every 2 minutes using archival tags. Data are arranged by \code{individual_id} and then \code{timestamp}.
#'
#' @format A dataframe with 75,000 observations and 5 variables:
#' \describe{
#'   \item{timestamp}{A POSIXct object which defines the time of each observation.}
#'   \item{individual_id}{A number which defines each individual (see \code{\link{dat_ids}}).}
#'   \item{transmitter_id}{The individual's acoustic transmitter (see \code{\link{dat_ids}}).}
#'   \item{archival_id}{The individual's archival tag (see \code{\link{dat_ids}}).}
#'   \item{depth}{A number which defines the depth (m) of the individual at each timestep}
#' }
#'
#' @source Data were collected by Marine Scotland Science and Scottish Natural Heritage. Data were processed by Edward Lavender.
"dat_archival"


#### End of code.
#####################################
#####################################
