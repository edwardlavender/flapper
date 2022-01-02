#####################################
#####################################
#### dat_ids

#' @title Example tagged individuals dataset
#' @description  A dataset containing the characteristics of a sample of flapper skate (\emph{Dipturus intermedius}) that were tagged as part of the Movement Ecology of Flapper Skate (MEFS) project. Data are arranged by \code{individual_id} (see below).
#'
#' @format A dataframe with 33 observations and 6 variables:
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
#' @source Data were collected by, and belong to, Marine Scotland Science and NatureScot. Data were processed by Edward Lavender. If you wish to use these data, please contact Marine Scotland Science and NatureScot for further information.
#'
#' @references Data collection and processing is described in Lavender (in prep). Fine-scale habitat use of the Critically Endangered flapper skate (\emph{Dipturus intermedius}). [Doctoral dissertation, University of St Andrews].
"dat_ids"


#####################################
#####################################
#### dat_moorings

#' @title Example passive acoustic telemetry receiver moorings dataset
#' @description A dataset containing a sample of passive acoustic telemetry receiver locations and associated information. Data are arranged by \code{receiver_id} (see below).
#'
#' @format A dataframe with 40 observations and 8 variables:
#' \describe{
#'   \item{receiver_id}{An integer which distinguishes each unique receiver deployment.}
#'   \item{receiver}{A character which defines each receiver.}
#'   \item{sentinel_id}{A character which defines each receiver's built-in acoustic transmitter (`sentinel tag'). Some receivers did not have sentinel tags; their \code{sentinel_id} is \code{NA}.}
#'   \item{receiver_start_date}{A date which defines the start date of each receiver's deployment.}
#'   \item{receiver_end_date}{A date which defines the end date of each receiver's deployment.}
#'   \item{receiver_long}{A number which defines the longitude (decimal degrees) of each receiver.}
#'   \item{receiver_lat}{A number which defines the latitude (decimal degrees) of each receiver.}
#'   \item{receiver_depth}{A number which defines the approximate depth of each receiver below the surface (m).}
#' }
#'
#' @source Data were collected by, and belong to, Marine Scotland Science and NatureScot. Data were processed by Edward Lavender. If you wish to use these data, please contact Marine Scotland Science and NatureScot for further information.
#'
#' @references Data collection and processing is described in Lavender (in prep). Fine-scale habitat use of the Critically Endangered flapper skate (\emph{Dipturus intermedius}). [Doctoral dissertation, University of St Andrews].
"dat_moorings"


#####################################
#####################################
#### dat_acoustics

#' @title Example passive acoustic telemetry detections dataset
#' @description A dataset containing a sample of processed flapper skate (\emph{Dipturus intermedius}) detection time series. Data are arranged by \code{individual_id}, \code{timestamp} and then \code{receiver_id} (see below).
#'
#' @format A dataframe with 39,242 observations and 9 variables:
#' \describe{
#'   \item{individual_id}{A unique identifier of the individual that was detected (see \code{\link{dat_ids}}).}
#'   \item{transmitter_id}{The acoustic transmitter that was detected (see \code{\link{dat_ids}}).}
#'   \item{index}{An integer individual-specific cumulative detection index.}
#'   \item{timestamp}{A POSIXct object which defines the time of each observation.}
#'   \item{receiver_id}{A unique identifier of the receiver at which the individual was detected (see \code{\link{dat_moorings}}).}
#'   \item{receiver}{The receiver at which the individual was detected (see \code{\link{dat_moorings}}).}
#'   \item{receiver_long}{The longitude of the receiver at which the individual was detected (see \code{\link{dat_moorings}}).}
#'   \item{receiver_lat}{The latitude of the receiver at which the individual was detected (see \code{\link{dat_moorings}}).}
#'   \item{receiver_depth}{The depth of the receiver at which the individual was detected (see \code{\link{dat_moorings}}).}
#' }
#'
#' @details Time stamps are expressed to the nearest two minutes. `Duplicate' detections (of the same individual at the same receiver in the same time step) are excluded. This supports the implementation of the acoustic-centroid (AC) and acoustic-centroid depth-contour (ACDC) algorithms (see \code{\link[flapper]{ac}} and \code{\link[flapper]{acdc}}).
#'
#' @source Data were collected by, and belong to, Marine Scotland Science and NatureScot. Data were processed by Edward Lavender. If you wish to use these data, please contact Marine Scotland Science and NatureScot for further information.
#'
#' @references Data collection and processing is described in Lavender (in prep). Fine-scale habitat use of the Critically Endangered flapper skate (\emph{Dipturus intermedius}). [Doctoral dissertation, University of St Andrews].
"dat_acoustics"


#####################################
#####################################
#### dat_sentinel

#' @title Example sentinel tag range testing dataset
#' @description A dataset containing a sample of transmissions and detections assembled from sentinel tags. Sentinel tags are built-in acoustic transmitters which release acoustic signals at programmed intervals. The receiver unit of each sentinel tag records the exact timing of these transmissions. At the same time, any nearby receivers may also record these transmissions as detections. After some processing, the result is a dataframe comprising transmissions and detections, which contains information on detection probability. Data are arranged by \code{source_id}, \code{timestamp} and then \code{sink_id}.
#'
#' @format A dataframe with 106,733 observations and 7 variables:
#' \describe{
#'   \item{timestamp}{A POSIXct object which defines the time of each observation.}
#'   \item{type}{A character which defines the type of observation: "transmission" or "detection".}
#'   \item{sink_id}{A unique identifier of the receiver which received the transmission (see \code{\link{dat_moorings}}).}
#'   \item{sink}{The receiver which received the transmission (see \code{\link{dat_moorings}}).}
#'   \item{source_id}{A unique identifier of the receiver which released the transmission (see \code{\link{dat_moorings}}).}
#'   \item{source}{The receiver which released the transmission (see \code{\link{dat_moorings}}).}
#'   \item{sentinel_id}{A unique identifier of the sentinel tag which released the transmission (see \code{\link{dat_moorings}}).}
#' }
#'
#' @source Data were collected by, and belong to, Marine Scotland Science and NatureScot. Data were processed by Edward Lavender. If you wish to use these data, please contact Marine Scotland Science and NatureScot for further information.
#'
#' @references Data collection and processing is described in Lavender (in prep). Fine-scale habitat use of the Critically Endangered flapper skate (\emph{Dipturus intermedius}). [Doctoral dissertation, University of St Andrews].
"dat_sentinel"


#####################################
#####################################
#### dat_archival

#' @title Example archival dataset
#' @description A dataset containing a sample of flapper skate (\emph{Dipturus intermedius}) depth (m) time series. Observations were sampled every 2 minutes using archival tags. Data are arranged by \code{individual_id} and then \code{timestamp}.
#'
#' @format A dataframe with 75,000 observations and 5 variables:
#' \describe{
#'   \item{timestamp}{A POSIXct object which defines the time of each observation.}
#'   \item{individual_id}{A number which defines each individual (see \code{\link{dat_ids}}).}
#'   \item{transmitter_id}{The individual's acoustic transmitter (see \code{\link{dat_ids}}).}
#'   \item{archival_id}{The individual's archival tag (see \code{\link{dat_ids}}).}
#'   \item{depth}{A number which defines the depth (m) of the individual at each time step.}
#' }
#'
#' @source Data were collected by, and belong to, Marine Scotland Science and NatureScot. Data were processed by Edward Lavender. If you wish to use these data, please contact Marine Scotland Science and NatureScot for further information.
#'
#' @references Data collection and processing is described in Lavender (in prep). Fine-scale habitat use of the Critically Endangered flapper skate (\emph{Dipturus intermedius}). [Doctoral dissertation, University of St Andrews].
"dat_archival"

#####################################
#####################################
#### dat_coast

#' @title The coastline around the MEFS Firth of Lorn acoustic array
#' @description A SpatialPolygonsDataFrame delineating the coastline around a subset of acoustic receivers set up by the Movement Ecology of Flapper Skate (MEFS) project in the Firth of Lorn, off the west coast of Scotland.
#'
#' @format A SpatialPolygonsDataFrame (see \code{\link[sp]{SpatialPolygonsDataFrame-class}}).
#'
#' @source Data were processed from https://biogeo.ucdavis.edu/data/gadm3.6/Rsp/gadm36_GBR_0_sp.rds.
"dat_coast"


#####################################
#####################################
#### dat_gebco

#' @title The bathymetry around the MEFS Firth of Lorn acoustic array
#' @description A dataset of the bathymetry (m) in an area around a subset of acoustic receivers set up by the Movement Ecology of Flapper Skate (MEFS) project in the Firth of Lorn, off the west coast of Scotland. Bathymetry data are provided by the General Bathymetric Chart of the Oceans (GEBCO).
#'
#' @format A \code{\link[raster]{raster}} with 36 rows, 36 columns and 1296 cells, with the following properties:
#' \describe{
#'   \item{dimensions}{57, 74, 4218 (nrow, ncol, ncell)}
#'   \item{resolution}{257, 463  (x, y)}
#'   \item{extent}{695492.1, 714510.1, 6246657, 6273048 (xmin, xmax, ymin, ymax)}
#'   \item{crs}{+proj=utm +zone=29 +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs}
#' }
#'
#' @source GEBCO Compilation Group (2019) GEBCO 2019 Grid (doi:10.5285/836f016a-33be-6ddc-e053-6c86abc0788e)
"dat_gebco"


#####################################
#####################################
#### dat_centroids

#' @title Example detection centroids from \code{\link[flapper]{acs_setup_centroids}}
#' @description A list of detection centroids created by \code{\link[flapper]{acs_setup_centroids}}. This is included principally to streamline function examples.
#'
#' @format A list with 57 elements, one for each receiver from 1:max(\code{\link[flapper]{dat_moorings}$receiver_id}). Each element contains a SpatialPolygonsDataFrame with the detection centroid for that receiver, under a detection range of 425 m and within the boundaries defined by \code{\link[flapper]{dat_coast}} and \code{\link[flapper]{dat_gebco}}.
"dat_centroids"


#####################################
#####################################
#### dat_dc

#' @title Example DC algorithm output
#' @description An object of class \code{\link[flapper]{acdc_archive-class}} from \code{\link[flapper]{dc}} created by a step-wise implementation of this function. This is included principally to streamline function examples.
#'
#' @format A named list with 4 elements:
#' \describe{
#'   \item{archive}{A list of results from internal function calls that contain essential results, with one element per chunk.}
#'   \item{ts_by_chunk}{A list of time series, with one element per chunk.}
#'   \item{time}{A dataframe that defines the times of sequential stages in the algorithm's progression.}
#'   \item{args}{A named list of user inputs that record the parameters used to generate the outputs.}
#' }
#'
#' @seealso See \code{\link[flapper]{acdc_archive-class}} for further information on this S3 class.
"dat_dc"


#####################################
#####################################
#### dat_acdc

#' @title Example ACDC algorithm output
#' @description An object of class \code{\link[flapper]{acdc_archive-class}} from \code{\link[flapper]{acdc}}, created by a (2-hour) chunk-wise implementation of this function. This is included principally to streamline function examples.
#'
#' @format A named list with 4 elements:
#' \describe{
#'   \item{archive}{A list of results from calls to \code{\link[flapper]{.acs}}, with one element per chunk.}
#'   \item{ts_by_chunk}{A list of time series, with one element per chunk.}
#'   \item{time}{A dataframe that defines the times of sequential stages in the algorithm's progression.}
#'   \item{args}{A named list of user inputs that record the parameters used to generate the outputs.}
#' }
#'
#' @seealso See \code{\link[flapper]{acdc_archive-class}} for further information on this S3 class.
"dat_acdc"


#####################################
#####################################
#### dat_dcpf_histories

#' @title Example DCPF algorithm output: particle histories
#' @description An example \code{\link[flapper]{pf_archive-class}} object. This comprises a named list of results from an implementation of the depth-contour particle filtering (DCPF) algorithm (see \code{\link[flapper]{dc}} and \code{\link[flapper]{pf}}). The algorithm was implemented for a series of depth observations collected from flapper skate (\emph{Dipturus intermedius}) off the west coast of Scotland, where the bathymetry is described by \code{\link[flapper]{dat_gebco}}. The outputs include algorithm parameters and a list of dataframes that record the particles sampled at each time step. This is included principally to streamline function examples.
#'
#' @format A named list, with two elements, which records the parameters used in the call to \code{\link[flapper]{pf}} and the particles sampled at each time step:
#' \describe{
#'   \item{args}{A named list that records the function arguments used to generate outputs.}
#'   \item{history}{A list of dataframes, one for each time step, that record particle samples. Each dataframe comprises \code{n = 10} rows and the following three columns:}
#'   \itemize{
#'     \item{id_current}{ An integer that uniquely defines each location on the bathymetry \code{\link[raster]{raster}} (\code{bathy}).}
#'     \item{pr_current} A double that defines the probability of movement into each cell.
#'     \item{timestep}{ An integer that defines each time step.}
#'   }
#' }
#'
"dat_dcpf_histories"


#####################################
#####################################
#### dat_dcpf_paths

#' @title Example DCPF algorithm output: reconstructed paths
#' @description An example \code{\link[flapper]{pf_path-class}} object. This comprises a dataframe of movement paths reconstructed by an implementation of the depth-contour particle filtering (DCPF) algorithm (see \code{\link[flapper]{dc}} and \code{\link[flapper]{pf}}). The algorithm was implemented for a series of depth observations collected from flapper skate (\emph{Dipturus intermedius}) off the west coast of Scotland, where the bathymetry is described by \code{\link[flapper]{dat_gebco}}. The resultant particle histories (\code{\link[flapper]{dat_dcpf_histories}}) were processed into this set of movement paths by \code{\link[flapper]{pf_simplify}}. This dataset is included principally to streamline function examples.
#'
#' @format A dataframe with seven columns and 180 observations that records the reconstructed paths. This contains the following columns:
#' \describe{
#'     \item{path_id}{ An integer that uniquely defines each path.}
#'     \item{timestep}{ An integer that defines each time step.}
#'     \item{cell_id}{ An integer that defines the cells ID of the surface raster (over which paths were reconstructed) retained by the algorithm.}
#'     \item{cell_x}{ A double that defines the cell x coordinate.}
#'     \item{cell_y}{ A double that defines the cell y coordinate.}
#'     \item{cell_z}{ A double that defines the value of the surface raster in each cell.}
#'     \item{cell_pr}{ A double that defines the probability of movement into each cell, given immediately preceding sampled locations (or the origin).}
#'   }
#'
"dat_dcpf_paths"


#### End of code.
#####################################
#####################################
