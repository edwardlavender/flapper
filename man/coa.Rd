% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/coa.R
\name{coa}
\alias{coa}
\title{The centres of activity (COA) algorithm}
\usage{
coa(
  mat,
  xy,
  cl = NULL,
  varlist = NULL,
  na_omit = TRUE,
  as_POSIXct = as.POSIXct,
  output = "matrix"
)
}
\arguments{
\item{mat}{A detection matrix, with one row for each time step and one column for each receiver, in which each cell defines the number of detections at each time step/receiver (for a particular individual) (see \code{\link[flapper]{make_matrix_detections}}). It is advisable that the rows and columns of this matrix are labelled by time stamp and receiver respectively, especially if there are any rows without detections (see \code{output} below).}

\item{xy}{A matrix that defines receiver locations (x, y). This should contain one row for each receiver (column) in \code{mat} (in the same order as in \code{mat}) and two columns for the coordinates. Planar coordinates (i.e., in Universal Transverse Mercator projection) are required for the averaging process.}

\item{cl, varlist}{(optional) Parallelisation options. \code{cl} is (a) a cluster object from \code{\link[parallel]{makeCluster}} or (b) an integer that defines the number of child processes. \code{varlist} is a character vector of variables for export (see \code{\link[flapper]{cl_export}}). Exported variables must be located in the global environment. If a cluster is supplied, the connection to the cluster is closed within the function (see \code{\link[flapper]{cl_stop}}). For further information, see \code{\link[flapper]{cl_lapply}} and \code{\link[flapper]{flapper-tips-parallel}}.}

\item{na_omit, as_POSIXct}{Processing options. \code{na_omit} is a logical variable that defines whether or not to omit NAs (i.e., rows in \code{mat} for which no detections were made and thus for which COAs cannot be calculated) from returned coordinates. If \code{output = "data.frame"} (see below), \code{as_POSIXct} is a function to convert time stamps, taken from the row names of \code{mat}, to a POSIXct vector. \code{as_POSIXct = NULL} suppresses this conversion.}

\item{output}{A character that defines the output format. Currently supported options are: \code{"matrix"}, which returns a matrix of the coordinates of COAs; and \code{"data.frame"}, which returns a dataframe with timestamps (taken from the row names of \code{mat}) and COA coordinates.}
}
\value{
The function returns a matrix or a dataframe, depending on the \code{output} argument, that represents a time series of COAs. If \code{na_omit = TRUE}, the time series may have `gaps' over which COAs could not be calculated due to the absence of detections.
}
\description{
This function implements the arithmetic mean-position algorithm to calculate an individual's centres of activity (COAs) though time from detections at passive acoustic telemetry receivers, as described by Simpfendorfer et al (2002). Under this approach, each COA is calculated as the arithmetic mean of the locations of receivers at which an individual was detected over a specified time interval, weighted by the frequency of detections at each of those receivers. To implement the function, a detection matrix that defines the number of detections of an individual along sequence of time steps at each receiver (e.g., from \code{\link[flapper]{make_matrix_detections}}) needs to be supplied, along a matrix of receiver locations in a planar coordinate system. For each time interval, the function calculates the centre of activity and returns a matrix or dataframe with this information.
}
\details{
Centres of activity (COA) are a widely used metric for the reconstruction of patterns of space use from passive acoustic telemetry detections. Several methods have been developed to calculate COAs, but the mean-position algorithm is the commonest. Under this approach, COAs are estimated as an average of the locations of receivers at which an individual is detected over a specified time interval, weighted by the frequency of detections at each of those receivers. Within \code{\link[flapper]{flapper}}, COAs are calculated in three stages by first exploring possible time intervals over which to calculate COAs with \code{\link[flapper]{coa_setup_delta_t}}; then summarising detections over those intervals with \code{\link[flapper]{make_matrix_detections}}; and finally passing the resultant detection matrix to the \code{\link[flapper]{coa}} function to calculate COAs. This implements the arithmetic version of the mean-position algorithm, calculating the arithmetic mean of the receiver locations, weighted by the frequency of detections at each receiver.

To generate estimates of space use, COAs are usually taken as point estimates from which utilisation distributions (typically kernel utilisation distributions, KUDs) are estimated. Thus, in the case of a coupled COA-KUD approach, usually the estimate of space use is a (kernel) utilisation distribution which describes the probability of relocating an individual in any given area at a randomly chosen time. Alternative methods of home range analysis, including those which incorporate time, such as dynamic Brownian bridge movement models, can be used to estimate the utilisation distribution. Generally, the COA approach is most suitable when detections are relatively frequent, and receivers are regularly distributed across an area. Under other conditions, its performance as a method for estimating space use has been subject to relatively limited evaluation but it can be problematic (e.g., in clustered arrays).
}
\examples{
#### Define data for the calculation of COAs

## (1) Define the period over which to calculate COAs
# ... by focusing on the time over which IDs were at liberty.
dat_ids$tag_start_date <- as.POSIXct(dat_ids$tag_start_date)
dat_ids$tag_end_date <- as.POSIXct("2017-06-02")

## (2) Define the receivers over which to calculate detections
# ... Here we use factors to ensure that the order of receiver coordinates
# ... (see below) and the order or receivers in the detection matrix
# ... matches.
dat_moorings$receiver_id <- factor(dat_moorings$receiver_id)
dat_acoustics$receiver_id <- factor(dat_acoustics$receiver_id,
  levels = levels(dat_moorings$receiver_id)
)

## (3) Define the detection matrix
# ... Here we simply create a detection matrix across all IDs and using
# ... a convenient (but unjustified) delta_t value. In reality, we would need
# ... to consider the time series for which to calculate
# ... COAs and the appropriate delta t value(s) more carefully.
detection_matrix_by_id <- make_matrix_detections(dat_acoustics,
  delta_t = "days",
  start = min(dat_ids$tag_start_date),
  end = max(dat_ids$tag_end_date)
)

## (4) Define receiver coordinates with UTM projection
proj_wgs84 <- sp::CRS(SRS_string = "EPSG:4326")
proj_utm <- sp::CRS(SRS_string = "EPSG:32629")
xy <- sp::SpatialPoints(
  dat_moorings[, c("receiver_long", "receiver_lat")],
  proj_wgs84
)
xy <- sp::spTransform(xy, proj_utm)
xy <- sp::coordinates(xy)

#### Example (1): Implement the COA algorithm for an example individual
coa_mat <- coa(mat = detection_matrix_by_id[[1]], xy = xy)
utils::str(coa_mat)

#### Example (2): Change the output format and coerce time stamps to POSIXct format
coa_dat <- coa(mat = detection_matrix_by_id[[1]], xy = xy, output = "data.frame")
utils::str(coa_dat)

#### Example (3): Implement the algorithm on a cluster
# This will only be faster for very large detection time series.
if (flapper_run_parallel) {
  coa_mat <- coa(
    mat = detection_matrix_by_id[[1]],
    xy = xy,
    cl = parallel::makeCluster(2L)
  )
}

}
\references{
Simpfendorfer, C. A., M. R. Heupel, and R. E. Hueter. 2002. Estimation of short-term centers of activity from an array of omnidirectional hydrophones and its use in studying animal movements. Canadian Journal of Fisheries and Aquatic Sciences 59:23-32.
}
\seealso{
\code{\link[flapper]{coa_setup_delta_t}} suggests suitable time intervals over which to calculate COAs. \code{\link[flapper]{make_matrix_detections}} makes the detection matrices from detection time series data required by this function. For data in the VEMCO Vue export format, the `COA' function in the VTrack package (https://github.com/RossDwyer/VTrack) can also be used to calculate centres of activity.
}
\author{
Edward Lavender
}
