######################################
######################################
#### acdc_record-class

#' @title "acdc_record" class
#' @description An S3 class that defines the object returned by an acoustic-container/depth-contour (AC/DC) algorithm (\code{\link[flapper]{ac}}, \code{\link[flapper]{dc}} or \code{\link[flapper]{acdc}}), after simplification via \code{\link[flapper]{acdc_simplify}} or directly from internal routines.

#' @return \subsection{A named list}{An "acdc_record" class object is a named list with the following elements: `map', `record', `time', `args', `chunks' and `simplify'. The main output element is the `map' RasterLayer that shows where the individual could have spent more or less time over the duration of the movement time series. The `record' element records time-specific maps of the possible locations of the individual, and can be used to plot maps of specific time points or to produce animations (for the time steps specified by the \code{save_record_spatial} argument). The `time' element is a dataframe that defines the times of sequential stages in the algorithm's progression, providing a record of computation time. The `args' element is a named list of user inputs that record the parameters used to generate the outputs (if \code{save_args = TRUE}, otherwise the `args' element is \code{NULL}). The `chunks' element is a list with chunk-specific information that is generated if the algorithm is implemented chunk-wise and then simplified via \code{\link[flapper]{acdc_simplify}} with the \code{keep_chunks = TRUE} argument. The \code{simplify} element is a logical value that defines whether or not the object was created from \code{\link[flapper]{ac}}/\code{\link[flapper]{dc}}/\code{\link[flapper]{acdc}} and \code{\link[flapper]{acdc_simplify}} or an internal routine. Below, more detail about the `map' and `record' elements is provided.}
#'
#' \subsection{(A) `map'}{This is a RasterLayer that defines the locations in which an individual could have been (or was not), given the movement time series. Each cell is a count of the number of times when the acoustic and/or archival data were compatible with the individual having been in that cell, given the algorithm (unless normalisation has been implemented via \code{normalise = TRUE} in the algorithm or in  post-processing via \code{\link[flapper]{acdc_simplify}}). This provides an overall measure of the locations in which an individual could have spent more or less time, but not where it was.}
#'
#' \subsection{(B) `record'}{This is a temporal record of the possible locations of each individual, with one element for each primary time step. For the AC* algorithms (\code{\link[flapper]{ac}} and \code{\link[flapper]{acdc}}), the `primary time step' refers to acoustic time steps; for the DC algorithm (\code{\link[flapper]{dc}}), the `primary time step' refers to archival time steps. Each element contains (a) a dataframe with the information for that time step and (b) a list of spatial objects, one for each `secondary time step' that record spatial outputs of the algorithm, as specified by \code{save_record_spatial} (see below). For the AC* algorithms, `secondary time steps' are the time steps in-between acoustic observations; for the DC algorithm, primary and secondary time steps are the same.}
#'
#' \subsection{(i) record$dat}{
#'
#' For the AC* algorithms, \code{record$dat} contains the following information:
#'
#' \itemize{
#'   \item `timestep_cumulative' is a cumulative vector that defines the time step of each observation (across all acoustic and archival time steps up to the current time step);
#'   \item `timestep_detection' defines the time step of the detection time series;
#'   \item `timestep_archival' defines the time steps between the current and next acoustic detection (if applicable);
#'   \item `receiver_1_id' and `receiver_2_id' define the receiver at which the individual was detected at the current and next acoustic time steps;
#'   \item `receiver_1_timestamp' and `receiver_2_timestamp' define the time stamps of these detections;
#'   \item `time_btw_dets' provides the time (s) between these acoustic detections;
#'   \item `archival_timestamp' defines the time stamp of archival observations between acoustic detections (if applicable);
#'   \item `archival_depth' defines the depth (m) of the individual each archival time step (if applicable);
#'   }
#'
#' For the DC algorithm, \code{record$dat} contains the \code{archival} data for the time step, as inputted, with the following additional columns:
#' \itemize{
#'   \item `timestep_cumulative' is a cumulative integer that defines the time step;
#'   \item `depth_lwr' is a number that defines the lower bound for the individual's possible depth on the bathymetry layer (\code{bathy}) given the depth error model (\code{calc_depth_error});
#'   \item `depth_upr' is a number that defines the upper bound for the individual's possible depth on the bathymetry layer (\code{bathy}) given the depth error model (\code{calc_depth_error});
#'   \item `availability' is a logical variable that defines whether or not there were any cells on the bathymetry layer (\code{bathy}) within the boundaries for the individual's depth at that time step;
#' }
#'
#' }
#'
#' \subsection{(ii) record$spatial}{If \code{save_record_spatial > 0}, then for those time steps specified by \code{save_record_spatial}, each element in `record' contains a `spatial' list, with one element for each secondary time stamp, that records spatial information pertaining to the possible locations of the individual at that time step. For all algorithms, each `spatial' list includes the following elements:
#' \itemize{
#'   \item `map_timestep' is a RasterLayer of all the positions that the individual could have occupied at that time step, given the algorithm;
#'   \item `map_cumulative' is a RasterLayer of the cumulative of the number of times when the movement data were compatible with the individual being in that cell, under the specified algorithm, from all previous time steps up to the current time step (i.e., the sum of `map_timestep' across all time steps from the first time step to the current time step) (unless \code{normalise = TRUE} in which case the interpretation differs);
#' }
#' For the AC* algorithms, this list also includes information on the acoustic containers and location probability:
#' \itemize{
#'   \item `container_ap` is \code{\link[sp]{SpatialPolygonsDataFrame}} that defines the boundaries of the individual's location from the perspective of its previous location;
#'   \item `container_an' is a \code{\link[sp]{SpatialPolygonsDataFrame}} that defines the boundaries of the individual's location from the perspective of the receiver that records a detection at the moment of detection (i.e., the detection container);
#'   \item `container_b` is \code{\link[sp]{SpatialPolygonsDataFrame}} that defines the boundaries of the individual's location from (A) the perspective of
#'   \item `container_c` is a \code{\link[sp]{SpatialPolygonsDataFrame}} that defines the boundaries of the individual's location at a given time step, accounting for the information provided by previous, current and future locations;
#'   \item `kernel` is a \code{\link[raster]{raster}} that defines location probability across the grid, accounting for the receiver(s) at which an individual was detected and the receiver(s) at which an individual was not detected and receiver arrangement;
#' }
#'
#' For \code{save_record_spatial = 0L} or elements outside of the specified input to \code{save_record_spatial}, `spatial' is simply an empty list.
#' }
#'
#'
#' @author Edward Lavender
#' @docType package
#' @name acdc_record-class
NULL


######################################
######################################
#### acdc_archive-class

#' @title "acdc-archive" class
#' @description An S3 class that defines the object returned by an acoustic-container/depth-contour (AC/DC) algorithm (\code{\link[flapper]{ac}}, \code{\link[flapper]{dc}} or \code{\link[flapper]{acdc}}).
#'
#' @return An \code{\link[flapper]{acdc_archive-class}} object is a named list that contains the output of a call to \code{\link[flapper]{ac}}, \code{\link[flapper]{dc}} or \code{\link[flapper]{acdc}}. This contains the following elements: `archive', `ts_by_chunk', `time' and `args'.
#'
#' The main output is the `archive' element. This contains a list of \code{\link[flapper]{acdc_record-class}} objects returned by the call(s) to the internal routines that implement the algorithm. If the algorithm is implemented step-wise, this list contains a single element; if the algorithm is implemented chunk-wise, there is one element for each chunk. Each element contains an \code{\link[flapper]{acdc_record-class}} object with results of the call to internal routines for that chunk. The results across chunks can be aggregated using \code{\link[flapper]{acdc_simplify}}.
#'
#' The `ts_by_chunk' element is a list, with one element for each chunk, that contains the movement time series that were used in that chunk. Each element contains `acoustics' and `archival' dataframes, though the former is \code{NULL} for the DC algorithm and the latter is \code{NULL} for the AC algorithm. For the AC* algorithms, if there is more than one chunk, the last observation of each acoustic chunk is the same as the first acoustic observation for the next chunk. This results from splitting chunks at acoustic observations, which enables the results of chunks that are computed independently to be simply aggregated without the loss of any information.
#'
#' The `time' element is a dataframe that defines the times of sequential stages in the algorithm's progression, providing a record of computation time.
#'
#' The `args' element is a named list of user inputs that records the parameters used to generate the outputs.
#'
#' @author Edward Lavender
#' @docType package
#' @name acdc_archive-class
NULL
