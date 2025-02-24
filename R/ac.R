#' @title The acoustic-container (AC) algorithm
#' @description This function implements the acoustic-container (AC) algorithm. To implement the function, a dataframe (or list) of passive acoustic telemetry detections is required (\code{acoustics}). At each time step, the algorithm integrates information from past and future acoustic detections in the form of acoustic containers to determine the possible locations of an individual in an area (see Details).
#'
#' Under the default options, the approach is implemented step-wise (i.e., step-by-step across the whole time series). The result is a named list of outputs, including a record of the results for each time step, as well as a cumulative map of the expected proportion of time spent in each part of the study area across the whole time series. Alternatively, the approach can be implemented chunk-wise, in which case the acoustic time series is split into chunks (e.g., hourly, daily, monthly segments) and the algorithm is implemented within each chunk step-by-step. The main benefits of this approach are that it can be used to reconstruct putative patterns in space use over biologically meaningful periods separately and/or the chunk-wise implementation can be parallelised, improving computation time. (Chunk-wise results results are easily combined across the duration of the original time series without the loss of information via \code{\link[flapper]{acdc_simplify}}.). This option is implemented if (a) a list, rather than a dataframe, of acoustic detections is provided (via \code{acoustics}); (b) the user specifies that the time series should be split into chunks of a particular duration before the algorithm is initiated (via the \code{split} argument); and/or (c) the algorithm is implemented in parallel via \code{cl}, in which case the acoustic time series is split (if necessary) into user-defined or automatically defined chunks prior to computation. In this case, the result is a named list of outputs, as described above, but in which the results for each chunk are returned separately. If the chunks have been implemented simply to improve computation time via parallelisation, then the maps of space use for each chunk can be combined easily to generate a single, overall map of space use via \code{\link[flapper]{acdc_simplify}}.
#'
#' @param acoustics A dataframe, or a list of dataframes, that contains passive acoustic telemetry detection time series (see \code{\link[flapper]{dat_acoustics}} for an example) for a single individual. Each dataframe should contain the following columns: an integer vector of receiver IDs, named `receiver_id'; an integer vector of detection indices, named `index'; and a POSIXct vector of time stamps when detections were made, named `timestamp'. If a list of dataframes is supplied, dataframes must be refer to the detections of a single individual and be ordered by time (e.g., in hourly chunks). In addition, sequential list elements must be linked by identical receiver pairs (i.e., the final receiver at which the individual was detected for any given chunk must be the same as the receiver at which the individual was next detected at the start of the next chunk) because it is only in this specific scenario that information does not need to be shared across time steps (see \code{split}). The algorithm will be implemented on each dataframe, termed `chunk', either in sequence or parallel. Any empty or \code{NULL} elements will be removed automatically.
#' @param step A number that defines the time step length (s). The time series is rounded to the nearest \code{step} to ensure alignment (see Details). `Duplicate' detections (of the same individual at the same receiver in the same step) are dropped.
#' @param plot_ts A logical input that defines whether or not to the plot detection  time series before the algorithm is initiated.
#' @param bathy A \code{\link[raster]{raster}} that defines a grid across the area within which the individual could have moved. The coordinate reference system should be the Universal Transverse Mercator system, with distances in metres (see also \code{\link[flapper]{acs_setup_containers}}).
#' @param detection_containers A list of detection containers, with one element for each number from \code{1:max(acoustics$receiver_id)}, from \code{\link[flapper]{acs_setup_containers}}.
#' @param detection_kernels A named list of detection probability kernels, from \code{\link[flapper]{acs_setup_detection_kernels}} and created using consistent parameters as specified for other \code{acs_setup_*} functions and here (i.e., see the \code{overlaps}, \code{calc_detection_pr} and \code{map} arguments in \code{\link[flapper]{acs_setup_detection_kernels}}).
#' @param detection_kernels_overlap (optional) A named list, from \code{\link[flapper]{get_detection_containers_overlap}}, that defines, for each receiver, for each day over its deployment period, whether or not its detection container overlapped with those of other receivers. If \code{detection_kernels_overlap} and \code{detection_time_window} (below) are supplied, the implementation of detection probability kernels when a detection is made accounts for overlaps in receivers' detection containers; if unsupplied, receiver detection probability kernels are assumed not to overlap.
#' @param detection_time_window (optional) A number that defines the maximum duration (s) between consecutive detections at different receivers such that they can be said to have occurred at `effectively the same time'. This indicates that the same transmission was detected by multiple receivers. If \code{detection_kernels_overlap} (above) and \code{detection_time_window} are supplied, the implementation of detection probability kernels when a detection is made accounts for overlaps in receivers' detection containers, by up-weighting overlapping areas between receivers that detected the transmission and down-weighting overlapping areas between receivers that did not detect the transmission (see Details in \code{\link[flapper]{acs_setup_detection_kernels}}). Note that the timing of detections is affected by \code{step} (see Details).
#' @param mobility A number that defines the (Euclidean) distance (m) that an individual could move in the time steps between sequential detections (see also \code{\link[flapper]{acs_setup_containers}}).
#' @param normalise A logical variable that defines whether or not to normalise the map of possible locations at each time step. (The cumulative surface can be normalised via \code{\link[flapper]{acdc_simplify}}).
#' @param save_record_spatial An integer vector that defines the time steps for which to save a record of the spatial information from each time step. \code{save_record_spatial = 0L} suppresses the return of this information and \code{save_record_spatial = NULL} returns this information for all time steps. If the algorithm is applied chunk-wise, this spatial information must be returned for at least the first time step (the default) to aggregate maps across chunks (see \code{\link[flapper]{acdc_simplify}}). This information can also be used to plot time-specific results of the algorithm using \code{\link[flapper]{acdc_plot_trace}}, \code{\link[flapper]{acdc_plot_record}} and \code{\link[flapper]{acdc_animate_record}}.
#' @param write_record_spatial_for_pf (optional) A named list, passed to \code{\link[raster]{writeRaster}}, to save the \code{\link[raster]{raster}} of the individual's possible positions at each time step to file. The `filename' argument should be the directory in which to save files. Files are named by acoustic and intermediate (archival) time steps. For example, the file for the first acoustic time step and the first archival time step is named acc_1_arc_1.
#' @param save_args A logical input that defines whether or not to save the list of function inputs in the returned object.
#' @param verbose A logical variable that defines whether or not to print messages to the console or to file to relay function progress. If \code{con = ""}, messages are printed to the console (which is only supported if the algorithm is not implemented in parallel: see below); otherwise, they are written to file (see below).
#' @param con If \code{verbose = TRUE}, \code{con} is character string that defines how messages relaying function progress are returned. If \code{con = ""}, messages are printed to the console (unless redirected by \code{\link[base]{sink}}), an approach that is only implemented if the function is not implemented in parallel. Otherwise, \code{con} defines the directory into which to write .txt files, into which messages are written to relay function progress. This approach, rather than printing to the console, is recommended for clarity, speed and debugging. If the algorithm is implemented step-wise, then a single file is written to the specified directory named acdc_log.txt. If the algorithm is implemented chunk-wise, then an additional file is written for each chunk (named dot_acdc_log_1.txt, dot_acdc_log_2.txt and so on), with the details for each chunk.
#' @param progress (optional) If the algorithm is implemented step-wise, \code{progress} is an integer (\code{1}, \code{2} or \code{3}) that defines whether or not to display a progress bar in the console as the algorithm moves over acoustic time steps (\code{1}), the archival time steps between each pair of acoustic detections (\code{2}) or both acoustic and archival time steps (\code{3}), in which case the overall acoustic progress bar is punctuated by an archival progress bar for each pair of acoustic detections. This option is useful if there is a large number of archival observations between acoustic detections. Any other input will suppress the progress bar. If the algorithm is implemented for chunks, inputs to \code{progress} are ignored and a single progress bar is shown of the progress across acoustic chunks.
#' @param split A character string that defines the (approximate) time unit used to split acoustic time series into chunks (e.g., \code{"12 hours"}). If provided, this must be supported by \code{\link[base]{cut.POSIXt}} (otherwise, a pre-defined list of acoustic time series can be passed to \code{acoustics}, e.g., specifying seasonal chunks). If \code{split = NULL} and a cluster has been specified (see \code{cl}) (and \code{acoustics} is a dataframe), then the acoustic time series is automatically split into chunks and the algorithm implemented for each chunk in parallel. In all cases, splitting is subject to the constraint that  chunks must join at identical receiver pairs (i.e., the last receiver at which the individual was detected on one chunk must match the first receiver at which the individual was next detected at the start of the next chunk): in these specific scenarios, information does not need to transfer from one time step to the next.
#' @param cl,varlist (optional) Parallelisation options. \code{cl} is (a) a cluster object from \code{\link[parallel]{makeCluster}} or (b) an integer that defines the number of child processes to implement the algorithm in parallel. If supplied, the algorithm is implemented for each chunk in a list of acoustic time series, either (a) as supplied by the user (if \code{acoustics} is a list), (b) as defined by the input to \code{split}, or (c) as defined automatically from the number of nodes in the cluster if \code{split = NULL}. If \code{cl} is supplied, \code{varlist} may also be required. This is a character vector of objects to export (see \code{\link[flapper]{cl_export}}). Exported variables must be located in the global environment. If a cluster is supplied, the connection to the cluster is closed within the function (see \code{\link[flapper]{cl_stop}}). For further information, see \code{\link[flapper]{cl_lapply}} and \code{\link[flapper]{flapper-tips-parallel}}.
#'
#' @details The acoustic-container (AC) algorithm is an approach which uses acoustic detections to infer the possible locations of tagged animals within an area over some time interval. The locational information provided by acoustic detections is represented by acoustic containers, which are areas around receivers that define where an individual could have been at each time point given the spatiotemporal pattern of detections at receivers, a model of detection probability and a movement parameter.
#'
#' In outline, the crux of the approach is the recognition that acoustic detections typically occur irregularly, but we can consider a sequence of regular time steps between any pair of detections. Each detection anchors our knowledge of the location of an individual around a particular receiver (assuming that all detections are true detections). As time passes between pairs of detections, our uncertainty in the geographical location of an individual expands around the receiver at which the individual was detected before shrinking towards the receiver at which it was next detected. The dynamics of this process are captured by the expansion, contraction and intersection of pairs of acoustic containers.
#'
#' More specifically, at each time step, we can consider the set of possible locations for the individual from perspective of (a) the receiver at which the individual was detected and (b) the receiver at which the individual is next detected. When an individual is detected, it must be within some radius---say 800 m---of that receiver termed the `detection container'. From the perspective of the receiver at which the individual is next detected, the set of possible locations of the individual is wider in line with the time between detections and the movement speed of the animal. The intersection of these two areas defines the set of possible locations for the individual. (In most cases, this is simply that defined by the detection container around the first receiver.) With a more-refined model of detection probability, it may be possible to predict more precisely where the individual is likely to have been within this area. (In situations with depth data, the ACDC algorithm further restricts the set of locations by incorporating information on the animal's depth via depth contours.) Moving forward in time, we can consider a number of regular time steps before the next detection. During this time, from the perspective of the receiver at which the individual was detected, the set of possible locations of the individual expands, because it could have moved further away from the receiver; meanwhile, from the perspective of the receiver at which the individual was next detected, the set of possible locations of the individual shrinks, as the individual must have been located within the detection container of that receiver by the time of the detection. This process is described by the expansion and contraction of `acoustic containers'. At each time step, the intersection of the two containers defines the set of possible locations of the individual, possibly weighted by a detection probability (given the lack of detections during this time). These dynamics recognise that as time passes between detections the individual could have moved away from the receiver at which it was last detected but only at a rate and in a direction that fits with the receiver at which the individual was next detected. Thus, when the individual is detected again, our uncertainty about where it could have been collapses to the detection container around the next receiver (and its intersection with the individual's previous (expanded) location and the acoustic container around the following receiver), possibly weighted by a model of detection probability. Throughout this process, the rate of change in container size depends a movement parameter that describes the maximum swimming speed.
#'
#' This discussion assumes that the timing of detections and intermediate time steps between detections are perfectly aligned. In reality, there is likely to be a mismatch between the timing of detections and the intermediate time steps between detections, which may be poorly approximated the constant expansion and contraction of acoustic containers. The simplest solution to this issue is to round the acoustic time steps to the resolution of the intermediate time steps (specified by \code{step}). This results in a small loss of precision, but assuming that \code{step} is relatively small, the effect should be negligible. In any case, clocks on different receivers are unlikely to be perfectly synced throughout a study. This solution is also computationally preferable to an alternative approach in which the expansion and contraction of containers varies through time, depending on the gaps between observations and other time-specific variables such as behavioural state. However, the latter approach may be implemented in due course.
#'
#' The end result is a map that shows the expected time spent in different parts of a study area. The main limitation of this approach is the simple treatment of movement, but particle filtering can be used to extent this approach via the incorporate a movement model (see \code{\link[flapper]{pf}}).
#'
#' @return The function returns an \code{\link[flapper]{acdc_archive-class}} object. If a connection to write files has also been specified, an overall log (acdc_log.txt) as well as chunk-specific logs from calls to \code{\link[flapper]{.acs}}, if applicable, are written to file.
#'
#' @seealso This function calls \code{\link[flapper]{.acs_pl}} and \code{\link[flapper]{.acs}} to implement the AC algorithm. \code{\link[flapper]{acs_setup_containers}} defines the detection containers required by this function. \code{\link[flapper]{acs_setup_mobility}} is used to examine the assumption of the constant `mobility' parameter. \code{\link[flapper]{acs_setup_detection_kernels}} produces detection probability kernels for incorporation into the function. \code{\link[flapper]{acdc_simplify}} simplifies the outputs and \code{\link[flapper]{acdc_plot_trace}}, \code{\link[flapper]{acdc_plot_record}} and \code{\link[flapper]{acdc_animate_record}} visualise the results. The AC algorithm can be extended to incorporate depth contours via \code{\link[flapper]{acdc}}. Particle filtering can be used to reconstruct movement path via \code{\link[flapper]{pf}}.
#'
#' @examples
#' #### Step (1) Implement setup_acdc_*() steps
#' # ... Define detection containers required for AC algorithm (see setup_acdc_containers())
#'
#' #### Step (2) Prepare movement time series for algorithm
#' # Focus on an example individual for speed
#' id <- 25
#' acc <- dat_acoustics[dat_acoustics$individual_id == id, ][1:25, ]
#' # Process time series (if necessary)
#' # ... Observations should be processed to the nearest time step
#' # ... Duplicate detections should be dropped
#'
#' #### Example (1) Implement AC algorithm with default arguments
#' # This implements the algorithm on a single core, printing messages
#' # ... to the console to monitor function progress.
#' out_ac <- ac(
#'   acoustics = acc,
#'   step = 120,
#'   bathy = dat_gebco,
#'   detection_containers = dat_containers,
#'   mobility = 200
#' )
#'
#' #### Subsequent examples follow the implementation given in acdc().
#'
#' @author Edward Lavender
#' @export

ac <- function(acoustics,
               step,
               bathy,
               plot_ts = TRUE,
               detection_containers,
               detection_kernels = NULL,
               detection_kernels_overlap = NULL,
               detection_time_window = 5,
               mobility,
               normalise = TRUE,
               save_record_spatial = 1L,
               write_record_spatial_for_pf = NULL,
               save_args = TRUE,
               verbose = TRUE,
               con = "",
               progress = 1L,
               split = NULL,
               cl = NULL,
               varlist = NULL) {
  # Initiate function
  t_onset <- Sys.time()
  message(paste0("flapper::ac() called (@ ", t_onset, ")..."))
  # Check for missing inputs
  for (arg in c(acoustics, step, bathy, detection_containers, mobility)) 1L
  # Pass inputs to .acs_pl()
  out <-
    .acs_pl(
      acoustics = acoustics,
      archival = NULL,
      step = step,
      plot_ts = plot_ts,
      bathy = bathy,
      detection_containers = detection_containers,
      detection_kernels = detection_kernels,
      detection_kernels_overlap = detection_kernels_overlap,
      detection_time_window = detection_time_window,
      mobility = mobility,
      calc_depth_error = NULL,
      normalise = normalise,
      save_record_spatial = save_record_spatial,
      write_record_spatial_for_pf = write_record_spatial_for_pf,
      verbose = verbose,
      con = con,
      progress = progress,
      split = split,
      cl = cl,
      varlist = varlist
    )
  if (save_args) {
    out$args <- list(
      acoustics = acoustics,
      step = step,
      bathy = bathy,
      plot_ts = plot_ts,
      detection_containers = detection_containers,
      detection_kernels = detection_kernels,
      detection_kernels_overlap = detection_kernels_overlap,
      detection_time_window = detection_time_window,
      mobility = mobility,
      normalise = normalise,
      save_record_spatial = save_record_spatial,
      write_record_spatial_for_pf = write_record_spatial_for_pf,
      save_args = save_args,
      verbose = verbose,
      con = con,
      progress = progress,
      split = split,
      cl = cl,
      varlist = varlist
    )
  }
  t_end <- Sys.time()
  message(paste0("flapper::ac() finished (@ ", t_end, ")..."))
  return(out)
}
