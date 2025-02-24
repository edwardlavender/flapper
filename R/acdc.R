######################################
######################################
#### acdc()

#' @title The acoustic-container depth-contour (ACDC) algorithm
#' @description This function implements the acoustic-container depth-contour (ACDC) algorithm. This is an extension of the AC algorithm implemented by \code{\link[flapper]{ac}} that that integrates acoustic detections and depth observations to the expected time spent in different parts of a study area over the period of observations. To implement the function, a dataframe (or list) of passive acoustic telemetry detections is required (\code{acoustics}), alongside a dataframe of depth observations (\code{archival}). At each time step, the algorithm integrates information from past and future acoustic detections in the form of acoustic containers and information from depth observations in the form of depth contours to determine the possible locations of an individual in an area (see Details). As for the AC algorithm (\code{\link[flapper]{ac}}), the function can be implemented step-wise or chunk-wise and the outputs can be processed via \code{\link[flapper]{acdc_simplify}}.
#'
#' @param acoustics A dataframe, or a list of dataframes, that contains passive acoustic telemetry detection time series (see \code{\link[flapper]{ac}}).
#' @param archival A dataframe that contains depth time series (see \code{\link[flapper]{dat_archival}} for an example) for a specific individual. This should contain the following columns: a numeric vector of observed depths, named `depth'; and a POSIXct vector of time stamps when observations were made, named `timestamp'. Depths should be recorded in the same units and with the same sign as the bathymetry data (see \code{bathy}). Absolute depths (m) are suggested. Unlike the detection time series, archival time stamps are assumed to have occurred at regular intervals. Both acoustic and archival time series are rounded to the resolution of \code{archival} observations to ensure alignment (see Details). `Duplicate' observations (of the same individual at the same receiver, if applicable, in the same time step) are dropped.
#' @param plot_ts A logical input that defines whether or not to plot the detection and depth time series before the algorithm is initiated. This provides a useful visualisation of the extent to which they overlap.
#' @param bathy A \code{\link[raster]{raster}} that defines the bathymetry across the area within which the individual could have moved. This must be recorded in the same units and with the same sign as the depth observations (see \code{archival}). The coordinate reference system should be the Universal Transverse Mercator system, with distances in metres (see also \code{\link[flapper]{acs_setup_containers}}).
#' @param detection_containers A list of detection containers, with one element for each number from \code{1:max(acoustics$receiver_id)}, from \code{\link[flapper]{acs_setup_containers}}.
#' @param detection_kernels A named list of detection probability kernels (see \code{\link[flapper]{ac}}).
#' @param detection_kernels_overlap (optional) A named list of detection probability kernel overlaps (see \code{\link[flapper]{ac}}).
#' @param detection_time_window A number that defines the detection time window (see \code{\link[flapper]{ac}}).
#' @param mobility The mobility parameter (see \code{\link[flapper]{ac}}).
#' @param calc_depth_error A function that returns the depth errors around a vector of depths. The function should accept vector of depths (from \code{archival$depth}) and return a matrix, with one row for each (lower and upper) error and one one column for each depth (if the error varies with depth). For each depth, the two numbers are added to the observed depth to define the range of depths on the bathymetry raster (\code{bathy}) that the individual could plausibly have occupied at any time. Since the depth errors are added to the individual's depth, the first number should be negative (i.e., the individual could have been shallower that observed) and the second positive (i.e., the individual could have been deeper than observed). The appropriate form for \code{calc_depth_error} depends on the species (pelagic versus demersal/benthic species), the measurement error for the depth observations in \code{archival} and bathymetry (\code{bathy}) data, as well as the tidal range (m) across the area (over the duration of observations). For example, for a pelagic species, the constant function \code{calc_depth_error = function(...) matrix(c(-2.5, Inf)} implies that the individual could have occupied bathymetric cells that are deeper than the observed depth + (-2.5) m and shallower than Inf m (i.e. the individual could have been in any location in which the depth was deeper than the shallow depth limit for the individual). In contrast, for a benthic species, the constant function \code{calc_depth_error = function(...) matrix(c(-2.5, 2.5), nrow = 2)} implies that the individual could have occupied bathymetric cells whose depth lies within the interval defined by the observed depth + (-2.5) and + (+2.5) m.
#'
#' @param normalise A logical input that defines whether or not to normalise maps (see \code{\link[flapper]{ac}}).
#' @param save_record_spatial An integer of the spatial layers to save (see \code{\link[flapper]{ac}}).

#' @param write_record_spatial_for_pf A named list used to write time step-specific maps to file (see \code{\link[flapper]{ac}}).
#' @param save_args A logical input that defines whether or not to save the list of function inputs (see \code{\link[flapper]{ac}}).
#' @param verbose A logical variable that defines whether or not to print messages to the console or to file to relay function progress (see \code{\link[flapper]{ac}}).
#' @param con If \code{verbose = TRUE}, \code{con} is character string defines how messages relaying function progress are returned (see \code{\link[flapper]{ac}}).
#' @param progress An integer controlling the progress bar (see \code{\link[flapper]{ac}}).
#' @param split A character string that defines the time unit used to split acoustic time series into chunks (see \code{\link[flapper]{ac}}).
#' @param cl,varlist (optional) Parallelisation options (see \code{\link[flapper]{ac}}).
#'
#' @details The acoustic-container depth-contour (ACDC) algorithm is an approach which integrates acoustic detections and depth observations to infer the possible locations of tagged animals within an area over some time interval. The locational information provided by acoustic detections is represented by acoustic containers, which are areas around receivers that define where an individual could have been at each time point given the spatiotemporal pattern of detections at receivers, a model of detection probability and a movement parameter (see \code{\link[flapper]{ac}}). The locational information provided by depth observations is represented by depth contours, which are areas that define where an individual could have been at each time point given its depth and the local bathymetry (see \code{\link[flapper]{dc}}).
#'
#' In outline, the crux of the approach is the recognition that acoustic detections typically occur irregularly, while archival observations occur at regular intervals. Each detection anchors our knowledge of the location of an individual around a particular receiver (assuming that all detections are true detections). As time passes between acoustic detections, our uncertainty about the geographical location of an individual expands around the receiver at which it was detected before shrinking towards the receiver at which it was next detected. (This is the AC algorithm.) During this time, regular depth observations restrict the number of possible locations in which the individual could have been located at each time step, assuming that the bathymetric landscape is non-homogeneous. This is true both for pelagic species, which must be in an area where the depth of the seabed is at least as deep as the observed depth, and for benthic/demersal species, which must be in an area where the depth of the seabed is relatively close to the observed depth. (This is the DC algorithm.) In reality, the timing of acoustic and archival observations may not align perfectly; this is currently resolved by rounding the acoustic time series to the temporal resolution of depth observations (see also \code{\link[flapper]{ac}}). In this way, the ACDC algorithm moves over acoustic and archival time steps, summing the positions in which individuals could have been located based on these two sources of data, to generate a map of that shows expected time spent in different parts of the study area over the time interval under consideration.

#' @return The function returns an \code{\link[flapper]{acdc_archive-class}} object. If a connection to write files has also been specified, an overall log (acdc_log.txt) as well as chunk-specific logs from calls to \code{\link[flapper]{.acs}}, if applicable, are written to file.
#'
#' @seealso This function calls \code{\link[flapper]{.acs_pl}} and \code{\link[flapper]{.acs}} to implement the ACDC algorithm. The AC component can be implemented via  \code{\link[flapper]{ac}} and the DC component via \code{\link[flapper]{dc}}. \code{\link[flapper]{acs_setup_containers}} defines the detection containers required by this function.  \code{\link[flapper]{acs_setup_mobility}} is used to examine the assumption of the constant `mobility' parameter. \code{\link[flapper]{acs_setup_detection_kernels}} produces detection probability kernels for incorporation into the function. \code{\link[flapper]{acdc_simplify}} simplifies the outputs and \code{\link[flapper]{acdc_plot_trace}}, \code{\link[flapper]{acdc_plot_record}} and \code{\link[flapper]{acdc_animate_record}} visualise the results. Particle filtering can be used to reconstruct movement paths.
#'
#' @examples
#' #### Step (1) Implement setup_acdc_*() steps
#' # ... Define detection containers required for ACDC algorithm (see setup_acdc_containers())
#'
#' #### Step (2) Prepare movement time series for algorithm
#' # Focus on an example individual
#' id <- 25
#' acc <- dat_acoustics[dat_acoustics$individual_id == id, ]
#' arc <- dat_archival[dat_archival$individual_id == id, ]
#' # Focus on the subset of data for which we have both acoustic and archival detections
#' acc <- acc[acc$timestamp >= min(arc$timestamp) - 2 * 60 &
#'   acc$timestamp <= max(arc$timestamp) + 2 * 60, ]
#' arc <- arc[arc$timestamp >= min(acc$timestamp) - 2 * 60 &
#'   arc$timestamp <= max(acc$timestamp) + 2 * 60, ]
#' # We'll focus on a one day period with overlapping detection/depth time series for speed
#' end <- as.POSIXct("2016-03-18")
#' acc <- acc[acc$timestamp <= end, ]
#' arc <- arc[arc$timestamp <= end, ]
#' arc <- arc[arc$timestamp >= min(acc$timestamp) - 2 * 60 &
#'   arc$timestamp <= max(acc$timestamp) + 2 * 60, ]
#' # Process time series (if necessary)
#' # ... Observations should be processed to the nearest time step
#' # ... Duplicate detections should be dropped
#'
#' #### Example (1) Implement ACDC algorithm with default arguments
#' # This implements the algorithm on a single core, printing messages
#' # ... to the console to monitor function progress.
#' out_acdc <- acdc(
#'   acoustics = acc,
#'   archival = arc,
#'   bathy = dat_gebco,
#'   detection_containers = dat_containers,
#'   mobility = 200,
#'   calc_depth_error = function(...) matrix(c(-2.5, 2.5), nrow = 2)
#' )
#' # The function returns a list with four elements
#' # ... archive contains the results of the algorithm, implemented by the back-end
#' # ... function .acs(). The other elements provide the time series
#' # ... for each chunk, the time of the algorithm and a list of user inputs
#' summary(out_acdc)
#'
#' #### Example (2): Write messages to file to monitor function progress via 'con'
#' out_acdc <- acdc(
#'   acoustics = acc,
#'   archival = arc,
#'   bathy = dat_gebco,
#'   detection_containers = dat_containers,
#'   mobility = 200,
#'   calc_depth_error = function(...) matrix(c(-2.5, 2.5), nrow = 2),
#'   con = tempdir()
#' )
#' acdc_log <- readLines(paste0(tempdir(), "/acdc_log.txt"))
#' utils::head(acdc_log, 10)
#' file.remove(paste0(tempdir(), "/acdc_log.txt"))
#'
#' #### Example (3): Implement the algorithm and return spatial information
#' # Specify save_record_spatial = NULL to include spatial information for all time steps
#' # ... or a vector to include this information for specific time steps
#' out_acdc <- acdc(
#'   acoustics = acc,
#'   archival = arc,
#'   bathy = dat_gebco,
#'   detection_containers = dat_containers,
#'   mobility = 200,
#'   calc_depth_error = function(...) matrix(c(-2.5, 2.5), nrow = 2),
#'   save_record_spatial = NULL
#' )
#'
#' #### Example (4): Implement the algorithm in parallel e.g. by supplying a cluster
#' # If verbose = TRUE (the default), it is necessary to specify a directory
#' # ... into which dot_acdc_log_*.txt files are saved (i.e., messages
#' # ... cannot be written to the console in parallel)
#' out_acdc <- acdc(
#'   acoustics = acc,
#'   archival = arc,
#'   bathy = dat_gebco,
#'   detection_containers = dat_containers,
#'   mobility = 200,
#'   calc_depth_error = function(...) matrix(c(-2.5, 2.5), nrow = 2),
#'   con = tempdir(),
#'   cl = parallel::makeCluster(2L)
#' )
#' ## Check logs
#' list.files(tempdir())
#' # "acdc_log.txt" contains the log for the overall function
#' acdc_log <- readLines(paste0(tempdir(), "/acdc_log.txt"))
#' head(acdc_log, 20)
#' # "acdc_log_1.txt", "acdc_log_2.txt" etc contain chunk-specific logs
#' acdc_log_1 <- readLines(paste0(tempdir(), "/dot_acdc_log_1.txt"))
#' utils::head(acdc_log_1)
#' utils::tail(acdc_log_1)
#' ## Examine outputs
#' # Note that there are now four elements in archive, one for each chunk
#' # Likewise, there are four elements in ts_by_chunk,
#' # ... containing the movement time series for each chunk.
#' summary(out_acdc)
#' # Note that the last observation of each time series overlaps with the
#' # ... first observation for the next chunk, to prevent loss of information
#' lapply(
#'   out_acdc$ts_by_chunk,
#'   function(chunk) chunk$acoustics[c(1, nrow(chunk$acoustics)), ]
#' )
#'
#' #### Example (5) Biologically meaningful chunks can be specified via
#' # .. the 'split' argument or by passing a list of acoustic time series
#' # .. already split by list of dataframes to 'acoustics'
#' ## Using the split argument:
#' out_acdc <- acdc(
#'   acoustics = acc,
#'   archival = arc,
#'   bathy = dat_gebco,
#'   detection_containers = dat_containers,
#'   mobility = 200,
#'   calc_depth_error = function(...) matrix(c(-2.5, 2.5), nrow = 2),
#'   con = tempdir(),
#'   cl = parallel::makeCluster(2L),
#'   split = "2 hours"
#' )
#' ## Passing a list of dataframes
#' # This option can provide more flexibility than split.
#' # ... This example could also be used using split, as described above,
#' # ... but this is not the case for all time categories (e.g., seasons).
#' acc$chunk <- cut(acc$timestamp, "2 hours")
#' acc_ls <- split(acc, acc$chunk)
#' out_acdc <- acdc(
#'   acoustics = acc_ls,
#'   archival = arc,
#'   bathy = dat_gebco,
#'   detection_containers = dat_containers,
#'   mobility = 200,
#'   calc_depth_error = function(...) matrix(c(-2.5, 2.5), nrow = 2),
#'   con = tempdir(),
#'   cl = parallel::makeCluster(2L)
#' )
#'
#' #### Example (5) Implement the algorithm for multiple individuals
#' # ... To do this, it is necessary to apply the function iteratively
#' # ... to each individual.
#' # Pre-processing to define computation time
#' # ... E.g., careful definition of time series
#' if (flapper_run_slow) {
#'   # Define cluster
#'   cluster <- FALSE
#'   if (cluster) {
#'     cl <- parallel::makeCluster(2L)
#'     parallel::clusterExport(cl = cl, varlist = c(
#'       "acdc",
#'       "dat_archival",
#'       "dat_gebco",
#'       "dat_containers"
#'     ))
#'   } else {
#'     cl <- NULL
#'   }
#'   # Implement algorithm for each individual
#'   acdc_out_by_id <-
#'     pbapply::pblapply(split(dat_acoustics, dat_acoustics$individual_id), cl = cl, function(acc) {
#'       # Define individual-specific folder in which to save function logs
#'       dir_global <- paste0(tempdir(), "/")
#'       dir_id <- paste0(dir_global, acc$individual_id[1], "/")
#'       if (!dir.exists(dir_id)) dir.create(dir_id)
#'       # Focus on a small sample of time series for speed
#'       acc <- acc[1:3, ]
#'       # Isolate archival data for individual
#'       arc <- dat_archival[dat_archival$individual_id == acc$id[1], ]
#'       # Implement algorithm
#'       acdc_out <- acdc(
#'         acoustics = acc,
#'         archival = dat_archival,
#'         bathy = dat_gebco,
#'         detection_containers = dat_containers,
#'         mobility = 200,
#'         calc_depth_error = function(...) matrix(c(-2.5, 2.5), nrow = 2),
#'         save_record_spatial = 1:10L,
#'         con = dir_id
#'       )
#'       # Include logs in output
#'       acdc_log <- lapply(list.files(dir_id, full.names = TRUE), readLines)
#'       # Simplify the results at this stage or outside of this loop
#'       # ... using acdc_simplify()
#'       # Return results for specified individual
#'       out <- list(acdc_id = acc$individual_id[1], acdc_out = acdc_out, acdc_log = acdc_log)
#'       return(out)
#'     })
#'   if (!is.null(cl)) parallel::stopCluster(cl)
#'   summary(acdc_out_by_id)
#' }
#'
#' #### Step (3) Simplify the function outputs
#' # This step aggregates information across chunks, which is necessary to
#' # ... plot information aggregated across all chunks (see below).
#'
#' #### Step (4) Examine function outputs, e.g., via plotting
#' # See acdc_plot_record() and acdc_animate_record() to visualise the results
#' # ... (either for a specific chunk or aggregated across all chunks
#' # ... using acdc_simplify() as described above).
#'
#' @author Edward Lavender
#' @export
#'

acdc <- function(acoustics,
                 archival,
                 plot_ts = TRUE,
                 bathy,
                 detection_containers,
                 detection_kernels = NULL, detection_kernels_overlap = NULL, detection_time_window = 5,
                 mobility,
                 calc_depth_error = function(...) matrix(c(-2.5, 2.5), nrow = 2),
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
  message(paste0("flapper::acdc() called (@ ", t_onset, ")..."))
  # Check for missing inputs
  for (arg in c(acoustics, archival, bathy, detection_containers, mobility)) 1L
  # Check archival names
  check_names(input = archival, req = "timestamp")
  out <-
    .acs_pl(
      acoustics = acoustics,
      archival = archival,
      step = as.numeric(difftime(archival$timestamp[2], archival$timestamp[1], units = "s")),
      plot_ts = plot_ts,
      bathy = bathy,
      detection_containers = detection_containers,
      detection_kernels = detection_kernels,
      detection_kernels_overlap = detection_kernels_overlap,
      detection_time_window = detection_time_window,
      mobility = mobility,
      calc_depth_error = calc_depth_error,
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
      archival = archival,
      plot_ts = plot_ts,
      bathy = bathy,
      detection_containers = detection_containers,
      detection_kernels = detection_kernels,
      detection_kernels_overlap = detection_kernels_overlap,
      detection_time_window = detection_time_window,
      mobility = mobility,
      calc_depth_error = calc_depth_error,
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
  }
  t_end <- Sys.time()
  message(paste0("flapper::acdc() finished (@ ", t_end, ")..."))
  return(out)
}
