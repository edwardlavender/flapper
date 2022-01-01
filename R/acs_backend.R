######################################
######################################
#### .acs()

#' @title Back-end implementation of the AC and ACDC algorithms
#' @description This function is the back-end of the acoustic-centroid (AC) and acoustic-centroid depth-contour (ACDC) algorithms.
#' @param acoustics A dataframe that contains passive acoustic telemetry detection time series for a specific individual (see \code{\link[flapper]{dat_acoustics}} for an example). This should contain the following columns: an integer vector of receiver IDs, named `receiver_id' (that must match that inputted to \code{\link[flapper]{acs_setup_centroids}}); an integer vector of detection indices, named `index'; a POSIXct vector of time stamps when detections were made, named `timestamp'; and a numeric vector of those time stamps, named `timestamp_num'.
#' @param archival For the ACDC algorithm, \code{archival} is a dataframe that contains depth time series for the same individual (see \code{\link[flapper]{dat_archival}} for an example). This should contain the following columns: a numeric vector of observed depths, named `depth'; a POSIXct vector of time stamps when observations were made, named `timestamp'; and a numeric vector of those time stamps, named `timestamp_num'. Depths should be recorded in the same units and with the same sign as the bathymetry data (see \code{bathy}). Absolute depths (m) are suggested. Unlike the detection time series, archival time stamps are assumed to have occurred at regular intervals.
#' @param step A number that defines the time step length (s) between consecutive detections. If \code{archival} is supplied, this is the resolution of the archival data (e.g., 120 s).
#' @param round_ts A logical input that defines whether or not to round the \code{acoustics} and/or \code{archival} time series to the nearest \code{step/60} minutes. This step is necessary to ensure appropriate alignment between the two time series given the assumption of a constant \code{mobility} parameter (see below). For implementations via \code{\link[flapper]{.acs_pl}}, this step is implemented in \code{\link[flapper]{.acs_pl}}.
#' @param plot_ts A logical input that defines whether or not to the plot movement series before the algorithm is initiated.
#' @param bathy A \code{\link[raster]{raster}} that defines the area (for the AC algorithm) or bathymetry (for the ACDC* algorithm) across the area within which the individual could have moved. For the ACDC algorithm, this must be recorded in the same units and with the same sign as the depth observations (see \code{archival}). The coordinate reference system should be the Universal Transverse Mercator system, with distances in metres (see also \code{\link[flapper]{acs_setup_centroids}}).
#' @param map (optional) A blank \code{\link[raster]{raster}}, with the same properties (i.e., dimensions, resolution, extent and coordinate reference system) as the area/bathymetry raster (see \code{bathy}), but in which all values are 0. If \code{NULL}, this is computed internally, but supplying a pre-defined raster can be more computationally efficient if the function is applied iteratively (e.g., over different time windows).
#' @param detection_kernels A named list of detection probability kernels, from \code{\link[flapper]{acs_setup_detection_kernels}} and created using consistent parameters as specified for other \code{acs_setup_*} functions and here (i.e., see the \code{overlaps}, \code{calc_detection_pr} and \code{map} arguments in \code{\link[flapper]{acs_setup_detection_kernels}}).
#' @param detection_kernels_overlap (optional) A named list (the `list_by_receiver' element from \code{\link[flapper]{get_detection_centroids_overlap}}), that defines, for each receiver, for each day over its deployment period, whether or not its detection centroid overlapped with those of other receivers. If \code{detection_kernels_overlap} and \code{detection_time_window} (below) are supplied, the implementation of detection probability kernels when a detection is made accounts for overlaps in receivers' detection centroids; if unsupplied, receiver detection probability kernels are assumed not to overlap.
#' @param detection_time_window (optional) A number that defines the maximum duration (s) between consecutive detections at different receivers such that they can be said to have occurred at `effectively the same time'. This indicates that the same transmission was detected by multiple receivers. If \code{detection_kernels_overlap} (above) and \code{detection_time_window} are supplied, the implementation of detection probability kernels when a detection is made accounts for overlaps in receivers' detection centroids, by up-weighting overlapping areas between receivers that detected the transmission and down-weighting overlapping areas between receivers that did not detect the transmission (see Details in \code{\link[flapper]{acs_setup_detection_kernels}}). Note that the timing of detections is affected by \code{step} (see Details).
#' @param detection_centroids A list of detection centroids, with one element for each number from \code{1:max(acoustics$receiver_id)}, from \code{\link[flapper]{acs_setup_centroids}}.
#' @param mobility A number that defines the distance (m) that an individual could move in the time steps between acoustic detections (see also \code{\link[flapper]{acs_setup_centroids}}).
#' @param calc_depth_error In the ACDC algorithm, \code{calc_depth_error} is function that returns the depth errors around a vector of depths. The function should accept vector of depths (from \code{archival$depth}) and return a matrix, with one row for each (lower and upper) error and one one column for each depth (if the error varies with depth). For each depth, the two numbers are added to the observed depth to define the range of depths on the bathymetry raster (\code{bathy}) that the individual could plausibly have occupied at any time. Since the depth errors are added to the individual's depth, the first number should be negative (i.e., the individual could have been shallower that observed) and the second positive (i.e., the individual could have been deeper than observed). The appropriate form for \code{calc_depth_error} depends on the species (pelagic versus demersal/benthic species), the measurement error for the depth observations in \code{archival} and bathymetry (\code{bathy}) data, as well as the tidal range (m) across the area (over the duration of observations). For example, for a pelagic species, the constant function \code{calc_depth_error = function(...) matrix(c(-2.5, Inf)} implies that the individual could have occupied bathymetric cells that are deeper than the observed depth + (-2.5) m and shallower than Inf m (i.e. the individual could have been in any location in which the depth was deeper than the shallow depth limit for the individual). In contrast, for a benthic species, the constant function \code{calc_depth_error = function(...) matrix(c(-2.5, 2.5), nrow = 2)} implies that the individual could have occupied bathymetric cells whose depth lies within the interval defined by the observed depth + (-2.5) and + (+2.5) m.
#' @param normalise A logical variable that defines whether or not to normalise the map of possible locations at each time step. In both cases, at each time step the possible locations of the individual are scaled so that the most probable locations have a score of 1 and other scores vary between 0--1. If \code{normalise = FALSE}, these scores are simply summed at each time step, in which case scores on the final map can be interpreted as the number of time steps when the individual could have been in any given location. In contrast, if \code{normalise = TRUE}, at each time step scores are normalised so that they sum to one; the consequence is that time steps with detections, when uncertainty in the individual's location concentrates in the detection centroid around a receiver, are weighted more strongly than time steps between detections, when the uncertainty in the individual's location is spread across a larger area. The final surface can be normalised within \code{\link[flapper]{acdc_simplify}}, with in each cell (0--1) providing a measure of the relative potential use of each location.
#' @param chunk An integer that defines the chunk ID (from \code{\link[flapper]{.acs_pl}}).
#' @param save_record_spatial An integer vector that defines the time steps for which to save a record of the spatial information from each time step. \code{save_record_spatial = 0} suppresses the return of this information and \code{save_record_spatial = NULL} returns this information for all time steps.
#' @param write_record_spatial_for_pf (optional) A named list, passed to \code{\link[raster]{writeRaster}}, to save a subset of the spatial record (specifically the \code{\link[raster]{raster}} of the individual's possible positions at each time step) to file. This forms the basis for extending maps of space use via particle filtering (see \code{\link[flapper]{pf}}.) The `filename' argument should be the directory in which to save files. Files are named by chunk ID, acoustic and internal (archival) time steps. For example, the file for the first chunk, the first acoustic time step and the first archival time step is named `chu_1_acc_1_arc_1'.
#' @param save_args A logical input that defines whether or not to include a list of function arguments in the outputs. This can be switched off if the function is applied iteratively.
#' @param verbose A logical variable that defines whether or not to print messages to the console or to file to relay function progress. If \code{con = ""}, messages are printed to the console; otherwise, they are written to file (see below).
#' @param con If \code{verbose = TRUE}, \code{con} is character string that defines the full pathway to a .txt file (which can be created on-the-fly) into which messages are written to relay function progress. This approach, rather than printing to the console, is recommended for clarity, speed and debugging.
#' @param progress (optional) If the algorithm is implemented step-wise, \code{progress} is an integer (\code{1}, \code{2} or \code{3}) that defines whether or not to display a progress bar in the console as the algorithm moves over acoustic time steps (\code{1}), the `archival' time steps between each pair of acoustic detections (\code{2}) or both acoustic and archival time steps (\code{3}), in which case the overall acoustic progress bar is punctuated by an archival progress bar for each pair of acoustic detections. This option is useful if there is a large number of archival observations between acoustic detections. Any other input will suppress the progress bar. If the algorithm is implemented for chunks, inputs to \code{progress} are ignored and a single progress bar is shown of the progress across acoustic chunks.

#' @param check A logical input that defines whether or not to check function inputs. This can be switched off to improve computation time when the function is applied iteratively or via a front-end function (e.g., \code{\link[flapper]{ac}} or \code{\link[flapper]{acdc}}).
#'
#' @return The function returns an \code{\link[flapper]{acdc_record-class}} object with the following elements: `map', `record', `time', `args', `chunks' and `simplify'. The main output element is the `map' RasterLayer that shows where the individual could have spent more or less time over the duration of the movement time series. The `record' element records time-specific information on the possible locations of the individual, and can be used to plot maps of specific time points or to produce animations (for the time steps specified by \code{save_record_spatial}). The `time' element is a dataframe that defines the times of sequential stages in the algorithm's progression, providing a record of computation time. The `args' element is a named list of user inputs that record the parameters used to generate the outputs (if \code{save_args = TRUE}, otherwise the `args' element is \code{NULL}).
#'
#' @seealso The front-end functions \code{\link[flapper]{ac}} and \code{\link[flapper]{acdc}} call \code{\link[flapper]{.acs_pl}} which in turn calls this function. \code{\link[flapper]{acs_setup_centroids}} defines the detection centroids required by this function.  \code{\link[flapper]{acs_setup_mobility}} is used to examine the assumption of the constant `mobility' parameter. \code{\link[flapper]{acs_setup_detection_kernels}} produces detection probability kernels for incorporation into the function. For calls via \code{\link[flapper]{ac}} and \code{\link[flapper]{acdc}}, \code{\link[flapper]{acdc_simplify}} simplifies the outputs and \code{\link[flapper]{acdc_plot_record}} and \code{\link[flapper]{acdc_animate_record}} visualise the results.
#'
#' @examples
#' #### Step (1) Prepare study site grid
#' # Grid resolution needs to be sufficiently high to capture detection probability/movement
#' # And sufficiently low to minimise computational costs
#' blank <- raster::raster(raster::extent(dat_gebco), res = c(75, 75))
#' gebco <- raster::resample(dat_gebco, blank)
#'
#' #### Step (2) Implement setup_acdc_*() steps
#' # ... Define detection centroids required for algorithm(s) (see acs_setup_centroids())
#'
#' #### Step (3) Prepare movement time series for algorithm(s)
#' # Add required columns to dataframes:
#' dat_acoustics$timestamp_num <- as.numeric(dat_acoustics$timestamp)
#' dat_archival$timestamp_num  <- as.numeric(dat_archival$timestamp)
#' # Focus on an example individual
#' id <- 25
#' acc <- dat_acoustics[dat_acoustics$individual_id == id, ]
#' arc <- dat_archival[dat_archival$individual_id == id, ]
#' # Focus on the subset of data for which we have both acoustic and archival detections
#' acc <- acc[acc$timestamp >= min(arc$timestamp) - 2*60 &
#'              acc$timestamp <= max(arc$timestamp) + 2*60, ]
#' arc <- arc[arc$timestamp >= min(acc$timestamp) - 2*60 &
#'              arc$timestamp <= max(acc$timestamp) + 2*60, ]
#' # We'll focus on a one day period with overlapping detection/depth time series for speed
#' end <- as.POSIXct("2016-03-18")
#' acc <- acc[acc$timestamp <= end, ]
#' arc <- arc[arc$timestamp <= end, ]
#' arc <- arc[arc$timestamp >= min(acc$timestamp) - 2*60 &
#'              arc$timestamp <= max(acc$timestamp) + 2*60, ]
#'
#' #### Example (1) Implement AC algorithm with default arguments
#' out_acdc_1 <- flapper:::.acs(acoustics = acc,
#'                              bathy = gebco,
#'                              detection_centroids = dat_centroids,
#'                              mobility = 200
#'                              )
#'
#' #### Example (2) Implement ACDC algorithm with default arguments
#' out_acdc_2 <- flapper:::.acs(acoustics = acc,
#'                              archival = arc,
#'                              bathy = gebco,
#'                              detection_centroids = dat_centroids,
#'                              mobility = 200,
#'                              calc_depth_error = function(...) matrix(c(-2.5, 2.5), nrow = 2)
#'                              )
#'
#' #### Example (3) Implement AC or ACDC algorithm with detection probability kernels
#'
#' ## (A) Get detection centroid overlaps
#' # Define receiver locations as a SpatialPointsDataFrame object with a UTM CRS
#' proj_wgs84 <- sp::CRS(SRS_string = "EPSG:4326")
#' proj_utm   <- sp::CRS(SRS_string = "EPSG:32629")
#' rownames(dat_moorings) <- dat_moorings$receiver_id
#' xy <- sp::SpatialPoints(dat_moorings[, c("receiver_long", "receiver_lat")],
#'                         proj_wgs84)
#' xy <- sp::spTransform(xy, proj_utm)
#' xy <- sp::SpatialPointsDataFrame(xy, dat_moorings[, c("receiver_id",
#'                                                       "receiver_start_date",
#'                                                       "receiver_end_date")])
#' # Get detection overlap(s) as a SpatialPolygonsDataFrame
#' centroids <- get_detection_centroids(xy = sp::SpatialPoints(xy),
#'                                      detection_range = 425,
#'                                      coastline = dat_coast,
#'                                      byid = TRUE)
#' centroids_df <- dat_moorings[, c("receiver_id",
#'                                  "receiver_start_date",
#'                                  "receiver_end_date")]
#' row.names(centroids_df) <- names(centroids)
#' centroids <- sp::SpatialPolygonsDataFrame(centroids, centroids_df)
#' overlaps <- get_detection_centroids_overlap(centroids =  centroids)
#'
#' ## (B) Define detection probability function based on distance and detection_range
#' calc_dpr <-
#'   function(x){
#'     ifelse(x <= 425, stats::plogis(2.5 + -0.02 * x), 0)
#'   }
#'
#' ## (C) Get detection kernels (a slow step)
#' kernels <- acs_setup_detection_kernels(xy = xy,
#'                                        centroids = dat_centroids,
#'                                        overlaps = overlaps,
#'                                        calc_detection_pr = calc_dpr,
#'                                        bathy = gebco)
#' ## (D) Implement algorithm
#' out_acdc_3 <- flapper:::.acs(acoustics = acc,
#'                              archival = arc,
#'                              bathy = gebco,
#'                              detection_centroids = dat_centroids,
#'                              detection_kernels = kernels,
#'                              detection_kernels_overlap = overlaps$list_by_receiver,
#'                              detection_time_window = 10,
#'                              mobility = 200,
#'                              calc_depth_error = function(...) matrix(c(-2.5, 2.5), nrow = 2)
#'                              )
#'
#' #### Example (4): Compare outputs with/without detection probability and normalisation
#' ## Without detection kernels
#' out_acdc_4a <- flapper:::.acs(acoustics = acc,
#'                               archival = arc,
#'                               bathy = gebco,
#'                               detection_centroids = dat_centroids,
#'                               mobility = 200,
#'                               calc_depth_error = function(...) matrix(c(-2.5, 2.5), nrow = 2)
#'                               )
#' out_acdc_4b <- flapper:::.acs(acoustics = acc,
#'                               archival = arc,
#'                               bathy = gebco,
#'                               detection_centroids = dat_centroids,
#'                               normalise = TRUE,
#'                               mobility = 200,
#'                               calc_depth_error = function(...) matrix(c(-2.5, 2.5), nrow = 2)
#'                               )
#' ## With detection kernels
#' out_acdc_4c <- flapper:::.acs(acoustics = acc,
#'                               archival = arc,
#'                               bathy = gebco,
#'                               detection_centroids = dat_centroids,
#'                               detection_kernels = kernels,
#'                               mobility = 200,
#'                               calc_depth_error = function(...) matrix(c(-2.5, 2.5), nrow = 2)
#'                               )
#' out_acdc_4d <- flapper:::.acs(acoustics = acc,
#'                               archival = arc,
#'                               bathy = gebco,
#'                                detection_centroids = dat_centroids,
#'                               detection_kernels = kernels, normalise = TRUE,
#'                               mobility = 200,
#'                               calc_depth_error = function(...) matrix(c(-2.5, 2.5), nrow = 2)
#'                               )
#' ## Comparison of final maps
#' pp <- par(mfrow = c(2, 2))
#' raster::plot(out_acdc_4a$map, main = "4a")
#' raster::plot(out_acdc_4b$map, main = "4b")
#' raster::plot(out_acdc_4c$map, main = "4c")
#' raster::plot(out_acdc_4d$map, main = "4d")
#' par(pp)
#'
#' #### Example (5): Implement AC or ACDC algorithm and write messages to file via 'con'
#' out_acdc_5 <- flapper:::.acs(acoustics = acc,
#'                              archival = arc,
#'                              bathy = gebco,
#'                              detection_centroids = dat_centroids,
#'                              mobility = 200,
#'                              calc_depth_error = function(...) matrix(c(-2.5, 2.5), nrow = 2),
#'                              verbose = TRUE,
#'                              con = paste0(tempdir(), "/", "acdc_log.txt")
#'                              )
#' # Check log
#' utils::head(readLines(paste0(tempdir(), "/", "acdc_log.txt")))
#' utils::tail(readLines(paste0(tempdir(), "/", "acdc_log.txt")))
#'
#' #### Example (6): Implement AC or ACDC algorithm and return spatial information
#' # Specify save_record_spatial = NULL to include spatial information for all time steps
#' # ... (used for plotting) or a vector to include this information for specific time steps
#' out_acdc_6 <- flapper:::.acs(acoustics = acc,
#'                              archival = arc,
#'                              bathy = gebco,
#'                              detection_centroids = dat_centroids,
#'                              mobility = 200,
#'                              calc_depth_error = function(...) matrix(c(-2.5, 2.5), nrow = 2),
#'                              save_record_spatial = NULL,
#'                              verbose = TRUE,
#'                              con = paste0(tempdir(), "/", "acdc_log.txt")
#'                              )
#'
#'
#' @author Edward Lavender
#' @keywords internal

.acs <-
  function(
    acoustics,
    archival = NULL,
    step = 120,
    round_ts = TRUE,
    plot_ts = TRUE,
    bathy,
    map = NULL,
    detection_centroids,
    detection_kernels = NULL, detection_kernels_overlap = NULL, detection_time_window = 5,
    mobility,
    calc_depth_error = function(...) matrix(c(-2.5, 2.5), nrow = 2),
    normalise = FALSE,
    chunk = 1L,
    save_record_spatial = 1L,
    write_record_spatial_for_pf = NULL,
    save_args = TRUE,
    verbose = TRUE,
    con = "",
    progress = 1L,
    check = TRUE
  ){


    ######################################
    ######################################
    #### Set up

    #### A list to store overall outputs:
    # This includes:
    # ... dataframe highlighting time steps etc,
    # ... any saved spatial data
    # ... the final space use raster
    out <- list(map = NULL, record = NULL, time = NULL, args = NULL, chunks = NULL, simplify = FALSE)
    if(save_args) {
      out$args <- list(acoustics = acoustics,
                       archival = archival,
                       step = step,
                       round_ts = round_ts,
                       plot_ts = plot_ts,
                       bathy = bathy,
                       map = map,
                       detection_centroids = detection_centroids,
                       detection_kernels = detection_kernels,
                       detection_kernels_overlap = detection_kernels_overlap,
                       detection_time_window = detection_time_window,
                       mobility = mobility,
                       calc_depth_error = calc_depth_error,
                       normalise = normalise,
                       chunk = chunk,
                       save_record_spatial = save_record_spatial,
                       write_record_spatial_for_pf = write_record_spatial_for_pf,
                       save_args = save_args,
                       verbose = verbose,
                       con = con,
                       progress = progress,
                       check = check)
    }
    out$record <- list()

    #### Define function for printing messages to file or console
    ## Check the connection for writing files, if applicable
    if(check & con != ""){
      if(!verbose) {
        message("Input to 'con' ignored since verbose = FALSE.")
      } else {
        # Check directory
        check_dir(input = dirname(con))
        # Write black file to directory if required
        if(!file.exists(con)){
          message("'con' file does not exist: attempting to write file in specified directory...")
          file.create(file1 = con)
          message("Blank 'con' successfully written to file.")
        }
      }
    }
    ## Define function
    append_messages <- ifelse(con == "", FALSE, TRUE)
    cat_to_cf <- function(..., message = verbose, file = con, append = append_messages){
      if(message) cat(paste(..., "\n"), file = con, append = append)
    }

    #### Checks
    ## Formally initiate function and implement remaining checks
    t_onset <- Sys.time()
    cat_to_cf(paste0("flapper::.acs() called (@ ", t_onset, ")..."))
    out$time <- data.frame(event = "onset", time = t_onset)
    if(check){
      cat_to_cf("... Checking user inputs...")
      # Check acoustics contains required column names and correct variable types
      check_names(input = acoustics,
                  req = c("index", "timestamp", "timestamp_num", "receiver_id"),
                  extract_names = colnames,
                  type = all)
      check_class(input = acoustics$timestamp, to_class = "POSIXct", type = "stop")
      check_class(input = acoustics$receiver_id, to_class = "integer", type = "stop")
      # Check archival contains required column names and correct variable types
      if(!is.null(archival)){
        check_names(input = archival,
                    req = c("timestamp", "timestamp_num", "depth"),
                    extract_names = colnames,
                    type = all)
        check_class(input = archival$timestamp, to_class = "POSIXct", type = "stop")
        check_class(input = archival$depth, to_class = "numeric", type = "stop")
        # Check data volume
        if(nrow(archival) <= 1) stop("'archival' dataframe only contains one or fewer rows.")
        # Check archival step length
        step_est <- as.numeric(difftime(archival$timestamp[2], archival$timestamp[1], units = "s"))
        if(!all.equal(step, step_est)){
          stop("'step' does not equal difftime(archival$timestamp[2], archival$timestamp[1], units = 's'.)")
        }
      }
      # Check detection centroids have been supplied as a list
      check_class(input = detection_centroids, to_class = "list", type = "stop")
      # Focus on the subset of data for which we have both acoustic and archival detections
      if(!is.null(archival)){
        nrw_acc_pre <- nrow(acoustics)
        nrw_arc_pre <- nrow(archival)
        acoustics <- acoustics[acoustics$timestamp >= min(archival$timestamp) - step &
                                 acoustics$timestamp <= max(archival$timestamp) + step, ]
        archival <- archival[archival$timestamp >= min(acoustics$timestamp) - step, ]
        nrw_acc_post <- nrow(acoustics)
        nrw_arc_post <- nrow(archival)
        nrw_acc_delta <- nrw_acc_pre - nrw_acc_post
        nrw_arc_delta <- nrw_arc_pre - nrw_arc_post
        if(nrw_acc_post == 0 | nrw_arc_post == 0) stop("No overlapping acoustic/archival observations to implement algorithm.")
        if(nrw_acc_delta != 0) message(paste(nrw_acc_delta, "acoustic observation(s) beyond the ranges of archival detections ignored."))
        if(nrw_arc_delta != 0) message(paste(nrw_arc_delta, "archival observation(s) before the start of (processed) acoustic detections ignored."))
      }
      # Check detection kernels
      if(!is.null(detection_kernels)) {
        # Check input is as expected
        check_names(input = detection_kernels,
                    req = c("receiver_specific_kernels",
                            "receiver_specific_inv_kernels",
                            "array_design_intervals",
                            "bkg_surface_by_design",
                            "bkg_inv_surface_by_design"),
                    type = all)
        # Check example kernel matches the properties of bathy
        index_of_kernel_to_check <- which.min(sapply(detection_kernels$receiver_specific_kernels, is.null))
        if(!is.null(bathy)){
          raster_comparison <-
            tryCatch(
              raster::compareRaster(bathy, detection_kernels$receiver_specific_kernels[[index_of_kernel_to_check]]),
              error = function(e) return(e)
            )
          if(inherits(raster_comparison, "error")){
            warning(paste0("Checked detection kernel (detection_kernels$receiver_specific_kernels[[",
                           index_of_kernel_to_check,
                           "]]) and 'bathy' have different properties."),
                    immediate. = TRUE, call. = FALSE)
            stop(raster_comparison)
          }
        }
      }
      # Check depth error
      if(!is.null(archival)){
        de <- calc_depth_error(archival$depth)
        if(inherits(de, "matrix")){
          if(nrow(de) == 2){
            if(ncol(de) == 1) { message("'calc_depth_error' function taken to be independent of depth.")
            } else {
              message("'calc_depth_error' taken to depend on depth.")
            }
            if(any(de[1, ] > 0) | any(de[2, ] < 0)) stop("'calc_depth_error' should be a function that returns a two-row matrix with lower (negative) adjustment(s) (top row) and upper (positive) adjustment(s) (bottom row).'", call. = FALSE)
          } else stop("'calc_depth_error' should return a two-row matrix.", call. = FALSE)
        } else stop("'calc_depth_error' should return a two-row matrix.", call. = FALSE)
      }
      # Check write opts
      if(!is.null(write_record_spatial_for_pf)){
        check_named_list(input = write_record_spatial_for_pf)
        check_names(input = write_record_spatial_for_pf, req = "filename")
        write_record_spatial_for_pf$filename <- check_dir(input = write_record_spatial_for_pf$filename, check_slash = TRUE)
        if(length(list.files(write_record_spatial_for_pf$filename)) != 0L) {
          warning("write_record_spatial_for_pf$filename' is not an empty directory.", immediate. = TRUE, call. = FALSE)
        }
      }
    }

    #### Process time series
    if(round_ts){
      acoustics$timestamp <- lubridate::round_date(acoustics$timestamp, paste0(step/60, "mins"))
      if(is.null(archival)) archival$timestamp <- lubridate::round_date(archival$timestamp, paste(step/60, "mins"))
    }

    #### Visualise time series
    if(plot_ts) {
      ## AC algorithm implementation
      if(is.null(archival)){
        prettyGraphics::pretty_line(acoustics$timestamp,
                                    pch = 21, col = "royalblue", bg = "royalblue")
        ## ACDC implementation
      } else {
        axis_ls <- prettyGraphics::pretty_plot(archival$timestamp, abs(archival$depth)*-1,
                                               pretty_axis_args = list(side = 3:2),
                                               xlab = "Timestamp", ylab = "Depth (m)",
                                               type = "l")
        prettyGraphics::pretty_line(acoustics$timestamp,
                                    pretty_axis_args = list(axis_ls = axis_ls),
                                    inherit = TRUE,
                                    replace_axis = list(side = 1, pos = axis_ls[[2]]$lim[1]),
                                    add = TRUE,
                                    pch = 21, col = "royalblue", bg = "royalblue")
      }
    }

    #### Define the starting number of archival time steps
    # ... This will be updated to become the number of time steps moved since the start of the algorithm
    timestep_cumulative <- 0

    #### Define space use raster that will be updated inside this function
    if(is.null(map)) {
      map <- raster::setValues(bathy, 0)
      map <- raster::mask(map, bathy)
    }
    map_cumulative <- map
    # Directory in which to save files
    write_record_spatial_for_pf_dir <- write_record_spatial_for_pf$filename

    ##### Define 'uniform' detection probability across study area if detection kernels unsupplied
    if(is.null(detection_kernels)){
      kernel <- raster::setValues(bathy, 1)
      kernel <- raster::mask(kernel, bathy)
    }

    #### Define depth errors (if applicable)
    if(!is.null(archival)){
      archival$depth_lwr <- archival$depth + calc_depth_error(archival$depth)[1, ]
      archival$depth_upr <- archival$depth + calc_depth_error(archival$depth)[2, ]
    }

    ##### Wipe timestep_archival and timestep_detection from memory for safety
    timestep_archival <- NULL; timestep_detection <- NULL
    rm(timestep_archival, timestep_detection)

    #### Define progress bar
    # this will indicate how far along acoustic time steps we are.
    pb1 <- pb2 <- NULL
    if(progress %in% c(1, 3)) {
      pb1 <- utils::txtProgressBar(min = 0, max = (nrow(acoustics)-1), style = 3)
    }


    ######################################
    ######################################
    #### Move over acoustic time steps

    # For each acoustic time step
    # ... (except the last one - because we can't calculate where the individual
    # ... is next if we're on the last acoustics df)
    # ... we'll implement the algorithm
    out$time <- rbind(out$time, data.frame(event = "algorithm_initiation", time = Sys.time()))
    cat_to_cf("... Initiating algorithm: moving over acoustic and internal ('archival') time steps...")
    for(timestep_detection in seq_along(1:(nrow(acoustics)-1))){


      ######################################
      #### Define the details of the current and next acoustic detection

      #### Get the time steps
      # timestep_detection <- 1
      cat_to_cf(paste0("... On acoustic time step ('timestep_detection') ", timestep_detection, "."))

      #### Obtain the details of the previous acoustic detection (if applicable)
      # This is only applicable if index > 1.
      # If we are on the first detection time stamp (for a given chunk),
      # ... then we know that the previous receiver is the same as the current receiver
      # ... because this is enforced during chunk definition
      # ... (i.e., the acoustics dataframe has been split to ensure
      # ... that the previous receiver at which the  individual was detected
      # ... is ALWAYS the same as the current receiver).
      # If we are on a later detection time stamp, then we can define the previous receiver
      # ... from the objects created at the previous time step.
      index <- acoustics$index[timestep_detection]
      if(index > 1L) {
        if(timestep_detection == 1L){
          receiver_0_id <- acoustics$receiver_id[timestep_detection]
        } else {
          receiver_0_id <- receiver_1_id
        }
      }

      #### Obtain details of current acoustic detection
      # Define the time of the current acoustic detection
      receiver_1_timestamp <- acoustics$timestamp[timestep_detection]
      receiver_1_timestamp_num <- acoustics$timestamp_num[timestep_detection]
      # Identify the receiver at which the individual was detected:
      receiver_1_id <- acoustics$receiver_id[timestep_detection]

      #### Obtain details of next acoustic detection
      # Define time of next acoustic detection
      receiver_2_timestamp <- acoustics$timestamp[timestep_detection+1]
      receiver_2_timestamp_num <- acoustics$timestamp_num[timestep_detection+1]
      # Identify the next receiver at which the individual was detected:
      receiver_2_id      <- acoustics$receiver_id[timestep_detection+1]

      #### Duration between detections
      time_btw_dets <- as.numeric(difftime(receiver_2_timestamp, receiver_1_timestamp, units = "s"))


      ######################################
      #### Identify the number of archival records between the current and next acoustic detections

      #### AC implementation
      if(is.null(archival)){
        # Define the sequence of time steps
        pos <- seq(receiver_1_timestamp, receiver_2_timestamp, step)
        # Number of steps
        lpos_loop <- lpos <- length(pos)
        # If the last step coincides with the detection at the next receiver, we will drop this step
        # ... (we will account for this at the next time step instead)
        if(lpos > 1 & max(pos) == receiver_2_timestamp) lpos_loop <- lpos - 1

        #### ACDC implementation
      } else {
        # Identify the position of archival record which is closest to the current acoustic detection
        archival_pos1 <- which.min(abs(archival$timestamp_num - receiver_1_timestamp_num))
        # Define the positions of archival records between the two detections
        pos <- seq(archival_pos1, (archival_pos1 + time_btw_dets/step), by = 1)
        # Define the number of archival records between acoustic detections
        lpos_loop <- lpos <- length(pos)
        if(length(pos) > 1 & max(archival$timestamp[pos]) == receiver_2_timestamp) lpos_loop <- lpos - 1
      }


      ######################################
      #### Loop over archival time steps

      # Define a blank list in which we'll store the outputs of
      # ... looping over every archival time step
      als <- list()

      # Define another blank list in which we'll store the any saved spatial objects.
      # For each archival value, we'll add an element to this list which contains
      # ... the spatial objects at the value (if we've specified save_record_spatial)
      spatial <- list()

      # Define progress bar for looping over archival time steps
      # This is useful if there are a large number of archival time steps between
      # ... acoustic detections:
      # NB: only bother having a second progress bar if there are more than 1 archival time steps to move through.
      if(progress %in% c(2, 3) & lpos > 1) pb2 <- utils::txtProgressBar(min = 0, max = lpos, style = 3)

      # For each archival time step between our acoustic detections...
      for(timestep_archival in 1:lpos_loop){


        ######################################
        #### Identify the depth at the current time step:

        # Print the archival time step we're on.
        # timestep_archival <- 1
        cat_to_cf(paste0("... ... On internal time step ('timestep_archival') ", timestep_archival, "."))

        # Define timestep_cumulative: this keeps track of the total number of time steps
        # ... moved over the course of the algorithm.
        timestep_cumulative <- timestep_cumulative + 1

        # Identify depth at that time
        if(!is.null(archival)) {
          arc_ind   <- pos[timestep_archival]
          depth     <- archival$depth[arc_ind]
          depth_lwr <- archival$depth_lwr[arc_ind]
          depth_upr <- archival$depth_upr[arc_ind]
        }


        ######################################
        #### Implement centroid expansion and contraction

        #### At the moment of detection, get the detection centroids
        if(timestep_archival == 1){

          ## Previous centroid (perspective Ap)
          # If the current detection is at a new receiver (compared to the previous time step),
          # ... get the (expanded) detection centroid from the previous time step
          centroid_ap <- NULL
          if(index > 1){
            if(receiver_0_id != receiver_1_id)
              centroid_ap <- rgeos::gBuffer(centroid_c, width = mobility, quadsegs = 1000)
          }

          ## Detection centroid (perspective A or An)
          # Get the detection centroid around the current receiver
          centroid_an <- detection_centroids[[receiver_1_id]]

          ## Next centroid (perspective B)
          # If the individual was next detected at a different receiver,
          # ... we also need the (expanded) centroid around that receiver
          if(receiver_1_id != receiver_2_id){
            centroid_b <- detection_centroids[[receiver_2_id]]
            centroid_b <- rgeos::gBuffer(centroid_b, width = mobility * (lpos - 1), quadsegs = 1000)
          } else centroid_b <- NULL

          ## Individual's centroid (perspective C)
          # (i.e., the boundaries of the individuals location, given perspectives A, Ap and B)
          centroid_c <- get_intersection(list(centroid_an, centroid_ap, centroid_b))

          #### For subsequent time steps in-between detections
        } else {

          #### Dynamics (1): If the current and next receiver are identical
          # .... we simply expand/contract the centroid around this receiver
          if(receiver_1_id == receiver_2_id){

            ## Define blank centroids
            centroid_an <- centroid_ap <- centroid_b <- NULL

            ## Option 1: as the time step increases until halfway through acoustic detections:
            # ... increase the size of the centroid by the value of mobility (m) at each time step
            # Note the use of ceiling: if lpos is odd, then this takes us to the middle timestep_archival;
            # ... if lpos is even, then there are two middle timestep_archivals, and this takes us to the first one.
            if(timestep_archival <= ceiling(lpos/2)){
              cat_to_cf(paste0("... ... ... Acoustic centroid is expanding..."))
              centroid_c <- rgeos::gBuffer(centroid_c, width = mobility, quadsegs = 1000)

              ## Option 2: as the time step increases, from beyond halfway through acoustic detections,
              # ... to time of the next acoustic detection, we decrease the centroid size
            } else if(timestep_archival > ceiling(lpos/2)){

              ## Retain or shrink the centroid as appropriate
              ## For sequences with an even number of steps, the first post-middle value the stays the same
              # For example, if lpos = 4, then
              # ... centroid 1 = (say) 800 m (at detection)
              # ... centroid 2 =       1000 m
              # ... centroid 3 =       1000 m [first post-middle value - constant]
              # ... centroid 4 =       800 m (at next detection)
              # In this scenario, we don't need to change centroid.
              ## For sequences with an odd number of steps, the fist post-middle value starts to  shrink
              # For example, if lpos = 5,  then
              # ... centroid 1 = 800 m
              # ... centroid 2 = 1000 m
              # ... centroid 3 = 1200 m
              # ... centroid 4 = 1000 m [first post-middle value - shrinks]
              # ... centroid 5 = 800 m
              if(timestep_archival == (lpos/2 + 1)){
                cat_to_cf(paste0("... ... ... Acoustic centroid is constant ..."))
              } else {
                cat_to_cf(paste0("... ... ... Acoustic centroid is shrinking ..."))
                centroid_c <- rgeos::gBuffer(centroid_c, width = -mobility, quadsegs = 1000)
              }
            }

            #### Dynamics (2): If the two receivers are different
            # ... we keep expanding the current centroid
            # ... and contract the centroid of the other receiver
            # ... note that a high 'quadsegs' here is required to ensure correct calculations.
          } else {
            # Blank centroid_an
            centroid_an <- NULL
            # Expand existing centroid
            centroid_ap <- rgeos::gBuffer(centroid_c, width = mobility, quadsegs = 1000)
            # Contract centroid around receiver_id_2
            centroid_b <- rgeos::gBuffer(centroid_b, width = -mobility, quadsegs = 1000)
            # Get the intersection between the expanding centroid for receiver_1_id and the contracting centroid for receiver_2_id
            centroid_c <- rgeos::gIntersection(centroid_ap, centroid_b)
          }
        }


        ######################################
        #### Update map based on detection kernels

        #### Incorporate detection kernels, if supplied
        if(!is.null(detection_kernels)){

          #### Detection kernel at the moment of detection (timestep_archival == 1)
          # (Because we only go to one archival step before the next detection, the individual is only detected at timestep_archival == 1)
          # Option (1): Receiver detection kernel does not overlap with any other kernels: standard detection probability kernel
          # Option (2): Receiver detection kernel overlaps with other kernels
          # ... ... i): ... If there is only a detection at the one receiver, we down-weight overlapping regions
          # ... ... ii): ... If there are detections are multiple receivers, we upweight/downweight overlapping regions as necessary
          if(timestep_archival == 1){

            #### Define which option to implement
            # We start by assuming option 1 (no overlapping receivers)
            use_detection_kernel_option_1 <- TRUE
            # If receivers are well-synchronised and there are overlapping receivers, we will check whether or not
            # ... there are overlapping receivers for the current receiver
            # ... at the current time step
            if(!is.null(detection_time_window) & !is.null(detection_kernels_overlap)){

              ## Identify the full set of receivers whose deployments overlapped in space/time with receiver_1_id
              # ... from 'detection_kernels_overlap$list_by_receiver' list, derived from get_detection_centroids_overlap().
              # ... This is a list, with one element for each receiver.
              # ... Each element is a time series dataframe, with columns for the time and each receiver, that defines
              # ... whether or not the each receiver was active over that receivers deployment period.
              receiver_id_at_time_all <- detection_kernels_overlap[[receiver_1_id]]
              receiver_id_at_time_all <-
                receiver_id_at_time_all[receiver_id_at_time_all$timestamp == as.Date(receiver_1_timestamp), ]
              receiver_id_at_time_all$timestamp <- NULL
              receiver_id_at_time_all$receiver_id <- NULL

              ## If there are overlapping receivers, we need to account for the overlap
              if(any(receiver_id_at_time_all == 1)){
                use_detection_kernel_option_1 <- FALSE
                receiver_id_at_time_all <- as.integer(colnames(receiver_id_at_time_all)[receiver_id_at_time_all == 1])
                receiver_id_at_time_all <- data.frame(receiver_id = receiver_id_at_time_all,
                                                      detection = 0)
                # Check whether any overlapping receivers made a detection
                # Only do this if the current and next detection are effectively at the same time (for speed)
                if(time_btw_dets <= detection_time_window){
                  # Isolate all detections that occurred within the clock drift
                  acc_at_time <- acoustics %>% dplyr::filter(.data$timestamp >= receiver_1_timestamp - detection_time_window &
                                                               .data$timestamp <= receiver_1_timestamp + detection_time_window)
                  # Identify unique receivers at which detections occurred within the clock drift
                  if(nrow(acc_at_time) > 0){
                    receiver_id_at_time_all$detection[which(receiver_id_at_time_all$receiver_id %in%
                                                              unique(acc_at_time$receiver_id))] <- 1
                  }
                }
              }
            }

            #### Implement option 1 (no overlapping receivers --> standard detection probability kernel)
            if(use_detection_kernel_option_1){
              kernel <- detection_kernels$receiver_specific_kernels[[receiver_1_id]]

              #### Implement option 2 (overlapping receivers --> upweight/downweight kernels)
            } else{

              kernel <- detection_kernels$receiver_specific_kernels[[receiver_1_id]]
              for(i in 1:nrow(receiver_id_at_time_all)){
                if(receiver_id_at_time_all$detection[i] == 0){
                  kernel <- kernel * detection_kernels$receiver_specific_inv_kernels[[receiver_id_at_time_all$receiver_id[i]]]
                } else {
                  kernel <- kernel * detection_kernels$receiver_specific_kernels[[receiver_id_at_time_all$receiver_id[i]]]
                }
              }
            }
            # Re-scale kernels so that the most probable areas have a score of one
            kernel <- kernel/raster::cellStats(kernel, "max")

            #### At other time steps (timestep_archival != 1 ) (in between detections),
            # ... we simply upweight areas away from receivers relative to those within detection centroids
            # ... by using a single inverse detection probability surface calculated across all receivers
          } else {

            # Define an index to select the right kernel
            kernel_index <- which(as.Date(receiver_1_timestamp + (timestep_archival - 1)*step) %within%
                                    detection_kernels$array_design_intervals$array_interval)
            kernel <- detection_kernels$bkg_inv_surface_by_design[[kernel_index]]

          }
        }

        #### Define uniform probabilities across study site, if detection kernels unsupplied
        # Defined as 'kernel' at the start of the function.
        # Areas beyond the current centroid are masked below.


        ######################################
        #### Update map based on the location and depth.

        #### AC algorithm implementation only depends on detection probability kernels
        if(is.null(archival)){
          map_timestep <- raster::mask(kernel, centroid_c)
          if(is.null(detection_kernels)) map_timestep <- raster::mask(map_timestep, bathy)

          #### ACDC algorithm implementation also incorporates depth
        } else {

          # Identify the area of the bathymetry raster which is contained within the
          # ... allowed polygon. This step can be slow.
          bathy_sbt <- raster::mask(bathy, centroid_c)
          if(is.null(detection_kernels)) bathy_sbt <- raster::mask(bathy_sbt, bathy)

          # Identify possible position of individual at this time step based on depth ± depth_error m:
          # this returns a raster with cells of value 0 (not depth constraint) or 1 (meets depth constraint)
          map_timestep <- bathy_sbt >= depth_lwr & bathy_sbt <= depth_upr

          # Weight by detection probability, if necessary
          if(!is.null(detection_kernels)) map_timestep <- map_timestep * kernel
        }

        # Normalise probabilities within centroid, if specified
        if(normalise) map_timestep <- map_timestep/raster::cellStats(map_timestep, stat = "sum")

        # Add these positions to the map raster
        map_cumulative <- sum(map_cumulative,  map_timestep, na.rm = TRUE)
        # map_cumulative <- raster::mask(map_cumulative, bathy)

        # If the user has specified spatial information to be recorded along the way...
        # ... (e.g. for illustrative plots and/or error checking):
        if(is.null(save_record_spatial) | timestep_cumulative %in% save_record_spatial){
          # Create a list, plot options (po), which contains relevant spatial information:
          po <- list(centroid_ap    = centroid_ap,
                     centroid_an     = centroid_an,
                     centroid_b     = centroid_b,
                     centroid_c     = centroid_c,
                     kernel         = kernel,
                     map_timestep   = map_timestep,
                     map_cumulative = map_cumulative)

          # If the user hasn't selected to save the spatial record, simply define po = NULL,
          # ... so we can proceed to add plot options to the list() below without errors:
        } else po <- NULL


        ######################################
        #### Add details to temporary dataframe

        # Define dataframe
        dat <- data.frame(
          timestep_cumulative  = timestep_cumulative,
          timestep_detection   = timestep_detection,
          timestep_archival    = timestep_archival,
          receiver_1_id        = receiver_1_id,
          receiver_2_id        = receiver_2_id,
          receiver_1_timestamp = receiver_1_timestamp,
          receiver_2_timestamp = receiver_2_timestamp,
          time_btw_dets        = time_btw_dets
        )
        if(!is.null(archival)){
          dat$archival_timestamp   <- archival$timestamp[pos[timestep_archival]]
          dat$archival_depth       <- depth
        }

        # Add the dataframe to als
        als[[timestep_archival]] <- dat

        # Add plot options to the spatial list:
        spatial[[timestep_archival]] <- po

        # Write to file
        if(!is.null(write_record_spatial_for_pf)){
          write_record_spatial_for_pf$x <- map_timestep
          write_record_spatial_for_pf$filename <- paste0(write_record_spatial_for_pf_dir, "chu_", chunk, "_acc_", timestep_detection, "_arc_", timestep_archival)
          do.call(raster::writeRaster, write_record_spatial_for_pf)
        }

        # Update progress bar describing moment over archival time steps
        # only define title on the first archival time step out of the sequence
        # (this stops the progress bar being replotted on every run of this loop with the same title)
        if(progress %in% c(2, 3) & lpos > 1){
          utils::setTxtProgressBar(pb2, timestep_archival)
        }

      } # close for(j in 1:lpos){ (looping over archival time steps)
      if(!is.null(pb2)) close(pb2)

      #### Update overall lists
      als_df <- dplyr::bind_rows(als)
      out$record[[timestep_detection]] <- list(dat = als_df, spatial = spatial)

      # Update progress bar:
      if(progress %in% c(1, 3)) {
        utils::setTxtProgressBar(pb1, timestep_detection)
      }

    } # close for(i in 1:nrow(acoustics)){ (looping over acoustic time steps)
    if(!is.null(pb1)) close(pb1)

    #### Return function outputs
    cat_to_cf("... Movement over acoustic and internal ('archival') time steps has been completed.")
    t_end <- Sys.time()
    out$map <- map_cumulative
    # timings
    out$time <- rbind(out$time, data.frame(event = "algorithm_competion", time = t_end))
    out$time$serial_duration <- Tools4ETS::serial_difference(out$time$time, units = "mins")
    out$time$total_duration <- NA
    total_duration <- sum(as.numeric(out$time$serial_duration), na.rm = TRUE)
    out$time$total_duration[nrow(out$time)] <- total_duration
    cat_to_cf(paste0("... flapper::.acs() call completed (@ ", t_end, ") after ~", round(total_duration, digits = 2), " minutes."))
    class(out) <- c(class(out), "acdc_record")
    return(out)

  }

