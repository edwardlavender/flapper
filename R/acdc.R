######################################
######################################
#### acdc_setup_centroids()

#' @title Setup the acoustic centroids required for the ACDC algorithm
#' @description This function produces the acoustic contours required by the acoustic-centroid depth-contour (ACDC) algorithm.
#' @param rs A integer vector of receiver IDs.
#' @param xy A \code{\link[sp]{SpatialPoints}} object that defines the locations of each receiver. The order of points in this object should match the order of receivers defined in \code{rs}. The coordinate reference system should be the Universal Transverse Mercator system with distances in metres (to match \code{detection_range}, see below).
#' @param detection_range A number that defines the maximum detection range (m) at which an individual could be detected from a receiver.
#' @param mobility A number that defines the distance that an individual could move in the time period between archival observations.
#' @param n_timesteps An integer that defines the the number of timesteps after a hypothetical detection for which centroids will be created, where the duration of each timestep is given by the duration between archival observations.
#' @param coastline (optional) A \code{\link[sp]{SpatialPolygonsDataFrame-class}} object that defines the coastline in an area. If provided, acoustic centroids are processed to remove any areas on land. Algorithm speed declines with the complexity of the coastline.
#' @param boundaries (optional) A \code{\link[raster]{extent}} object that defines the boundaries of an area within which individuals are assumed to have remained. If provided, acoustic centroids are processed to remain within this area.
#' @param plot A logical input that defines whether or not to produce a plot of the area, including receivers, the coastline and the area boundaries (if provided), and acoustic centroids. This is useful for checking purposes but it can reduce algorithm speed.
#' @param cl (optional) A cluster object created by \code{\link[parallel]{makeCluster}} to implement the algorithm in parallel. The connection to the cluster is closed within the function.
#' @param verbose A logical input that defines whether or not to print messages to the console to relay function progress.
#'
#' @details Given detection at a particular receiver at a particular time, the detection range of the receiver and the movement speed of the animal, an acoustic centroid defines the possible locations of a detected individual at that time or a subsequent time, given only this information. More specifically, when an individual is located at a receiver, its location must be within some radius of a receiver defined by the maximum detection distance. This radius expands with the duration since detection. This function defines, for each receiver, a list of acoustic centroids that reflect the possible locations of an individual were it to have been detected from 0 to \code{n_timesteps} ago at that receiver, accounting for the coastline and within a defined area if necessary. Using the observed detection data, the ACDC and ACDCMP algorithms pull the relevant centroids out of this list, which substantially saves computation time because acoustic centroids are not computed on-the-fly. These centroids are processed within these algorithms (e.g., if an individual is detected at two different receivers, then at the halfway point between detections its location must be within the intersection of the relevant centroids for this two receivers) and combine them with other information to reconstruct where an individual could have been over time.
#'
#' @return The function returns a list of \code{\link[sp]{SpatialPolygonsDataFrame-class}} objects, with one element for all numbers from 1 to the maximum receiver number (\code{rx}). Any list elements that do not correspond to receivers contain a \code{NULL} element. List elements that correspond to receivers contain a \code{\link[sp]{SpatialPolygonsDataFrame-class}} object containing all the centroids for that receiver.
#'
#' @examples
#' #### Define data for acdc_setup_centroids()
#' ## Define coordinates of receivers as SpatialPoints with UTM CRS
#' # CRS of receiver locations as recorded in dat_moorings
#' proj_wgs84 <- sp::CRS("+init=epsg:4326")
#' # CRS of receiver locations required
#' proj_utm <- sp::CRS(paste("+proj=utm +zone=29 +datum=WGS84",
#'                           "+units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))
#' # Define SpatialPoints object
#' xy_wgs84 <- sp::SpatialPoints(dat_moorings[, c("receiver_long", "receiver_lat")], proj_wgs84)
#' xy_utm <- sp::spTransform(xy_wgs84, proj_utm)
#'
#' #### Example (1): Define a list of centroids with specified parameters
#' # ... (Argument values are small to reduce computation time for examples)
#' centroids <- acdc_setup_centroids(rs = dat_moorings$receiver_id,
#'                         xy = xy_utm,
#'                         detection_range = 500,
#'                         mobility = 250,
#'                         n_timesteps = 3
#'                         )
#' # A list of SpatialPolygonsDataFrames is returned with elements from 1:max(rs)
#' # NULL elements correspond to numbers in this sequence that do not refer to receivers
#' # Otherwise a SpatialPolygonsDataFrame is returned with all the centroids for that receiver
#' centroids
#'
#' #### Example (2): Visualise the acoustic centroids produced via plot = TRUE
#' centroids <- acdc_setup_centroids(rs = dat_moorings$receiver_id,
#'                         xy = xy_utm,
#'                         detection_range = 500,
#'                         mobility = 250,
#'                         n_timesteps = 3,
#'                         plot = TRUE
#'                         )
#'
#' #### Example (3): Remove areas of the centroids that overlap with coastline
#' centroids <- acdc_setup_centroids(rs = dat_moorings$receiver_id,
#'                         xy = xy_utm,
#'                         detection_range = 500,
#'                         mobility = 250,
#'                         n_timesteps = 3,
#'                         plot = TRUE,
#'                         coastline = dat_coast
#'                         )
#'
#' #### Example (4): Remove areas of the centroids beyond a boundary
#' xy_utm_coords <- sp::coordinates(xy_utm)
#' boundaries <- raster::extent(min(xy_utm_coords[, 1]),
#'                              max(xy_utm_coords[, 1]),
#'                              min(xy_utm_coords[, 2]),
#'                              max(xy_utm_coords[, 2])
#'                         )
#' centroids <- acdc_setup_centroids(rs = dat_moorings$receiver_id,
#'                         xy = xy_utm,
#'                         detection_range = 500,
#'                         mobility = 250,
#'                         n_timesteps = 3,
#'                         plot = TRUE,
#'                         coastline = dat_coast,
#'                         boundaries = boundaries
#'                         )
#'
#' #### Example (5): Implement the algorithm in parallel
#' centroids <- acdc_setup_centroids(rs = dat_moorings$receiver_id,
#'                         xy = xy_utm,
#'                         detection_range = 500,
#'                         mobility = 250,
#'                         n_timesteps = 3,
#'                         plot = TRUE,
#'                         coastline = dat_coast,
#'                         boundaries = boundaries,
#'                         cl = parallel::makeCluster(2L)
#'                         )
#'
#' #### Example (6): Acoustic centroids can be saved to file using rlist::list.save()
#' # rlist::list.save(centroids, paste0(tempdir(), "/centroids.RData"))

#' @author Edward Lavender
#' @export
#'

acdc_setup_centroids <- function(
  rs,
  xy,
  detection_range,
  mobility,
  n_timesteps = 250,
  coastline = NULL,
  boundaries = NULL,
  plot = FALSE,
  cl = NULL,
  verbose = TRUE
){

  #### Initiate function
  cat_to_console <- function(..., show = verbose) if(show) cat(paste(..., "\n"))
  cat_to_console("flapper::acdc_setup_centroids() called...")
  if(is.numeric(rs)) rs <- as.integer(rs)
  if(!is.integer(rs)) stop(paste("Argument 'rs' must be of class 'integer', not class(es):"), class(rs))
  if(any(rs <= 0)) stop("Argument 'rs' cannot contain receiver IDs <= 0.")
  if(any(duplicated(rs))){
    message("Argument 'rs' contains duplicate elements. rs has been simplified to unique(rs).")
    rs <- unique(rs)
  }
  if(!is.null(coastline) & !is.null(boundaries)) {
    coastline <- raster::crop(coastline, boundaries)
    if(is.null(coastline)) message("No coastline within defined boundaries. \n")
  }
  if(plot){
    cat_to_console("... Plotting background map of area...")
    if(!is.null(coastline)) {
      raster::plot(coastline, col = "lightgreen", border = "darkgreen", lwd = 1.5)
      graphics::points(xy, pch = 21, col = "royalblue", bg = "royalblue")
    } else {
      raster::plot(xy, pch = 21, col = "royalblue", bg = "royalblue")
    }
    if(!is.null(boundaries)) raster::lines(boundaries, col = "red", lty = 3)
  }


  #### Define a list of acoustic centroids for each receiver
  cat_to_console("... Building a nested list of acoustic centroids. This is the slow step...")

  ## Define a sequence of centroid radiuss
  # Around each receiver, we'll create a polygon of this radius
  radius_seq <- seq(detection_range, length.out = n_timesteps, by = mobility)

  ## Define a list of receivers, with a list of centroids for each receiver
  bathy_ls <- pbapply::pblapply(1:length(rs), cl = cl, function(i){

      centroids_ls <- lapply(radius_seq, function(radius){

        # Define a buffer around the current receiver of appropriate radius
        bathy_poly <- rgeos::gBuffer(xy[i], width = radius)

        # Reduce the size of the polygon by overlapping to remove areas on land
        # This keeps the polygons as small as possible which is important for ACDC/MP algorithm computation efficiency.
        if(!is.null(coastline)) {
          bathy_poly <- rgeos::gDifference(bathy_poly, coastline, byid = FALSE)
        }

        # Remove any areas beyond specified boundaries
        # Again this keeps polygon size to a minimum
        if(!is.null(boundaries)) {
          bathy_poly <- raster::crop(bathy_poly, boundaries)
        }

        # Return acoustic centroid
        return(bathy_poly)

      })

    # Define names of the rasters forn receiver, i, based on radius
    names(centroids_ls) <- paste0("s_", radius_seq)
    return(centroids_ls)

    })
  if(!is.null(cl)) parallel::stopCluster(cl)
  names(bathy_ls) <- rs

  #### Add NULL elements to the list for any receivers in the range 1:max(rs) that are not in rs
  # This means we can use receiver numbers to go straight to the correct element in the list in ACDC/MP algorithms.
  bathy_ls <- lapply(as.integer(1:max(rs)), function(i){
    if(i %in% rs){
      return(bathy_ls[[as.character(i)]])
    } else{
      return(NULL)
    }
  })

  #### Convert nested list of polygons to a SpatialPolygonsDataFrame
  # ... with one element for each receiver and each dataframe containing all the polygons for that receiver
  cat_to_console("... Converting the nested list of acoustic centroids to a SpatialPolygonsDataFrame...")
  if(is.null(bathy_ls)) {
    stop("There are no acoustic centroids within defined spatial boundaries.")
  }
  spdf_ls <- pbapply::pblapply(bathy_ls, cl = NULL, function(element){
    if(!is.null(element)){
      # bind all the sub-elements in each element together into a single spatial polygon
      sp <- raster::bind(element)
      # convert this into a spatial polygons dataframe
      spdf <- sp::SpatialPolygonsDataFrame(sp, data.frame(id = 1:length(sp)))
      return(spdf)
    }})

  #### Visualise centroids
  if(plot){
    cat_to_console("... Plotting centroids on map...")
    pbapply::pblapply(spdf_ls, function(spdf) if(!is.null(spdf)) raster::lines(spdf, col = "dimgrey", lwd = 0.75))
  }

  #### Return list of SpatialPolygonsDataFrame
  return(spdf_ls)

}



######################################
######################################
#### .acdc()

#' @title Back end implementation of the ACDC algorithm
#' @description This function is the back end of the ACDC algorithm.
#' @param acoustics A dataframe that contains passive acoustic telemetry detection time-series (see \code{\link[flapper]{dat_acoustics}} for an example). This should contain the following columns: an integer vector of receiver IDs, named 'receiver_id'; a POSIXct vector of timestamps when detections were made, named 'timestamp'; and a numeric vector of those timestamps, named 'timestamp_num'.
#' @param archival A dataframe that contains depth time-series (see \code{\link[flapper]{dat_archival}} for an example). This should contain the following columns: a numeric vector of observed depths, named 'depth'; a POSIXct vector of timestamps when observations were made, named 'timestamp'; and a numeric vector of those timestamps, named 'timestamp_num'. Depths should be recorded in the same units and with the same sign as the bathymetry data (see \code{bathy}). Absolute depths (m) are suggested. Unlike the detection time-series, archival timestamps are assumed to have occurred at regular intervals. Two-minute intervals are currently assumed.
#' @param bathy A \code{\link[raster]{raster}} that defines the bathymetry across the area within which the individual could have moved. This must be recorded in the same units and with the same sign as the depth observations (see \code{archival}). The coordinate reference system should be the Universal Transeverse Mercator system, with distances in metres (see also \code{\link[flapper]{acdc_setup_centroids}}).
#' @param map (optional) A blank \code{\link[raster]{raster}}, with the same properties (i.e., dimensions, resolution, extent and coordinate reference system) as the bathymetry raster (see \code{bathy}), but in which all values are 0. If \code{NULL}, this is computed internally, but supplying a pre-defined raster can be more computationally efficient if the function is applied iteratively (e.g., over different time windows).
#' @param detection_range A number that defines the maximum detection range (m) at which an individual could be detected from a receiver (see also \code{\link[flapper]{acdc_setup_centroids}}).
#' @param mobility A number that defines the distance (m) that an individual could move in the time period between archival observations (see also \code{\link[flapper]{acdc_setup_centroids}}).
#' @param depth_error A number that defines the interval around each depth (m) observation that defines the range of depths on the bathymetry raster (see \code{bathy}) that the individual could plausibly have occupied at that time. For example, \code{depth_error = 2.5} m implies that the individual could have occupied bathymetric cells whose depth lies within the interval defined by the observed depth (m) +/- 2.5 m. The appropriate value for \code{depth_error} depends on measurement error for the archival and bathymetry data, as well as the tidal range (m) across the area.
#' @param acc_centroids A list of acoustic centroids, with one element for each number from \code{1:max(acoustics$receiver_id)}, from \code{\link[flapper]{acdc_setup_centroids}}.
#' @param plot An integer vector that defines the time steps for which to return the necessary spatial information required to plot the plausible locations of the individual, given detection and depth time-series. \code{plot = 0} suppresses the return of this information and \code{plot = NULL} returns this information for all time steps. This spatial information can be used to plot time-specific results of the algorithm using \code{\link[flapper]{acdc_plot}}.
#' @param plot_ts A logical input that defines whether or not to the plot detection and depth time-series before the algorithm is initiated. This provides a useful visualisation of the extent to which they overlap.
#' @param verbose A logical variable that defines whether or not to print messages to the console or to file to relay function progress. If \code{con = ""}, messages are printed to the console; otherwise, they are written to file (see below).
#' @param con If \code{verbose = TRUE}, \code{con} is character string that defines the full pathway to a .txt file into which messages are written to relay function progress. This is approach, rather than printing to the console, is recommended for clarity, speed and debugging.
#' @param progress An integer (\code{1}, \code{2} or \code{3}) that defines whether or not to display a progress bar in the console as the algorithm moves over acoustic timesteps (\code{1}), the archival timesteps between each pair of acoustic detections (\code{2}) or both acoustic and archival timesteps (\code{3}), in which case the overall acoustic progress bar is punctuated by an archival progress bar for each pair of acoustic detections. This option is useful if there is a large number of archival observations between acoustic detections. Any other input will suppress the progress bar.
#' @param check A logical input that defines whether or not to check function inputs. This can be switched off to improve computation time when the function is applied iteratively.
#' @param ... Additional arguments (none implemented).
#'
#' @return The function returns a named list with the following elements: ‘map’, ‘record’, ‘time’ and ‘args’. The main output of the function is the ‘map’ RasterLayer that shows where the individual could have spent more or less time over the duration of the movement time-series. The ‘record’ element records time-specific maps of the possible locations of the individual, and can be used to plot maps of specific time points or to produce animations (for the time steps specified by \code{plot}). The ‘time’ element is a dataframe that defines the times of sequential stages in the algorithm's progression, providing a record of computation time; and the ‘args’ element is a named list of user inputs that record the parameters used to generate the outputs.
#'
#' @seealso \code{\link[flapper]{acdc_setup_centroids}} defines the acoustic centroids required by this function. \code{\link[flapper]{acdc_plot}} and \code{\link[flapper]{acdc_animate}} visualise the results.
#'
#' @examples
#' #### Step (1) Implement setup_acdc_*() steps
#' # ... Define acoustic centroids required for ACDC algorithm (see setup_acdc_centroids())
#'
#' #### Step (2) Prepare movement time-series for algorithm
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
#' #### Example (1) Implement ACDC algorithm with one default arguments
#' out_acdc <- .acdc(acoustics = acc,
#'                   archival = arc,
#'                   bathy = dat_gebco,
#'                   detection_range = 425,
#'                   mobility = 200,
#'                   depth_error = 2.5,
#'                   acc_centroids = dat_centroids,
#'                   )
#'
#' #### Example (2): Implement algorithm and write messages to file via 'con'
#' \dontrun{
#' out_acdc <- .acdc(acoustics = acc,
#'                   archival = arc,
#'                   bathy = dat_gebco,
#'                   space_use = NULL,
#'                   detection_range = 425,
#'                   mobility = 200,
#'                   depth_error = 2.5,
#'                   acc_centroids = dat_centroids,
#'                   verbose = TRUE,
#'                   con = paste0(tempdir(), "/", "acdc_log.txt")
#'                  )
#' # Check log
#' utils::head(readLines(paste0(tempdir(), "/", "acdc_log.txt")))
#' utils::tail(readLines(paste0(tempdir(), "/", "acdc_log.txt")))
#' }
#'
#' #### Example (3): Implement algorithm and return plotting information
#' # Specify plot = NULL to include plotting information for all time steps
#' # ... or a vector to include this information for specific time steps
#' \dontrun{
#' out_acdc <- .acdc(acoustics = acc,
#'                   archival = arc,
#'                   bathy = dat_gebco,
#'                   space_use = NULL,
#'                   detection_range = 425,
#'                   mobility = 200,
#'                   depth_error = 2.5,
#'                   acc_centroids = dat_centroids,
#'                   plot = NULL,
#'                   png_param = list(),
#'                   verbose = TRUE,
#'                   con = paste0(tempdir(), "/", "acdc_log.txt")
#'                   )
#' }
#'
#' @author Edward Lavender
#' @export
#'

.acdc <-
  function(
    acoustics,
    archival,
    bathy,
    map = NULL,
    detection_range,
    mobility,
    depth_error = 2.5,
    acc_centroids,
    plot = 0L,
    plot_ts = TRUE,
    verbose = TRUE,
    con = "",
    progress = 1L,
    check = TRUE,...
    ){


    ######################################
    ######################################
    #### Set up

    #### A list to store overall outputs:
    # This includes:
    # ... dataframe highlighting timesteps etc,
    # ... any saved spatial data
    # ... the final space use raster
    out <- list(map = NULL, record = NULL, time = NULL, args = NULL)
    out$args <- list(acoustics = acoustics,
                     archival = archival,
                     bathy = bathy,
                     map = map,
                     detection_range = detection_range,
                     mobility = mobility,
                     depth_error = depth_error,
                     acc_centroids = acc_centroids,
                     plot = plot,
                     verbose = verbose,
                     con = con,
                     progress = progress,
                     check = check,
                     dots = list(...))
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
    cat_to_cf(paste0("flapper::.acdc() called (@ ", t_onset, ")..."))
    out$time <- data.frame(event = "onset", time = t_onset)
    if(check){
      cat_to_cf("... Checking user inputs...")
      # Check acoustics contains required column names and correct variable types
      check_names(input = acoustics,
                  req = c("timestamp", "timestamp_num", "receiver_id"),
                  extract_names = colnames,
                  type = all)
      check_class(input = acoustics$timestamp, to_class = "POSIXct", type = "stop")
      check_class(input = acoustics$receiver_id, to_class = "integer", type = "stop")
      # Check archival contains required column names and correct variable types
      check_names(input = archival,
                  req = c("timestamp", "timestamp_num", "depth"),
                  extract_names = colnames,
                  type = all)
      check_class(input = archival$timestamp, to_class = "POSIXct", type = "stop")
      check_class(input = archival$depth, to_class = "numeric", type = "stop")
      # Check acoustic centroids have been supplied as a list
      check_class(input = acc_centroids, to_class = "list", type = "stop")
      # Focus on the subset of data for which we have both acoustic and archival detections
      nrw_acc_pre <- nrow(acoustics)
      nrw_arc_pre <- nrow(archival)
      acoustics <- acoustics[acoustics$timestamp >= min(archival$timestamp) - 2*60 &
                               acoustics$timestamp <= max(archival$timestamp) + 2*60, ]
      archival <- archival[archival$timestamp >= min(acoustics$timestamp) - 2*60, ]
      nrw_acc_post <- nrow(acoustics)
      nrw_arc_post <- nrow(archival)
      nrw_acc_delta <- nrw_acc_pre - nrw_acc_post
      nrw_arc_delta <- nrw_arc_pre - nrw_arc_post
      if(nrw_acc_post == 0 | nrw_arc_post == 0) stop("No overlapping acoustic/archival observations to implement algorithm.")
      if(nrw_acc_delta != 0) message(paste(nrw_acc_delta, "acoustic observation(s) beyond the ranges of archival detections ignored."))
      if(nrw_arc_delta != 0) message(paste(nrw_arc_delta, "archival observation(s) before the start of (processed) acoustic detections ignored."))
    }

    #### Visualise time-series
    if(plot_ts) {
      axis_ls <- prettyGraphics::pretty_plot(archival$timestamp, abs(archival$depth)*-1,
                                             pretty_axis_args = list(side = 3:2),
                                             xlab = "Timestamp", ylab = "Depth (m)",
                                             type = "l",
                                             return_list = TRUE)
      prettyGraphics::pretty_line(acoustics$timestamp,
                                  pretty_axis_args = list(axis_ls = axis_ls),
                                  inherit = TRUE,
                                  replace_axis = list(side = 1, pos = axis_ls[[2]]$lim[1]),
                                  add = TRUE,
                                  pch = 21, col = "royalblue", bg = "royalblue")
    }

    #### Define the starting number of archival timesteps
    # ... This will be updated to become the number of timesteps moved since the start of the algorithm
    timestep_cumulative <- 0

    #### Define space use raster that will be updated inside this function
    if(is.null(map)) {
      map <- bathy
      map <- raster::setValues(map, 0)
    }
    map_cumulative <- map

    #### Define radius_seq (the sequence of radius sizes that polygons take in acc_centroids)
    # Select the first element in the list which is not NULL:
    sele <- which(lapply(acc_centroids, function(element){!is.null(element)}) %>% unlist() %>% as.vector())[1]
    # Define the length of that element
    selel <- length(acc_centroids[[sele]])
    # Define radius seq using this information:
    radius_seq <- seq(detection_range, length.out = selel, by = mobility)
    # Define max radius:
    max_radius_seq <- max(radius_seq)

    ##### Wipe timestep_archival and timestep_detection from memory for safety
    timestep_archival <- NULL; timestep_detection <- NULL
    rm(timestep_archival, timestep_detection)

    #### Define progress bar
    # this will indicate how far along acoustic timesteps we are.
    if(progress %in% c(1, 3)) {
      pb1 <- utils::txtProgressBar(min = 0, max = (nrow(acoustics)-1), style = 3)
    }


    ######################################
    ######################################
    #### Move over acoustic timesteps

    # For each acoustic timestep
    # ... (except the last one - because we can't calculate where the individual
    # ... is next if we're on the last acoustics df
    # ... we'll implement the algorithm
    out$time <- rbind(out$time, data.frame(event = "algorithm_initiation", time = Sys.time()))
    cat_to_cf("... Initiating algorithm: moving over acoustic and archival timesteps...")
    for(timestep_detection in seq_along(1:(nrow(acoustics)-1))){


      ######################################
      #### Define the details of the current and next acoustic detection

      #### Print the acoustic timestep
      cat_to_cf(paste0("... On acoustic timestep ('timestep_detection') ", timestep_detection, "."))

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

      #### Duration between detections to the nearest 2 minutes (floored):
      time_btw_dets <- plyr::round_any(as.numeric(difftime(receiver_2_timestamp, receiver_1_timestamp, units = "mins")), 2)


      ######################################
      #### Identify the number of archival records between the current and next acoustic detections

      # Identify the position of archival record which is closest to
      # ... the current acoustic detection
      archival_pos1 <- which.min(abs(archival$timestamp_num - receiver_1_timestamp_num))

      # Define the positions of archival records between the two detections
      pos <- seq(archival_pos1, (archival_pos1 + time_btw_dets/2), by = 1) # assumes 2-minute timesteps

      # Define the number of archival records between acoustic detections
      lpos <- length(pos)


      ######################################
      #### Loop over archival timesteps

      # Define a blank list in which we'll store the outputs of
      # ... looping over every archival timestep.
      als <- list()

      # Define another blank list in which we'll store the any saved
      # ... spatial objects.
      # For each archival value, we'll add an element to this list which contains
      # ... the spatial objects at the value (if we've selected plot = TRUE)
      # plot options list... :
      spatial <- list()

      # Define progress bar for looping over archival timesteps
      # This is useful if there are a large number of archival timesteps between
      # ... acoustic detections:
      # NB: only bother having a second progress bar if there are more than 1 archival timesteps to move through.
      if(progress %in% c(2, 3) & lpos > 1) pb2 <- utils::txtProgressBar(min = 0, max = lpos, style = 3)

      # For each archival timestep between our acoustic detections...
      for(timestep_archival in 1:lpos){


        ######################################
        #### Identify the depth at the current timestep:

        # Print the archival timestep we're on.
        cat_to_cf(paste0("... ... On archival timestep ('timestep_archival') ", timestep_archival, "."))

        # Define timestep_cumulative: this keeps track of the total number of timesteps
        # ... moved over the course of the algorithm.
        timestep_cumulative <- timestep_cumulative + 1

        # Identify depth at that time
        depth <- archival$depth[pos[timestep_archival]]


        ######################################
        #### Identify the area in which the individual could be located at each timestep
        # ... based on time between the current and next location,
        # ... and the location of those two locations:

        #### Option 1: as the timestep increases until halfway through acoustic detections:
        # ... increase the size of the radius by the value of mobility (m) at each timestep
        # ... but keep the location of the current receiver as the centre of the radius in which we search.
        # Note the use of ceiling: if lpos is odd, then this takes us to the middle timestep_archival;
        # ...if lpos is even, then there are two middle t_dsts, and this takes us to the first one.
        if(timestep_archival == 1 | timestep_archival <= ceiling(lpos/2)){ # less than or equal to (more conservative):

          # The radius increases by detection_range + mobility for each timestep
          # ... that the individual is not detected by a receiver (until half way)
          radius <- detection_range + mobility * (timestep_archival - 1)
          if(radius > max_radius_seq) radius <- max_radius_seq
          # Print statement which explains the acoustic radius is constant or increasing...
          cat_to_cf(paste0("... ... ... Acoustic radius is constant or increasing (radius = ", radius, ")."))

          # Define approximate location in which individual must be located:
          # Determine position of radius in radius_seq
          radius_pos <- which(radius_seq == radius)
          # Extract appropriate acoustic centroid (polygon) from list of shapefiles based on receiver_1_id (!)
          centroid <- acc_centroids[[receiver_1_id]][radius_pos, ]


        #### Option 2: as the timestep increases, from beyond halfway through acoustic detections,
        # ... to time of the next acoustic detection:
        # ... We will shift the centre of the radius to be at the location of the receiver the individual was next detected at
        # ... And we'll shrink the radius gradually around this receiver until we reach ca. 800 m
        # ... (max detection range) at the time when the individual was actually detected
        } else if(timestep_archival > ceiling(lpos/2)){

          # If lpos is odd (e.g. lpos = 3), then if we're on lpos/2 + 1, we're on the second
          # ... halfway value and we'll keep the same radius, otherwise, we'll decrease the radius
          if(!(is.integer(lpos/2) & timestep_archival == ((lpos/2) + 1))){
            # radius is detection_range + mobility * multiplier which shrinks
            # ... as we move beyond halfway and towards the next acoustic detection
            radius <- detection_range + mobility * (lpos - timestep_archival)
            if(radius > max_radius_seq) radius <- max_radius_seq
          }
          cat_to_cf(paste0("... ... ... Acoustic radius is constant or decreasing (radius = ", radius, ")."))

          # If the receiver at which the individual has been detected is different from the one at
          # ... which it is next detected, then some adjustments are going to be necessary (see below).
          if(receiver_1_id != receiver_2_id & (timestep_archival == (ceiling(lpos/2) + 1))){
            # We need to copy the centroid of the previous timestep (centroid) into a new object (centroid_previous)
            # ... before this is replaced below for the current timestep.
            # (But we only need to do this if we're halfway between detections...)
            centroid_previous <- centroid
          }

          # Define approximate location in which individual must be located:
          # Determine position of radius in radius_seq
          radius_pos <- which(radius_seq == radius)

          # Extract appropriate spatial polygon from list of shapefiles
          # based on receiver_2_id (!), rather than receiver_1_id as above.
          # ... If the current and next receivers are identical, this doesnt make any difference
          # ... But if the next receiver is different, the centre of the acoustic area
          # ... shifts to the new location, and then starts to shrink around that location
          centroid <- acc_centroids[[receiver_2_id]][radius_pos,]
          centroid_next <- centroid

          # The code above derives the centroid within which the individual must be located based on acoustic data
          # ... Once we're halfway between two acoustic detections,
          # ... some additional steps are needed here if receiver_1_id != receiver_2_id.
          # ... Specifically, when we're half way between detections, the individual must be within the intersection
          # ... of the two centroid at the previous and next archival timestep (which may be less than the area of centroid)
          # ... At latter timesteps, this intersection area needs to grow by 200 m each timestep,
          # ... because that is how far the individual can move. But it can actually only be
          # ... in a portion of this area (the portion within the shrinking area centroid which defines
          # ... the maximum distance the individual can be from the next receiver in order to
          # ... 'get there in time' to be detected):
          if(receiver_1_id != receiver_2_id){
            # If we are halfway between acoustic detections....
            if(timestep_archival == (ceiling(lpos/2) + 1)){
              # Determine the overlap between centroid_previous and centroid: this is where the individual can be.
              centroid_overlap <- rgeos::gIntersection(centroid, centroid_previous)
              if(is.null(centroid_overlap)) {
                message("The algorithm is about to break due to an issue during the intersection of acoustic centroids.")
                message("Plotting the acoustic centroids to aid diagnosis before the function stops...")
                raster::plot(bathy)
                raster::plot(centroid, lwd = 2, add = TRUE)
                raster::plot(centroid_previous, lwd = 2, add = TRUE)
                message("Returning the outputs up to the previous timestamp before stopping...")
                return(out)
                msg <- paste0("The algorithm is halfway between two acoustic detections at two different receivers (",  receiver_1_id, " and ", receiver_2_id, ").",
                              "The acoustic centroid around receiver", receiver_1_id, " does not intersect with the centroid around receiver ", receiver_2_id, ", ",
                              "yet it is in the intersecting region that the individual is assumed to be at the halfway point between detections.",
                              "Either the maximum centroid size in 'acc_acoustics' is too small or the mobility parameter may be too low.")
                stop(msg)
              }
              # Copy centroid_overlap into centroid_overlap_expanded (see below):
              centroid_overlap_expanded <- centroid_overlap
              # Redefine centroid:
              centroid <- centroid_overlap
            # Else, if we're beyond halfway...
            } else if(timestep_archival > (ceiling(lpos/2) + 1)){
              # We need to grow centroid_overlap_expanded by c. 200 m at each timestep
              # ... compared to the last timestep (i.e. last saved version of centroid_overlap_expanded)
              centroid_overlap_expanded <- rgeos::gBuffer(centroid_overlap_expanded, width = mobility)
              # Now we need to work out the overlap between centroid_overlap_expanded and centroid
              # ... because the individual must be within the latter.
              # Adjust centroid accordingly:
              centroid <- rgeos::gIntersection(centroid_overlap_expanded, centroid)
            }
          }
        }


        ######################################
        #### Update map based on the location and depth.

        # Identify the area of the bathymetry raster which is contained within the
        # ... allowed polygon. This step can be slow.
        bathy_sbt <- raster::mask(bathy, centroid)

        # Identify possible position of individual at this time step based on depth ± depth_error m:
        # this returns a raster with cells of value 0 (not depth constraint) or 1 (meets depth constraint)
        map_timestep <- bathy_sbt >= (depth - depth_error) & bathy_sbt <= (depth + depth_error)

        # Add these positions to the map raster
        map_cumulative <- sum(map_cumulative,  map_timestep, na.rm = T)

        # If the user has specified spatial information to be recorded along the way...
        # ... (e.g. for illustrative plots and/or error checking):
        if(is.null(plot) | timestep_cumulative %in% plot){
          # Create a list, plot options (po), which contains relevant spatial information:
          po <- list(centroid = centroid,
                     map_timestep = map_timestep,
                     map_cumulative = map_cumulative)
          # Add other spatial information that may have been created along the way:
          if(receiver_1_id != receiver_2_id && timestep_archival >= (ceiling(lpos/2) + 1)){
            # add to po:
            po <- append(po,
                         list(
                           # previous centroid around receiver:
                           centroid_previous = centroid_previous,
                           # the centroid around the new receiver:
                           centroid_next = centroid_next,
                           # overlap between centroid_previous and centroid/centroid_next:
                           centroid_overlap = centroid_overlap,
                           # the grown overlap between centroids
                           centroid_overlap_expanded = centroid_overlap_expanded
                         )
            )
          }

        # If the user hasn't selected to produce plots, simply define po = NULL,
        # ... so we can proceed to add plot options to the list() below without errors:
        } else{
          po <- NULL
        }


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
          time_btw_dets        = time_btw_dets,
          archival_timestamp   = archival$timestamp[pos[timestep_archival]],
          archival_depth       = depth,
          centroid_radius      = radius
          )

        # Add the dataframe to als
        als[[timestep_archival]] <- dat

        # Add plot options to the spatial list:
        spatial[[timestep_archival]] <- po

        # Update progress bar describing moment over archival timesteps
        # only define title on the first archival timestep out of the sequence
        # (this stops the progress bar being replotted on every run of this loop with the same title)
        if(progress %in% c(2, 3) & lpos > 1){
          utils::setTxtProgressBar(pb2, timestep_archival)
        }

      } # close for(j in 1:lpos){ (looping over archival timesteps)

      #### Update overall lists
      als_df <- dplyr::bind_rows(als)
      out$record[[timestep_detection]] <- list(dat = als_df, spatial = spatial)

      # Update progress bar:
      if(progress %in% c(1, 3)) {
        utils::setTxtProgressBar(pb1, timestep_detection)
      }

    } # close for(i in 1:nrow(acoustics)){ (looping over acoustic timesteps)


    #### Return function outputs
    cat_to_cf("... Movement over acoustic and archival timesteps has been completed.")
    t_end <- Sys.time()
    out$map <- map_cumulative
    out$time <- rbind(out$time, data.frame(event = "algorithm_competion", time = t_end))
    out$time$serial_duration <- Tools4ETS::serial_difference(out$time$time, units = "mins")
    out$time$total_duration <- NA
    total_duration <- sum(as.numeric(out$time$serial_duration), na.rm = TRUE)
    out$time$total_duration[nrow(out$time)] <- total_duration
    cat_to_cf(paste0("... flapper::.acdc() call completed (@ ", t_end, ") after ~", round(total_duration, digits = 2), " minutes."))
    return(out)

  }


######################################
######################################
#### acdc()




######################################
######################################
#### acdc_plot()

#' @title Plot the results of the ACDC algorithm
#' @description This function is used to plot the results of the ACDC algorithm. To implement the function, a named list from \code{\link[flapper]{.acdc}} must be supplied, from which the results can be extracted and plotted. For each specified time-step, the function extracts the necessary information; sets up a blank background plot using \code{\link[raster]{plot}} and \code{\link[prettyGraphics]{pretty_axis}} and then adds requested spatial layers to this plot. Depending on user-inputs, this will usually show a cumulative map of where the individual could have spent more or less time, summed from the start of the algorithm to each time point. Coastline, receivers and acoustic centroids can be added and customised and the finalised plots can be returned or saved to file.
#' @param acdc A named list from \code{\link[flapper]{.acdc}}.
#' @param plot An integer vector that defines the time steps for which to make plots. If \code{plot = NULL}, the function will make a plot for all time steps for which the necessary information is available in \code{acdc}.
#' @param add_coastline (optional) A named list of arguments, passed to \code{\link[raster]{plot}}, to add a polygon (i.e., of the coastline), to the plot. If provided, this must contain an 'x' element that contains the coastline as a spatial object (e.g., a SpatialPolygonsDataFrame: see \code{\link[flapper]{dat_coast}} for an example).
#' @param add_receivers (optional) A named list of arguments, passed to \code{\link[graphics]{points}}, to add points (i.e., receivers) to the plot. If provided, this must contain an 'x' element that is a SpatialPoints object that specifies receiver locations (in the same coordinate reference system as other spatial data).
#' @param add_raster (optional) A named list of arguments, passed to \code{\link[fields]{image.plot}}, to plot the RasterLayer of possible locations that is extracted from \code{acdc}.
#' @param add_centroids (optional) A named list of arguments, passed to \code{\link[raster]{plot}}, to add the acoustic centroid to the plot.
#' @param add_additional (optional) A stand-alone function, to be executed after the background plot has been made and any specified spatial layers have been added to this, to customise the result (see Examples).
#' @param crop_spatial A logical variable that defines whether or not to crop spatial data to lie within the axis limits.
#' @param xlim,ylim,fix_zlim,pretty_axis_args Axis control arguments. \code{xlim} and \code{ylim} control the axis limits, following the rules of the 'lim' argument in \code{\link[prettyGraphics]{pretty_axis}}. \code{fix_zlim} is a logical input that defines whether or not to fix z axis limits across all plots (to facilitate comparisons), or a vector of two numbers that define a custom range for the z axis which is fixed across all plots. \code{fix_zlim = FALSE} produces plots in which the z axis is allowed to vary flexibly between time units. Other axis options supported by \code{\link[prettyGraphics]{pretty_axis}} are implemented by passing a named list of arguments to this function via \code{pretty_axis_args}.
#' @param par_param (optional) A named list of arguments, passed to \code{\link[graphics]{par}} to control the plotting window. This is executed before plotting is initiated and therefore affects all plots.
#' @param png_param (optional) A named list of arguments, passed to \code{\link[grDevices]{png}}, to save plots to file. If supplied, the plot for each time step is saved separately. The 'filename' argument should be the directory in which plots are saved. Plots are then saved as "1.png", "2.png" and so on.
#' @param cl A cluster object created by \code{\link[parallel]{makeCluster}}. If supplied, the function loops over specified time steps in parallel to make plots. This is only implemented if plots are saved to file (i.e., \code{png_param} is supplied). If supplied, the connection to the cluster is closed within the function.
#' @param verbose A logical variable that defines whether or not relay messages to the console to monitor function progress.
#' @param check A logical variable that defines whether or not to check user inputs to the function before its initation.
#' @param ... Additional arguments, passed to \code{\link[raster]{plot}}, to customise the blank background plot onto which spatial layers are added, such as \code{xlab}, \code{ylab} and \code{main}.
#'
#' @return The function plots the results of the ACDC algorithm at specified time steps, with one plot per time step. Plots are saved to file if \code{png_param} is supplied.
#' @examples
#' #### Example (1): The default options simply plot the first surface
#' acdc_plot(acdc = dat_acdc)
#'
#' #### Example (2): Define the number of plots to be produced and control the plotting window
#' acdc_plot(acdc = dat_acdc,
#'           plot = 1:2,
#'           par_param = list(mfrow = c(1, 2), mar = c(8, 8, 8, 8)))
#'
#' #### Example (3): Add and customise spatial information via add_* args
#' ## Define a SpatialPoints object of receiver locations
#' proj_wgs84 <- sp::CRS("+init=epsg:4326")
#' proj_utm <- sp::CRS(paste("+proj=utm +zone=29 +datum=WGS84",
#'                           "+units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))
#' rsp <- sp::SpatialPoints(dat_moorings[, c("receiver_long", "receiver_lat")], proj_wgs84)
#' rsp <- sp::spTransform(rsp, proj_utm)
#' ## Plot with receiver locations and coastline, customise the centroids and the raster
#' acdc_plot(acdc = dat_acdc,
#'           add_coastline = list(x = dat_coast, col = "darkgreen"),
#'           add_receivers = list(x = rsp, pch = 4, col = "royalblue"),
#'           add_centroids = list(col = "red"),
#'           add_raster = list(col = rev(topo.colors(100)))
#'           )
#'
#' #### Example (4): Control axis properties
#' # ... via smallplot argument for raster, pretty_axis_args, xlim, ylim and fix_zlim
#' # ... set crop_spatial = TRUE to crop spatial data within adjusted limits
#' acdc_plot(acdc = dat_acdc,
#'           add_coastline = list(x = dat_coast, col = "darkgreen"),
#'           add_receivers = list(x = rsp, pch = 4, col = "royalblue"),
#'           add_centroids = list(col = "red"),
#'           add_raster = list(smallplot= c(0.85, 0.9, 0.25, 0.75)),
#'           crop_spatial = TRUE,
#'           pretty_axis_args = list(side = 1:4,
#'                                   control_sci_notation = list(magnitude = 16L, digits = 0)),
#'           xlim = raster::extent(dat_coast)[1:2],
#'           ylim = raster::extent(dat_coast)[3:4],
#'           fix_zlim = c(0, 1)
#'           )
#'
#' #### Example (5): Modify each plot after it is produced via add_additional
#' # Specify a function to add titles to a plot
#' add_titles <- function(){
#'   mtext(side = 1, "x (UTM)", line = 2)
#'   mtext(side = 2, "y (UTM)", line = -8)
#' }
#' # Make plots with added titles
#' acdc_plot(acdc = dat_acdc,
#'           plot = 1:2,
#'           par_param = list(mfrow = c(1, 2), mar = c(8, 8, 8, 8)),
#'           add_coastline = list(x = dat_coast, col = "darkgreen"),
#'           add_receivers = list(x = rsp, pch = 4, col = "royalblue"),
#'           add_centroids = list(col = "red"),
#'           add_raster = list(),
#'           crop_spatial = TRUE,
#'           xlim = raster::extent(dat_coast)[1:2],
#'           ylim = raster::extent(dat_coast)[3:4],
#'           add_additional = add_titles
#'           )
#'
#' #### Example (6) Save plots via png_param
#' list.files(tempdir())
#' acdc_plot(acdc = dat_acdc,
#'           plot = 1:2,
#'           png_param = list(filename = tempdir())
#'           )
#' list.files(tempdir())
#' @author Edward Lavender
#' @export
#'

acdc_plot <- function(acdc,
                      plot = 1,
                      add_coastline = NULL,
                      add_receivers = NULL,
                      add_raster = list(col = rev(grDevices::terrain.colors(255))),
                      add_centroids = list(),
                      add_additional = NULL,
                      crop_spatial = FALSE,
                      xlim = NULL, ylim = NULL, fix_zlim = FALSE,
                      pretty_axis_args = list(side = 1:4,
                                              axis = list(list(),
                                                          list(),
                                                          list(labels = FALSE),
                                                          list(labels = FALSE)),
                                              control_axis = list(las = TRUE),
                                              control_sci_notation = list(magnitude = 16L, digits = 0)),
                      par_param = list(),
                      png_param = list(),
                      cl = NULL,
                      verbose = TRUE,
                      check = TRUE,...){

  #### Checks
  cat_to_console <- function(..., show = verbose) if(show) cat(paste(..., "\n"))
  cat_to_console("flapper::acdc_plot() called...")
  if(check){
    cat_to_console("... Checking function inputs...")
    ## Check plots to be produced
    if(any(plot <= 0L)) stop("Input to 'plot' must be > 0.")
    ## Check spatial data have been provided correctly
    if(!is.null(add_coastline)) {
      check_named_list(input = add_coastline)
      check_names(input = add_coastline, req = "x")
    }
    if(!is.null(add_receivers)) {
      check_named_list(input = add_receivers)
      check_names(input = add_receivers, req = "x")
      check_class(input = add_receivers$x, to_class = "SpatialPoints", type = "stop")
    }
    if(!is.null(add_raster)) check_named_list(input = add_raster)
    if(!is.null(add_centroids)) check_named_list(input = add_centroids)
    ## Check plotting window param
    check_named_list(input = par_param, ignore_empty = TRUE)
    ## Check png_param, if provided
    if(length(png_param) > 0){
      check_named_list(input = png_param)
      check_names(input = png_param, req = c("filename"), extract_names = names, type = all)
      png_param$filename <- check_dir(input = png_param$filename, check_slash = TRUE)
      save_png <- TRUE
    } else save_png <- FALSE
    ## Check cluster
    if(!is.null(cl) & is.null(png_param$filename)){
      message("Input to 'cl' ignored as png_param$filename is unspecified.")
      cl <- NULL
    }
    ## Check spatial information
    if(all(c(is.null(add_raster), is.null(add_coastline), is.null(add_receivers)))){
      stop("At least one argument out of 'add_raster', 'add_coastline' and 'add_receivers' should be provided.")
    }
    ## Other plotting param
    if(!is.null(pretty_axis_args$side)) {
      if(length(pretty_axis_args$side) == 1) stop("At least two sides in pretty_axis_args$side should be specified for a map.")
    }
    # Check dots
    check...("zlim",...)
  }

  #### Define data for background plot
  cat_to_console("... Defining data for background plot...")
  ## Unpack information for plotting and isolate relevant plots
  acdc_plot <- lapply(acdc$record, function(elm) elm$spatial)
  acdc_plot <- lapply(acdc_plot, function(x) if(length(x) > 0) return(x))
  acdc_plot <- purrr::flatten(acdc_plot)
  if(!is.null(plot)) acdc_plot <- acdc_plot[plot]
  acdc_plot <- acdc_plot[which(!sapply(acdc_plot, is.null))]
  if(is.null(plot)) plot <- 1:length(acdc_plot)
  if(check) {
    if(length(acdc_plot) <= 0) stop("No plotting data available for selected plot(s).")
    if(any(plot > length(acdc_plot))) stop("Some inputs to 'plot' larger than the number of available plots.")
  }

  ## Extent of area
  if(!is.null(add_raster)) {
    first_raster <- acdc_plot[[1]]$map_cumulative
    ext_ras <- raster::extent(first_raster)
    x_ras   <- ext_ras[1:2]
    y_ras   <- ext_ras[3:4]
  } else{
    x_ras <- NULL
    y_ras <- NULL
  }
  ## Extent of coastline provided
  if(!is.null(add_coastline)) {
    ext_coastline <- raster::extent(add_coastline$x)
    x_coastline   <- ext_coastline[1:2]
    y_coastline   <- ext_coastline[3:4]
  } else {
    x_coastline <- NULL
    y_coastline <- NULL
  }
  ## Extent of receiver locations
  if(!is.null(add_receivers)) {
    rxy <- add_receivers$x
    rxy <- sp::coordinates(rxy)
    x_rxy <- range(rxy[, 1])
    y_rxy <- range(rxy[, 2])
  } else {
    x_rxy <- NULL
    y_rxy <- NULL
  }
  ## Define x and and y values to the extremes of these possible limits
  x <- range(x_ras, x_coastline, x_rxy)
  y <- range(y_ras, y_coastline, y_rxy)
  ## Define axis limits
  axis_param <- prettyGraphics::implement_pretty_axis_args(x = list(x, y),
                                                           pretty_axis_args = pretty_axis_args,
                                                           xlim = xlim,
                                                           ylim = ylim)
  xlim <- axis_param[[1]]$lim
  ylim <- axis_param[[2]]$lim
  ext <- raster::extent(xlim, ylim)
  ## Crop spatial data within limits
  if(!is.null(add_coastline)) if(crop_spatial) add_coastline$x <- raster::crop(add_coastline$x, ext)

  ## Define zlim, if requested
  if(!is.null(add_raster)) {
    if(is.logical(fix_zlim)) {
      if(fix_zlim){
        range_use <- lapply(acdc_plot, function(map_info){
          min_use <- raster::minValue(map_info$map_cumulative)
          max_use <- raster::maxValue(map_info$map_cumulative)
          return(c(min_use, max_use))
        })
        range_use <- do.call(rbind, range_use)
        min_use <- min(range_use[, 1])
        max_use <- max(range_use[, 2])
        zlim <- c(min_use, max_use)
      }
    }
  }


  #### Define plotting window
  cat_to_console("... Setting plotting window...")
  pp <- do.call(graphics::par, par_param)

  #### Loop over every detection
  cat_to_console("... Making plots for each timestep ...")
  pbapply::pblapply(1:length(acdc_plot), cl = cl, function(i){

    #### Set up image to save
    if(save_png){
      title <- paste0(i, ".png")
      png_param_tmp <- png_param
      png_param_tmp$filename <- paste0(png_param_tmp$filename, title)
      do.call(grDevices::png, png_param_tmp)
    }

    #### Define background plot
    map_info <- acdc_plot[[i]]
    map      <- map_info$map_cumulative
    area <- sp::SpatialPoints(raster::coordinates(ext), raster::crs(first_raster))
    raster::plot(x = area, xlim = xlim, ylim = ylim,...)

    #### Define time-specific zlim, if requested
    if(!is.null(add_raster)) {
      if(is.logical(fix_zlim)){
        if(!fix_zlim) {
          min_use <- raster::minValue(map)
          max_use <- raster::maxValue(map)
          zlim <- c(min_use, max_use)
        }
      } else {
        zlim <- fix_zlim
      }
    }

    #### Add spatial objects
    # Add spatial use surface
    if(!is.null(add_raster)) {
      add_raster$x <-  map
      add_raster$zlim <- zlim
      if(crop_spatial) add_raster$x <- raster::crop(add_raster$x, ext)
      add_raster$add <- TRUE
      do.call(fields::image.plot, add_raster)
    }
    # Add acoustic centroid
    if(!is.null(add_centroids)) {
      add_centroids$x   <- map_info$centroid
      if(crop_spatial) add_centroids$x <- raster::crop(add_centroids$x, ext)
      add_centroids$add <- TRUE
      do.call(raster::plot, add_centroids)
    }
    # Add the coastline (note that the coastline has already been cropped, if necessary)
    if(!is.null(add_coastline)) {
      add_coastline$add <- TRUE
      do.call(raster::plot, add_coastline)
    }
    # Add receivers
    if(!is.null(add_receivers)) {
      do.call(graphics::points, add_receivers)
    }
    # Add additional
    if(!is.null(add_additional)) {
      add_additional()
    }

    #### Add axes back at end
    prettyGraphics::pretty_axis(axis_ls = axis_param, add = TRUE)

    #### Save fig
    if(save_png) grDevices::dev.off()
  })
  if(!is.null(cl)) parallel::stopCluster(cl)

  invisible()

}


######################################
######################################
#### acdc_animate()

#' @title Create a html animation of the ACDC algorithm
#' @description This function is a simple wrapper for \code{\link[flapper]{acdc_plot}} and \code{\link[animation]{saveHTML}} which creates an animation of the ACDC algorithm over time. To implement this function, a named list of arguments for \code{\link[flapper]{acdc_plot}}, which creates the plots, must be supplied. This is embedded within \code{\link[animation]{saveHTML}}, which creates a folder in the working directory named 'images' that contains a .png file for each time step and an animation as a .html file.
#' @param expr_param A named list of arguments, passed to \code{\link[flapper]{acdc_plot}}, to create plots.
#' @param html_name A string that defines the name of the html file (see 'htmlfile' argument in \code{\link[animation]{saveHTML}}).
#' @param image_name A string that defines the names of the individual .png files creates (see 'img.name' argument in \code{\link[animation]{saveHTML}}).
#' @param html_title,html_description Character strings that provide a title and a description that are displayed within the html animation (see 'title' and 'description' arguments in \code{\link[animation]{saveHTML}}).
#' @param navigator A logical variable that defines whether or not to add a navigator panel to the animation (see 'navigator' argument in \code{\link[animation]{saveHTML}}).
#' @param ani_height,ani_width,ani_res Numbers that define the size and the resolution of the animation (see 'ani.height' 'ani.width' and 'ani.res' arguments in \code{\link[animation]{ani.options}}).
#' @param interval A number that defines the time interval between sequential frames (see 'interval' argument in \code{\link[animation]{ani.options}}).
#' @param verbose A logical or character variable that defines whether or not, or what, to write as a footer to the html animation (see 'verbose' argument in \code{\link[animation]{ani.options}}).
#' @param ... Additional arguments passed to \code{\link[animation]{ani.options}}.
#'
#' @return The function produces an animation in .html format in the working directory (or a sub-directory of this). A folder 'images' is also produced which contains the images for each time step. A 'css' and 'js' folder are also produced by \code{\link[animation]{saveHTML}} which creates the animation.
#'
#' @examples
#' dir_current <- getwd()
#' setwd(tempdir())
#' acdc_animate(expr_param = list(acdc = dat_acdc,
#'                                add_coastline = list(x = dat_coast, col = "darkgreen"),
#'                                plot = 1:5,
#'                                fix_zlim = FALSE)
#' )
#' setwd(dir_current)
#' @details This function requires the \code{\link[animation]{animation}} package.
#' @author Edward Lavender
#' @export
#'

acdc_animate <-
  function(expr_param,
           html_name = "ACDC_algorithm_demo.html",
           image_name = "ACDC",
           html_title = "Demonstration of the ACDC Algorithm",
           html_description = "",
           navigator = FALSE,
           ani_height = 800,
           ani_width = 800,
           ani_res = 1200,
           interval = 0.1,
           verbose = FALSE,
           ...){

    if (!requireNamespace("pbapply", quietly = TRUE)) {
      stop("This function requires the 'animation' package. Please install it before continuing with install.packages('animation').")
    }
    animation::saveHTML({
      do.call(acdc_plot, expr_param)
    },
    htmlfile = html_name,
    img.name = image_name,
    title = html_title,
    description = html_description,
    navigator = navigator,
    ani.height = ani_height,
    ani.width = ani_width,
    ani.res = ani_res,
    interval = interval,
    verbose = verbose,...
    )
  }
