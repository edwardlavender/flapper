######################################
######################################
#### acdc_setup_mobility()

#' @title Examine the constant `mobility' assumption of the AC* algorithm(s)
#' @description In the simplest and fastest (and only) version of the acoustic-centroid (AC) and acoustic-centroid depth-contour (ACDC) algorithms currently implemented by \code{\link[flapper]{flapper}}, the rate at which acoustic centroids expand and contract depends on a single `mobility' parameter that describes how the uncertainty in an individual's location changes through time according to the passive acoustic telemetry data. These changes essentially reflect the maximum horizontal distance that an individual can move in the regularised time steps between sequential detections. However, in some situations, a fixed parameter for horizontal movement may not be well supported. For instance, in cases with archival data (the ACDC algorithm), animals changes dramatically in depth and for which changes in depth are likely to curtail the extent of horizontal movement*. Thus, this function investigates the extent to which the horizontal distance an animal could travel changes through time if the mobility parameter is instead conceptualised as the diagonal distance between sequential depth observations.
#'
#' @param depth A vector of depth observations, whose units match the \code{mobility} parameter, to be incorporated into the ACDC algorithm. (Note that the ACDC algorithm may drop depth observations internally depending on their alignment with the acoustic time series and so, ideally, pre-processed time series should be passed to \code{depth} to ensure inferences correspond directly to the time series modelled by \code{\link[flapper]{acdc}}.) Depth observations should be regularly spaced in time (i.e., represent a time series for a single individual, without gaps).
#' @param mobility A number, in the same units as \code{depth}, that defines the maximum horizontal distance that an individual could move in the time period between archival observations (see \code{\link[flapper]{acdc_setup_centroids}}).
#' @param plot A logical variable that defines whether or not to plot the distribution of horizontal distances that the individual could travel, given its depth time series and \code{mobility}, as a histogram.
#' @param add_mobility (optional) If \code{plot = TRUE}, \code{add_mobility} is a named list of arguments, passed to \code{\link[graphics]{abline}}, to add the \code{mobility} parameter as a vertical line to the histogram as a reference point. \code{add_mobility = NULL} suppresses this line.
#' @param add_rug (optional) If \code{plot = TRUE}, \code{add_rug} is a named list of arguments, passed to \code{\link[graphics]{rug}}, to add the estimated distances to the histogram as a rug. \code{add_rug = NULL} suppresses this rug.
#' @param xlab,... Additional arguments, passed to \code{\link[prettyGraphics]{pretty_hist}}, to customise the histogram.
#'
#' @details The function uses a Pythagorean approximation to estimate the \eqn{distance} that an individual could travel in each time step (\eqn{t}), given the maximum horizontal distance that the individual could move (\eqn{mobility}) and sequential changes in \eqn{depth}, according to the equation:
#'
#' \deqn{distance_{t + 1} = \sqrt{mobility^2 + (depth_{t + 1} - depth_{t})^2},}
#'
#' where \eqn{depth_{t + 1}} for the final (\eqn{n^{th}}) observation is defined as \eqn{depth_{t = n}} and thus \eqn{distance_{t = n} = mobility}.
#'
#' *If the horizontal distances that an individual could travel are not very variable, then the benefits of a single mobility parameter are likely to outweigh the costs. On the other hand, substantial variation in the horizontal distances that an individual could travel may suggest that pre-processed acoustic centroids are inappropriate; in this situation, the correct centroids could be computed on-the-fly within the ACDC algorithm, but this is not currently implemented. However, the particle filtering algorithms can account for this by the incorporation of movement models within/between acoustic centroids.
#'
#' @return The function returns a numeric vector of distances and, if \code{plot = TRUE}, a histogram of those distances.
#'
#' @examples
#' #### Example (1) Explore mobility for single individual
#' mob <-
#'   acdc_setup_mobility(depth = dat_archival$depth[dat_archival$individual_id == 25],
#'                       mobility = 200)
#'
#' #### Example (2) Customise histogram
#' # ... suppress plot with plot =  FALSE
#' # ... suppress add_* lists with NULL or customise
#' # ... pass args to prettyGraphics::pretty_hist() via ...
#' mob <-
#'   acdc_setup_mobility(depth = dat_archival$depth[dat_archival$individual_id == 25],
#'                       mobility = 200,
#'                       add_mobility = NULL,
#'                       add_rug = list(side = 3, pos = 25000,
#'                                      col = "darkred", lwd = 0.25, ticksize = 0.01),
#'                       breaks = 25)
#'
#' #### Example (3) Explore mobility for mulitple individuals
#' pp <- graphics::par(mfrow = c(2, 2))
#' mob_ls <- lapply(unique(dat_archival$individual_id), function(id){
#'   acdc_setup_mobility(depth = dat_archival$depth[dat_archival$individual_id == id],
#'                       mobility = 200)
#' })
#' graphics::par(pp)
#'
#' #### Results
#' # For these sample time series, even though the animals change depth
#' # ... (quite substantially) though time, the sequential changes in depth
#' # ... are small and the influence on the horizontal distances that they
#' # ... are assumed to be able to travel in the time gap between archival
#' # ... observations is so small that a constant mobility parameter
#' # ... is reasonable without further information (e.g., models of the underlying
#' # ... behavioural state of the animal). However, it is still likely to be
#' # ... beneficial to include a movement model to join locations within/
#' # ... between acoustic centroids for some applications via particle filtering.
#'
#' @seealso \code{\link[flapper]{acdc_setup_mobility}}, \code{\link[flapper]{acdc_setup_n_centroids}},  \code{\link[flapper]{acdc_setup_centroids}} and \code{\link[flapper]{acdc_setup_detection_kernels}} are used to set up the AC and ACDC algorithms as implemented by \code{\link[flapper]{ac}} and \code{\link[flapper]{acdc}}.
#'
#' @author Edward Lavender
#' @export
#'

acdc_setup_mobility <- function(depth,
                                mobility,
                                plot = TRUE,
                                add_mobility = list(col = "royalblue", lwd = 2),
                                add_rug = list(),
                                xlab,...){

  #### Calculate mobility
  # depths at each time step
  dat <- data.frame(depth_t1 = depth)
  # depths at next time step
  dat$depth_t2 <- dplyr::lead(dat$depth_t1)
  # force the last depth to be the same as the previous depth (i.e., fall back on mobility value)
  dat$depth_t2[nrow(dat)] <- dat$depth_t1[nrow(dat)]
  # calculate mobility via Pythagoras' Theorem
  dat$depth_delta_sq <- (dat$depth_t2 - dat$depth_t1)^2
  mobility_sq <- mobility^2
  dat$mobility <- sqrt(mobility_sq + dat$depth_delta_sq)

  #### Make plot
  if(plot) {
    # Histogram
    if(missing(xlab)) xlab <- "Mobility (m)"
    prettyGraphics::pretty_hist(x = dat$mobility, xlab = xlab,...)
    # Add line for mobility
    if(!is.null(add_mobility)) {
      add_mobility$v <- mobility
      do.call(graphics::abline, add_mobility)
    }
    # Add rug
    if(!is.null(add_rug)) {
      add_rug$x <- dat$mobility
      if(is.null(add_rug$pos)) add_rug$pos <- 0
      do.call(graphics::rug, add_rug)
    }
  }

  #### Return mobility vector
  return(dat$mobility)
}


######################################
######################################
#### acdc_setup_n_centroids()

#' @title Suggest the number of centroids for the AC* algorithm(s)
#' @description The acoustic-centroid (AC) and acoustic-centroid depth-contour (ACDC) algorithms require a list of acoustic centroids from \code{\link[flapper]{acdc_setup_centroids}}. The number of centroids that is required depends on (a) the duration between detections; (b) the distances among receivers; (c) the area of interest; and (d) other considerations such as the spatial extent of the study area. This function implements three methods to facilitate a sensible choice of the number of centroids for the AC* algorithm(s) (see Details). This ensures that the number of centroids is sufficient to cover requirements but not so large that the centroids become slow to compute and unwieldy. The function requires a time series of detections, the time step length, the deployment details of passive acoustic telemetry receivers, a mobility parameter and a Spatial* or Raster* layer that defines the boundaries of the study site. Given these inputs, the function returns suggested bounds for the number of centroids.
#'
#' @param detections For Method One, \code{detections} is a POSIXct vector of time stamps when detections were made (for a particular individual).
#' @param step For Method One, \code{step} is the duration (s) between sequential time steps inbetween acoustic detections. For the ACDC algorithm, this is the duration between sequential archival observations.
#' @param hist,... For Method One, \code{hist} is a logical value that defines whether or not to plot the distribution of gaps across the detection time series as a histogram and \code{...} includes any additional arguments passed to \code{\link[prettyGraphics]{pretty_hist}}.
#' @param moorings For Method Two, \code{moorings} is a dataframe that defines passive acoustic telemetry receiver locations and deployment periods. This must contain the following columns: `receiver_id', a unique identifier of each receiver; `receiver_lat' and `receiver_long', the latitude and longitude of each receiver in decimal degrees; and `receiver_start_date' and `receiver_end_date', the start and end time of each receiver's deployment period (see \code{\link[flapper]{dat_moorings}} for an example).
#' @param mobility For Methods Two and Three, \code{mobility} is a number that defines the distance (m) that an individual could move in the time period in the time steps between acoustic observations (see also \code{\link[flapper]{acdc_setup_centroids}}).
#' @param double For Method Two, \code{double} is a logical variable that defines whether or not to double the minimum number of centroids (see Details).
#' @param boundaries For Method Three, a Spatial* or Raster* object that defines the study site, from which area boundaries can be extracted (via \code{\link[raster]{extent}}).
#'
#' @details
#' \subsection{Method (1)}{The first method provides a reasonable upper bound for the number of centroids. This method is based on the gaps between detections. During this time, acoustic centroids increase in size, reflecting the increasing uncertainty in the location of an individual, until the half way point between detections. At this point, uncertainty in the individual's location is maximised and the acoustic centroids reach their maximum size. Thereafter, uncertainty in the individual's location and the size of the acoustic centroids decrease towards the receiver at which the individual was next detected. Under this perspective, the minimum number of centroids is the number required such that the centroids continue to increase in size until the halfway point between the two most temporally distance detections. For example, if the time steps between acoustic detections are two-minutes in duration, the minimum number of centroids is one quarter of the duration (minutes) of the longest gap between acoustic detections. To implement this method, a time series of \code{detections} (for the individual and time period for which the AC* algorithm(s) will be implemented) and the time \code{step} length need to be provided. This method provides a sensible upper bound for the number of centroids, since it allows centroids to continue to expand to their maximum possible size, given the data and the assumptions of the algorithm. However, even with modest gaps between detections and a modest movement capacity, this may suggest a very large area, and it may take a long time to compute the requisite number of centroids using \code{\link[flapper]{acdc_setup_centroids}}. Moreover, in practice, the area of interest may be smaller.}
#'
#' \subsection{Method (2)}{The second method provides a reasonable lower bound for the number of centroids. This method is based on the locations of receivers and the assumption that, when an individual is detected by two different receivers, at the halfway point between detections, its potential location is described by the intersection of the acoustic centroid around the receiver at which it was previously detected and the centroid around the receiver at which it is next detected (evaluated at the halfway point between detections). Under this perspective, the minimum centroid size is half of the distance between the furthest two receivers that were operational at the same time and the minimum centroid number is this size over the animal's mobility†. (This could be restricted to the subset of receivers at which the individual was detected sequentially but this is not implemented.) In practice, it is advisable to double this minimum number to ensure a reasonable degree of overlap between the centroids of the two receivers. To implement this method, a dataframe that contains receiver locations and deployment periods (\code{moorings}) must be supplied, along with the \code{mobility} parameter that describes how far an individual could move in any archival time step. This method provides a sensible lower bound for the number of centroids that is array-specific. However, it does not account for the full span of movements that are possible under the AC* algorithm(s).}
#'
#' \subsection{Method (3)}{The final method provides an alternative suggestion for the number of centroids based on the minimum number that is required for the largest centroid around each receiver to fill the study area. This method requires the locations of receivers and a Spatial* or Raster* layer that defines the study area, the boundaries of which are approximated as a bounding box with four vertices and transformed to degrees latitude/longitude. For each receiver, the function calculates the distance to each of the four boundary coordinates. The longest distance between a receiver and a boundary of the study area is used to define the number of centroids required for the largest centroid to fill the study area†. In many settings, this method may be most suitable because it accounts for the full span of movements that are possible under the AC* algorithm(s) while focusing on a specific area.}
#'
#' \subsection{Other considerations}{In practice, the most appropriate number of centroids is likely to be a compromise between these minimum and maximum values that depends on other considerations, particularly the spatial scale of the study, the effect of boundaries on the final centroid size and, perhaps in some cases, the geographical range of the species. For acoustic time series that are followed by a long tail of archival observations after final acoustic detection, both methods may underestimate the number of centroids required to capture the increasing uncertainty in the individual's location over this time. However, in this case, it is advisable to use the depth-contour (DC) algorithm (see \code{\link[flapper]{dc}}) at some point after the last acoustic detection since the influence of that observation on putative patterns of space use will decay through time and the DC algorithm is more computationally efficient. If in doubt, opt to make more, rather than fewer centroids: if there are too few centroids, the AC* algorithm(s) can fail to converge if centroids around receivers with sequential detections do not overlap.}
#'
#' †Both Methods Two and Three assume that the detection range is the same as the \code{mobility}. In practice, this is unlikely to be the case, but if the detection range is even approximately equal to \code{mobility} then the difference is negligible.
#'
#' @return The function prints a table with the suggested number of centroids and the maximum radius size from each method. Invisibly, the function returns an integer vector with three suggested values for the number of centroids from methods (1), (2) and (3). The parameters used to generate these suggestions (i.e., \code{detections}, \code{moorings}, \code{mobility}, \code{double} and \code{boundaries}) are also included in a `param' attribute.
#'
#' @seealso This function is designed to facilitate an informed choice for the `n_timesteps' argument in \code{\link[flapper]{acdc_setup_centroids}}. \code{\link[flapper]{acdc_setup_mobility}} can also guide the implementation of this function. This underpins the AC and ACDC algorithms, which are implemented by \code{\link[flapper]{ac}} and \code{\link[flapper]{acdc}}.
#'
#' @examples
#' n_timesteps <-
#'   acdc_setup_n_centroids(detections = dat_acoustics$timestamp[dat_acoustics$individual_id == 25],
#'                          step = 120,
#'                          moorings = dat_moorings,
#'                          mobility = 200,
#'                          double = TRUE,
#'                          boundaries = dat_gebco)
#' utils::str(n_timesteps)
#'
#' @author Edward Lavender
#' @export

acdc_setup_n_centroids <- function(detections,
                                   step,
                                   hist = TRUE,...,
                                   moorings,
                                   mobility,
                                   double = TRUE,
                                   boundaries){

  #### Checks
  # verbose = TRUE
  # t_onset <- Sys.time()
  # cat_to_console <- function(..., show = verbose) if(show) cat(paste(..., "\n"))
  # cat_to_console(paste0("flapper::acdc_setup_n_centroids() called (@ ", t_onset, ")..."))
  # Check moorings contains required information
  check_names(input = moorings,
              req = c("receiver_id",
                      "receiver_start_date", "receiver_end_date",
                      "receiver_lat", "receiver_long"),
              extract_names = colnames,
              type = all)

  #### Method (1): Determine the maximum gap between detections
  # cat_to_console("... Implementing Mpproach (1)...")
  # Calculate the duration of gaps
  gaps <- Tools4ETS::serial_difference(detections, units = "s")
  gaps <- as.numeric(gaps)
  max_gap <- max(gaps, na.rm = TRUE)
  # Visualise the distribution of gaps
  if(hist) {
    prettyGraphics::pretty_hist(gaps,...)
    graphics::abline(v = max_gap, col = "red", lty = 3)
  }
  # The minimum number of time steps under this method is max_gap/2/step
  minimum_n_timesteps_1 <- ceiling(max_gap/2/step)

  #### Method (2): Determine the number of timesteps to overlap centroids between receivers
  # cat_to_console("... Implementing Method (2)...")

  ## Identify time intervals over which receivers were deployed
  moorings$interval <- lubridate::interval(moorings$receiver_start_date, moorings$receiver_end_date)

  ## Define the minimum number of timesteps for each combination of receivers
  # ... that was deployed for an overlapping interval
  minimum_n_timesteps_2 <-
    sapply(unique(moorings$interval), function(interval){
      # Identify receivers whose deployment periods overlapped with this interval
      moorings_tmp <- moorings
      moorings_tmp$overlap <- lubridate::int_overlaps(interval, moorings_tmp$interval)
      moorings_tmp <- moorings_tmp[which(moorings_tmp$overlap), ]
      dist_btw_receivers_m <- dist_btw_receivers(moorings_tmp, f = function(x) x*1000)
      # Calculate the minimum number of timesteps for intervals to overlap
      minimum_n_timesteps <- max(dist_btw_receivers_m$dist/mobility)
      # Double the overlap so that the centroids would fully overlap
      if(double) minimum_n_timesteps <- minimum_n_timesteps*2
      return(minimum_n_timesteps)
    })
  # Across all receivers deployment combinations, calculate the minimum overlap
  # ... which the the maximum of the minimum overlaps across all receivers
  minimum_n_timesteps_2 <- max(minimum_n_timesteps_2)

  #### Method (3): The number of time steps for the centroids to fill the area
  # cat_to_console("... Implementing Method (3)...")
  rxy <- as.matrix(moorings[, c("receiver_long", "receiver_lat")])
  ext <- raster::extent(boundaries)
  if(is.na(raster::crs(boundaries))) stop("The 'boundaries' CRS is NA.")
  ext <- sp::SpatialPoints(ext, proj4string = raster::crs(boundaries))
  ext <- sp::spTransform(ext, sp::CRS("+proj=longlat +datum=WGS84 +no_defs"))
  dist_from_rxy_to_boundaries <- raster::pointDistance(ext, rxy, lonlat = TRUE)
  minimum_n_timesteps_3 <- ceiling(max(dist_from_rxy_to_boundaries)/mobility)

  #### Return suggestions
  # Collect suggestions
  minimum_n_timesteps <- c(minimum_n_timesteps_1, minimum_n_timesteps_2, minimum_n_timesteps_3)
  minimum_n_timesteps <- as.integer(ceiling(minimum_n_timesteps))
  attributes(minimum_n_timesteps)$method <- 1:3
  attributes(minimum_n_timesteps)$param <- list(detections = detections,
                                                moorings = moorings,
                                                mobility = mobility,
                                                double = double)
  # Print results
  cat("Results -------------------------------\n")
  results <- data.frame(method = c(1, 2, 3),
                        n_centroids = minimum_n_timesteps,
                        radius = minimum_n_timesteps * mobility)
  rownames(results) <- NULL
  print(results)
  cat("---------------------------------------\n")
  # Close function and return outputs
  t_end <- Sys.time()
  # duration <- difftime(t_end, t_onset, units = "mins")
  # cat_to_console(paste0("... flapper::acdc_setup_n_centroids() call completed (@ ", t_end, ") after ~", round(duration, digits = 2), " minutes."))
  return(invisible(minimum_n_timesteps))
}


######################################
######################################
#### acdc_setup_centroids()

#' @title Setup the acoustic centroids required for the AC* algorithm(s)
#' @description This function produces the acoustic centroids required by the acoustic-centroid (AC) and acoustic-centroid depth-contour (ACDC) algorithms.
#' @param xy A \code{\link[sp]{SpatialPointsDataFrame}} object that defines receiver IDs and locations. The \code{data} slot must include a dataframe with a column that provides a unique, integer identifier for each receiver (`receiver_id'). The coordinate reference system should be the Universal Transverse Mercator system with distances in metres (to match \code{detection_range}, see below).
#' @param detection_range A number that defines the maximum detection range (m) at which an individual could be detected from a receiver.
#' @param mobility A number that defines the distance that an individual could move in the time period between archival observations.
#' @param n_timesteps An integer that defines the number of timesteps after a hypothetical detection for which centroids will be created, where the duration of each timestep is given by the duration between archival observations (see \code{\link[flapper]{acdc_setup_n_centroids}}).
#' @param coastline (optional) A \code{\link[sp]{SpatialPolygonsDataFrame-class}} object that defines the coastline in an area. If provided, acoustic centroids are processed to remove any areas on land. Algorithm speed declines with the complexity of the coastline.
#' @param boundaries (optional) An \code{\link[raster]{extent}} object that defines the boundaries of an area within which individuals are assumed to have remained. If provided, acoustic centroids are processed to remain within this area.
#' @param plot A logical input that defines whether or not to produce a plot of the area, including receivers, the coastline and the area boundaries (if provided), and acoustic centroids. This is useful for checking purposes but it can reduce algorithm speed.
#' @param cl (optional) A cluster object created by \code{\link[parallel]{makeCluster}} to implement the algorithm in parallel. The connection to the cluster is closed within the function.
#' @param verbose A logical input that defines whether or not to print messages to the console to relay function progress.
#'
#' @details Given only a detection at a particular receiver at a particular time, the detection range of the receiver and the movement speed of the animal, an acoustic centroid defines the possible locations of a detected individual at that time or a subsequent time. More specifically, when an individual is located at a receiver, its location must be within some radius of that receiver defined by the maximum detection distance. This radius expands with the duration since detection. This function defines, for each receiver, a list of acoustic centroids that reflect the possible locations of an individual were it to have been detected from 0 to \code{n_timesteps} ago at that receiver, accounting for the coastline and within a defined area if necessary. Using the observed detection data, the AC (\code{\link[flapper]{ac}}) and ACDC (\code{\link[flapper]{acdc}}) algorithms pull the relevant centroids out of this list, which substantially saves computation time because acoustic centroids are not computed on-the-fly. These centroids are processed within these algorithms (e.g., if an individual is detected at two different receivers, then at the halfway point between detections its location must be within the intersection of the relevant centroids for these two receivers) and combine them with other information to reconstruct where an individual could have been over time.
#'
#' @return The function returns a list of \code{\link[sp]{SpatialPolygonsDataFrame-class}} objects, with one element for all numbers from 1 to the maximum receiver number (\code{xy$receiver_id}). Any list elements that do not correspond to receivers contain a \code{NULL} element. List elements that correspond to receivers contain a \code{\link[sp]{SpatialPolygonsDataFrame-class}} object containing all the centroids for that receiver.
#'
#' @examples
#' #### Define data for acdc_setup_centroids()
#' ## Define coordinates of receivers as SpatialPointsDataFrame with UTM CRS
#' # CRS of receiver locations as recorded in dat_moorings
#' proj_wgs84 <- sp::CRS("+init=epsg:4326")
#' # CRS of receiver locations required
#' proj_utm <- sp::CRS(paste("+proj=utm +zone=29 +datum=WGS84",
#'                           "+units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))
#' # Define SpatialPoints object
#' xy_wgs84 <- sp::SpatialPoints(dat_moorings[, c("receiver_long", "receiver_lat")], proj_wgs84)
#' xy_utm <- sp::spTransform(xy_wgs84, proj_utm)
#' # Link with receiver IDs to define a SpatialPointsDataFrame
#' xy_utm <-
#'  sp::SpatialPointsDataFrame(xy_utm,
#'                             dat_moorings[, "receiver_id", drop = FALSE])
#'
#' #### Example (1): Define a list of centroids with specified parameters
#' # ... (Argument values are small to reduce computation time for examples)
#' centroids <- acdc_setup_centroids(xy = xy_utm,
#'                                   detection_range = 500,
#'                                   mobility = 250,
#'                                   n_timesteps = 3
#'                                  )
#' # A list of SpatialPolygonsDataFrames is returned
#' # with elements from 1:max(xy_utm$receiver_id)
#' # NULL elements correspond to numbers in this sequence that do not refer to receivers
#' # Otherwise a SpatialPolygonsDataFrame is returned with all the centroids for that receiver
#' centroids
#'
#' #### Example (2): Visualise the acoustic centroids produced via plot = TRUE
#' centroids <- acdc_setup_centroids(xy = xy_utm,
#'                                   detection_range = 500,
#'                                   mobility = 250,
#'                                   n_timesteps = 3,
#'                                   plot = TRUE
#'                                   )
#'
#' #### Example (3): Remove areas of the centroids that overlap with coastline
#' centroids <- acdc_setup_centroids(xy = xy_utm,
#'                                   detection_range = 500,
#'                                   mobility = 250,
#'                                   n_timesteps = 3,
#'                                   plot = TRUE,
#'                                   coastline = dat_coast
#'                                   )
#'
#' #### Example (4): Remove areas of the centroids beyond a boundary
#' xy_utm_coords <- sp::coordinates(xy_utm)
#' boundaries <- raster::extent(min(xy_utm_coords[, 1]),
#'                              max(xy_utm_coords[, 1]),
#'                              min(xy_utm_coords[, 2]),
#'                              max(xy_utm_coords[, 2])
#'                         )
#' centroids <- acdc_setup_centroids(xy = xy_utm,
#'                                   detection_range = 500,
#'                                   mobility = 250,
#'                                   n_timesteps = 3,
#'                                   plot = TRUE,
#'                                   coastline = dat_coast,
#'                                   boundaries = boundaries
#'                                   )
#'
#' #### Example (5): Implement the algorithm in parallel
#' centroids <- acdc_setup_centroids(xy = xy_utm,
#'                                   detection_range = 500,
#'                                   mobility = 250,
#'                                   n_timesteps = 3,
#'                                   plot = TRUE,
#'                                   coastline = dat_coast,
#'                                   boundaries = boundaries,
#'                                   cl = parallel::makeCluster(2L)
#'                                   )
#'
#' #### Example (6): Acoustic centroids can be saved to file using rlist::list.save()
#' # rlist::list.save(centroids, paste0(tempdir(), "/centroids.RData"))

#' @author Edward Lavender
#' @export
#'

acdc_setup_centroids <- function(
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
  if(!inherits(xy, "SpatialPointsDataFrame"))
    stop(paste("Argument 'xy' must be of class 'SpatialPointsDataFrame', not class(es):"), class(xy))
  rs <- xy$receiver_id
  if(is.numeric(rs)) rs <- as.integer(rs)
  if(!is.integer(rs))
    stop(paste("Argument 'xy$receiver_id' must be of class 'integer', not class(es):"), class(rs))
  if(any(rs <= 0))
    stop("Argument 'xy$receiver_id' cannot contain receiver IDs <= 0.")
  if(any(duplicated(rs))){
    message("Argument 'xy$receiver_id' contains duplicate elements. 'xy' has been simplified to contain only unique receiver IDs.")
    xy <- xy[!duplicated(xy$receiver_id), ]
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

  ## Define a sequence of centroid radius
  # Around each receiver, we'll create a polygon of this radius
  radius_seq <- seq(detection_range, length.out = n_timesteps, by = mobility)

  ## Define a list of receivers, with a list of centroids for each receiver
  bathy_ls <- pbapply::pblapply(1:length(rs), cl = cl, function(i){

      centroids_ls <- lapply(radius_seq, function(radius){

        # Define a buffer around the current receiver of appropriate radius
        bathy_poly <- rgeos::gBuffer(xy[i, ], width = radius)

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
#### acdc_setup_detection_kernels()

#' @title Setup detection probability kernels for the AC* algorithm(s)
#' @description This function produces detection probability kernels for incorporation into the acoustic-contour* (AC*) algorithms. Within acoustic centroids, the incorporation of detection probability kernels reduces uncertainty and increases precision in simulated patterns of space use by up-weighting areas nearer to receivers when an individual is detected and down-weighing areas nearer to receivers when an individual is not detected (see Details).
#'
#' To implement the function, a \code{\link[sp]{SpatialPointsDataFrame}} that defines receiver IDs, locations and deployment dates must be supplied (via \code{xy}). A record of servicing events for receivers can also be supplied (via \code{services}). Detection probability kernels are calculated around each receiver, using a user-defined function based on Euclidean distances (\code{calc_detection_pr}) across a \code{\link[raster]{raster}} (\code{map}). Kernels can be restricted by the \code{coastline} and \code{boundaries} of an area if applicable. These kernels are used to weight possible locations around a receiver when an individual is detected.
#'
#' For each unique array design (i.e. set of active receivers, given receiver deployment dates and servicing events, if applicable), a detection probability surface across the whole area is also created, which is used to weight possible locations of the individual in the time steps between detections (up-weighting locations away from receivers). By default, these calculations account for any areas of overlap in the detection probability kernels of multiple receivers. This step is computationally demanding, but it can be suppressed or sped-up via the \code{overlaps} argument. Outputs are returned in a named list that is designed to be incorporated into the AC* algorithm(s).
#'
#' @param xy A \code{\link[sp]{SpatialPointsDataFrame}} that defines receiver IDs, locations and deployment dates. The \code{data} slot must include a dataframe with the following columns: an unique, integer identifier for each receiver (`receiver_id') and receiver deployment \code{\link[base]{Dates}} (`receiver_start_date' and `receiver_end_date'). For receiver locations, the coordinate reference system should be the Universal Transverse Mercator system with distances in metres (as for \code{map}, below).
#' @param services (optional) A dataframe that defines receiver IDs and servicing \code{\link[base]{Dates}} (times during the deployment period of a receiver when it was not active due to servicing). If provided, this must contain the following columns: an integer identifier for serviced receivers (named ‘receiver_id’) and two columns that define the time of the service(s) (‘service_start_date’ and ‘service_end_date’) (see \code{\link[flapper]{make_matrix_receivers}}).
#' @param centroids The list of acoustic centroids, with one element for each number from \code{1:max(xy$receiver_id)}, from \code{\link[flapper]{acdc_setup_centroids}}.
#' @param overlaps (optional) A named list, from \code{\link[flapper]{get_detection_centroids_overlap}}, that defines, for each receiver, for each day over its deployment period, whether or not its detection centroid overlapped with those of other receivers. If provided, this speeds up detection probability calculations in overlapping regions by focusing on only the subset of receivers with overlapping detection probability kernels. If there are no overlapping receivers, \code{FALSE} can be supplied instead to suppress these calculations.
#' @param calc_detection_pr A function that takes in a vector of distances and returns a vector of detection probabilities (around a receiver). Detection probability should decline to 0 after the \code{detection_range} distance from a receiver (see \code{\link[flapper]{acdc_setup_centroids}}).
#' @param map A (blank) \code{\link[raster]{raster}} on which receiver locations and detection probability kernels are represented. As for \code{xy}, the coordinate reference system should be the Universal Transverse Mercator system with distances in metres. Resolution needs to be sufficiently high such that detection probability can be represented effectively, given \code{calc_detection_pr}, but sufficiently low for convenient run times. For the \code{\link[flapper]{acdc}} algorithm, this should have the same properties as \code{bathy}.
#' @param coastline (optional) A \code{\link[sp]{SpatialPolygonsDataFrame-class}} mask applied to detection probability kernels (e.g. on land, as in  \code{\link[flapper]{acdc_setup_centroids}}). Algorithm speed declines with the complexity of \code{coastline}.
#' @param boundaries An \code{\link[raster]{extent}} object that defines the boundaries of an area within which individuals are assumed to have remained (as in \code{\link[flapper]{acdc_setup_centroids}}), if these are different from the extent of \code{map}. If provided, detection probability kernels are processed to remain within this area.
#' @param verbose A logical input that defines whether or not to print messages to the console to relay function progress.
#' @param ... Additional arguments (none implemented).
#'
#' @details A detection probability kernel is a bivariate probability density function that describes detection probability around a receiver. Typically, this takes the shape of a dome whereby detection probability declines uniformly around a receiver with increasing distance. Accordingly, this function assumes that detection probability kernels only depend on distance (via \code{calc_detection_pr}) and are constant in space and time. Spatially variable detection probability kernels can be incorporated just as easily into the AC* algorithm(s) but need to be created manually. However, temporally variable detection probability kernels are not currently supported by the AC* algorithm(s).
#'
#' The purpose of detection probability kernels within the AC* algorithm(s) is to reduce the uncertainty of the possible positions of an individual within the acoustic centroids, both when an individual is detected and when it is not. When an individual is detected at a receiver, a central assumption of the AC* algorithm(s) is that the individual is within a finite radius of that receiver in an area defined as the ‘detection centroid’. Under the simplest implementation of the algorithm, this represents a threshold detection probability model in which the probability of detection is certain within the centroid and zero outside. Under this model, the individual is equally likely to have been in any location within the centroid and all possible locations of the individual are weighted equally. This approach is suitable for the most conservative analyses, since even unlikely positions of the individual receive equal weighting. However, detection probability typically declines with distance around a receiver and the incorporation of this information in the form of kernels around the receiver at which is detected (and any receivers with overlapping kernels) improve precision by increasing the weighting for some areas over others. The way in which this increase in precision is realised over space depends on whether or not detection probability kernels and the timing of detections at multiple receivers overlap.
#'
#' In the simplest scenario, an individual is detected at a receiver whose centroid does not overlap with any other receiver and the kernel simply increases the weight of locations nearest to the receiver (according to a user-specified detection probability function). For some array designs, an alternative possibility is that an individual is detected at a receiver whose centroid overlaps with other receiver(s). In this case, detection or lack of detection at the same moment at the receivers with overlapping centroids provides further information on the location of the individual. One the one hand, if the individual is not detected at the other receiver(s), then the probability that the individual is in the overlapping region(s) is reduced in line with the overlap in detection probability, which decreases the weight of potential locations in these areas. On the other hand, if the individual is detected at effectively the same time at multiple receiver(s), then the individual is more likely to be within the overlapping parts of their centroids, and the weight of possible locations here is correspondingly elevated. The definition of ‘effectively’ is context specific but depends on the accuracy with which the clocks of different receivers are synchronised. (For example, ± 5 s might be reasonable.) In these situations, the detection centroid essentially just acts as a computational device by restricting calculations to the centroid(s) and ‘cutting’ the probability of detection to zero beyond this area.
#'
#' This process follows the standard rules of probability, with the probability of any location (xy) when the individual is detected given by the product of the detection probability at all the receivers at which the individual was detected (\eqn{\prod Pr(det)_{xy, detection = 1}}, which up-weights areas that intersect between receivers that recorded detections), multiplied by the product of not being detected at all the receivers at which it was not detected (\eqn{\prod (1 - Pr(det)_{xy, detection = 0}}, which down-weights areas that intersect with receivers that did not record detections); i.e.,
#'
#' \deqn{\prod Pr(det)_{xy, detection = 1} \times \prod (1 - Pr(det)_{xy, detection = 0}}
#'
#' By way of illustration, consider a simple array comprising three equidistant receivers with equally overlapping detection centroids. An individual is detected at two receivers but not the third. In this scenario, the individual must be located in the intersection between the centroids of the two receivers at which it was detected, but it is more likely to be located in the part of this region that does not intersect with the third receiver (call this area A and the intersecting area for all three receivers B). To assign some numbers to this example, consider the overlap of a single detection probability, say the contour \eqn{Pr(det) = 0.2}. In this case, the probability of the individual being located in area A is the probability of being detected at receiver 1 (0.2) and receiver 2 (0.2) but not receiver (3) (\eqn{1 - 0}), which equals 0.04. In comparison, the probability of the individual being located in area B is the probability of being detected at receiver 1 (0.2) and receiver 2 (0.2), but not receiver 3 (\eqn{1 - 0.2}), which equals 0.032. (In the AC* algorithms(s), these probabilities are re-scaled, but the point is detection probability kernels up-weight some areas and down-weight others, using the rules of probability, in line with intuitive expectations.)
#'
#' When an individual is not detected, the detection centroid grows into a set of ‘acoustic centroids’ that describe our increasing uncertainty in the location of the individual. As they grow, they may encompass other detection centroids before they shrink towards the receiver at which the individual is next detected. During this time, the AC* algorithm(s) identify the possible locations of the individual within these areas. Under the most conservative approach, at each time step all positions are treated as equally likely (although normalisation within the algorithm(s) can be implemented to down-weight time steps in which the location of the individual is more uncertain). This includes any positions within the detection centroids of other receivers since, under realistic conditions, there is usually a non-zero probability that an individual can be near a receiver and yet remain undetected. However, this is typically unlikely and when the goal of the analysis is to create more precise estimates of space use, the incorporation of detection probability kernel(s) around receivers effectively reduces the probability that the individual is within their detection centroids during this time. Again, to incorporate detection probability kernels in this way, it is necessary to account for overlapping detection ranges, where the probability of detection is higher and, therefore, the probability of a possible location in such an area is lower given the absence of a detection. As above, this process follows the laws of probability. In any given location, the probability of no detection(s) is given by the product of receivers' inverse detection probability kernels. This is used to down-weight possible locations near receivers, especially in areas with overlapping kernels, effectively up-weighting possible locations further away, when an individual is not detected.
#'
#' In summary, in the AC* algorithm(s), detection and acoustic centroids describe the spatial extent of our uncertainty when an individual is detected and in the time between detections respectively. The purpose of this function is to pre-process and package the information provided by detection probability kernels in such a way as to facilitate its incorporation into the AC* algorithm(s). This improves the precision of simulated patterns of space use by down- or up-weighting areas according to a model of detection probability. This provides a reasonable overall assessment of the places in which an individual could have spend more or less time over a period of study.
#'
#' However, it is worth noting that, within acoustic centroids, unrealistic areas may be highlighted in which an individual could not have been located because of movement constraints which are not captured by centroid-level expansion and contraction. For example, there may be parts of a centroid that are highly unlikely given a detection at a nearby receiver because they are too far away from the receiver at which the individual was subsequently detected. In the AC*PF algorithm(s), the incorporation of a movement model further reduces uncertainty by accounting for the constraints imposed by an individual’s current position on its next position.
#'
#' @return The function returns a named list with five elements:
#' \enumerate{
#'   \item \strong{receiver_specific_kernels.} A list, with one element for all integers from 1 to the maximum receiver number. Any elements that do not correspond to receivers contain a \code{NULL} element. List elements that correspond to receivers contain a \code{\link[raster]{raster}} of the detection probability kernel around the relevant receiver. Cells values define the detection probability around a receiver, given \code{calc_detection_pr}. In the AC* algorithm(s), these kernels are used to up-weight location probabilities near to a receiver when it is detected (following modification to account for overlapping areas, if necessary).
#'   \item \strong{receiver_specific_inv_kernels.} A list, as for \code{receiver_specific_kernels}, but in which elements contain the inverse detection probability kernels (i.e., 1 - detection probability). In the AC* algorithm(s), these is used to down-weight-weight location probabilities in the overlapping regions between a receiver that recorded detections and others that did not at the same time.
#'   \item \strong{array_design_intervals.} A dataframe that defines the number and deployment times of each unique array design, resulting from receiver deployment, servicing and removal. In the times between detections, this is used to select the appropriate `background' detection probability surface (see below). This contains the following columns:
#'     \itemize{
#'       \item \code{array_id} An integer vector that defines each unique array design.
#'       \item \code{array_start_date} A Date that defines the start date of each array design.
#'       \item \code{array_end_date} A Date that defines the end date of each array design.
#'       \item \code{array_interval} An \code{\link[lubridate]{Interval-class}} vector that defines the deployment period of each array design.
#'     }
#'   \item \strong{bkg_surface_by_design.} A list, with one element for each array design, that defines the detection probability surface across all receivers deployed in that phase of the study. In areas that are covered by the detection probability kernel of a single receiver, the detection probability depends only on distance to that receiver (via \code{calc_detection_pr}). In areas covered by multiple, overlapping kernels, detection probability represents the combined detection probability across all overlapping kernels (see Details).
#'   \item \strong{bkg_inv_surface_by_design.} A list, as above for \code{bkg_surface_by_design}, but which contains the inverse detection probability surface (i.e., 1 - \code{bkg_surface_by_design}). In the AC* algorithm(s), this is used to up-weight areas away from receivers (or, equivalently, down-weight areas near to receivers) in the time steps between detections.
#' }
#'
#' @examples
#' #### Set up data for examples
#'
#' ## Define receiver IDs, locations and deployment dates
#' # Focus on a subset of receivers for example speed
#' moorings <- dat_moorings[1:5, ]
#' # Define receiver locations as a SpatialPoints object with a UTM CRS
#' proj_wgs84 <- sp::CRS("+init=epsg:4326")
#' proj_utm <- sp::CRS(paste("+proj=utm +zone=29 +datum=WGS84",
#'                           "+units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))
#' xy <- sp::SpatialPoints(moorings[, c("receiver_long", "receiver_lat")],
#'                         proj_wgs84)
#' xy <- sp::spTransform(xy, proj_utm)
#' # Link receiver IDs, locations and deployment dates to form a SpatialPointsDataFrame
#' # ... Note required column names and class types.
#' xy <- sp::SpatialPointsDataFrame(xy, moorings[, c("receiver_id",
#'                                                   "receiver_start_date",
#'                                                   "receiver_end_date")])
#'
#' ## Define blank map of area for which AC* algorithm(s) will be implemented
#' # ... The resolution must be sufficiently high
#' # ... such that there are areas with non zero detection probability.
#' # ... However, function speed will fall with large, high resolution rasters.
#' # ... Here, we set a low resolution for example speed.
#' map_blank <- raster::raster(raster::extent(dat_gebco), res = c(50, 50))
#' map_blank <- raster::setValues(map_blank, 0)
#'
#' ## Define detection probability function
#' # This should depend on distance alone
#' # Detection probability should decline to 0 after detection_range
#' # ... (defined in flapper::acdc_setup_centroids()).
#' # ... Here, we assume detection_range = 425 m.
#' calc_dpr <-
#'   function(x){
#'     ifelse(x <= 425, stats::plogis(2.5 + -0.02 * x), 0)
#'   }
#' plot(0:1000, calc_dpr(0:1000), type = "l")
#'
#' ## Get detection centroids and, if applicable, information on their overlap(s)
#' # We'll use the example centroids provided in dat_centroids
#' # We'll get their overlaps via flapper::get_detection_centroids_overlap
#' overlaps <-
#'   get_detection_centroids_overlap(
#'     centroids = get_detection_centroids(xy = xy, byid = TRUE)
#'   )
#'
#' #### Example (1): Implement function using default options
#' kernels <- acdc_setup_detection_kernels(xy = xy,
#'                                         centroids = dat_centroids,
#'                                         overlaps = overlaps,
#'                                         calc_detection_pr = calc_dpr,
#'                                         map = map_blank)
#'
#' # Examine list elements
#' summary(kernels)
#'
#' # Examine example receiver-specific kernels
#' pp <- graphics::par(mfrow = c(1, 2))
#' raster::plot(kernels$receiver_specific_kernels[[3]])
#' points(xy[xy$receiver_id == 3, ], cex = 2)
#' raster::plot(kernels$receiver_specific_kernels[[4]])
#' points(xy[xy$receiver_id == 4, ], cex = 2)
#' graphics::par(pp)
#'
#' # Examine example receiver-specific inverse kernels
#' pp <- graphics::par(mfrow = c(1, 2))
#' raster::plot(kernels$receiver_specific_inv_kernels[[3]])
#' points(xy[xy$receiver_id == 3, ], cex = 2)
#' raster::plot(kernels$receiver_specific_inv_kernels[[4]])
#' points(xy[xy$receiver_id == 4, ], cex = 2)
#' graphics::par(pp)
#'
#' # Examine background detection Pr surfaces
#' # (for each unique combination of receivers that were deployed)
#' pp <- graphics::par(mfrow = c(2, 3), mar = c(0, 0, 0, 0))
#' lapply(kernels$bkg_surface_by_design, function(bkg) {
#'   raster::plot(bkg, axes = FALSE)
#'   graphics::box()
#' })
#' graphics::par(pp)
#'
#' # Examine background inverse detection Pr surfaces
#' pp <- graphics::par(mfrow = c(2, 3), mar = c(0, 0, 0, 0))
#' lapply(kernels$bkg_inv_surface_by_design, function(bkg) {
#'   raster::plot(bkg, axes = FALSE)
#'   graphics::box()
#' })
#' graphics::par(pp)
#'
#' #### Example (2): Incorporate spatial information (coastline, area boundaries)
#'
#' # Define spatial mask for the land and check this via a plot
#' sea <- invert_poly(dat_coast)
#' area <- raster::setValues(map_blank, 1)
#' raster::plot(raster::mask(area, sea))
#'
#' # Implement algorithm
#' kernels <- acdc_setup_detection_kernels(xy = xy,
#'                                         centroids = dat_centroids,
#'                                         overlaps = overlaps,
#'                                         calc_detection_pr = calc_dpr,
#'                                         map = map_blank,
#'                                         coastline = sea,
#'                                         boundaries = update_extent(map_blank, -1000))
#'
#' # Examine example receiver-specific kernels
#' pp <- graphics::par(mfrow = c(1, 2))
#' raster::plot(kernels$receiver_specific_kernels[[3]])
#' points(xy[xy$receiver_id == 3, ], cex = 2)
#' raster::plot(kernels$receiver_specific_kernels[[4]])
#' points(xy[xy$receiver_id == 4, ], cex = 2)
#' graphics::par(pp)
#'
#' # Examine example receiver-specific inverse kernels
#' pp <- graphics::par(mfrow = c(1, 2))
#' raster::plot(kernels$receiver_specific_inv_kernels[[3]])
#' points(xy[xy$receiver_id == 3, ], cex = 2)
#' raster::plot(kernels$receiver_specific_inv_kernels[[4]])
#' points(xy[xy$receiver_id == 4, ], cex = 2)
#' graphics::par(pp)
#'
#' # Examine background detection Pr surfaces
#' # ... (for each unique combination of receivers that were deployed)
#' pp <- graphics::par(mfrow = c(2, 3), mar = c(0, 0, 0, 0))
#' lapply(kernels$bkg_surface_by_design, function(bkg) {
#'   raster::plot(bkg, axes = FALSE)
#'   graphics::box()
#' })
#' graphics::par(pp)
#'
#' # Examine background inverse detection Pr surfaces
#' pp <- graphics::par(mfrow = c(2, 3), mar = c(0, 0, 0, 0))
#' lapply(kernels$bkg_inv_surface_by_design, function(bkg) {
#'   raster::plot(bkg, axes = FALSE)
#'   graphics::box()
#' })
#' graphics::par(pp)
#'
#' @seealso This is one of a number of functions used to set up the AC and ACDC algorithms implemented by \code{\link[flapper]{ac}} and \code{\link[flapper]{acdc}}: \code{\link[flapper]{acdc_setup_mobility}}, \code{\link[flapper]{acdc_setup_n_centroids}}, \code{\link[flapper]{acdc_setup_centroids}} and \code{\link[flapper]{acdc_setup_detection_kernels}}. This function is supported by \code{\link[flapper]{make_matrix_receivers}}, which defines receiver activity statuses; \code{\link[flapper]{acdc_setup_centroids}}, which defines acoustic centroids; and \code{\link[flapper]{get_detection_centroids_overlap}}, which defines detection centroid overlaps
#' @author Edward Lavender
#' @export
#'

acdc_setup_detection_kernels <-
  function(xy, services = NULL,
           centroids, overlaps = NULL,
           calc_detection_pr,
           map,
           coastline = NULL, boundaries = NULL,
           verbose = TRUE,...
  ){


    ######################################
    #### Initiation, set up and checks

    #### Initiation
    cat_to_console <- function(..., show = verbose){
      if(show) cat(paste(..., "\n"))
    }
    t_onset <- Sys.time()
    cat_to_console(paste0("flapper::acdc_setup_detection_kernels() called (@ ", t_onset, ")..."))
    cat_to_console("... Setting up function...")

    #### Checks
    ## xy
    if(!inherits(xy, "SpatialPointsDataFrame")) stop("'xy' must be a SpatialPointsDataFrame.")
    ## moorings
    moorings <- data.frame(xy)
    check_names(input = moorings, req = c("receiver_id", "receiver_start_date", "receiver_end_date"),
                extract_names = colnames, type = all)
    if(is.numeric(moorings$receiver_id)) moorings$receiver_id <- as.integer(moorings$receiver_id)
    if(!is.integer(moorings$receiver_id))
      stop(paste("Argument 'xy$receiver_id' must be of class 'integer', not class(es):"), class(moorings$receiver_id))
    if(any(moorings$receiver_id <= 0))
      stop("Argument 'xy$receiver_id' cannot contain receiver IDs <= 0.")
    if(any(duplicated(moorings$receiver_id )))
      stop("Argument 'xy$receiver_id' contains duplicate elements.")
    ## services
    if(!is.null(services)){
      check_names(input = services, req = c("receiver_id", "service_start_date", "service_end_date"),
                  extract_names = colnames, type = all)
      if(is.numeric(services$receiver_id)) services$receiver_id <- as.integer(services$receiver_id)
      if(!is.integer(services$receiver_id))
        stop(paste("Argument 'services$receiver_id' must be of class 'integer', not class(es):"), class(services$receiver_id))
      if(!all(unique(services$receiver_id) %in% unique(moorings$receiver_id))){
        message("Not all receivers in services$receiver_id are in moorings$receiver_id.")
      }
    }
    ## overlaps
    if(is.logical(overlaps)) {
      if(overlaps) stop("'overlaps' must be a named list or FALSE.")
    } else {
      check_named_list(input = overlaps, ignore_empty = FALSE)
      check_names(input = overlaps, req = "list_by_receiver", type = all)
    }
    ## spatial data
    # check spatial data
    if(!is.null(coastline) & !is.null(boundaries)) {
      coastline <- raster::crop(coastline, boundaries)
      if(is.null(coastline)) message("No coastline within defined boundaries. \n")
    }
    # process map/area
    map <- raster::setValues(map, 0)
    if(!is.null(boundaries)) map <- raster::crop(map, boundaries)
    if(!is.null(coastline)) map <- raster::mask(map, coastline)
    # pull out detection centroids (for loop speed)
    centroids <- lapply(centroids, function(elm) if(!is.null(elm)) elm[1, ])
    # cluster: not currently implemented due to issues with raster temporary files
    # if(is.null(cl) & !is.null(varlist)) {
    #   warning("cl = NULL but varlist is not NULL: varlist ignored.",
    #           immediate. = TRUE, call. = TRUE)
    #   varlist <- NULL
    # }


    ######################################
    ####  Receiver-specific kernels (for detection)

    #### Calculate detection Pr around each receiver
    # (used to up-weight areas around a receiver with a detection)
    cat_to_console("... Getting receiver-specific kernels (for detection)...")
    xy_ls <- lapply(1:length(xy), function(i) xy[i, ])
    # if(!is.null(cl) & !is.null(varlist)) parallel::clusterExport(cl = cl, varlist = varlist)
    detection_kernels_by_xy <-
      pbapply::pblapply(xy_ls, cl = NULL, function(xyi){
        ## Focus on area within centroid, for speed
        cat_to_console(paste0("\n... ... For receiver ", xyi$receiver_id, " ..."))
        cat_to_console("... ... ... Isolating detection centroid ...")
        xyi_centroid <- centroids[[xyi$receiver_id]]
        xyi_map <- raster::crop(map, xyi_centroid)
        ## Calc distance from receiver
        cat_to_console("... ... ... Calculating distances from the receiver ...")
        dist_around_xyi <- raster::distanceFromPoints(xyi_map, xyi)
        # raster::plot(dist_around_xyi)
        ## Calc detection probability around receiver
        cat_to_console("... ... ... Calculating detection probability ...")
        det_pr_around_xyi <- raster::calc(dist_around_xyi, fun = calc_detection_pr)
        # raster::plot(det_pr_around_xyi)
        ## Process kernels
        cat_to_console("... ... ... Processing kernel ...")
        if(!is.null(boundaries)) raster::crop(det_pr_around_xyi, boundaries)
        det_pr_around_xyi <- raster::extend(det_pr_around_xyi, raster::extent(map), value = 0)
        if(!is.null(coastline))  det_pr_around_xyi <- raster::mask(det_pr_around_xyi, coastline)
        dpr_at_xyi <- raster::extract(det_pr_around_xyi, xyi)
        if(is.na(dpr_at_xyi)) {
          warning("Detection probability is NA at receiver ", xyi$receiver_id, ".",
                  immediate. = TRUE, call. = FALSE)
          if(!is.null(coastline)) message("... Is the 'coastline' mask the right way around (masking the land and not the sea)?")
        } else if(dpr_at_xyi == 0) warning("Detection probability = NA at receiver ", xyi$receiver_id,
                                           immediate. = TRUE, call. = FALSE)
        ## Create modified kernels
        # ... (a) restricted to where detection Pr is positive (detection centroid)
        # ... (b) that show where detection Pr is positive
        det_pr_around_xyi_only <- det_pr_around_xyi
        det_pr_around_xyi_only[det_pr_around_xyi_only == 0] <- NA
        det_pr_around_xyi_positive <- det_pr_around_xyi > 0
        ## Return detection probability kernels
        out <- list(det_pr_around_xyi = det_pr_around_xyi,
                    det_pr_around_xyi_only = det_pr_around_xyi_only,
                    det_pr_around_xyi_positive = det_pr_around_xyi_positive)
        return(out)
      })
    # if(!is.null(cl)) parallel::stopCluster(cl = cl)
    names(detection_kernels_by_xy)   <- as.character(xy$receiver_id)
    receiver_specific_kernels        <- lapply(detection_kernels_by_xy, function(elm) elm$det_pr_around_xyi)
    names(receiver_specific_kernels) <- names(detection_kernels_by_xy)

    #### Calculate inverse detection Pr around each receiver
    # (used in calculations to down-weight areas, around a receiver that recorded detections,
    # ... that overlap with receivers that didn't record a detection)
    cat_to_console("... Getting receiver-specific inverse kernels...")
    receiver_specific_inv_kernels        <- lapply(receiver_specific_kernels, function(k) 1 - k)
    names(receiver_specific_inv_kernels) <- names(receiver_specific_kernels)


    ######################################
    #### Area-wide kernel (for non detection)

    #### Get dates of changes in array design
    cat_to_console("... Getting area-wide kernels (for non-detection)...")
    cat_to_console("... ... Get unique array designs...")
    # Get receiver status matrix from moorings and services (in units of days, time unit is Date)
    rs_mat <- make_matrix_receivers(moorings = moorings,
                                    services = services,
                                    delta_t = "days",
                                    as_POSIXct = NULL)
    # Get receiver status change points (dates when the array design changed)
    rs_mat_cp <- unique(rs_mat)
    # Define the time interval for each array design
    array_design_intervals <- data.frame(array_id = 1:nrow(rs_mat_cp),
                                         array_start_date = as.Date(rownames(rs_mat_cp)))
    array_design_intervals$array_end_date <-
      dplyr::lead(array_design_intervals$array_start_date) - 1
    array_design_intervals$array_end_date[nrow(array_design_intervals)] <-
      max(moorings$receiver_end_date)
    array_design_intervals$array_interval <-
      lubridate::interval(array_design_intervals$array_start_date,
                          array_design_intervals$array_end_date)

    #### For each unique array design, create the area-wide kernel surface that represents detection Pr/inverse detection Pr
    cat_to_console("... ... Get area wide kernels for each array design...")
    bkgs_by_design <-
      pbapply::pblapply(1:nrow(rs_mat_cp), function(icp){

        #### Identify active receivers on that date
        # icp <- 1
        cat_to_console(paste0("\n... ... ... For design ", icp, "/", nrow(rs_mat_cp), "..."))
        cp <- rs_mat_cp[icp, , drop = FALSE]
        rs_active <- colnames(cp)[which(cp == 1)]

        #### Pull out necessary kernels for active receivers from detection_kernels_by_xy into a list
        cat_to_console("... ... ... ... Extract detection probability kernels for active receivers...")
        detection_kernels_inv_by_rs_active <-
          lapply(rs_active, function(ra) receiver_specific_inv_kernels[[ra]])

        #### Calculate the probability of not being detected in each cell
        cat_to_console("... ... ... ... Combining detection kernels to calculate the background detection probability surfaces (this is a slow step)...")
        if(length(rs_active) == 1){
          bkg_inv <- detection_kernels_inv_by_rs_active[[1]]
        } else {
          detection_kernels_inv_for_rs_active <- raster::stack(detection_kernels_inv_by_rs_active)
          bkg_inv <- prod(detection_kernels_inv_for_rs_active)
        }
        # Get the probability of at least one detection in each grid cell:
        bkg <- 1 - bkg_inv

        #### Return surfaces
        return(list(bkg = bkg, bkg_inv = bkg_inv))
      })
    bkg_by_design     <- lapply(bkgs_by_design, function(elm) elm$bkg)
    bkg_inv_by_design <- lapply(bkgs_by_design, function(elm) elm$bkg_inv)

    #### Process outputs
    cat_to_console("... Process detection probability kernels ...")
    # For receiver-specific detection probability kernel lists,
    # add NULL elements to the list for any receivers
    # in the range 1:max(rs) that are not in rs.
    # This means we can use receiver numbers to go straight
    # ... to the correct element in the list from the integer receiver ID.
    receiver_specific_kernels <- lapply(as.integer(1:max(moorings$receiver_id)), function(i){
      if(i %in% moorings$receiver_id) return(receiver_specific_kernels[[as.character(i)]]) else return(NULL)
    })
    receiver_specific_inv_kernels <- lapply(as.integer(1:max(moorings$receiver_id)), function(i){
      if(i %in% moorings$receiver_id) return(receiver_specific_inv_kernels[[as.character(i)]]) else return(NULL)
    })

    #### Return outputs
    # Define outputs
    out <- list()
    out$receiver_specific_kernels     <- receiver_specific_kernels
    out$receiver_specific_inv_kernels <- receiver_specific_inv_kernels
    out$array_design_intervals        <- array_design_intervals
    out$bkg_surface_by_design         <- bkg_by_design
    out$bkg_inv_surface_by_design     <- bkg_inv_by_design
    # Check function duration
    t_end <- Sys.time()
    total_duration <- difftime(t_end, t_onset, units = "mins")
    cat_to_console(paste0("... flapper::acdc_setup_detection_centroids() call completed (@ ", t_end, ") after ~", round(total_duration, digits = 2), " minutes."))
    # Return outputs
    return(out)

  }


######################################
######################################
#### acdc_simplify()

#' @title Simplify the outputs of the AC* algorithms
#' @description This function simplifies the output of \code{\link[flapper]{ac}} and \code{\link[flapper]{acdc}}, by processing information from the '.acdc' elements of a \code{\link[flapper]{acdc-class}} object that hold the results of calls to the workhorse function \code{\link[flapper]{.acdc}}. This is especially useful if the AC* algorithm(s) have been applied chunk-wise, in which case the results for each acoustic chunk are returned in a list. The function aggregates information across chunks to generate a continuous time series of results and a map of where the individual could have spend more or less time over the entire time series.
#' @param acdc An \code{\link[flapper]{acdc-class}} object returned by \code{\link[flapper]{ac}} or \code{\link[flapper]{acdc}}.
#' @param mask (optional) A spatial mask (e.g., the argument passed to \code{bathy} in \code{\link[flapper]{ac}} or \code{\link[flapper]{acdc}}) to mask areas (e.g., land) from the overall map. If implemented, cells in masked areas are assigned NAs rather than a score of 0.
#' @param normalise A logical input that defines whether or not to normalise the overall map so that cell scores sum to one.
#' @param keep_chunks A logical variable that defines whether or not to retain all chunk-specific information.
#' @param ... Additional arguments (none implemented).
#' @return The function returns an object of class \code{\link[flapper]{.acdc-class}}.
#' @details If the \code{\link[flapper]{ac}} or \code{\link[flapper]{acdc}} function was implemented step-wise, this function simply extracts the necessary information and re-packages it into an \code{\link[flapper]{.acdc-class}} object. For a chunk-wise implementation, the function (a) computes the map of where the individual could have spent more or less time by aggregating the chunk-specific maps, accounting for the overlap between chunks; (b) simplifies chunk-specific records into a single contiguous time series, with re-defined time stamps from the start to the end of the time series to (c) return a \code{\link[flapper]{.acdc-class}} object.
#' @seealso The AC and ACDC algorithms are implemented by \code{\link[flapper]{ac}} and \code{\link[flapper]{acdc}}, via the \code{\link[flapper]{.acdc_pl}} and ultimately \code{\link[flapper]{.acdc}}. After simplification, \code{\link[flapper]{acdc_plot}} and \code{\link[flapper]{acdc_animate}} can be implemented to visualise time-specific results.
#' @author Edward Lavender
#' @export
#'

acdc_simplify <- function(acdc, mask = NULL, normalise = FALSE, keep_chunks = FALSE,...) {

  #### Checks
  if(!(inherits(acdc, "acdc") | !inherits(acdc, ".acdc"))){
    stop("Object of class 'acdc' expected.")
  }
  if(inherits(acdc, ".acdc")) {
    "class(acdc) == '.acdc': 'acdc' returned unchanged."
    return(acdc)
  }

  #### Set up
  out <- list(map = NULL, record = NULL, time = acdc$time, args = acdc$args, chunks = NULL, simplify = TRUE)

  #### Simplify extract outputs the algorithm has only been implemented for a single chunk
  if("map" %in% names(acdc$.acdc)) {
    out$map <- acdc$.acdc$map
    out$record <- acdc$.acdc$record

    #### Otherwise aggregate information across chunks
  } else{

    #### Compute final map from the sum of the maps for each chunk
    # Get the first map of each chunk
    maps_first <- lapply(acdc$.acdc, function(chunk) {
      map_1 <- chunk$record[[1]]$spatial[[1]]$map_timestep
      if(is.null(map_1)) {
        stop("chunk$record[[1]]$spatial[[1]]$map_timestep is NULL. In flapper::acdc(), save_spatial_record = 1L (or greater/NULL) is required to return the necessary spatial information to correct for overlapping detection time series across chunks in the summation of chunk-specific maps.")
      }
      return(map_1)
    })
    # Get the last map of each chunk
    maps_last <- lapply(acdc$.acdc, function(chunk) chunk$map)
    # Correct for the repeated influence of the first map due to the overlapping detection time series
    maps <- lapply(2:length(maps_last), function(i){
      map <- maps_last[[i]] - maps_first[[i]]
      return(map)
    })
    # Sum the adjusted maps across chunks
    out$map <- raster::brick(maps)
    out$map <- raster::calc(out$map, sum, na.rm = TRUE)

    #### Simplify records
    ## Add chunk-specific records
    out$record <- lapply(acdc$.acdc, function(chunk) chunk$record)
    ## Delete the last element of each chunk (except the last chunk) since chunks are overlapping
    out$record <- lapply(out$record, function(chunk) chunk[1:(length(chunk)-1)])
    ## Define a dataframe to adjust the time stamps recorded for each chunks
    # For chunks 2:n_chunks, we will add the time stamps reached by the previous chunk
    # ... up to the current chunk
    adjust_timestep <- lapply(out$record, function(chunk_record){
      # chunk_record <- out$record[[1]]
      dat <- chunk_record[[length(chunk_record)]]$dat
      adjustment <- dat[nrow(dat), c("timestep_cumulative", "timestep_detection")]
      return(adjustment)
    })
    adjust_timestep <- do.call(rbind, adjust_timestep)
    adjust_timestep$timestep_cumulative <- cumsum(adjust_timestep$timestep_cumulative)
    adjust_timestep$timestep_detection <- cumsum(adjust_timestep$timestep_detection)
    ## Adjust time stamps and add the chunk to the dataframe for each time stamp
    out$record <- lapply(1:length(out$record), function(i) {
      chunk_record <- out$record[[i]]
      if(i == 1) {
        adjustment <- data.frame(timestep_cumulative = 0, timestep_detection = 0)
      } else{
        adjustment <- adjust_timestep[i-1, ]
      }
      chunk_record <- lapply(chunk_record, function(t){
        t$dat$timestep_cumulative <- t$dat$timestep_cumulative + adjustment$timestep_cumulative
        t$dat$timestep_detection <- t$dat$timestep_detection + adjustment$timestep_detection
        t$dat$chunk <- i
        return(t)
      })
      return(chunk_record)
    })
    ## Flatten record list across chunks
    out$record <- purrr::flatten(out$record)

  }

  #### Normalise the final map
  if(!is.null(mask)) out$map <- raster::mask(out$map, mask)
  if(normalise) out$map <- out$map/raster::cellStats(out$map, "sum")

  #### Keep chunk-specific information, if requested
  if(keep_chunks) out$chunks <- acdc$.acdc

  #### Return outputs
  class(out) <- c(class(out), ".acdc")
  return(out)

}

######################################
######################################
#### acdc_plot()

#' @title Plot time-specific maps from the AC* algorithm(s)
#' @description This function is used to plot time-specific maps from the AC* algorithm(s). To implement the function, a \code{\link[flapper]{.acdc-class}} list from \code{\link[flapper]{ac}} or \code{\link[flapper]{acdc}} and \code{\link[flapper]{acdc_simplify}} (or \code{\link[flapper]{.acdc}}) must be supplied, from which the results can be extracted and plotted for specified time steps. For each time step, the function extracts the necessary information; sets up a blank background plot using \code{\link[raster]{plot}} and \code{\link[prettyGraphics]{pretty_axis}} and then adds requested spatial layers to this plot. Depending on user-inputs, this will usually show a cumulative map of where the individual could have spent more or less time, summed from the start of the algorithm to each time point. Coastline, receivers and acoustic centroids can be added and customised and the finalised plots can be returned or saved to file.
#' @param acdc An \code{\link[flapper]{.acdc-class}} object.
#' @param plot An integer vector that defines the time steps for which to make plots. If \code{plot = NULL}, the function will make a plot for all time steps for which the necessary information is available in \code{acdc}.
#' @param add_coastline (optional) A named list of arguments, passed to \code{\link[raster]{plot}}, to add a polygon (i.e., of the coastline), to the plot. If provided, this must contain an `x' element that contains the coastline as a spatial object (e.g., a SpatialPolygonsDataFrame: see \code{\link[flapper]{dat_coast}} for an example).
#' @param add_receivers (optional) A named list of arguments, passed to \code{\link[graphics]{points}}, to add points (i.e., receivers) to the plot. If provided, this must contain an `x' element that is a SpatialPoints object that specifies receiver locations (in the same coordinate reference system as other spatial data).
#' @param add_raster (optional) A named list of arguments, passed to \code{\link[fields]{image.plot}}, to plot the RasterLayer of possible locations that is extracted from \code{acdc}.
#' @param add_centroids (optional) A named list of arguments, passed to \code{\link[raster]{plot}}, to add the acoustic centroid to the plot.
#' @param add_additional (optional) A stand-alone function, to be executed after the background plot has been made and any specified spatial layers have been added to this, to customise the result (see Examples).
#' @param crop_spatial A logical variable that defines whether or not to crop spatial data to lie within the axis limits.
#' @param xlim,ylim,fix_zlim,pretty_axis_args Axis control arguments. \code{xlim} and \code{ylim} control the axis limits, following the rules of the \code{lim} argument in \code{\link[prettyGraphics]{pretty_axis}}. \code{fix_zlim} is a logical input that defines whether or not to fix z axis limits across all plots (to facilitate comparisons), or a vector of two numbers that define a custom range for the z axis which is fixed across all plots. \code{fix_zlim = FALSE} produces plots in which the z axis is allowed to vary flexibly between time units. Other axis options supported by \code{\link[prettyGraphics]{pretty_axis}} are implemented by passing a named list of arguments to this function via \code{pretty_axis_args}.
#' @param par_param (optional) A named list of arguments, passed to \code{\link[graphics]{par}}, to control the plotting window. This is executed before plotting is initiated and therefore affects all plots.
#' @param png_param (optional) A named list of arguments, passed to \code{\link[grDevices]{png}}, to save plots to file. If supplied, the plot for each time step is saved separately. The `filename' argument should be the directory in which plots are saved. Plots are then saved as "1.png", "2.png" and so on.
#' @param cl A cluster object created by \code{\link[parallel]{makeCluster}}. If supplied, the function loops over specified time steps in parallel to make plots. This is only implemented if plots are saved to file (i.e., \code{png_param} is supplied). If supplied, the connection to the cluster is closed within the function.
#' @param verbose A logical variable that defines whether or not relay messages to the console to monitor function progress.
#' @param check A logical variable that defines whether or not to check user inputs to the function before its initiation.
#' @param ... Additional arguments, passed to \code{\link[raster]{plot}}, to customise the blank background plot onto which spatial layers are added, such as \code{xlab}, \code{ylab} and \code{main}.
#'
#' @return The function plots the results of the AC* algorithm(s) at specified time steps, with one plot per time step. Plots are saved to file if \code{png_param} is supplied.
#' @examples
#' #### Step (1): Implement AC* algorithm(s)
#' # ... see examples via ac() and acdc()
#'
#' #### Step (2): Simplify outputs of the AC* algorithm(s)
#' dat_acdc <- acdc_simplify(dat_acdc)
#'
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
#'
#' #### Example (7) To plot the overall map, you can also just use a
#' # ... a raster plotting function like prettyGraphics::pretty_map()
#' ext <- update_extent(raster::extent(dat_coast), -1000)
#' prettyGraphics::pretty_map(x = ext,
#'                            add_rasters = list(x = dat_acdc$map),
#'                            add_points = list(x = rsp, pch = "*", col = "red"),
#'                            add_polys = list(x = dat_coast, col = "lightgreen"),
#'                            crop_spatial = TRUE,
#'                            xlab = "Easting", ylab = "Northing"
#'                            )
#'
#' @seealso This function is typically used following calls to \code{\link[flapper]{acdc}} and \code{\link[flapper]{acdc_simplify}}.
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
    ## Check object class
    if(!inherits(acdc, ".acdc")){
      if(inherits(acdc, "acdc")) {
        stop("'acdc' must be converted from 'acdc' to '.acdc. via flapper::acdc_simplify() before implementing this function.")
      } else{
        stop("'acdc' must be of class '.acdc'.")
      }
    }
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
  cat_to_console("... Making plots for each time step ...")
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

#' @title Create a html animation of the AC* algorithm(s)
#' @description This function is a simple wrapper for \code{\link[flapper]{acdc_plot}} and \code{\link[animation]{saveHTML}} which creates an animation of the AC* algorithm(s) over time. To implement this function, a named list of arguments for \code{\link[flapper]{acdc_plot}}, which creates the plots, must be supplied. This is embedded within \code{\link[animation]{saveHTML}}, which creates a folder in the working directory named `images' that contains a .png file for each time step and an animation as a .html file.
#' @param expr_param A named list of arguments, passed to \code{\link[flapper]{acdc_plot}}, to create plots.
#' @param html_name A string that defines the name of the html file (see `htmlfile' argument in \code{\link[animation]{saveHTML}}).
#' @param image_name A string that defines the names of the individual .png files creates (see `img.name' argument in \code{\link[animation]{saveHTML}}).
#' @param html_title,html_description Character strings that provide a title and a description that are displayed within the html animation (see `title' and `description' arguments in \code{\link[animation]{saveHTML}}).
#' @param navigator A logical variable that defines whether or not to add a navigator panel to the animation (see `navigator' argument in \code{\link[animation]{saveHTML}}).
#' @param ani_height,ani_width,ani_res Numbers that define the size and the resolution of the animation (see `ani.height' `ani.width' and `ani.res' arguments in \code{\link[animation]{ani.options}}).
#' @param interval A number that defines the time interval between sequential frames (see `interval' argument in \code{\link[animation]{ani.options}}).
#' @param verbose A logical or character variable that defines whether or not, or what, to write as a footer to the html animation (see `verbose' argument in \code{\link[animation]{ani.options}}).
#' @param ... Additional arguments passed to \code{\link[animation]{ani.options}}.
#'
#' @return The function produces an animation in .html format in the working directory (or a sub-directory of this). A folder named `images' is also produced which contains the images for each time step. The `css' and `js' folders are also produced by \code{\link[animation]{saveHTML}} which creates the animation.
#'
#' @examples
#' dir_current <- getwd()
#' setwd(tempdir())
#' dat_acdc <- acdc_simplify(dat_acdc)
#' acdc_animate(expr_param = list(acdc = dat_acdc,
#'                                add_coastline = list(x = dat_coast, col = "darkgreen"),
#'                                plot = 1:5,
#'                                fix_zlim = FALSE)
#'                                )
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
    #### Checks
    ## animation package
    if (!requireNamespace("animation", quietly = TRUE)) {
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
