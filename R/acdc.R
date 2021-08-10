######################################
######################################
#### acdc_setup_mobility()

#' @title Examine the constant `mobility' assumption of the ACDC algorithm
#' @description In the simplest and fastest (and only) version of the ACDC algorithm currently implemented by \code{\link[flapper]{flapper}}, the rate at which acoustic centroids expand and contract depends on a single `mobility' parameter that describes how the uncertainty in an individual's location changes through time according to the passive acoustic telemetry data. These changes essentially reflect the maximum horizontal distance that an individual can move in the time period between archival observations. However, in some situations, a fixed parameter for horizontal movement may not be well supported by archival data, particularly for animals that change dramatically in depth and for which changes in depth are likely to curtail the extent of horizontal movement*. Thus, this function investigates the extent to which the horizontal distance an animal could travel changes through time if the mobility parameter is instead conceptualised as the diagonal distance between sequential depth observations.
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
#' *If the horizontal distances that an individual could travel are not very variable, then the benefits of a single mobility parameter are likely to outweigh the costs. On the other hand, substantial variation in the horizontal distances that an individual could travel may suggest that pre-processed acoustic centroids are inappropriate; in this situation, the correct centroids could be computed on-the-fly within the ACDC algorithm, but this is not currently implemented. However, the ACDCPF/ACDCMP algorithms can account for this by the incorporation of movement models within/between acoustic centroids. (Indeed, even if the horizontal mobility remains approximately constant, the integration of a movement model via ACDCPF/MP can still be beneficial.)
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
#' # ... between acoustic centroids for some applications via the ACDCPF/MP
#' # ... algorithms.
#'
#' @seealso \code{\link[flapper]{acdc_setup_mobility}}, \code{\link[flapper]{acdc_setup_n_centroids}} and \code{\link[flapper]{acdc_setup_centroids}} are used to set up the ACDC algorithm as implemented by \code{\link[flapper]{acdc}}.
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

#' @title Suggest the number of centroids for the ACDC algorithm
#' @description The ACDC algorithm requires a list of acoustic centroids from \code{\link[flapper]{acdc_setup_centroids}}. The number of centroids that is required depends on (a) the duration between detections; (b) the distances among receivers; (c) the area of interest; and (d) other considerations such as the spatial extent and resolution of the bathymetry data that is required by the ACDC algorithm and presence or absence of a long tail of archival observations after the final acoustic detection. This function implements two approaches to facilitate a sensible choice of the number of centroids for the ACDC algorithm (see Details). This ensures that the number of centroids is sufficient to cover requirements but not so large that the centroids become slow to compute and unwieldy. The function requires a time series of detections, the deployment details of passive acoustic telemetry receivers and a mobility parameter. Given these inputs, the function returns a suggested upper and lower bound for the number of centroids.
#'
#' @param detections A POSIXct vector of time stamps when detections were made (for a particular individual).
#' @param moorings A dataframe that defines passive acoustic telemetry receiver locations and deployment periods. This must contain the following columns: `receiver_id', a unique identifier of each receiver; `receiver_lat' and `receiver_long', the latitude and longitude of each receiver in decimal degrees; and `receiver_start_date' and `receiver_end_date', the start and end time of each receiver's deployment period (see \code{\link[flapper]{dat_moorings}} for an example).
#' @param mobility A number that defines the distance (m) that an individual could move in the time period between archival observations (see also \code{\link[flapper]{acdc_setup_centroids}}).
#' @param double A logical variable that defines whether or not to double the minimum number of centroids given by approach two (see Details).
#' @param hist A logical value that defines whether or not to plot the distribution of gaps across the detection time series as a histogram.
#' @param ... Additional arguments passed to \code{\link[prettyGraphics]{pretty_hist}}.
#'
#' @details
#' \subsection{Method (1)}{The first approach provides a reasonable upper bound for the number of centroids. This approach is based on the gaps between detections. During this time, acoustic centroids increase in size, reflecting the increasing uncertainty in the location of an individual, until the half way point between detections. At this point, uncertainty in the individual's location is maximised and the acoustic centroids reach their maximum size. Thereafter, uncertainty in the individual's location and the size of the acoustic centroids decrease towards the receiver at which the individual was next detected. Under this perspective, the minimum number of centroids is the number required such that the centroids continue to increase in size until the halfway point between the two most temporally distance detections. More specifically, assuming that archival time steps are two-minutes in duration, the minimum number of centroids is one quarter of the duration (minutes) of the longest gap between acoustic detections. To implement this approach, a time series of \code{detections} (for the individual and time period for which the ACDC algorithm will be implemented) needs to be provided. This approach provides a sensible upper bound for the number of centroids, since it allows centroids to continue to expand to their maximum possible size, given the data and the assumptions of the algorithm. However, even with modest gaps between detections and a modest movement capacity, this may suggest a very large area, and it may take a long time to compute the requisite number of centroids using \code{\link[flapper]{acdc_setup_centroids}}. Moreover, in practice, the area of interest may be smaller.}
#'
#' \subsection{Method (2)}{The second approach provides a reasonable lower bound for the number of centroids. This approach is based on the locations of receivers and the assumption that, when an individual is detected by two different receivers, at the halfway point between detections, its potential location is described by the intersection of the acoustic centroid around the receiver at which it was previously detected and the centroid around the receiver at which it is next detected (evaluated at the halfway point between detections). Under this perspective, the minimum centroid size is half of the distance between the furthest two receivers that were operational at the same time. (This could be restricted to the subset of receivers at which the individual was detected sequentially but this is not implemented.) In practice, it is advisable to double this minimum number to ensure a reasonable degree of overlap between the centroids of the two receivers. To implement this approach, a dataframe that contains receiver locations and deployment periods (\code{moorings}) must be supplied, along with the \code{mobility} parameter that describes how far an individual could move in any archival time step. This approach provides a sensible lower bound for the number of centroids that is array-specific. However, it does not account for the full span of movements that are possible under the ACDC algorithm.}
#'
#' \subsection{Other considerations}{In practice, the most appropriate number of centroids is likely to be a compromise between these minimum and maximum values that depends on other considerations, particularly the spatial scale of the study, the effect of boundaries on the final centroid size and, perhaps in some cases, the geographical range of the species. For acoustic time series that are followed by a long tail of archival observations after final acoustic detection, both approaches may underestimate the number of centroids required to capture the increasing uncertainty in the individual's location over this time. However, in this case, it is advisable to use the depth-contour (DC) algorithm (see \code{\link[flapper]{dc}}) at some point after the last acoustic detection since the influence of that observation on putative patterns of space use will decay through time and the DC algorithm is more computationally efficient.}
#'
#' @return The function returns an integer vector with the upper and lower suggested value for the number of centroids from methods (1) and (2). The parameters used to generate these suggestions (i.e., \code{detections}, \code{moorings}, \code{mobility} and \code{double}) are also included in a `param' attribute.
#'
#' @seealso This function is designed to facilitate an informed choice for the `n_timesteps' argument in \code{\link[flapper]{acdc_setup_centroids}}. \code{\link[flapper]{acdc_setup_mobility}} can also guide the implementation of this function. This underpins the ACDC algorithm, which is implemented by \code{\link[flapper]{acdc}}.
#'
#' @examples
#' n_timesteps <-
#'   acdc_setup_n_centroids(dat_acoustics$timestamp[dat_acoustics$individual_id == 25],
#'                          dat_moorings,
#'                          mobility = 200,
#'                          double = TRUE)
#' utils::str(n_timesteps)
#'
#' @author Edward Lavender
#' @export

acdc_setup_n_centroids <- function(detections, moorings, mobility, double = TRUE, hist = TRUE,...){

  #### Checks
  # Check moorings contains required information
  check_names(input = moorings,
              req = c("receiver_id",
                      "receiver_start_date", "receiver_end_date",
                      "receiver_lat", "receiver_long"),
              extract_names = colnames,
              type = all)

  #### Approach (1): Determine the maximum gap between detections
  # Calculate the duration of gaps
  gaps <- Tools4ETS::serial_difference(detections, units = "mins")
  gaps <- as.numeric(gaps)
  max_gap <- max(gaps, na.rm = TRUE)
  # Visualise the distribution of gaps
  if(hist) {
    prettyGraphics::pretty_hist(gaps,...)
    graphics::abline(v = max_gap, col = "red", lty = 3)
  }
  # The minimum number of time steps under this approach is max_gap/2/2
  # ... assuming two minute time stamps
  # ... and since after half way the centroid will start to shrink
  minimum_n_timesteps_1 <- max_gap/2/2

  #### Approach (2): Determine the number of timesteps to overlap centroids between receivers

  ## Identify time intervals over which receivers were deployed
  moorings$interval <- lubridate::interval(moorings$receiver_start_date, moorings$receiver_end_date)

  ## Define the minimum number of timesteps for each combination of receivers
  # ... that was deployed for an overlapping interval
  minimum_n_timesteps_2 <-
    pbapply::pbsapply(unique(moorings$interval), function(interval){
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

  #### Return suggestions
  minimum_n_timesteps <- c(minimum_n_timesteps_1, minimum_n_timesteps_2)
  minimum_n_timesteps <- as.integer(ceiling(minimum_n_timesteps))
  attributes(minimum_n_timesteps)$method <- 1:2
  attributes(minimum_n_timesteps)$param <- list(detections = detections,
                                                moorings = moorings,
                                                mobility = mobility,
                                                double = double)
  return(minimum_n_timesteps)
}


######################################
######################################
#### acdc_setup_centroids()

#' @title Setup the acoustic centroids required for the ACDC algorithm
#' @description This function produces the acoustic centroids required by the acoustic-centroid depth-contour (ACDC) algorithm.
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
#' @details Given only a detection at a particular receiver at a particular time, the detection range of the receiver and the movement speed of the animal, an acoustic centroid defines the possible locations of a detected individual at that time or a subsequent time. More specifically, when an individual is located at a receiver, its location must be within some radius of that receiver defined by the maximum detection distance. This radius expands with the duration since detection. This function defines, for each receiver, a list of acoustic centroids that reflect the possible locations of an individual were it to have been detected from 0 to \code{n_timesteps} ago at that receiver, accounting for the coastline and within a defined area if necessary. Using the observed detection data, the ACDC (\code{\link[flapper]{acdc}}) and ACDCPF/MP algorithms pull the relevant centroids out of this list, which substantially saves computation time because acoustic centroids are not computed on-the-fly. These centroids are processed within these algorithms (e.g., if an individual is detected at two different receivers, then at the halfway point between detections its location must be within the intersection of the relevant centroids for these two receivers) and combine them with other information to reconstruct where an individual could have been over time.
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
#' # We'll get their overlaps via flapper::get_detection_centroid_overlaps
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
#### acdc()

#' @title The acoustic-centroid depth-contour (ACDC) algorithm
#' @description This function implements the acoustic-centroid depth-contour (ACDC) algorithm. This is an approach that integrates acoustic detections and depth observations to infer where benthic/demersal animals could have spent more or less time over the period of observations.
#'
#' To implement the function, a dataframe (or list) of passive acoustic telemetry detections is required (\code{acoustics}), alongside a dataframe of depth observations (\code{archival}). At each time step, the algorithm integrates information from past and future acoustic detections in the form of acoustic centroids and information from depth observations in the form of depth contours to determine the possible locations of an individual in an area (see Details).
#'
#' Under the default options, the approach is implemented step-wise (i.e., step-by-step across the whole time series). The result is a named list of outputs, including a record of the results for each time step, as well as a cumulative map of where the individual could have spent more or less time summed across the whole time series. Alternatively, the approach can be implemented chunk-wise, in which case the acoustic time series is split into chunks (e.g., hourly, daily, monthly segments) and the algorithm is implemented within each chunk step-by-step. The main benefits of this approach are that it can be used to reconstruct putative patterns in space use over biologically meaningful periods separately and/or the chunk-wise implementation can be parallelised, improving computation time. (Chunk-wise results results are easily combined across the duration of the original time series without the loss of information via \code{\link[flapper]{acdc_simplify}}.). This option is implemented if (a) a list, rather than a dataframe, of acoustic detections is provided (via \code{acoustics}); (b) the user specifies that the time series should be split into chunks of a particular duration before the algorithm is initiated (via the \code{split} argument); and/or (c) the algorithm is implemented on a cluster via \code{cl}, in which case the acoustic time series is split (if necessary) into user-defined or automatically defined chunks prior to computation. In this case, the result is a named list of outputs, as described above, but in which the results for each chunk are returned separately. If the chunks have been implemented simply to improve computation time via parallelisation, then the maps of space use for each chunk can be combined easily to generate a single, overall map of space use via \code{\link[flapper]{acdc_simplify}}.
#'
#' @param acoustics A dataframe, or a list of dataframes, that contains passive acoustic telemetry detection time series (see \code{\link[flapper]{dat_acoustics}} for an example) for a single individual. Each dataframe should contain the following columns: an integer vector of receiver IDs, named `receiver_id'; a POSIXct vector of time stamps when detections were made, named `timestamp'; and a numeric vector of those time stamps, named `timestamp_num'. If a list of dataframes is supplied, dataframes must be refer to the detections of a single individual and be ordered by time (e.g., in hourly chunks). The algorithm will be implemented on each dataframe, termed `chunk', either in sequence or parallel. Any empty or NULL elements will be removed automatically.
#' @param archival A dataframe that contains depth time series (see \code{\link[flapper]{dat_archival}} for an example). This should contain the following columns: a numeric vector of observed depths, named `depth'; a POSIXct vector of time stamps when observations were made, named `timestamp'; and a numeric vector of those time stamps, named `timestamp_num'. Depths should be recorded in the same units and with the same sign as the bathymetry data (see \code{bathy}). Absolute depths (m) are suggested. Unlike the detection time series, archival time stamps are assumed to have occurred at regular intervals. Two-minute intervals are currently assumed.
#' @param bathy A \code{\link[raster]{raster}} that defines the bathymetry across the area within which the individual could have moved. This must be recorded in the same units and with the same sign as the depth observations (see \code{archival}). The coordinate reference system should be the Universal Transverse Mercator system, with distances in metres (see also \code{\link[flapper]{acdc_setup_centroids}}).
#' @param map (optional) A blank \code{\link[raster]{raster}}, with the same properties (i.e., dimensions, resolution, extent and coordinate reference system) as the bathymetry raster (see \code{bathy}), but in which all values are 0. If \code{NULL}, this is computed internally, but supplying a pre-defined raster can be more computationally efficient if the function is applied iteratively (e.g., over multiple individuals).
#' @param detection_range A number that defines the maximum detection range (m) at which an individual could be detected from a receiver (see also \code{\link[flapper]{acdc_setup_centroids}}).
#' @param mobility A number that defines the distance (m) that an individual could move in the time period between archival observations (see also \code{\link[flapper]{acdc_setup_centroids}}).
#' @param calc_depth_error A function that returns the depth error around a given depth. This should accept a single depth value (from \code{archival$depth}) and return two numbers that, when added to that depth, define the range of depths on the bathymetry raster (\code{bathy}) that the individual could plausibly have occupied at any time, given its depth. Since the depth errors are added to the individual's depth, the first number should be negative (i.e., the individual could have been slightly shallower that observed) and the second positive (i.e., the individual could have been slightly deeper than observed). For example, the constant function \code{calc_depth_error = function(...) c(-2.5, 2.5)} implies that the individual could have occupied bathymetric cells whose depth lies within the interval defined by the observed depth + (-2.5) and + (+2.5) m. The appropriate form for \code{calc_depth_error} depends on measurement error for the depth observations in \code{archival} and bathymetry (\code{bathy}) data, as well as the tidal range (m) across the area (over the duration of observations), but this implementation allows the depth error to depend on depth and for the lower and upper error around an observation to differ.
#' @param acc_centroids A list of acoustic centroids, with one element for each number from \code{1:max(acoustics$receiver_id)}, from \code{\link[flapper]{acdc_setup_centroids}}.
#' @param plot An integer vector that defines the time steps for which to return the necessary spatial information required to plot the plausible locations of the individual, given detection and depth time series. \code{plot = 0L} suppresses the return of this information and \code{plot = NULL} returns this information for all time steps. If the algorithm is applied chunk-wise, this spatial information must be returned for at least the first time step (the default) to aggregate maps across chunks (see \code{\link[flapper]{acdc_simplify}}). This information can also be used to plot time-specific results of the algorithm using \code{\link[flapper]{acdc_plot}} and \code{\link[flapper]{acdc_animate}}.
#' @param plot_ts A logical input that defines whether or not to the plot detection and depth time series before the algorithm is initiated. This provides a useful visualisation of the extent to which they overlap.
#' @param verbose A logical variable that defines whether or not to print messages to the console or to file to relay function progress. If \code{con = ""}, messages are printed to the console (which is only supported if the algorithm is not implemented in parallel: see below); otherwise, they are written to file (see below).
#' @param con If \code{verbose = TRUE}, \code{con} is character string defines how messages relaying function progress are returned. If \code{con = ""}, messages are printed to the console (unless redirected by \code{\link[base]{sink}}), an approach that is only implemented if the function is not implemented in parallel. Otherwise, \code{con} defines the directory into which to write .txt files, into which messages are written to relay function progress. This approach, rather than printing to the console, is recommended for clarity, speed and debugging. If the algorithm is implemented step-wise, then a single file is written to the specified directory named acdc_log.txt. If the algorithm is implemented chunk-wise, then an additional file is written for each chunk (named dot_acdc_log_1.txt, dot_acdc_log_2.txt and so on), with the details for each chunk.
#' @param progress (optional) If the algorithm is implemented step-wise, \code{progress} is an integer (\code{1}, \code{2} or \code{3}) that defines whether or not to display a progress bar in the console as the algorithm moves over acoustic time steps (\code{1}), the archival time steps between each pair of acoustic detections (\code{2}) or both acoustic and archival time steps (\code{3}), in which case the overall acoustic progress bar is punctuated by an archival progress bar for each pair of acoustic detections. This option is useful if there is a large number of archival observations between acoustic detections. Any other input will suppress the progress bar. If the algorithm is implemented for chunks, inputs to \code{progress} are ignored and a single progress bar is shown of the progress across acoustic chunks.
#' @param split A character string that defines the time unit used to split acoustic time series into chunks (e.g., \code{"12 hours"}). If provided, this must be supported by \code{\link[lubridate]{floor_date}} (otherwise, a pre-defined list of acoustic time series can be passed to \code{acoustics}, e.g., specifying seasonal chunks). If \code{split = NULL} and a cluster has been specified (see \code{cl}) (and \code{acoustics} is a dataframe), then the acoustic time series is automatically split into chunks and the algorithm implemented for each chunk in parallel.
#' @param cl A cluster object created by \code{\link[parallel]{makeCluster}} to implement the algorithm in parallel. If supplied, the algorithm is implemented for each chunk in a list of acoustic time series as supplied by the user (if \code{acoustics} is a list) or of the time units specified via \code{split} by the user or defined automatically based on the number of nodes in the cluster if \code{split = NULL}.
#' @param ... Additional arguments (none implemented).
#'
#' @details The acoustic-centroid depth-contour (ACDC) algorithm is an approach which integrates acoustic detections and depth observations to infer the possible locations of benthic or demersal animals within an area over some time interval. The locational information provided by acoustic detections is represented by acoustic centroids, which are areas around receivers that define where an individual could have been at each timepoint given the spatiotemporal pattern of detections at receivers, a model of detection probability and a movement parameter. The locational information provided by depth observations is represented by depth contours, which are areas that define where an individual could have been at each time point given its depth and the local bathymetry.
#'
#' In outline, the crux of the approach is the recognition that acoustic detections typically occur irregularly, while archival observations occur at regular intervals. Each detection anchors our knowledge of the location of an individual around a particular receiver (assuming that all detections are true detections). As time passes between acoustic detections, our uncertainty about the geographical location of an individual expands around the receiver at which it was detected before shrinking towards the receiver at which it was next detected. During this time, regular depth observations restrict the number of possible locations in which the individual could have been located at each time step.
#'
#' More specifically, when an individual is detected, it must be within some radius---say 800 m---of that receiver. This is the starting acoustic centroid. With a more-refined model of detection probability, it may be possible to predict more precisely where the individual is likely to have been within this centroid (but this approach is not yet implemented). The observed depth at this time further restricts the positions in which the individual could have been, assuming a benthic/demersal lifestyle and a non-homogeneous bathymetric landscape. Moving forward in time, a number of depth records may be made before another acoustic detection. During this time, our uncertainty about where the individual could have been gets larger, because it could have moved further away from the receiver, so the acoustic centroids that define this uncertainty expand to a maximum size at the halfway point between acoustic detections. After that, the individual must have started to move towards the receiver at which it was next detected, so these acoustic centroids start to shrink towards that receiver. If the individual was detected by different receivers, the overlap between the centroids of these two receivers at the halfway point defines the set of positions in which the individual could have been at this time. Thereafter, our uncertainty in the individual's location is given by the overlap between the expansion of this centroid region and the contraction of the centroid around the receiver at which it was next detected. Thus, when the individual is detected again, our uncertainty about where it could have been collapses to the detection radius around the next receiver, possibly weighted by a model of detection probability around this receiver (although that is not yet implemented). The rate of change in centroid size depends a movement parameter that describes an average swimming speed, which will depend on some underlying estimated behavioural state (although that is not yet implemented).
#'
#' The result is a map that shows where the individual could have spent more or less (or no) time over the time interval under construction. The main limitation of this approach is that reconstructs where the individual could have been, but not where it was. In reality, the individual's current position constrains where it can go next. The ACDCPF/MP algorithms are extensions of this approach that incorporate a movement model for this reason.
#'
#' @return The function returns a \code{\link[flapper]{acdc-class}} object. If a connection to write files has also been specified, an overall log (acdc_log.txt) as well as chunk-specific logs from calls to \code{\link[flapper]{.acdc}}, if applicable, are written to file.
#'
#' @seealso The `depth-contour' component of the algorithm can be implemented via \code{\link[flapper]{dc}}. For more information on the ACDC algorithm, see \code{\link[flapper]{acdc_setup_n_centroids}} and \code{\link[flapper]{acdc_setup_centroids}} which defined the acoustic centroids required by this function and \code{\link[flapper]{acdc_setup_mobility}} which examines the assumption of a constant `mobility' parameter; \code{\link[flapper]{.acdc}}, the back-end workhorse of this function; \code{\link[flapper]{acdc_simplify}} which simplifies the outputs of the \code{\link[flapper]{acdc}} function; and \code{\link[flapper]{acdc_plot}} and \code{\link[flapper]{acdc_animate}}, which visualise the results.
#'
#' @examples
#' #### Step (1) Implement setup_acdc_*() steps
#' # ... Define acoustic centroids required for ACDC algorithm (see setup_acdc_centroids())
#'
#' #### Step (2) Prepare movement time series for algorithm
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
#' #### Example (1) Implement ACDC algorithm with default arguments
#' # This implements the algorithm on a single core, printing messages
#' # ... to the console to monitor function progress.
#' out_acdc <- acdc(acoustics = acc,
#'                  archival = arc,
#'                  bathy = dat_gebco,
#'                  detection_range = 425,
#'                  mobility = 200,
#'                  calc_depth_error = function(...) c(-2.5, 2.5),
#'                  acc_centroids = dat_centroids
#'                  )
#' # The function returns a list with four elements
#' # ... .acdc contains the results of the algorithm, implemented by the back-end
#' # ... function .acdc(). The other elements provide the time series
#' # ... for each chunk, the time of the algorithm and a list of user inputs
#' summary(out_acdc)
#'
#' #### Example (2): Write messages to file to monitor function progress via 'con'
#' out_acdc <- acdc(acoustics = acc,
#'                  archival = arc,
#'                  bathy = dat_gebco,
#'                  detection_range = 425,
#'                  mobility = 200,
#'                  calc_depth_error = function(...) c(-2.5, 2.5),
#'                  acc_centroids = dat_centroids,
#'                  con = tempdir()
#'                  )
#' acdc_log <- readLines(paste0(tempdir(), "/acdc_log.txt"))
#' utils::head(acdc_log, 10)
#' file.remove(paste0(tempdir(), "/acdc_log.txt"))
#'
#' #### Example (3): Implement the algorithm and return plotting information
#' # Specify plot = NULL to include plotting information for all time steps
#' # ... or a vector to include this information for specific time steps
#' out_acdc <- acdc(acoustics = acc,
#'                  archival = arc,
#'                  bathy = dat_gebco,
#'                  detection_range = 425,
#'                  mobility = 200,
#'                  calc_depth_error = function(...) c(-2.5, 2.5),
#'                  acc_centroids = dat_centroids,
#'                  plot = NULL
#'                  )
#'
#' #### Example (4): Implement the algorithm in parallel by supplying a cluster
#' # If verbose = TRUE (the default), it is necessary to specify a directory
#' # ... into which dot_acdc_log_*.txt files are saved (i.e., messages
#' # ... cannot be written to the console in parallel)
#' out_acdc <- acdc(acoustics = acc,
#'                  archival = arc,
#'                  bathy = dat_gebco,
#'                  detection_range = 425,
#'                  mobility = 200,
#'                  calc_depth_error = function(...) c(-2.5, 2.5),
#'                  acc_centroids = dat_centroids,
#'                  con = tempdir(),
#'                  cl = parallel::makeCluster(2L)
#'                  )
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
#' # Note that there are now four elements in .acdc, one for each chunk
#' # Likewise, there are four elements in ts_by_chunk,
#' # ... containing the movement time series for each chunk.
#' summary(out_acdc)
#' # Note that the last observation of each time series overlaps with the
#' # ... first observation for the next chunk, to prevent loss of information
#' lapply(out_acdc$ts_by_chunk,
#'   function(chunk) chunk$acoustics[c(1, nrow(chunk$acoustics)), ])
#'
#' #### Example (5) Biologically meaningful chunks can be specified via
#' # .. the 'split' argument or by passing a list of acoustic time series
#' # .. already split by list of dataframes to 'acoustics'
#' ## Using the split argument:
#' out_acdc <- acdc(acoustics = acc,
#'                  archival = arc,
#'                  bathy = dat_gebco,
#'                  detection_range = 425,
#'                  mobility = 200,
#'                  calc_depth_error = function(...) c(-2.5, 2.5),
#'                  acc_centroids = dat_centroids,
#'                  con = tempdir(),
#'                  cl = parallel::makeCluster(2L),
#'                  split = "2 hours"
#'                  )
#' ## Passing a list of dataframes
#' # This option can provide more flexibility than split, which only
#' # ... understands time categories supported by lubridate::floor_date()
#' # ... This example could also be used using split, as described above,
#' # ... but this is not the case for all time categories (e.g., seasons).
#' acc$chunk <- cut(acc$timestamp, "2 hours")
#' acc_ls <- split(acc, acc$chunk)
#' out_acdc <- acdc(acoustics = acc_ls,
#'                  archival = arc,
#'                  bathy = dat_gebco,
#'                  detection_range = 425,
#'                  mobility = 200,
#'                  calc_depth_error = function(...) c(-2.5, 2.5),
#'                  acc_centroids = dat_centroids,
#'                  con = tempdir(),
#'                  cl = parallel::makeCluster(2L)
#'                  )
#'
#' #### Example (5) Implement the algorithm for multiple individuals
#' # ... To do this, it is necessary to apply the function iteratively
#' # ... to each individual.
#' # Pre-processing to define computation time
#' # ... E.g., careful definition of time series
#' # ... E.g., define 'map' argument
#' map <- dat_gebco
#' map <- raster::setValues(map, 0)
#' # Define cluster
#' cluster <- FALSE
#' if(cluster){
#'   cl <- parallel::makeCluster(2L)
#'   parallel::clusterExport(cl = cl, varlist = c("acdc",
#'                                                "dat_archival",
#'                                                "dat_gebco",
#'                                                "map",
#'                                                "dat_centroids"
#'                                                 ))
#' } else cl<- NULL
#' # Implement algorithm for each individual
#' acdc_out_by_id <-
#'   pbapply::pblapply(split(dat_acoustics, dat_acoustics$individual_id), cl = cl, function(acc){
#'     # Define individual-specific folder in which to save function logs
#'     dir_global <- paste0(tempdir(), "/")
#'     dir_id     <- paste0(dir_global, acc$individual_id[1], "/")
#'     if(!dir.exists(dir_id)) dir.create(dir_id)
#'     # Focus on a small sample of time series for speed
#'     acc <- acc[1:3, ]
#'     # Isolate archival data for individual
#'     arc <- dat_archival[dat_archival$individual_id == acc$id[1], ]
#'     # Implement algorithm
#'     acdc_out <- acdc(acoustics = acc,
#'                      archival = dat_archival,
#'                      bathy = dat_gebco,
#'                      map = map,
#'                      detection_range = 425,
#'                      mobility = 200,
#'                      calc_depth_error = function(...) c(-2.5, 2.5),
#'                      acc_centroids = dat_centroids,
#'                      plot = 1:10L,
#'                      con = dir_id
#'                      )
#'     # Include logs in output
#'     acdc_log <- lapply(list.files(dir_id, full.names = TRUE), readLines)
#'     # Simplify the results at this stage or outside of this loop
#'     # ... using acdc_simplify()
#'     # Return results for specified individual
#'     out <- list(acdc_id = acc$individual_id[1], acdc_out = acdc_out, acdc_log = acdc_log)
#'     return(out)
#'   })
#' if(!is.null(cl)) parallel::stopCluster(cl)
#' summary(acdc_out_by_id)
#'
#' #### Step (3) Simplify the function outputs
#' # This step aggregates information across chunks, which is necessary to
#' # ... plot information aggregated across all chunks (see below).
#'
#' #### Step (4) Examine function outputs, e.g., via plotting
#' # See acdc_plot() and acdc_animate() to visualise the results
#' # ... (either for a specific chunk or aggregated across all chunks
#' # ... using acdc_simplify() as described above).
#'
#' @author Edward Lavender
#' @export
#'

acdc <- function(
  acoustics,
  archival,
  bathy,
  map = NULL,
  detection_range,
  mobility,
  calc_depth_error = function(...) c(-2.5, 2.5),
  acc_centroids,
  plot = 1L,
  plot_ts = TRUE,
  verbose = TRUE,
  con = "",
  progress = 1L,
  split = NULL,
  cl = NULL,...
  ){


  ######################################
  #### Set up

  #### Initiate function
  t_onset <- Sys.time()
  message(paste0("flapper::acdc() called (@ ", t_onset, ")..."))

  #### A list to store overall outputs
  out <- list(.acdc = NULL, ts_by_chunk = NULL, time = NULL, args = NULL)
  out$args <- list(acoustics = acoustics,
                   archival = archival,
                   bathy = bathy,
                   map = map,
                   detection_range = detection_range,
                   mobility = mobility,
                   calc_depth_error = calc_depth_error,
                   acc_centroids = acc_centroids,
                   plot = plot,
                   verbose = verbose,
                   con = con,
                   progress = progress,
                   split = split,
                   cl = cl,
                   dots = list(...))

  #### Check parallelisation options
  if(is.null(cl)) n_cores <- 1 else n_cores <- length(cl)
  # if(n_cores == 1 & !is.null(split)) {
  #  message("Input to 'split' is ignored since cl = NULL.")
  #  split <- NULL
  # }
  if(inherits(acoustics, "list") & !is.null(split)) message("Input to 'split' ignored since inherits(acoustics, 'list') == TRUE.")

  #### Define function for printing messages to file or console
  ## Check the connection for writing files, if applicable

  if(con != ""){
    if(!verbose) {
      message("Input to 'con' ignored since verbose = FALSE.")
    } else {
      # Check directory
      con <- check_dir(input = con, check_slash = TRUE)
      con_dir <- con
      # Define file
      con <- paste0(con_dir, "acdc_log.txt")
      if(!file.exists(con)){
        message(paste0(con, " does not exist: attempting to write file in specified directory..."))
        file.create(file1 = con)
        message("... Blank file successfully written to file.")
      }
    }
  } else{
    if(n_cores > 1) stop("con = '' is not implemented in parallel (!is.null(cl)). Please supply a directory.")
  }
  ## Define function
  append_messages <- ifelse(con == "", FALSE, TRUE)
  cat_to_cf <- function(..., message = verbose, file = con, append = append_messages){
    if(message) cat(paste(..., "\n"), file = con, append = append)
  }

  #### Checks
  ## Formally initiate function and implement remaining checks
  cat_to_cf(paste0("flapper::acdc() called (@ ", t_onset, ")..."))
  out$time <- data.frame(event = "onset", time = t_onset)
  cat_to_cf("... Checking user inputs...")
  # Check acoustics contains required column names and correct variable types
  if(!inherits(acoustics, "list")) acoustics_tmp <- list(acoustics) else acoustics_tmp <- acoustics
  length_acoustics_tmp <- length(acoustics_tmp)
  lapply(acoustics_tmp, function(acc) {
    check_names(arg = "acoustics",
                input = acc,
                req = c("timestamp", "timestamp_num", "receiver_id"),
                extract_names = colnames,
                type = all)
    check_class(input = acc$timestamp, to_class = "POSIXct", type = "stop")
    check_class(input = acc$receiver_id, to_class = "integer", type = "stop")
  })
  # Check archival contains required column names and correct variable types
  check_names(input = archival,
              req = c("timestamp", "timestamp_num", "depth"),
              extract_names = colnames,
              type = all)
  check_class(input = archival$timestamp, to_class = "POSIXct", type = "stop")
  check_class(input = archival$depth, to_class = "numeric", type = "stop")
  # Check acoustic centroids have been supplied as a list
  check_class(input = acc_centroids, to_class = "list", type = "stop")
  out$time <- rbind(out$time, data.frame(event = "initial_checks_passed", time = Sys.time()))
  # Check depth error
  de_1 <- calc_depth_error(archival$depth[1])
  if(length(de_1) != 2){
    stop("'calc_depth_error' should be a function that returns a numeric vector of length two (i.e., a lower and upper depth adjustment).")
  }
  if(de_1[1] > 0 | de_1[2] < 0){
    stop("'calc_depth_error' should return a negative and a postive adjustment (in that order).")
  }

  #### Study site rasters
  ## Blank map for space use over the study area
  if(is.null(map)) {
    map <- raster::setValues(bathy, 0)
    map <- raster::mask(map, bathy)
  }


  ######################################
  #### Implement splitting (if necessary)

  #### Define a list of dataframes
  # .. If the algorithm is to be implemented in parallel
  if(inherits(acoustics, "list") | n_cores > 1 | !is.null(split)){

    #### Implement splitting
    if(length_acoustics_tmp == 1){

      ## Define split if not provided
      # If the split hasn't been specified, then the user doesn't care
      # ... about splitting the outputs into biologically interpretable time intervals
      # ... However, we"ll still define a split factor, to be based on computational perspectives
      if(is.null(split)) {
        cat_to_cf("... Splitting 'acoustics' into chunks...")
        chunks <- seq(min(acoustics$timestamp), max(acoustics$timestamp), length.out = n_cores+1)
        dft <- difftime(chunks[2], chunks[1])
        dft_num <- as.numeric(dft)
        dft_num <- floor(dft_num)
        dft_units <- attr(dft, "units")
        message(paste("'acoustics' dataframe split into chunks of ~", dft_num, dft_units, "across", n_cores, "core(s)."))
        split <- paste(dft_num, dft_units)
      }

      ## Split dataframe
      acoustics$split <- lubridate::floor_date(acoustics$timestamp, unit = split)
      acoustics_ls <- split(acoustics, f = acoustics$split)
    } else{
      acoustics_ls <- acoustics
    }

    #### Process split dataframes
    cat_to_cf("... Processing acoustics chunks...")

    ## Remove NULL/length 0 elements
    cat_to_cf("... ... Checking for NULL/empty chunks...")
    empty_elms <- sapply(acoustics_ls, function(x) is.null(x) | nrow(x) == 0)
    if(any(empty_elms)) {
      msg <- paste0("acoustics_ls[c(", paste0(which(empty_elms), collapse = ","), ")] chunks are empty/NULL and will be removed...")
      message(msg)
      cat_to_cf(paste("... ... ...", msg))
      acoustics_ls <- acoustics_ls[which(!empty_elms)]
    }

    ## Force overlapping time series
    # If we naively split the dataframe into number of different windows,
    # ... on every run, we have to stop before the last acoustic reading
    # ... (because we can't identify the next receiver - their isn't one in the split dataframe)
    # ...which means we're not including some information when we estimate space use.
    # To get around this, in the list dataframes, we need to add the first line of every dataframe
    # ... to the previous dataframe. Then, when we split the dataframe, we won't be loosing information
    # ...because we've copied the last line.
    cat_to_cf("... ... Overlapping chunks...")
    acoustics_ls_wth_overlap <-
      lapply(2:(length(acoustics_ls)), function(i){
               # define an adjusted dataframe, binds the previous dataframe
               # ... with the first row of the dataframe in question:
               adj <- rbind(acoustics_ls[[i-1]], acoustics_ls[[i]][1, ])
               return(adj)
             })
    # Add back the final element:
    acoustics_ls_wth_overlap[[length(acoustics_ls)]] <- acoustics_ls[[length(acoustics_ls)]]
    names(acoustics_ls_wth_overlap) <- names(acoustics_ls)

    #### Additional checks
    # Check the number of rows in acoustics_ls_wth_overlap. This cannot be less than two.
    # ... If there are no rows or only one row,
    # ... then we can't calculate where the individual was
    # ... next detected, which will cause problems.
    l <- length(acoustics_ls_wth_overlap)
    lapply(1:l, function(i){
      nrw <- nrow(acoustics_ls_wth_overlap[[i]])
      if(nrw < 2){
        stop(paste("acoustics_ls_wth_overlap[[", i, "]] has less than two rows. This is not allowed."))
      }
    })

    out$time <- rbind(out$time, data.frame(event = "acoustics_chunks_defined", time = Sys.time()))

  } else {
    acoustics_ls_wth_overlap <- acoustics_tmp
  }


  ######################################
  #### Visualise time series

  #### Focus on the data for which we have both acoustic and archival observations
  cat_to_cf("... Processing acoustic and archival time series...")
  movement_ts <- lapply(1:length(acoustics_ls_wth_overlap), function(i){
    acc <- acoustics_ls_wth_overlap[[i]]
    nrw_acc_pre <- nrow(acc)
    nrw_arc_pre <- nrow(archival)
    acc <- acc[acc$timestamp >= min(archival$timestamp) - 2*60 &
                 acc$timestamp <= max(archival$timestamp) + 2*60, ]
    arc <- archival[archival$timestamp >= min(acc$timestamp) - 2*60, ]
    if(i < length(acoustics_ls_wth_overlap)){
      arc <- arc[arc$timestamp <= max(acc$timestamp) + 2*60, ]
    }
    nrw_acc_post <- nrow(acc)
    nrw_arc_post <- nrow(arc)
    nrw_acc_delta <- nrw_acc_pre - nrw_acc_post
    nrw_arc_delta <- nrw_arc_pre - nrw_arc_post
    if(nrw_acc_post == 0 | nrw_arc_post == 0) stop("No overlapping acoustic/archival observations to implement algorithm.")
    if(nrw_acc_delta != 0) {
      cat_to_cf(paste("... ...  Chunk", i, ":", nrw_acc_delta, "acoustic observation(s) beyond the ranges of archival observations ignored."))
      message(paste("Chunk", i, ":", nrw_acc_delta, "acoustic observation(s) beyond the ranges of archival observations ignored."))
      }
    if(nrw_arc_delta != 0) {
      cat_to_cf(paste("... ... Chunk", i, ":", nrw_arc_delta, "archival observation(s) beyond the ranges of (processed) acoustic detections ignored."))
      message(paste("Chunk", i, ":", nrw_arc_delta, "archival observation(s) beyond the ranges of (processed) acoustic detections ignored."))
    }
    ls <- list(acoustics = acc, archival = arc)
    return(ls)
  })
  out$ts_by_chunk <- movement_ts
  out$time <- rbind(out$time, data.frame(event = "movement_time_series_processed", time = Sys.time()))

  #### Visualise processed time series
  if(plot_ts) {
    cat_to_cf("... Plotting acoustic and archival time series (for each chunk)...")
    if(length(movement_ts) < 25) pp <- graphics::par(mfrow = prettyGraphics::par_mf(length(movement_ts)))
    lapply(movement_ts, function(move){
      acoustics <- move$acoustics
      archival  <- move$archival
      axis_ls <- prettyGraphics::pretty_plot(archival$timestamp, abs(archival$depth)*-1,
                                             pretty_axis_args = list(side = 3:2,
                                                                     axis = list(list(format = "%H:%M:%S %d-%m-%y"),
                                                                                 list()
                                                                                 )
                                                                     ),
                                             xlab = "Time stamp", ylab = "Depth (m)",
                                             type = "l")
      prettyGraphics::pretty_line(acoustics$timestamp,
                                  pretty_axis_args = list(axis_ls = axis_ls),
                                  inherit = TRUE,
                                  replace_axis = list(side = 1, pos = axis_ls[[2]]$lim[1]),
                                  add = TRUE,
                                  pch = 21, col = "royalblue", bg = "royalblue")
    })
    if(length(movement_ts) < 25) graphics::par(pp)
    out$time <- rbind(out$time, data.frame(event = "time_series_plotted", time = Sys.time()))
  }


  ######################################
  #### Implement ACDC algorithm

  #### Checks
  n_chunks <- length(acoustics_ls_wth_overlap)
  # Define a list of files, one for each chunk
  if(verbose & con != "") {
    con_ls <- lapply(1:n_chunks, function(i) {
      file <- paste0(con_dir, "dot_acdc_log_", i, ".txt")
      return(file)
    })
  } else {
    con_ls <- lapply(1:n_chunks, function(i) {
      return("")
    })
  }

  # Write blank files to directory if required
  if(verbose & con != "" & n_chunks > 1) {
    cat_to_cf("... Defining chunk-specific log files as dot_acdc_log_1.txt, dot_acdc_log_2.txt etc...")
    lapply(con_ls, function(file) {
      if(!file.exists(file)){
        msg1 <- paste(file, "does not exist: attempting to write file in specified directory...")
        cat_to_cf(paste("... ...", msg1))
        message(msg1)
        file.create(file1 = file)
        cat_to_cf("... ... ... Blank file successfully written to file.")
        message("... Blank file successfully written to file.")
      }
    })
  }

  #### Implement ACDC algorithm directly via .acdc back-end
  if(length(acoustics_ls_wth_overlap) == 1) {

    #### Implement algorithm
    cat_to_cf("... Calling .acdc() to implement ACDC algorithm on one chunk...")
    out$time <- rbind(out$time, data.frame(event = "calling_.acdc", time = Sys.time()))
    .out <- .acdc(acoustics = movement_ts[[1]]$acoustics,
                  archival = movement_ts[[1]]$archival,
                  bathy = bathy,
                  map = map,
                  detection_range = detection_range,
                  mobility = mobility,
                  calc_depth_error = calc_depth_error,
                  acc_centroids = acc_centroids,
                  plot = plot,
                  plot_ts = FALSE,
                  verbose = verbose,
                  con = con,
                  progress = progress,
                  check = FALSE,
                  keep_args = FALSE)

  } else {

    #### Implement algorithm in parallel
    cat_to_cf(paste("... Calling .acdc() to implement ACDC algorithm on", length(acoustics_ls_wth_overlap), "chunks, using", n_cores, "cores..."))
    out$time <- rbind(out$time, data.frame(event = "calling_.acdc", time = Sys.time()))
    .out <- pbapply::pblapply(1:length(acoustics_ls_wth_overlap), cl = cl, function(i){

      #### Implement algorithm
      .out <- .acdc(acoustics = movement_ts[[i]]$acoustics,
                    archival = movement_ts[[i]]$archival,
                    bathy = bathy,
                    map = map,
                    detection_range = detection_range,
                    mobility = mobility,
                    calc_depth_error = calc_depth_error,
                    acc_centroids = acc_centroids,
                    plot = plot,
                    plot_ts = FALSE,
                    verbose = verbose,
                    con = con_ls[[i]],
                    progress = 0L,
                    check = FALSE,
                    keep_args = FALSE)
      return(.out)
    })
    if(!is.null(cl)) parallel::stopCluster(cl = cl)
  }

  #### Return outputs
  out$.acdc <- .out
  t_end <- Sys.time()
  out$time <- rbind(out$time, data.frame(event = "algorithm_competion", time = t_end))
  out$time$serial_duration <- Tools4ETS::serial_difference(out$time$time, units = "mins")
  out$time$total_duration <- NA
  total_duration <- sum(as.numeric(out$time$serial_duration), na.rm = TRUE)
  out$time$total_duration[nrow(out$time)] <- total_duration
  cat_to_cf(paste0("... flapper::acdc() call completed (@ ", t_end, ") after ~", round(total_duration, digits = 2), " minutes."))
  class(out) <- c(class(out), "acdc")
  return(out)

}


######################################
######################################
#### acdc_simplify()

#' @title Simplify the outputs of \code{\link[flapper]{acdc}}
#' @description This function simplifies the output of \code{\link[flapper]{acdc}}, by processing information from the '.acdc' elements of a \code{\link[flapper]{acdc-class}} object that hold the results of calls to the workhorse function \code{\link[flapper]{.acdc}}. This is especially useful if the ACDC algorithm has been applied chunk-wise, in which case the results for each acoustic chunk are returned in a list. The function aggregates information across chunks to generate a continuous time series of results and a map of where the individual could have spend more or less time over the entire time series.
#' @param acdc An \code{\link[flapper]{acdc-class}} object returned by \code{\link[flapper]{acdc}}.
#' @param keep_chunks A logical variable that defines whether or not to retain all chunk-specific information.
#' @param ... Additional arguments (none implemented).
#' @return The function returns an object of class \code{\link[flapper]{.acdc-class}}.
#' @details If the \code{\link[flapper]{acdc}} function was implemented step-wise, this function simply extracts the necessary information and re-packages it into an \code{\link[flapper]{.acdc-class}} object. For a chunk-wise implementation, the function (a) computes the map of where the individual could have spent more or less time by aggregating the chunk-specific maps, accounting for the overlap between chunks; (b) simplifies chunk-specific records into a single contiguous time series, with re-defined time stamps from the start to the end of the time series to (c) return a \code{\link[flapper]{.acdc-class}} object.
#' @seealso The ACDC algorithm is implemented by \code{\link[flapper]{acdc}}, via the back-end function \code{\link[flapper]{.acdc}}. After simplification, \code{\link[flapper]{acdc_plot}} and \code{\link[flapper]{acdc_animate}} can be implemented to visualise time-specific results.
#' @author Edward Lavender
#' @export
#'

acdc_simplify <- function(acdc, keep_chunks = FALSE,...) {

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
        stop("chunk$record[[1]]$spatial[[1]]$map_timestep is NULL. In flapper::acdc(), plot = 1L (or greater/NULL) is required to return the necessary spatial information to correct for overlapping detection time series across chunks in the summation of chunk-specific maps.")
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
    out$map <- raster::calc(out$map, sum)

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

  #### Keep chunk-specific information, if requested
  if(keep_chunks) out$chunks <- acdc$.acdc

  #### Return outputs
  class(out) <- c(class(out), ".acdc")
  return(out)

}

######################################
######################################
#### acdc_plot()

#' @title Plot time-specific maps from the ACDC algorithm
#' @description This function is used to plot time-specific maps from the ACDC algorithm. To implement the function, a \code{\link[flapper]{.acdc-class}} list from \code{\link[flapper]{acdc}} and \code{\link[flapper]{acdc_simplify}} (or \code{\link[flapper]{.acdc}}) must be supplied, from which the results can be extracted and plotted for specified time steps. For each time step, the function extracts the necessary information; sets up a blank background plot using \code{\link[raster]{plot}} and \code{\link[prettyGraphics]{pretty_axis}} and then adds requested spatial layers to this plot. Depending on user-inputs, this will usually show a cumulative map of where the individual could have spent more or less time, summed from the start of the algorithm to each time point. Coastline, receivers and acoustic centroids can be added and customised and the finalised plots can be returned or saved to file.
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
#' @return The function plots the results of the ACDC algorithm at specified time steps, with one plot per time step. Plots are saved to file if \code{png_param} is supplied.
#' @examples
#' #### Step (1): Implement ACDC algorithm
#' # ... see examples via acdc()
#'
#' #### Step (2): Simplify outputs of the ACDC algorithm
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
        stop("'acdc' must be converted from 'acdc' to '.acdc. via flapper::acdc_simpify() before implementing this function.")
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

#' @title Create a html animation of the ACDC algorithm
#' @description This function is a simple wrapper for \code{\link[flapper]{acdc_plot}} and \code{\link[animation]{saveHTML}} which creates an animation of the ACDC algorithm over time. To implement this function, a named list of arguments for \code{\link[flapper]{acdc_plot}}, which creates the plots, must be supplied. This is embedded within \code{\link[animation]{saveHTML}}, which creates a folder in the working directory named `images' that contains a .png file for each time step and an animation as a .html file.
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
