######################################
######################################
#### acs_setup_mobility()

#' @title Examine the constant `mobility' assumption of the AC* algorithm(s)
#' @description In the simplest and fastest (and only) version of the acoustic-container (AC) and acoustic-container depth-contour (ACDC) algorithms currently implemented by \code{\link[flapper]{flapper}}, the rate at which acoustic containers expand and contract depends on a single `mobility' parameter that describes how the set of possible locations for an individual changes through time according to the passive acoustic telemetry data. These changes essentially reflect the maximum horizontal distance that an individual can move in the regularised time steps between sequential detections. However, in some situations, a fixed parameter for horizontal movement may not be well supported. For instance, in cases with archival data (the ACDC algorithm), animals changes dramatically in depth and for which changes in depth are likely to curtail the extent of horizontal movement*. Thus, this function investigates the extent to which the horizontal distance an animal could travel changes through time if the mobility parameter is instead conceptualised as the diagonal distance between sequential depth observations.
#'
#' @param depth A vector of depth observations, whose units match the \code{mobility} parameter, to be incorporated into the ACDC algorithm. (Note that the ACDC algorithm may drop depth observations internally depending on their alignment with the acoustic time series and so, ideally, pre-processed time series should be passed to \code{depth} to ensure inferences correspond directly to the time series modelled by \code{\link[flapper]{acdc}}.) Depth observations should be regularly spaced in time (i.e., represent a time series for a single individual, without gaps).
#' @param mobility A number, in the same units as \code{depth}, that defines the maximum horizontal distance that an individual could move in the time period between archival observations (see \code{\link[flapper]{acs_setup_containers}}).
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
#' *If the horizontal distances that an individual could travel are not very variable, then the benefits of a single mobility parameter are likely to outweigh the costs. On the other hand, substantial variation in the horizontal distances that an individual could travel may suggest that a constant mobility parameter is inappropriate; in this situation, the correct containers could be computed on-the-fly within the ACDC algorithm, but this is not currently implemented. However, the particle filtering algorithms can account for this by the incorporation of movement models within/between acoustic containers.
#'
#' @return The function returns a numeric vector of distances and, if \code{plot = TRUE}, a histogram of those distances.
#'
#' @examples
#' #### Example (1) Explore mobility for single individual
#' mob <-
#'   acs_setup_mobility(depth = dat_archival$depth[dat_archival$individual_id == 25],
#'                      mobility = 200)
#'
#' #### Example (2) Customise histogram
#' # ... suppress plot with plot =  FALSE
#' # ... suppress add_* lists with NULL or customise
#' # ... pass args to prettyGraphics::pretty_hist() via ...
#' mob <-
#'   acs_setup_mobility(depth = dat_archival$depth[dat_archival$individual_id == 25],
#'                       mobility = 200,
#'                       add_mobility = NULL,
#'                       add_rug = list(side = 3, pos = 25000,
#'                                      col = "darkred", lwd = 0.25, ticksize = 0.01),
#'                       breaks = 25)
#'
#' #### Example (3) Explore mobility for mulitple individuals
#' pp <- graphics::par(mfrow = c(2, 2))
#' mob_ls <- lapply(unique(dat_archival$individual_id), function(id){
#'   acs_setup_mobility(depth = dat_archival$depth[dat_archival$individual_id == id],
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
#' # ... between acoustic containers for some applications via particle filtering.
#'
#' @seealso \code{\link[flapper]{acs_setup_mobility}}, \code{\link[flapper]{acs_setup_containers}} and \code{\link[flapper]{acs_setup_detection_kernels}} are used to set up the AC and ACDC algorithms as implemented by \code{\link[flapper]{ac}} and \code{\link[flapper]{acdc}}.
#'
#' @author Edward Lavender
#' @export
#'

acs_setup_mobility <- function(depth,
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
#### acs_setup_containers()

#' @title Setup the detection containers required for the AC* algorithm(s)
#' @description This function produces the detection containers required by the acoustic-container (AC) and acoustic-container depth-contour (ACDC) algorithms.
#' @param xy A \code{\link[sp]{SpatialPointsDataFrame}} object that defines receiver IDs and locations. The \code{data} slot must include a dataframe with a column that provides a unique, integer identifier for each receiver (`receiver_id'). The coordinate reference system should be the Universal Transverse Mercator system with distances in metres (to match \code{detection_range}, see below).
#' @param detection_range A number that defines the maximum detection range (m) at which an individual could be detected from a receiver.
#' @param coastline (optional) A \code{\link[sp]{SpatialPolygonsDataFrame-class}} object that defines the coastline in an area. If provided, detection containers are processed to remove any areas on land. Algorithm speed declines with the complexity of the coastline.
#' @param boundaries (optional) An \code{\link[raster]{extent}} object that defines the boundaries of an area within which individuals are assumed to have remained. If provided, acoustic containers are processed to remain within this area.
#' @param ... Additional arguments passed to \code{\link[flapper]{get_detection_containers}} and, ultimately, \code{\link[rgeos]{gBuffer}}, except \code{byid} which is necessarily \code{TRUE}.
#' @param plot A logical input that defines whether or not to produce a plot of the area, including receivers, the coastline and the area boundaries (if provided), and acoustic containers. This is useful for checking purposes but it can reduce algorithm speed.
#' @param verbose A logical input that defines whether or not to print messages to the console to relay function progress.
#'
#' @details Given a detection at a particular receiver at a particular time, the detection container defines the boundaries of the area around a receiver within which the individual must have been located (from the perspective of that receiver).
#'
#' For the AC* algorithms, note that in some coastal settings the representation of detection containers as \code{\link[sp]{SpatialPolygonsDataFrame-class}} objects may cause a mismatch with detection kernels and the bathymetry \code{\link[raster]{raster}} in terms of what is defined as water versus land. At the time of writing, in the AC* algorithms the detection kernels/bathymetry data take precedence, with any grid cells that have a value of \code{NA} masked, even if within `detection containers'. In the future, these disparities should be resolved by redefining detection containers on the bathymetry grid too.
#'
#' @return The function returns a list of \code{\link[sp]{SpatialPolygonsDataFrame-class}} objects, with one element for all numbers from 1 to the maximum receiver number (\code{xy$receiver_id}). Any list elements that do not correspond to receivers contain a \code{NULL} element. List elements that correspond to receivers contain a \code{\link[sp]{SpatialPolygonsDataFrame-class}} object containing the detection container for that receiver.
#'
#' @examples
#' #### Define data for acs_setup_containers()
#' ## Define coordinates of receivers as SpatialPointsDataFrame with UTM CRS
#' # CRS of receiver locations as recorded in dat_moorings
#' proj_wgs84 <- sp::CRS(SRS_string = "EPSG:4326")
#' # CRS of receiver locations required
#' proj_utm   <- sp::CRS(SRS_string = "EPSG:32629")
#' # Define SpatialPoints object
#' xy_wgs84 <- sp::SpatialPoints(dat_moorings[, c("receiver_long", "receiver_lat")], proj_wgs84)
#' xy_utm <- sp::spTransform(xy_wgs84, proj_utm)
#' # Link with receiver IDs to define a SpatialPointsDataFrame
#' xy_utm <-
#'  sp::SpatialPointsDataFrame(xy_utm,
#'                             dat_moorings[, "receiver_id", drop = FALSE])
#'
#' #### Example (1): Define a list of containers with specified parameters
#' # ... (Argument values are small to reduce computation time for examples)
#' containers <- acs_setup_containers(xy = xy_utm,
#'                                  detection_range = 500
#'                                  )
#' # A list of SpatialPolygonsDataFrames is returned
#' # with elements from 1:max(xy_utm$receiver_id)
#' # NULL elements correspond to numbers in this sequence that do not refer to receivers
#' # Otherwise a SpatialPolygonsDataFrame is returned with all the containers for that receiver
#' containers
#'
#' #### Example (2): Visualise the containers produced via plot = TRUE
#' containers <- acs_setup_containers(xy = xy_utm,
#'                                   detection_range = 500,
#'                                   plot = TRUE
#'                                   )
#'
#' #### Example (3): Remove areas of the containers that overlap with coastline
#' containers <- acs_setup_containers(xy = xy_utm,
#'                                   detection_range = 500,
#'                                   plot = TRUE,
#'                                   coastline = dat_coast
#'                                   )
#'
#' #### Example (4): Remove areas of the containers beyond a boundary
#' xy_utm_coords <- sp::coordinates(xy_utm)
#' boundaries <- raster::extent(min(xy_utm_coords[, 1]),
#'                              max(xy_utm_coords[, 1]),
#'                              min(xy_utm_coords[, 2]),
#'                              max(xy_utm_coords[, 2])
#'                         )
#' containers <- acs_setup_containers(xy = xy_utm,
#'                                   detection_range = 500,
#'                                   plot = TRUE,
#'                                   coastline = dat_coast,
#'                                   boundaries = boundaries
#'                                   )
#'
#' @author Edward Lavender
#' @export
#'

acs_setup_containers <- function(xy,
                                detection_range,
                                coastline = NULL,
                                boundaries = NULL,
                                plot = FALSE,
                                verbose = TRUE,...){

  #### Initiate function
  t_onset <- Sys.time()
  cat_to_console <- function(..., show = verbose) if(show) cat(paste(..., "\n"))
  cat_to_console(paste0("flapper::acs_setup_detection_containers() called (@ ", t_onset, ")..."))

  #### Function checks
  cat_to_console("... Checking user inputs...")
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
  if(!all(rownames(xy@data) == as.character(rs))){
    warning("'rownames(xy@data) overwritten with xy$receiver_id.", immediate. = TRUE, call. = FALSE)
  }
  if(!is.null(coastline)){
    check_class(input = coastline, to_class = "SpatialPolygonsDataFrame")
    if(length(coastline) != 1)
      stop("'coastline' has multiple features which need to be dissolved into a single feature SpatialPolgyonsDataFrame.")
  }
  if(!is.null(coastline) & !is.null(boundaries)) {
    coastline <- raster::crop(coastline, boundaries)
    if(is.null(coastline)) message("No coastline within defined boundaries. \n")
  }
  check...("byid",...)

  #### Plot map of area
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

  #### Make containers
  cat_to_console("... Making containers...")
  containers <- get_detection_containers(xy = xy,
                                       detection_range = detection_range,
                                       boundaries = boundaries,
                                       coastline = coastline,
                                       plot = FALSE,
                                       byid = TRUE,...)
  containers <- sp::spChFIDs(containers, as.character(rs))
  xyd <- data.frame(xy)
  rownames(xyd) <- as.character(xy$receiver_id)
  containers <- sp::SpatialPolygonsDataFrame(containers, xyd)

  #### Add containers to map
  if(plot){
    cat_to_console("... Plotting containers on map...")
    raster::lines(containers, col = "royalblue")
  }

  #### Define a list, with one element for each container
  cat_to_console("... Processing containers...")
  containers <- lapply(1:max(rs), function(i){
    if(i %in% rs){
      out <- containers[containers$receiver_id == i, ]
    } else {
      out <- NULL
    }
    return(out)
  })

  #### Return containers
  t_end <- Sys.time()
  total_duration <- difftime(t_end, t_onset, units = "mins")
  cat_to_console(paste0("... flapper::acs_setup_detection_containers() call completed (@ ", t_end, ") after ~", round(total_duration, digits = 2), " minutes."))
  return(containers)
}


######################################
######################################
#### acs_setup_detection_kernels()

#' @title Setup detection probability kernels for the AC* algorithm(s)
#' @description This function produces detection probability kernels for incorporation into the acoustic-container* (AC*) algorithms. Within acoustic containers, the incorporation of detection probability kernels reduces uncertainty and increases precision in simulated patterns of space use by up-weighting areas nearer to receivers when an individual is detected and down-weighing areas nearer to receivers when an individual is not detected (see Details).
#'
#' To implement the function, a \code{\link[sp]{SpatialPointsDataFrame}} that defines receiver IDs, locations and deployment dates must be supplied (via \code{xy}). A record of servicing events for receivers can also be supplied (via \code{services}). Detection probability kernels are calculated around each receiver, using a user-defined function based on Euclidean distances (\code{calc_detection_pr}) across a \code{\link[raster]{raster}} (\code{bathy}). Kernels can be restricted by barriers to movement as defined by \code{NAs} in \code{bathy} and the boundaries of the area. These kernels are used to weight possible locations around a receiver when an individual is detected.
#'
#' For each unique array design (i.e. set of active receivers, given receiver deployment dates and servicing events, if applicable), a detection probability surface across the whole area is also created, which is used to weight possible locations of the individual in the time steps between detections (up-weighting locations away from receivers). By default, these calculations account for any areas of overlap in the detection probability kernels of multiple receivers. This step is computationally demanding, but it can be suppressed or sped-up via the \code{overlaps} argument. Outputs are returned in a named list that is designed to be incorporated into the AC* algorithm(s).
#'
#' @param xy A \code{\link[sp]{SpatialPointsDataFrame}} that defines receiver IDs, locations and deployment dates. The \code{data} slot must include a dataframe with the following columns: an unique, integer identifier for each receiver (`receiver_id') and receiver deployment \code{\link[base]{Dates}} (`receiver_start_date' and `receiver_end_date'). For receiver locations, the coordinate reference system should be the Universal Transverse Mercator system with distances in metres (as for \code{bathy}, below).
#' @param services (optional) A dataframe that defines receiver IDs and servicing \code{\link[base]{Dates}} (times during the deployment period of a receiver when it was not active due to servicing). If provided, this must contain the following columns: an integer identifier for serviced receivers (named ‘receiver_id’) and two columns that define the time of the service(s) (‘service_start_date’ and ‘service_end_date’) (see \code{\link[flapper]{make_matrix_receivers}}).
#' @param containers The list of detection containers, with one element for each number from \code{1:max(xy$receiver_id)}, from \code{\link[flapper]{acs_setup_containers}}.
#' @param overlaps (optional) A named list, from \code{\link[flapper]{get_detection_containers_overlap}}, that defines, for each receiver, for each day over its deployment period, whether or not its detection container overlapped with those of other receivers. If provided, this speeds up detection probability calculations in overlapping regions by focusing on only the subset of receivers with overlapping detection probability kernels. If there are no overlapping receivers, \code{FALSE} can be supplied instead to suppress these calculations.
#' @param calc_detection_pr A function that takes in a vector of distances and returns a vector of detection probabilities (around a receiver). Detection probability should decline to 0 after the \code{detection_range} distance from a receiver (see \code{\link[flapper]{acs_setup_containers}}).
#' @param bathy A \code{\link[raster]{raster}} of the bathymetry across an area (see \code{\link[flapper]{ac}}). Receiver locations and detection probability kernels are represented across this \code{\link[raster]{raster}}. As for \code{xy}, the coordinate reference system should be the Universal Transverse Mercator system with distances in metres. Resolution needs to be sufficiently high such that detection probability can be represented effectively, given \code{calc_detection_pr}, but sufficiently low for convenient run times. If \code{bathy} contains NAs, it is also taken as a mask that is applied to detection probability kernels to remove `impossible` areas (e.g., on land).
#' @param verbose A logical input that defines whether or not to print messages to the console to relay function progress.
#'
#' @details A detection probability kernel is a bivariate probability density function that describes detection probability around a receiver. Typically, this takes the shape of a dome whereby detection probability declines uniformly around a receiver with increasing distance. Accordingly, this function assumes that detection probability kernels only depend on distance (via \code{calc_detection_pr}) and are constant in space and time. Spatially variable detection probability kernels can be incorporated just as easily into the AC* algorithm(s) but need to be created manually. For example, in areas of complex coastline, narrow peninsulas that punctuate detection containers may effectively block transmissions from the outer regions of detection containers from being detected at receivers, in which case a model that incorporates the `line of sight' between receivers and the surrounding regions may be appropriate. However, temporally variable detection probability kernels are not currently supported by the AC* algorithm(s).
#'
#' The purpose of detection probability kernels within the AC* algorithm(s) is to reduce the uncertainty of the possible positions of an individual within the acoustic containers, both when an individual is detected and when it is not. When an individual is detected at a receiver, a central assumption of the AC* algorithm(s) is that the individual is within a finite radius of that receiver in an area defined as the ‘detection container’. Under the simplest implementation of the algorithm, this represents a threshold detection probability model in which the probability of detection is certain within the container and zero outside. Under this model, the individual is equally likely to have been in any location within the container and all possible locations of the individual are weighted equally. This approach is suitable for the most conservative analyses, since even unlikely positions of the individual receive equal weighting. However, detection probability typically declines with distance around a receiver and the incorporation of this information in the form of kernels around the receiver at which is detected (and any receivers with overlapping kernels) improve precision by increasing the weighting for some areas over others. The way in which this increase in precision is realised over space depends on whether or not detection probability kernels and the timing of detections at multiple receivers overlap.
#'
#' In the simplest scenario, an individual is detected at a receiver whose container does not overlap with any other receiver and the kernel simply increases the weight of locations nearest to the receiver (according to a user-specified detection probability function). For some array designs, an alternative possibility is that an individual is detected at a receiver whose container overlaps with other receiver(s). In this case, detection or lack of detection at the same moment at the receivers with overlapping containers provides further information on the location of the individual. One the one hand, if the individual is not detected at the other receiver(s), then the probability that the individual is in the overlapping region(s) is reduced in line with the overlap in detection probability, which decreases the weight of potential locations in these areas. On the other hand, if the individual is detected at effectively the same time at multiple receiver(s), then the individual is more likely to be within the overlapping parts of their containers, and the weight of possible locations here is correspondingly elevated. The definition of ‘effectively’ is context specific but depends on the accuracy with which the clocks of different receivers are synchronised. (For example, ± 5 s might be reasonable.) In these situations, the detection container essentially just acts as a computational device by restricting calculations to the container(s) and ‘cutting’ the probability of detection to zero beyond this area. This process follows the standard rules of probability.
#'
#' When an individual is not detected, the detection container grows into a set of ‘acoustic containers’ that contain the set of possible locations for the individual. As they grow, they may encompass other detection containers before they shrink towards the receiver at which the individual is next detected. During this time, the AC* algorithm(s) identify the possible locations of the individual within these areas. Under the most conservative approach, at each time step all positions are treated as equally likely (although normalisation within the algorithm(s) effectively down-weights time steps in which the location of the individual is more uncertain). This includes any positions within the detection containers of other receivers since, under realistic conditions, there is usually a non-zero probability that an individual can be near a receiver and yet remain undetected. However, this is typically unlikely and when the goal of the analysis is to create more precise estimates of space use, the incorporation of detection probability kernel(s) around receivers effectively reduces the probability that the individual is within their detection containers during this time. Again, to incorporate detection probability kernels in this way, it is necessary to account for overlapping detection ranges, where the probability of detection is higher and, therefore, the probability of a possible location in such an area is lower given the absence of a detection. As above, this process follows the laws of probability.
#'
#' In summary, in the AC* algorithm(s), detection and acoustic containers describe the spatial extent of our uncertainty when an individual is detected and in the time between detections respectively. The purpose of this function is to pre-process and package the information provided by detection probability kernels in such a way as to facilitate its incorporation into the AC* algorithm(s). This improves the precision of simulated patterns of space use by down- or up-weighting areas according to a model of detection probability. This provides a reasonable overall assessment of the places in which an individual could have spend more or less time over a period of study.
#'
#' However, it is worth noting that, within acoustic containers, unrealistic areas may be highlighted in which an individual could not have been located because of movement constraints which are not captured by container-level expansion and contraction. For example, there may be parts of a container that are highly unlikely given a detection at a nearby receiver because they are too far away from the receiver at which the individual was subsequently detected. In the AC*PF algorithm(s), the incorporation of a movement model further reduces uncertainty by accounting for the constraints imposed by an individual’s current position on its next position.
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
#' proj_wgs84 <- sp::CRS(SRS_string = "EPSG:4326")
#' proj_utm   <- sp::CRS(SRS_string = "EPSG:32629")
#' xy <- sp::SpatialPoints(moorings[, c("receiver_long", "receiver_lat")],
#'                         proj_wgs84)
#' xy <- sp::spTransform(xy, proj_utm)
#' # Link receiver IDs, locations and deployment dates to form a SpatialPointsDataFrame
#' # ... Note required column names and class types.
#' xy <- sp::SpatialPointsDataFrame(xy, moorings[, c("receiver_id",
#'                                                   "receiver_start_date",
#'                                                   "receiver_end_date")])
#'
#' ## Define bathymetric map of area for which AC* algorithm(s) will be implemented
#' # ... The resolution must be sufficiently high
#' # ... such that there are areas with non zero detection probability.
#' # ... However, function speed will fall with large, high resolution rasters.
#' # ... Here, we set a low resolution for example speed.
#' surface <- raster::raster(raster::extent(dat_gebco), res = c(50, 50))
#' surface <- raster::resample(dat_gebco, surface)
#'
#' ## Define detection probability function
#' # This should depend on distance alone
#' # Detection probability should decline to 0 after detection_range
#' # ... (defined in flapper::acs_setup_containers()).
#' # ... Here, we assume detection_range = 425 m.
#' calc_dpr <-
#'   function(x){
#'     ifelse(x <= 425, stats::plogis(2.5 + -0.02 * x), 0)
#'   }
#' plot(0:1000, calc_dpr(0:1000), type = "l")
#'
#' ## Get detection containers and, if applicable, information on their overlap(s)
#' # We'll use the example containers provided in dat_containers
#' # We'll get their overlaps via flapper::get_detection_containers_overlap
#' overlaps <-
#'   get_detection_containers_overlap(
#'     containers = get_detection_containers(xy = xy, byid = TRUE)
#'   )
#'
#' #### Example (1): Implement function using default options
#' kernels <- acs_setup_detection_kernels(xy = xy,
#'                                        containers = dat_containers,
#'                                        overlaps = overlaps,
#'                                        calc_detection_pr = calc_dpr,
#'                                        bathy = surface)
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
#' @seealso This is one of a number of functions used to set up the AC and ACDC algorithms implemented by \code{\link[flapper]{ac}} and \code{\link[flapper]{acdc}}: \code{\link[flapper]{acs_setup_mobility}}, \code{\link[flapper]{acs_setup_containers}} and \code{\link[flapper]{acs_setup_detection_kernels}}. This function is supported by \code{\link[flapper]{make_matrix_receivers}}, which defines receiver activity statuses; \code{\link[flapper]{acs_setup_containers}}, which defines acoustic containers; and \code{\link[flapper]{get_detection_containers_overlap}}, which defines detection container overlaps
#' @author Edward Lavender
#' @export
#'

acs_setup_detection_kernels <-
  function(xy, services = NULL,
           containers,
           overlaps = NULL,
           calc_detection_pr,
           bathy,
           verbose = TRUE
  ){


    ######################################
    #### Initiation, set up and checks

    #### Initiation
    cat_to_console <- function(..., show = verbose){
      if(show) cat(paste(..., "\n"))
    }
    t_onset <- Sys.time()
    cat_to_console(paste0("flapper::acs_setup_detection_kernels() called (@ ", t_onset, ")..."))
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
    if(any(duplicated(moorings$receiver_id)))
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
    map           <- raster::setValues(bathy, 0)
    if(length(raster::Which(is.na(bathy), cells = TRUE)) > 0L)
      mask_layer <- TRUE
    else
      mask_layer <- FALSE


    ######################################
    ####  Receiver-specific kernels (for detection)

    #### Calculate detection Pr around each receiver
    # (used to up-weight areas around a receiver with a detection)
    cat_to_console("... Getting receiver-specific kernels (for detection)...")
    xy_ls <- lapply(1:length(xy), function(i) xy[i, ])
    detection_kernels_by_xy <-
      pbapply::pblapply(xy_ls, cl = NULL, function(xyi){
        ## Focus on area within container, for speed
        # xyi <- xy_ls[[1]]
        cat_to_console(paste0("\n... ... For receiver ", xyi$receiver_id, " ..."))
        cat_to_console("... ... ... Isolating detection container ...")
        xyi_container <- containers[[xyi$receiver_id]]
        xyi_map <- raster::crop(map, xyi_container)
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
        det_pr_around_xyi <- raster::extend(det_pr_around_xyi, raster::extent(map), value = 0)
        if(mask_layer) det_pr_around_xyi <- raster::mask(det_pr_around_xyi, bathy)
        dpr_at_xyi <- raster::extract(det_pr_around_xyi, xyi)
        if(is.na(dpr_at_xyi)) {
          warning("Detection probability is NA at receiver ", xyi$receiver_id, ".",
                  immediate. = TRUE, call. = FALSE)
        } else if(dpr_at_xyi == 0) warning("Detection probability = NA at receiver ", xyi$receiver_id,
                                           immediate. = TRUE, call. = FALSE)
        ## Create modified kernels
        # ... (a) restricted to where detection Pr is positive (detection container)
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
    cat_to_console(paste0("... flapper::acs_setup_detection_kernels() call completed (@ ", t_end, ") after ~", round(total_duration, digits = 2), " minutes."))
    # Return outputs
    return(out)

  }
