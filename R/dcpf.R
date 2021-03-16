########################################
########################################
#### dcpf_setup_movement_pr()

#' @title A simple movement model dependent on distance
#' @description This function provides a simple movement model that calculates the probability of movement between two locations according to the distance between them, using an logistic equation with pre-defined parameters.
#'
#' @param distance A numeric vector of distances (m).
#'
#' @details Under this model, for distance(s) \eqn{ \leq 500 } m, \eqn{Pr(distance) = logistic(10 + distance  -0.05)}; otherwise, \eqn{Pr(distance) = 0}. This particular model is designed for flapper skate (\emph{Dipturus intermedius}) and represents a reasonable model for the probability of moving a given distance in a two-minute period (in the absence of additional information).
#'
#' @return The function returns a numeric vector of probabilities that represent the probability of movement between two or more areas given the distances between them.
#'
#' @examples
#' pr <- dcpf_setup_movement_pr(1:1000)
#' plot(pr, type = "l", xlab = "Distance (m)", ylab = "Pr(distance)")
#' @seealso This function is used as the default movement model in \code{\link[flapper]{dcpf}}.
#' @author Edward Lavender
#' @export

dcpf_setup_movement_pr <- function(distance) {
  pr <- stats::plogis(10 + distance * -0.05)
  pr[distance > 500] <- 0
  return(pr)
}


########################################
########################################
#### dcpf()

#' @title The depth-contour particle filtering (DCPF) algorithm
#' @description This function implements the depth-contour particle filtering (DCPF) algorithm. This is an extension of the DC algorithm (\code{\link[flapper]{dc}}) that implements particle filtering to reconstruct possible movement paths of an individual (i.e., a benthic animal) over a surface (i.e., the seabed). As in the DC algorithm, at each time step the possible locations of an individual are determined from its depth and the bathymetric landscape (plus a measurement error term). The extension is the incorporation of a  movement model, via a simulation-based particle filtering process, that connects a subset of these locations between time steps into movement paths.
#'
#' To implement this approach, a numeric vector of depth observations resulting from movement over a surface (\code{archival}), as well as a \code{\link[raster]{raster}} of the surface over which movement occurred (\code{bathy}) and a measurement error parameter (\code{depth_error}), must be supplied. A starting location (\code{origin}) can be supplied to constrain the set of possible locations of the individual within this area. At each time step, \code{n} possible locations (`particles') are sampled (with replacement) from the set of possible locations. For each (\code{1:n}) particle, a movement model is used to simulate where the individual could have moved to at the next time step, if it was in any of those locations. In the current framework, the probability of movement into surrounding cells depends on the distance to those cells, which can be represented as using Euclidean or least-cost distances depending on the distance method (\code{calc_distance}), and a user-defined movement model (\code{calc_movement_pr}) that links distances to movement probabilities.
#'
#' At each subsequent time step, this process repeats, with \code{n} possible locations of the individual sampled according to the probability that the individual could have been in that cell, given a previously sampled location, its depth and the bathymetry. The result is a set of paths over the surface that are consistent with the data and model parameters.
#'
#' @param archival A numeric vector of depth (m) observations. Depth should be recorded using absolute values in the same units as the bathymetry (\code{bathy}, see below). Observations are assumed to have been made at regular time intervals.
#' @param bathy A \code{\link[raster]{raster}} of the bathymetry in the area within which the individual was located over the study. Bathymetry values should be recorded as absolute values and in the same units (m) as for depths (see \code{archival}). The coordinate reference system should be the Universal Transverse Mercator projection. The resolution of this layer needs to be sufficiently high such that an individual could transition between cells in the duration between \code{archival} observations (see \code{calc_movement_pr}). If the `shortest distances' method is used for distance calculations (i.e., \code{calc_distance = "lcp"}, see below), then the resolution of the surface in x and y directions should also be equal (for surface's with unequal horizontal resolution, \code{\link[raster]{resample}} can be used to equalise resolution: see \code{\link[flapper]{lcp_over_surface}}). For computational efficiency, it is beneficial if \code{bathy} is cropped to the smallest possible area, with any areas in which movement is impossible (e.g., on land for benthic animals) set to NA (see \code{\link[raster]{crop}}, \code{\link[raster]{mask}} and the processing implemented by \code{\link[flapper]{lcp_over_surface}}).
#' @param origin (optional) A matrix that defines the coordinates (x, y) of the individual's initial location. Coordinates should follow the restrictions for \code{bathy} and lie on a plane (i.e., Universal Transverse Mercator projection).
#' @param depth_error A number that defines the interval around each depth (m) observation that defines the range of depths on the bathymetry raster (\code{bathy}) that the individual could plausibly have occupied at that time. For example, \code{depth_error = 2.5} m implies that the individual could have occupied bathymetric cells whose depth lies within the interval defined by the observed depth (m) +/- 2.5 m. The appropriate value for depth_error depends on measurement error for the \code{archival} and bathymetry (\code{bathy}) data, as well as the tidal range (m) across the area (over the duration of observations).
#' @param calc_distance A character that defines the method used to calculate distances between a point (i.e., a sampled location) and the surrounding cells. This drives the probability of movement into those cells via a movement model (see \code{calc_movement_pr}). Currently supported options are Euclidean distances (\code{"euclid"}) or shortest (least-cost) distances (\code{"lcp"}), which represent the shortest distance that an individual would have to move over the surface (\code{bathy}) to traverse between locations (accounting for both planar and vertical distances). Note that this option requires that resolution of \code{bathy} in the x and y directions is equal. At small spatial scales, this option provides more realistic distances in hilly landscapes, but it is more computationally expensive. At larger scales, horizontal distances tend to overwhelm vertical distances, so Euclidean distances may be acceptable. A pragmatic option is to implement the algorithm (possibly for a subset of the \code{archival} time series) using the Euclidean distances method and then interpolate least-cost paths between sequential positions returned by this approach (see \code{\link[flapper]{lcp_interp}}). This approach will demonstrate whether sequential positions are plausible (i.e., not too far apart) once the bathymetry is taken into account. If so, the shortest-distances derived using this method can then be used for post-hoc adjustment of movement probabilities. Alternatively, this approach may demonstrate that the algorithm should be re-implemented using the shortest distances method (see \code{\link[flapper]{lcp_interp}}).
#' @param calc_movement_pr The movement model. Currently, the only supported option is a function that calculates the probability of movement between two locations in the time between depth observations, given the distance between them. The default option is a declining logistic curve, designed for flapper skate (\emph{Dipturus intermedius}), representing a high probability of movement between nearby locations and a lower probability of movement between distant locations (see \code{\link[flapper]{dcpf_setup_movement_pr}}). For computational efficiency, it is beneficial if the probability of movement is set to zero beyond a distance considered to be highly unlikely, because such locations are excluded from subsequent calculations (see \code{\link[flapper]{dcpf_setup_movement_pr}} and the \code{mobility} argument, below).
#' @param mobility (optional) A number that defines the maximum horizontal distance (m) that the individual could travel in the time period between \code{archival} observations. While this is optional, it is usually computationally beneficial to define \code{mobility} because this restricts distance and movement probability calculations at each time step within the smallest appropriate range (rather than across the full surface).
#' @param n An integer that defines the number of particles (i.e., the number of locations sampled at each time step from the set of possible locations at that time step).
#' @param cl,use_all_cores Parallelisation options. The algorithm can be parallelised within time steps over (1) paths via \code{cl} (a cluster object created by \code{\link[parallel]{makeCluster}}) or (2) within paths for the calculation of shortest distances (if \code{calc_distance = "lcp"}) via a logical input (\code{TRUE}) to \code{use_all_cores}. The most efficient solution is context-specific.
#' @param seed (optional) An integer to define the seed for reproducible simulations (see \code{\link[base]{set.seed}}).
#' @param verbose A logical variable that defines whether or not to print messages to the console to relay function progress
#' @param ... Additional arguments (none implemented).
#'
#' @details
#' \subsection{Background}{The DCPF algorithm simulates possible movement paths of a benthic animal over the seabed, given a regular sequence of depth observations (\code{archival}), the bathymetry (\code{bathy}) over which movement occurred and a movement model (\code{calc_movement_pr}) that specifies the probability of movement from a given location to any other, given the distance between them. The function was motivated by small scale applications concerning the reconstruction of possible movement paths of flapper skate (\emph{Dipturus intermedius}) tagged with archival tags, following capture and release in a given location, for short-periods of time post-release.}
#'
#' \subsection{Methods}{At the first time step, the function identifies all of the locations in which the animal could have been located, based on its depth, the bathymetry and some measurement error (\code{depth_error}) that depends on the accuracy of the depth observations, the bathymetry data and the magnitude of the tidal range over the period of observations. From this set of possible locations, \code{n} starting points (particles) are selected. If an \code{origin} is specified, this selection can be biased towards cells near the origin by the movement model. (Location probability could also be weighted by the proximity between the observed depth and the depths of cells within the range defined by the observed depth and the measurement error (i.e., \code{archival[1] - depth_error, archival[1] + depth error}), but this is not currently implemented.)
#'
#' From each starting position, the Euclidean or shortest distances to cells of the appropriate depth at the next time step are calculated and passed to a movement model that assigns movement probabilities to each cell. While movement probabilities are likely to be behaviourally dependent, time-varying movement parameters are not currently learnt from data or implemented. As such, the movement model will typically depend on the maximum distance that an individual could travel within the time period between depth observations. Across all particles, cells within the required depth range with a movement probability of more than zero from the set of locations from which \code{n} particles are re-sampled, with replacement and according to their probability, at the next time step. This process repeats until the end of the time series. However, note that in the current implementation of the algorithm, unlike for the start of this process, there is no constraint that forces the individual to return to a specified geographical location at the end of the time series. Indeed, other than an (optional) \code{origin} and the information provided by the depth time series (\code{archival}), geographic restrictions (e.g., from acoustic detections) on the location of the animal over the period of observations cannot be incorporated. Therefore, this algorithm is best-suited to small scale applications (e.g., to examine the movements of individuals tagged with archival for short periods of time immediately post-release). If geographical observations (i.e., detections at acoustic receivers) are available for an individual over its time at liberty, the ACDCPF algorithm (currently unavailable) is required to integrate this information into the construction of movement paths.
#'
#' The result is a set of simulated pathways over a surface that are consistent with the data and the model parameters. While the number of particles is predetermined by \code{n}, more than \code{n} possible pathways may be defined by all of the combinations of sequential particle movements. Indeed, if there are vast numbers of possible paths through the landscape captured by particle movements, the function can run into vector memory limitations when assembling paths. Therefore, it is advisable to initiate the algorithm with a small number of paths to get a sense of how many possible paths there are, before increasing the number of paths incrementally. In the future, this may be resolved with a separate routine for processing paths (e.g., \code{dcpf_simplify}, like \code{\link[flapper]{acdc_simplify}}), but this is not currently implemented.}
#'
#' \subsection{Convergence}{Algorithm convergence is not guaranteed. There are four main circumstances in which the algorithm may fail to return any paths that span the start to the end of the depth time series:
#'  \enumerate{
#'    \item \strong{Chance.} All \code{n} paths may be `dead ends'. This possibility can be mitigated by increasing \code{n}.
#'    \item \strong{Movement model.} The movement model may be too limiting. This possibility can be mitigated by ensuring that the movement model realistically represents the probability that an individual can transition between cells given the distance between them.
#'    \item \strong{Depth error.} The depth error may be too restrictive, given the accuracy of the depth observations, the bathymetry data and the tidal height across an area. This possibility can be mitigated by ensuring that the depth error is appropriate.
#'    \item \strong{Other assumptions.} Other assumptions (e.g., strict benthic habit) may be violated.
#'  }
#' In these scenarios, the function returns a message that it is about to fail and a list, with one element for each time step, that records the sampled locations and their associated probabilities for each time step from the start of the algorithm until the current time step, before stopping.
#' }
#'
#' @examples
#' #### Define data
#'
#' ## Sample species
#' # In this example, we consider flapper skate (Dipturus intermedius)
#' # ... off the west coast of Scotland.
#'
#' ## Define a starting location (optional) in UTM coordinates
#' xy <- matrix(c(708886.3, 6254404), ncol = 2)
#' ## Define 'observed' depth time series using absolute values
#' # Imagine these are observations made every two minutes
#' depth <- c(163.06, 159.71, 153.49, 147.04, 139.86, 127.19, 114.75,
#'            99.44,  87.01,  78.16,  70.03,  60.23,  49.96,  35.39,
#'            27.75,  20.13,  12.73,  11.32)
#'
#' ## Define surface over which movement occurred
#' # We will use the example dat_gebco bathymetry dataset
#' # This is relatively coarse in resolution, so we need to re-sample
#' # ... the raster to generate a finer-resolution raster such that
#' # ... our animal can transition between cells
#' # ... in the duration between depth observations. For speed,
#' # ... we will focus on a small area around the origin. We could
#' # ... also process the raster in other ways (e.g., mask any areas of land)
#' # ... to improve efficiency.
#' surface    <- dat_gebco
#' boundaries <- raster::extent(707884.6, 709884.6, 6253404, 6255404)
#' blank      <- raster::raster(boundaries, res = c(5, 5))
#' surface    <- raster::resample(surface, blank)
#'
#' ## Define movement model
#' # The default movement model is suitable, with skate moving typically
#' # ... less than 200 m in a two-minute period.
#'
#' ## Visualise movement surface, with starting location overlaid
#' prettyGraphics::pretty_map(add_rasters = list(x = surface),
#'                            add_points = list(x = xy),
#'                            verbose = FALSE)
#'
#' #### Example (1): Implement algorithm using default options
#' # Note that because the bathymetry data is very coarse, and the bathymetry is
#' # ... actually very complex in this region of Scotland, we have to
#' # ... force the depth_error to be high in this example.
#' paths <- dcpf(archival = depth,
#'               bathy = surface,
#'               origin  = xy,
#'               depth_error = 30,
#'               calc_distance = "euclid",
#'               calc_movement_pr = dcpf_setup_movement_pr,
#'               n = 10L,
#'               seed = 1)
#' # The function returns a dataframe with information for each path
#' utils::str(paths)
#' # Algorithm duration during testing ~0.21 minutes
#'
#' \dontrun{
#'
#' #### Example (2): Implement a blanket mobility restriction
#' # This can improve computational efficiency but offers no improvement
#' # ... in this example (e.g., due to relatively small raster, so the cost of
#' # ... focusing in on a specific area matches the speed gains of
#' # ... implementing the distance/movement
#' # ... pr calculations over a smaller area at each time step)
#' paths <- dcpf(archival = depth,
#'               bathy = surface,
#'               origin  = xy,
#'               depth_error = 30,
#'               calc_distance = "euclid",
#'               calc_movement_pr = dcpf_setup_movement_pr,
#'               mobility = 200,
#'               n = 10L,
#'               seed = 1)
#' # Algorithm duration during testing ~0.22 minutes
#'
#' #### Example (3): Implement algorithm using shortest distances
#' # Note the need for a surface with equal resolution if this option is implemented.
#' paths <- dcpf(archival = depth,
#'               bathy = surface,
#'               origin  = xy,
#'               depth_error = 30,
#'               calc_distance = "lcp",
#'               calc_movement_pr = dcpf_setup_movement_pr,
#'               n = 10L,
#'               seed = 1)
#' # This option is slower: algorithm duration during testing ~0.73 minutes
#'
#' #### Example (4): Implement algorithm using shortest distances with mobility restriction
#' paths <- dcpf(archival = depth,
#'               bathy = surface,
#'               origin  = xy,
#'               depth_error = 30,
#'               calc_distance = "lcp",
#'               calc_movement_pr = dcpf_setup_movement_pr,
#'               mobility = 200,
#'               n = 10L,
#'               seed = 1)
#' # Algorithm duration during testing ~0.63 minutes
#' # With shortest distances, the mobility restriction makes more difference
#' # ... in improving the computation time, because calculations are more involved.
#'
#' #### Example (5): Parallelisation for Euclidean distances is via cl argument
#' # ... which implements parallelisation across paths
#' paths <- dcpf(archival = depth,
#'               bathy = surface,
#'               origin  = xy,
#'               depth_error = 30,
#'               calc_distance = "euclid",
#'               calc_movement_pr = dcpf_setup_movement_pr,
#'               mobility = 200,
#'               n = 10L,
#'               cl = parallel::makeCluster(2L),
#'               seed = 1)
#'
#' #### Example (6): Parallelisation for shortest distances is usually best
#' # ... via use_all_cores = TRUE
#' paths <- dcpf(archival = depth,
#'               bathy = surface,
#'               origin  = xy,
#'               depth_error = 30,
#'               calc_distance = "lcp",
#'               calc_movement_pr = dcpf_setup_movement_pr,
#'               mobility = 200,
#'               n = 10L,
#'               use_all_cores = TRUE,
#'               seed = 1)
#' # But the speed benefits in this case are minimal.
#' # Algorithm duration during testing ~0.61 minutes.
#'
#' }
#'
#' @return The function returns a named list that records the parameters used to generate function outputs (`args'), in a named list, and a dataframe (`dcpf') that records possible movement paths over a surface. The dataframe includes the following columns: a unique identifier for each path (`path_id'); the time step (`time_step'), the location (cell ID and three-dimensional coordinates) on \code{bathy} (`cell_id', `cell_x', `cell_y' and `cell_z') and the probability associated with that cell, given movement from the previous cell (`cell_pr'). Rows are ordered by path and then time step. \code{\link[flapper]{dat_dcpf}} provides an example.
#'
#' @seealso For the movement model, Euclidean distances are obtained from \code{\link[raster]{distanceFromPoints}} or shortest distances are obtained from \code{\link[flapper]{lcp_costs}} and \code{\link[flapper]{lcp_from_point}}. The default movement model applied to these distances is \code{\link[flapper]{dcpf_setup_movement_pr}}. Paths between sequential locations can be interpolated via \code{\link[flapper]{lcp_interp}}, which is a wrapper for the \code{\link[flapper]{lcp_over_surface}} routine. This can be useful for checking whether the faster Euclidean distances method is acceptable and, if so, for post-hoc adjustments of movement probabilities based on shortest distances (see \code{\link[flapper]{lcp_interp}}). The results of the algorithm can be visualised with \code{\link[flapper]{dcpf_plot_1d}}, \code{\link[flapper]{dcpf_plot_2d}} and \code{\link[flapper]{dcpf_plot_3d}}. The log-likelihood of the paths, given the movement model, can be calculated via \code{\link[flapper]{dcpf_loglik}}. In terms of related algorithms, the DC algorithm is implemented via \code{\link[flapper]{dc}}. This is extended by the ACDC algorithm, through the integration of locational information from acoustic detections, via \code{\link[flapper]{acdc}}. For depth time series accompanied by acoustic detections, movement paths can be reconstructed via the ACDCPF algorithm (currently unavailable).
#'
#' @author Edward Lavender
#' @export

dcpf <- function(archival,
                 bathy,
                 origin = NULL,
                 depth_error = 2.5,
                 calc_distance = c("euclid", "lcp"),
                 calc_movement_pr = dcpf_setup_movement_pr,
                 mobility = NULL,
                 n = 10L,
                 cl = NULL, use_all_cores = FALSE,
                 seed = NULL,
                 verbose = TRUE,...){

  #### Set up function
  # Function onset
  cat_to_cf <- function(..., show = verbose){
    if(show) cat(paste(..., "\n"))
  }
  t_onset <- Sys.time()
  cat_to_cf(paste0("flapper::dcpf() called (@ ", t_onset, ")..."))
  cat_to_cf("... Setting up function...")
  # List for outputs
  out_dcpf <- list(dcpf = NULL,
                   args = list(archival= archival,
                               bathy = bathy,
                               origin = origin,
                               depth_error = depth_error,
                               calc_distance = calc_distance,
                               calc_movement_pr = calc_movement_pr,
                               mobility = mobility,
                               n = n,
                               cl = cl,
                               use_all_cores = use_all_cores,
                               seed = NULL,
                               verbose = TRUE,
                               dots = ...))
  # Seed
  if(!is.null(seed)) set.seed(seed)
  # Blank list to store path histories
  history <- list()
  # Bathy param
  n_cell  <- raster::ncell(bathy)
  proj    <- raster::crs(bathy)
  mask_na <- is.na(bathy)
  boundaries <- raster::extent(bathy)
  if(is.null(mobility)) bathy_sbt <- bathy
  # Set up distance calculations
  calc_distance <- match.arg(calc_distance)
  if(calc_distance == "lcp") {
    cat_to_cf("... Setting up cost-surface for calc_distance = 'lcp'...")
    if(!all.equal(raster::res(bathy)[1], raster::res(bathy)[2])){
      stop("raster::res(bathy)[1] must equal raster::res(bathy)[2] for this option.")
    }
    costs <- lcp_costs(bathy, verbose = verbose)
    cost  <- costs$dist_total
    graph <- lcp_graph_surface(surface = bathy, cost = cost, verbose = verbose)
  } else {
    if(use_all_cores) {
      warning("use_all_cores = TRUE ignored for calc_distance = 'euclid'.",
              call. = FALSE, immediate. = TRUE)
      use_all_cores <- FALSE
    }
  }
  # Other cluster options
  if(!is.null(cl) & use_all_cores) {
    warning("'cl' and 'use_all_cores' cannot both be specified: setting cl = NULL.",
            call. = FALSE, immediate. = TRUE)
    cl <- NULL
  }
  if(!is.null(cl)) .cl <- cl else .cl <- NULL

  #### Get possible cells the individual could have occupied at each time step based on depth
  cat_to_cf("... Determining possible cells the individual could have occupied based on depth at each time step...")
  cells_by_time <- lapply(archival, function(y){
    # y <- archival[1]
    cells <- flapper::cells_from_val(x = bathy,
                                     y = c(y - depth_error, y + depth_error))
    return(cells)
  })
  # Check there are possible cells at all time steps
  lapply(1:length(cells_by_time), function(t){
    cells_at_time <- cells_by_time[[t]]
    if(length(cells_at_time) < 1) {
      y1 <- archival[t] - depth_error
      y2 <- archival[t] + depth_error
      stop("There are no cells within the required depth range {", y1, ", ", y2, "} at time ", t, ".")
    }
  })

  #### For first time step define set of cells
  cat_to_cf("... Determining the set of possible starting locations (t = 1)...")
  cells_at_time_current     <- cells_by_time[[1]]
  cells_at_time_current     <- data.frame(id_current = cells_at_time_current,
                                          pr_current = 1)
  # Adjust cell probabilities by distance from origin, if applicable, using mobility model
  if(!is.null(origin)){
    # Crop raster to focus on area within mobility for speed
    if(!is.null(mobility)){
      origin_sp <- sp::SpatialPoints(origin, proj4string = proj)
      buf       <- rgeos::gBuffer(origin_sp, width = mobility)
      bathy_sbt <- raster::crop(bathy, buf)
    }
    # Get distances
    if(calc_distance == "euclid"){
      dist_1 <- raster::distanceFromPoints(bathy_sbt, origin)
    } else if(calc_distance == "lcp"){
      get_destination_index <- function(x) !is.na(x) & x >= archival[1] - depth_error & x <= archival[1] + depth_error
      dist_1 <- lcp_from_point(origin = origin,
                               destination = get_destination_index,
                               surface = bathy_sbt,
                               graph = graph,
                               use_all_cores = use_all_cores,
                               verbose = verbose)
    }
    pr_1 <- raster::calc(dist_1, calc_movement_pr)
    pr_1[is.na(pr_1)] <- 0
    pr_1[bathy < archival[1] - depth_error] <- 0
    pr_1[bathy > archival[1] + depth_error] <- 0
    pr_1[mask_na] <- 0
    if(!is.null(mobility)) pr_1 <- raster::extend(pr_1, boundaries, 0)
    cells_at_time_current$pr_current  <- raster::extract(pr_1, cells_at_time_current$id)
  }

  #### Implement algorithm iteratively
  cat_to_cf("... Implementing algorithm iteratively over time steps...")
  for(t in 1:(length(archival))){

    #### For the current time step, select likely cells
    ## Check there are available cells at the current time step; if not
    # ... movement model is insufficient
    # ... algorithm 'unlucky' in that all n paths up to this point are 'dead ends'
    # ... depth_error is too small
    cat_to_cf(paste0("... ... Time = ", t, "..."))
    if(all(cells_at_time_current$pr_current == 0)) {
      message("The probability of all cells at time ", t, " is 0. Either (1) the algorithm has been 'unlucky' and all n = ", n, " paths up to this point are 'dead ends', (2) the mobility model ('mobility') is too limiting and/or (3) the depth error parameter ('depth_error') is too small. The function will now stop, returning outputs up until this point.")
      return(history)
    }
    ## select cells
    cat_to_cf("... ... ... Selecting candidate starting positions for the current time step...")
    cells_at_time_current_sbt <- cells_at_time_current[sample(x = 1:length(cells_at_time_current$id_current),
                                                              size = n,
                                                              prob = cells_at_time_current$pr_current,
                                                              replace = TRUE), ]
    cells_at_time_current_sbt$timestep <- t
    rownames(cells_at_time_current_sbt) <- NULL
    history[[t]] <- cells_at_time_current_sbt

    #### For each particle, identify possible next locations
    if(t <= (length(archival) - 1)) {
      cat_to_cf("... ... ... For each particle, identifying possible positions for the next time step...")
      cells_at_time_next <- cells_by_time[[t + 1]]
      mask_1 <- bathy < archival[t + 1] - depth_error
      mask_2 <- bathy > archival[t + 1] + depth_error
      if(calc_distance == "lcp"){
        get_destination_index <-
          function(x) !is.na(x) & x >= archival[t + 1] - depth_error & x <= archival[t + 1] + depth_error
      }
      cells_from_current_to_next <- pbapply::pblapply(1:n, cl = .cl, function(j){
        # Define location
        # j <- 1
        cell_j <- cells_at_time_current_sbt[j, ]
        cell_j_xy <- raster::xyFromCell(bathy, cell_j$id_current)
        # Subset bathy around location for speed
        if(!is.null(mobility)){
          cell_j_sp <- sp::SpatialPoints(cell_j_xy, proj4string = proj)
          buf       <- rgeos::gBuffer(cell_j_sp, width = mobility)
          bathy_sbt <- raster::crop(bathy, buf)
        }
        if(calc_distance == "euclid"){
          dist_j <- raster::distanceFromPoints(bathy_sbt, cell_j_xy)
        } else if(calc_distance == "lcp"){
          dist_j <- lcp_from_point(origin = cell_j_xy,
                                   destination = get_destination_index,
                                   surface = bathy_sbt,
                                   graph = graph,
                                   use_all_cores = use_all_cores,
                                   verbose = FALSE)
        }
        pr_j <- raster::calc(dist_j, calc_movement_pr)
        pr_j[is.na(pr_j)] <- 0
        pr_j[mask_1]  <- 0
        pr_j[mask_2]  <- 0
        pr_j[mask_na] <- 0
        if(!is.null(mobility)) pr_j <- raster::extend(pr_j, boundaries, 0)
        pr_j <- data.frame(id_current = cell_j$id_current, pr_current = cell_j$pr_current, id_next = 1:n_cell, pr_next = as.vector(pr_j))
        pr_j <- pr_j[pr_j$pr_next > 0, ]
        if(nrow(pr_j) > 0) out <- NULL else out <- pr_j
        return(pr_j)
      })
      cells_from_current_to_next <- plyr::compact(cells_from_current_to_next)
      cells_from_current_to_next <- do.call(rbind, cells_from_current_to_next)
      if(length(cells_from_current_to_next) == 0L){
        message("Unable to sample any cells at time ", t, " for the next location for any of the ", n, " paths. Either (1) the algorithm has been 'unlucky' and all n = ", n, " paths up to this point are 'dead ends', (2) the mobility model ('mobility') is too limiting and/or (3) the depth error parameter ('depth_error') is too small. The function will now stop, returning outputs up until this point.")
        return(history)
      }
      cells_at_time_current <- cells_from_current_to_next
      colnames(cells_at_time_current) <- c("id_previous", "pr_previous", "id_current", "pr_current")
    }
  }
  if(!is.null(cl)) parallel::stopCluster(cl = .cl)

  #### Process paths from cell pairs (id_current, id_next)
  cat_to_cf("... Processing movement paths...")
  path <- list()
  path[[1]] <- history[[1]]
  path[[1]]$id_1 <- path[[1]]$id_current
  path[[1]]$pr_1 <- path[[1]]$pr_current
  path[[1]] <- path[[1]][, c("id_1", "pr_1", "id_current")]
  for(t in 1:(length(archival) - 1)){
    history_for_pair <- dplyr::right_join(path[[t]], history[[t + 1]], by = c("id_current" = "id_previous"))
    history_for_pair[, paste0("id_", t+1)]  <- history_for_pair$id_current.y
    history_for_pair[, paste0("pr_", t+1)]  <- history_for_pair$pr_current
    history_for_pair[, "id_current"] <- history_for_pair$id_current.y
    keep <- c(colnames(path[[t]])[!(colnames(path[[t]]) %in% "id_current")],
              paste0("id_", t+1), paste0("pr_", t+1),
              "id_current")
    history_for_pair <- history_for_pair[, keep]
    path[[t+1]] <- history_for_pair
  }
  # Isolate the last element of the list, which contains all paths
  paths <- path[[length(archival)]]
  paths$id_current <- NULL
  # Reformat paths into a list of dataframes, one for each path, containing the path id, cell id and cell pr
  path_ls <- lapply(1:nrow(paths), function(i){
    d <- paths[i, ]
    dat_for_path <- data.frame(path_id = i,
                               cell_id = as.integer(d[seq(1, ncol(d)-1, by = 2)]),
                               cell_pr = as.numeric(d[seq(2, ncol(d), by = 2)])
    )
    return(dat_for_path)
  })
  # Define a dataframe with each path, including cell coordinates and depths
  path_df <- do.call(rbind, path_ls)
  path_df[, c("cell_x", "cell_y")] <- raster::xyFromCell(bathy, path_df$cell_id)
  path_df[, "cell_z"] <- raster::extract(bathy, path_df[, c("cell_x", "cell_y")])
  path_df[, "timestep"] <- rep(1:length(archival), length(unique(path_df$path_id)))
  path_df <- path_df[, c("path_id", "timestep", "cell_id", "cell_x", "cell_y", "cell_z", "cell_pr")]
  out_dcpf$dcpf <- path_df

  #### Return outputs
  if(!is.null(seed)) set.seed(NULL)
  t_end <- Sys.time()
  total_duration <- difftime(t_end, t_onset, units = "mins")
  cat_to_cf(paste0("... flapper::dcpf() call completed (@ ", t_end, ") after ~", round(total_duration, digits = 2), " minutes."))
  return(out_dcpf)

}


########################################
########################################
#### dcpf_loglik()

#' @title Calculate the log-likelihood of movement paths from the DCPF algorithm
#' @importFrom rlang .data
#' @description This function calculates the total log-likelihood of each movement path reconstructed by the depth-contour particle filtering (DCPF) algorithm.
#' @param paths A dataframe containing movement paths from \code{\link[flapper]{dcpf}}. At a minimum, this should contain a unique identifier for each path (named `path_id') and the probability associated with each cell along each path (`cell_pr').
#' @param ... Additional arguments (none implemented).
#' @details For each path, at each time step the probability associated with the sampled location depends on a user-defined movement model that is driven by the distance between the sampled locations for the individual at the previous and current time steps. This function simply sums the logarithms of these probabilities for each path as a measure of their relative likelihood, given the movement model.
#' @examples
#' # An example with the example DCPF outputs included in flapper
#' dcpf_loglik(dat_dcpf$dcpf)
#' @return The function returns a dataframe with the log likelihood (`loglik') of each path (`path_id'). Rows are ordered by log-likelihood and a `delta' column is provided with the differences in log-likelihood between the most likely path and every other path.
#' @author Edward Lavender
#' @export
#'

dcpf_loglik <- function(paths,...){
  check_names(input = paths, req = c("path_id", "cell_pr"))
  op <- options()
  options(dplyr.summarise.inform = FALSE)
  paths <-
    paths %>%
    dplyr::group_by(.data$path_id) %>%
    dplyr::mutate(log_pr = log(.data$cell_pr)) %>%
    dplyr::summarise(loglik = sum(.data$log_pr)) %>%
    dplyr::mutate(delta = max(.data$loglik) - .data$loglik) %>%
    dplyr::arrange(dplyr::desc(.data$loglik)) %>%
    as.data.frame()
  options(op)
  return(paths)
}


########################################
########################################
#### dcpf_plot_1d()

#' @title Plot one-dimensional depth time series from the DCPF algorithm
#' @description This function plots the observed depth time series and the depth time series associated with each path reconstructed by the depth-contour particle filtering (DCPF) algorithm.
#' @param archival A numeric vector of depth (m) observations, as used by \code{\link[flapper]{dcpf}}.
#' @param paths A dataframe containing reconstructed movement path(s) from \code{\link[flapper]{dcpf}}. At a minimum, this should contain a unique identifier for each path (named `path_id’) and the depth associated with each cell along each path (`cell_z').
#' @param scale A number that vertically scales the depth time series for the observations and the reconstructed path(s). By default, absolute values for depth are assumed and negated for ease of visualisation.
#' @param pretty_axis_args,xlab,ylab,type,... Plot customisation arguments passed to \code{\link[prettyGraphics]{pretty_plot}}.
#' @param add_lines A named list, passed to \code{\link[graphics]{lines}}, to customise the appearance of the depth time series for reconstructed path(s).
#' @param prompt A logical input that defines whether or not plot the observed depth time series with each reconstructed depth time series on a separate plot, sequentially, with a pause between plots (\code{prompt = TRUE}), or with all reconstructed time series on a single plot (\code{prompt = FALSE}).
#' @details Observed and reconstructed depth time series can differ due to measurement error, which is controlled via the \code{depth_error} parameter in the DCPF algorithm (see \code{\link[flapper]{dcpf}}).
#' @examples
#' #### Implement dcpf() algorithm
#' # Here, we extract the necessary information from pre-defined outputs for speed
#' archival <- dat_dcpf$args$archival
#' paths    <- dat_dcpf$dcpf
#'
#' #### Example (1): The default implementation
#' dcpf_plot_1d(archival, paths)
#'
#' #### Example (2): Plot customisation options, e.g.:
#' dcpf_plot_1d(archival, paths, scale = 1, pretty_axis_args = list(side = 1:2))
#' dcpf_plot_1d(archival, paths, type = "l")
#' dcpf_plot_1d(archival, paths, add_lines = list(col = "red", lwd = 0.5))
#'
#' #### Example (3): Plot individual comparisons
#' if(interactive()){
#'   pp <- par(mfrow = c(3, 4))
#'   dcpf_plot_1d(depth, paths, prompt = TRUE)
#'   par(pp)
#' }
#' @return The function returns a plot of the observed and reconstructed depth time series, either for all paths at once (if \code{prompt = FALSE}) or each path separately (if \code{prompt = TRUE}).
#' @seealso \code{\link[flapper]{dcpf}} implements the DCPF algorithm. \code{\link[flapper]{dcpf_plot_1d}}, \code{\link[flapper]{dcpf_plot_2d}} and \code{\link[flapper]{dcpf_plot_3d}} provide plotting routines. \code{\link[flapper]{dcpf_loglik}} calculates the log-likelihood of each path.
#' @author Edward Lavender
#' @export

dcpf_plot_1d <- function(archival,
                         paths,
                         scale = -1,
                         pretty_axis_args = list(side = 3:2),
                         xlab = "Time (index)", ylab = "Depth (m)", type = "b",
                         add_lines = list(col = "royalblue", type = "b"),
                         prompt = FALSE,
                         ...){
  check_names(input = paths, req = c("path_id", "cell_z"))
  if(!prompt) prettyGraphics::pretty_plot(archival*scale,
                                          pretty_axis_args = pretty_axis_args,
                                          xlab = xlab, ylab = ylab,...)
  lapply(split(paths, paths$path_id), function(d){
    if(prompt) prettyGraphics::pretty_plot(archival*scale,
                                           pretty_axis_args = pretty_axis_args,
                                           xlab = xlab, ylab = ylab,...)
    add_lines$x <- d$cell_z*scale
    do.call(graphics::lines, add_lines)
    if(prompt) readline(prompt = "Press [enter] to continue...")
  })
  return(invisible())
}


########################################
########################################
#### dcpf_plot_2d()

#' @title Map two-dimensional paths from the DCPF algorithm
#' @description This function is a simple wrapper for \code{\link[prettyGraphics]{pretty_map}} that maps the paths reconstructed by the DCPF algorithm over a surface.
#' @param paths A dataframe containing reconstructed movement path(s) from \code{\link[flapper]{dcpf}}. At a minimum, this should contain a unique identifier for each path (named `path_id') and the x and y coordinates that define the trajectory of each path (`cell_x' and `cell_y').
#' @param bathy A \code{\link[raster]{raster}} of the bathymetry over which movement was reconstructed.
#' @param add_paths A named list, passed to \code{\link[prettyGraphics]{add_sp_path}}, to customise the appearance of the paths.
#' @param prompt A logical input that defines whether or not plot each path on a separate plot, sequentially, with a pause between plots (\code{prompt = TRUE}), or all paths on a single plot (\code{prompt = FALSE}).
#' @param ... Additional arguments, passed to \code{\link[prettyGraphics]{pretty_map}}, for plot customisation.
#' @return The function maps the trajectories of reconstructed paths across the bathymetry surface, returning a single map if \code{prompt = FALSE} or one map for each path if \code{prompt = TRUE}.
#' @examples
#' #### Implement dcpf() algorithm
#' # Here, we extract the necessary information from pre-defined outputs for speed
#' paths    <- dat_dcpf$dcpf
#' bathy    <- dat_dcpf$args$bathy
#'
#' #### Example (1): The default implementation
#' dcpf_plot_2d(paths, bathy)
#'
#' #### Example (2): Plot customisation options
#' # Customise the appearance of the path(s)
#' dcpf_plot_2d(paths, bathy,
#'              add_paths = list(length = 0.075, col = viridis::viridis(100)))
#' # Pass arguments to prettyGraphics::pretty_map() via ... , e.g.:
#' dcpf_plot_2d(paths, bathy, xlab = "Easting (UTM)", ylab = "Northing (UTM)")
#'
#' #### Example (3): Plot individual paths separately
#' if(interactive()){
#'   pp <- par(mfrow = c(3, 4))
#'   dcpf_plot_2d(paths, bathy, add_paths = list(length = 0.01),
#'                prompt = TRUE, verbose = FALSE)
#'   par(pp)
#' }
#'
#' @seealso \code{\link[flapper]{dcpf}} implements the DCPF algorithm. \code{\link[flapper]{dcpf_plot_1d}}, \code{\link[flapper]{dcpf_plot_2d}} and \code{\link[flapper]{dcpf_plot_3d}} provide plotting routines. For mapping, it can be useful to interpolate shortest (least-cost) paths between sequential locations via \code{\link[flapper]{lcp_interp}}. \code{\link[flapper]{dcpf_loglik}} calculates the log-likelihood of each path.
#' @author Edward Lavender
#' @export
#'

dcpf_plot_2d <- function(paths,
                         bathy,
                         add_paths = list(),
                         prompt = FALSE,...){
  check_names(input = paths, req = c("path_id", "cell_z"))
  if(!prompt) prettyGraphics::pretty_map(add_rasters = list(x = bathy),...)
  lapply(split(paths, paths$path_id), function(d){
    if(prompt) prettyGraphics::pretty_map(add_rasters = list(x = bathy),...)
    add_paths$x <- d$cell_x
    add_paths$y <- d$cell_y
    do.call(prettyGraphics::add_sp_path, add_paths)
    if(prompt) readline(prompt = "Press [enter] to continue...")
  })
  return(invisible())
}


########################################
########################################
#### dcpf_plot_3d()

#' @title Map three-dimensional paths from the DCPF algorithm
#' @description This function is a simple wrapper for \code{\link[prettyGraphics]{pretty_scape_3d}} that maps the paths reconstructed by the DCPF algorithm over a surface in three dimensions.
#' @param paths A dataframe containing reconstructed movement path(s) from \code{\link[flapper]{dcpf}}. At a minimum, this should contain a unique identifier for each path (named `path_id') and the x, y and z coordinates that define the trajectory of each path (`cell_x', `cell_y' and `cell_z').
#' @param bathy A \code{\link[raster]{raster}} of the bathymetry over which movement was reconstructed.
#' @param add_paths A named list, passed to \code{\link[plotly]{add_paths}}, to customise the appearance of the paths.
#' @param shift A number that vertically shifts the paths above the surface (\code{bathy}). The default is \code{shift = 5}, which shifts paths 5 m above the surface. This helps to ensure that paths are visible on interactive, three-dimensional \code{\link[plotly]{plotly}} plots.
#' @param stretch A number that vertically stretches the height of the surface (see \code{\link[prettyGraphics]{pretty_scape_3d}}). The default is \code{-5} which negates the bathymetry and stretches it five-fold.
#' @param aspectmode A character that defines the shape of the plot: \code{"cube"} produces a cube; \code{"data"} produces a plot whether the size of the x, y and z axes is scaled according to the data (see \code{\link[prettyGraphics]{pretty_scape_3d}}).
#' @param prompt A logical input that defines whether or not plot each path on a separate plot, sequentially, with a pause between plots (\code{prompt = TRUE}), or all paths on a single plot (\code{prompt = FALSE}).
#' @param ... Additional arguments, passed to \code{\link[prettyGraphics]{pretty_scape_3d}}, for plot customisation.
#'
#' @return The function maps the trajectories of reconstructed paths across the bathymetry surface in three-dimensions, returning a single map if \code{prompt = FALSE} or one map for each path if \code{prompt = TRUE}. The function also invisibly returns the plot object, if \code{prompt = TRUE}, or a list of plot objects, if \code{prompt = FALSE} (with one element for each path), to facilitate further modification.
#'
#' @examples
#' #### Implement dcpf() algorithm
#' # Here, we extract the necessary information from pre-defined outputs for speed
#' # Note that it may be beneficial to interpolate paths between points
#' # ... e.g., via lcp_interp() prior to plotting, but we will not do that here.
#' paths <- dat_dcpf$dcpf
#' bathy <- dat_dcpf$args$bathy
#'
#' #### Example (1): Visualise paths using the default options
#' dcpf_plot_3d(paths, bathy)
#'
#' #### Example (2): Customise the plot
#' # Customise via add_paths() list
#' dcpf_plot_3d(paths, bathy,
#'              add_paths = list(line = list(color = "black", width = 10),
#'                               marker = list(color = "blue", size = 10)))
#' # Adjust shift, stretch or aspectmode
#' dcpf_plot_3d(paths, bathy, shift = 200, stretch = -10)
#' # Customise via ... e.g., add coastline:
#' coast <- raster::crop(dat_coast, bathy)
#' dcpf_plot_3d(paths, bathy, coastline = coast)
#' # The returned plot objects can also be used for further customisation.
#'
#' #### Example (3): Plot individual paths separately
#' if(interactive()) {
#'   dcpf_plot_3d(paths, bathy, prompt = TRUE)
#' }
#'
#' @details This function requires the \code{\link[plotly]{plotly}} package.
#'
#' @seealso \code{\link[flapper]{dcpf}} implements the DCPF algorithm. \code{\link[flapper]{dcpf_plot_1d}}, \code{\link[flapper]{dcpf_plot_2d}} and \code{\link[flapper]{dcpf_plot_3d}} provide plotting routines. For mapping, it can be useful to interpolate shortest (least-cost) paths between sequential locations via \code{\link[flapper]{lcp_interp}}. \code{\link[flapper]{dcpf_loglik}} calculates the log-likelihood of each path.
#'
#' @author Edward Lavender
#' @export

dcpf_plot_3d <- function(paths,
                         bathy,
                         add_paths = list(line = list(width = 10)),
                         shift = 5,
                         stretch = -5,
                         aspectmode = "data",
                         prompt = FALSE,...){
  # Check for plotly
  if(!requireNamespace("plotly", quietly = TRUE)) stop("This function requires the 'plotly' package. Please install it with `install.packages('plotly')` first.")
  # Define a list of outputs
  p_ls <- list()
  # Plot the surface
  if(!prompt) p <- prettyGraphics::pretty_scape_3d(r = bathy,
                                                   stretch = stretch,
                                                   aspectmode = aspectmode,...)
  # Add paths sequentially to the surface
  paths$cell_z <- paths$cell_z * stretch + shift
  paths_ls <- split(paths, paths$path_id)
  for(i in 1:length(paths_ls)){
    if(prompt) p <- prettyGraphics::pretty_scape_3d(r = bathy,
                                                    stretch = stretch,
                                                    aspectmode = aspectmode,...)
    d <- paths_ls[[i]]
    add_paths$p <- p
    add_paths$x <- d$cell_x
    add_paths$y <- d$cell_y
    add_paths$z <- d$cell_z
    p <- do.call(plotly::add_paths, add_paths)
    if(prompt) {
      print(p)
      readline(prompt = "Press [enter] to continue...")
      p_ls[[i]] <- p
    }
  }
  if(!prompt) {
    print(p)
    p_ls <- p
  }
  # Return invisible plot (list)
  return(invisible(p_ls))
}
