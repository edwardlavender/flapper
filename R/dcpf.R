########################################
########################################
#### dcpf_setup_movement_pr()

#' @title A simple movement model dependent on distance
#' @description This function provides a simple movement model that calculates the probability of movement between two locations according to the distance between them, using an logistic equation with pre-defined parameters.
#'
#' @param distance A numeric vector of distances (m).
#' @param ... Additional arguments (none implemented).
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

dcpf_setup_movement_pr <- function(distance,...) {
  pr <- stats::plogis(10 + distance * -0.05)
  pr[distance > 500] <- 0
  return(pr)
}


########################################
########################################
#### dcpf_setup_cells_by_time()

#' @title Define the \code{cells_by_time} list for the DCPF algorithm
#' @description For a time series of depth observations for a benthic animal (\code{archival}), this function defines a list of integer vectors, one for each time step, that defines the IDs of all of the cells on a bathymetry \code{\link[raster]{raster}} (\code{bathy}) that the individual could have occupied at that time step, given its depth and a measurement error parameter. This procedure is necessary for the depth-contour particle filtering (DCPF) algorithm (see \code{\link[flapper]{dcpf}}).
#'
#' @param archival A dataframe that defines the depth time series for \code{\link[flapper]{dcpf}}. This must contain a column named `depth' with depth observations.
#' @param bathy A \code{\link[raster]{raster}} of the bathymetry for \code{\link[flapper]{dcpf}}.
#' @param calc_depth_error The measurement error function in \code{\link[flapper]{dcpf}}.
#' @param cl (optional) A cluster object created by \code{\link[parallel]{makeCluster}}. If supplied, the connection to the cluster is stopped within the function.
#' @param varlist (optional) A character vector of names of objects to export, to be passed to the \code{varlist} argument of \code{\link[parallel]{clusterExport}}. This may be required if \code{cl} is supplied. Exported objects must be located in the global environment.
#' @param verbose A logical input that defines whether or not to print messages to the console to monitor function progress.
#'
#' @details While this routine can be implemented internally within the DCPF algorithm via \code{\link[flapper]{dcpf}}, it can be beneficial to implement the routine beforehand and pass the resultant list to \code{\link[flapper]{dcpf}}. This is especially the case for large bathymetry surfaces for which computations can take some time: this procedure can be parallelised here, while it is not within \code{\link[flapper]{dcpf}}.
#'
#' @return The function returns a list, with one element for each time step (i.e., observation in \code{archival}) that defines the IDs of all the cells in \code{bathy} that the individual could have occupied at that time step, given the data and the \code{calc_depth_error} function.
#' @examples
#' #### Define a dataframe of depth observations for DCPF algorithm (see ?dcpf)
#' depth <- c(163.06, 159.71, 153.49, 147.04, 139.86, 127.19, 114.75,
#'            99.44,  87.01,  78.16,  70.03,  60.23,  49.96,  35.39,
#'            27.75,  20.13,  12.73,  11.32)
#' depth <- data.frame(depth = depth)
#'
#' #### Define bathymetry surface for DCPF algorithm (see ?dcpf)
#' surface    <- dat_gebco
#' boundaries <- raster::extent(707884.6, 709884.6, 6253404, 6255404)
#' blank      <- raster::raster(boundaries, res = c(5, 5))
#' surface    <- raster::resample(surface, blank)
#'
#' #### Example (1): Implement the function in series
#' cells_by_time <- dcpf_setup_cells_by_time(archival = depth,
#'                                           bathy = surface,
#'                                           calc_depth_error = function(...) c(-30, 30))
#' # The function returns a list of vectors that define all of the cells the individual
#' # ... could possibly have occupied at each time step, based only on depth, the
#' # ... bathymetry and the depth error
#' utils::str(cells_by_time)
#'
#' #### Example (2): Implement the function in parallel
#' # This is only beneficial for large areas/long time series
#' cells_by_time <- dcpf_setup_cells_by_time(archival = depth,
#'                                           bathy = surface,
#'                                           calc_depth_error = function(...) c(-30, 30),
#'                                           cl = parallel::makeCluster(2L),
#'                                           varlist = "surface")
#'
#' #### Example (3): The function returns an error if there are any
#' # ... time points with no possible positions for the individual, given the
#' # ... depth, bathymetry and depth error. This may indicate that the depth error
#' # ... is too small, the area of the bathymetry is too small (i.e., the individual
#' # ... may have moved beyond this area) or other model assumptions (e.g.,
#' # ... benthic habit) are violated.
#' \dontrun{
#'   cells_by_time <- dcpf_setup_cells_by_time(archival = depth,
#'                                             bathy = surface,
#'                                             calc_depth_error = function(...) c(-0.1, 0.1))
#' }
#' @seealso The DCPF algorithm is implemented by \code{\link[flapper]{dcpf}}.
#' @author Edward Lavender
#' @export
#'

dcpf_setup_cells_by_time <- function(archival,
                                     bathy,
                                     calc_depth_error = function(...) c(-2.5, 2.5),
                                     cl = NULL,
                                     varlist = NULL,
                                     verbose = TRUE){

  #### Function set up
  cat_to_console <- function(..., show = verbose){
    if(show) cat(paste(..., "\n"))
  }
  t_onset <- Sys.time()
  cat_to_console(paste0("flapper::dcpf_setup_cells_by_time() called (@ ", t_onset, ")..."))

  #### Checks
  if(!inherits(archival, "data.frame")) stop("'archival' must be a data.frame")
  check_names(input = archival, req = "depth", extract_names = colnames, type = all)
  de_1 <- calc_depth_error(archival$depth[1])
  if(length(de_1) != 2){
    stop("'calc_depth_error' should be a function that returns a numeric vector of length two (i.e., a lower and upper depth adjustment).")
  }
  if(de_1[1] > 0 | de_1[2] < 0){
    stop("'calc_depth_error' should return a negative and a postive adjustment (in that order).")
  }

  #### Get possible cells the individual could have occupied at each time step
  cat_to_console("... Step 1/2: Determining possible cells the individual could have occupied based on depth at each time step...")
  if(!is.null(cl) & !is.null(varlist)) parallel::clusterExport(cl = cl, varlist = varlist)
  cells_by_time <- pbapply::pblapply(archival$depth, cl = cl, function(y){
    depth_lwr <- y + calc_depth_error(y)[1]
    depth_upr <- y + calc_depth_error(y)[2]
    cells <- flapper::cells_from_val(x = bathy, y = c(depth_lwr, depth_upr))
    return(cells)
  })
  if(!is.null(cl)) parallel::stopCluster(cl)

  #### Check there are possible cells at all time steps
  cat_to_console("... Step 2/2: Checking there are possible cells at each time step...")
  pbapply::pblapply(1:length(cells_by_time), function(t){
    cells_at_time <- cells_by_time[[t]]
    if(length(cells_at_time) < 1) {
      y <- archival$depth[t]
      depth_lwr <- y + calc_depth_error(y)[1]
      depth_upr <- y + calc_depth_error(y)[2]
      stop("There are no cells within the required depth range {", depth_lwr, ", ", depth_upr, "} at time ", t, ".", call. = FALSE)
    }
  })
  cat_to_console("... .... Checks passed.")

  #### Return outputs
  t_end <- Sys.time()
  total_duration <- difftime(t_end, t_onset, units = "mins")
  cat_to_console(paste0("... flapper::dcpf_setup_cells_by_time() call completed (@ ", t_end, ") after ~", round(total_duration, digits = 2), " minutes."))
  return(cells_by_time)
}


########################################
########################################
#### dcpf()

#' @title The depth-contour particle filtering (DCPF) algorithm
#' @description This function implements the depth-contour particle filtering (DCPF) algorithm. This is an extension of the DC algorithm (\code{\link[flapper]{dc}}) that implements particle filtering to reconstruct possible movement paths of an individual (i.e., a benthic animal) over a surface (i.e., the seabed). As in the DC algorithm, at each time step the possible locations of an individual are determined from its depth and the bathymetric landscape (plus some measurement error). The extension is the incorporation of a movement model, via a simulation-based particle filtering process, that connects a subset of these locations between time steps into movement paths.
#'
#' To implement this approach, a dataframe of depth observations resulting from movement over a surface (\code{archival}), as well as a \code{\link[raster]{raster}} of the surface over which movement occurred (\code{bathy}) and a function that defines the measurement error at a given depth (\code{calc_depth_error}), must be supplied. A starting location (\code{origin}) can be supplied to constrain the initial set of sampled locations of the individual. At each time step, \code{n} possible locations (`particles') are sampled (with replacement) from the set of possible locations. For each (\code{1:n}) particle, a movement model is used to simulate where the individual could have moved to at the next time step, if it was in any of those locations. In the current framework, the probability of movement into surrounding cells depends on the distance to those cells, which can be represented as using Euclidean or least-cost distances depending on the distance method (\code{calc_distance}), and user-defined movement models (\code{calc_movement_pr_from_origin} and \code{calc_movement_pr}) that link distances to movement probabilities at each time step.
#'
#' At each subsequent time step, this process repeats, with \code{n} possible locations of the individual sampled according to the probability that the individual could have been in that cell, given a previously sampled location, its depth and the bathymetry. The result is a set of locations on the surface at each time step that are consistent with the data and model parameters. Sampled locations can be connected into movement paths via \code{\link[flapper]{dcpf_simplify}}.
#'
#' @param archival A dataframe of depth time series (for a single individual). At a minimum, this should contain a column named `depth' with depth observations. Depth should be recorded using absolute values in the same units as the bathymetry (\code{bathy}, see below). Observations are assumed to have been made at regular time intervals. Other columns can be included for variables that affect movement probabilities (see \code{calc_movement_pr}).
#' @param bathy A \code{\link[raster]{raster}} of the bathymetry in the area within which the individual was located over the study. Bathymetry values should be recorded as absolute values and in the same units (m) as for depths (see \code{archival}). The coordinate reference system should be the Universal Transverse Mercator projection. The resolution of this layer needs to be sufficiently high such that an individual could transition between cells in the duration between \code{archival} observations (see \code{calc_movement_pr}). If the `shortest distances' method is used for distance calculations (i.e., \code{calc_distance = "lcp"}, see below), then the resolution of the surface in x and y directions should also be equal (for surface's with unequal horizontal resolution, \code{\link[raster]{resample}} can be used to equalise resolution: see \code{\link[flapper]{lcp_over_surface}}). For computational efficiency, it is beneficial if \code{bathy} is cropped to the smallest possible area, with any areas in which movement is impossible (e.g., on land for benthic animals) set to NA (see \code{\link[raster]{crop}}, \code{\link[raster]{mask}} and the processing implemented by \code{\link[flapper]{lcp_over_surface}}). It may also be desirable to aggregate high-resolution \code{\link[raster]{raster}}s (see \code{\link[flapper]{process_surface}}).
#' @param cells_by_time (optional) A list, with one element per time step, that defines the IDs of all the cells in \code{bathy} that the individual could have occupied at each time step based on its depth, \code{bathy} and the measurement error (see \code{calc_depth_error}, below). This can be computed via \code{\link[flapper]{dcpf_setup_cells_by_time}}. If un-supplied, it is computed internally.
#' @param origin (optional) A matrix that defines the coordinates (x, y) of the individual's initial location. Coordinates should follow the restrictions for \code{bathy} and lie on a plane (i.e., Universal Transverse Mercator projection).
#' @param calc_depth_error A function that returns the depth error around a given depth. This should accept a single depth value (from \code{archival$depth}) and return two numbers that, when added to that depth, define the range of depths on the bathymetry raster (\code{bathy}) that the individual could plausibly have occupied at any time, given its depth. Since the depth errors are added to the individual's depth, the first number should be negative (i.e., the individual could have been slightly shallower that observed) and the second positive (i.e., the individual could have been slightly deeper than observed). For example, the constant function \code{calc_depth_error = function(...) c(-2.5, 2.5)} implies that the individual could have occupied bathymetric cells whose depth lies within the interval defined by the observed depth + (-2.5) and + (+2.5) m. The appropriate form for \code{calc_depth_error} depends on measurement error for the depth observations in \code{archival} and bathymetry (\code{bathy}) data, as well as the tidal range (m) across the area (over the duration of observations), but this implementation allows the depth error to depend on depth and for the lower and upper error around an observation to differ.
#' @param calc_distance A character that defines the method used to calculate distances between a point (i.e., a sampled location) and the surrounding cells. This drives the probability of movement into those cells via a movement model (see \code{calc_movement_pr} and \code{calc_movement_pr_from_origin}). Currently supported options are Euclidean distances (\code{"euclid"}) or shortest (least-cost) distances (\code{"lcp"}), which represent the shortest distance that an individual would have to move over the surface (\code{bathy}) to traverse between locations (accounting for both planar and vertical distances). Note that this option requires that resolution of \code{bathy} in the x and y directions is equal. At small spatial scales, this option provides more realistic distances in hilly landscapes, but it is more computationally expensive. At larger scales, horizontal distances tend to overwhelm vertical distances, so Euclidean distances may be acceptable. A pragmatic option is to implement the algorithm (possibly for a subset of the \code{archival} time series) using the Euclidean distances method and then interpolate least-cost paths between (a) the subset of sampled locations returned by this approach via \code{\link[flapper]{dcpf_simplify}} or (b) (b) the subset of paths returned by \code{\link[flapper]{dcpf_simplify}} via \code{\link[flapper]{lcp_interp}}. This two-step approach will demonstrate whether sequential positions are plausible (i.e., not too far apart) once the bathymetry is taken into account. If so, the shortest-distances derived using this method can then be used for post-hoc adjustment of movement probabilities. Alternatively, this approach may demonstrate that the algorithm should be re-implemented using the shortest distances method (see \code{\link[flapper]{lcp_interp}}).
#' @param calc_distance_euclid_fast If \code{calc_distance = "euclid"}, \code{calc_distance_euclid_fast} is a logical input that defines whether or not to use the `fast' method for distance and probability calculations. Under the `slow' method (\code{calc_distance_euclid_fast = FALSE}), at each time step the algorithm iterates over each particle, calculates the distance around that particle to neighbouring cells, converts distances to movement probabilities, accounting for the depth, and samples future locations accordingly. Since each particle is considered in turn, each particle has a `history' that is `remembered', which makes path assembly relatively straightforward. In contrast, the faster implementation considers all particles at the same time: a single distance/movement probability surface is calculated around sampled particles, from which future particles are sampled. This approach really excels as the number of particles increases (see \code{n}). The disadvantage is that particles do not `remember' where they have come from and as a result movement path assembly is more involved (see \code{\link[flapper]{dcpf_simplify}}). However, in most cases \code{calc_distance_euclid_fast} is probably desirable, since \code{\link[flapper]{dcpf_simplify}} takes care of movement path assembly.
#' @param calc_distance_graph (optional) If \code{calc_distance = "lcp"}, \code{calc_distance_graph} is a graph object from \code{\link[flapper]{lcp_graph_surface}} that defines the distances between connected cells in \code{bathy}. This is useful in iterative applications in the same environment, especially for large bathymetry rasters (see \code{bathy}) because the calculations of movement costs between adjacent cells and graph construction - both required for the calculation of shortest paths - can be skipped.
#' @param calc_movement_pr The movement model. This must be a function that calculates the probability of movement between two locations in the time between depth observations, given the distance between them (the first argument) and (optionally) any other pertinent information in the archival dataframe for the relevant time step (the second argument). For example, behavioural states could be defined in a column in the \code{archival} dataframe and then included as part of the movement model. Both arguments must be accepted, even if the second is ignored. The default option is a declining logistic curve, designed for flapper skate (\emph{Dipturus intermedius}), representing a high probability of movement between nearby locations and a lower probability of movement between distant locations (see \code{\link[flapper]{dcpf_setup_movement_pr}}). For computational efficiency, it is beneficial if the probability of movement is set to zero beyond a distance considered to be highly unlikely, because such locations are excluded from subsequent calculations (see \code{\link[flapper]{dcpf_setup_movement_pr}} and the \code{mobility} argument, below). Currently, the movement model cannot incorporate the individual's location (for spatially variable fields such as water currents).
#' @param calc_movement_pr_from_origin (optional) If an \code{origin} is supplied, \code{calc_movement_pr_from_origin} can be supplied to define the probability of sampling possible starting locations depending on their distance from the \code{origin}. By default, \code{calc_movement_pr_from_origin = calc_movement_pr} and the same guidance applies to both arguments. Specifying \code{calc_movement_pr_from_origin} specifically may be necessary if the duration between the time at which the \code{origin} was observed and the first \code{archival} observation differs from the regular interval between \code{archival} observations. This allows the effect of the \code{origin} to be weaker or stronger than the effect of locations at later time steps on sampled locations, if necessary.
#' @param mobility (optional) A number that defines the maximum horizontal distance (m) that the individual could travel in the time period between \code{archival} observations. While this is optional, it is usually computationally beneficial to define \code{mobility} because this restricts distance and movement probability calculations at each time step within the smallest appropriate range (rather than across the full surface).
#' @param mobility_from_origin (optional) As above for \code{mobility}, but for the first time step (i.e., the horizontal movement the individual could have moved from the \code{origin}) (see \code{calc_movement_pr_from_origin}).
#' @param n An integer that defines the number of particles (i.e., the number of locations sampled at each time step from the set of possible locations at that time step).
#' @param resample (optional) An integer that defines the minimum number of unique cells that should be sampled, given the movement model (\code{calc_movement_pr} or \code{calc_movement_pr_from_origin}). If supplied, if fewer than \code{resample} unique cells are sampled at a time step, \code{n} particles are re-sampled with replacement with equal probability from all of the cells with non-zero probability. This may facilitate algorithm convergence if there are some cells that are overwhelmingly more probable (given the movement model) are `dead ends': re-sampling all possible cells with equal probability allows the algorithm to explore less likely routes more easily when the number of routes becomes low. \code{resample} must be less than or equal to \code{n}.
#' @param update_history (optional) A list that defines particle histories from an earlier implementation of \code{\link[flapper]{dcpf}} (i.e., the `history' element of a \code{\link[flapper]{.dcpf-class}} object). If provided, the function attempts to continue paths from an earlier time step (see \code{update_history_from_time_step}). This can be useful if the algorithm fails to converge on its initial implementation because the algorithm can be re-started part-way through the time series.
#' @param update_history_from_time_step If \code{update_history} is provided, \code{update_history_from_time_step} is an integer that defines the time step (i.e., element in \code{update_history}) from which to restart the algorithm. If provided, the algorithm continues from this point, by taking the starting positions for \code{n} particles from \code{update_history[[update_history_from_time_step]]$id_current}.
#' @param cl,use_all_cores Parallelisation options. These can be implemented for the approaches that consider particles iteratively (i.e., \code{calc_distance = "euclid"} with \code{calc_distance_euclid_fast = FALSE} or \code{calc_distance = "lcp"}. The algorithm can be parallelised within time steps over (1) paths via \code{cl} (an integer defining the number of child processes (ignored on Windows) or a cluster object created by \code{\link[parallel]{makeCluster}} (see \code{\link[pbapply]{pblapply}})) or (2) within paths for the calculation of shortest distances (if \code{calc_distance = "lcp"}) via a logical input (\code{TRUE}) to \code{use_all_cores}. For \code{calc_distance = "euclid"}, parallelisation is typically only beneficial for relatively large numbers of particles, because the substantial computation overhead associated with parallelisation across paths at each time step is substantial. For \code{calc_distance = "lcp"}, \code{use_all_cores = TRUE} is typically a better option for parallelisation; however, this may only be beneficial if \code{mobility} relatively large. At present, parallelisation is not very effective at present, especially if a cluster is supplied, and may be slower.
#' @param seed (optional) An integer to define the seed for reproducible simulations (see \code{\link[base]{set.seed}}).
#' @param verbose A logical variable that defines whether or not to print messages to the console to relay function progress
#' @param ... Additional arguments (none implemented).
#'
#' @details
#' \subsection{Background}{The DCPF algorithm simulates possible movement paths of a benthic animal over the seabed, given a regular sequence of depth observations (\code{archival}), the bathymetry (\code{bathy}) over which movement occurred and a movement model (\code{calc_movement_pr}) that specifies the probability of movement from a given location to any other, given the distance between them (and any other pre-defined time-dependent parameters in \code{archival}). The function was motivated by small scale applications concerning the reconstruction of possible movement paths of flapper skate (\emph{Dipturus intermedius}) tagged with archival tags, following capture and release in a given location, for short-periods of time post-release.}
#'
#' \subsection{Methods}{At the first time step, the function identifies all of the locations in which the animal could have been located, based on its depth, the bathymetry and some measurement error (determined by \code{calc_depth_error}) that depends on the accuracy of the depth observations, the bathymetry data and the magnitude of the tidal range over the period of observations. From this set of possible locations, \code{n} starting points (`particles') are selected. If an \code{origin} is specified, this selection can be biased towards cells near the origin by the movement model. (Location probability could also be weighted by the proximity between the observed depth and the depths of cells within the range defined by the observed depth and the measurement error (i.e., \code{archival$depth[1] + calc_depth_error(archival$depth[1])[1], archival$depth[1] + calc_depth_error(archival$depth[1])[2]}), but this is not currently implemented.)
#'
#' From each starting position, the Euclidean or shortest distances to cells of the appropriate depth at the next time step are calculated and passed to a movement model that assigns movement probabilities to each cell. Since movement probabilities are likely to be behaviourally dependent, the movement model can also depend on any other relevant information in \code{archival}. However, currently, the model cannot depend on sampled locations; therefore, at least under some conditions, the movement model will need to reflect the maximum distance that an individual could travel within the time period between depth observations, accounting for the possible effects of water currents and any other influences on swimming speed. From the set of cells that are within the required depth range with a movement probability of more than zero, \code{n} particles are sampled, with replacement and according to their probability, and taken as possible starting positions at the next time step. This process repeats until the end of the time series. However, note that in the current implementation of the algorithm, unlike for the start of this process, there is no constraint that forces the individual to return to a specified geographical location at the end of the time series. Indeed, other than an (optional) \code{origin} and the information provided by the depth time series (\code{archival}), geographic restrictions (e.g., from acoustic detections) on the location of the animal over the period of observations cannot be incorporated. Therefore, this algorithm is best-suited to small scale applications (e.g., to examine the movements of individuals tagged with archival for short periods of time immediately post-release). If geographical observations (i.e., detections at acoustic receivers) are available for an individual over its time at liberty, the ACDCPF algorithm (currently unavailable) is required to integrate this information into the construction of movement paths.
#'
#' The result is a set of simulated particles that represent the possible locations of the individual at each time step. This can be assembled into a set of movement paths over a surface that is consistent with the data and the model parameters via \code{\link[flapper]{dcpf_simplify}}. While the number of particles is predetermined by \code{n}, more than \code{n} possible pathways may be defined by all of the combinations of sequential particle movements.}
#'
#' \subsection{Convergence}{Algorithm convergence is not guaranteed. There are four main circumstances in which the algorithm may fail to return any paths that span the start to the end of the depth time series:
#'  \enumerate{
#'    \item \strong{Chance.} All \code{n} paths may be `dead ends'. This possibility can be mitigated by increasing \code{n}.
#'    \item \strong{Movement model.} The movement model may be too limiting. This possibility can be mitigated by ensuring that the movement model realistically represents the probability that an individual can transition between cells given the distance between them. This may be guided by data on the study species or similar species. The movement model may need to account for the effect of water currents, which may increase maximum `swimming' speeds in some directions. (Unfortunately, spatially variable swimming speeds are not currently implemented.) If maximum swimming speeds are uncertain, implementing the algorithm over longer time series (e.g., every \eqn{2^{nd}} observation in \code{archival}), with a suitably relaxed movement model, may facilitate convergence if maximum speeds are unlikely to be maintained for long periods.
#'    \item \strong{Depth error.} The depth error may be too restrictive, given the accuracy of the depth observations, the bathymetry data and the tidal height across an area. This possibility can be mitigated by ensuring that the depth error is appropriate.
#'    \item \strong{Other assumptions.} Other assumptions (e.g., benthic habit) may be violated.
#'  }
#' In these scenarios, the function returns a message that it is about to fail and the results from the start of the algorithm until the current time step, before stopping.
#' }
#'
#' \subsection{Computational considerations}{ This algorithm is computationally intensive. It is advisable that it is run initially with a small time series in a small area and a small number of particles. For larger datasets, there are some tricks that can improve computation time.
#'   \itemize{
#'     \item \strong{Temporal resolution.} Reduce the temporal resolution of the \code{archival} time series so that there are fewer time steps.
#'     \item \strong{Bathymetric resolution.} Reduce the resolution of \code{bathy}, propagating the additional error induced by this process via \code{calc_depth_error} (e.g., see \code{\link[flapper]{process_surface}}).
#'     \item \strong{Mobility limits.} Check whether or not setting mobility limits (\code{mobility}, \code{mobility_from_origin}) improves speed.
#'     \item \strong{Distance calculations}. This step is particularly slow because the distance from each location to many or all surrounding locations are calculated. To speed up this step, implement \code{calc_distance = "euclid"} with \code{calc_distance_euclid_fast = TRUE}. If necessary, interpolate least-cost paths after algorithm completion within \code{\link[flapper]{dcpf_simplify}} or \code{\link[flapper]{lcp_interp}}. In the future, a Markov chain Monte Carlo style approach may be implemented in which distances (and probabilities) for randomly selected `proposal' cells are calculated, with those cells then rejected or retained, rather than calculating distances to many or all surrounding cells, but this is unlikely to be faster in many settings.
#'     \item \strong{Parallelisation.} If the fast Euclidean distances method is not used, test alternative options for parallelisation. For the \code{cl} argument, specifying an integer on non-Windows platforms may be faster than a cluster from \code{\link[parallel]{makeCluster}} (see \code{\link[pbapply]{pblapply}}). Parallelisation may be slower in some circumstances.
#'   }
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
#' depth <- data.frame(depth = depth)
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
#' ## Define depth error function
#' # Because the bathymetry data is very coarse, and the bathymetry is
#' # ... very complex in this region of Scotland, we have to
#' # ... force a high depth error to be high in this example.
#' # The calc_depth_error() function can depend on depth, but in this example
#' # ... we assume the depth error is independent of depth.
#' cde <- function(...) c(-30, 30)
#'
#' ## Check that there are locations on the surface within the requisite depth range
#' # ... at each time step via dcpf_setup_cells_by_time(). This is optional, but
#' # ... will speed up the initial stages of the algorithm because it does not
#' # ... have to compute the cells_by_time list.
#' cells_by_time <- dcpf_setup_cells_by_time(depth, surface, calc_depth_error = cde)
#' # ... Checks passed.
#'
#' ## Define movement model
#' # The default movement model is suitable, with skate moving typically
#' # ... less than 200 m in a two-minute period.
#' # You could use a separate movement model (and mobility restriction) for the origin
#' # ... if necessary, but for brevity we don't implement that here.
#'
#' ## Visualise movement surface, with starting location overlaid
#' prettyGraphics::pretty_map(add_rasters = list(x = surface),
#'                            add_points = list(x = xy),
#'                            verbose = FALSE)
#'
#' #### Example (1): Implement algorithm using default options
#' out_1 <- dcpf(archival = depth,
#'               bathy = surface,
#'               cells_by_time = cells_by_time,
#'               origin  = xy,
#'               calc_depth_error = function(...) c(-30, 30),
#'               calc_distance = "euclid",
#'               calc_movement_pr = dcpf_setup_movement_pr,
#'               n = 10L,
#'               seed = 1)
#' # The function returns a .dcpf-class object
#' class(out_1)
#' utils::str(out_1)
#' # Algorithm duration during testing ~0.03 minutes
#' # ... (vs ~0.21 minutes with calc_distance_euclid_fast = FALSE)
#'
#' #### Example (2): Implement a blanket mobility restriction
#' # This can improve computational efficiency but offers no improvement
#' # ... in this example (e.g., due to relatively small raster, so the cost of
#' # ... focusing in on a specific area matches the speed gains of
#' # ... implementing the distance/movement
#' # ... pr calculations over a smaller area at each time step)
#' out_2 <- dcpf(archival = depth,
#'               bathy = surface,
#'               origin  = xy,
#'               calc_depth_error = function(...) c(-30, 30),
#'               calc_distance = "euclid",
#'               calc_movement_pr = dcpf_setup_movement_pr,
#'               mobility = 250,
#'               n = 10L,
#'               seed = 1)
#' # Algorithm duration during testing ~0.04 minutes
#'
#' \dontrun{
#'
#' #### Example (3): Implement algorithm using shortest distances
#' # Note the need for a surface with equal resolution if this option is implemented.
#' # To speed up the initial stages of the algorithm, you can supply the graph
#' # ... required for least-cost calculations via calc_distance_graph, but for
#' # ... brevity we don't implement that here.
#' out_3 <- dcpf(archival = depth,
#'               bathy = surface,
#'               origin  = xy,
#'               calc_depth_error = function(...) c(-30, 30),
#'               calc_distance = "lcp",
#'               calc_movement_pr = dcpf_setup_movement_pr,
#'               n = 10L,
#'               seed = 1)
#' # This option is slower: algorithm duration during testing ~0.73 minutes
#'
#' #### Example (4): Implement algorithm using shortest distances with mobility restriction
#' out_4 <- dcpf(archival = depth,
#'               bathy = surface,
#'               origin  = xy,
#'               calc_depth_error = function(...) c(-30, 30),
#'               calc_distance = "lcp",
#'               calc_movement_pr = dcpf_setup_movement_pr,
#'               mobility = 250,
#'               n = 10L,
#'               seed = 1)
#' # Algorithm duration during testing ~0.63 minutes
#' # With shortest distances, the mobility restriction makes more difference
#' # ... in improving the computation time, because calculations are more involved.
#'
#' #### Example (5): Parallelisation for Euclidean distances is via cl argument
#' # ... which implements parallelisation across paths. This is only implemented
#' # ... for calc_distance_euclid_fast = FALSE, which is rarely desirable.
#' out_5 <- dcpf(archival = depth,
#'               bathy = surface,
#'               origin  = xy,
#'               calc_depth_error = function(...) c(-30, 30),
#'               calc_distance = "euclid",
#'               calc_distance_euclid_fast = FALSE,
#'               calc_movement_pr = dcpf_setup_movement_pr,
#'               mobility = 250,
#'               n = 10L,
#'               cl = parallel::makeCluster(2L),
#'               seed = 1)
#'
#' #### Example (6): Parallelisation for shortest distances is usually best
#' # ... via use_all_cores = TRUE
#' # However, the benefits of parallelisation depend on the number of least-cost
#' # ... paths calculations that need to be performed, which depends on the size
#' # ... of bathy, and may be minimal (or negative) for small areas.
#' out_6 <- dcpf(archival = depth,
#'               bathy = surface,
#'               origin  = xy,
#'               calc_depth_error = function(...) c(-30, 30),
#'               calc_distance = "lcp",
#'               calc_movement_pr = dcpf_setup_movement_pr,
#'               mobility = 250,
#'               n = 10L,
#'               use_all_cores = TRUE,
#'               seed = 1)
#' # But the speed benefits in this case are minimal.
#' # Algorithm duration during testing ~0.61 minutes.
#'
#' #### Example (7): Extend the movement model via variables in archival
#' # For example, you could assign behavioural 'states', such as resting versus
#' # ... non resting. Here, we use information on the change in depth to restrict
#' # ... movement, based on the idea that large changes in depth correlate
#' # ... with shorter horizontal movements. This is appropriate for calc_distance =
#' # ... 'euclid', which does not account for movement over the bathymetry and
#' # ... may be quicker than implementing least-cost distances.
#' ## Define absolute vertical activity (VA)
#' # For the last observation, we need to assign a VA of 0 to ensure that
#' # ... we do not include NAs in the calculations.
#' depth$va_abs <- abs(depth$depth - dplyr::lead(depth$depth))
#' depth$va_abs[nrow(depth)] <- 0
#' ## Define movement model, depending on distance and a dataframe with other information
#' setup_movement_pr <- function(distance, data){
#'   beta <- -0.05 + -0.005 * data$va_abs
#'   pr <- stats::plogis(10 + distance * beta)
#'   pr[distance > 500] <- 0
#'   return(pr)
#' }
#' ## Examine the movement model with distance and VA
#' pr <- setup_movement_pr(1:1000, data.frame(va_abs = 0))
#' prettyGraphics::pretty_plot(pr, type = "n")
#' for(va_abs in seq(min(depth$va_abs), max(depth$va_abs), length.out = 5)){
#'   lines(1:1000, setup_movement_pr(1:1000, data.frame(va_abs = va_abs)), lwd = 0.25 + va_abs/5)
#' }
#' ## Implement algorithm with adjusted movement model
#' out_7 <- dcpf(archival = depth,
#'               bathy = surface,
#'               origin  = xy,
#'               calc_depth_error = function(...) c(-30, 30),
#'               calc_distance = "euclid",
#'               calc_movement_pr = setup_movement_pr,
#'               mobility = 250,
#'               n = 10L,
#'               seed = 1)
#' # In this case, the algorithm has failed to reach the end of the time series,
#' # ... which provides a useful example to introduce strategies for managing
#' # ... convergence failures (see below). Options include,
#' # ... ... Updating the outputs from an earlier time step
#' # ... ... Increase the number of particles
#' # ... ... Resampling
#' # ... ... Reduce the number of temporal resolution
#' # ... ... Tweak movement model, depth error, mobility etc.
#'
#' #### Example (8): Update particle histories from an earlier time step
#' # ... via update_history and update_history_from_time_step
#' out_8 <- dcpf(archival = depth,
#'               bathy = surface,
#'               origin  = xy,
#'               calc_depth_error = function(...) c(-30, 30),
#'               calc_distance = "euclid",
#'               calc_movement_pr = setup_movement_pr,
#'               mobility = 250,
#'               n = 10L,
#'               seed = 1,
#'               update_history = out_7$history,
#'               update_history_from_time_step = 5)
#'
#' #### Example (9): Implement re-sampling via resample
#' # Here, for demonstration purposes, we implement the algorithm from
#' # ... scratch with the original movement model, this time re-sampling
#' # ... possible positions with equal probability if there are
#' # ... fewer than 10 unique positions.
#' # ... (Normally resample would be less than n.)
#' # In this example, the function re-samples candidate starting positions at t = 9
#' out_9 <- dcpf(archival = depth,
#'               bathy = surface,
#'               origin  = xy,
#'               calc_depth_error = function(...) c(-30, 30),
#'               calc_distance = "euclid",
#'               calc_movement_pr = dcpf_setup_movement_pr,
#'               mobility = 250,
#'               n = 10L,
#'               resample = 10,
#'               seed = 1)
#'
#' #### Example (10): A simulation workflow
#' # This example provides a simulation workflow for comparing simulated and
#' # ... reconstructed paths. Specifically, we simulate movement in an area,
#' # ... so that we know the 'true' path. We then implement the dcpf() algorithm
#' # ... to reconstruct possible paths
#' # ... and compare the true and reconstructed paths.
#'
#' ## (A) Set seed for reproducibility
#' seed <- 2021
#' set.seed(seed)
#'
#' ## (B) Define area for simulation
#' # Here, we use the sample area as above, but reduce the resolution
#' # ... for example speed.
#' dat_gebco_planar <- raster::raster(crs = raster::crs(dat_gebco),
#'                                    ext = raster::extent(dat_gebco),
#'                                    resolution = 25)
#' dat_gebco_planar <- raster::resample(dat_gebco, dat_gebco_planar, method = "bilinear")
#' # Define 'sea' for movement simulation
#' dat_coast <- raster::crop(dat_coast, raster::extent(dat_gebco))
#' dat_sea <- invert_poly(dat_coast)
#' # Visualise area
#' prettyGraphics::pretty_map(dat_gebco_planar,
#'                            add_rasters = list(x = dat_gebco_planar),
#'                            add_polys = list(x = dat_sea))
#'
#' ## (C) Simulate movement path
#' # Define movement parameters
#' sim_steps <- function(...) stats::rgamma(1, shape = 15, scale = 5)
#' prettyGraphics::pretty_hist(stats::rgamma(10000, shape = 15, scale = 5), breaks = 100)
#' # Define animal's origin
#' origin_sim <- sp::spsample(dat_sea, n = 1, type = "random")
#' origin_sim <- sp::coordinates(origin_sim)
#' # Simulate path
#' path_sim <- sim_path_sa(n = 10,
#'                         p_1 = origin_sim,
#'                         area = dat_sea,
#'                         sim_step = sim_steps,
#'                         add_rasters = list(x = dat_gebco_planar),
#'                         seed = seed)
#' # Get resultant depth time series
#' path_sim <- path_sim$xy_mat
#' path_sim <- data.frame(path_id = 1,
#'                        cell_id = raster::cellFromXY(dat_gebco_planar, path_sim),
#'                        cell_x = path_sim[, 1],
#'                        cell_y = path_sim[, 2],
#'                        timestep = 1:nrow(path_sim))
#' path_sim$cell_z <- raster::extract(dat_gebco_planar, path_sim$cell_id)
#' prettyGraphics::pretty_plot(path_sim$cell_z, type = "l")
#' # Simulate 'observed' depth time series given some error
#' # ... For illustration, we will make the error smaller in this example
#' cde <- function(...) c(-2.5, 2.5)
#' depth_obs <- runif(length(path_sim$cell_z),
#'                    path_sim$cell_z + cde(path_sim$cell_z)[1],
#'                    path_sim$cell_z + cde(path_sim$cell_z)[2])
#' depth_obs <- data.frame(depth = depth_obs)
#' # Compare 'observed' and 'true' depth time series
#' dcpf_plot_1d(depth_obs, path_sim,
#'              type = "b", cex = 0.5,
#'              add_lines = list(col = "royalblue", type = "l"))
#'
#' ## Implement dcpf() on 'observed' time series
#' # We will assume that the origin was known.
#' # ... We will mostly use the default options.
#' history_dcpf <- dcpf(archival = depth_obs,
#'                      bathy = dat_gebco_planar,
#'                      origin = origin_sim,
#'                      calc_depth_error = cde,
#'                      calc_distance = "euclid",
#'                      mobility = 200,
#'                      n = 10L
#' )
#'
#' ## Visualise particle histories in relation to simulated paths
#' # Here, each plot shows the particles sampled at a particular time step
#' # ... The green area shows areas of the requisite depth at that time step and the
#' # ... particles show sampled locations at that time step. The simulated path
#' # ... is shown in black.
#' pp <- graphics::par(mfrow = c(3, 4))
#' dcpf_plot_history(history_dcpf,
#'                   add_particles = list(pch = 21),
#'                   add_paths = list(x = path_sim$cell_x, path_sim$cell_y, length = 0.05),
#'                   xlim = range(path_sim$cell_x), ylim = range(path_sim$cell_y),
#'                   crop_spatial = TRUE,
#'                   prompt = FALSE)
#' graphics::par(pp)
#'
#' ## Assemble paths
#' path_dcpf <- dcpf_simplify(history_dcpf)
#'
#' ## Compare 'observed' and reconstructed depth time series
#' dcpf_plot_1d(depth_obs, path_dcpf)
#'
#' ## Show that the distances between sequential positions are within the restrictions
#' # ... of the movement model
#' require(rlang)
#' path_dcpf <-
#'   path_dcpf %>% dplyr::group_by(.data$path_id) %>%
#'   dplyr::mutate(cell_x2 = dplyr::lead(.data$cell_x),
#'                 cell_y2 = dplyr::lead(.data$cell_y),
#'                 dist_1 = sqrt((.data$cell_x2 - .data$cell_x)^2 +
#'                                 (.data$cell_y2 -.data$ cell_y)^2))
#' range(path_dcpf$dist_1, na.rm =TRUE)
#'
#' ## Visualise paths
#' # Zoom around path
#' xlim <- range(c(path_sim$cell_x, path_dcpf$cell_x), na.rm = TRUE)
#' ylim <- range(c(path_sim$cell_y, path_dcpf$cell_y), na.rm = TRUE)
#' boundaries <- raster::extent(xlim, ylim)
#' area <- raster::crop(dat_gebco_planar, boundaries)
#' # Define function to add simulated path for comparison
#' add_paths_sim <-
#'   function() prettyGraphics::add_sp_path(path_sim$cell_x, path_sim$cell_y,
#'                                          lwd = 2, length = 0.01)
#' # Make plots
#' if(interactive()){
#'   dcpf_plot_2d(path_dcpf, area, add_paths = list(length = 0.05),
#'                add_additional = add_paths_sim,
#'                prompt = TRUE)
#' }
#'
#' }
#'
#'
#' @return The function returns a \code{\link[flapper]{.dcpf-class}} object. This is a named list that includes the parameters used to generate function outputs (`args') and the particles sampled at each time step (`history'). The latter can be assembled into a dataframe of movement paths via \code{\link[flapper]{dcpf_simplify}}.
#'
#' @seealso For the movement model, Euclidean distances are obtained from \code{\link[raster]{distanceFromPoints}} or shortest distances are obtained from \code{\link[flapper]{lcp_from_point}} (via \code{\link[flapper]{lcp_costs}} and \code{\link[flapper]{lcp_graph_surface}}, unless \code{calc_distance_graph} is supplied). The default movement model applied to these distances is \code{\link[flapper]{dcpf_setup_movement_pr}}. \code{\link[flapper]{process_behav_rest}} provides a means to assign resting/non-resting behaviour to \code{archival} time steps for behaviourally dependent movement models. Particle histories can be visualised with \code{\link[flapper]{dcpf_plot_history}} and joined into paths via \code{\link[flapper]{dcpf_simplify}}. For processed paths, shortest distances/paths between sequential locations can be interpolated via \code{\link[flapper]{lcp_interp}}, which is a wrapper for the \code{\link[flapper]{lcp_over_surface}} routine. This can be useful for checking whether the faster Euclidean distances method is acceptable and, if so, for post-hoc adjustments of movement probabilities based on shortest distances (see \code{\link[flapper]{lcp_interp}}). Paths can be visualised with \code{\link[flapper]{dcpf_plot_1d}}, \code{\link[flapper]{dcpf_plot_2d}} and \code{\link[flapper]{dcpf_plot_3d}}. The log-likelihood of the paths, given the movement model, can be calculated via \code{\link[flapper]{dcpf_loglik}}. In terms of related algorithms, the DC algorithm is implemented via \code{\link[flapper]{dc}}. This is extended by the ACDC algorithm, through the integration of locational information from acoustic detections, via \code{\link[flapper]{acdc}}. For depth time series accompanied by acoustic detections, movement paths can be reconstructed via the ACDCPF algorithm (currently unavailable).
#'
#' @author Edward Lavender
#' @export

dcpf <- function(archival,
                 bathy,
                 cells_by_time = NULL,
                 origin = NULL,
                 calc_depth_error = function(...) c(-2.5, 2.5),
                 calc_distance = c("euclid", "lcp"),
                 calc_distance_euclid_fast = TRUE,
                 calc_distance_graph = NULL,
                 calc_movement_pr = dcpf_setup_movement_pr,
                 calc_movement_pr_from_origin = calc_movement_pr,
                 mobility = NULL,
                 mobility_from_origin = mobility,
                 n = 10L,
                 resample = NULL,
                 update_history = NULL,
                 update_history_from_time_step = 1,
                 cl = NULL, use_all_cores = FALSE,
                 seed = NULL,
                 verbose = TRUE,...){

  #### Set up function
  # Function onset
  cat_to_console <- function(..., show = verbose){
    if(show) cat(paste(..., "\n"))
  }
  t_onset <- Sys.time()
  cat_to_console(paste0("flapper::dcpf() called (@ ", t_onset, ")..."))
  cat_to_console("... Setting up function...")

  #### Initial checks
  if(!is.null(names(list(...)))){
    warning(paste0("The following argument(s) passed via ... are not supported: ",
                paste(names(list(...)), collapse = ", "), "."),
            call. = FALSE, immediate. = TRUE)
  }
  if(!inherits(archival, "data.frame")) stop("'archival' must be a data.frame")
  check_names(input = archival, req = "depth", extract_names = colnames, type = all)
  if(any(is.na(archival$depth))) stop("'archival$depth' contains NAs.")
  de_1 <- calc_depth_error(archival$depth[1])
  if(length(de_1) != 2){
    stop("'calc_depth_error' should be a function that returns a numeric vector of length two (i.e., a lower and upper depth adjustment).")
  }
  if(de_1[1] > 0 | de_1[2] < 0){
    stop("'calc_depth_error' should return a negative and a postive adjustment (in that order).")
  }
  if(!is.null(resample)){
    if(resample > n) stop("'resample' must be <= 'n'.")
  }

  #### Define necessary objects
  # List for outputs
  out_dcpf <- list(history = NULL,
                   args = list(archival = archival,
                               bathy = bathy,
                               cells_by_time = cells_by_time,
                               origin = origin,
                               calc_depth_error = calc_depth_error,
                               calc_distance = calc_distance,
                               calc_distance_euclid_fast = calc_distance_euclid_fast,
                               calc_distance_graph = calc_distance_graph,
                               calc_movement_pr = calc_movement_pr,
                               calc_movement_pr_from_origin = calc_movement_pr_from_origin,
                               mobility = mobility,
                               mobility_from_origin = mobility_from_origin,
                               n = n,
                               resample = resample,
                               update_history = update_history,
                               update_history_from_time_step = update_history_from_time_step,
                               cl = cl,
                               use_all_cores = use_all_cores,
                               seed = seed,
                               verbose = verbose,
                               dots = list(...)
                   )
  )
  # Seed
  if(!is.null(seed)) set.seed(seed)
  # Blank list to store path histories
  if(is.null(update_history)) history <- list() else history <- update_history
  # Bathy param
  n_cell  <- raster::ncell(bathy)
  proj    <- raster::crs(bathy)
  mask_na <- is.na(bathy)
  boundaries <- raster::extent(bathy)
  if(is.null(mobility_from_origin)) bathy_sbt_1 <- bathy
  if(is.null(mobility)) bathy_sbt <- bathy
  # Origin
  if(!is.null(origin)) if(!inherits(origin, "matrix")) stop("'origin' coordinates must be supplied as a matrix.")
  # Set up distance calculations
  calc_distance <- match.arg(calc_distance)
  if(calc_distance == "lcp") {
    cat_to_console("... Setting up cost-surface for calc_distance = 'lcp'...")
    if(!all.equal(raster::res(bathy)[1], raster::res(bathy)[2])){
      stop("raster::res(bathy)[1] must equal raster::res(bathy)[2] for this option.")
    }
    if(is.null(calc_distance_graph)){
      costs <- lcp_costs(bathy, verbose = verbose)
      cost  <- costs$dist_total
      calc_distance_graph <- lcp_graph_surface(surface = bathy, cost = cost, verbose = verbose)
      out_dcpf$args$calc_distance_graph <- calc_distance_graph
    }
  } else {
    if(use_all_cores){
      warning("use_all_cores = TRUE ignored for calc_distance = 'euclid'.",
              call. = FALSE, immediate. = TRUE)
      use_all_cores <- FALSE
    }
    if(calc_distance_euclid_fast){
      if(!is.null(cl)){
        warning("'cl' argument ignored for calc_distance_euclid_fast = TRUE.",
                call. = FALSE, immediate. = TRUE)
        cl <- NULL
      }
    }
    if(!is.null(calc_distance_graph)){
      warning("'calc_distance_graph' ignored for calc_distance = 'euclid'.",
              call. = FALSE, immediate = TRUE)
      calc_distance_graph <- NULL
    }
  }
  # Other cluster options
  if(!is.null(cl) & use_all_cores){
    warning("'cl' and 'use_all_cores' cannot both be specified: setting cl = NULL.",
            call. = FALSE, immediate. = TRUE)
    cl <- NULL
  }
  .cl <- cl

  #### Get possible cells the individual could have occupied at each time step based on depth
  if(is.null(cells_by_time)){
    cells_by_time <- dcpf_setup_cells_by_time(archival = archival,
                                              bathy = bathy,
                                              calc_depth_error = calc_depth_error,
                                              cl = NULL,
                                              verbose = verbose)
    out_dcpf$args$cells_by_time <- cells_by_time
  }

  #### For first time step define set of cells
  if(is.null(update_history)) {
    cat_to_console("... Determining the set of possible starting locations (t = 1)...")
    arc_1 <- archival[1, , drop = FALSE]
    depth_lwr <- arc_1$depth + calc_depth_error(arc_1$depth)[1]
    depth_upr <- arc_1$depth + calc_depth_error(arc_1$depth)[2]
    cells_at_time_current     <- cells_by_time[[1]]
    cells_at_time_current     <- data.frame(id_current = cells_at_time_current,
                                            pr_current = 1)
    # Adjust cell probabilities by distance from origin, if applicable, using mobility model
    if(!is.null(origin)){
      # Re-define origin on bathy grid for consistency
      # ... This is necessary so that Euclidean and LCP distances are calculated in the same way
      origin_cell_id <- raster::cellFromXY(bathy, origin)
      origin <- raster::xyFromCell(bathy, origin_cell_id)
      # Crop raster to focus on area within mobility_from_origin for speed
      if(!is.null(mobility_from_origin)){
        origin_sp <- sp::SpatialPoints(origin, proj4string = proj)
        buf       <- rgeos::gBuffer(origin_sp, width = mobility_from_origin)
      }
      # Get distances
      if(calc_distance == "euclid"){
        if(!is.null(mobility_from_origin)) bathy_sbt_1 <- raster::crop(bathy, buf)
        dist_1 <- raster::distanceFromPoints(bathy_sbt_1, origin)
        if(!is.null(mobility_from_origin)) dist_1 <- raster::extend(dist_1, boundaries, NA)
      } else if(calc_distance == "lcp"){
        if(!is.null(mobility_from_origin)) bathy_sbt_1 <- raster::mask(bathy, buf)
        get_destination_index <- function(x) !is.na(x) & x >= depth_lwr & x <= depth_upr
        dist_1 <- lcp_from_point(origin = origin,
                                 destination = get_destination_index,
                                 surface = bathy_sbt_1,
                                 graph = calc_distance_graph,
                                 use_all_cores = use_all_cores,
                                 verbose = verbose)
      }
      pr_1 <- raster::calc(dist_1, function(x) calc_movement_pr_from_origin(x, arc_1))
      pr_1[is.na(pr_1)] <- 0
      pr_1[bathy < depth_lwr] <- 0
      pr_1[bathy > depth_upr] <- 0
      pr_1[mask_na] <- 0
      # if(!is.null(mobility_from_origin)) pr_1 <- raster::extend(pr_1, boundaries, 0)
      cells_at_time_current$pr_current  <- raster::extract(pr_1, cells_at_time_current$id)
    }
  } else {
    cat_to_console("... Using update_history...")
    cat_to_console(paste0("... Determining the set of possible starting locations for update (t = ",
                          update_history_from_time_step, ")..."))
    cells_at_time_current <- history[[update_history_from_time_step]]
    print(cells_at_time_current)
  }


  #### Implement algorithm iteratively
  cat_to_console("... Implementing algorithm iteratively over time steps...")
  for(t in update_history_from_time_step:nrow(archival)){

    #### For the current time step, select likely cells
    ## Check there are available cells at the current time step; if not
    # ... movement model is insufficient
    # ... algorithm 'unlucky' in that all n paths up to this point are 'dead ends'
    # ... calc_depth_error is too small
    cat_to_console(paste0("... ... Time = ", t, "..."))
    if(all(cells_at_time_current$pr_current == 0)) {
      message("The probability of all cells at time ", t, " is 0. Either (1) the algorithm has been 'unlucky' and all n = ", n, " particles up to this point have led to 'dead ends', (2) the mobility model ('mobility') is too limiting and/or (3) the depth error parameter from 'calc_depth_error' is too small. The function will now stop, returning outputs up until this point.")
      out_dcpf$history <- history
      class(out_dcpf) <- c(class(out_dcpf), ".dcpf")
      return(out_dcpf)
    }
    ## Select cells
    # Initial selection
    cat_to_console("... ... ... Selecting candidate starting positions for the current time step...")
    cells_at_time_current_sbt <- cells_at_time_current[sample(x = 1:length(cells_at_time_current$id_current),
                                                              size = n,
                                                              prob = cells_at_time_current$pr_current,
                                                              replace = TRUE), ]
    # Re-implement selection, equally across all cells with non 0 Pr, if fewer than a threshold number of cells have been selected
    if(!is.null(resample)) {
      if(length(unique(cells_at_time_current_sbt$id_current)) < resample){
        cat_to_console("... ... ... ... Resampling candidate starting positions for the current time step...")
        cells_at_time_current$pr_current_adj <- cells_at_time_current$pr_current
        cells_at_time_current$pr_current_adj[cells_at_time_current$pr_current_adj != 0] <- 1
        cells_at_time_current_sbt <- cells_at_time_current[sample(x = 1:length(cells_at_time_current$id_current),
                                                                  size = n,
                                                                  prob = cells_at_time_current$pr_current_adj,
                                                                  replace = TRUE), ]
        cells_at_time_current$pr_current_adj <- NULL
      }
    }
    cells_at_time_current_sbt$timestep <- t
    rownames(cells_at_time_current_sbt) <- NULL
    history[[t]] <- cells_at_time_current_sbt

    #### For each particle, identify possible next locations
    if(t <= (nrow(archival) - 1)) {

      ## Set up
      cat_to_console("... ... ... For each particle, identifying possible positions for the next time step...")
      arc_t_next <- archival[t + 1, , drop = FALSE]
      depth_lwr <- arc_t_next$depth + calc_depth_error(arc_t_next$depth)[1]
      depth_upr <- arc_t_next$depth + calc_depth_error(arc_t_next$depth)[2]
      cells_at_time_next <- cells_by_time[[t + 1]]
      mask_1 <- bathy < depth_lwr
      mask_2 <- bathy > depth_upr
      if(calc_distance == "lcp"){
        get_destination_index <-
          function(x) !is.na(x) & x >= depth_lwr & x <= depth_upr
      }

      ## Fast euclidean distances method
      if(calc_distance == "euclid" & calc_distance_euclid_fast){
        # Get cell IDs and coordinates for the 'current' time step
        cell_all_xy <- raster::xyFromCell(bathy, cells_at_time_current_sbt$id_current)
        if(!is.null(mobility)){
          cell_all_sp <- sp::SpatialPoints(cell_all_xy, proj4string = proj)
          buf       <- rgeos::gBuffer(cell_all_sp, width = mobility)
          bathy_sbt <- raster::crop(bathy, buf)
        }
        # Get distances around cells
        dist_all <- raster::distanceFromPoints(bathy_sbt, cell_all_xy)
        if(!is.null(mobility)) dist_all <- raster::extend(dist_all, boundaries, NA)
        # Get probabilities around cells
        pr_all <- raster::calc(dist_all, function(x) calc_movement_pr(x, arc_t_next))
        if(!is.null(mobility)) pr_all <- raster::mask(pr_all, buf, updatevalue = 0)
        pr_all[is.na(pr_all)] <- 0
        pr_all[mask_1]  <- 0
        pr_all[mask_2]  <- 0
        pr_all[mask_na] <- 0
        # Define a dataframe of probabilities
        pr_all <- data.frame(id_previous = NA,
                             pr_previous = NA,
                             id_current = 1:n_cell,
                             pr_current = as.vector(pr_all))
        cells_from_current_to_next <- pr_all[pr_all$pr_current > 0, ]

      ## Other distance methods (point-by-point)
      } else {
        cells_from_current_to_next <- pbapply::pblapply(1:n, cl = .cl, function(j){
          # Define location
          # j <- 1
          cell_j <- cells_at_time_current_sbt[j, ]
          cell_j_xy <- raster::xyFromCell(bathy, cell_j$id_current)
          # Subset bathy around location for speed
          if(!is.null(mobility)){
            cell_j_sp <- sp::SpatialPoints(cell_j_xy, proj4string = proj)
            buf       <- rgeos::gBuffer(cell_j_sp, width = mobility)
            # bathy_sbt <- raster::crop(bathy, buf)
          }
          if(calc_distance == "euclid"){
            if(!is.null(mobility)) bathy_sbt <- raster::crop(bathy, buf)
            dist_j <- raster::distanceFromPoints(bathy_sbt, cell_j_xy)
            if(!is.null(mobility)) dist_j <- raster::extend(dist_j, boundaries, NA)
          } else if(calc_distance == "lcp"){
            if(!is.null(mobility)) bathy_sbt <- raster::mask(bathy, buf)
            dist_j <- lcp_from_point(origin = cell_j_xy,
                                     destination = get_destination_index,
                                     surface = bathy_sbt,
                                     graph = calc_distance_graph,
                                     use_all_cores = use_all_cores,
                                     verbose = FALSE)
          }
          pr_j <- raster::calc(dist_j, function(x) calc_movement_pr(x, arc_t_next))
          pr_j[is.na(pr_j)] <- 0
          pr_j[mask_1]  <- 0
          pr_j[mask_2]  <- 0
          pr_j[mask_na] <- 0
          # if(!is.null(mobility)) pr_j <- raster::extend(pr_j, boundaries, 0)
          pr_j <- data.frame(id_current = cell_j$id_current, pr_current = cell_j$pr_current, id_next = 1:n_cell, pr_next = as.vector(pr_j))
          pr_j <- pr_j[pr_j$pr_next > 0, ]
          if(nrow(pr_j) > 0) out <- NULL else out <- pr_j
          return(pr_j)
        })
        cells_from_current_to_next <- plyr::compact(cells_from_current_to_next)
        cells_from_current_to_next <- do.call(rbind, cells_from_current_to_next)
      }

      ## Processing
      if(length(cells_from_current_to_next) == 0L){
        message("Unable to sample any cells at time ", t, " for the next location for any of the ", n, " particles. Either (1) the algorithm has been 'unlucky' and all n = ", n, " particles up to this point have led to 'dead ends', (2) the mobility model ('mobility') is too limiting and/or (3) the depth error parameter from 'calc_depth_error' is too small. The function will now stop, returning outputs up until this point.")
        out_dcpf$history <- history
        class(out_dcpf) <- c(class(out_dcpf), ".dcpf")
        return(out_dcpf)
      }
      cells_at_time_current <- cells_from_current_to_next
      colnames(cells_at_time_current) <- c("id_previous", "pr_previous", "id_current", "pr_current")
    }
  }
  if(!is.null(.cl) & inherits(.cl, "cluster")) parallel::stopCluster(cl = .cl)

  #### Define outputs
  out_dcpf$history <- history
  class(out_dcpf) <- c(class(out_dcpf), ".dcpf")

  #### Return outputs
  if(!is.null(seed)) set.seed(NULL)
  t_end <- Sys.time()
  total_duration <- difftime(t_end, t_onset, units = "mins")
  cat_to_console(paste0("... flapper::dcpf() call completed (@ ", t_end, ") after ~", round(total_duration, digits = 2), " minutes."))
  return(out_dcpf)

}


########################################
########################################
#### dcpf_simplify()

#' @title Convert particle histories from \code{\link[flapper]{dcpf}} into movement paths
#' @description This function is designed to simplify the \code{\link[flapper]{.dcpf-class}} object from \code{\link[flapper]{dcpf}} that defines sampled particle histories into a set of movement paths. The function identifies pairs of cells between which movement may have occurred at each time step (if necessary), (re)calculates distances and probabilities between connected cell pairs and then links pairwise movements between cells into a set of possible movement paths.
#' @param record A \code{\link[flapper]{.dcpf-class}} object from \code{\link[flapper]{dcpf}}.
#' @param calc_distance A character that defines the method used to calculate distances between sequential combinations of particles (see \code{\link[flapper]{dcpf}}). Currently supported options are Euclidean distances (\code{"euclid"}) or shortest (least-cost) distances ("lcp"). Note that \code{calc_distance} does not need to be the same method as used for \code{\link[flapper]{dcpf}}: it is often computationally beneficial to implement \code{\link[flapper]{dcpf}} using Euclidean distances and then, for the subset of sampled particles, implement in \code{\link[flapper]{dcpf_simplify}} with \code{calc_distance = "lcp"} to re-compute distances using the shortest-distances algorithm, along with the adjusted probabilities. However, for large paths, the quickest option is to implement both functions using \code{calc_distance = "euclid"} and then interpolate shortest paths only for the set of returned paths (see \code{\link[flapper]{lcp_interp}}). If \code{calc_distance = NULL}, the method saved in \code{record} is used.
#' @param calc_distance_graph (optional) If \code{calc_distance = "lcp"}, \code{calc_distance_graph} is a graph object that defines the distances between connected cells in the bathymetry layer (stored in \code{record$args$bathy}). If un-supplied, this is taken from \code{record$args$calc_distance_graph}, if available, or computed via \code{\link[flapper]{lcp_graph_surface}}.
#' @param max_n_copies (optional) An integer that specifies the maximum number of copies of a sampled cell that are retained at each time stamp. Each copy represents a different route to that cell. By default, all copies (i.e. routes to that cell are retained) via \code{max_n_copies = NULL}. However, in cases where there are a large number of paths through a landscape, the function can run into vector memory limitations during path assembly, so \code{max_n_copies} may need to be set. In this case, at each time step, if there are more than \code{max_n_copies} paths to a given cell, then a subset of these (\code{max_n_copies}) are sampled, according to the \code{sample_method} argument.
#' @param sample_method (optional) If \code{max_n_copies} is supplied, \code{sample_method} is a character that defines the sampling method. Currently supported options are: \code{"random"}, which implements random sampling; \code{"weighted"}, which implements weighted sampling, with random samples taken according to their probability at the current time step; and \code{"max"}, which selects for the top \code{max_n_copies} most likely copies of a given cell according to the probability associated with movement into that cell from the previous location.
#' @param add_origin A logical input that defines whether or not to include the origin in the returned dataframe.
#' @param verbose A logical input that defines whether or not to print messages to the console to monitor function progress.
#' @param ... Additional arguments (none implemented).
#'
#' @details The implementation of this function depends on how \code{\link[flapper]{dcpf}} has been implemented. Under the default options in \code{\link[flapper]{dcpf}}, the fast Euclidean distances method is used to sample sequential particle positions, in which case the history of each particle through the landscape is not retained and has to be assembled afterwards. In this case, \code{\link[flapper]{dcpf_simplify}} calculates the distances between all combinations of cells at each time step, using either a Euclidean distances or shortest distances algorithm according to the input to \code{calc_distance}. Distances are converted to probabilities using the movement models retained in \code{record} from the call to \code{\link[flapper]{dcpf}} to identify possible movement paths between cells at each time step. Pairwise cell movements are then assembled into complete movement paths. If the fast Euclidean distances method has not been used, then pairwise cell movements are retained by  \code{\link[flapper]{dcpf}}. In this case, the function simply recalculates distances between sequential cell pairs and the associated movement probabilities, which are used to assemble a set of movement paths.
#'
#' @return The function returns a \code{\link[flapper]{dcpf-class}} object, which is a dataframe that defines the movement paths.
#'
#' @examples
#' #### Example particle histories
#' # In these examples, we will use the example particle histories included in flapper
#' summary(dat_dcpf_histories)
#'
#' #### Example (1): The default implementation
#' paths_1 <- dcpf_simplify(dat_dcpf_histories)
#'
#' ## Demonstration that the distance and probabilities calculations are correct
#' # This method works for Euclidean distances and when calc_movement_pr applies to
#' # ... all time stamps.
#' require(magrittr)
#' require(rlang)
#' paths_1 <-
#'   paths_1 %>% dplyr::group_by(.data$path_id) %>%
#'   dplyr::mutate(cell_xp = dplyr::lag(.data$cell_x),
#'                 cell_yp = dplyr::lag(.data$cell_y),
#'                 cell_dist_chk = sqrt((.data$cell_xp - .data$cell_x)^2 +
#'                                        (.data$cell_yp - .data$cell_y)^2),
#'                 cell_pr_chk = dat_dcpf_histories$args$calc_movement_pr(.data$cell_dist_chk),
#'                 dist_equal = .data$cell_dist_chk == .data$cell_dist_chk,
#'                 pr_equal = .data$cell_pr == .data$cell_pr_chk) %>%
#'   data.frame()
#' utils::head(paths_1)
#'
#' ## Demonstration that the depths of sampled cells are correct
#' paths_1$cell_z_chk <- raster::extract(dat_dcpf_histories$args$bathy,
#'                                       paths_1$cell_id)
#' all.equal(paths_1$cell_z, paths_1$cell_z_chk)
#'
#' ## Compare depth time series
#' # There is a relatively large degree of mismatch here, which reflects
#' # ... the low resolution bathymetry data used for the algorithm.
#' dcpf_plot_1d(dat_dcpf_histories$args$archival, paths_1)
#'
#' ## Examine paths
#' # Log likelihood
#' dcpf_loglik(paths_1)
#' # 2-d visualisation
#' dcpf_plot_2d(paths_1, dat_dcpf_histories$args$bathy,
#'              add_paths = list(length = 0.05))
#' # 3-d visualisation
#' dcpf_plot_3d(paths_1, dat_dcpf_histories$args$bathy)
#'
#' #### Example (2): Re-calculate distances using another method
#' # Use shortest distances:
#' paths_2a <- dcpf_simplify(dat_dcpf_histories, calc_distance = "lcp")
#' # Speed up shortest distance calculations by supplying the graph object:
#' costs <- lcp_costs(dat_dcpf_histories$args$bathy)
#' graph <- lcp_graph_surface(dat_dcpf_histories$args$bathy, costs$dist_total)
#' paths_2b <- dcpf_simplify(dat_dcpf_histories,
#'                           calc_distance = "lcp",
#'                           calc_distance_graph = graph)
#'
#' ## Demonstrate the LCP calculations are correct
#' paths_2d_lcps <- lcp_interp(paths_2b,
#'                             dat_dcpf_histories$args$bathy,
#'                             calc_distance = TRUE)
#' head(cbind(paths_2b$dist, paths_2d_lcps$dist_lcp$dist))
#'
#' #### Example (3): Restrict the number of routes to each cell at each time step
#' # Implement approach for different numbers of copies
#' # Since we only have sampled a small number of particles for this simulation
#' # ... this does not make any difference here, but it can dramatically reduce
#' # ... the time taken to assemble paths and prevent vector memory issues.
#' paths_3a <- dcpf_simplify(dat_dcpf_histories, max_n_copies = 1)
#' paths_3b <- dcpf_simplify(dat_dcpf_histories, max_n_copies = 5)
#' paths_3c <- dcpf_simplify(dat_dcpf_histories, max_n_copies = 7)
#' # Compare the number of paths retained
#' unique(paths_3a$path_id)
#' unique(paths_3b$path_id)
#' unique(paths_3c$path_id)
#'
#' #### Example (4): Change the sampling method used to retain paths
#' # Again, this doesn't make a difference here, but it can when there are
#' # ... more particles.
#' paths_4a <- dcpf_simplify(dat_dcpf_histories,
#'                           max_n_copies = 5,
#'                           sample_method = "random")
#' paths_4b <- dcpf_simplify(dat_dcpf_histories,
#'                           max_n_copies = 5,
#'                           sample_method = "weighted")
#' paths_4c <- dcpf_simplify(dat_dcpf_histories,
#'                           max_n_copies = 5,
#'                           sample_method = "max")
#' # Compare retained paths
#' dcpf_loglik(paths_3a)
#' dcpf_loglik(paths_3b)
#' dcpf_loglik(paths_3c)
#'
#' #### Example (5): Retain/drop the origin, if specified
#' # For the example particle histories, an origin was specified
#' dat_dcpf_histories$args$origin
#' # This is included as 'timestep = 0' in the returned dataframe
#' # ... with the coordinates re-defined on bathy:
#' paths_5a <- dcpf_simplify(dat_dcpf_histories)
#' paths_5a[1, c("cell_x", "cell_y")]
#' raster::xyFromCell(dat_dcpf_histories$args$bathy,
#'                    raster::cellFromXY(dat_dcpf_histories$args$bathy,
#'                                       dat_dcpf_histories$args$origin))
#' head(paths_5a)
#' # If specified, the origin is dropped with add_origin = FALS
#' paths_5b <- dcpf_simplify(dat_dcpf_histories, add_origin = FALSE)
#' head(paths_5b)
#'
#' @author Edward Lavender
#' @export

dcpf_simplify <- function(record,
                          calc_distance = NULL,
                          calc_distance_graph = NULL,
                          max_n_copies = NULL,
                          sample_method = c("random", "weighted", "max"),
                          add_origin = TRUE,
                          verbose = TRUE,...){


  ########################################
  #### Function set up

  #### Set up
  cat_to_console <- function(..., show = verbose){
    if(show) cat(paste(..., "\n"))
  }
  t_onset <- Sys.time()
  cat_to_console(paste0("flapper::dcpf_simplify() called (@ ", t_onset, ")..."))
  if(!inherits(record, ".dcpf")) stop("'record' must be a '.dcpf-class' object.")
  history  <- record$history
  archival <- record$args$archival
  bathy    <- record$args$bathy
  origin   <- record$args$origin
  if(!is.null(origin)) {
    origin_cell_id <- raster::cellFromXY(bathy, origin)
    origin         <- raster::xyFromCell(bathy, origin_cell_id)
  }
  if(is.null(calc_distance)) calc_distance <- record$args$calc_distance
  calc_distance                <- match.arg(calc_distance, c("euclid", "lcp"))
  calc_movement_pr_from_origin <- record$args$calc_movement_pr_from_origin
  calc_movement_pr             <- record$args$calc_movement_pr
  mobility                     <- record$args$mobility
  mobility_from_origin         <- record$args$mobility_from_origin
  sample_method                <- match.arg(sample_method)


  ########################################
  #### Identify sequential, pairwise connections between cells

  ## Implement setup for LCP distance calculations, if necessary.
  cat_to_console(paste("... Getting pairwise cell movements based on calc_distance = ", calc_distance, "..."))
  if(calc_distance == "lcp"){
    cat_to_console("... Setting up LCP calculations...")
    if(is.null(calc_distance_graph)) calc_distance_graph <- record$args$calc_distance_graph
    if(is.null(calc_distance_graph)){
      cat_to_console("... ... Setting up cost-surface for calc_distance = 'lcp'...")
      costs <- lcp_costs(bathy, verbose = verbose)
      cost  <- costs$dist_total
      calc_distance_graph <- lcp_graph_surface(surface = bathy, cost = cost, verbose = verbose)
    }
  }

  #### Calculate distances and probabilities between cell pairs
  # ... (1) Method for DCPF outputs derived via calc_distance_euclid_fast
  cat_to_console("... ... Stepping through time steps to join coordinate pairs...")
  if(record$args$calc_distance == "euclid" & record$args$calc_distance_euclid_fast){

    #### Add element to history for movement from time 0 to time 1
    if(!is.null(origin)){
      origin_dat <- data.frame(id_current = origin_cell_id,
                               pr_current = 1,
                               timestep = 0)
      history <- append(list(origin_dat), history)
    } else {
      history <- append(history[1], history)
      history[[1]]$timestep <- 0
    }

    #### Re-define history with cell -pairs
    history <- pbapply::pblapply(2:length(history), function(t){

      ## Get retained cells from the previous and current step and coordinates
      z1 <- history[[t - 1]]
      z1_xy <- raster::xyFromCell(bathy, z1$id_current)
      z2 <- history[[t]]
      z2_xy <- raster::xyFromCell(bathy, z2$id_current)

      ## Calculate Euclidean distances between all combinations of cells to identify possible movement paths between cell pairs
      # For both the origin (if applicable) and subsequent locations, we will calculate distances from cells IDs on bathy
      # ... For the origin, this is may result in slightly less accurate distances,
      # ... but it ensures consistency across cells and methods (euclidean and lcps)
      # This returns a matrix
      # ... First argument forms rows of output
      # ... Second argument forms columns of output
      if(calc_distance == "euclid"){
        dist_btw_cells <- sp::spDists(z1_xy, z2_xy, longlat = FALSE, segments = FALSE, diagonal = FALSE)
      } else if(calc_distance == "lcp"){
        dist_btw_cells <- cppRouting::get_distance_matrix(Graph = calc_distance_graph,
                                                          from = z1$id_current,
                                                          to = z2$id_current,...)
      }
      # Adjust distances according to mobility_from_origin or mobility parameters
      if((t - 1) == 1){
        if(!is.null(origin)){
          if(is.null(mobility_from_origin)) dist_btw_cells[dist_btw_cells > mobility_from_origin] <- NA
        } else {
          dist_btw_cells <- matrix(c(rep(0, nrow(z2))), nrow = 1)
        }
      } else {
        if(!is.null(mobility)) dist_btw_cells[dist_btw_cells > mobility] <- NA
      }

      ## Calculate probabilities of movement FROM previous cell INTO current cell
      # ... I.e., probability associated with CURRENT cells
      # ... For the first time step, probabilities depend on the origin or are simply assigned 1
      # ... For later time steps, probabilities depend on distances (and any other parameters in archival)
      if((t - 1) == 1){
        if(!is.null(origin)){
          pr_btw_cells <- calc_movement_pr_from_origin(dist_btw_cells, archival[t-1, ])
        } else {
          pr_btw_cells <- matrix(c(rep(1, nrow(z2))), nrow = 1)
        }
      } else {
        pr_btw_cells <- calc_movement_pr(dist_btw_cells, archival[t-1, ])
      }

      ## Define dataframe from matrix that includes all possible movements between sampled locations at t-1 and t
      tmp <- data.frame(id_previous = NA,
                        pr_previous = NA,
                        row = as.vector(row(pr_btw_cells)),
                        col = as.vector(col(pr_btw_cells)),
                        id_current = NA,
                        pr_current = as.vector(pr_btw_cells),
                        dist_current = as.vector(dist_btw_cells),
                        timestep = z1$timestep[1] + 1
                        )
      tmp <- tmp[!is.na(tmp$pr_current), ]
      tmp$id_previous <- z1$id_current[tmp$row]
      if(is.null(origin) & (t - 1) == 1) tmp$id_previous <- NA
      tmp$id_current  <- z2$id_current[tmp$col]
      tmp$row <- NULL
      tmp$col <- NULL
      return(tmp)
    })

  # ... (2) Method for DCPF outputs not derived via calc_distance_euclid_fast
  } else {

    #### Add origin to history, if necessary
    if(!is.null(origin)){
      history[[1]]$timestep      <- 0
      history[[1]]$id_previous   <- origin_cell_id
      history[[1]]$pr_previous   <- 1
      # history[[1]]$id_previous_x <- origin[1]
      # history[[1]]$id_previous_y <- origin[2]
      # history[[1]][, c("id_current_x", "id_current_y")]  <- raster::xyFromCell(bathy, history[[1]]$id_current)
    }

    #### Update history with distances and probabilities of movement between connected cells
    history <- lapply(1:length(history), function(t){

      #### Get history for time t
      d <- history[[t]]
      d$timestep <- t

      #### Calculate distances between connected cells
      if(!rlang::has_name(d, "id_previous")){
        d$dist_current <- NA
        d$pr_current   <- 1
      } else{
        d[, c("id_previous_x", "id_previous_y")]  <- raster::xyFromCell(bathy, d$id_previous)
        d[, c("id_current_x", "id_current_y")]    <- raster::xyFromCell(bathy, d$id_current)
        if(calc_distance == "euclid"){
          d$dist_current <-
            raster::pointDistance(d[, c("id_previous_x", "id_previous_y")],
                                  d[, c("id_current_x", "id_current_y")], lonlat = FALSE)
        } else if(calc_distance == "lcp"){
          d$dist_current <- cppRouting::get_distance_pair(Graph = calc_distance_graph,
                                                          from = d$id_previous,
                                                          to = d$id_current,...)
        }
      }

      #### Re-calculate probabilities from distances between connected cells
      if(t == 1){
        if(is.null(origin)){
          d$pr_current <- 1
        } else {
          d$pr_current <- calc_movement_pr_from_origin(d$dist_current, archival[1, ])
        }
      } else {
        d$pr_current <- calc_movement_pr(d$dist_current, archival[t, ])
      }
      return(d)
    })
  }


  ########################################
  #### Assemble and format paths

  #### Assemble paths
  cat_to_console("... Assembling paths...")
  path <- list()
  path[[1]] <- history[[1]]
  path[[1]]$id_1 <- path[[1]]$id_current
  path[[1]]$pr_1 <- path[[1]]$pr_current
  path[[1]]$dist_1 <- path[[1]]$dist_current
  path[[1]] <- path[[1]][, c("id_1", "pr_1", "dist_1", "id_current")]
  for(t in 1:(length(history) - 1)){
    history_for_pair <- dplyr::inner_join(path[[t]], history[[t + 1]], by = c("id_current" = "id_previous"))
    history_for_pair <- dplyr::distinct(history_for_pair)
    if(!is.null(max_n_copies)){
      history_for_pair <-
        history_for_pair %>%
        dplyr::group_by(.data$id_current.y) %>%
        dplyr::arrange(dplyr::desc(.data$pr_current))
      if(sample_method == "random"){
        history_for_pair <- history_for_pair %>% dplyr::slice_sample(n = max_n_copies)
      } else if(sample_method == "weighted"){
        history_for_pair <- history_for_pair %>% dplyr::slice_sample(n = max_n_copies, weight_by = .data$pr_current)
      } else if(sample_method == "max"){
        history_for_pair <- history_for_pair %>% dplyr::slice(1:max_n_copies)
      }
    }
    history_for_pair[, paste0("id_", t+1)]  <- history_for_pair$id_current.y
    history_for_pair[, paste0("pr_", t+1)]  <- history_for_pair$pr_current
    history_for_pair[, paste0("dist_", t+1)]  <- history_for_pair$dist_current
    history_for_pair[, "id_current"] <- history_for_pair$id_current.y
    keep <- c(colnames(path[[t]])[!(colnames(path[[t]]) %in% c("id_current", "dist"))],
              paste0("id_", t+1), paste0("pr_", t+1), paste0("dist_", t+1),
              "id_current")
    history_for_pair <- history_for_pair[, keep]
    path[[t+1]] <- history_for_pair
  }
  # Isolate the last element of the list, which contains all paths
  paths <- path[[length(history)]]
  paths$id_current <- NULL

  #### Reformat paths
  # Reformat paths into a list of dataframes, one for each path, containing the path id, cell id and cell pr
  cat_to_console("... Formatting paths...")
  path_ls <- pbapply::pblapply(1:nrow(paths), function(i){
    d <- paths[i, ]
    # colnames(d)[seq(1, ncol(d), by = 3)]
    # colnames(d)[seq(2, ncol(d), by = 3)]
    # colnames(d)[seq(3, ncol(d), by = 3)]
    dat_for_path <- data.frame(path_id = i,
                               timestep = 1:(ncol(d)/3),
                               cell_id = as.integer(d[seq(1, ncol(d), by = 3)]),
                               cell_pr = as.numeric(d[seq(2, ncol(d), by = 3)]),
                               dist = as.numeric(d[seq(3, ncol(d), by = 3)])
                               )
    if(!is.null(origin) & add_origin){
      d_origin <- data.frame(path_id = i,
                             timestep = 0,
                             cell_id = origin_cell_id,
                             cell_pr = 1,
                             dist = NA)
      dat_for_path <- rbind(d_origin, dat_for_path)
    }
    return(dat_for_path)
  })

  #### Define a dataframe with each path, including cell coordinates and depths
  cat_to_console("... Adding cell coordinates and depths...")
  path_df <- do.call(rbind, path_ls)
  path_df[, c("cell_x", "cell_y")] <- raster::xyFromCell(bathy, path_df$cell_id)
  if(!is.null(origin) & add_origin){
    path_df[path_df$timestep == 0, "cell_x"] <- origin[1]
    path_df[path_df$timestep == 0, "cell_y"] <- origin[2]
  }
  path_df[, "cell_z"] <- raster::extract(bathy, path_df[, c("cell_x", "cell_y")])
  path_df <- path_df[, c("path_id", "timestep", "cell_id", "cell_x", "cell_y", "cell_z", "cell_pr", "dist")]
  class(path_df) <- c(class(path_df), "dcpf")

  #### Return outputs
  t_end <- Sys.time()
  total_duration <- difftime(t_end, t_onset, units = "mins")
  cat_to_console(paste0("... flapper::dcpf_simplify() call completed (@ ", t_end, ") after ~", round(total_duration, digits = 2), " minutes."))
  return(path_df)
}


########################################
########################################
#### dcpf_loglik()

#' @title Calculate the log-likelihood of movement paths from the DCPF algorithm
#' @importFrom rlang .data
#' @description This function calculates the total log-likelihood of each movement path reconstructed by the depth-contour particle filtering (DCPF) algorithm.
#' @param paths A dataframe containing movement paths from \code{\link[flapper]{dcpf}} plus \code{\link[flapper]{dcpf_simplify}} (see \code{\link[flapper]{dcpf-class}}). At a minimum, this should contain a unique identifier for each path (named `path_id') and the probability associated with each cell along each path (`cell_pr').
#' @param ... Additional arguments (none implemented).
#' @details For each path, at each time step the probability associated with the sampled location depends on a user-defined movement model that is driven by the distance between the sampled locations for the individual at the previous and current time steps (and other user-defined parameters). This function simply sums the logarithms of these probabilities for each path as a measure of their relative likelihood, given the movement model.
#' @examples
#' # An example with the example DCPF paths dataset included in flapper
#' dcpf_loglik(dat_dcpf_paths)
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
#### dcpf_plot_history()

#' @title Plot particle histories from the DCPF algorithm
#' @description This function plots the spatiotemporal particle histories from the depth-contour particle filtering (DCPF) algorithm. This produces, for each time step, a map of the individual's possible locations, solely based on its depth (including the bathymetry and the measurement error), with sampled locations overlaid.
#' @param record A \code{\link[flapper]{.dcpf-class}} object from \code{\link[flapper]{dcpf}} that contains particle histories.
#' @param time_steps An integer vector that defines the time steps for which to plot particle histories.
#' @param add_bathy A named list, passed to \code{\link[prettyGraphics]{pretty_map}}, to customise the appearance of the bathymetry surface on each map.
#' @param add_particles A named list, passed to \code{\link[prettyGraphics]{pretty_map}}, to customise the appearance of the particles on each map.
#' @param forwards A logical variable that defines whether or not create plots forwards (i.e., from the first to the last \code{time_steps}) or backwards (i.e., from the last to the first \code{time_steps}).
#' @param prompt A logical input that defines whether or not to pause between plots (\code{prompt = TRUE}).
#' @param ... Plot customisation arguments passed to \code{\link[prettyGraphics]{pretty_map}}.
#'
#' @examples
#' #### Implement dcpf() algorithm
#' # Here, we use pre-defined outputs for speed
#'
#' #### Example (1): The default implementation
#' dcpf_plot_history(dat_dcpf_histories, time_steps = 1)
#'
#' #### Example (2): Plot customisation options, e.g.:
#' # Customise bathy via add_bathy()
#' dcpf_plot_history(dat_dcpf_histories,
#'                   time_steps = 1,
#'                   add_bathy = list(col = c(grDevices::topo.colors(2))))
#' # Customise particles via add_particles
#' dcpf_plot_history(dat_dcpf_histories,
#'                   time_steps = 1,
#'                   add_particles = list(col = "red"))
#' # Pass other arguments to prettyGraphics::pretty_map() via ...
#' dcpf_plot_history(dat_dcpf_histories,
#'                   time_steps = 1,
#'                   add_polys = list(x = dat_coast, col = "brown"),
#'                   crop_spatial = TRUE)
#'
#' #### Example (3): Plot multiple time steps
#' pp <- graphics::par(mfrow = c(2, 2))
#' dcpf_plot_history(dat_dcpf_histories, time_steps = 1:4, prompt = FALSE)
#' graphics::par(pp)
#'
#' @return The function returns a plot, for each time step, of all the possible locations of the individual, given the \code{archival}, \code{bathy} and \code{calc_depth_error} arguments to code{\link[flapper]{dcpf}} (marked [1]), with sampled locations overlaid.
#'
#' @seealso \code{\link[flapper]{dcpf}} implements the DCPF algorithm. \code{\link[flapper]{dcpf_simplify}} assembles paths from particle histories. \code{\link[flapper]{dcpf_plot_1d}}, \code{\link[flapper]{dcpf_plot_2d}} and \code{\link[flapper]{dcpf_plot_3d}} provide plotting routines for paths. \code{\link[flapper]{dcpf_loglik}} calculates the log-likelihood of each path.
#'
#' @author Edward Lavender
#' @export

dcpf_plot_history <- function(record,
                              time_steps = 1:length(history),
                              add_bathy = list(),
                              add_particles = list(pch = "."),
                              forwards = TRUE,
                              prompt = TRUE,...){
  if(!inherits(record, ".dcpf")) stop("'record' must be a '.dcpf' class object.")
  depth            <- record$args$archival$depth
  bathy            <- record$args$bathy
  calc_depth_error <- record$args$calc_depth_error
  history          <- record$history
  time_steps       <- sort(time_steps)
  if(!forwards) time_steps <- rev(time_steps)
  lapply(time_steps, function(t){
    depth_lwr <- depth[t] + calc_depth_error(depth[t])[1]
    depth_upr <- depth[t] + calc_depth_error(depth[t])[2]
    title <- paste0("Time ", t, " [", round(depth_lwr, 2), ": ", round(depth_upr, 2), " (m)]")
    r <- bathy >= depth_lwr & bathy <= depth_upr
    add_bathy$x <- r
    xy_t <- raster::xyFromCell(bathy, history[[t]]$id_current)
    add_particles$x <- xy_t[, 1]
    add_particles$y <- xy_t[, 2]
    prettyGraphics::pretty_map(r,
                               add_rasters = add_bathy,
                               add_points = add_particles,
                               main = title,
                               verbose = FALSE,...)
    if(prompt * length(time_steps) > 1) readline(prompt = "Press [enter] to continue or [Esc] to exit...")
  })
  return(invisible())
}


########################################
########################################
#### dcpf_plot_1d()

#' @title Plot one-dimensional depth time series from the DCPF algorithm
#' @description This function plots the observed depth time series and the depth time series associated with each path reconstructed by the depth-contour particle filtering (DCPF) algorithm.
#' @param archival A dataframe of depth (m) observations named `depth', as used by \code{\link[flapper]{dcpf}}.
#' @param paths A dataframe containing reconstructed movement path(s) from \code{\link[flapper]{dcpf}} via \code{\link[flapper]{dcpf_simplify}} (see \code{\link[flapper]{dcpf-class}}). At a minimum, this should contain a unique identifier for each path (named `path_id’), timesteps (`timestep') and the depth associated with each cell along each path (`cell_z').
#' @param scale A number that vertically scales the depth time series for the observations and the reconstructed path(s). By default, absolute values for depth are assumed and negated for ease of visualisation.
#' @param pretty_axis_args,xlab,ylab,type,... Plot customisation arguments passed to \code{\link[prettyGraphics]{pretty_plot}}.
#' @param add_lines A named list, passed to \code{\link[graphics]{lines}}, to customise the appearance of the depth time series for reconstructed path(s).
#' @param prompt A logical input that defines whether or not plot the observed depth time series with each reconstructed depth time series on a separate plot, sequentially, with a pause between plots (\code{prompt = TRUE}), or with all reconstructed time series on a single plot (\code{prompt = FALSE}).
#' @details Observed and reconstructed depth time series can differ due to measurement error, which is controlled via the \code{calc_depth_error} function in the DCPF algorithm (see \code{\link[flapper]{dcpf}}).
#' @examples
#' #### Implement dcpf() algorithm
#' # Here, we use pre-defined outputs for speed
#' archival <- dat_dcpf_histories$args$archival
#' paths    <- dat_dcpf_paths
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
#'   pp <- graphics::par(mfrow = c(3, 4))
#'   dcpf_plot_1d(depth, paths, prompt = TRUE)
#'   graphics::par(pp)
#' }
#'
#' @return The function returns a plot of the observed and reconstructed depth time series, either for all paths at once (if \code{prompt = FALSE}) or each path separately (if \code{prompt = TRUE}).
#' @seealso \code{\link[flapper]{dcpf}} implements the DCPF algorithm. \code{\link[flapper]{dcpf_plot_history}} visualises particle histories and \code{\link[flapper]{dcpf_simplify}} processes the outputs into a dataframe of movement paths. \code{\link[flapper]{dcpf_plot_1d}}, \code{\link[flapper]{dcpf_plot_2d}} and \code{\link[flapper]{dcpf_plot_3d}} provide plotting routines for paths. \code{\link[flapper]{dcpf_loglik}} calculates the log-likelihood of each path.
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
  # Checks
  check_names(input = archival, req = c("depth"))
  check_names(input = paths, req = c("path_id", "timestep", "cell_z"), type = all)
  # Drop origin, if supplied, to avoid alignment issues
  # paths <- paths[paths$timestep != 0, ]
  if(nrow(paths) < 1) stop("'paths' does not contain any timesteps post-origin.")
  # Make plots
  dots <- list(...)
  if(!prompt) prettyGraphics::pretty_plot(1:nrow(archival), archival$depth*scale,
                                          pretty_axis_args = pretty_axis_args,
                                          xlab = xlab, ylab = ylab,
                                          type = type,...)
  lapply(split(paths, paths$path_id), function(d){
    if(prompt) prettyGraphics::pretty_plot(archival$depth*scale,
                                           pretty_axis_args = pretty_axis_args,
                                           xlab = xlab, ylab = ylab,
                                           type = type,...)
    add_lines$x <- d$timestep
    add_lines$y <- d$cell_z*scale
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
#' @param paths A dataframe containing reconstructed movement path(s) from \code{\link[flapper]{dcpf}} via \code{\link[flapper]{dcpf_simplify}} (see \code{\link[flapper]{dcpf-class}}). At a minimum, this should contain a unique identifier for each path (named `path_id') and the x and y coordinates that define the trajectory of each path (`cell_x' and `cell_y').
#' @param bathy A \code{\link[raster]{raster}} of the bathymetry over which movement was reconstructed.
#' @param add_paths A named list, passed to \code{\link[prettyGraphics]{add_sp_path}}, to customise the appearance of the paths.
#' @param prompt A logical input that defines whether or not plot each path on a separate plot, sequentially, with a pause between plots (\code{prompt = TRUE}), or all paths on a single plot (\code{prompt = FALSE}).
#' @param ... Additional arguments, passed to \code{\link[prettyGraphics]{pretty_map}}, for plot customisation.
#' @return The function maps the trajectories of reconstructed paths across the bathymetry surface, returning a single map if \code{prompt = FALSE} or one map for each path if \code{prompt = TRUE}.
#' @examples
#' #### Implement dcpf() algorithm
#' # Here, we use pre-defined outputs for speed
#' bathy <- dat_dcpf_histories$args$bathy
#' paths <- dat_dcpf_paths
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
#'   pp <- graphics::par(mfrow = c(3, 4))
#'   dcpf_plot_2d(paths, bathy, add_paths = list(length = 0.01),
#'                prompt = TRUE, verbose = FALSE)
#'   graphics::par(pp)
#' }
#'
#' @seealso \code{\link[flapper]{dcpf}} implements the DCPF algorithm. \code{\link[flapper]{dcpf_plot_history}} visualises particle histories and \code{\link[flapper]{dcpf_simplify}} processes these into a dataframe of movement paths. \code{\link[flapper]{dcpf_plot_1d}}, \code{\link[flapper]{dcpf_plot_2d}} and \code{\link[flapper]{dcpf_plot_3d}} provide plotting routines for paths. For mapping, it can be useful to interpolate shortest (least-cost) paths between sequential locations via \code{\link[flapper]{lcp_interp}}. \code{\link[flapper]{dcpf_loglik}} calculates the log-likelihood of each path.
#' @author Edward Lavender
#' @export
#'

dcpf_plot_2d <- function(paths,
                         bathy,
                         add_paths = list(),
                         prompt = FALSE,...){
  check_names(input = paths, req = c("path_id", "cell_x", "cell_y"))
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
#' @param paths A dataframe containing reconstructed movement path(s) from \code{\link[flapper]{dcpf}} via \code{\link[flapper]{dcpf_simplify}} (see \code{\link[flapper]{dcpf-class}}). At a minimum, this should contain a unique identifier for each path (named `path_id') and the x, y and z coordinates that define the trajectory of each path (`cell_x', `cell_y' and `cell_z').
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
#' # Here, we use pre-defined outputs for speed
#' # Note that it may be beneficial to interpolate paths between points
#' # ... e.g., via lcp_interp() prior to plotting, but we will not do that here.
#' bathy <- dat_dcpf_histories$args$bathy
#' paths <- dat_dcpf_paths
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
#' @seealso \code{\link[flapper]{dcpf}} implements the DCPF algorithm. \code{\link[flapper]{dcpf_plot_history}} visualises particle histories and \code{\link[flapper]{dcpf_simplify}} processes these into a dataframe of movement paths. \code{\link[flapper]{dcpf_plot_1d}}, \code{\link[flapper]{dcpf_plot_2d}} and \code{\link[flapper]{dcpf_plot_3d}} provide plotting routines for paths. For mapping, it can be useful to interpolate shortest (least-cost) paths between sequential locations via \code{\link[flapper]{lcp_interp}}. \code{\link[flapper]{dcpf_loglik}} calculates the log-likelihood of each path.
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
  # Checks
  if(!requireNamespace("plotly", quietly = TRUE)) stop("This function requires the 'plotly' package. Please install it with `install.packages('plotly')` first.")
  check_names(input = paths, req = c("path_id", "cell_x", "cell_y", "cell_z"))
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
