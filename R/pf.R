########################################
########################################
#### pf()

#' @title The particle filtering routine
#' @importFrom data.table :=
#' @description This function implements a simulation-based particle filtering process for the reconstruction of animal movement paths. This extends the AC, DC and ACDC algorithms in \code{\link[flapper]{flapper}} (\code{\link[flapper]{ac}}, \code{\link[flapper]{dc}} and \code{\link[flapper]{acdc}}), which determine the possible locations of an individual through time based on acoustic detections and/or depth contours, by the incorporation of a movement model that connects a subset of the animal's possible locations between time steps into movement paths.
#'
#' To implement this approach, a list of \code{\link[raster]{raster}} layers that record the possible locations of the individual at each time step (or a list pointers to these layers) needs to be supplied. A starting location (\code{origin}) can be supplied to constrain the initial set of sampled locations of the individual. At each time step, \code{n} possible locations (`particles') are sampled (with replacement) from the set of possible locations. For each (\code{1:n}) particle, a movement model is used to simulate where the individual could have moved to at the next time step, if it was in any of those locations. In the current framework, the probability of movement into surrounding cells depends on the distance to those cells, which can be represented as using Euclidean or least-cost distances depending on the distance method (\code{calc_distance}), and user-defined movement models (\code{calc_movement_pr_from_origin} and \code{calc_movement_pr}) that link distances to movement probabilities at each time step.
#'
#' At each subsequent time step, this process repeats, with \code{n} possible locations of the individual sampled according to the probability that the individual could have been in that cell, given a previously sampled location and the set of possible locations at the next time step. The result is a set of locations at each time step that are consistent with the data and model parameters. Sampled locations can be connected into movement paths via \code{\link[flapper]{pf_simplify}}.
#'
#' @param record A list of \code{\link[raster]{raster}} layers, or a named list of source files that can be sequentially loaded with \code{\link[raster]{writeRaster}}, that define the set of possible locations of the individual at each time step. This can be derived from the `record' elements in \code{\link[flapper]{ac}}, \code{\link[flapper]{dc}} or \code{\link[flapper]{acdc}}. As in these algorithms, for each \code{\link[raster]{raster}}, the coordinate reference system should be the Universal Transverse Mercator projection. Raster resolution needs to be sufficiently high such that an individual could transition between cells in the time steps between sequential layers (see \code{calc_movement_pr}). If the `shortest distances' method is used for distance calculations (i.e., \code{calc_distance = "lcp"}, see below), then the resolution of the surface in x and y directions should also be equal (for surface's with unequal horizontal resolution, \code{\link[raster]{resample}} can be used to equalise resolution: see \code{\link[flapper]{lcp_over_surface}}). For computational efficiency, it is beneficial if inputs are cropped to the smallest possible area, with any areas in which movement is impossible (e.g., on land for benthic animals) set to NA (see \code{\link[raster]{crop}}, \code{\link[raster]{mask}} and the processing implemented by \code{\link[flapper]{lcp_over_surface}}). It may also be desirable to aggregate high-resolution \code{\link[raster]{raster}}s (see \code{\link[flapper]{process_surface}}).
#' @param data (optional) A dataframe, with one row for each time step (i.e. element in \code{record}), and columns for any time-specific variables included in the movement models (see \code{calc_movement_pr} and \code{calc_movement_pr_from_origin}). If unsupplied, movement probabilities are constant through time.
#' @param origin (optional) A matrix that defines the coordinates (x, y) of the individual's initial location. Coordinates should follow the restrictions for \code{record} and lie on a plane (i.e., Universal Transverse Mercator projection).
#' @param calc_distance A character that defines the method used to calculate distances between a point (i.e., a sampled location) and the surrounding cells. This drives the probability of movement into those cells via a movement model (see \code{calc_movement_pr} and \code{calc_movement_pr_from_origin}). Currently supported options are Euclidean distances (\code{"euclid"}) or shortest (least-cost) distances (\code{"lcp"}), which represent the shortest distance that an individual would have to move over the surface (\code{bathy}, see below) to traverse between locations (accounting for both planar and vertical distances). Note that this option requires that resolution of \code{bathy} in the x and y directions is equal. At small spatial scales, this option provides more realistic distances in hilly landscapes, but it is more computationally expensive. At larger scales, horizontal distances tend to overwhelm vertical distances, so Euclidean distances may be acceptable. A pragmatic option is to implement the algorithm (possibly for a subset of the \code{record} time series) using the Euclidean distances method and then interpolate least-cost paths between (a) the subset of sampled locations returned by this approach via \code{\link[flapper]{pf_simplify}} or (b) (b) the subset of paths returned by \code{\link[flapper]{pf_simplify}} via \code{\link[flapper]{lcp_interp}}. This two-step approach will demonstrate whether sequential positions are plausible (i.e., not too far apart) once the bathymetry is taken into account. If so, the shortest-distances derived using this method can be used for post-hoc adjustment of movement probabilities. Alternatively, this approach may demonstrate that the algorithm should be re-implemented using the shortest distances method (see \code{\link[flapper]{lcp_interp}}).
#' @param calc_distance_euclid_fast If \code{calc_distance = "euclid"}, \code{calc_distance_euclid_fast} is a logical input that defines whether or not to use the `fast' method for distance and probability calculations. Under the `slow' method (\code{calc_distance_euclid_fast = FALSE}), at each time step the algorithm iterates over each particle, calculates the distance around that particle to neighbouring cells, converts distances to movement probabilities, combines these with the intrinsic probabilities associated with each cell (assigned by the AC, DC or ACDC algorithm) and samples future locations accordingly. Since each particle is considered in turn, each particle has a `history' that is `remembered', which makes path assembly relatively straightforward. In contrast, the faster implementation considers all particles at the same time: a single distance/movement probability surface is calculated around sampled particles, from which future particles are sampled. This approach really excels as the number of particles increases (see \code{n}). The disadvantage is that particles do not `remember' where they have come from and as a result movement path assembly is more involved (see \code{\link[flapper]{pf_simplify}}). However, in most cases \code{calc_distance_euclid_fast} is probably desirable, since \code{\link[flapper]{pf_simplify}} takes care of movement path assembly.
#' @param bathy (optional) If \code{calc_distance = "lcp"}, \code{bathy} is \code{\link[raster]{raster}} that defines the bathymetry across the area within which the individual could have moved. \code{bathy} must be planar (i.e., Universal Transverse Mercator projection) with units of metres in x, y and z directions (m). The surface's resolution is taken to define the distance between horizontally and vertically connected cells and must be the same in both x and y directions. Any cells with NA values (e.g., due to missing data) are treated as `impossible' to move though by the algorithm (see \code{\link[flapper]{lcp_over_surface}}).
#' @param calc_distance_graph (optional) If \code{calc_distance = "lcp"}, \code{calc_distance_graph} is a graph object from \code{\link[flapper]{lcp_graph_surface}} that defines the distances between connected cells in \code{bathy}. This is useful in iterative applications in the same environment, especially for large bathymetry rasters (see \code{bathy}) because the calculations of movement costs between adjacent cells and graph construction---both required for the calculation of shortest paths---can be skipped.
#' @param calc_movement_pr The movement model. This must be a function that calculates the probability of movement between two locations in the time between sequential \code{record}s, given the distance between them (the first argument) and (optionally) any other pertinent information in the \code{data} for the relevant time step (the second argument). For example, behavioural states could be defined in a column in the \code{data} and then included as part of the movement model. Both arguments must be accepted, even if the second is ignored. The default option is a declining logistic curve, designed for flapper skate (\emph{Dipturus intermedius}), representing a high probability of movement between nearby locations and a lower probability of movement between distant locations (see \code{\link[flapper]{pf_setup_movement_pr}}). For computational efficiency, it is beneficial if the probability of movement is set to zero beyond a distance considered to be highly unlikely, because such locations are excluded from subsequent calculations (see \code{\link[flapper]{pf_setup_movement_pr}} and the \code{mobility} argument, below). Currently, the movement model cannot incorporate the individual's location (for spatially variable fields such as water currents).
#' @param calc_movement_pr_from_origin (optional) If an \code{origin} is supplied, \code{calc_movement_pr_from_origin} can be supplied to define the probability of sampling possible starting locations depending on their distance from the \code{origin}. By default, \code{calc_movement_pr_from_origin = calc_movement_pr} and the same guidance applies to both arguments. Specifying \code{calc_movement_pr_from_origin} specifically may be necessary if the duration between the time at which the \code{origin} was observed and the first \code{record} differs from the regular interval between subsequent \code{record}s. This allows the effect of the \code{origin} to be weaker or stronger than the effect of locations at later time steps on sampled locations, if necessary.
#' @param mobility (optional) A number that defines the maximum horizontal distance (m) that the individual could travel in the time period between sequential \code{\link[raster]{raster}} layers. While this is optional, it is often computationally beneficial to define \code{mobility} because this restricts distance and probability calculations at each time step within the smallest appropriate range (rather than across the full surface).
#' @param mobility_from_origin (optional) As above for \code{mobility}, but for the first time step (i.e., the horizontal movement the individual could have moved from the \code{origin}) (see \code{calc_movement_pr_from_origin}).
#' @param n An integer that defines the number of particles (i.e., the number of locations sampled at each time step from the set of possible locations at that time step).
#' @param resample (optional) An integer that defines the minimum number of unique cells that should be sampled, given the movement model (\code{calc_movement_pr} or \code{calc_movement_pr_from_origin}). If supplied, if fewer than \code{resample} unique cells are sampled at a time step, \code{n} particles are re-sampled with replacement with equal probability from all of the cells with non-zero probability. This may facilitate algorithm convergence if there are some cells that are overwhelmingly more probable (given the movement model) are `dead ends': re-sampling all possible cells with equal probability allows the algorithm to explore less likely routes more easily when the number of routes becomes low. \code{resample} must be less than or equal to \code{n}.
#' @param update_history (optional) A list that defines particle histories from an earlier implementation of \code{\link[flapper]{pf}} (i.e., the `history' element of a \code{\link[flapper]{pf_archive-class}} object). If provided, the function attempts to continue paths from an earlier time step (see \code{update_history_from_time_step}). This can be useful if the algorithm fails to converge on its initial implementation because the algorithm can be re-started part-way through the time series.
#' @param update_history_from_time_step If \code{update_history} is provided, \code{update_history_from_time_step} is an integer that defines the time step (i.e., element in \code{update_history}) from which to restart the algorithm. If provided, the algorithm continues from this point, by taking the starting positions for \code{n} particles from \code{update_history[[update_history_from_time_step]]$id_current}.
#' @param write_history (optional) A named list, passed to \code{\link[base]{saveRDS}}, to save a dataframe of the sampled particles at each time step to file. The `file' argument should be the directory in which to save files. Files are named by time steps as `pf_1', `pf_2' and so on.
#' @param cl,use_all_cores (optional) Parallelisation options. These can be implemented for the approaches that consider particles iteratively (i.e., \code{calc_distance = "euclid"} with \code{calc_distance_euclid_fast = FALSE}, or \code{calc_distance = "lcp"}). The algorithm can be parallelised within time steps over (1) particles via \code{cl} (an integer defining the number of child processes (ignored on Windows) or a cluster object created by \code{\link[parallel]{makeCluster}} (see \code{\link[pbapply]{pblapply}})) or (2) within particles for the calculation of shortest distances (if \code{calc_distance = "lcp"}) via a logical input (\code{TRUE}) to \code{use_all_cores}. For \code{calc_distance = "euclid"}, parallelisation is typically only beneficial for relatively large numbers of particles, because the substantial computation overhead associated with parallelisation across paths at each time step is substantial. For \code{calc_distance = "lcp"}, \code{use_all_cores = TRUE} is typically a better option for parallelisation; however, this may only be beneficial if \code{mobility} relatively large. At present, parallelisation is not very effective, especially if a cluster is supplied, and may be slower.
#' @param seed (optional) An integer to define the seed for reproducible simulations (see \code{\link[base]{set.seed}}).
#' @param verbose A logical variable that defines whether or not to print messages to the console or to file to relay function progress. If \code{con = ""}, messages are printed to the console; otherwise, they are written to file (see below).
#' @param con If \code{verbose = TRUE}, \code{con} is character string that defines how messages relaying function progress are returned. If \code{con = ""}, messages are printed to the console (unless redirected by \code{\link[base]{sink}}). Otherwise, \code{con} defines the full pathway to a .txt file (which can be created on-the-fly) into which messages are written to relay function progress. This approach, rather than printing to the console, is recommended for clarity, speed and debugging.
#' @param optimisers A named list of optimisation controls from \code{\link[flapper]{pf_setup_optimisers}}.
#'
#' @details
#' \subsection{Background}{This function implements a widely applicable particle simulation and filtering based approach to refine maps of possible locations of an individual through time via the incorporation of a movement model that facilitates the reconstruction of movement paths. Within \code{\link[flapper]{flapper}}, the acoustic-container (AC), depth-contour (DC) and acoustic-container depth-contour (ACDC) algorithms, which define the possible locations of an individual through time based on acoustic containers and/or depth contours, can be passed through a this process, resulting in the DCPF, ACPF and ACDCPF algorithms.
#' \itemize{
#'   \item \strong{ACPF}. The ACPF algorithm combines the AC algorithm with particle filtering. This is designed to simulate possible movement paths and emergent patterns of space use in passive acoustic telemetry arrays. This algorithm is widely applicable.
#'   \item \strong{DCPF}. The DCPF algorithm combines the DC algorithm with particle filtering. This is designed to simulate possible movement paths of a tagged animal over the seabed, given a regular sequence of depth observations (\code{archival}), the bathymetry (\code{bathy}) over which movement occurred and a movement model (\code{calc_movement_pr}) that specifies the probability of movement from a given location to any other, given the distance between them (and any other pre-defined time-dependent parameters in \code{archival}). The function was motivated by small scale applications concerning the reconstruction of possible movement paths of flapper skate (\emph{Dipturus intermedius}) tagged with archival tags, following capture and release in a given location, for short-periods of time post-release.
#'   \item \strong{ACDCPF}. The ACDCPF algorithm combines the ACDC algorithm with particle filtering. For tagged animals with acoustic and archival data, this algorithm is designed to reconstruct fine-scale movement paths and emergent patterns of space use across an area.
#'   }
#' }
#'
#' \subsection{Methods}{At the first time step, \code{n} starting points (`particles') are selected from the set of possible locations of the individual. If an \code{origin} is specified, this selection can be biased towards cells near the origin by the movement model. From each starting position, the Euclidean or shortest distances to cells in which the individual could have been located at the next time step are calculated and passed to a movement model that assigns movement probabilities to each cell. Since movement probabilities are likely to be behaviourally dependent, the movement model can depend on time step-specific information specified in \code{data}. However, currently, the model cannot depend on sampled locations; therefore, at least under some conditions, the movement model will need to reflect the maximum distance that an individual could travel within the time period between depth observations, accounting for the possible effects of water currents and any other influences on swimming speed. Movement probabilities are combined with the 'intrinsic' probability associated with each location (as defined in the \code{record}), giving a holistic measure of the probability of movement into each cell. From the set of cells in which the individual could have been located with a probability of more than zero, \code{n} particles are sampled, with replacement and according to their probability, and taken as possible starting positions at the next time step. This process repeats until the end of the time series. Since locations are sampled from the set of possible locations at each time step, any restrictions on the individual's movement (e.g., from acoustic detections) incorporated within \code{record} are directly incorporated. The result is a set of simulated particles that represent possible locations of the individual at each time step, accounting for movement restrictions. This can be assembled into a set of movement paths over a surface that is consistent with the data and the model parameters via \code{\link[flapper]{pf_simplify}}. While the number of particles is predetermined by \code{n}, more than \code{n} possible pathways may be defined by all of the combinations of sequential particle movements.}
#'
#' \subsection{Convergence}{Algorithm convergence is not guaranteed. There are three main circumstances in which the algorithm may fail to return any paths that span the start to the end of the depth time series:
#'  \enumerate{
#'    \item \strong{Chance.} All \code{n} paths may be `dead ends'. This possibility can be mitigated by increasing \code{n}.
#'    \item \strong{Movement model.} The movement model may be too limiting. This possibility can be mitigated by ensuring that the movement model realistically represents the probability that an individual can transition between cells given the distance between them. This may be guided by data on the study species or similar species. The movement model may need to account for the effect of water currents, which may increase maximum `swimming' speeds in some directions. (Unfortunately, spatially variable swimming speeds are not currently implemented.) If maximum swimming speeds are uncertain, implementing the algorithm over longer time series (e.g., every \eqn{2^{nd}} \code{record}), with a suitably relaxed movement model, may facilitate convergence if maximum speeds are unlikely to be maintained for long periods.
#'    \item \strong{Other assumptions.} The particle filtering process is based on pre-defined surfaces of the possible locations of the individual at each time step and the assumptions in the computation of these surfaces may be violated. For example, in the DC and ACDC algorithms. the depth error may be too restrictive, given the accuracy of the depth observations, the bathymetry data and the tidal height across an area.
#'  }
#' In these scenarios, the function returns a message that it is about to fail and the results from the start of the algorithm until the current time step, before stopping.
#' }
#'
#' \subsection{Computational considerations}{ This algorithm is computationally intensive. It is advisable that it is run initially with a small time series in a small area and a small number of particles. For larger datasets, there are some tricks that can improve computation time.
#'   \itemize{
#'     \item \strong{Temporal resolution.} Reduce the temporal resolution of the \code{record} time series so that there are fewer time steps.
#'     \item \strong{Bathymetric resolution.} Reduce the resolution of \code{record}s. For the DC or ACDC algorithms, this will require re-implementing \code{\link[flapper]{dc}} or \code{\link[flapper]{acdc}} with a lower resolution surface and propagating the additional error induced by this process via \code{calc_depth_error} (e.g., see \code{\link[flapper]{process_surface}}) and then using the recomputed \code{record}s in this function.
#'     \item \strong{Mobility limits.} Check whether or not setting mobility limits (\code{mobility}, \code{mobility_from_origin}) improves speed.
#'     \item \strong{Distance calculations}. This step is particularly slow because the distance from each location to many or all surrounding locations are calculated. To speed up this step, implement \code{calc_distance = "euclid"} with \code{calc_distance_euclid_fast = TRUE}. If necessary, interpolate least-cost paths after algorithm completion within \code{\link[flapper]{pf_simplify}} or \code{\link[flapper]{lcp_interp}}. In the future, a Markov chain Monte Carlo style approach may be implemented in which distances (and probabilities) for randomly selected `proposal' cells are calculated, with those cells then rejected or retained, rather than calculating distances to many or all surrounding cells, but this is unlikely to be faster in many settings.
#'     \item \strong{Parallelisation.} If the fast Euclidean distances method is not used, test alternative options for parallelisation. For the \code{cl} argument, specifying an integer on non-Windows platforms may be faster than a cluster from \code{\link[parallel]{makeCluster}} (see \code{\link[pbapply]{pblapply}}). Parallelisation may be slower in some circumstances.
#'   }
#' }
#'
#' @examples
#' #### Summary
#' # In these examples, we consider an example depth time series, which we generate.
#' # We use the DC algorithm to reconstruct the possible locations of the individual through time.
#' # We then use particle filtering (the DCPF algorithm in this case) to refine maps of the
#' # ... individuals possible location via the incorporation of a movement model. The
#' # ... same principles apply to other algorithms e.g., AC and ACDC.
#'
#' #### Step (1): Generate some example movement time series over a surface
#'
#' ## Sample species
#' # In this example, we consider flapper skate (Dipturus intermedius)
#' # ... off the west coast of Scotland.
#'
#' ## Define a starting location (optional) in UTM coordinates
#' xy <- matrix(c(708886.3, 6254404), ncol = 2)
#' ## Define 'observed' depth time series using absolute values
#' # Imagine these are observations made every two minutes
#' depth <- c(
#'   163.06, 159.71, 153.49, 147.04, 139.86, 127.19, 114.75,
#'   99.44, 87.01, 78.16, 70.03, 60.23, 49.96, 35.39,
#'   27.75, 20.13, 12.73, 11.32
#' )
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
#' surface <- dat_gebco
#' boundaries <- raster::extent(707884.6, 709884.6, 6253404, 6255404)
#' blank <- raster::raster(boundaries, res = c(5, 5))
#' surface <- raster::resample(surface, blank)
#'
#' ## Define depth error function
#' # Because the bathymetry data is very coarse, and the bathymetry is
#' # ... very complex in this region of Scotland, we have to
#' # ... force a high depth error to be high in this example.
#' # The calc_depth_error() function can depend on depth, but in this example
#' # ... we assume the depth error is independent of depth.
#' cde <- function(...) matrix(c(-30, 30), nrow = 2)
#'
#' ## Define movement model
#' # The default movement model is suitable, with skate moving typically
#' # ... less than 200 m in a two-minute period.
#' # You could use a separate movement model (and mobility restriction) for the origin
#' # ... if necessary, but for brevity we don't implement that here.
#'
#' ## Visualise movement surface, with starting location overlaid
#' prettyGraphics::pretty_map(
#'   add_rasters = list(x = surface),
#'   add_points = list(x = xy),
#'   verbose = FALSE
#' )
#'
#' #### Step (2): Use the DC algorithm to get individual's the possible locations
#' dc_out <- dc(
#'   archival = depth,
#'   bathy = surface,
#'   calc_depth_error = cde,
#'   save_record_spatial = NULL
#' )
#' # Extract time-specific maps
#' # ... Either directly for single chunk implementations of dc()
#' record <- dc_out$archive$record
#' record <- lapply(record, function(r) r$spatial[[1]]$map_timestep)
#' # ... Or indirectly via acdc_simplify() in general
#' dc_out_summary <- acdc_simplify(dc_out, type = "dc")
#' record <- lapply(dc_out_summary$record, function(r) r$spatial[[1]]$map_timestep)
#' # Plot maps
#' lapply(record, function(r) raster::plot(r))
#'
#' #### Example (1): Implement algorithm using default options
#' out_1 <- pf(
#'   record = record,
#'   origin = xy,
#'   calc_distance = "euclid",
#'   calc_movement_pr = pf_setup_movement_pr,
#'   n = 10L,
#'   seed = 1
#' )
#' # The function returns a pf_archive-class object
#' class(out_1)
#' utils::str(out_1)
#' # Algorithm duration during testing ~0.03 minutes
#' # ... (vs ~0.21 minutes with calc_distance_euclid_fast = FALSE)
#'
#' #### Example (2): Implement algorithm reading rasters from file on the fly
#' # Save files (this can be done directly within dc())
#' tmp_root <- paste0(tempdir(), "/dc_files/")
#' dir.create(tmp_root)
#' lapply(1:length(record), function(i) raster::writeRaster(record[[i]], paste0(tmp_root, i, ".tif")))
#' # Create pointer to files (in the correct order by time step)
#' record_pointer <- list.files(tmp_root)
#' record_pointer <- substr(record_pointer, 1, nchar(record_pointer) - 4)
#' record_pointer <- data.frame(index = 1:length(record_pointer), name = record_pointer)
#' record_pointer <- record_pointer[order(as.integer(record_pointer$name)), ]
#' record_pointer <- list.files(tmp_root, full.names = TRUE)[record_pointer$index]
#' out_2 <- pf(
#'   record = record_pointer,
#'   origin = xy,
#'   calc_distance = "euclid",
#'   calc_movement_pr = pf_setup_movement_pr,
#'   n = 10L,
#'   seed = 1
#' )
#'
#' #### Example (3): Implement a blanket mobility restriction
#' # This can improve computational efficiency but offers no improvement
#' # ... in this example (e.g., due to relatively small raster, so the cost of
#' # ... focusing in on a specific area matches the speed gains of
#' # ... implementing the distance/movement
#' # ... pr calculations over a smaller area at each time step)
#' out_3 <- pf(
#'   record = record,
#'   origin = xy,
#'   calc_distance = "euclid",
#'   calc_movement_pr = pf_setup_movement_pr,
#'   mobility = 250,
#'   n = 10L,
#'   seed = 1
#' )
#' # Algorithm duration during testing ~0.04 minutes
#'
#' #### Example (4): Implement algorithm using shortest distances
#' # Note the need for a surface with equal resolution if this option is implemented.
#' # To speed up the initial stages of the algorithm, you can supply the graph
#' # ... required for least-cost calculations via calc_distance_graph, but for
#' # ... brevity we don't implement that here.
#' out_4 <- pf(
#'   record = record,
#'   origin = xy,
#'   calc_distance = "lcp",
#'   bathy = surface,
#'   calc_movement_pr = pf_setup_movement_pr,
#'   n = 10L,
#'   seed = 1
#' )
#' # This option is slower: algorithm duration during testing ~0.73 minutes
#'
#' #### Example (5): Implement algorithm using shortest distances with mobility restriction
#' out_5 <- pf(
#'   record = record,
#'   origin = xy,
#'   calc_distance = "lcp",
#'   bathy = surface,
#'   calc_movement_pr = pf_setup_movement_pr,
#'   mobility = 250,
#'   n = 10L,
#'   seed = 1
#' )
#' # Algorithm duration during testing ~0.63 minutes
#' # With shortest distances, the mobility restriction makes more difference
#' # ... in improving the computation time, because calculations are more involved.
#'
#' #### Example (6): Parallelisation for Euclidean distances is via cl argument
#' # ... which implements parallelisation across paths. This is only implemented
#' # ... for calc_distance_euclid_fast = FALSE, which is rarely desirable.
#' out_6 <- pf(
#'   record = record,
#'   origin = xy,
#'   calc_distance = "euclid",
#'   calc_distance_euclid_fast = FALSE,
#'   calc_movement_pr = pf_setup_movement_pr,
#'   mobility = 250,
#'   n = 10L,
#'   cl = parallel::makeCluster(2L),
#'   seed = 1
#' )
#'
#' #### Example (7): Parallelisation for shortest distances is usually best
#' # ... via use_all_cores = TRUE
#' # However, the benefits of parallelisation depend on the number of least-cost
#' # ... paths calculations that need to be performed, which depends on the size
#' # ... of bathy, and may be minimal (or negative) for small areas.
#' out_7 <- pf(
#'   record = record,
#'   origin = xy,
#'   calc_distance = "lcp",
#'   bathy = surface,
#'   calc_movement_pr = pf_setup_movement_pr,
#'   mobility = 250,
#'   n = 10L,
#'   use_all_cores = TRUE,
#'   seed = 1
#' )
#' # But the speed benefits in this case are minimal.
#' # Algorithm duration during testing ~0.61 minutes.
#'
#' #### Example (8): Extend the movement model via variables in data
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
#' setup_movement_pr <- function(distance, data) {
#'   beta <- -0.05 + -0.005 * data$va_abs
#'   pr <- stats::plogis(10 + distance * beta)
#'   pr[distance > 500] <- 0
#'   return(pr)
#' }
#' ## Examine the movement model with distance and VA
#' pr <- setup_movement_pr(1:1000, data.frame(va_abs = 0))
#' prettyGraphics::pretty_plot(pr, type = "n")
#' for (va_abs in seq(min(depth$va_abs), max(depth$va_abs), length.out = 5)) {
#'   lines(1:1000, setup_movement_pr(1:1000, data.frame(va_abs = va_abs)), lwd = 0.25 + va_abs / 5)
#' }
#' ## Implement algorithm with adjusted movement model
#' out_8 <- pf(
#'   record = record,
#'   data = depth,
#'   origin = xy,
#'   calc_distance = "euclid",
#'   calc_movement_pr = setup_movement_pr,
#'   mobility = 250,
#'   n = 10L,
#'   seed = 1
#' )
#'
#' #### Dealing with convergence failures.
#' # The algorithm can fail to converge. There are a variety of options
#' # ... that can be trialled to deal with this:
#' # ... ... Updating the outputs from an earlier time step
#' # ... ... Increase the number of particles
#' # ... ... Resampling
#' # ... ... Reduce the number of temporal resolution
#' # ... ... Tweak movement model, depth error, mobility etc.
#'
#' #### Example (9): Update particle histories from an earlier time step
#' # ... via update_history and update_history_from_time_step
#' out_9 <- pf(
#'   record = record,
#'   data = depth,
#'   origin = xy,
#'   calc_distance = "euclid",
#'   calc_movement_pr = setup_movement_pr,
#'   mobility = 250,
#'   n = 10L,
#'   seed = 1,
#'   update_history = out_8$history,
#'   update_history_from_time_step = 5
#' )
#'
#' #### Example (10): Implement re-sampling via resample
#' # Here, for demonstration purposes, we implement the algorithm from
#' # ... scratch with the original movement model, this time re-sampling
#' # ... possible positions with equal probability if there are
#' # ... fewer than 10 unique positions.
#' # ... (Normally resample would be less than n.)
#' # In this example, the function re-samples candidate starting positions at t = 9
#' out_10 <- pf(
#'   record = record,
#'   data = depth,
#'   origin = xy,
#'   calc_distance = "euclid",
#'   calc_movement_pr = pf_setup_movement_pr,
#'   mobility = 250,
#'   n = 10L,
#'   resample = 10,
#'   seed = 1
#' )
#'
#' #### Example (10): A simulation workflow
#' # This example provides a simulation workflow for comparing simulated and
#' # ... reconstructed paths. Specifically, we simulate movement in an area,
#' # ... so that we know the 'true' path. We then implement the DCPF algorithm
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
#' dat_gebco_planar <- raster::raster(
#'   crs = raster::crs(dat_gebco),
#'   ext = raster::extent(dat_gebco),
#'   resolution = 25
#' )
#' dat_gebco_planar <- raster::resample(dat_gebco, dat_gebco_planar, method = "bilinear")
#' # Define 'sea' for movement simulation
#' dat_coast <- raster::crop(dat_coast, raster::extent(dat_gebco))
#' dat_sea <- invert_poly(dat_coast)
#' # Visualise area
#' prettyGraphics::pretty_map(dat_gebco_planar,
#'   add_rasters = list(x = dat_gebco_planar),
#'   add_polys = list(x = dat_sea)
#' )
#'
#' ## (C) Simulate movement path
#' # Define movement parameters
#' sim_steps <- function(...) stats::rgamma(1, shape = 15, scale = 5)
#' prettyGraphics::pretty_hist(stats::rgamma(10000, shape = 15, scale = 5), breaks = 100)
#' # Define animal's origin
#' origin_sim <- sp::spsample(dat_sea, n = 1, type = "random")
#' origin_sim <- sp::coordinates(origin_sim)
#' # Simulate path
#' path_sim <- sim_path_sa(
#'   n = 10,
#'   p_1 = origin_sim,
#'   area = dat_sea,
#'   sim_step = sim_steps,
#'   add_rasters = list(x = dat_gebco_planar),
#'   seed = seed
#' )
#' # Get resultant depth time series
#' path_sim <- path_sim$xy_mat
#' path_sim <- data.frame(
#'   path_id = 1,
#'   cell_id = raster::cellFromXY(dat_gebco_planar, path_sim),
#'   cell_x = path_sim[, 1],
#'   cell_y = path_sim[, 2],
#'   timestep = 1:nrow(path_sim)
#' )
#' path_sim$cell_z <- raster::extract(dat_gebco_planar, path_sim$cell_id)
#' prettyGraphics::pretty_plot(path_sim$cell_z, type = "l")
#' # Check simulated movements on the scale of the grid
#' sp::spDists(raster::xyFromCell(dat_gebco_planar, path_sim$cell_id),
#'   segments = TRUE
#' )
#' # Simulate 'observed' depth time series given some error
#' # ... For illustration, we will make the error smaller in this example
#' cde <- function(...) matrix(c(-2.5, 2.5), nrow = 2)
#' depth_obs <- runif(
#'   length(path_sim$cell_z),
#'   path_sim$cell_z + cde(path_sim$cell_z)[1],
#'   path_sim$cell_z + cde(path_sim$cell_z)[2]
#' )
#' depth_obs <- data.frame(depth = depth_obs)
#' # Compare 'true' and 'observed' depth time series
#' pf_plot_1d(path_sim, depth_obs,
#'   type = "b", cex = 0.5,
#'   add_lines = list(col = "royalblue", type = "l")
#' )
#'
#' ## Implement dc() on 'observed' time series
#' dc_out <- dc(
#'   archival = depth_obs,
#'   bathy = dat_gebco_planar,
#'   calc_depth_error = cde,
#'   save_record_spatial = NULL
#' )
#'
#' ## Implement pf() on the results from dc()
#' # We will assume that the origin was known.
#' # ... We will mostly use the default options.
#' history_dcpf <-
#'   pf(
#'     record = acdc_access_maps(acdc_simplify(dc_out, type = "dc"),
#'       type = "map_timestep"
#'     ),
#'     origin = origin_sim,
#'     calc_distance = "euclid",
#'     mobility = 200,
#'     n = 10L
#'   )
#'
#' ## Visualise particle histories in relation to simulated paths
#' # Here, each plot shows the particles sampled at a particular time step
#' # ... The green area shows areas of the requisite depth at that time step and the
#' # ... particles show sampled locations at that time step. The simulated path
#' # ... is shown in black.
#' pp <- graphics::par(mfrow = c(3, 4))
#' pf_plot_history(history_dcpf,
#'   add_particles = list(pch = 21),
#'   add_paths = list(x = path_sim$cell_x, path_sim$cell_y, length = 0.05),
#'   xlim = range(path_sim$cell_x), ylim = range(path_sim$cell_y),
#'   crop_spatial = TRUE,
#'   prompt = FALSE
#' )
#' graphics::par(pp)
#'
#' ## Assemble paths
#' path_dcpf <- pf_simplify(history_dcpf, bathy = dat_gebco_planar)
#'
#' ## Compare reconstructed versus observed depth time series
#' pf_plot_1d(path_dcpf, depth_obs)
#'
#' ## Show that the distances between sequential positions are within the restrictions
#' # ... of the movement model
#' require(rlang)
#' path_dcpf <-
#'   path_dcpf %>%
#'   dplyr::group_by(.data$path_id) %>%
#'   dplyr::mutate(
#'     cell_x2 = dplyr::lag(.data$cell_x),
#'     cell_y2 = dplyr::lag(.data$cell_y),
#'     dist_1 = sqrt((.data$cell_x2 - .data$cell_x)^2 +
#'       (.data$cell_y2 - .data$ cell_y)^2)
#'   )
#' path_dcpf
#' range(path_dcpf$dist_1, na.rm = TRUE)
#'
#' ## Visualise paths
#' # Zoom around path
#' xlim <- range(c(path_sim$cell_x, path_dcpf$cell_x), na.rm = TRUE)
#' ylim <- range(c(path_sim$cell_y, path_dcpf$cell_y), na.rm = TRUE)
#' boundaries <- raster::extent(xlim, ylim)
#' area <- raster::crop(dat_gebco_planar, boundaries)
#' # Define function to add simulated path for comparison
#' add_paths_sim <-
#'   function() {
#'     prettyGraphics::add_sp_path(path_sim$cell_x, path_sim$cell_y,
#'       lwd = 2, length = 0.01
#'     )
#'   }
#' # Make plots
#' if (interactive()) {
#'   pf_plot_2d(path_dcpf, area,
#'     add_paths = list(length = 0.05),
#'     add_additional = add_paths_sim,
#'     prompt = TRUE
#'   )
#' }
#'
#' #### Example (11): Write a dataframe of sampled particles to file at each time step
#' # Define directory in which to save files
#' root <- paste0(tempdir(), "/pf/")
#' dir.create(root)
#' # Implement PF, writing the history of particle samples at each time step
#' out_11 <- pf(
#'   record = record,
#'   origin = xy,
#'   write_history = list(file = root)
#' )
#' # Read the record into R
#' history_from_file <- lapply(pf_access_history_files(root), readRDS)
#' # Show that the history recorded in 'out_11' is the same as that recorded by the saved files
#' length(out_11$history)
#' length(history_from_file)
#' identical(out_11$history, history_from_file)
#' cbind(
#'   out_11$history[[1]],
#'   history_from_file[[1]]
#' )
#'
#' @return The function returns a \code{\link[flapper]{pf_archive-class}} object. This is a named list that includes the parameters used to generate function outputs (`args') and the particles sampled at each time step (`history'). The latter can be assembled into a dataframe of movement paths via \code{\link[flapper]{pf_simplify}}.
#'
#' @seealso This routine builds on the AC (\code{\link[flapper]{ac}}), (\code{\link[flapper]{dc}}) and (\code{\link[flapper]{acdc}}) algorithms. For the movement model, Euclidean distances are obtained from \code{\link[raster]{distanceFromPoints}} or shortest distances are obtained from \code{\link[flapper]{lcp_from_point}} (via \code{\link[flapper]{lcp_costs}} and \code{\link[flapper]{lcp_graph_surface}}, unless \code{calc_distance_graph} is supplied). The default movement model applied to these distances is \code{\link[flapper]{pf_setup_movement_pr}}. \code{\link[flapper]{get_mvt_resting}} provides a means to assign resting/non-resting behaviour to time steps (in \code{data}) for behaviourally dependent movement models. Particle histories can be visualised with \code{\link[flapper]{pf_plot_history}} and joined into paths via \code{\link[flapper]{pf_simplify}}. For processed paths, shortest distances/paths between sequential locations can be interpolated via \code{\link[flapper]{lcp_interp}}, which is a wrapper for the \code{\link[flapper]{lcp_over_surface}} routine. This can be useful for checking whether the faster Euclidean distances method is acceptable and, if so, for post-hoc adjustments of movement probabilities based on shortest distances (see \code{\link[flapper]{lcp_interp}}). Paths can be visualised with \code{\link[flapper]{pf_plot_1d}}, \code{\link[flapper]{pf_plot_2d}} and \code{\link[flapper]{pf_plot_3d}}. The log-probability of the paths, given the movement model, can be calculated via \code{\link[flapper]{pf_loglik}}.
#'
#' @author Edward Lavender
#' @export

pf <- function(record,
               data = NULL,
               origin = NULL,
               calc_distance = c("euclid", "lcp"),
               calc_distance_euclid_fast = TRUE,
               bathy = NULL,
               calc_distance_graph = NULL,
               calc_movement_pr = pf_setup_movement_pr,
               calc_movement_pr_from_origin = calc_movement_pr,
               mobility = NULL,
               mobility_from_origin = mobility,
               n = 10L,
               resample = NULL,
               update_history = NULL,
               update_history_from_time_step = 1L,
               write_history = NULL,
               cl = NULL, use_all_cores = FALSE,
               seed = NULL,
               verbose = TRUE, con = "",
               optimisers = pf_setup_optimisers()) {
  #### Set up function
  t_onset <- Sys.time()
  if (con != "") {
    if (!verbose) {
      message("Input to 'con' ignored since verbose = FALSE.")
    } else {
      # Check directory
      check_dir(input = dirname(con))
      # Write black file to directory if required
      if (!file.exists(con)) {
        message(paste0(con, " does not exist: attempting to write file in specified directory..."))
        file.create(file1 = con)
        message("... Blank file successfully written to file.")
      }
    }
  }
  ## Define function
  append_messages <- ifelse(con == "", FALSE, TRUE)
  cat_to_cf <- function(..., message = verbose, file = con, append = append_messages) {
    if (message) cat(paste(..., "\n"), file = con, append = append)
  }
  cat_to_cf(paste0("flapper::pf() called (@ ", t_onset, ")..."))
  cat_to_cf("... Setting up function...")

  #### Define storage container for outputs
  if (!is.null(seed)) set.seed(seed)
  out_pf <- list(
    history = NULL,
    method = "pf",
    args = list(
      record = record,
      data = data,
      origin = origin,
      calc_distance = calc_distance,
      calc_distance_euclid_fast = calc_distance_euclid_fast,
      bathy = bathy,
      calc_distance_graph = calc_distance_graph,
      calc_movement_pr = calc_movement_pr,
      calc_movement_pr_from_origin = calc_movement_pr_from_origin,
      mobility = mobility,
      mobility_from_origin = mobility_from_origin,
      n = n,
      resample = resample,
      update_history = update_history,
      update_history_from_time_step = update_history_from_time_step,
      write_history = write_history,
      cl = cl,
      use_all_cores = use_all_cores,
      seed = seed,
      verbose = verbose,
      optimisers = optimisers
    )
  )

  #### Checks

  ## Check optimisers
  check_class(input = optimisers, to_class = "pf_optimiser", type = "stop")

  ## Check data and define data_1 and data_t_next if unsupplied
  if (!is.null(data)) {
    if (!inherits(data, "data.frame")) stop("'data' must be a data.frame")
    if (length(record) != nrow(data)) stop("The number of records (`length(record)`) does not equal the number of observations in 'data' (`nrow(data)`).")
    data <- data.table::as.data.table(data)
  } else {
    data_1 <- data_t_next <- NULL
  }

  ## Define blank list to store path histories
  if (is.null(update_history)) history <- list() else history <- update_history
  ## Check record and define record raster param
  record_1 <- record[[1]]
  read_records <- FALSE
  if (inherits(record_1, "character")) {
    read_records <- TRUE
    if (!file.exists(record_1)) stop(paste0("record[[1]] ('", record_1, "') does not exist."))
    record_1 <- raster::raster(record_1)
  }
  n_cell <- raster::ncell(record_1)
  all_cells <- 1:n_cell
  proj <- raster::crs(record_1)
  boundaries <- raster::extent(record_1)
  record_1_sbt <- record_1
  ## Check origin
  if (!is.null(origin)) if (!inherits(origin, "matrix")) stop("'origin' coordinates must be supplied as a matrix.")

  ## Check re-sample
  if (!is.null(resample)) {
    if (resample > n) stop("'resample' must be <= 'n'.")
  }

  ## Check and set up distance calculations and associated cluster options
  calc_distance <- match.arg(calc_distance)
  out_pf$args$calc_distance <- calc_distance
  if (calc_distance == "lcp") {
    if (is.null(bathy)) stop("'bathy' must be supplied if calc_distance = 'lcp'. ")
    if (!all.equal(raster::res(bathy)[1], raster::res(bathy)[2])) {
      stop("raster::res(bathy)[1] must equal raster::res(bathy)[2] for this option.")
    }
    raster_comparison <- tryCatch(raster::compareRaster(record_1, bathy), error = function(e) {
      return(e)
    })
    if (inherits(raster_comparison, "error")) {
      warning("record[[1]] and 'bathy' have different properties",
        immediate. = TRUE, call. = FALSE
      )
      stop(raster_comparison)
    }
    if (is.null(calc_distance_graph)) {
      cat_to_cf("... Setting up cost-surface for calc_distance = 'lcp'...")
      costs <- lcp_costs(bathy, verbose = verbose)
      cost <- costs$dist_total
      calc_distance_graph <- lcp_graph_surface(surface = bathy, cost = cost, verbose = verbose)
      out_pf$args$calc_distance_graph <- calc_distance_graph
    }
    if (is.null(mobility_from_origin)) bathy_sbt_1 <- bathy
    if (is.null(mobility)) bathy_sbt <- bathy
  } else {
    if (optimisers$use_calc_distance_euclid_backend_grass) {
      if (!requireNamespace("fasterRaster", quietly = TRUE)) {
        stop("'The 'fasterRaster' package is required if optimisers$use_calc_distance_euclid_backend_grass = TRUE.")
      }
      if (is.null(optimisers$use_grass_dir)) stop("For 'optimisers$use_calc_distance_euclid_backend_grass' distances method, 'optimisers$use_grass_dir' must be specified.")
    }
    if (use_all_cores) {
      warning("use_all_cores = TRUE ignored for calc_distance = 'euclid'.",
        call. = FALSE, immediate. = TRUE
      )
      use_all_cores <- FALSE
    }
    if (calc_distance_euclid_fast) {
      if (!is.null(cl)) {
        warning("'cl' argument ignored for calc_distance_euclid_fast = TRUE.",
          call. = FALSE, immediate. = TRUE
        )
        cl <- NULL
      }
    }
    if (!is.null(calc_distance_graph)) {
      warning("'calc_distance_graph' ignored for calc_distance = 'euclid'.",
        call. = FALSE, immediate = TRUE
      )
      calc_distance_graph <- NULL
    }
  }
  # Other cluster options
  if (!is.null(cl) & use_all_cores) {
    warning("'cl' and 'use_all_cores' cannot both be specified: setting cl = NULL.",
      call. = FALSE, immediate. = TRUE
    )
    cl <- NULL
  }
  .cl <- cl
  # Write history to file
  if (!is.null(write_history)) {
    check_named_list(input = write_history)
    check_names(input = write_history, req = "file")
    write_history$file <- check_dir(input = write_history$file, check_slash = TRUE)
    write_history_dir <- write_history$file
  }
  # Global variables for data.table
  id_current <- dist <- pr_current <- pr_current_adj <- NULL


  #### For first time step define set of cells that the individual could have occupied
  if (is.null(update_history)) {
    cat_to_cf("... Determining the set of possible starting locations (t = 1)...")
    if (optimisers$use_raster_operations) {
      cells_at_time_current <- raster::Which(x = record_1 > 0, cells = TRUE, na.rm = TRUE)
      cells_at_time_current <-
        data.table::data.table(
          id_current = cells_at_time_current,
          pr_current = raster::extract(record_1, cells_at_time_current)
        )
    } else {
      cells_at_time_current <-
        data.table::data.table(
          id_current = all_cells,
          pr_current = as.vector(record_1)
        )
    }

    # Adjust cell probabilities by distance from origin, if applicable, using mobility model
    if (!is.null(origin)) {
      # Re-define origin on record_1 grid for consistency
      # ... This is necessary so that Euclidean and LCP distances are calculated in the same way
      origin_cell_id <- raster::cellFromXY(record_1, origin)
      origin <- raster::xyFromCell(record_1, origin_cell_id)
      # Crop raster to focus on area within mobility_from_origin for speed
      if (!is.null(mobility_from_origin)) {
        origin_sp <- sp::SpatialPoints(origin, proj4string = proj)
        buf <- rgeos::gBuffer(origin_sp, width = mobility_from_origin)
      }
      # Get distances
      if (calc_distance == "euclid") {
        if (!is.null(mobility_from_origin)) {
          record_1_sbt <- raster::crop(record_1, buf)
          record_1_sbt <- raster::mask(record_1, buf)
        }
        if (!optimisers$use_calc_distance_euclid_backend_grass) {
          dist_1 <- raster::distanceFromPoints(record_1_sbt, origin)
        } else {
          dist_1 <- fasterRaster::fasterVectToRastDistance(
            rast = record_1_sbt,
            vect = sp::SpatialPointsDataFrame(
              coords = origin,
              data = data.frame(id = 1:nrow(origin)),
              proj4string = proj
            ),
            grassDir = optimisers$use_grass_dir
          )
        }
        # if(!is.null(mobility_from_origin)) dist_1 <- raster::extend(dist_1, boundaries, NA)
      } else if (calc_distance == "lcp") {
        if (!is.null(mobility_from_origin)) bathy_sbt_1 <- raster::mask(bathy, buf)
        dist_1 <- lcp_from_point(
          origin = origin,
          destination = raster::Which(x = record_1 > 0, cells = TRUE, na.rm = TRUE),
          surface = bathy_sbt_1,
          graph = calc_distance_graph,
          use_all_cores = use_all_cores,
          verbose = verbose
        )
      }
      # Convert distances to movement probabilities and combine with location probabilities
      if (!is.null(data)) data_1 <- data[1, , drop = FALSE]
      if (optimisers$use_raster_operations) {
        pr_1 <- raster::calc(dist_1, function(x) calc_movement_pr_from_origin(x, data_1))
        pr_1 <- pr_1 * record_1_sbt
        if (calc_distance == "euclid" & !is.null(mobility_from_origin)) pr_1 <- raster::extend(pr_1, boundaries, 0)
        cells_at_time_current[, pr_current := raster::extract(pr_1, cells_at_time_current$id)]
      } else {
        if (calc_distance == "euclid" & !is.null(mobility_from_origin)) dist_1 <- raster::extend(dist_1, boundaries, 0)
        cells_at_time_current[, dist := as.vector(dist_1)]
        cells_at_time_current[, pr_current := calc_movement_pr_from_origin(pr_current, data_1)]
        cells_at_time_current[, pr_current := pr_current * as.vector(record_1)]
      }
    }
    cells_at_time_current[is.na(pr_current), pr_current := 0]
  } else {
    cat_to_cf("... Using update_history...")
    cat_to_cf(paste0(
      "... Determining the set of possible starting locations for update (t = ",
      update_history_from_time_step, ")..."
    ))
    cells_at_time_current <- history[[update_history_from_time_step]]
  }


  #### Implement algorithm iteratively
  cat_to_cf("... Implementing algorithm iteratively over time steps...")
  for (t in update_history_from_time_step:length(record)) {
    #### For the current time step, select likely cells
    ## Check there are available cells at the current time step; if not
    # ... movement model is insufficient
    # ... algorithm 'unlucky' in that all n paths up to this point are 'dead ends'
    cat_to_cf(paste0("... ... Time = ", t, "..."))
    if (all(cells_at_time_current$pr_current == 0)) {
      message("The probability of all cells at time ", t, " is 0. Either (1) the algorithm has been 'unlucky' and all n = ", n, " particles up to this point have led to 'dead ends', (2) the mobility model ('mobility') is too limiting and/or (3) other model assumptions have been violated. The function will now stop, returning outputs up until this point.")
      out_pf$history <- history
      class(out_pf) <- c(class(out_pf), "pf_archive")
      return(out_pf)
    }
    ## Select cells
    # Initial selection
    cat_to_cf("... ... ... Selecting candidate starting positions for the current time step...")
    cells_at_time_current_sbt <- cells_at_time_current[sample(
      x = 1:length(cells_at_time_current$id_current),
      size = n,
      prob = cells_at_time_current$pr_current,
      replace = TRUE
    ), ]
    # Re-implement selection, equally across all cells with non 0 Pr, if fewer than a threshold number of cells have been selected
    if (!is.null(resample)) {
      if (length(unique(cells_at_time_current_sbt$id_current)) < resample) {
        cat_to_cf("... ... ... ... Resampling candidate starting positions for the current time step...")
        cells_at_time_current[, pr_current_adj := pr_current]
        cells_at_time_current[pr_current_adj != 0, pr_current_adj := 1]
        cells_at_time_current_sbt <- cells_at_time_current[sample(
          x = 1:length(cells_at_time_current$id_current),
          size = n,
          prob = cells_at_time_current$pr_current_adj,
          replace = TRUE
        ), ]
        cells_at_time_current[, pr_current_adj := NULL]
      }
    }
    cells_at_time_current_sbt$timestep <- t
    rownames(cells_at_time_current_sbt) <- NULL
    history[[t]] <- as.data.frame(cells_at_time_current_sbt)
    # Write history to file (if specified)
    if (!is.null(write_history)) {
      write_history$object <- history[[t]]
      write_history$file <- paste0(write_history_dir, "pf_", t, ".rds")
      do.call(saveRDS, write_history)
    }

    #### For each particle, identify possible next locations
    if (t <= (length(record) - 1)) {
      ## Set up
      cat_to_cf("... ... ... For each particle, getting the possible positions for the next time step...")

      # Get record for next time step
      record_at_time_next <- record[[t + 1]]
      if (read_records) record_at_time_next <- raster::raster(record_at_time_next)
      record_sbt <- record_at_time_next
      if (!is.null(data)) data_t_next <- data[t + 1, , drop = FALSE]

      ## Fast euclidean distances method
      if (calc_distance == "euclid" & calc_distance_euclid_fast) {
        # Get cell IDs and coordinates for the 'current' time step
        # ... (using record_at_time_next and cells_at_time_current_sbt)
        cell_all_xy <- raster::xyFromCell(record_at_time_next, cells_at_time_current_sbt$id_current)
        if (!is.null(mobility)) {
          cell_all_sp <- sp::SpatialPoints(cell_all_xy, proj4string = proj)
          buf <- rgeos::gBuffer(cell_all_sp, width = mobility)
          record_sbt <- raster::crop(record_at_time_next, buf)
          record_sbt <- raster::mask(record_sbt, buf)
        }
        # Get distances around cells
        if (!optimisers$use_calc_distance_euclid_backend_grass) {
          dist_all <- raster::distanceFromPoints(record_sbt, cell_all_xy)
        } else {
          dist_all <- fasterRaster::fasterVectToRastDistance(
            rast = record_sbt,
            vect = sp::SpatialPointsDataFrame(
              coords = cell_all_xy,
              data = data.frame(id = 1:nrow(cell_all_xy)),
              proj4string = proj
            ),
            grassDir = optimisers$use_grass_dir
          )
        }
        # if(!is.null(mobility)) dist_all <- raster::extend(dist_all, boundaries, NA)
        # Get probabilities based on distance and combine with layer of possible locations for next time step
        if (optimisers$use_raster_operations) {
          pr_all <- raster::calc(dist_all, function(x) calc_movement_pr(x, data_t_next))
          pr_all <- pr_all * record_sbt
          if (calc_distance == "euclid" & !is.null(mobility)) pr_all <- raster::extend(pr_all, boundaries, 0)
          pr_all <- data.table::data.table(
            id_previous = NA,
            pr_previous = NA,
            id_current = all_cells,
            pr_current = as.vector(pr_all)
          )
        } else {
          if (calc_distance == "euclid" & !is.null(mobility)) dist_all <- raster::extend(dist_all, boundaries, 0)
          pr_all <- data.table::data.table(
            id_previous = NA,
            pr_previous = NA,
            id_current = all_cells,
            dist = as.vector(dist_all)
          )
          pr_all[, pr_current := calc_movement_pr(dist, data_t_next)]
          pr_all[, pr_current := pr_current * as.vector(record_at_time_next)]
          pr_all[, dist := NULL]
        }
        cells_from_current_to_next <- pr_all[which(pr_all$pr_current > 0), ]


        ## Other distance methods (point-by-point)
      } else {
        cells_from_current_to_next <- pbapply::pblapply(1:n, cl = .cl, function(j) {
          # Define location
          # j <- 1
          cell_j <- cells_at_time_current_sbt[j, ]
          cell_j_xy <- raster::xyFromCell(record_at_time_next, cell_j$id_current)
          # Subset record_1 around location for speed
          if (!is.null(mobility)) {
            cell_j_sp <- sp::SpatialPoints(cell_j_xy, proj4string = proj)
            buf <- rgeos::gBuffer(cell_j_sp, width = mobility)
          }
          if (calc_distance == "euclid") {
            if (!is.null(mobility)) {
              record_sbt <- raster::crop(record_at_time_next, buf)
              record_sbt <- raster::mask(record_sbt, buf)
            }
            if (!optimisers$use_calc_distance_euclid_backend_grass) {
              dist_j <- raster::distanceFromPoints(record_sbt, cell_j_xy)
            } else {
              dist_j <- fasterRaster::fasterVectToRastDistance(
                rast = record_sbt,
                vect = sp::SpatialPointsDataFrame(
                  coords = cell_j_xy,
                  data = data.frame(id = 1L),
                  proj4string = proj
                ),
                grassDir = optimisers$use_grass_dir
              )
            }
          } else if (calc_distance == "lcp") {
            if (!is.null(mobility)) bathy_sbt <- raster::mask(bathy, buf)
            dist_j <- lcp_from_point(
              origin = cell_j_xy,
              destination = raster::Which(record_at_time_next > 0, cells = TRUE, na.rm = TRUE),
              surface = bathy_sbt,
              graph = calc_distance_graph,
              use_all_cores = use_all_cores,
              verbose = FALSE
            )
          }
          pr_j <- raster::calc(dist_j, function(x) calc_movement_pr(x, data_t_next))
          pr_j <- pr_j * record_sbt
          if (calc_distance == "euclid" & !is.null(mobility)) pr_j <- raster::extend(pr_j, boundaries, 0)
          pr_j <- data.table::data.table(
            id_current = cell_j$id_current,
            pr_current = cell_j$pr_current,
            id_next = 1:n_cell,
            pr_next = as.vector(pr_j)
          )
          pr_j <- pr_j[which(pr_j$pr_next > 0), ]
          if (nrow(pr_j) == 0) pr_j <- NULL
          return(pr_j)
        })
        cells_from_current_to_next <- compact(cells_from_current_to_next)
        cells_from_current_to_next <- do.call(rbind, cells_from_current_to_next)
        if (is.null(cells_from_current_to_next)) {
          cells_from_current_to_next <- data.table::data.table(
            id_current = integer(),
            pr_current = numeric(),
            id_next = integer(),
            pr_next = numeric()
          )
        }
      }

      ## Processing
      if (nrow(cells_from_current_to_next) == 0L) {
        message("Unable to sample any cells at time ", t, " for the next location for any of the ", n, " particles. Either (1) the algorithm has been 'unlucky' and all n = ", n, " particles up to this point have led to 'dead ends', (2) the mobility model ('mobility') is too limiting and/or (3) the depth error parameter from 'calc_depth_error' is too small. The function will now stop, returning outputs up until this point.")
        out_pf$history <- history
        class(out_pf) <- c(class(out_pf), "pf_archive")
        return(out_pf)
      }
      cells_at_time_current <- cells_from_current_to_next
      colnames(cells_at_time_current) <- c("id_previous", "pr_previous", "id_current", "pr_current")
    }
  }
  cl_stop(.cl)

  #### Define outputs
  out_pf$history <- history
  class(out_pf) <- c(class(out_pf), "pf_archive")

  #### Return outputs
  if (!is.null(seed)) set.seed(NULL)
  t_end <- Sys.time()
  total_duration <- difftime(t_end, t_onset, units = "mins")
  cat_to_cf(paste0("... flapper::pf() call completed (@ ", t_end, ") after ~", round(total_duration, digits = 2), " minutes."))
  return(out_pf)
}
