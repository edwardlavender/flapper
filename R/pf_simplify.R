########################################
########################################
#### pf_simplify()

#' @title Convert particle histories from \code{\link[flapper]{pf}} into movement paths
#' @description This function is designed to simplify the \code{\link[flapper]{pf_archive-class}} object from \code{\link[flapper]{pf}} that defines sampled particle histories into a set of movement paths. The function identifies pairs of cells between which movement may have occurred at each time step (if necessary), (re)calculates distances and probabilities between connected cell pairs and then, if specified, links pairwise movements between cells into a set of possible movement paths.
#' @param archive A \code{\link[flapper]{pf_archive-class}} object from \code{\link[flapper]{pf}}.
#' @param max_n_particles (optional) An integer that defines the maximum number of particles to selected at each time step. If supplied, particle samples are thinned, with \code{max_n_particles} retained at each time step (in a method specified by \code{max_n_particles_sampler}), prior to processing. This reduces computation time but may lead to failures in path building (if the subset of sample particles cannot be connected into any movement paths).
#' @param max_n_particles_sampler If \code{max_n_particles} is supplied, \code{max_n_particles_sampler} is a character that defines the sampling method. Currently supported options are: \code{"random"}, which implements random sampling; \code{"weighted"}, which implements weighted sampling, with random samples taken according to their probability at the current time step; and \code{"max"}, which selects the top \code{max_n_particles} most likely particles.
#' @param bathy A \code{\link[raster]{raster}} that defines the grid across the area within which the individual could have moved. \code{bathy} must be planar (i.e., Universal Transverse Mercator projection) with units of metres in x and y directions. If \code{calc_distance = "lcp"}, the surface's resolution is taken to define the distance between horizontally and vertically connected cells and must be the same in both x and y directions. Any cells with NA values (e.g., due to missing data) are treated as `impossible' to move though by the algorithm (see \code{\link[flapper]{lcp_over_surface}}). If unsupplied, \code{bathy} can be extracted from \code{archive}, if available.
#' @param calc_distance A character that defines the method used to calculate distances between sequential combinations of particles (see \code{\link[flapper]{pf}}). Currently supported options are Euclidean distances (\code{"euclid"}) or shortest (least-cost) distances ("lcp"). In practice, Euclidean distances are calculated and then, for the subset of connections that meet specified criteria (see \code{calc_distance_limit}, \code{calc_distance_barrier}, \code{mobility_from_origin} and \code{mobility}), Euclidean distances are updated to shortest distances, if specified. Note that \code{calc_distance} does not need to be the same method as used for \code{\link[flapper]{pf}}: it is often computationally beneficial to implement \code{\link[flapper]{pf}} using Euclidean distances and then, for the subset of sampled particles, implement \code{\link[flapper]{pf_simplify}} with \code{calc_distance = "lcp"} to re-compute distances using the shortest-distances algorithm, along with the adjusted probabilities. However, for large paths, the quickest option is to implement both functions using \code{calc_distance = "euclid"} and then interpolate shortest paths only for the set of returned paths (see \code{\link[flapper]{lcp_interp}}). If \code{calc_distance = NULL}, the method saved in \code{archive} is used.
#' @param calc_distance_lcp_fast (optional) If \code{calc_distance = "lcp"}, \code{calc_distance_lcp_fast} is a function that predicts the shortest distance. This function must accept two arguments (even if they are ignored): (a) a numeric vector of Euclidean distances and (b) an integer vector of barrier overlaps that defines whether (1) or not (1) each transect crosses a barrier (see \code{\link[flapper]{lcp_comp}}). This option avoids the need to implement shortest-distance algorithms.
#' @param calc_distance_graph (optional) If \code{calc_distance = "lcp"} and \code{calc_distance_lcp_fast = NULL}, \code{calc_distance_graph} is a graph object that defines the distances between connected cells in \code{bathy}. If unsupplied, this is taken from \code{archive$args$calc_distance_graph}, if available, or computed via \code{\link[flapper]{lcp_graph_surface}}.
#' @param calc_distance_limit (optional) If \code{calc_distance = "lcp"} and \code{calc_distance_lcp_fast = NULL}, \code{calc_distance_limit} is a number that defines the lower Euclidean distance limit for shortest-distances calculations. If supplied, shortest distances are only calculated for cell connections that are more than a Euclidean distance of \code{calc_distance_limit} apart. In other words, if supplied, it is assumed that there exists a valid shortest path (shorter than the maximum distance imposed by the animal's mobility constraints) if the Euclidean distance between two points is less than \code{calc_distance_limit} (unless that segment crosses a barrier: see \code{calc_distance_barrier}). This option can improve the speed of distance calculations. However, if supplied, note that Euclidean and shortest distances (and resultant probabilities) may be mixed in function outputs. This argument is currently only implemented for \code{\link[flapper]{pf_archive-class}} objects derived with \code{calc_distance = TRUE} and \code{calc_distance_euclid_fast = TRUE} via \code{\link[flapper]{pf}} for which connected cell pairs need to be derived (i.e., not following a previous implementation of \code{\link[flapper]{pf_simplify}}).
#' @param calc_distance_barrier (optional) If \code{calc_distance = "lcp"}, \code{calc_distance_barrier} is a simple feature geometry that defines a barrier, such as the coastline, to movement (see \code{\link[flapper]{segments_cross_barrier}}). The coordinate reference system for this object must match \code{bathy}. If supplied, shortest distances are only calculated for segments that cross a barrier (or exceed \code{calc_distance_limit}). This option can improve the speed of distance calculations. However, if supplied, note that Euclidean and shortest distances (and resultant probabilities) may be mixed in function outputs. All \code{calc_distance_barrier*} arguments are currently only implemented for \code{\link[flapper]{pf_archive-class}} objects derived with \code{calc_distance = TRUE} and \code{calc_distance_euclid_fast = TRUE} via \code{\link[flapper]{pf}} for which connected cell pairs need to be derived (i.e., not following a previous implementation of \code{\link[flapper]{pf_simplify}}).
#' @param calc_distance_barrier_limit (optional) If \code{calc_distance_barrier} is supplied, \code{calc_distance_barrier_limit} is the lower Euclidean limit for determining barrier overlaps. If supplied, barrier overlaps are only determined for cell connections that are more than \code{calc_distance_barrier_limit} apart. This option can reduce the number of cell connections for which barrier overlaps need to be determined (and ultimately the speed of distance calculations).
#' @param calc_distance_barrier_grid (optional) If \code{calc_distance_barrier} is supplied, \code{calc_distance_barrier_grid} is a \code{\link[raster]{raster}} that defines the distance to the barrier (see \code{\link[flapper]{segments_cross_barrier}}). The coordinate reference system must match \code{bathy}. If supplied, \code{mobility_from_origin} and \code{mobility} are also required. \code{calc_distance_barrier_grid} is used to reduce the wall time of the barrier-overlap routine (see \code{\link[flapper]{segments_cross_barrier}}).
#' @param calc_distance_restrict (optional) If and \code{calc_distance_lcp_fast = NULL} and \code{calc_distance_limit} and/or \code{calc_distance_barrier} are supplied, \code{calc_distance_restrict} is a logical variable that defines whether (\code{TRUE}) or not (\code{FALSE}) to restrict further the particle pairs for which shortest distances are calculated. If \code{TRUE}, the subset of particles flagged by \code{calc_distance_limit} and \code{calc_distance_barrier} for shortest-distance calculations is further restricted to consider only those pairs of particles for which there is not a valid connection to or from those particles under the Euclidean distance metric. Either way, there is no reduction in the set of particles, only in the number of connections evaluated between those particles. This can improve the speed of distance calculations.
#' @param calc_distance_algorithm,calc_distance_constant Additional shortest-distance calculation options if \code{calc_distance_lcp_fast = NULL} (see \code{\link[cppRouting]{get_distance_pair}}). \code{calc_distance_algorithm} is a character that defines the algorithm: \code{"Dijkstra"}, \code{"bi"}, \code{"A*"} or \code{"NBA"} are supported. \code{calc_distance_constant} is the numeric constant required to maintain the heuristic function admissible in the A* and NBA algorithms. For shortest distances (based on costs derived via \code{\link[flapper]{lcp_costs}}), the default (\code{calc_distance_constant = 1}) is appropriate.
#' @param mobility,mobility_from_origin (optional) The mobility parameters (see \code{\link[flapper]{pf}}). If unsupplied, these can be extracted from \code{archive}, if available. However, even if \code{\link[flapper]{pf}} was implemented without these options, it is beneficial to specify mobility limits here (especially if \code{calc_distance = "lcp"}) because they restrict the number of calculations that are required (for example, for shortest distances, at each time step, distances are only calculated for the subset of particle connections below \code{mobility_from_origin} or \code{mobility} in distance).
#' @param write_history A named list of arguments, passed to \code{\link[base]{saveRDS}}, to write `intermediate' files. The `file' argument should be the directory in which to write files. This argument is currently only implemented for \code{\link[flapper]{pf_archive-class}} objects derived with \code{calc_distance = TRUE} and \code{calc_distance_euclid_fast = TRUE} via \code{\link[flapper]{pf}} for which connected cell pairs need to be derived (i.e., not following a previous implementation of \code{\link[flapper]{pf_simplify}}). If supplied, two directories are created in `file' (1/ and 2/), in which dataframes of the pairwise distances between connected cells and the subset of those  that formed continuous paths from the start to the end of the time series are written, respectively. Files are named by time steps as `pf_1', `pf_2' and so on. Files for each time step are written and re-read from this directory during particle processing. This helps to minimise vector memory requirements because the information for all time steps does not have to be retained in memory at once.
#' @param cl,varlist,use_all_cores (optional) Parallelisation options for the first stage of the algorithm, which identifies connected cell pairs, associated distances and movement probabilities. The first parallelisation option is to parallelise the algorithm over time steps via \code{cl}. \code{cl} is (a) a cluster object from \code{\link[parallel]{makeCluster}} or (b) an integer that defines the number of child processes. If \code{cl} is supplied, \code{varlist} may be required. This is a character vector of variables for export (see \code{\link[flapper]{cl_export}}). Exported variables must be located in the global environment. If a cluster is supplied, the connection to the cluster is closed within the function (see \code{\link[flapper]{cl_stop}}). For further information, see \code{\link[flapper]{cl_lapply}} and \code{\link[flapper]{flapper-tips-parallel}}.The second parallelisation option is to parallelise shortest-distance calculations within time steps via a logical input (\code{TRUE}) to \code{use_all_cores} that is passed to \code{\link[cppRouting]{get_distance_matrix}}. This option is only implemented for \code{calc_distance = "lcp"}.
#' @param return A character (\code{return = "path"} or \code{return = "archive"}) that defines the type of object that is returned (see Details).
#' @param summarise_pr (optional) For \code{return = "archive"}, \code{summarise_pf} is a function or a logical input that defines whether or not (and how) to summarise the probabilities of duplicate cell records for each time step. If a function is supplied, only one record of each sampled cell is returned per time step, with the associated probability calculated from the probabilities of each sample of that cell using the supplied function. For example, \code{summarise_pr = max} returns the most probable sample. Alternatively, if a logical input (\code{summarise_pr = TRUE}) is supplied, only one record of each sampled cell is returned per time step, with the associated probability calculated as the sum of the normalised probabilities of all samples for that cell, rescaled so that the maximum probability takes a score of one. Specifying \code{summarise_pr} is useful for deriving maps of the `probability of use' across an area based on particle histories because it ensures that `probability of use' scores depend on the number of time steps during which an individual could have occupied a location, rather than the total number of samples of that location (see \code{\link[flapper]{pf_plot_map}}). Both \code{summarise_pr = NULL} and \code{summarise_pr = FALSE} suppress this argument.
#' @param max_n_copies (optional) For \code{return = "path"}, \code{max_n_copies} is an integer that specifies the maximum number of copies of a sampled cell that are retained at each time stamp. Each copy represents a different route to that cell. By default, all copies (i.e. routes to that cell are retained) via \code{max_n_copies = NULL}. However, in cases where there are a large number of paths through a landscape, the function can run into vector memory limitations during path assembly, so \code{max_n_copies} may need to be set. In this case, at each time step, if there are more than \code{max_n_copies} paths to a given cell, then a subset of these (\code{max_n_copies}) are sampled, according to the \code{max_n_copies_sampler} argument.
#' @param max_n_copies_sampler (optional) For \code{return = "path"}, if \code{max_n_copies} is supplied, \code{max_n_copies_sampler} is a character that defines the sampling method. Currently supported options are: \code{"random"}, which implements random sampling; \code{"weighted"}, which implements weighted sampling, with random samples taken according to their probability at the current time step; and \code{"max"}, which selects for the top \code{max_n_copies} most likely copies of a given cell according to the probability associated with movement into that cell from the previous location.
#' @param max_n_paths (optional) For \code{return = "path"}, \code{max_n_paths} is an integer that specifies the maximum number of paths to be reconstructed. During path assembly, following the implementation of \code{max_n_copies} (if provided), \code{max_n_paths} are selected at random at each time step. This option is provided to improve the speed of path assembly in situations with large numbers of paths. \code{max_n_paths = NULL} computes all paths (if possible).
#' @param add_origin For \code{return = "path"}, \code{add_origin} is a logical input that defines whether or not to include the origin in the returned dataframe.
#' @param verbose A logical input that defines whether or not to print messages to the console to monitor function progress.
#'
#' @details The implementation of this function depends on how \code{\link[flapper]{pf}} has been implemented and the \code{return} argument. Under the default options in \code{\link[flapper]{pf}}, the fast Euclidean distances method is used to sample sequential particle positions, in which case the history of each particle through the landscape is not retained and has to be assembled afterwards. In this case, \code{\link[flapper]{pf_simplify}} calculates the distances between all combinations of cells at each time step, using either a Euclidean distances or shortest distances algorithm according to the input to \code{calc_distance}. Distances are converted to probabilities using the `intrinsic' probabilities associated with each location and the movement models retained in \code{archive} from the call to \code{\link[flapper]{pf}} to identify possible movement paths between cells at each time step. If the fast Euclidean distances method has not been used, then pairwise cell movements are retained by \code{\link[flapper]{pf}}. In this case, the function simply recalculates distances between sequential cell pairs and the associated cell probabilities, which are then processed according to the \code{return} argument.
#'
#' Following the identification of pairwise cell movements, if \code{return = "archive"}, the function selects all of the unique cells at each time step that were connected to cells at the next time step. (For cells that were selected multiple times at a given time step, due to sampling with replacement in \code{\link[flapper]{pf}}, if \code{summarise_pr} is supplied, only one sample is retained: in maps of the `probability of use' across an area (see \code{\link[flapper]{pf_plot_map}}), this ensures that cell scores depend on the number of time steps when the individual could have occupied a given cell, rather than the total number of samples of a location.) Otherwise, if \code{return = "path"}, pairwise cell movements are assembled into complete movement paths.
#'
#' @return If \code{return = "archive"}, the function returns a \code{\link[flapper]{pf_archive-class}} object, as inputted, but in which only the most likely record of each cell that was connected to cells at the next time step is retained and with the \code{method = "pf_simplify"} flag. If \code{return = "path"}, the function returns a \code{\link[flapper]{pf_path-class}} object, which is a dataframe that defines the movement paths.
#'
#' @examples
#' #### Example particle histories
#' # In these examples, we will use the example particle histories included in flapper
#' summary(dat_dcpf_histories)
#'
#' #### Example (1): The default implementation
#' paths_1   <- pf_simplify(dat_dcpf_histories)
#'
#' ## Demonstration that the distance and probabilities calculations are correct
#' # The simple method below works if three conditions are met:
#' # ... The 'intrinsic' probability associated with each cell is the same (as for DC algorithm);
#' # ... Paths have been reconstructed via pf_simplify() using Euclidean distances;
#' # ... The calc_movement_pr() movement model applies to all time steps;
#' require(magrittr)
#' require(rlang)
#' paths_1 <-
#'   paths_1 %>%
#'   dplyr::group_by(.data$path_id) %>%
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
#' pf_plot_1d(paths_1, dat_dc$args$archival)
#'
#' ## Examine paths
#' # Log likelihood
#' pf_loglik(paths_1)
#' # 2-d visualisation
#' pf_plot_2d(paths_1, dat_dcpf_histories$args$bathy,
#'            add_paths = list(length = 0.05))
#' # 3-d visualisation
#' pf_plot_3d(paths_1, dat_dcpf_histories$args$bathy)
#'
#' #### Example (2): Re-calculate distances as shortest distances
#'
#' ## Implement flapper::pf()
#' # For this example, we need to increase the number of particles
#' # ... for Euclidean-based sampling to generate viable paths
#' # ... when we consider shortest distances
#' set.seed(1)
#' dcpf_args <- dat_dcpf_histories$args
#' dcpf_args$calc_distance_euclid_fast <- TRUE
#' dcpf_args$n <- 50
#' out_dcpf_2 <- do.call(pf, dcpf_args)
#'
#' ## Implement pf_simplify() using shortest distances
#' paths_2 <- pf_simplify(out_dcpf_2, calc_distance = "lcp")
#' # ... Duration: ~ 0.655 s
#' system.time(
#'   invisible(utils::capture.output(
#'     pf_simplify(out_dcpf_2, calc_distance = "lcp")
#'     ))
#'   )
#'
#' ## Demonstrate the LCP calculations are correct
#' paths_2_lcps <- lcp_interp(paths_2,
#'                             out_dcpf_2$args$bathy,
#'                             calc_distance = TRUE)
#' head(cbind(paths_2$dist, paths_2_lcps$dist_lcp$dist))
#'
#' ## Trial options for increasing speed of shortest-distance calculations
#' # Speed up shortest-distance calculations via (a) the graph:
#' # ... Duration: ~0.495 s
#' # ... Note that you may achieve further speed improvements via
#' # ... a simplified/contracted graph
#' # ... ... see cppRouting::cpp_simplify()
#' # ... ... see cppRouting::cpp_contract()
#' costs <- lcp_costs(out_dcpf_2$args$bathy)
#' graph <- lcp_graph_surface(out_dcpf_2$args$bathy, costs$dist_total)
#' system.time(
#'   invisible(utils::capture.output(
#'     pf_simplify(out_dcpf_2,
#'                 calc_distance = "lcp",
#'                 calc_distance_graph = graph)
#'   ))
#' )
#' # Speed up shortest-distance calculations via (b) the lower Euclid dist limit
#' # ... Duration: ~0.493 s
#' costs <- lcp_costs(out_dcpf_2$args$bathy)
#' graph <- lcp_graph_surface(out_dcpf_2$args$bathy, costs$dist_total)
#' system.time(
#'   invisible(utils::capture.output(
#'     pf_simplify(out_dcpf_2,
#'                 calc_distance = "lcp",
#'                 calc_distance_graph = graph,
#'                 calc_distance_limit = 100)
#'   ))
#' )
#' # Speed up shortest-distance calculations via (c) the barrier
#' # ... Duration: ~1.411 s (much slower in this example)
#' coastline <- sf::st_as_sf(dat_coast)
#' sf::st_crs(coastline) <- NA
#' system.time(
#'   invisible(utils::capture.output(
#'     pf_simplify(out_dcpf_2,
#'                 calc_distance = "lcp",
#'                 calc_distance_graph = graph,
#'                 calc_distance_limit = 100,
#'                 calc_distance_barrier = coastline)
#'   ))
#' )
#' # Speed up calculations via (d) mobility limits
#' # ... (In the examples above, the mobility parameters
#' # ... can be extracted from out_dcpf_2,
#' # ... so specifying them directly here in this example makes
#' # ... no material difference, but this is not necessarily the case
#' # ... if pf() has been implemented without mobility parameters).
#' system.time(
#'   invisible(utils::capture.output(
#'     pf_simplify(out_dcpf_2,
#'                 calc_distance = "lcp",
#'                 calc_distance_graph = graph,
#'                 calc_distance_limit = 100,
#'                 mobility = 200,
#'                 mobility_from_origin = 200)
#'   ))
#' )
#' # Speed up calculations via (e) parallelisation
#' # ... see the details in the documentation.
#'
#' #### Example (3): Restrict the number of routes to each cell at each time step
#' # Implement approach for different numbers of copies
#' # Since we only have sampled a small number of particles for this simulation
#' # ... this does not make any difference here, but it can dramatically reduce
#' # ... the time taken to assemble paths and prevent vector memory issues.
#' paths_3a <- pf_simplify(dat_dcpf_histories, max_n_copies = 1)
#' paths_3b <- pf_simplify(dat_dcpf_histories, max_n_copies = 5)
#' paths_3c <- pf_simplify(dat_dcpf_histories, max_n_copies = 7)
#' # Compare the number of paths retained
#' unique(paths_3a$path_id)
#' unique(paths_3b$path_id)
#' unique(paths_3c$path_id)
#'
#' #### Example (4): Change the sampling method used to retain paths
#' # Again, this doesn't make a difference here, but it can when there are
#' # ... more particles.
#' paths_4a <- pf_simplify(dat_dcpf_histories,
#'                           max_n_copies = 5,
#'                           max_n_copies_sampler = "random")
#' paths_4b <- pf_simplify(dat_dcpf_histories,
#'                           max_n_copies = 5,
#'                           max_n_copies_sampler = "weighted")
#' paths_4c <- pf_simplify(dat_dcpf_histories,
#'                           max_n_copies = 5,
#'                           max_n_copies_sampler = "max")
#' # Compare retained paths
#' pf_loglik(paths_3a)
#' pf_loglik(paths_3b)
#' pf_loglik(paths_3c)
#'
#' #### Example (5): Set the maximum number of paths for reconstruction (for speed)
#' # Reconstruct all paths (note you may experience vector memory limitations)
#' # unique(pf_simplify(dat_dcpf_histories, max_n_paths = NULL)$path_id)
#' # Reconstruct one path
#' unique(pf_simplify(dat_dcpf_histories, max_n_paths = 1)$path_id)
#' # Reconstruct (at most) five paths
#' unique(pf_simplify(dat_dcpf_histories, max_n_paths = 5)$path_id)
#'
#' #### Example (6): Retain/drop the origin, if specified
#' # For the example particle histories, an origin was specified
#' dat_dcpf_histories$args$origin
#' # This is included as 'timestep = 0' in the returned dataframe
#' # ... with the coordinates re-defined on bathy:
#' paths_5a <- pf_simplify(dat_dcpf_histories)
#' paths_5a[1, c("cell_x", "cell_y")]
#' raster::xyFromCell(dat_dcpf_histories$args$bathy,
#'                    raster::cellFromXY(dat_dcpf_histories$args$bathy,
#'                                       dat_dcpf_histories$args$origin))
#' head(paths_5a)
#' # If specified, the origin is dropped with add_origin = FALSE
#' paths_5b <- pf_simplify(dat_dcpf_histories, add_origin = FALSE)
#' head(paths_5b)
#'
#' #### Example (6) Get particle samples for connected particles
#' ## Implement DCPF with more particles for demonstration purposes
#' set.seed(1)
#' dcpf_args <- dat_dcpf_histories$args
#' dcpf_args$calc_distance_euclid_fast <- TRUE
#' dcpf_args$n <- 250
#' out_dcpf_6a <- do.call(pf, dcpf_args)
#' head(out_dcpf_6a$history[[1]])
#' ## Extract particle samples for connected particles
#' # There may be multiple records of any given cell at any given time step
#' # ... due to sampling with replacement.
#' out_dcpf_6b <- pf_simplify(out_dcpf_6a, return = "archive")
#' head(out_dcpf_6b$history[[1]])
#' table(duplicated(out_dcpf_6b$history[[1]]$id_current))
#' ## Extract particle samples for connected particles,
#' # ... with only the most likely record of each particle returned.
#' # We can implement the approach using out_dcpf_6b
#' # ... to skip distance calculations.
#' # Now, there is only one (the most likely) record of sampled cells
#' # ... at each time step.
#' out_dcpf_6c <- pf_simplify(out_dcpf_6b,
#'                            summarise_pr = TRUE,
#'                            return = "archive")
#' head(out_dcpf_6c$history[[1]])
#' table(duplicated(out_dcpf_6c$history[[1]]$id_current))
#' ## Make movement paths
#' # Again, we use out_dcpf_6b to skip distance calculations.
#' out_dcpf_6d <- pf_simplify(out_dcpf_6b,
#'                            max_n_copies = 2L,
#'                            return = "path")
#' ## Compare resultant maps
#' # The map for all particles is influenced by particles that were 'dead ends',
#' # ... which isn't ideal for a map of space use.
#' # The map for connected samples deals with this problem, but is influenced by
#' # ... the total number of samples of each cell, rather than the number of time steps
#' # ... in which the individual could have been located in a given cell.
#' # The map for unique, connected samples deals with this issue, so that scores
#' # ... represent the number of time steps in which the individual could have occupied
#' # ... a given cell, over the length of the time series.
#' # The map for the paths are sparser because paths have only been reconstructed
#' # ... for a sample of sampled particles.
#' pp <- par(mfrow = c(2, 2), oma = c(2, 2, 2, 2), mar = c(2, 4, 2, 4))
#' paa <- list(side = 1:4, axis = list(labels = FALSE))
#' transform = NULL
#' m_1 <- pf_plot_map(out_dcpf_6a, dcpf_args$bathy,
#'                    transform = transform,
#'                    pretty_axis_args = paa, main = "all samples")
#' m_2 <- pf_plot_map(out_dcpf_6b, dcpf_args$bathy,
#'                    transform = transform,
#'                    pretty_axis_args = paa, main = "connected samples")
#' m_3 <- pf_plot_map(out_dcpf_6c, dcpf_args$bathy,
#'                    transform = transform,
#'                    pretty_axis_args = paa, main = "unique, connected samples")
#' m_4 <- pf_plot_map(out_dcpf_6d, dcpf_args$bathy,
#'                    transform = transform,
#'                    pretty_axis_args = paa, main = "paths")
#' par(pp)
#' # Note that all locations in reconstructed paths are derived from PF samples
#' all(out_dcpf_6d$cell_id[out_dcpf_6d$timestep != 0] %in%
#'       do.call(rbind, out_dcpf_6c$history)$id_current)
#' # But the paths only contain a subset of sampled particles
#' h_6a <- lapply(out_dcpf_6a$history, function(elm) elm[, "id_current"])
#' table(unique(unlist(h_6a)) %in% out_dcpf_6d$cell_id[out_dcpf_6d$timestep != 0])
#' h_6b <- lapply(out_dcpf_6b$history, function(elm) elm[, "id_current"])
#' table(unique(unlist(h_6b)) %in% out_dcpf_6d$cell_id[out_dcpf_6d$timestep != 0])
#'
#' @author Edward Lavender
#' @export

pf_simplify <- function(archive,
                        max_n_particles = NULL,
                        max_n_particles_sampler = c("random", "weighted", "max"),
                        bathy = NULL,
                        calc_distance = NULL,
                        calc_distance_lcp_fast = NULL,
                        calc_distance_graph = NULL,
                        calc_distance_limit = NULL,
                        calc_distance_barrier = NULL,
                        calc_distance_barrier_limit = NULL,
                        calc_distance_barrier_grid = NULL,
                        calc_distance_restrict = FALSE,
                        calc_distance_algorithm = "bi", calc_distance_constant = 1,
                        mobility = NULL,
                        mobility_from_origin = mobility,
                        write_history = NULL,
                        cl = NULL, varlist = NULL, use_all_cores = FALSE,
                        return = c("path", "archive"),
                        summarise_pr = FALSE,
                        max_n_copies = NULL,
                        max_n_copies_sampler = c("random", "weighted", "max"),
                        max_n_paths = 100L,
                        add_origin = TRUE,
                        verbose = TRUE){


  ########################################
  #### Function set up

  #### Set up
  cat_to_console <- function(..., show = verbose){
    if(show) cat(paste(..., "\n"))
  }
  t_onset <- Sys.time()
  cat_to_console(paste0("flapper::pf_simplify() called (@ ", t_onset, ")..."))
  if(!inherits(archive, "pf_archive")) stop("'archive' must be a 'pf_archive-class' object.", call. = FALSE)
  history <- archive$history
  if(is.null(bathy)) bathy <- archive$args$bathy
  if(is.null(bathy)) stop("'bathy' must be supplied via 'bathy' or 'archive$args$bathy' for this function.", call. = FALSE)
  bathy_xy <- raster::coordinates(bathy)
  layers   <- archive$args$record
  layers_1 <- layers[[1]]
  read_layers <- FALSE
  if(inherits(layers_1, "character")){
    read_layers <- TRUE
    if(!file.exists(layers_1)) stop(paste0("archive$args$record[[1]] ('", layers_1, "') does not exist."), call. = FALSE)
    layers_1    <- raster::raster(layers_1)
  }
  raster_comparison <- tryCatch(raster::compareRaster(layers_1, bathy), error = function(e) return(e))
  if(inherits(raster_comparison, "error")){
    warning("archive$args$record[[1]] and 'bathy' have different properties.",
            immediate. = TRUE, call. = FALSE)
    stop(raster_comparison, call. = FALSE)
  }
  data     <- archive$args$data
  origin   <- archive$args$origin
  if(!is.null(origin)) {
    origin_cell_id <- raster::cellFromXY(bathy, origin)
    origin         <- bathy_xy[origin_cell_id, , drop = FALSE]
  }
  if(is.null(calc_distance)) calc_distance <- archive$args$calc_distance
  calc_distance <- match.arg(calc_distance, c("euclid", "lcp"))
  if(!is.null(calc_distance_lcp_fast) & !inherits(calc_distance_lcp_fast, "function")){
    stop("'calc_distance_lcp_fast' expects NULL or a function.", call. = FALSE)
  }
  if(calc_distance == "lcp"){
    if(!is.null(calc_distance_barrier)) {
      check_crs(bathy, calc_distance_barrier, calc_distance_barrier_grid)
    } else{
      if(!is.null(calc_distance_barrier_limit) | !is.null(calc_distance_barrier_grid)){
        calc_distance_barrier_limit <- calc_distance_barrier_grid <- NULL
        warning("'calc_distance_barrier' is NULL; other 'calc_distance_barrier_*' arguments are ignored.",
                immediate. = TRUE, call. = FALSE)
      }
    }
    if(!inherits(calc_distance_lcp_fast, "function")){
      if(!is.null(calc_distance_limit) & !is.null(calc_distance_barrier_limit)){
        if(calc_distance_barrier_limit >= calc_distance_limit){
          calc_distance_barrier_limit <- NULL
          warning("'calc_distance_barrier_limit' >= 'calc_distance_limit': 'calc_distance_barrier_limit' ignored.",
                  immediate. = TRUE, call. = FALSE)
        }
      }
      if(calc_distance_restrict & is.null(calc_distance_limit) & is.null(calc_distance_barrier)){
        warning("'calc_distance_restrict' is only implemented if 'calc_distance_limit' and/or 'calc_distance_barrier' are supplied.", immediate. = TRUE, call. = FALSE)
        calc_distance_restrict <- FALSE
      }
    } else {
      if(!all(is.null(c(calc_distance_graph, calc_distance_limit))) | calc_distance_restrict){
        calc_distance_graph <- calc_distance_limit <- NULL
        calc_distance_restrict <- FALSE
        warning("'calc_distance_graph', 'calc_distance_limit' and/or 'calc_distance_restrict' arguments are ignored with 'calc_distance_lcp_fast'.",
                immediate. = TRUE, call. = FALSE)
      }
    }
  } else {
    if(!all(is.null(c(calc_distance_lcp_fast,
                      calc_distance_graph,
                      calc_distance_limit,
                      calc_distance_barrier,
                      calc_distance_barrier_limit,
                      calc_distance_barrier_grid))) | calc_distance_restrict){
      calc_distance_lcp_fast <-
        calc_distance_graph <-
        calc_distance_limit <-
        calc_distance_barrier <-
        calc_distance_barrier_limit <-
        calc_distance_barrier_grid <- NULL
      calc_distance_restrict <- FALSE
      warning("All other calc_distance_* arguments are ignored for calc_distance = 'euclid'.",
              immediate. = TRUE, call. = FALSE)
    }
  }
  calc_movement_pr_from_origin <- archive$args$calc_movement_pr_from_origin
  calc_movement_pr             <- archive$args$calc_movement_pr
  if(is.null(mobility)) mobility <- archive$args$mobility
  if(is.null(mobility_from_origin)) mobility_from_origin <- archive$args$mobility_from_origin
  if(is.null(mobility) | is.null(mobility_from_origin))
    message("'mobility' and/or 'mobility_from_origin' taken as NULL.")
  if(!is.null(calc_distance_barrier_grid) & (is.null(mobility_from_origin) | is.null(mobility)))
    stop("'mobility_from_origin' and 'mobility' are required if 'calc_distance_barrier_grid' is specified.", call. = FALSE)
  if(!is.null(write_history)){
    check_named_list(input = write_history)
    check_names(input = write_history, req = "file")
    write_history$file  <- check_dir(input = write_history$file, check_slash = TRUE)
    write_history_dir   <- write_history$file
    write_history_dir_1 <- paste0(write_history_dir, "1/")
    write_history_dir_2 <- paste0(write_history_dir, "2/")
    dir.create(write_history_dir_1)
    dir.create(write_history_dir_2)
    # Stop if write_history_dir_1/write_history_dir_2 are not empty
    # ... because all files are read from these directories later
    # ... which can lead to lists of the wrong length/contents if they
    # ... are already populated (e.g., with testing files).
    if(length(list.files(write_history_dir_1)) != 0L)
      stop(paste0("'", write_history_dir_1, "' is not empty."), call. = FALSE)
    if(length(list.files(write_history_dir_2)) != 0L)
      stop(paste0("'", write_history_dir_2, "' is not empty."), call. = FALSE)
  }
  # Match samplers
  max_n_particles_sampler <- match.arg(max_n_particles_sampler)
  max_n_copies_sampler    <- match.arg(max_n_copies_sampler)
  # Check cluster
  if(!is.null(cl) && use_all_cores){
    warning("Both 'cl' and 'use_all_cores' supplied: 'use_all_cores' ignored.",
            immediate. = TRUE, call. = FALSE)
    use_all_cores <- FALSE
  }
  if(use_all_cores && calc_distance != "lcp") {
    warning("'use_all_cores' ignored: calc_distance != 'lcp'.",
            immediate. = TRUE, call. = FALSE)
    use_all_cores <- FALSE
  }
  # Match return
  return <- match.arg(return)


  ########################################
  #### Thin particle samples

  if(!is.null(max_n_particles)){
    history <- pbapply::pblapply(history, function(history_for_t){
      if(max_n_particles_sampler == "random"){
        history_for_t <- history_for_t %>% dplyr::slice_sample(n = max_n_particles)
      } else if(max_n_particles_sampler == "weighted"){
        history_for_t <- history_for_t %>% dplyr::slice_sample(n = max_n_particles, weight_by = .data$pr_current)
      } else if(max_n_particles_sampler == "max"){
        history_for_t <- history_for_t %>% dplyr::arrange(.data$pr_current) %>% dplyr::slice(1:max_n_particles)
      }
      return(history_for_t)
    })
  }


  ########################################
  #### Identify sequential, pairwise connections between cells

  # This is necessary following all implementations of pf().
  # But it may not be necessary if pf_simplify() has previously been called with return = 'archive',
  # ... in which case this step has already been done and we can skip distance calculations.

  if(!rlang::has_name(history[[1]], "dist_current")){

    ## Implement setup for LCP distance calculations, if necessary.
    cat_to_console(paste0("... Getting pairwise cell movements based on calc_distance = '", calc_distance, "'..."))
    if(calc_distance == "lcp" & !inherits(calc_distance_lcp_fast, "function")){
      cat_to_console("... Setting up LCP calculations...")
      if(!all.equal(raster::res(bathy)[1], raster::res(bathy)[2])){
        stop("raster::res(bathy)[1] must equal raster::res(bathy)[2] for calc_distance = 'lcp'.", call. = FALSE)
      }
      if(is.null(calc_distance_graph)) calc_distance_graph <- archive$args$calc_distance_graph
      if(is.null(calc_distance_graph)){
        cat_to_console("... ... Setting up cost-surface for calc_distance = 'lcp'...")
        costs <- lcp_costs(bathy, verbose = verbose)
        cost  <- costs$dist_total
        calc_distance_graph <- lcp_graph_surface(surface = bathy, cost = cost, verbose = verbose)
      }
    }

    #### Calculate distances and probabilities between cell pairs
    # ... (1) Method for pf outputs derived via calc_distance_euclid_fast
    cat_to_console("... ... Stepping through time steps to join coordinate pairs...")
    if(archive$args$calc_distance == "euclid" & archive$args$calc_distance_euclid_fast){

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
      layers <- append(raster::setValues(layers_1, 1), layers)

      #### Re-define history with cell pairs
      # Define cluster/chunks
      if(is.null(cl))
        chunks <- 1L
      else
        if(inherits(cl, "cluster")) chunks <- length(cl) else chunks <- cl
      t_vec <- 2:length(history)
      t_ind_by_chunk <- parallel::splitIndices(length(t_vec), chunks)
      cl_export(cl, varlist)
      history_by_chunk <- pbapply::pblapply(t_ind_by_chunk, cl = cl, function(t_ind_for_chunk){
        ## Define time step indices for chunk
        t_vec_for_chunk <- t_vec[t_ind_for_chunk]
        history_for_chunk <- lapply(t_vec_for_chunk, function(t){
          ## Get full suite of allowed cells from current time step
          # (This is combined with movement probabilities to generate adjusted location probabilities below.)
          layers_2 <- layers[[t]]
          if(read_layers) layers_2 <- raster::raster(layers_2)

          ## Get retained cells from the previous and current step and coordinates
          tms <- history[[t - 1]]$timestep[1] + 1
          z1_id_current_unq <- unique(history[[t - 1]]$id_current)
          z2_id_current_unq <- unique(history[[t]]$id_current)

          ## Calculate Euclidean distances between all combinations of cells to identify possible movement paths between cell pairs
          # For both the origin (if applicable) and subsequent locations, we will calculate distances from cells IDs
          # ... For the origin, this is may result in slightly less accurate distances,
          # ... but it ensures consistency across cells and methods (euclidean and lcps).
          # Following the calculation of Euclidean distances (which is fast), we will drop any pairwise combinations that
          # ... exceed mobility_from_origin or mobility. We will then proceed to calculate LCPs for the remaining
          # ... combinations (if necessary), because this is a slower step.

          if((t - 1) == 1 & is.null(origin)){
            tmp <- data.frame(timestep = tms,
                              id_previous = NA,
                              id_current = z2_id_current_unq,
                              dist_current = 0,
                              pr_current = 1
            )

          } else {
            # Calculate Euclidean distances (using internal raster function for speed)
            # ... This returns a matrix where the first arg forms the rows
            # ... and the second arg forms the columns.
            dist_btw_cells <- .planedist2(bathy_xy[z1_id_current_unq, , drop = FALSE],
                                          bathy_xy[z2_id_current_unq, , drop = FALSE])
            # rownames(dist_btw_cells) <- z1_id_current_unq
            # colnames(dist_btw_cells) <- z2_id_current_unq
            # Extract distances into dataframe
            # ... Note that the matrix is non symmetrical (so we need to extract all distances).
            tmp <- data.frame(timestep     = tms,
                              id_previous  = z1_id_current_unq[as.vector(row(dist_btw_cells))],
                              id_current   = z2_id_current_unq[as.vector(col(dist_btw_cells))],
                              dist_current = as.vector(dist_btw_cells)
            )
            # Adjust distances according to mobility_from_origin or mobility
            if((t - 1) == 1)  mob <- mobility_from_origin else mob <- mobility
            if(!is.null(mob)) tmp <- tmp %>% dplyr::filter(.data$dist_current <= mob)

            # Calculate LCPs (if specified)
            if(calc_distance == "lcp"){
              calc_distance_lcp <- TRUE

              # ... Identify segments that exceed the lower distance limit
              if(!is.null(calc_distance_limit)) {
                index_in_limit <- which(tmp$dist_current > calc_distance_limit)
              } else index_in_limit <- integer(0)

              # ... Identify segments that cross a barrier
              # ... ... This step can be slow for large numbers of paths,
              # ... ... so we will focus on only those paths that
              # ... ... (a) were not flagged by the check above (if applicable), and/or
              # ... ... (b) meet the barrier distance threshold (if applicable)
              index_in_barrier <- integer(0)
              if(!is.null(calc_distance_barrier)){
                if(length(index_in_limit) > 0L){
                  if(!is.null(calc_distance_barrier_limit)){
                    index_for_barrier <- which(tmp$dist_current > calc_distance_barrier_limit &
                                                 tmp$dist_current <= calc_distance_limit)

                  } else {
                    index_for_barrier <- seq_len(nrow(tmp))
                    index_for_barrier <- index_for_barrier[!(index_for_barrier %in% index_in_limit)]
                  }
                } else {
                  if(!is.null(calc_distance_barrier_limit)){
                    index_for_barrier <- which(tmp$dist_current > calc_distance_barrier_limit)
                  } else {
                    index_for_barrier <- seq_len(nrow(tmp))
                  }
                }
                if(length(index_for_barrier) > 0L){
                  index_in_barrier <-
                    which(segments_cross_barrier(start = bathy_xy[tmp$id_previous[index_for_barrier], , drop = FALSE],
                                                 end = bathy_xy[tmp$id_current[index_for_barrier], , drop = FALSE],
                                                 barrier = calc_distance_barrier,
                                                 distance = calc_distance_barrier_grid,
                                                 mobility = mob))
                  index_in_barrier <- index_for_barrier[index_in_barrier]
                }
              }
              # ... Define index of segments for LCP calculations
              index_in_cond <- c(index_in_limit, index_in_barrier)
              if(length(index_in_cond) > 0L) index_in_cond <- unique(index_in_cond)
              # ... Drop elements for which the previous and current cells are
              # ... ... already represented in 'tmp' with valid paths (that passed the checks above)
              if(calc_distance_restrict){
                if(length(index_in_cond) > 0L){
                  tmp$recalc                <- FALSE
                  tmp$recalc[index_in_cond] <- TRUE
                  index_out_cond            <- !tmp$recalc
                  # This indexing approach is faster than filtering and Rcpp pf_contains()
                  # Referring to all cells is also much faster than unique(tmp$id_previous[index_out_cond]))
                  tmp$recalc[index_in_cond][
                    (tmp$id_previous[index_in_cond] %in% tmp$id_previous[index_out_cond]) &
                      (tmp$id_current[index_in_cond] %in% tmp$id_current[index_out_cond])
                  ] <- FALSE
                  index_in_cond <- which(tmp$recalc)
                }
              }
              # ... Implement LCP calculations across selected or all segments as required
              if(inherits(calc_distance_lcp_fast, "function")){
                tmp$barrier                   <- 0L
                tmp$barrier[index_in_barrier] <- 1L
                tmp$dist_current <- calc_distance_lcp_fast(tmp$dist_current, tmp$barrier)
              } else {
                if(!is.null(calc_distance_limit) | !is.null(calc_distance_barrier) | calc_distance_restrict){
                  if(length(index_in_cond) == 0L) calc_distance_lcp <- FALSE
                }
                if(calc_distance_lcp){
                  if(length(index_in_cond) > 0L){
                    tmp$dist_current[index_in_cond] <-
                      cppRouting::get_distance_pair(Graph = calc_distance_graph,
                                                    from = tmp$id_previous[index_in_cond],
                                                    to = tmp$id_current[index_in_cond],
                                                    algorithm = calc_distance_algorithm,
                                                    constant = calc_distance_constant,
                                                    allcores = use_all_cores)
                  } else {
                    tmp$dist_current <-
                      cppRouting::get_distance_pair(Graph = calc_distance_graph,
                                                    from = tmp$id_previous,
                                                    to = tmp$id_current,
                                                    algorithm = calc_distance_algorithm,
                                                    constant = calc_distance_constant,
                                                    allcores = use_all_cores)
                  }
                }
              }


              # Repeat filtration based on mobility_from_origin or mobility (copied from above)
              if(!is.null(mob)) tmp <- tmp %>% dplyr::filter(.data$dist_current <= mob)

              # Check there are remaining cells
              # ... If Euclidean distances have been used for sampling
              # ... it is not guaranteed that any particles will meet mobility_from_origin/mobility
              # ... criteria under if Euclidean distances used for sampling and LCPs used here
              if((t - 1) == 1){
                if(!is.null(mobility_from_origin))
                  if(nrow(tmp) == 0)
                    stop(paste0("No possible pairwise connections at time = ", t,
                                " under shortest distances given 'mobility_from_origin'."),
                         call. = FALSE)
              } else {
                if(!is.null(mobility))
                  if(nrow(tmp) == 0)
                    stop(paste0("No possible pairwise connections at time = ", t,
                                " under shortest distances given 'mobility'."),
                         call. = FALSE)
              }
            }

            ## Calculate probabilities of movement FROM previous cell INTO current cell
            # ... I.e., probability associated with CURRENT cells
            # ... For the first time step, probabilities depend on
            # ... ... If an origin has been specified then (a) distance from origin and (b) intrinsic probabilities
            # ... ... If an origin has not been specified then just the (b) intrinsic probabilities
            # ... For later time steps, probabilities depend on
            # ... ... (a) distances (and any other parameters in archival)
            # ... ... (b) intrinsic probabilities associated with cells (e.g., due to detection pr)
            # (a) Get probabilities based on movement
            if((t - 1) == 1){
              if(!is.null(mobility_from_origin))
                tmp$pr_current <- calc_movement_pr_from_origin(tmp$dist_current, data[t-1, ])
              else tmp$pr_current <- 1L
            } else {
              tmp$pr_current   <- calc_movement_pr(tmp$dist_current, data[t-1, ])
            }

            # (b) Combine with 'intrinsic' probabilities in each cell
            tmp$pr_current_intrinsic <- raster::extract(layers_2, tmp$id_current)
            tmp$pr_current           <- tmp$pr_current * tmp$pr_current_intrinsic
            tmp$pr_current_intrinsic <- NULL
          }

          ## Filter impossible movements
          tmp <- tmp %>% dplyr::filter(.data$pr_current > 0)

          ## Write dataframe to file (minimise memory requirements) or return
          if(!is.null(write_history)){
            write_history$object <- tmp
            write_history$file   <- paste0(write_history_dir_1, "pf_", t, ".rds")
            do.call(saveRDS, write_history)
            return(NULL)
          } else return(tmp)
        })
        return(history_for_chunk)
      })
      if(return == "path") cl_stop(cl)
      history <- purrr::flatten(history_by_chunk)

      # ... (2) Method for pf outputs not derived via calc_distance_euclid_fast
    } else {

      #### write_history is not currently implemented for this option
      if(!is.null(write_history))
        warning("write_history is not currently implemented for pf_archive objects not derived via the 'fast Euclidean distances' method.",
                immediate. = TRUE, call. = FALSE)
      write_history <- NULL

      #### Add origin to history, if necessary
      if(!is.null(origin)){
        history[[1]]$timestep      <- 0
        history[[1]]$id_previous   <- origin_cell_id
        history[[1]]$pr_previous   <- 1
      }

      #### Update history with distances and probabilities of movement between connected cells
      cl_export(cl, varlist)
      history <- pbapply::pblapply(1:length(history), cl = cl, function(t){

        #### Get full suite of allowed cells from current time step
        layers_2 <- layers[[t]]
        if(read_layers) layers_2 <- raster::raster(layers_2)

        #### Get history for time t
        d <- history[[t]]
        d$timestep <- t

        #### Calculate distances between connected cells
        if(!rlang::has_name(d, "id_previous")){
          d$dist_current <- NA
          d$pr_current   <- 1
        } else {
          d[, c("id_previous_x", "id_previous_y")]  <- raster::xyFromCell(layers_2, d$id_previous)
          d[, c("id_current_x", "id_current_y")]    <- raster::xyFromCell(layers_2, d$id_current)
          if(calc_distance == "euclid"){
            d$dist_current <-
              raster::pointDistance(d[, c("id_previous_x", "id_previous_y")],
                                    d[, c("id_current_x", "id_current_y")], lonlat = FALSE)
          } else if(calc_distance == "lcp"){
            # Calculate LCPs
            d$dist_current <- cppRouting::get_distance_pair(Graph = calc_distance_graph,
                                                            from = d$id_previous,
                                                            to = d$id_current,
                                                            algorithm = calc_distance_algorithm,
                                                            constant = calc_distance_constant,
                                                            allcores = use_all_cores)
            # Filter movements that exceed mobility_from_origin/mobility
            if(t == 1) {
              if(!is.null(mobility_from_origin)) d <- d %>% dplyr::filter(.data$dist_current <= mobility_from_origin)
            } else {
              if(!is.null(mobility)) d <- d %>% dplyr::filter(.data$dist_current <= mobility)
            }
          }
        }

        #### Re-calculate probabilities associated with movement into cells
        ## Get movement probabilities from distances between connected cells
        if(t == 1){
          if(is.null(origin)){
            d$pr_current <- 1
          } else {
            d$pr_current <- calc_movement_pr_from_origin(d$dist_current, data[1, ])
          }
        } else {
          d$pr_current <- calc_movement_pr(d$dist_current, data[t, ])
        }
        d$pr_current_intrinsic <- raster::extract(layers_2, d$id_current)
        d$pr_current <- d$pr_current * d$pr_current_intrinsic
        d$pr_current_intrinsic <- NULL
        d <- d %>% dplyr::filter(.data$pr_current > 0)
        return(d)
      })
      if(return == "path") cl_stop(cl)
    }


    ########################################
    #### Select connected cells

    #### Select the subset of cells at each time step that are connected to past positions
    cat_to_console("... ... Identifying connected cells...")
    if(!is.null(write_history)) {
      history_files <- pf_access_history_files(write_history_dir_1)
      history_files <- history <- lapply(history_files, function(x) x)
    }
    len_history <- length(history)
    if(!is.null(write_history)) history[[len_history]] <- readRDS(history_files[[len_history]])
    if(len_history > 5) pb <- pbapply::startpb(min = 2, max = len_history - 2)
    for(t in len_history:2){
      # t = len_history
      if(!is.null(write_history)) history[[t - 1]] <- readRDS(history_files[[t - 1]])
      history[[t - 1]] <-
        history[[t - 1]] %>%
        dplyr::filter(.data$id_current %in% history[[t]]$id_previous)
      if(nrow(history[[t - 1]]) == 0) stop(paste0("There are no connected cells at time = ", t-1,
                                                  " (given cells at time = ", t, ")."),
                                           call. = FALSE)
      if(!is.null(write_history)) {
        write_history$object <- history[[t - 1]]
        write_history$file   <- paste0(write_history_dir_2, "pf_", t-1, ".rds")
        do.call(saveRDS, write_history)
        history[[t + 1]] <- NA
      }
      if(len_history > 5) pbapply::setpb(pb, len_history - t)
    }
    if(len_history > 5) pbapply::closepb(pb)

  } else {
    if(!is.null(write_history)){
      warning("'write_history' is not currently implemented for this option.",
              immediate. = TRUE, call. = FALSE)
      write_history <- NULL
    }
  }

  #### List updated history files for processing (for return = 'archive' or 'path')
  if(!is.null(write_history)) {
    history_files <- pf_access_history_files(write_history_dir_2)
    history_files <- history <- lapply(history_files, function(x) x)
  }


  ########################################
  #### Process selected particles (if specified)

  if(return == "archive"){
    cat_to_console("... ... Processing connected cells for return = 'archive'...")
    if(is.null(cl)) {
      chunks <- 1L
    } else {
      if(inherits(cl, "cluster")) chunks <- length(cl) else chunks <- cl
    }
    t_vec <- seq_len(length(history))
    t_ind_by_chunk <- parallel::splitIndices(length(t_vec), chunks)
    history_by_chunk <- pbapply::pblapply(t_ind_by_chunk, cl = cl, function(t_ind_for_chunk){
      t_vec_for_chunk <- t_vec[t_ind_for_chunk]
      history_for_chunk <-
        lapply(t_vec_for_chunk, function(t){
          history_for_t <- history[[t]]
          if(!is.null(write_history)) history_for_t <- readRDS(history_for_t)
          history_for_t <-
            history_for_t %>%
            dplyr::group_by(.data$id_current) %>%
            dplyr::arrange(.data$id_current, dplyr::desc(.data$pr_current))
          if(!is.null(summarise_pr)){
            if(inherits(summarise_pr, "function")){
              history_for_t <-
                history_for_t %>%
                dplyr::mutate(pr_current = summarise_pr(.data$pr_current)) %>%
                dplyr::slice(1L)
            } else if(inherits(summarise_pr, "logical")){
              if(summarise_pr){
                denom <- sum(history_for_t$pr_current)
                history_for_t <-
                  history_for_t %>%
                  dplyr::mutate(pr_current = .data$pr_current/denom) %>%
                  dplyr::mutate(pr_current = sum(.data$pr_current)) %>%
                  dplyr::ungroup() %>%
                  dplyr::group_by(.data$id_current) %>%
                  dplyr::slice(1L)
              }
            } else stop("Implementation of 'summarise_pr' unrecognised: only functions or logical inputs are allowed.",
                        call. = FALSE)
          }
          history_for_t <-
            history_for_t %>%
            dplyr::ungroup() %>%
            data.frame()
          history_for_t[, colnames(history_for_t)[colnames(history_for_t) %in% c("id_previous", "pr_previous",
                                                                                 "id_current", "pr_current",
                                                                                 "timestep",
                                                                                 "dist_current")]]
          return(history_for_t)
        })
      return(history_for_chunk)
    })
    cl_stop(cl)
    history <- purrr::flatten(history_by_chunk)
    archive_for_connected_cells <- list(history = history,
                                        method = "pf_simplify",
                                        args = archive$args
    )
    class(archive_for_connected_cells) <- c(class(archive_for_connected_cells), "pf_archive")
    out <- archive_for_connected_cells


    ########################################
    #### Assemble and format paths (if specified)

  } else if(return == "path"){

    #### Assemble paths
    cat_to_console("... Assembling paths...")
    path <- list()
    if(!is.null(write_history)) history[[1]] <- readRDS(history_files[[1]])
    path[[1]] <- history[[1]]
    path[[1]]$id_1 <- path[[1]]$id_current
    path[[1]]$pr_1 <- path[[1]]$pr_current
    path[[1]]$dist_1 <- path[[1]]$dist_current
    path[[1]] <- path[[1]][, c("id_1", "pr_1", "dist_1", "id_current")]
    pb <- pbapply::startpb(min = 2, max = length(history) - 1)
    for(t in 1:(length(history) - 1)){
      if(!is.null(write_history)) history[[t + 1]] <- readRDS(history_files[[t + 1]])
      history_for_pair <- dplyr::inner_join(path[[t]], history[[t + 1]], by = c("id_current" = "id_previous"))
      history_for_pair <- dplyr::distinct(history_for_pair)
      if(!is.null(max_n_copies)){
        history_for_pair <-
          history_for_pair %>%
          dplyr::group_by(.data$id_current.y) %>%
          dplyr::arrange(dplyr::desc(.data$pr_current))
        if(max_n_copies_sampler == "random"){
          history_for_pair <- history_for_pair %>% dplyr::slice_sample(n = max_n_copies)
        } else if(max_n_copies_sampler == "weighted"){
          history_for_pair <- history_for_pair %>% dplyr::slice_sample(n = max_n_copies, weight_by = .data$pr_current)
        } else if(max_n_copies_sampler == "max"){
          history_for_pair <- history_for_pair %>% dplyr::slice(1:max_n_copies)
        }
      }
      if(!is.null(max_n_paths)){
        if(nrow(history_for_pair) > max_n_paths){
          history_for_pair <- history_for_pair %>% dplyr::slice_sample(n = max_n_paths)
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
      if(!is.null(write_history)){
        if(t > 2){
          path[[t-2]]    <- NA
          history[[t-2]] <- NA
        }
      }
      pbapply::setpb(pb, t)
    }
    pbapply::closepb(pb)
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
    path_df[, c("cell_x", "cell_y")] <- raster::xyFromCell(layers_1, path_df$cell_id)
    if(!is.null(origin) & add_origin){
      path_df[path_df$timestep == 0, "cell_x"] <- origin[1]
      path_df[path_df$timestep == 0, "cell_y"] <- origin[2]
    }
    if(!is.null(bathy)) {
      path_df[, "cell_z"] <- raster::extract(bathy, path_df[, c("cell_x", "cell_y")])
    } else path_df[, "cell_z"] <- NA
    path_df <- path_df[, c("path_id", "timestep", "cell_id", "cell_x", "cell_y", "cell_z", "cell_pr", "dist")]
    class(path_df) <- c(class(path_df), "pf_path")
    out <- path_df
  }

  #### Return outputs
  t_end <- Sys.time()
  total_duration <- difftime(t_end, t_onset, units = "mins")
  cat_to_console(paste0("... flapper::pf_simplify() call completed (@ ", t_end, ") after ~", round(total_duration, digits = 2), " minutes."))
  return(out)
}
