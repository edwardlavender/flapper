########################################
########################################
#### pf_setup_movement_pr()

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
#' pr <- pf_setup_movement_pr(1:1000)
#' plot(pr, type = "l", xlab = "Distance (m)", ylab = "Pr(distance)")
#' @seealso This function is used as the default movement model in \code{\link[flapper]{pf}}.
#' @author Edward Lavender
#' @export

pf_setup_movement_pr <- function(distance,...) {
  pr <- stats::plogis(10 + distance * -0.05)
  pr[distance > 500] <- 0
  return(pr)
}


########################################
########################################
#### pf_plot_history()

#' @title Plot particle histories from a PF algorithm
#' @description This function plots the spatiotemporal particle histories from a particle filtering (PF) algorithm (the acoustic-centroid PF, the depth-contour PF or the acoustic-centroid depth-contour PF). This produces, for each time step, a map of the individual's possible locations (from the AC, DC or ACDC algorithm), with sampled locations (derived via the particle filtering routine) overlaid.
#' @param record A \code{\link[flapper]{.pf-class}} object from \code{\link[flapper]{pf}} that contains particle histories.
#' @param time_steps An integer vector that defines the time steps for which to plot particle histories.
#' @param add_surface A named list, passed to \code{\link[prettyGraphics]{pretty_map}}, to customise the appearance of the surface, which shows the set of possible positions that the individual could have occupied at a given time step (from \code{\link[flapper]{ac}}, \code{\link[flapper]{dc}} and \code{\link[flapper]{acdc}}), on each map.
#' @param add_particles A named list, passed to \code{\link[prettyGraphics]{pretty_map}}, to customise the appearance of the particles on each map.
#' @param forwards A logical variable that defines whether or not create plots forwards (i.e., from the first to the last \code{time_steps}) or backwards (i.e., from the last to the first \code{time_steps}).
#' @param prompt A logical input that defines whether or not to pause between plots (\code{prompt = TRUE}).
#' @param ... Plot customisation arguments passed to \code{\link[prettyGraphics]{pretty_map}}.
#'
#' @examples
#' #### Implement pf() algorithm
#' # Here, we use pre-defined outputs for speed
#'
#' #### Example (1): The default implementation
#' pf_plot_history(dat_dcpf_histories, time_steps = 1)
#'
#' #### Example (2): Plot customisation options, e.g.:
#' # Customise bathy via add_bathy()
#' pf_plot_history(dat_dcpf_histories,
#'                 time_steps = 1,
#'                 add_surface = list(col = c(grDevices::topo.colors(2))))
#' # Customise particles via add_particles
#' pf_plot_history(dat_dcpf_histories,
#'                 time_steps = 1,
#'                 add_particles = list(col = "red"))
#' # Pass other arguments to prettyGraphics::pretty_map() via ...
#' pf_plot_history(dat_dcpf_histories,
#'                 time_steps = 1,
#'                 add_polys = list(x = dat_coast, col = "brown"),
#'                 crop_spatial = TRUE)
#'
#' #### Example (3): Plot multiple time steps
#' pp <- graphics::par(mfrow = c(2, 2))
#' pf_plot_history(dat_dcpf_histories, time_steps = 1:4, prompt = FALSE)
#' graphics::par(pp)
#'
#' @return The function returns a plot, for each time step, of all the possible locations of the individual, with sampled locations overlaid.
#'
#' @seealso \code{\link[flapper]{pf}} implements PF. \code{\link[flapper]{pf_simplify}} assembles paths from particle histories. \code{\link[flapper]{pf_plot_map}} creates an overall `probability of use' map from particle histories. \code{\link[flapper]{pf_plot_1d}}, \code{\link[flapper]{pf_plot_2d}} and \code{\link[flapper]{pf_plot_3d}} provide plotting routines for paths. \code{\link[flapper]{pf_loglik}} calculates the log-likelihood of each path.
#'
#' @seealso \code{\link[flapper]{pf}} implements PF. \code{\link[flapper]{pf_simplify}} assembles paths from particle histories. \code{\link[flapper]{pf_plot_history}} visualises particle histories. \code{\link[flapper]{pf_plot_1d}}, \code{\link[flapper]{pf_plot_2d}} and \code{\link[flapper]{pf_plot_3d}} provide plotting routines for paths. \code{\link[flapper]{pf_loglik}} calculates the log-likelihood of each path.
#' @author Edward Lavender
#' @export

pf_plot_history <- function(record,
                            time_steps = 1:length(history),
                            add_surface = list(),
                            add_particles = list(pch = "."),
                            forwards = TRUE,
                            prompt = TRUE,...){
  if(!inherits(record, ".pf")) stop("'record' must be a '.pf' class object.")
  layers           <- record$args$record
  history          <- record$history
  time_steps       <- sort(time_steps)
  if(!forwards) time_steps <- rev(time_steps)
  lapply(time_steps, function(t){
    title <- paste0("Time ", t)
    r <- layers[[t]]
    if(inherits(r, "character")) r <- raster::raster(r)
    add_surface$x <- r
    xy_t <- raster::xyFromCell(r, history[[t]]$id_current)
    add_particles$x <- xy_t[, 1]
    add_particles$y <- xy_t[, 2]
    prettyGraphics::pretty_map(r,
                               add_rasters = add_surface,
                               add_points = add_particles,
                               main = title,
                               verbose = FALSE,...)
    if(prompt * length(time_steps) > 1) readline(prompt = "Press [enter] to continue or [Esc] to exit...")
  })
  return(invisible())
}


######################################
######################################
#### pf_plot_map()

#' @title Plot `probability of use' from a PF algorithm
#' @description This function creates a plot of the `probability of use' across an area based on particles sampled by a particle filtering (PF) algorithm. To implement the function, a \code{\link[flapper]{.pf-class}} object that contains particles (locations) sampled by \code{\link[flapper]{pf}} must be supplied. The function extracts all sampled locations and, for each location, calculates `the probability of use' for that location over the time series. This is returned (invisibly) as a \code{\link[raster]{raster}} and plotted.
#' @param pf A \code{\link[flapper]{.pf-class}} object (from \code{\link[flapper]{pf}}).
#' @param map A \code{\link[raster]{raster}} that defines a grid across the area of interest.
#' @param scale A character that defines how \code{\link[raster]{raster}} values are scaled: \code{"original"} uses the original values; \code{"max"} scales values by the maximum value so that they lie between zero and one; and \code{"sum"} scales values by their sum so that they sum to one.
#' @param add_rasters A named list, passed to \code{\link[prettyGraphics]{pretty_map}}, to customise the appearance of the plotted surface.
#' @param ... Additional arguments passed to \code{\link[prettyGraphics]{pretty_map}}.
#'
#' @details For each location, the 'probability of use' is calculated as the sum of the number of times that the location was sampled, weighted by the associated probabilities of each sample.
#'
#' @return The function invisibly returns a \code{\link[raster]{raster}}, in which each cell contains the `probability of use' score and produces a plot of this surface.
#'
#' @examples
#' #### Example (1): Implement the function with default options
#' # using the example 'dat_dcpf_histories' data
#' pf_plot_map(dat_dcpf_histories, map = dat_dc$args$bathy)
#'
#' #### Example (2): Re-scale the map
#' pf_plot_map(dat_dcpf_histories, map = dat_dc$args$bathy, scale = "max")
#' pf_plot_map(dat_dcpf_histories, map = dat_dc$args$bathy, scale = "sum")
#'
#' #### Example (3): Customise the map
#' pf_plot_map(dat_dcpf_histories, map = dat_dc$args$bathy,
#'             add_rasters = list(col = grDevices::grey.colors(n = 100)),
#'             xlab = "x", ylab = "y")
#'
#' @seealso \code{\link[flapper]{pf}} implements PF. \code{\link[flapper]{pf_simplify}} assembles paths from particle histories. \code{\link[flapper]{pf_plot_history}} visualises particle histories. \code{\link[flapper]{pf_plot_1d}}, \code{\link[flapper]{pf_plot_2d}} and \code{\link[flapper]{pf_plot_3d}} provide plotting routines for paths. \code{\link[flapper]{pf_loglik}} calculates the log-likelihood of each path.
#' @author Edward Lavender
#' @export

pf_plot_map <- function(pf,
                        map,
                        scale = c("original", "max", "sum"),
                        add_rasters = list(),...){
  # Check inputs
  check_class(input = pf, to_class = ".pf")
  scale <- match.arg(scale)
  # Extract particle histories as a single dataframe
  pf_particle_histories <- lapply(pf$history, function(elm) elm[, c("id_current", "pr_current"), drop = FALSE])
  pf_particle_histories <- do.call(rbind, pf_particle_histories)
  # Calculated the frequency with which each cell was sampled (weighted by the probability)
  wt_freq <- stats::aggregate(x = list("wt" = pf_particle_histories$pr_current),
                              by = list("id" = pf_particle_histories$id_current),
                              FUN = sum)
  # Re-scale weighted frequencies e.g. so that the maximum value has a score of one
  if(scale == "max"){
    wt_freq$wt <- wt_freq$wt/max(wt_freq$wt)
  } else if(scale == "sum"){
    wt_freq$wt <- wt_freq$wt/sum(wt_freq$wt)
  }
  # Assign scores to map
  p <- raster::setValues(map, 0)
  p[wt_freq$id] <- wt_freq$wt
  p <- raster::mask(p, map)
  if(!is.null(add_rasters)) add_rasters$x <- p
  prettyGraphics::pretty_map(x = p,
                             add_rasters = add_rasters,...)
  return(invisible(p))
}


########################################
########################################
#### pf_simplify()

#' @title Convert particle histories from \code{\link[flapper]{pf}} into movement paths
#' @description This function is designed to simplify the \code{\link[flapper]{.pf-class}} object from \code{\link[flapper]{pf}} that defines sampled particle histories into a set of movement paths. The function identifies pairs of cells between which movement may have occurred at each time step (if necessary), (re)calculates distances and probabilities between connected cell pairs and then links pairwise movements between cells into a set of possible movement paths.
#' @param record A \code{\link[flapper]{.pf-class}} object from \code{\link[flapper]{pf}}.
#' @param calc_distance A character that defines the method used to calculate distances between sequential combinations of particles (see \code{\link[flapper]{pf}}). Currently supported options are Euclidean distances (\code{"euclid"}) or shortest (least-cost) distances ("lcp"). Note that \code{calc_distance} does not need to be the same method as used for \code{\link[flapper]{pf}}: it is often computationally beneficial to implement \code{\link[flapper]{pf}} using Euclidean distances and then, for the subset of sampled particles, implement \code{\link[flapper]{pf_simplify}} with \code{calc_distance = "lcp"} to re-compute distances using the shortest-distances algorithm, along with the adjusted probabilities. However, for large paths, the quickest option is to implement both functions using \code{calc_distance = "euclid"} and then interpolate shortest paths only for the set of returned paths (see \code{\link[flapper]{lcp_interp}}). If \code{calc_distance = NULL}, the method saved in \code{record} is used.
#' @param bathy (optional) If \code{calc_distance = "lcp"}, \code{bathy} is \code{\link[raster]{raster}} that defines the bathymetry across the area within which the individual could have moved. \code{bathy} must be planar (i.e., Universal Transverse Mercator projection) with units of metres in x, y and z directions (m). The surface's resolution is taken to define the distance between horizontally and vertically connected cells and must be the same in both x and y directions. Any cells with NA values (e.g., due to missing data) are treated as `impossible' to move though by the algorithm (see \code{\link[flapper]{lcp_over_surface}}).
#' @param calc_distance_graph (optional) If \code{calc_distance = "lcp"}, \code{calc_distance_graph} is a graph object that defines the distances between connected cells in \code{bathy}. If unsupplied, this is taken from \code{record$args$calc_distance_graph}, if available, or computed via \code{\link[flapper]{lcp_graph_surface}}.
#' @param max_n_copies (optional) An integer that specifies the maximum number of copies of a sampled cell that are retained at each time stamp. Each copy represents a different route to that cell. By default, all copies (i.e. routes to that cell are retained) via \code{max_n_copies = NULL}. However, in cases where there are a large number of paths through a landscape, the function can run into vector memory limitations during path assembly, so \code{max_n_copies} may need to be set. In this case, at each time step, if there are more than \code{max_n_copies} paths to a given cell, then a subset of these (\code{max_n_copies}) are sampled, according to the \code{sample_method} argument.
#' @param sample_method (optional) If \code{max_n_copies} is supplied, \code{sample_method} is a character that defines the sampling method. Currently supported options are: \code{"random"}, which implements random sampling; \code{"weighted"}, which implements weighted sampling, with random samples taken according to their probability at the current time step; and \code{"max"}, which selects for the top \code{max_n_copies} most likely copies of a given cell according to the probability associated with movement into that cell from the previous location.
#' @param add_origin A logical input that defines whether or not to include the origin in the returned dataframe.
#' @param verbose A logical input that defines whether or not to print messages to the console to monitor function progress.
#' @param ... Additional arguments (none implemented).
#'
#' @details The implementation of this function depends on how \code{\link[flapper]{pf}} has been implemented. Under the default options in \code{\link[flapper]{pf}}, the fast Euclidean distances method is used to sample sequential particle positions, in which case the history of each particle through the landscape is not retained and has to be assembled afterwards. In this case, \code{\link[flapper]{pf_simplify}} calculates the distances between all combinations of cells at each time step, using either a Euclidean distances or shortest distances algorithm according to the input to \code{calc_distance}. Distances are converted to probabilities using the `intrinsic' probabilities associated with each location and the movement models retained in \code{record} from the call to \code{\link[flapper]{pf}} to identify possible movement paths between cells at each time step. Pairwise cell movements are then assembled into complete movement paths. If the fast Euclidean distances method has not been used, then pairwise cell movements are retained by  \code{\link[flapper]{pf}}. In this case, the function simply recalculates distances between sequential cell pairs and the associated cell probabilities, which are used to assemble a set of movement paths.
#'
#' @return The function returns a \code{\link[flapper]{pf-class}} object, which is a dataframe that defines the movement paths.
#'
#' @examples
#' #### Example particle histories
#' # In these examples, we will use the example particle histories included in flapper
#' summary(dat_dcpf_histories)
#'
#' #### Example (1): The default implementation
#' paths_1 <- pf_simplify(dat_dcpf_histories)
#'
#' ## Demonstration that the distance and probabilities calculations are correct
#' # The simple method below works if three conditions are met:
#' # ... The 'intrinsic' probability associated with each cell is the same (as for DC algorithm);
#' # ... Paths have been reconstructed via pf_simplify() using Euclidean distances;
#' # ... The calc_movement_pr() movement model applies to all time steps;
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
#' pf_plot_1d(dat_dc$args$archival, paths_1)
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
#' #### Example (2): Re-calculate distances using another method
#' # Use shortest distances:
#' paths_2a <- pf_simplify(dat_dcpf_histories, calc_distance = "lcp")
#' # Speed up shortest distance calculations by supplying the graph object:
#' costs <- lcp_costs(dat_dcpf_histories$args$bathy)
#' graph <- lcp_graph_surface(dat_dcpf_histories$args$bathy, costs$dist_total)
#' paths_2b <- pf_simplify(dat_dcpf_histories,
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
#'                           sample_method = "random")
#' paths_4b <- pf_simplify(dat_dcpf_histories,
#'                           max_n_copies = 5,
#'                           sample_method = "weighted")
#' paths_4c <- pf_simplify(dat_dcpf_histories,
#'                           max_n_copies = 5,
#'                           sample_method = "max")
#' # Compare retained paths
#' pf_loglik(paths_3a)
#' pf_loglik(paths_3b)
#' pf_loglik(paths_3c)
#'
#' #### Example (5): Retain/drop the origin, if specified
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
#' # If specified, the origin is dropped with add_origin = FALS
#' paths_5b <- pf_simplify(dat_dcpf_histories, add_origin = FALSE)
#' head(paths_5b)
#'
#' @author Edward Lavender
#' @export

pf_simplify <- function(record,
                        calc_distance = NULL,
                        bathy = NULL,
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
  cat_to_console(paste0("flapper::pf_simplify() called (@ ", t_onset, ")..."))
  if(!inherits(record, ".pf")) stop("'record' must be a '.pf-class' object.")
  history <- record$history
  layers  <- record$args$record
  layers_1   <- layers[[1]]
  read_layers <- FALSE
  if(inherits(layers_1, "character")){
    read_layers <- TRUE
    if(!file.exists(layers_1)) stop(paste0("record[[1]] ('", layers_1, "') does not exist."))
    layers_1     <- raster::raster(layers_1)
  }
  data    <- record$args$data
  origin   <- record$args$origin
  if(!is.null(origin)) {
    origin_cell_id <- raster::cellFromXY(layers_1, origin)
    origin         <- raster::xyFromCell(layers_1, origin_cell_id)
  }
  if(is.null(calc_distance)) calc_distance <- record$args$calc_distance
  calc_distance                <- match.arg(calc_distance, c("euclid", "lcp"))
  calc_movement_pr_from_origin <- record$args$calc_movement_pr_from_origin
  calc_movement_pr             <- record$args$calc_movement_pr
  mobility                     <- record$args$mobility
  mobility_from_origin         <- record$args$mobility_from_origin
  sample_method                <- match.arg(sample_method)
  if(is.null(bathy)) bathy <- record$args$bathy


  ########################################
  #### Identify sequential, pairwise connections between cells

  ## Implement setup for LCP distance calculations, if necessary.
  cat_to_console(paste("... Getting pairwise cell movements based on calc_distance = ", calc_distance, "..."))
  if(calc_distance == "lcp"){
    cat_to_console("... Setting up LCP calculations...")
    if(is.null(bathy)) stop("'bathy' must be supplied if calc_distance = 'lcp'. ")
    if(!all.equal(raster::res(bathy)[1], raster::res(bathy)[2])){
      stop("raster::res(bathy)[1] must equal raster::res(bathy)[2] for this option.")
    }
    raster_comparison <- tryCatch(raster::compareRaster(layers_1, bathy), error = function(e) return(e))
    if(inherits(raster_comparison, "error")){
      warning("record$args$record[[1]] and 'bathy' have different properties",
              immediate. = TRUE, call. = FALSE)
      stop(raster_comparison)
    }
    if(is.null(calc_distance_graph)) calc_distance_graph <- record$args$calc_distance_graph
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
    layers <- append(raster::setValues(layers_1, 1), layers)

    #### Re-define history with cell -pairs
    history <- pbapply::pblapply(2:length(history), function(t){

      ## Get full suite of allowed cells from current time step
      layers_2 <- layers[[t]]
      if(read_layers) layers_2 <- raster::raster(layers_2)

      ## Get retained cells from the previous and current step and coordinates
      z1 <- history[[t - 1]]
      z1_xy <- raster::xyFromCell(layers_2, z1$id_current)
      z2 <- history[[t]]
      z2_xy <- raster::xyFromCell(layers_2, z2$id_current)

      ## Calculate Euclidean distances between all combinations of cells to identify possible movement paths between cell pairs
      # For both the origin (if applicable) and subsequent locations, we will calculate distances from cells IDs
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
      # ... For the first time step, probabilities depend on
      # ... ... If an origin has been specified then (a) distance from origin and (b) intrinsic probabilities
      # ... ... If an origin has not been specified then just the (b) intrinsic probabilities
      # ... For later time steps, probabilities depend on
      # ... ... (a) distances (and any other parameters in archival)
      # ... ... (b) intrinsic probabilities associated with cells (e.g., due to detection pr)
      # (a) Get probabilities based on movement
      if((t - 1) == 1){
        if(!is.null(origin)){
          pr_btw_cells <- calc_movement_pr_from_origin(dist_btw_cells, data[t-1, ])
        } else {
          pr_btw_cells <- matrix(c(rep(1, nrow(z2))), nrow = 1)
        }
      } else {
        pr_btw_cells <- calc_movement_pr(dist_btw_cells, data[t-1, ])
      }
      # (b) Combine with 'intrinsic' probabilities in each cell
      pr_intrinsic_cells <- raster::extract(layers_2, z2$id_current)
      pr_intrinsic_cells <- matrix(pr_intrinsic_cells, nrow = 1)
      pr_intrinsic_cells <- pr_intrinsic_cells[rep(1:nrow(pr_intrinsic_cells), times = nrow(pr_btw_cells)), ]
      pr_of_cells        <- pr_btw_cells * pr_intrinsic_cells

      ## Define dataframe from matrix that includes all possible movements between sampled locations at t-1 and t
      tmp <- data.frame(id_previous = NA,
                        pr_previous = NA,
                        row = as.vector(row(pr_of_cells)),
                        col = as.vector(col(pr_of_cells)),
                        id_current = NA,
                        pr_current = as.vector(pr_of_cells),
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

    # ... (2) Method for pf outputs not derived via calc_distance_euclid_fast
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
      } else{
        d[, c("id_previous_x", "id_previous_y")]  <- raster::xyFromCell(layers_2, d$id_previous)
        d[, c("id_current_x", "id_current_y")]    <- raster::xyFromCell(layers_2, d$id_current)
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
  path_df[, c("cell_x", "cell_y")] <- raster::xyFromCell(layers_1, path_df$cell_id)
  if(!is.null(origin) & add_origin){
    path_df[path_df$timestep == 0, "cell_x"] <- origin[1]
    path_df[path_df$timestep == 0, "cell_y"] <- origin[2]
  }
  if(!is.null(bathy)) {
    path_df[, "cell_z"] <- raster::extract(bathy, path_df[, c("cell_x", "cell_y")])
  } else path_df[, "cell_z"] <- NA
  path_df <- path_df[, c("path_id", "timestep", "cell_id", "cell_x", "cell_y", "cell_z", "cell_pr", "dist")]
  class(path_df) <- c(class(path_df), "pf")

  #### Return outputs
  t_end <- Sys.time()
  total_duration <- difftime(t_end, t_onset, units = "mins")
  cat_to_console(paste0("... flapper::pf_simplify() call completed (@ ", t_end, ") after ~", round(total_duration, digits = 2), " minutes."))
  return(path_df)
}


########################################
########################################
#### pf_loglik()

#' @title Calculate the log-likelihood of movement paths from a PF algorithm
#' @importFrom rlang .data
#' @description This function calculates the total log-likelihood of each movement path reconstructed by a particle filtering (PF) algorithm, including the acoustic-centroid (AC), depth-contour (DC) or acoustic-centroid depth-contour (ACDC) algorithms.
#' @param paths A dataframe containing movement paths from \code{\link[flapper]{pf}} plus \code{\link[flapper]{pf_simplify}} (see \code{\link[flapper]{pf-class}}). At a minimum, this should contain a unique identifier for each path (named `path_id') and the probability associated with each cell along each path (`cell_pr').
#' @param ... Additional arguments (none implemented).
#' @details For each path, at each time step the probability associated with the sampled location depends on (a) the `intrinsic' probability associated with each cell (assigned by the AC, DC or ACDC algorithm) and (b) a user-defined movement model that is driven by the distance between the sampled locations for the individual at the previous and current time steps (and other user-defined parameters). This function simply sums the logarithms of these probabilities for each path as a measure of their relative likelihood, given the movement model.
#' @examples
#' # An example with the DCPF paths dataset included in flapper
#' pf_loglik(dat_dcpf_paths)
#' @return The function returns a dataframe with the log likelihood (`loglik') of each path (`path_id'). Rows are ordered by log-likelihood and a `delta' column is provided with the differences in log-likelihood between the most likely path and every other path.
#' @author Edward Lavender
#' @export
#'

pf_loglik <- function(paths,...){
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
#### pf_plot_1d()

#' @title Plot one-dimensional depth time series from a PF algorithm
#' @description This function plots the observed depth time series and the depth time series associated with each path reconstructed by the depth-contour particle filtering (DCPF) or acoustic-centroid depth-contour particle filtering (ACDCPF) algorithm.
#' @param archival A dataframe of depth (m) observations named `depth', as used by \code{\link[flapper]{dc}} and \code{\link[flapper]{acdc}}.
#' @param paths A dataframe containing reconstructed movement path(s) from \code{\link[flapper]{pf}} via \code{\link[flapper]{pf_simplify}} (see \code{\link[flapper]{pf-class}}). At a minimum, this should contain a unique identifier for each path (named `path_id’), timesteps (`timestep') and the depth associated with each cell along each path (`cell_z').
#' @param scale A number that vertically scales the depth time series for the observations and the reconstructed path(s). By default, absolute values for depth are assumed and negated for ease of visualisation.
#' @param pretty_axis_args,xlab,ylab,type,... Plot customisation arguments passed to \code{\link[prettyGraphics]{pretty_plot}}.
#' @param add_lines A named list, passed to \code{\link[graphics]{lines}}, to customise the appearance of the depth time series for reconstructed path(s).
#' @param prompt A logical input that defines whether or not plot the observed depth time series with each reconstructed depth time series on a separate plot, sequentially, with a pause between plots (\code{prompt = TRUE}), or with all reconstructed time series on a single plot (\code{prompt = FALSE}).
#' @details Observed and reconstructed depth time series can differ due to measurement error, which is controlled via the \code{calc_depth_error} function in the DC and ACDC algorithms (see \code{\link[flapper]{dc}} and \code{\link[flapper]{acdc}}).
#' @examples
#' #### Implement pf() algorithm
#' # Here, we use pre-defined outputs for speed
#' archival <- dat_dc$args$archival
#' paths    <- dat_dcpf_paths
#'
#' #### Example (1): The default implementation
#' pf_plot_1d(archival, paths)
#'
#' #### Example (2): Plot customisation options, e.g.:
#' pf_plot_1d(archival, paths, scale = 1, pretty_axis_args = list(side = 1:2))
#' pf_plot_1d(archival, paths, type = "l")
#' pf_plot_1d(archival, paths, add_lines = list(col = "red", lwd = 0.5))
#'
#' #### Example (3): Plot individual comparisons
#' if(interactive()){
#'   pp <- graphics::par(mfrow = c(3, 4))
#'   pf_plot_1d(depth, paths, prompt = TRUE)
#'   graphics::par(pp)
#' }
#'
#' @return The function returns a plot of the observed and reconstructed depth time series, either for all paths at once (if \code{prompt = FALSE}) or each path separately (if \code{prompt = TRUE}).
#' @seealso \code{\link[flapper]{pf}} implements the pf algorithm. \code{\link[flapper]{pf_plot_history}} visualises particle histories, \code{\link[flapper]{pf_plot_map}} creates an overall `probability of use' map from particle histories and \code{\link[flapper]{pf_simplify}} processes the outputs into a dataframe of movement paths. \code{\link[flapper]{pf_plot_1d}}, \code{\link[flapper]{pf_plot_2d}} and \code{\link[flapper]{pf_plot_3d}} provide plotting routines for paths. \code{\link[flapper]{pf_loglik}} calculates the log-likelihood of each path.
#' @author Edward Lavender
#' @export

pf_plot_1d <- function(archival,
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
  if(any(is.na(paths$cell_z))) stop("paths$cell_z contains NAs.")
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
#### pf_plot_2d()

#' @title Map two-dimensional paths from a PF algorithm
#' @description This function is a simple wrapper for \code{\link[prettyGraphics]{pretty_map}} that maps the paths reconstructed by a particle filtering (PF) algorithm over a surface.
#' @param paths A dataframe containing reconstructed movement path(s) from \code{\link[flapper]{pf}} via \code{\link[flapper]{pf_simplify}} (see \code{\link[flapper]{pf-class}}). At a minimum, this should contain a unique identifier for each path (named `path_id') and the x and y coordinates that define the trajectory of each path (`cell_x' and `cell_y').
#' @param bathy A \code{\link[raster]{raster}} of the surface over which movement was reconstructed.
#' @param add_paths A named list, passed to \code{\link[prettyGraphics]{add_sp_path}}, to customise the appearance of the paths.
#' @param prompt A logical input that defines whether or not plot each path on a separate plot, sequentially, with a pause between plots (\code{prompt = TRUE}), or all paths on a single plot (\code{prompt = FALSE}).
#' @param ... Additional arguments, passed to \code{\link[prettyGraphics]{pretty_map}}, for plot customisation.
#' @return The function maps the trajectories of reconstructed paths across the surface, returning a single map if \code{prompt = FALSE} or one map for each path if \code{prompt = TRUE}.
#' @examples
#' #### Implement pf() algorithm
#' # Here, we use pre-defined outputs for speed
#' bathy <- dat_dcpf_histories$args$bathy
#' paths <- dat_dcpf_paths
#'
#' #### Example (1): The default implementation
#' pf_plot_2d(paths, bathy)
#'
#' #### Example (2): Plot customisation options
#' # Customise the appearance of the path(s)
#' pf_plot_2d(paths, bathy,
#'              add_paths = list(length = 0.075, col = viridis::viridis(100)))
#' # Pass arguments to prettyGraphics::pretty_map() via ... , e.g.:
#' pf_plot_2d(paths, bathy, xlab = "Easting (UTM)", ylab = "Northing (UTM)")
#'
#' #### Example (3): Plot individual paths separately
#' if(interactive()){
#'   pp <- graphics::par(mfrow = c(3, 4))
#'   pf_plot_2d(paths, bathy, add_paths = list(length = 0.01),
#'                prompt = TRUE, verbose = FALSE)
#'   graphics::par(pp)
#' }
#'
#' @seealso \code{\link[flapper]{pf}} implements the pf algorithm. \code{\link[flapper]{pf_plot_history}} visualises particle histories, \code{\link[flapper]{pf_plot_map}} creates an overall `probability of use' map from particle histories and \code{\link[flapper]{pf_simplify}} processes these into a dataframe of movement paths. \code{\link[flapper]{pf_plot_1d}}, \code{\link[flapper]{pf_plot_2d}} and \code{\link[flapper]{pf_plot_3d}} provide plotting routines for paths. For mapping, it can be useful to interpolate shortest (least-cost) paths between sequential locations via \code{\link[flapper]{lcp_interp}}. \code{\link[flapper]{pf_loglik}} calculates the log-likelihood of each path.
#' @author Edward Lavender
#' @export
#'

pf_plot_2d <- function(paths,
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
#### pf_plot_3d()

#' @title Map three-dimensional paths from a PF algorithm
#' @description This function is a simple wrapper for \code{\link[prettyGraphics]{pretty_scape_3d}} that maps the paths reconstructed by the depth-contour or acoustic-centroid depth-contour particle filtering algorithms (DCPF and ACDCPF) over a surface in three dimensions.
#' @param paths A dataframe containing reconstructed movement path(s) from \code{\link[flapper]{pf}} via \code{\link[flapper]{pf_simplify}} (see \code{\link[flapper]{pf-class}}). At a minimum, this should contain a unique identifier for each path (named `path_id') and the x, y and z coordinates that define the trajectory of each path (`cell_x', `cell_y' and `cell_z').
#' @param bathy A \code{\link[raster]{raster}} of the bathymetry surface over which movement was reconstructed.
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
#' #### Implement pf() algorithm
#' # Here, we use pre-defined outputs for speed
#' # Note that it may be beneficial to interpolate paths between points
#' # ... e.g., via lcp_interp() prior to plotting, but we will not do that here.
#' bathy <- dat_dcpf_histories$args$bathy
#' paths <- dat_dcpf_paths
#'
#' #### Example (1): Visualise paths using the default options
#' pf_plot_3d(paths, bathy)
#'
#' #### Example (2): Customise the plot
#' # Customise via add_paths() list
#' pf_plot_3d(paths, bathy,
#'            add_paths = list(line = list(color = "black", width = 10),
#'                             marker = list(color = "blue", size = 10)))
#' # Adjust shift, stretch or aspectmode
#' pf_plot_3d(paths, bathy, shift = 200, stretch = -10)
#' # Customise via ... e.g., add coastline:
#' coast <- raster::crop(dat_coast, bathy)
#' pf_plot_3d(paths, bathy, coastline = coast)
#' # The returned plot objects can also be used for further customisation.
#'
#' #### Example (3): Plot individual paths separately
#' if(interactive()) {
#'   pf_plot_3d(paths, bathy, prompt = TRUE)
#' }
#'
#' @details This function requires the \code{\link[plotly]{plotly}} package.
#'
#' @seealso \code{\link[flapper]{pf}} implements the pf algorithm. \code{\link[flapper]{pf_plot_history}} visualises particle histories, \code{\link[flapper]{pf_plot_map}} creates an overall `probability of use' map from particle histories and \code{\link[flapper]{pf_simplify}} processes these into a dataframe of movement paths. \code{\link[flapper]{pf_plot_1d}}, \code{\link[flapper]{pf_plot_2d}} and \code{\link[flapper]{pf_plot_3d}} provide plotting routines for paths. For mapping, it can be useful to interpolate shortest (least-cost) paths between sequential locations via \code{\link[flapper]{lcp_interp}}. \code{\link[flapper]{pf_loglik}} calculates the log-likelihood of each path.
#'
#' @author Edward Lavender
#' @export

pf_plot_3d <- function(paths,
                       bathy,
                       add_paths = list(line = list(width = 10)),
                       shift = 5,
                       stretch = -5,
                       aspectmode = "data",
                       prompt = FALSE,...){
  # Checks
  if(!requireNamespace("plotly", quietly = TRUE)) stop("This function requires the 'plotly' package. Please install it with `install.packages('plotly')` first.")
  check_names(input = paths, req = c("path_id", "cell_x", "cell_y", "cell_z"))
  if(any(is.na(paths$cell_z))) stop("paths$cell_z contains NAs.")
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
