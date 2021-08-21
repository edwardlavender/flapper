########################################
########################################
#### pf_simplify()

#' @title Convert particle histories from \code{\link[flapper]{pf}} into movement paths
#' @description This function is designed to simplify the \code{\link[flapper]{pf_archive-class}} object from \code{\link[flapper]{pf}} that defines sampled particle histories into a set of movement paths. The function identifies pairs of cells between which movement may have occurred at each time step (if necessary), (re)calculates distances and probabilities between connected cell pairs and then links pairwise movements between cells into a set of possible movement paths.
#' @param archive A \code{\link[flapper]{pf_archive-class}} object from \code{\link[flapper]{pf}}.
#' @param calc_distance A character that defines the method used to calculate distances between sequential combinations of particles (see \code{\link[flapper]{pf}}). Currently supported options are Euclidean distances (\code{"euclid"}) or shortest (least-cost) distances ("lcp"). Note that \code{calc_distance} does not need to be the same method as used for \code{\link[flapper]{pf}}: it is often computationally beneficial to implement \code{\link[flapper]{pf}} using Euclidean distances and then, for the subset of sampled particles, implement \code{\link[flapper]{pf_simplify}} with \code{calc_distance = "lcp"} to re-compute distances using the shortest-distances algorithm, along with the adjusted probabilities. However, for large paths, the quickest option is to implement both functions using \code{calc_distance = "euclid"} and then interpolate shortest paths only for the set of returned paths (see \code{\link[flapper]{lcp_interp}}). If \code{calc_distance = NULL}, the method saved in \code{archive} is used.
#' @param bathy (optional) If \code{calc_distance = "lcp"}, \code{bathy} is \code{\link[raster]{raster}} that defines the bathymetry across the area within which the individual could have moved. \code{bathy} must be planar (i.e., Universal Transverse Mercator projection) with units of metres in x, y and z directions (m). The surface's resolution is taken to define the distance between horizontally and vertically connected cells and must be the same in both x and y directions. Any cells with NA values (e.g., due to missing data) are treated as `impossible' to move though by the algorithm (see \code{\link[flapper]{lcp_over_surface}}).
#' @param calc_distance_graph (optional) If \code{calc_distance = "lcp"}, \code{calc_distance_graph} is a graph object that defines the distances between connected cells in \code{bathy}. If unsupplied, this is taken from \code{archive$args$calc_distance_graph}, if available, or computed via \code{\link[flapper]{lcp_graph_surface}}.
#' @param max_n_copies (optional) An integer that specifies the maximum number of copies of a sampled cell that are retained at each time stamp. Each copy represents a different route to that cell. By default, all copies (i.e. routes to that cell are retained) via \code{max_n_copies = NULL}. However, in cases where there are a large number of paths through a landscape, the function can run into vector memory limitations during path assembly, so \code{max_n_copies} may need to be set. In this case, at each time step, if there are more than \code{max_n_copies} paths to a given cell, then a subset of these (\code{max_n_copies}) are sampled, according to the \code{sample_method} argument.
#' @param sample_method (optional) If \code{max_n_copies} is supplied, \code{sample_method} is a character that defines the sampling method. Currently supported options are: \code{"random"}, which implements random sampling; \code{"weighted"}, which implements weighted sampling, with random samples taken according to their probability at the current time step; and \code{"max"}, which selects for the top \code{max_n_copies} most likely copies of a given cell according to the probability associated with movement into that cell from the previous location.
#' @param add_origin A logical input that defines whether or not to include the origin in the returned dataframe.
#' @param verbose A logical input that defines whether or not to print messages to the console to monitor function progress.
#' @param ... Additional arguments (none implemented).
#'
#' @details The implementation of this function depends on how \code{\link[flapper]{pf}} has been implemented. Under the default options in \code{\link[flapper]{pf}}, the fast Euclidean distances method is used to sample sequential particle positions, in which case the history of each particle through the landscape is not retained and has to be assembled afterwards. In this case, \code{\link[flapper]{pf_simplify}} calculates the distances between all combinations of cells at each time step, using either a Euclidean distances or shortest distances algorithm according to the input to \code{calc_distance}. Distances are converted to probabilities using the `intrinsic' probabilities associated with each location and the movement models retained in \code{archive} from the call to \code{\link[flapper]{pf}} to identify possible movement paths between cells at each time step. Pairwise cell movements are then assembled into complete movement paths. If the fast Euclidean distances method has not been used, then pairwise cell movements are retained by  \code{\link[flapper]{pf}}. In this case, the function simply recalculates distances between sequential cell pairs and the associated cell probabilities, which are used to assemble a set of movement paths.
#'
#' @return The function returns a \code{\link[flapper]{pf_path-class}} object, which is a dataframe that defines the movement paths.
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

pf_simplify <- function(archive,
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
  if(!inherits(archive, "pf_archive")) stop("'archive' must be a 'pf_archive-class' object.")
  history <- archive$history
  layers  <- archive$args$record
  layers_1   <- layers[[1]]
  read_layers <- FALSE
  if(inherits(layers_1, "character")){
    read_layers <- TRUE
    if(!file.exists(layers_1)) stop(paste0("archive$args$record[[1]] ('", layers_1, "') does not exist."))
    layers_1     <- raster::raster(layers_1)
  }
  data     <- archive$args$data
  origin   <- archive$args$origin
  if(!is.null(origin)) {
    origin_cell_id <- raster::cellFromXY(layers_1, origin)
    origin         <- raster::xyFromCell(layers_1, origin_cell_id)
  }
  if(is.null(calc_distance)) calc_distance <- archive$args$calc_distance
  calc_distance                <- match.arg(calc_distance, c("euclid", "lcp"))
  calc_movement_pr_from_origin <- archive$args$calc_movement_pr_from_origin
  calc_movement_pr             <- archive$args$calc_movement_pr
  mobility                     <- archive$args$mobility
  mobility_from_origin         <- archive$args$mobility_from_origin
  sample_method                <- match.arg(sample_method)
  if(is.null(bathy)) bathy <- archive$args$bathy


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
      warning("archive$args$record[[1]] and 'bathy' have different properties",
              immediate. = TRUE, call. = FALSE)
      stop(raster_comparison)
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
  class(path_df) <- c(class(path_df), "pf_path")

  #### Return outputs
  t_end <- Sys.time()
  total_duration <- difftime(t_end, t_onset, units = "mins")
  cat_to_console(paste0("... flapper::pf_simplify() call completed (@ ", t_end, ") after ~", round(total_duration, digits = 2), " minutes."))
  return(path_df)
}
