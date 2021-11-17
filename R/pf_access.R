######################################
######################################
#### pf_access_history_from_file()

#' @title List `history' files from a PF algorithm
#' @description This function creates an ordered vector (or list) of `history' files derived from the particle filtering (PF) algorithm (\code{\link[flapper]{pf}}). This is applicable if \code{\link[flapper]{pf}} is implemented with the \code{write_history} argument specified.
#'
#' @param root A string that defines the directory in which files are located.
#' @param use_absolute_paths A logical variable that defines whether to return relative paths (\code{FALSE}) or absolute paths (\code{TRUE}) (see \code{\link[tools]{file_path_as_absolute}}).
#' @param use_list A logical variable that defines whether or not return a vector (\code{use_list = FALSE}) or a list (\code{use_list = TRUE}).
#' @param ... Additional arguments passed to \code{\link[base]{list.files}} (excluding \code{full.names}).
#'
#' @details This function requires the \code{\link[stringr]{stringr}} package.
#'
#' @return The function returns an ordered list of file paths.
#'
#' @examples
#' #### Example (1): Example with default arguments
#' # Define a directory in which to save files from PF
#' root <- paste0(tempdir(), "/pf/")
#' dir.create(root)
#' # Implement the PF algorithm with write_history specified
#' # ... For speed, we will implement the algorithm using pre-defined data
#' pf_args <- dat_dcpf_histories$args
#' pf_args$calc_distance_euclid_fast <- TRUE
#' pf_args$write_history             <- list(file = root)
#' do.call(pf, pf_args)
#' # List the files
#' files <- pf_access_history_files(root)
#' utils::head(files)
#'
#' @seealso This function is designed to list outputs from \code{\link[flapper]{pf}} (see the \code{write_history} argument).
#' @author Edward Lavender
#' @export

pf_access_history_files <- function(root, use_absolute_paths = FALSE, use_list = FALSE,...){
  if(!requireNamespace("stringr", quietly = TRUE)){
    stop("This function requires the 'stringr' package. Please install it before continuing with install.packages('stringr').")
  }
  check...("full.names",...)
  check_dir(input = root)
  files <- list.files(root,...)
  if(!grepl("pf_", files[1], fixed = TRUE)){
    stop("File naming structure is unrecognised.", immediate. = TRUE)
  }
  files <- data.frame(index = 1:length(files), name = files)
  files$pf_id <- stringr::str_split_fixed(files$name, "_", 2)[, 2]
  files$pf_id <- substr(files$pf_id, 1, nchar(files$pf_id) - 4)
  files$pf_id <- as.integer(as.character(files$pf_id))
  files <- files %>% dplyr::arrange(.data$pf_id)
  files <- list.files(root, full.names = TRUE,...)[files$index]
  if(use_absolute_paths) {
    files <- sapply(files, function(f) tools::file_path_as_absolute(f))
    names(files) <- NULL
  }
  if(use_list) files <- as.list(files)
  return(files)
}


########################################
########################################
#### pf_access_history

#' @title Access the `history' element of a \code{\link[flapper]{pf_archive-class}} object
#' @description This function accesses and simplifies the `history' list in a \code{\link[flapper]{pf_archive-class}} object.
#' @param archive A \code{\link[flapper]{pf_archive-class}} object.
#' @param bathy (optional) A \code{\link[raster]{raster}} that defines the grid across the area over which particle filtering was applied. If unsupplied, this is extracted from \code{archive} if available.
#' @details From the `history' element of a \code{\link[flapper]{pf_archive-class}} object, this function extracts particle samples as a dataframe with columns for time steps, cell IDs, cell probabilities and coordinates (if \code{bathy} is available).
#' @return The function returns a dataframe that defines, for each time step (`timestep'), particle samples (`cell_id'), associated probabilities (`cell_pr') and, if \code{bathy} is available, cell coordinates (`cell_x', `cell_y' and `cell_z').
#' @examples
#' pf_access_history(dat_dcpf_histories)
#' @author Edward Lavender
#' @export

pf_access_history <- function(archive,
                              bathy = NULL
){
  check_class(input = archive, to_class = "pf_archive")
  if(is.null(bathy)) bathy <- archive$args$bathy
  history <- lapply(1:length(archive$history), function(t){
    elm <- archive$history[[t]]
    if(!rlang::has_name(elm, "timestep")) elm$timestep <- t
    elm <- elm[, c("timestep", "id_current", "pr_current")]
  })
  history <- do.call(rbind, history)
  colnames(history) <- c("timestep", "cell_id", "cell_pr")
  if(!is.null(bathy)){
    history[, c("cell_x", "cell_y")] <- raster::extract(bathy, history$cell_id)
    history$cell_z <- raster::extract(bathy, history$cell_id)
  }
  cols <- c("timestep", "cell_id", "cell_x", "cell_y", "cell_z", "cell_pr")
  history[, cols[cols %in% colnames(history)]]
  history <-
    history %>%
    dplyr::arrange(.data$timestep, .data$cell_id, .data$cell_pr)
  return(history)
}


########################################
########################################
#### pf_access_particles_unique()

#' @title Access the cells sampled by PF
#' @description Given a list of particle histories (or a list of file paths), this function accesses the unique particles (cells) sampled by a particle filtering (PF) algorithm (\code{\link[flapper]{pf}}).
#'
#' @param history A \code{\link[flapper]{pf_archive-class}} class object from \code{\link[flapper]{pf}}, the list of particle histories (the `history' element of a \code{\link[flapper]{pf_archive-class}} object) from \code{\link[flapper]{pf}} or a list of file paths to particle histories.
#' @param use_memory_safe If \code{history} is a record of file paths, \code{use_memory_safe} is a logical variable that defines whether or not to use the `memory-safe(r)' method to access unique cell samples. If specified, the function sequentially loads each file and re-defines the vector of unique particles at each time step as the unique combination of previous (unique) samples and the samples from the current time step. This may be slow. Alternatively, under the default \code{use_memory_safe = FALSE} option, the function loads each file (in parallel if specified), retaining all sampled particles, before selecting the unique particles (once) at the end of this process. This option should be faster.
#' @param cl,varlist Parallelisation options implemented if (a) particle histories are contained in memory or (b) particle histories are supplied as a list of file paths with \code{use_memory_safe = FALSE}. \code{cl} is a cluster object created by \code{\link[parallel]{makeCluster}}. If supplied, the connection to the cluster is closed within the function. \code{varlist} is a character vector of names of objects to export that is passed to the \code{varlist} argument of \code{\link[parallel]{clusterExport}}. Exported objects must be located in the global environment.
#'
#' @return The function returns a vector of the unique particles sampled by the PF algorithm.
#'
#' @examples
#' #### Example (1): Access unique particles when 'history' exists in memory
#' # Access unique particles from a pf_archive object
#' pf_access_particles_unique(dat_dcpf_histories)
#' # Access unique particles from a list of particle histories
#' pf_access_particles_unique(dat_dcpf_histories$history)
#' # Supply a cluster to speed up the algorithm (for very large lists)
#' pf_access_particles_unique(dat_dcpf_histories$history,
#'                            cl = parallel::makeCluster(2L))
#'
#' #### Example (2): Access unique particles when 'history' is a list of file paths
#'
#' ## Write example particle histories to file (to load)
#' root <- paste0(tempdir(), "/pf/")
#' dir.create(root)
#' pf_args <- dat_dcpf_histories$args
#' pf_args$calc_distance_euclid_fast <- TRUE
#' pf_args$write_history <- list(file = root)
#' out_pf <- do.call(pf, pf_args)
#'
#' ## Access particle histories using default options (use_memory_safe = FALSE)
#' # Access particle histories via pf_access_history_files()
#' pf_access_particles_unique(pf_access_history_files(root, use_list = TRUE))
#' # Supply a cluster to speed up the algorithm (for very large lists)
#' pf_access_particles_unique(pf_access_history_files(root, use_list = TRUE),
#'                            cl = parallel::makeCluster(2L))
#'
#' ## Access particle histories using the 'memory_safe' option
#' # For large lists, this is likely to be slower
#' # ... but it may be the only option in some cases.
#' pf_access_particles_unique(pf_access_history_files(root, use_list = TRUE),
#'                            use_memory_safe = TRUE)
#'
#' @seealso \code{\link[flapper]{pf}} implements particle filtering.
#' @author Edward Lavender
#' @export

pf_access_particles_unique <- function(history,
                                       use_memory_safe = FALSE,
                                       cl = NULL, varlist = NULL){

  #### Setup
  check_class(input = history, to_class = "list")
  if(inherits(history, "pf_archive")) history <- history$history
  history_1   <- history[[1]]
  read_history <- FALSE
  if(inherits(history_1, "character")){
    read_history <- TRUE
    if(!file.exists(history_1)) stop(paste0("history[[1]] ('", history_1, "') does not exist."))
    history_1 <- readRDS(history_1)
  }
  if(is.null(cl) & !is.null(varlist)) {
    warning("'cl' is NULL but 'varlist' is not: 'varlist' ignored.",
            immediate. = TRUE, call. = FALSE)
    varlist <- FALSE
  }

  #### Access unique cells (from file)
  if(read_history){

    ## Memory-safe(r) option: load files sequentially, selecting unique cells at each step
    if(use_memory_safe){
      if(!is.null(cl)) {
        warning("'cl' is not implemented for loading files when use_memory_safe = TRUE.",
                immediate. = TRUE, call. = FALSE)
        cl <- varlist <- NULL
      }
      cells <- unique(history_1$id_current)
      for(i in 2:length(history)){
        cells <- unique(c(cells, readRDS(history[[i]])$id_current))
      }

      ## Faster option: load all cells for each time step, selecting unique cells at the end
    } else {
      if(!is.null(cl) & is.null(varlist)) parallel::clusterExport(cl = cl, varlist = varlist)
      cells_by_time <-
        pbapply::pblapply(history,
                          cl = cl,
                          function(f) readRDS(f)$id_current)
      if(!is.null(cl)) parallel::stopCluster(cl = cl)
      cells <- unlist(cells_by_time)
      cells <- unique(cells)

    }

    #### Access unique cells (in memory)
  } else {

    if(!is.null(cl) & is.null(varlist)) parallel::clusterExport(cl = cl, varlist = varlist)
    cells <- pbapply::pblapply(history, cl = cl, function(d) d$id_current)
    if(!is.null(cl)) parallel::stopCluster(cl = cl)
    cells <- unlist(cells)
    cells <- unique(cells)

  }

  #### Return unique cells
  return(cells)
}

######################################
######################################
#### pf_access_distance_matrix()

#' @title Calculate distances between particles from PF
#' @description This function calculates the Euclidean or shortest distances between unique particle samples from a particle filtering (PF) algorithm (\code{\link[flapper]{pf}}). This is designed to support particle processing and the reconstruction of movement paths via \code{\link[flapper]{pf_simplify}} (see Details).
#'
#' @param archive A \code{\link[flapper]{pf_archive-class}} object from \code{\link[flapper]{pf}}.
#' @param bathy (optional) The bathymetry \code{\link[raster]{raster}} (see \code{\link[flapper]{pf}}). If unsupplied, \code{bathy} can be extracted from \code{archive$args$bathy}.
#' @param mobility (optional) The mobility parameter (see \code{\link[flapper]{pf}}). If supplied, distances that exceed \code{mobility} are set to zero in the sparse \code{\link[Matrix]{dgCMatrix-class}} that is returned (see Details). Particularly for shortest-distances implementations (see \code{calc_distance}), supplying \code{mobility} should reduce memory requirements and computation time (see Details).
#' @param calc_distance A character that defines the method used to calculate distances between sequential combinations of particles (see \code{\link[flapper]{pf}}). Currently supported options are Euclidean distances ("euclid") or shortest (least-cost) distances ("lcp"). Note that \code{calc_distance} does not need to be the same method as used for \code{\link[flapper]{pf}}. If \code{calc_distance = NULL}, the method saved in \code{archive} is used.
#' @param calc_distance_graph (optional) If \code{calc_distance = "lcp"}, \code{calc_distance_graph} is a graph object that defines the distances between connected cells in \code{bathy}. If unsupplied, this is taken from \code{archive$args$calc_distance_graph}, if available, or computed via \code{\link[flapper]{lcp_graph_surface}}.
#' @param ... If \code{calc_distance = "lcp"}, \code{...} are additional arguments passed to \code{\link[cppRouting]{get_distance_pair}} to calculate shortest distances. Allowed arguments are \code{algorithm}, \code{constant} and \code{allcores}.
#' @param verbose A logical input that defines whether or not to print messages to the console to monitor function progress.
#'
#' @details This function implements a stepwise approach to calculate distances between all unique pairwise combinations of sampled particles:
#' \enumerate{
#'   \item All unique particle samples are extracted from \code{archive} (via \code{\link[flapper]{pf_access_particles_unique}}).
#'   \item The Euclidean distances between all unique pairwise combinations are calculated (via \code{\link[raster]{pointDistance}}).
#'   \item If \code{calc_distance = "lcp"}, for all pairwise particle combinations for which Euclidean distances are less than \code{mobility} (if supplied), shortest distances are calculated (via \code{\link[cppRouting]{get_distance_pair}}).
#'   \item All pairwise particle combinations for which calculated distances are less than \code{mobility} are represented on a sparse \code{\link[Matrix]{dgCMatrix-class}} matrix, with one row and one column for each cell in \code{bathy}. Pairwise distances between unsampled cells and between sampled cells for which the distance exceeds \code{mobility} are represented as zeros.
#' }
#'
#' The resultant matrix can be passed to \code{\link[flapper]{pf_simplify}} to streamline distance calculations and/or integrate shortest distances in particle processing/movement path reconstruction. In both cases, the stepwise approach implemented by this function helps to minimise the computation time required for distance calculations because (a) only distances between relevant locations are calculated (once) and (b) any Euclidean distances that exceed \code{mobility} can be ignored when it comes to shortest-distance calculations.
#'
#'
#' @return The function returns a \code{\link[Matrix]{dgCMatrix-class}} sparse matrix with one row and one column for each pairwise cell combination in \code{bathy}. By default cells are empty (value zero). For the unique, pairwise cells in \code{archive} (that are less than \code{mobility} apart, if applicable), cell values represent the Euclidean or shortest distance between those cells in metres.
#'
#' @examples
#' #### Example (1): Implement function with the default options
#' # Implement function
#' bathy <- dat_dcpf_histories$args$bathy
#' mat <- pf_access_distance_matrix(dat_dcpf_histories, bathy = bathy)
#' # The function returns a dgCMatrix-class matrix
#' str(mat)
#' # There is one row and one column for each cell in bathy
#' mat@Dim == raster::ncell(bathy)
#' # For cells that were not sampled & for cells that are more than mobility apart
#' # ... (if applicable), distances have not been calculated and are given as zero.
#' # ... For example:
#' mat[4, 5]
#' # For cells that were sampled, either Euclidean or shortest distances are given
#' # ... depending on calc_distance. For example, here is the distance from an example
#' # ... sampled cell to every other cell (where zero is the default 'non calculated')
#' # ... distance:
#' cells <- pf_access_particles_unique(dat_dcpf_histories)
#' mat[cells[1], ] # or equivalently mat[, cells[1]]
#' # Here is the distance between two sampled cells, which we can check as follows:
#' mat[cells[1], cells[2]]
#' mat[cells[2], cells[1]]
#' raster::pointDistance(raster::xyFromCell(bathy, cells[1]),
#'                       raster::xyFromCell(bathy, cells[2]),
#'                       lonlat = FALSE)
#'
#' #### Example (2): Implement the function using shortest distances
#' # Implement algorithm
#' mat <- pf_access_distance_matrix(dat_dcpf_histories,
#'                                  bathy = bathy,
#'                                  calc_distance = "lcp")
#' # Check shortest distances for our example cell pair
#' mat[cells[1], cells[2]]
#' mat[cells[2], cells[1]]
#' lcp_over_surface(origin = raster::xyFromCell(bathy, cells[1]),
#'                  destination = raster::xyFromCell(bathy, cells[2]),
#'                  surface = bathy,
#'                  goal = 1L)$dist_lcp
#'
#' #### Example (3): For shortest distances, supply additional args via ...
#' mat <- pf_access_distance_matrix(dat_dcpf_histories,
#'                                  bathy = bathy,
#'                                  calc_distance = "lcp",
#'                                  allcores = TRUE)
#'
#' @seealso \code{\link[flapper]{pf}}, \code{\link[flapper]{pf_simplify}}
#' @author Edward Lavender
#' @export

pf_access_distance_matrix <- function(archive,
                                      bathy = NULL,
                                      mobility = NULL,
                                      calc_distance = NULL,
                                      calc_distance_graph = NULL,
                                      ...,
                                      verbose = TRUE){

  #### Set up
  t_onset <- Sys.time()
  cat_to_cf <- function(..., show = verbose){
    if(show) cat(paste(..., "\n"))
  }
  cat_to_cf(paste0("flapper::pf_access_distance_matrix() called (@ ", t_onset, ")..."))

  #### Get distance parameters
  check_class(input = archive, to_class = "pf_archive")
  cat_to_cf("... Getting parameters")
  if(is.null(bathy)) bathy <- archive$args$bathy
  if(is.null(bathy)) stop("'bathy' must be supplied via the 'bathy' argument or archive$args$bathy.")
  bathy_xy     <- raster::coordinates(bathy)
  bathy_n_cell <- raster::ncell(bathy)
  if(is.null(mobility)) {
    mobility  <- archive$args$mobility
    if(is.null(mobility)){
      message("'mobility' taken as NULL.")
    } else {
      message("'mobility' taken as ", mobility, " m.")
    }
  }
  if(is.null(calc_distance)){
    calc_distance <- archive$args$calc_distance
  } else {
    calc_distance <- match.arg(calc_distance, c("euclid", "lcp"))
  }
  message("'calc_distance' taken as '", calc_distance, "'.")
  if(is.null(calc_distance_graph)) calc_distance_graph <- archive$args$calc_distance_graph

  #### Get unique cells
  cat_to_cf("... Getting unique particle samples...")
  cells <- pf_access_particles_unique(archive)

  #### Get (unique) cell pairs
  cat_to_cf("... Getting (unique) pairwise particle combinations...")
  pairs <- combn(cells, 2)
  pairs <- data.frame(cell_1 = pairs[1, ],
                      cell_2 = pairs[2, ])

  #### Get Euclidean distances between unique cells
  cat_to_cf("... Calculating Euclidean distances...")
  xy_1 <- bathy_xy[pairs$cell_1, ]
  xy_2 <- bathy_xy[pairs$cell_2, ]
  pairs$dist <- pairs$dist_euclid <- raster::pointDistance(xy_1, xy_2, lonlat = FALSE)
  if(!is.null(mobility)) {
    cat_to_cf("... Filtering by 'mobility'...")
    pairs$dist_euclid[pairs$dist_euclid > mobility] <- NA

  }

  #### Get shortest distances for cells less than mobility apart
  if(calc_distance == "lcp" & !all(is.na(pairs$dist_euclid))){
    ## Set up LCP calculations
    cat_to_cf("... Getting shortest distances...")
    if(is.null(calc_distance_graph)){
      cat_to_cf("... ... Calculating cost-surface for calc_distance = 'lcp'...")
      costs <- lcp_costs(bathy, verbose = verbose)
      cost  <- costs$dist_total
      cat_to_cf("... ... Assembling graph for calc_distance = 'lcp'...")
      calc_distance_graph <- lcp_graph_surface(surface = bathy, cost = cost, verbose = verbose)
    }
    cat_to_cf("... ... Calling cppRouting::get_distance_pair() for LCP calculations...")
    pairs_for_lcp <- pairs %>% dplyr::filter(!is.na(.data$dist_euclid))
    ## Get LCPs
    pairs_for_lcp$dist_lcp <- cppRouting::get_distance_pair(Graph = calc_distance_graph,
                                                            from = pairs_for_lcp$cell_1,
                                                            to = pairs_for_lcp$cell_2)
    if(!is.null(mobility)) pairs_for_lcp$dist_lcp[pairs_for_lcp$dist_lcp > mobility] <- NA
    ## Update pairs
    pairs$dist[pairs_for_lcp$index] <- pairs_for_lcp$dist_lcp
    cat_to_cf("... Filtering by 'mobility'...")

  }

  #### Make (sparse) matrix
  pairs <- pairs %>% dplyr::filter(!is.na(.data$dist))
  mat <- Matrix::Matrix(nrow = bathy_n_cell,
                        ncol = bathy_n_cell,
                        data = 0,
                        sparse = TRUE)
  if(nrow(pairs) > 0){
    mat[cbind(pairs$cell_1, pairs$cell_2)] <- pairs$dist
    mat[cbind(pairs$cell_2, pairs$cell_1)] <- pairs$dist
  } else {
    message("All pairwise cell distances exceed 'mobility':empty sparse matrix returned.")
  }

  #### Return outputs
  t_end <- Sys.time()
  total_duration <- difftime(t_end, t_onset, units = "mins")
  cat_to_cf(paste0("... flapper::pf_access_distance_matrix() call completed (@ ", t_end, ") after ~", round(total_duration, digits = 2), " minutes."))
  return(mat)

}
