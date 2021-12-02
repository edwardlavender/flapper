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
#' @param cl,varlist (optional) Parallelisation options implemented if (a) particle histories are contained in memory or (b) particle histories are supplied as a list of file paths with \code{use_memory_safe = FALSE}. \code{cl} is (a) a cluster object from \code{\link[parallel]{makeCluster}} or (b) an integer that defines the number of child processes. \code{varlist} is a character vector of variables for export (see \code{\link[flapper]{cl_export}}). Exported variables must be located in the global environment. If a cluster is supplied, the connection to the cluster is closed within the function (see \code{\link[flapper]{cl_stop}}). For further information, see \code{\link[flapper]{cl_lapply}} and \code{\link[flapper]{flapper-tips-parallel}}.
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
      cells_by_time <- cl_lapply(history,
                                 cl = cl, varlist = varlist,
                                 function(f) readRDS(f)$id_current)
      cells <- unlist(cells_by_time)
      cells <- unique(cells)

    }

    #### Access unique cells (in memory)
  } else {

    cells <- cl_lapply(history,
                       cl = cl, varlist = varlist,
                       fun = function(d) d$id_current)
    cells <- unlist(cells)
    cells <- unique(cells)

  }

  #### Return unique cells
  return(cells)
}

