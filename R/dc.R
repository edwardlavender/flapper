#' @title The depth-contour (DC) algorithm
#' @description This function implements the depth-contour (DC) algorithm. Under the assumption that individuals are benthic/demersal, this algorithm relates one-dimensional depth time series to a two-dimensional bathymetry surface to determine the extent to which different parts of an area might have (or have not) been used, or effectively represent occupied depths, over time. Given a sequence of depth observations (\code{archival}) from a benthic animal and a measurement error parameter (\code{calc_depth_error}), at each time step the function determines the cells on a bathymetry \code{\link[raster]{raster}} (\code{bathy}) that match the observed depth. Across all time steps, matches are summed to produce a single map representing the number of occasions when the depth in each cell matched the observed depth.
#'
#' @param archival A dataframe of depth time series (for a single individual). At a minimum, this should contain a column named `depth' with depth observations. Depth should be recorded using absolute values in the same units as the bathymetry (\code{bathy}, see below).
#' @param bathy A \code{\link[raster]{raster}} of the bathymetry in an area within which the animal is likely to have been located over the study. Bathymetry values should be recorded as absolute values and in the same units as for depths (see \code{archival}).
#' @param calc_depth_error A function that returns the depth error around a given depth. This should accept a single depth value (from \code{archival$depth}) and return two numbers that, when added to that depth, define the range of depths on the bathymetry raster (\code{bathy}) that the individual could plausibly have occupied at any time, given its depth. Since the depth errors are added to the individual's depth, the first number should be negative (i.e., the individual could have been slightly shallower that observed) and the second positive (i.e., the individual could have been slightly deeper than observed). For example, the constant function \code{calc_depth_error = function(...) c(-2.5, 2.5)} implies that the individual could have occupied bathymetric cells whose depth lies within the interval defined by the observed depth + (-2.5) and + (+2.5) m. The appropriate form for \code{calc_depth_error} depends on measurement error for the depth observations in \code{archival} and bathymetry (\code{bathy}) data, as well as the tidal range (m) across the area (over the duration of observations), but this implementation allows the depth error to depend on depth and for the lower and upper error around an observation to differ.
#' @param check_availability A logical input that defines whether or not to record explicitly, for each time step, whether or not there were any cells on \code{bathy} that matched the observed depth.
#' @param plot An integer vector that defines the time steps for which to return time step-specific and cumulative maps of the individual's possible locations. \code{plot = 0} suppresses the return of this information and \code{plot = NULL} returns this information for all time steps.
#' @param write_history (optional) A named list, passed to \code{\link[raster]{writeRaster}}, to save the \code{\link[raster]{raster}} of the individual's possible positions at each time step to file. The `filename' argument should be the directory in which to save files. Files are named by archival time steps as 'arc_1', 'arc_2' and so on.
#' @param split,cl,varlist (optional) Parallelisation arguments. \code{split} is an integer which, if supplied, splits the \code{archival} dataframe every \code{n}th row into chunks*. The algorithm is applied sequentially within each chunk (if applicable) and chunk-wise maps are summed afterwards to create a single map of space use. The advantage of this approach is that chunks can be analysed in parallel, via \code{cl} and \code{varlist}, while memory use is minimised. \code{cl} is a cluster object created by \code{\link[parallel]{makeCluster}} to implement the algorithm in parallel. \code{varlist} is a character vector of names of objects to export, to be passed to the \code{varlist} argument of \code{\link[parallel]{clusterExport}}. This may be required if \code{cl} is supplied. Exported objects must be located in the global environment.
#' @param verbose A logical variable that defines whether or not to print messages to the console to relay function progress.
#' @param ... Additional arguments (none implemented).
#'
#' @details *Under the default options (\code{split = NULL}), the function starts with a blank map of the area and iterates over each time step, adding the `possible positions' of the individual to the map at each step. By continuously updating a single map, this approach is slow but minimises memory requirements. An alternative approach is to split the time series into chunks, implement an iterative approach within each chunk, and then maps for each chunk. This is implemented by \code{split}.
#'
#' @return The function returns a named list with the following elements. `args' is a named list of the arguments used to call the function. `archival' is the \code{archival} dataframe. This includes two new columns that define the lower and upper bounds for the possible depth of the individual on \code{bathy} at each time step (`depth_lwr' and `depth_upper') derived from \code{calc_depth_error} and, if \code{calc_availability = TRUE}, a logical vector that defines whether or not there are any cells on \code{bathy} of the required depth range at each time step. `spatial' is a named list, if \code{plot != 0}, which contains time step-specific (`map_timestep') and cumulative maps (`map_cumulative') of the individual's possible locations for each specified time step. `dc' is a \code{\link[raster]{raster}}, with the same properties as \code{bathy}, in which the value of each cell is the number of times that the depth in that cell overlapped with the individual's depth.
#'
#' @examples
#' #### Define depth time series for examples
#' # We will use a sample depth time series for one individual
#' # We will select a small sample of observations for example speed
#' depth <- dat_archival[dat_archival$individual_id == 25, ][1:100, ]
#'
#' #### Example (1): Implement algorithm with default options
#' dc_out <- dc(archival = depth,
#'              bathy = dat_gebco)
#' # Examine map
#' # Each cell shows the number of time steps when the bathymetry data in each cell
#' # ... matches the archival data
#' prettyGraphics::pretty_map(add_rasters = list(x = dc_out$dc),
#'                            add_polys = list(x = dat_coast))
#' # Convert counts on map to percentages
#' prettyGraphics::pretty_map(add_rasters = list(x = dc_out$dc/nrow(depth) * 100,
#'                                               zlim = c(0, 100)),
#'                            add_polys = list(x = dat_coast))
#' # Check for occasions when the individual's depth was not consistent
#' # ... with the depth data for the area and the depth error e.g., possibly
#' # ... due to movement beyond this area:
#' any(dc_out$archival$availability == FALSE)
#'
#' #### Example (2): Implement the algorithm in parallel
#' # Trial different options for 'split' and compare speed
#' system.time(
#'   dc_out <- dc(archival = depth,
#'                bathy = dat_gebco,
#'                split = 1,
#'                cl = parallel::makeCluster(2L)
#'   )
#' )
#' system.time(
#'   dc_out <- dc(archival = depth,
#'                bathy = dat_gebco,
#'                split = 5,
#'                cl = parallel::makeCluster(2L)
#'   )
#' )
#'
#' @seealso \code{\link[flapper]{dcq}} implements a faster version of this algorithm termed the `quick depth-contour' (DCQ) algorithm. Rather than considering the depth interval that the individual could have occupied at each time step, the DCQ algorithm considers a sequence of depth bins (e.g., 10 m bins), isolates these on the bathymetry \code{\link[raster]{raster}} (\code{bathy}) and counts the number of matches in each cell. The DCPF algorithm (see \code{\link[flapper]{dcpf}}) extends the DC algorithm via particle filtering to reconstruct possible movement paths over \code{bathy}. The ACDC algorithm (see \code{\link[flapper]{acdc}}) extends the depth-contour algorithm by integrating information from acoustic detections of individuals at each time step to restrict the locations in which depth contours are identified.
#'
#' @author Edward Lavender
#' @export

dc <- function(archival,
               bathy,
               calc_depth_error = function(...) c(-2.5, 2.5),
               check_availability = TRUE,
               plot = 1L,
               write_history = NULL,
               split = NULL,
               cl = NULL, varlist = NULL,
               verbose = TRUE,...){

  #### Set up function
  # Function onset
  cat_to_console <- function(..., show = verbose){
    if(show) cat(paste(..., "\n"))
  }
  t_onset <- Sys.time()
  cat_to_console(paste0("flapper::dc() called (@ ", t_onset, ")..."))
  cat_to_console("... Setting up function...")

  #### Define storage container
  out <- list(args = NULL, archival = NULL, spatial = NULL, dc = NULL)
  out$args <- list(archival = archival,
                   bathy = bathy,
                   calc_depth_error = calc_depth_error,
                   split = split,
                   cl = cl, varlist = varlist,
                   verbose = verbose,
                   dots = list(...))

  #### Checks
  # Check arguments passed via ...
  if(!is.null(names(list(...)))){
    warning(paste0("The following argument(s) passed via ... are not supported: ",
                   paste(names(list(...)), collapse = ", "), "."),
            call. = FALSE, immediate. = TRUE)
  }
  # Check archival dataframe
  if(!inherits(archival, "data.frame")) stop("'archival' must be a data.frame")
  check_names(input = archival, req = "depth", extract_names = colnames, type = all)
  if(any(is.na(archival$depth))) stop("'archival$depth' contains NAs.")
  # Check calc_depth_error
  de_1 <- calc_depth_error(archival$depth[1])
  if(length(de_1) != 2){
    stop("'calc_depth_error' should be a function that returns a numeric vector of length two (i.e., a lower and upper depth adjustment).")
  }
  if(de_1[1] > 0 | de_1[2] < 0){
    stop("'calc_depth_error' should return a negative and a postive adjustment (in that order).")
  }
  # Check cluster
  if(is.null(split) & !is.null(cl)) {
    warning("'cl' argument ignored unless 'split' is supplied.",
            immediate. = TRUE, call. = TRUE)
    cl <- NULL
    varlist <- NULL
  }
  if(is.null(cl) & !is.null(varlist)) {
    warning("'varlist' is supplied but 'cl' is NULL.",
            immediate. = TRUE, call. = TRUE)
    varlist <- NULL
  }
  # Define an index for saving/writing files
  drop_index_col <- TRUE
  if(rlang::has_name(archival, "index")) {
    drop_index_col <- FALSE
    warning("'archival$index overwritten.", call. = FALSE)
  }
  archival$index <- 1:nrow(archival)
  if(is.null(plot)) plot <- archival$index
  # Check write_history inputs
  if(!is.null(write_history)){
    # Check directory
    check_named_list(input = write_history)
    check_names(input = write_history, req = "filename")
    write_history$filename <- check_dir(input = write_history$filename, check_slash = TRUE)
    write_history_dir <- write_history$filename
  }

  #### Implement calc_depth_error()
  cat_to_console("... Implementing calc_depth_error()...")
  archival$depth_lwr <- archival$depth + calc_depth_error(archival$depth)[1]
  archival$depth_upr <- archival$depth + calc_depth_error(archival$depth)[2]

  #### Implement algorithm

  # Define a blank map
  blank <- bathy
  blank <- raster::setValues(blank, 0)
  # Define 'availability'
  if(check_availability) archival$availability <- NA

  # Implement algorithm
  if(is.null(split)){
    cat_to_console("... Implementing algorithm over time steps...")
    use <- blank
    if(verbose) pb <- utils::txtProgressBar(min = 0, max = nrow(archival), style = 3)
    for(i in 1:nrow(archival)){
      avail <- bathy >= archival$depth_lwr[i] & bathy <= archival$depth_upr[i]
      if(!is.null(write_history)){
        write_history$x <- avail
        write_history$filename <- paste0(write_history_dir, "arc_", archival$index[i])
        do.call(raster::writeRaster, write_history)
      }
      if(check_availability) archival$availability[i] <- raster::maxValue(avail) == 1
      use <- use + avail
      if(plot != 0) {
        if(archival$index[i] %in% plot){
          out$spatial[[i]] <- list(map_timestep = avail, map_cumulative = use)
        }
      }
      if(verbose) utils::setTxtProgressBar(pb, i)
    }
    if(verbose) close(pb)

  } else {
    cat_to_console("... Implementing algorithm over chunks...")
    archival_ls <- split(archival, rep(1:ceiling(nrow(archival)/split), each = split)[1:nrow(archival)])
    if(!is.null(cl) & !is.null(varlist)) parallel::clusterExport(cl = cl, varlist = varlist)
    use_and_avail_by_chunk <- pbapply::pblapply(archival_ls, cl = cl, function(d){
      ret <- list(spatial = NULL)
      blank <- raster::setValues(bathy, 0)
      use <- blank
      for(i in 1:nrow(d)){
        avail <- bathy >= d$depth_lwr[i] & bathy <= d$depth_upr[i]
        if(!is.null(write_history)){
          write_history$x <- avail
          write_history$filename <- paste0(write_history_dir, "arc_", d$index[i])
          do.call(raster::writeRaster, write_history)
        }
        if(check_availability) d$availability[i] <- raster::maxValue(avail) == 1
        use <- use + avail
        if(plot != 0) {
          if(d$index[i] %in% plot){
            ret$spatial[[i]] <- list(map_timestep = avail, map_cumulative = use)
          }
        }
      }
      ret$avail <- d
      ret$use   <- use
      return(ret)
    })
    if(!is.null(cl)) parallel::stopCluster(cl)
    use_by_chunk     <- lapply(use_and_avail_by_chunk, function(elm) elm$use)
    avail_by_chunk   <- lapply(use_and_avail_by_chunk, function(elm) elm$avail)
    if(plot != 0) {
      spatial_by_chunk <- lapply(use_and_avail_by_chunk, function(elm) elm$spatial)
      out$spatial      <- purrr::flatten(spatial_by_chunk)
      names(out$spatial) <- as.character(plot)
    }
    archival <- do.call(rbind, avail_by_chunk)
    if(drop_index_col) archival$index <- NULL
    cat_to_console("... Stacking chunk usage maps...")
    use <- raster::stack(use_by_chunk)
    cat_to_console("... Summing chunk usage maps to generate a single map of space use...")
    use <- sum(use)
  }

  #### Return outputs
  out$archival <- archival
  out$dc <- use
  t_end <- Sys.time()
  total_duration <- difftime(t_end, t_onset, units = "mins")
  cat_to_console(paste0("... flapper::dc() call completed (@ ", t_end, ") after ~", round(total_duration, digits = 2), " minutes."))
  return(out)

}




