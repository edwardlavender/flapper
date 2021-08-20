#' @title The depth-contour (DC) algorithm
#' @description This function implements the depth-contour (DC) algorithm. Under the assumption that individuals are benthic/demersal, this algorithm relates one-dimensional depth time series to a two-dimensional bathymetry surface to determine the extent to which different parts of an area might have (or have not) been used, or effectively represent occupied depths, over time. Given a sequence of depth observations (\code{archival}) from a benthic animal and a measurement error parameter (\code{calc_depth_error}), at each time step the function determines the cells on a bathymetry \code{\link[raster]{raster}} (\code{bathy}) that match the observed depth*. Across all time steps, matches are summed to produce a single map representing the number of occasions when the depth in each cell matched the observed depth.
#'
#' @param archival A dataframe of depth time series (for a single individual). At a minimum, this should contain a column named `depth' with depth observations. Depth should be recorded using absolute values in the same units as the bathymetry (\code{bathy}, see below).
#' @param bathy A \code{\link[raster]{raster}} of the bathymetry in an area within which the animal is likely to have been located over the study. Bathymetry values should be recorded as absolute values and in the same units as for depths (see \code{archival}).
#' @param plot_ts A logical input that defines whether or not to the depth time series before the algorithm is initiated.
#' @param calc_depth_error A function that returns the depth errors around a vector of depths. The function should accept vector of depths (from \code{archival$depth}) and return a matrix, with one row for each (lower and upper) error and one one column for each depth (if the error varies with depth). For each depth, the two numbers are added to the observed depth to define the range of depths on the bathymetry raster (\code{bathy}) that the individual could plausibly have occupied at any time. Since the depth errors are added to the individual's depth, the first number should be negative (i.e., the individual could have been slightly shallower that observed) and the second positive (i.e., the individual could have been slightly deeper than observed). For example, the constant function \code{calc_depth_error = function(...) matrix(c(-2.5, 2.5), nrow = 2)} implies that the individual could have occupied bathymetric cells whose depth lies within the interval defined by the observed depth + (-2.5) and + (+2.5) m. The appropriate form for \code{calc_depth_error} depends on measurement error for the depth observations in \code{archival} and bathymetry (\code{bathy}) data, as well as the tidal range (m) across the area (over the duration of observations), but this implementation allows the depth error to depend on depth and for the lower and upper error around an observation to differ.
#' @param check_availability A logical input that defines whether or not to record explicitly, for each time step, whether or not there were any cells on \code{bathy} that matched the observed depth (within the bounds defined by \code{calc_depth_error}).
#' @param normalise A logical variable that defines whether or not to normalise the map of possible locations at each time step so that their scores sum to one.
#' @param save_record_spatial An integer vector that defines the time steps for which to return time step-specific and cumulative maps of the individual's possible locations. \code{save_record_spatial = 0} suppresses the return of this information and \code{save_record_spatial = NULL} returns this information for all time steps.
#' @param write_record_spatial_for_pf (optional) A named list, passed to \code{\link[raster]{writeRaster}}, to save the \code{\link[raster]{raster}} of the individual's possible positions at each time step to file. The `filename' argument should be the directory in which to save files. Files are named by archival time steps as `arc_1', `arc_2' and so on.
#' @param save_args A logical input that defines whether or not to save the list of function inputs in the returned object.
#' @param verbose A logical variable that defines whether or not to print messages to the console or to file to relay function progress. If \code{con = ""}, messages are printed to the console; otherwise, they are written to file (see below).
#' @param con If \code{verbose = TRUE}, \code{con} is character string that defines how messages relaying function progress are returned. If \code{con = ""}, messages are printed to the console (unless redirected by \code{\link[base]{sink}}). Otherwise, \code{con} defines the full pathway to a .txt file (which can be created on-the-fly) into which messages are written to relay function progress.
#' @param split,cl,varlist (optional) Parallelisation arguments. \code{split} is an integer which, if supplied, splits the \code{archival} dataframe every \code{n}th row into chunks†. The algorithm is applied sequentially within each chunk (if applicable) and chunk-wise maps can be summed afterwards to create a single map of space use. The advantage of this approach is that chunks can be analysed in parallel, via \code{cl} and \code{varlist}, while memory use is minimised. \code{cl} is a cluster object created by \code{\link[parallel]{makeCluster}} to implement the algorithm in parallel. \code{varlist} is a character vector of names of objects to export, to be passed to the \code{varlist} argument of \code{\link[parallel]{clusterExport}}. This may be required if \code{cl} is supplied. Exported objects must be located in the global environment.

#' @param ... Additional arguments (none implemented).
#'
#' @details
#'
#' *Location probability could be weighted by the proximity between the observed depth and the depths of cells within the range defined by the observed depth and the measurement error (i.e., \code{archival$depth[1] + calc_depth_error(archival$depth[1])[1], archival$depth[1] + calc_depth_error(archival$depth[1])[2]}), but this is not currently implemented.
#'
#' †Under the default options (\code{split = NULL}), the function starts with a blank map of the area and iterates over each time step, adding the `possible positions' of the individual to the map at each step. By continuously updating a single map, this approach is slow but minimises memory requirements. An alternative approach is to split the time series into chunks, implement an iterative approach within each chunk, and then join the maps for each chunk. This is implemented by \code{split} (plus \code{\link[flapper]{acdc_simplify}}).
#'
#' @return  The function returns a \code{\link[flapper]{acdc-class}} object. If a connection to write files has also been specified, messages are also written to file relaying function progress.
#'
#'
#' @examples
#' #### Define depth time series for examples
#' # We will use a sample depth time series for one individual
#' # We will select a small sample of observations for example speed
#' depth <- dat_archival[dat_archival$individual_id == 25, ][1:100, ]
#'
#' #### Example (1): Implement algorithm with default options
#' ## Implement algorithm
#' out_dc <- dc(archival = depth, bathy = dat_gebco)
#' ## dc() returns a named list that is an 'acdc' class object
#' summary(out_dc)
#' # ... .acdc elements contain the key information from the algorithm
#' # ... ts_by_chunk contains the time series for each chunk
#' # ... time contains a dataframe of the algorithm's progression
#' # ... args contains a list of user inputs
#' ## Simplify outputs via acdc_simplify()
#' dc_summary <- acdc_simplify(out_dc, type = "dc")
#' summary(dc_summary)
#' ## Examine time-specific maps via acdc_plot_record()
#' acdc_plot_record(dc_summary)
#' ## Examine overall map via a raster* plotting function
#' # Each cell shows the number of time steps when the bathymetry data
#' # ... in each cell matches the archival data
#' prettyGraphics::pretty_map(add_rasters = list(x = dc_summary$map),
#'                            add_polys = list(x = dat_coast))
#' # Convert counts on map to percentages
#' prettyGraphics::pretty_map(add_rasters =
#'                              list(x = dc_summary$map/nrow(depth) * 100,
#'                                   zlim = c(0, 100)),
#'                            add_polys = list(x = dat_coast))
#' # Check for occasions when the individual's depth was not consistent
#' # ... with the depth data for the area and the depth error e.g., possibly
#' # ... due to movement beyond this area:
#' any(do.call(rbind, lapply(dc_summary$record, function(elm) elm$dat))$availability == FALSE)
#'
#' #### Example (2): Implement depth error functions that depend on depth
#' ## Here, we will define a calc_depth_error function that depends on:
#' # ... tag error
#' # ... tidal range
#' # ... a depth-dependent bathymetry error
#' cde <- function(depth){
#'   e <- 4.77 + 2.5 + sqrt(0.5 ^2 + (0.013 * depth)^2)
#'   e <- matrix(c(-e, e), nrow = 2)
#'   return(e)
#'   }
#' # Vectorise function over depths (essential)
#' cde <- Vectorize(cde)
#' ## Implement algorithm with  depth-dependent error
#' out_dc <- dc(archival = depth, bathy = dat_gebco, calc_depth_error = cde)
#'
#' #### Example (3): Write timestep-specific maps of allowed positions to file
#' out_dc <- dc(archival = depth, bathy = dat_gebco,
#'              write_record_spatial_for_pf = list(filename = tempdir()))
#' list.files(tempdir())
#'
#' #### Example (4): Implement the algorithm in parallel
#' ## Trial different options for 'split' and compare speed
#' # Approach (1)
#' at1 <- Sys.time()
#' cl <- parallel::makeCluster(2L)
#' parallel::clusterEvalQ(cl = cl, library(raster))
#' out_dc <- dc(archival = depth,
#'              bathy = dat_gebco,
#'              plot_ts = FALSE,
#'              split = 1,
#'              cl = cl)
#' at2 <- Sys.time()
#' # Approach (2)
#' bt1 <- Sys.time()
#' cl <- parallel::makeCluster(2L)
#' parallel::clusterEvalQ(cl = cl, library(raster))
#' out_dc <- dc(archival = depth,
#'              bathy = dat_gebco,
#'              plot_ts = FALSE,
#'              split = 5,
#'              cl = cl)
#' bt2 <- Sys.time()
#' # Compare timings
#' difftime(at2, at1)
#' difftime(bt2, bt1)
#'
#' #### Example (5): Write messages to file via con
#' out_dc <- dc(archival = depth, bathy = dat_gebco,
#'              con = paste0(tempdir(), "/dc_log.txt"))
#' readLines(paste0(tempdir(), "/dc_log.txt"))
#'
#' #### Example (6) Compare an automated vs. manual implementation of DC
#'
#' ## (A) Compare step-wise results for randomly select time steps
#' # Implement algorithm (using default calc_depth_error)
#' out_dc <- dc(archival = depth,
#'              bathy = dat_gebco,
#'              save_record_spatial = NULL)
#' # Compare results for randomly selected time steps to manual implementation
#' # ... with the same simple calc_depth_error model
#' pp <- graphics::par(mfrow = c(5, 2))
#' for(i in sample(1:nrow(depth), 5)){
#'   # Extract map created via dc()
#'   raster::plot(out_dc$.acdc$record[[i]]$spatial[[1]]$map_timestep)
#'   # Compare to manual creation of map
#'   raster::plot(dat_gebco >= depth$depth[i] - 2.5 & dat_gebco <= depth$depth[i] + 2.5)
#' }
#' graphics::par(pp)
#'
#' ## (B) Compare the chunk-wise results and manual implementation
#' # Implement algorithm chunk-wise
#' out_dc <- dc(archival = depth,
#'              bathy = dat_gebco,
#'              save_record_spatial = NULL,
#'              split = 10L)
#' # Get overall map via acdc_simplify()
#' map_1 <- acdc_simplify(out_dc, type = "dc", mask = dat_gebco)$map
#' # Get overall map via a simple custom approach
#' map_2 <- raster::setValues(dat_gebco, 0)
#' for(i in 1:nrow(depth)){
#'   map_2 <- map_2 + (dat_gebco >= depth$depth[i] - 2.5 & dat_gebco <= depth$depth[i] + 2.5)
#' }
#' map_2 <- raster::mask(map_2, dat_gebco)
#' # Show that the two overall maps are identical
#' pp <- par(mfrow = c(1, 3))
#' raster::plot(map_1, main = "Automated")
#' raster::plot(map_2, main = "Custom")
#' raster::plot(map_1 - map_2, main = "Difference")
#' par(pp)
#'
#' @seealso \code{\link[flapper]{acdc_simplify}} simplifies the outputs of the algorithm. \code{\link[flapper]{acdc_plot_record}} and \code{\link[flapper]{acdc_animate_record}} provide plotting routines. \code{\link[flapper]{dcq}} implements a faster version of this algorithm termed the `quick depth-contour' (DCQ) algorithm. Rather than considering the depth interval that the individual could have occupied at each time step, the DCQ algorithm considers a sequence of depth bins (e.g., 10 m bins), isolates these on the bathymetry \code{\link[raster]{raster}} (\code{bathy}) and counts the number of matches in each cell. The DCPF algorithm (see \code{\link[flapper]{pf}}) extends the DC algorithm via particle filtering to reconstruct possible movement paths over \code{bathy}. The ACDC algorithm (see \code{\link[flapper]{acdc}}) extends the depth-contour algorithm by integrating information from acoustic detections of individuals at each time step to restrict the locations in which depth contours are identified.
#'
#' @author Edward Lavender
#' @export

dc <- function(archival,
               bathy,
               plot_ts = TRUE,
               calc_depth_error = function(...) matrix(c(-2.5, 2.5), nrow = 2),
               check_availability = TRUE,
               normalise = FALSE,
               save_record_spatial = 1L,
               write_record_spatial_for_pf = NULL,
               save_args = TRUE,
               verbose = TRUE,
               con = "",
               split = NULL,
               cl = NULL,
               varlist = NULL
               ){

  #### Set up function
  t_onset <- Sys.time()
  if(con != ""){
    if(!verbose) {
      message("Input to 'con' ignored since verbose = FALSE.")
    } else {
      # Check directory
      check_dir(input = dirname(con))
      # Write black file to directory if required
      if(!file.exists(con)){
        message(paste0(con, " does not exist: attempting to write file in specified directory..."))
        file.create(file1 = con)
        message("... Blank file successfully written to file.")
      }
    }
  }
  append_messages <- ifelse(con == "", FALSE, TRUE)
  cat_to_cf <- function(..., message = verbose, file = con, append = append_messages){
    if(message) cat(paste(..., "\n"), file = con, append = append)
  }

  #### Initiate function with storage container for outputs
  cat_to_cf(paste0("flapper::dc() called (@ ", t_onset, ")..."))
  cat_to_cf("... Setting up function...")
  out <- list(.acdc = NULL, ts_by_chunk = NULL, time = NULL, args = NULL)
  out$time <- data.frame(event = "onset", time = t_onset)
  if(save_args){
    out$args <- list(archival = archival,
                     bathy = bathy,
                     plot_ts = plot_ts,
                     calc_depth_error = calc_depth_error,
                     check_availability = check_availability,
                     normalise = normalise,
                     save_record_spatial = save_record_spatial,
                     write_record_spatial_for_pf = write_record_spatial_for_pf,
                     save_args = save_args,
                     verbose = verbose,
                     con = con,
                     split = split,
                     cl = cl,
                     varlist = varlist
                     )
    }

  #### Checks
  ## Check archival dataframe
  if(!inherits(archival, "data.frame")) stop("'archival' must be a data.frame")
  check_names(input = archival, req = "depth", extract_names = colnames, type = all)
  if(any(is.na(archival$depth))) stop("'archival$depth' contains NAs.")
  ## Check calc_depth_error
  de <- calc_depth_error(archival$depth)
  if(inherits(de, "matrix")){
    if(nrow(de) == 2){
      if(ncol(de) == 1) { message("'calc_depth_error' function taken to be independent of depth.")
      } else {
        message("'calc_depth_error' taken to depend on depth.")
      }
      if(any(de[1, ] > 0) | any(de[2, ] < 0)) stop("'calc_depth_error' should be a function that returns a two-row matrix with lower (negative) adjustment(s) (top row) and upper (positive) adjustment(s) (bottom row).'", call. = FALSE)
    } else stop("'calc_depth_error' should return a two-row matrix.", call. = FALSE)
  } else stop("'calc_depth_error' should return a two-row matrix.", call. = FALSE)
  ## Check cluster
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
  ## Check write_record_spatial_for_pf inputs
  if(!is.null(write_record_spatial_for_pf)){
    check_named_list(input = write_record_spatial_for_pf)
    check_names(input = write_record_spatial_for_pf, req = "filename")
    write_record_spatial_for_pf$filename <- check_dir(input = write_record_spatial_for_pf$filename, check_slash = TRUE)
    write_record_spatial_for_pf_dir <- write_record_spatial_for_pf$filename
  }

  #### Set up objects for algorithm
  ## Archival index (for saving the spatial record)
  if(rlang::has_name(archival, "index")) {
    warning("'archival$index overwritten.", immediate. = TRUE, call. = FALSE)
  }
  archival$index <- 1:nrow(archival)
  if(is.null(save_record_spatial)) save_record_spatial <- archival$index
  ## Archival time series with depth error
  cat_to_cf("... Implementing calc_depth_error()...")
  archival$depth_lwr <- archival$depth + calc_depth_error(archival$depth)[1, ]
  archival$depth_upr <- archival$depth + calc_depth_error(archival$depth)[2, ]
  ## Archival time series with availability
  if(check_availability) archival$availability <- NA
  ## Archival time series as list (by chunk)
  if(is.null(split)){
    archival_ls <- list(archival)
  } else {
    archival_ls <- split(archival, rep(1:ceiling(nrow(archival)/split), each = split)[1:nrow(archival)])
  }
  out$ts_by_chunk <- lapply(archival_ls, function(d) list(acoustics = NULL, archival = d))
  ## Blank map for updating
  blank <- raster::setValues(bathy, 0)

  #### Plot time series (for each chunk)
  if(plot_ts){
    cat_to_cf("... Plotting movement time series (for each chunk)...")
    if(length(archival_ls) < 25) pp <- graphics::par(mfrow = prettyGraphics::par_mf(length(archival_ls)))
    lapply(archival_ls, function(d){
      if(nrow(d) == 1) pt <- "p" else pt <- "l"
      prettyGraphics::pretty_plot(1:nrow(d), abs(d$depth) * - 1,
                                  pretty_axis_args = list(side = 3:2),
                                  xlab = "Time (index)", ylab = "Depth (m)",
                                  type = pt)
      graphics::lines(1:nrow(d), abs(d$depth_lwr) * -1, lty = 3, col = "royalblue", type = pt)
      graphics::lines(1:nrow(d), abs(d$depth_upr) * -1, lty = 3, col = "royalblue", type = pt)
    })
    if(length(archival_ls) < 25) graphics::par(pp)
  }

  #### Implement algorithm
  ## Initiation
  out$time <- rbind(out$time, data.frame(event = "algorithm_initiation", time = Sys.time()))
  if(is.null(split)){
    cat_to_cf("... Implementing algorithm over time steps...")
  } else {
    cat_to_cf("... Implementing algorithm over chunks...")
  }
  if(!is.null(cl) & !is.null(varlist)) parallel::clusterExport(cl = cl, varlist = varlist)
  ## Implement .dc() algorithm over chunks
  .out_by_chunk <- pbapply::pblapply(archival_ls, cl = cl, function(d){
    # Define storage container for chunk-specific outputs
    .out <- list(map = NULL,
                 record = list(),
                 time = NULL,
                 args = NULL,
                 chunks = NULL,
                 simplify = FALSE)
    use <- blank
    for(i in 1:nrow(d)){
      avail <- bathy >= d$depth_lwr[i] & bathy <= d$depth_upr[i]
      if(check_availability) d$availability[i] <- raster::maxValue(avail) == 1
      if(normalise) avail <- avail/raster::cellStats(avail, "sum")
      if(!is.null(write_record_spatial_for_pf)){
        write_record_spatial_for_pf$x <- avail
        write_record_spatial_for_pf$filename <- paste0(write_record_spatial_for_pf_dir, "arc_", d$index[i])
        do.call(raster::writeRaster, write_record_spatial_for_pf)
      }
      use <- use + avail
      .out$record[[i]] <- list(dat = d[i, , drop = FALSE],
                               spatial = list())
      if(save_record_spatial[1] != 0) {
        if(d$index[i] %in% save_record_spatial){
          .out$record[[i]]$spatial[[1]] <- list(map_timestep = avail, map_cumulative = use)
        }
      }
    }
    # Update and return outputs
    .out$map <- use
    class(.out) <- c(class(.out), ".acdc")
    return(.out)
  })
  if(!is.null(cl)) parallel::stopCluster(cl)

  #### Finalise algorithm
  if(length(archival_ls) == 1){
    out$.acdc <- .out_by_chunk[[1]]
  } else {
    out$.acdc <- .out_by_chunk
  }
  t_end <- Sys.time()
  out$time <- rbind(out$time, data.frame(event = "algorithm_competion", time = t_end))
  out$time$serial_duration <- Tools4ETS::serial_difference(out$time$time, units = "mins")
  out$time$total_duration <- NA
  total_duration <- sum(as.numeric(out$time$serial_duration), na.rm = TRUE)
  out$time$total_duration[nrow(out$time)] <- total_duration
  class(out) <- c(class(out), "acdc")
  cat_to_cf(paste0("... flapper::dc() call completed (@ ", t_end, ") after ~", round(total_duration, digits = 2), " minutes."))
  return(out)
}
