######################################
######################################
#### acdc_simplify()

#' @title Simplify the outputs of the AC/DC algorithms
#' @description This function simplifies the output of \code{\link[flapper]{ac}}, \code{\link[flapper]{dc}} and \code{\link[flapper]{acdc}}, by processing information from the 'archive' elements of an \code{\link[flapper]{acdc_archive-class}} object that hold the results of calls to the workhorse routines. This is especially useful if the algorithm(s) have been applied chunk-wise, in which case the results for each chunk are returned in a list. The function aggregates information across chunks to generate a continuous time series of results and a map of the expected proportion of time steps spent in each grid cell.
#' @param archive An \code{\link[flapper]{acdc_archive-class}} object returned by \code{\link[flapper]{ac}}, \code{\link[flapper]{dc}} or \code{\link[flapper]{acdc}}.
#' @param type A character that defines whether the function should be implemented for the outputs of a call to an AC* algorithm (\code{\link[flapper]{ac}} or \code{\link[flapper]{acdc}}), in which case \code{type = "acs"}, or the DC algorithm (\code{\link[flapper]{dc}}), in which case \code{type = "dc"}.
#' @param mask (optional) A spatial mask (e.g., the argument passed to \code{bathy} in \code{\link[flapper]{ac}}, \code{\link[flapper]{dc}} or \code{\link[flapper]{acdc}}) to mask areas (e.g., land) from the overall map. If implemented, cells in masked areas are assigned NAs rather than a score of 0.
#' @param normalise A logical input that defines whether or not to normalise the overall map so that cell scores sum to one. If \code{normalise = FALSE}, the overall map represents the expected number of time steps spent in each grid cell; if \code{normalise = TRUE}, the overall map represents the expected proportion of time steps spent in each grid cell.
#' @param keep_chunks A logical variable that defines whether or not to retain all chunk-specific information.
#' @return The function returns an object of class \code{\link[flapper]{acdc_record-class}}.
#' @details If the \code{\link[flapper]{ac}}, \code{\link[flapper]{dc}} or \code{\link[flapper]{acdc}} function was implemented step-wise, this function simply extracts the necessary information and re-packages it into an \code{\link[flapper]{acdc_record-class}} object. For a chunk-wise implementation, the function (a) computes the map of where the individual could have spent more or less time by aggregating the chunk-specific maps) and (b) simplifies chunk-specific records into a single contiguous time series, with re-defined time stamps from the start to the end of the time series (for AC* algorithm(s)) to return an \code{\link[flapper]{acdc_record-class}} object.
#' @seealso The AC, DC and ACDC algorithms are implemented by \code{\link[flapper]{ac}}, \code{\link[flapper]{dc}} and \code{\link[flapper]{acdc}}. After simplification, \code{\link[flapper]{acdc_plot_trace}}, \code{\link[flapper]{acdc_plot_record}} and \code{\link[flapper]{acdc_animate_record}} can be implemented to visualise time-specific results.
#' @author Edward Lavender
#' @export
#'

acdc_simplify <- function(archive,
                          type = c("acs", "dc"),
                          mask = NULL,
                          normalise = TRUE,
                          keep_chunks = FALSE) {
  #### Checks
  if (!(inherits(archive, "acdc_archive") | !inherits(archive, "acdc_record"))) {
    stop("Object of class 'acdc_archive' expected.")
  }
  if (inherits(archive, "acdc_record")) {
    message("class(archive) == 'acdc_record': 'archive' returned unchanged.")
    return(archive)
  }
  type <- match.arg(type)
  message("acdc_simplify() implemented for type = '", type, "'.")

  #### Define container for outputs
  out <- list(map = NULL, record = NULL, time = archive$time, args = archive$args, chunks = NULL, simplify = TRUE)

  #### Keep chunk-specific information (unchanged), if requested
  if (keep_chunks) out$chunks <- archive$archive

  #### Simplify extract outputs the algorithm has only been implemented for a single chunk
  if (length(archive$archive) == 1) {
    out$map <- archive$archive[[1]]$map
    out$record <- archive$archive[[1]]$record

    #### Otherwise aggregate information across chunks
  } else {
    #### Get a list of the cumulative maps from each chunk (to be summed below)
    maps <- lapply(archive$archive, function(chunk) chunk$map)

    #### Process spatial elements so that 'map_cumulative' elements are carried forward (summed) across chunks, if necessary
    try_update_spatial <- TRUE
    if (!raster::inMemory(maps[[1]])) {
      if (!file.exists(maps[[1]]@file@name)) try_update_spatial <- FALSE
    }
    if (!is.null(archive$args)) {
      if (isTRUE(archive$args$save_record_spatial == 0)) try_update_spatial <- FALSE
    }
    if (try_update_spatial) {
      archive$archive <-
        lapply(1:length(archive$archive), function(chunk_id) {
          # chunk_id <- 2
          folder <- archive$archive[[chunk_id]]
          if (chunk_id > 1) {
            if (chunk_id == 2) {
              maps_for_previous_chunks <- maps[[1]]
            } else {
              maps_for_previous_chunks <- maps[1:(chunk_id - 1)]
              maps_for_previous_chunks <- raster::brick(maps_for_previous_chunks)
              maps_for_previous_chunks <- raster::calc(maps_for_previous_chunks, sum, na.rm = TRUE)
            }
            folder$record <-
              lapply(folder$record, function(record_elm) {
                record_elm$spatial <-
                  lapply(record_elm$spatial, function(spatial_elm) {
                    if (rlang::has_name(spatial_elm, "map_cumulative")) {
                      spatial_elm$map_cumulative <- sum(spatial_elm$map_cumulative, maps_for_previous_chunks, na.rm = TRUE)
                    }
                    return(spatial_elm)
                  })
                return(record_elm)
              })
          }
          return(folder)
        })
    }

    #### Simplify records
    out$record <- lapply(archive$archive, function(chunk) chunk$record)

    #### Process record time stamps, if necessary
    if (type == "acs") {
      ## Define a dataframe to adjust the time stamps recorded for each chunks
      # For chunks 2:n_chunks, we will add the time stamps reached by the previous chunk
      # ... up to the current chunk
      adjust_timestep <- lapply(out$record, function(chunk_record) {
        # chunk_record <- out$record[[1]]
        dat <- chunk_record[[length(chunk_record)]]$dat
        adjustment <- dat[nrow(dat), c("timestep_cumulative", "timestep_detection")]
        return(adjustment)
      })
      adjust_timestep <- do.call(rbind, adjust_timestep)
      adjust_timestep$timestep_cumulative <- cumsum(adjust_timestep$timestep_cumulative)
      adjust_timestep$timestep_detection <- cumsum(adjust_timestep$timestep_detection)
      ## Adjust time stamps and add the chunk to the dataframe for each time stamp
      out$record <- lapply(1:length(out$record), function(i) {
        chunk_record <- out$record[[i]]
        if (i == 1) {
          adjustment <- data.frame(timestep_cumulative = 0, timestep_detection = 0)
        } else {
          adjustment <- adjust_timestep[i - 1, ]
        }
        chunk_record <- lapply(chunk_record, function(t) {
          t$dat$timestep_cumulative <- t$dat$timestep_cumulative + adjustment$timestep_cumulative
          t$dat$timestep_detection <- t$dat$timestep_detection + adjustment$timestep_detection
          t$dat$chunk <- i
          return(t)
        })
        return(chunk_record)
      })
    }

    #### Sum chunk-specific maps across chunks
    if (raster::inMemory(maps[[1]]) | file.exists(maps[[1]]@file@name)) {
      out$map <- raster::brick(maps)
      out$map <- raster::calc(out$map, sum, na.rm = TRUE)
    }

    #### Flatten record list across chunks
    out$record <- purrr::flatten(out$record)
  }

  #### Mask and normalise the final map
  if (!is.null(out$map)) {
    if (raster::inMemory(out$map) | file.exists(out$map@file@name)) {
      if (!is.null(mask)) out$map <- raster::mask(out$map, mask)
      if (normalise) out$map <- out$map / raster::cellStats(out$map, "sum")
    }
  }
  if (is.null(out$map)) warning("out$map could not be processed.", call. = FALSE, immediate. = TRUE)

  #### Return outputs
  class(out) <- c(class(out), "acdc_record")
  return(out)
}
