########################################
########################################
#### pf_plot_history()

#' @title Plot particle histories from a PF algorithm
#' @description This function plots the spatiotemporal particle histories from a particle filtering (PF) algorithm (the acoustic-centroid PF, the depth-contour PF or the acoustic-centroid depth-contour PF). This produces, for each time step, a map of the individual's possible locations (from the AC, DC or ACDC algorithm), with sampled locations (derived via the particle filtering routine) overlaid.
#' @param archive A \code{\link[flapper]{pf_archive-class}} object from \code{\link[flapper]{pf}} that contains particle histories.
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

pf_plot_history <- function(archive,
                            time_steps = 1:length(history),
                            add_surface = list(),
                            add_particles = list(pch = "."),
                            forwards = TRUE,
                            prompt = TRUE,...){
  if(!inherits(archive, "pf_archive")) stop("'archive' must be a 'pf_archive' class object.")
  layers           <- archive$args$record
  history          <- archive$history
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
#' @description This function creates a plot of the `probability of use' across an area based on particles sampled by a particle filtering (PF) algorithm. To implement the function, a \code{\link[flapper]{pf_archive-class}} object that contains particles (locations) sampled by \code{\link[flapper]{pf}} must be supplied. The function extracts all sampled locations and, for each location, calculates `the probability of use' for that location over the time series. This is returned (invisibly) as a \code{\link[raster]{raster}} and plotted.
#' @param archive A \code{\link[flapper]{pf_archive-class}} object (from \code{\link[flapper]{pf}}).
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

pf_plot_map <- function(archive,
                        map,
                        scale = c("original", "max", "sum"),
                        add_rasters = list(),...){
  # Check inputs
  check_class(input = archive, to_class = "pf_archive")
  scale <- match.arg(scale)
  # Extract particle histories as a single dataframe
  pf_particle_histories <- lapply(archive$history, function(elm) elm[, c("id_current", "pr_current"), drop = FALSE])
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
