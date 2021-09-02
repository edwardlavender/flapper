########################################
########################################
#### pf_loglik()

#' @title Calculate the log-likelihood of movement paths from a PF algorithm
#' @importFrom rlang .data
#' @description This function calculates the total log-likelihood of each movement path reconstructed by a particle filtering (PF) algorithm, including the acoustic-centroid (AC), depth-contour (DC) or acoustic-centroid depth-contour (ACDC) algorithms.
#' @param paths A dataframe containing movement paths from \code{\link[flapper]{pf}} plus \code{\link[flapper]{pf_simplify}} (see \code{\link[flapper]{pf_path-class}}). At a minimum, this should contain a unique identifier for each path (named `path_id') and the probability associated with each cell along each path (`cell_pr').
#'
#' @details For each path, at each time step the probability associated with the sampled location depends on (a) the `intrinsic' probability associated with each cell (assigned by the AC, DC or ACDC algorithm) and (b) a user-defined movement model that is driven by the distance between the sampled locations for the individual at the previous and current time steps (and other user-defined parameters). This function simply sums the logarithms of these probabilities for each path as a measure of their relative likelihood, given the movement model.
#' @examples
#' # An example with the DCPF paths dataset included in flapper
#' pf_loglik(dat_dcpf_paths)
#' @return The function returns a dataframe with the log likelihood (`loglik') of each path (`path_id'). Rows are ordered by log-likelihood and a `delta' column is provided with the differences in log-likelihood between the most likely path and every other path.
#' @author Edward Lavender
#' @export
#'

pf_loglik <- function(paths){
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
#' @param paths A dataframe containing reconstructed movement path(s) from \code{\link[flapper]{pf}} via \code{\link[flapper]{pf_simplify}} (see \code{\link[flapper]{pf_path-class}}). At a minimum, this should contain a unique identifier for each path (named `path_id’), timesteps (`timestep') and the depth associated with each cell along each path (`cell_z').
#' @param archival A dataframe of depth (m) observations named `depth', as used by \code{\link[flapper]{dc}} and \code{\link[flapper]{acdc}}.
#' @param scale A number that vertically scales the depth time series for the observations and the reconstructed path(s). By default, absolute values for depth are assumed and negated for ease of visualisation.
#' @param pretty_axis_args,xlab,ylab,type,... Plot customisation arguments passed to \code{\link[prettyGraphics]{pretty_plot}}.
#' @param add_lines A named list, passed to \code{\link[graphics]{lines}}, to customise the appearance of the depth time series for reconstructed path(s).
#' @param prompt A logical input that defines whether or not plot the observed depth time series with each reconstructed depth time series on a separate plot, sequentially, with a pause between plots (\code{prompt = TRUE}), or with all reconstructed time series on a single plot (\code{prompt = FALSE}).
#' @details Observed and reconstructed depth time series can differ due to measurement error, which is controlled via the \code{calc_depth_error} function in the DC and ACDC algorithms (see \code{\link[flapper]{dc}} and \code{\link[flapper]{acdc}}).
#' @examples
#' #### Implement pf() algorithm
#' # Here, we use pre-defined outputs for speed
#' paths    <- dat_dcpf_paths
#' archival <- dat_dc$args$archival
#'
#' #### Example (1): The default implementation
#' pf_plot_1d(paths, archival)
#'
#' #### Example (2): Plot customisation options, e.g.:
#' pf_plot_1d(paths, archival, scale = 1, pretty_axis_args = list(side = 1:2))
#' pf_plot_1d(paths, archival, type = "l")
#' pf_plot_1d(paths, archival, add_lines = list(col = "red", lwd = 0.5))
#'
#' #### Example (3): Plot individual comparisons
#' if(interactive()){
#'   pp <- graphics::par(mfrow = c(3, 4))
#'   pf_plot_1d(paths, depth, prompt = TRUE)
#'   graphics::par(pp)
#' }
#'
#' @return The function returns a plot of the observed and reconstructed depth time series, either for all paths at once (if \code{prompt = FALSE}) or each path separately (if \code{prompt = TRUE}).
#' @seealso \code{\link[flapper]{pf}} implements the pf algorithm. \code{\link[flapper]{pf_plot_history}} visualises particle histories, \code{\link[flapper]{pf_plot_map}} creates an overall `probability of use' map from particle histories and \code{\link[flapper]{pf_simplify}} processes the outputs into a dataframe of movement paths. \code{\link[flapper]{pf_plot_1d}}, \code{\link[flapper]{pf_plot_2d}} and \code{\link[flapper]{pf_plot_3d}} provide plotting routines for paths. \code{\link[flapper]{pf_loglik}} calculates the log-likelihood of each path.
#' @author Edward Lavender
#' @export

pf_plot_1d <- function(paths,
                       archival,
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
#' @param paths A dataframe containing reconstructed movement path(s) from \code{\link[flapper]{pf}} via \code{\link[flapper]{pf_simplify}} (see \code{\link[flapper]{pf_path-class}}). At a minimum, this should contain a unique identifier for each path (named `path_id') and the x and y coordinates that define the trajectory of each path (`cell_x' and `cell_y').
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
#' @param paths A dataframe containing reconstructed movement path(s) from \code{\link[flapper]{pf}} via \code{\link[flapper]{pf_simplify}} (see \code{\link[flapper]{pf_path-class}}). At a minimum, this should contain a unique identifier for each path (named `path_id') and the x, y and z coordinates that define the trajectory of each path (`cell_x', `cell_y' and `cell_z').
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
