########################################
########################################
#### pf_plot_history()

#' @title Plot particle histories from a PF algorithm
#' @description This function plots the spatiotemporal particle histories from a particle filtering (PF) algorithm (the acoustic-centroid PF, the depth-contour PF or the acoustic-centroid depth-contour PF). This produces, for each time step, a map of the individual's possible locations (from the AC, DC or ACDC algorithm), with sampled locations (derived via the particle filtering routine) overlaid.
#' @param archive A \code{\link[flapper]{pf_archive-class}} object from \code{\link[flapper]{pf}}, or \code{\link[flapper]{pf}} plus \code{\link[flapper]{pf_simplify}} with the \code{record = "archive"} argument, that contains particle histories.
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
#' #### Example (4): Compare outputs for sampled versus connected particles
#' dat_dcpf_histories_connected <-
#'   pf_simplify(dat_dcpf_histories, return = "archive")
#' pp <- graphics::par(mfcol = c(2, 4))
#' pf_plot_history(dat_dcpf_histories, time_steps = 1:4,
#'                 add_particles = list(pch = 21, bg = "black"),
#'                 prompt  = FALSE)
#' pf_plot_history(dat_dcpf_histories_connected, time_steps = 1:4,
#'                 add_particles = list(pch = 21, bg = "black"),
#'                 prompt = FALSE)
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
#### pf_animate_history()

#' @title Create a html animation of the PF algorithm(s)
#' @description This function is a simple wrapper for \code{\link[flapper]{pf_plot_history}} and \code{\link[animation]{saveHTML}} which creates an animation of the particle filtering (PF) algorithm(s) over time. To implement this function, a named list of arguments for \code{\link[flapper]{pf_plot_history}}, which creates the plots, must be supplied. This is embedded within \code{\link[animation]{saveHTML}}, which creates a folder in the specified directory named `images' that contains a .png file for each time step and an animation as a .html file.
#' @param expr_param A named list of arguments, passed to \code{\link[flapper]{pf_plot_history}}, to create plots.
#' @param dir (optional) A string that defines the directory in which to save files. If unsupplied, if available, \code{dir} is taken from \code{html_name} using \code{\link[base]{dirname}}.
#' @param html_name A string that defines the name of the html file (see `htmlfile' argument in \code{\link[animation]{saveHTML}}).
#' @param image_name A string that defines the names of the individual .png files creates (see `img.name' argument in \code{\link[animation]{saveHTML}}).
#' @param html_title,html_description Character strings that provide a title and a description that are displayed within the html animation (see `title' and `description' arguments in \code{\link[animation]{saveHTML}}).
#' @param navigator A logical variable that defines whether or not to add a navigator panel to the animation (see `navigator' argument in \code{\link[animation]{saveHTML}}).
#' @param ani_height,ani_width,ani_res Numbers that define the size and the resolution of the animation (see `ani.height' `ani.width' and `ani.res' arguments in \code{\link[animation]{ani.options}}).
#' @param interval A number that defines the time interval between sequential frames (see `interval' argument in \code{\link[animation]{ani.options}}).
#' @param verbose A logical or character variable that defines whether or not, or what, to write as a footer to the html animation (see `verbose' argument in \code{\link[animation]{ani.options}}).
#' @param ... Additional arguments passed to \code{\link[animation]{ani.options}}.
#'
#' @return The function produces an animation in .html format in the specified directory. A folder named `images' is also produced which contains the images for each time step. The `css' and `js' folders are also produced by \code{\link[animation]{saveHTML}} which creates the animation.
#'
#' @examples
#' #### Example (1): Create a zoomed-in animation
#' pf_animate_history(
#'   expr_param = list(archive = dat_dcpf_histories,
#'                     add_particles = list(cex = 2.5, pch = 21,
#'                                          col = "black", bg = "black"),
#'                     prompt = FALSE),
#'   dir = tempdir(),
#'   interval = 0.25)
#'
#' #### Example (2): Create a wider scale animation
#' boundaries <- raster::extent(dat_coast)
#' pf_animate_history(
#'   expr_param = list(archive = dat_dcpf_histories,
#'                     add_particles = list(cex = 0.5, pch = 21,
#'                                           col = "black", bg = "black"),
#'                     add_polys = list(x = dat_coast, col = "brown"),
#'                     xlim = boundaries[1:2], ylim = boundaries[3:4],
#'                     prompt = FALSE),
#'   dir = tempdir())
#'
#' @details This function requires the \code{\link[animation]{animation}} package.
#' @author Edward Lavender
#' @export
#'

pf_animate_history <-
  function(expr_param,
           dir = NULL,
           html_name = "PF_algorithm_demo.html",
           image_name = "PF",
           html_title = "Demonstration of PF",
           html_description = "",
           navigator = FALSE,
           ani_height = 800,
           ani_width = 800,
           ani_res = 1200,
           interval = 0.1,
           verbose = FALSE,
           ...){
    #### Checks
    ## animation package
    if (!requireNamespace("animation", quietly = TRUE)) {
      stop("This function requires the 'animation' package. Please install it before continuing with install.packages('animation').")
    }
    #### Set directory
    if(is.null(dir)) dir <- dirname(html_name)
    wd <- getwd()
    check_dir(input = dir)
    setwd(dir)
    html_name <- basename(html_name)
    on.exit(setwd(wd), add = TRUE)
    #### Make plot
    animation::saveHTML({
      do.call(pf_plot_history, expr_param)
    },
    htmlfile = html_name,
    img.name = image_name,
    title = html_title,
    description = html_description,
    navigator = navigator,
    ani.height = ani_height,
    ani.width = ani_width,
    ani.res = ani_res,
    interval = interval,
    verbose = verbose,...
    )
    return(invisible())
  }


######################################
######################################
#### pf_plot_map()

#' @title Plot `probability of use' from a PF algorithm
#' @description This function creates a plot of the `probability of use' across an area based on particles sampled by a particle filtering (PF) algorithm. To implement the function, a \code{\link[flapper]{pf_archive-class}} object that contains connected particles (locations) sampled by \code{\link[flapper]{pf}} and processed by \code{\link[flapper]{pf_simplify}} must be supplied. The function extracts sampled locations and, for each location, calculates `the probability of use' for that location over the time series. This is returned (invisibly) as a \code{\link[raster]{raster}} and plotted.
#' @param archive A \code{\link[flapper]{pf_archive-class}} object (from \code{\link[flapper]{pf}} plus \code{\link[flapper]{pf_simplify}} with \code{record = "archive"}).
#' @param map A \code{\link[raster]{raster}} that defines a grid across the area of interest.
#' @param transform (optional) A function to transform cell weights (e.g, \code{\link[base]{log}}).
#' @param scale A character that defines how \code{\link[raster]{raster}} values are scaled: \code{"original"} uses the original values; \code{"max"} scales values by the maximum value (so that, if \code{transform = NULL}, they lie between zero and one; and \code{"sum"} scales values by their sum so that they sum to one.
#' @param add_rasters A named list, passed to \code{\link[prettyGraphics]{pretty_map}}, to customise the appearance of the plotted surface.
#' @param ... Additional arguments passed to \code{\link[prettyGraphics]{pretty_map}}.
#'
#' @details For each location, the 'probability of use' is calculated as the sum of the number of times (time steps) that the location was sampled, weighted by the associated probabilities of each sample.
#'
#' @return The function invisibly returns a \code{\link[raster]{raster}}, in which each cell contains the `probability of use' score and produces a plot of this surface.
#'
#' @examples
#' #### Prepare data
#' # The example data 'dat_dcpf_histories' contains all particles sampled
#' # ... by an implementation of the DCPF algorithm. However, not all particles
#' # ... that were sampled at one time step may have been near to particles sampled
#' # ... at the next time step. In addition, some particles may have been sampled
#' # ... multiple times at one time step, but our maps of space use should reflect
#' # ... the number of time steps that the individual could have occupied a location,
#' # ... rather than the total number of samples of a location. Hence, to map
#' # ... space use, we should focus on the subset of particles that were connected
#' # ... between time steps and only retain one record of each particle at each time step
#' # ... using pf_simplify() with record = "archive"
#' dat_dcpf_histories_connected <- pf_simplify(dat_dcpf_histories, return = "archive")
#'
#' #### Example (1): Implement the function with default options
#' pf_plot_map(dat_dcpf_histories_connected, map = dat_dc$args$bathy)
#'
#' #### Example (2): Re-scale the map
#' pf_plot_map(dat_dcpf_histories_connected, map = dat_dc$args$bathy, scale = "max")
#' pf_plot_map(dat_dcpf_histories_connected, map = dat_dc$args$bathy, scale = "sum")
#'
#' #### Example (3): Customise the map
#' pf_plot_map(dat_dcpf_histories_connected, map = dat_dc$args$bathy,
#'             add_rasters = list(col = grDevices::grey.colors(n = 100)),
#'             xlab = "x", ylab = "y")
#'
#' @seealso \code{\link[flapper]{pf}} implements PF. \code{\link[flapper]{pf_simplify}} assembles paths from particle histories. \code{\link[flapper]{pf_plot_history}} visualises particle histories. \code{\link[flapper]{pf_plot_1d}}, \code{\link[flapper]{pf_plot_2d}} and \code{\link[flapper]{pf_plot_3d}} provide plotting routines for paths. \code{\link[flapper]{pf_loglik}} calculates the log-likelihood of each path.
#' @author Edward Lavender
#' @export

pf_plot_map <- function(archive,
                        map,
                        transform = NULL,
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
  # Transform weighted frequencies
  if(!is.null(transform)) {
    wt_freq$wt <- transform(wt_freq$wt)
    if(any(is.na(wt_freq$wt))) stop("'transform' function has created NAs.")
  }
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
