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


######################################
######################################
#### pf_kud()

#' @title Apply kernel smoothing to particle samples from a PF algorithm
#' @description This function is a wrapper designed to apply kernel utilisation distribution (KUD) estimation to the outputs of a particle filtering (PF) algorithm. To implement the approach, a \code{\link[flapper]{pf_archive-class}} object from \code{\link[flapper]{pf}} (plus \code{\link[flapper]{pf_simplify}} with the \code{return = "archive"} argument) containing particle histories for connected particles must be supplied. Using all, or a subset of these particles, the function applies a KUD smoother to sampled particles via a user-supplied estimation routine (i.e., \code{\link[adehabitatHR]{kernelUD}} or \code{\link[flapper]{kud_around_coastline}}). The function extracts the utilisation distribution as a \code{\link[raster]{raster}}, applies a spatial mask (e.g.  coastline) and plots the distribution, if specified, before returning a \code{\link[raster]{raster}} of the processed KUD.
#'
#' @param archive A \code{\link[flapper]{pf_archive-class}} object (from \code{\link[flapper]{pf}} plus \code{\link[flapper]{pf_simplify}} with \code{record = "archive"}).
#' @param bathy A \code{\link[raster]{raster}} that defines the grid across which \code{\link[flapper]{pf}} was applied.
#' @param sample_size (optional) An integer that defines the number of particles to sample from \code{archive} for the estimation. If \code{sample_size = NULL}, all particles are used. If specified, \code{sample_size} particles are sampled without replacement in line with their probability. Sampling is implemented to improve estimation speed.
#' @param estimate_ud A function (either \code{\link[adehabitatHR]{kernelUD}} or \code{\link[flapper]{kud_around_coastline}}) that estimates kernel utilisation distributions. The latter option is used when kernels need to be processed to account for barriers to movement that cannot be modelled via \code{\link[adehabitatHR]{kernelUD}}.
#' @param grid,... Arguments passed to \code{estimate_ud} (and ultimately \code{\link[adehabitatHR]{kernelUD}}, where they are defined) to estimate the kernel utilisation distribution. If  \code{\link[flapper]{kud_around_coastline}} is supplied to \code{estimate_ud}, then \code{grid} must be a \code{\link[sp]{SpatialPixelsDataFrame}}.
#' @param mask (optional) A spatial mask (see \code{\link[raster]{mask}}).
#' @param plot A logical input that defines whether or not to plot the KUD.
#'
#' @return The function returns a \code{\link[raster]{raster}} of the KUD.
#'
#' @examples
#' #### Define a grid across which to implement estimation
#' # This grid takes values of 0 on land and values of 1 in the sea
#' bathy <- dat_dcpf_histories$args$bathy
#' grid <- raster::raster(raster::extent(bathy), nrows = 100, ncols = 100)
#' raster::values(grid) <- 0
#' grid <- raster::mask(grid, dat_coast, updatevalue = 1)
#' grid <- methods::as(grid, "SpatialPixelsDataFrame")
#'
#' #### Example (1): Implement function using default options
#' pf_kud(pf_simplify(dat_dcpf_histories, return = "archive"),
#'        bathy = bathy,
#'        estimate_ud = kud_around_coastline, grid = grid)
#'
#' #### Example (2): Implement function using random sampling (for speed)
#' pf_kud(pf_simplify(dat_dcpf_histories, return = "archive"),
#'        bathy = bathy,
#'        sample_size = 100,
#'        estimate_ud = kud_around_coastline, grid = grid)
#'
#' @seealso \code{\link[flapper]{pf}}, \code{\link[flapper]{pf_simplify}}, \code{\link[flapper]{pf_plot_map}}, \code{\link[adehabitatHR]{kernelUD}}, \code{\link[flapper]{kud_around_coastline}}, \code{\link[flapper]{eval_by_kud}}
#' @export
#' @author Edward Lavender
#'

pf_kud <- function(archive,
                   bathy,
                   sample_size = 10,
                   estimate_ud = adehabitatHR::kernelUD,
                   grid, ...,
                   mask = NULL,
                   plot = TRUE){

  #### Checks
  check_class(input = archive, to_class = "pf_archive")
  if(archive$method != "pf_simplify"){
    warning("archive$method != 'pf_simplify'", immediate. = TRUE, call. = FALSE)
  }

  #### Extract particle coordinates
  pairs_df <-
    lapply(archive$history, function(elm) {
      elm[, c("id_current", "pr_current"), drop = FALSE]
    })
  pairs_df <- do.call(rbind, pairs_df)
  pairs_xy <- raster::xyFromCell(bathy, pairs_df$id_current)

  #### Sample locations according to their probability
  if(!is.null(sample_size)){
    if(sample_size > nrow(pairs_xy)){
      warning(paste0("'sample_size' (n = ", sample_size,
                     ") exceeds the number of sampled locations (n = ", nrow(pairs_xy),
                     "): implementing estimation using all sampled locations."),
              immediate. = TRUE)
      sample_size <- NULL
    }
    if(!is.null(sample_size)){
      pairs_xy <-
        pairs_xy[sample(x = 1:nrow(pairs_xy),
                        size = sample_size,
                        prob = pairs_df$pr_current), ]
    }
  }

  #### Implement KUD estimation
  pairs_xy_spdf <- sp::SpatialPointsDataFrame(
    pairs_xy,
    data = data.frame(ID = factor(rep(1, nrow(pairs_xy)))),
    proj4string = raster::crs(bathy))
  pairs_ud <- estimate_ud(xy = pairs_xy_spdf, grid = grid,...)

  #### Process KUD
  pairs_ud <- raster::raster(pairs_ud[[1]])
  if(!is.null(mask)) {
    if(raster::crs(bathy) != raster::crs(mask)){
      warning(paste0("'bathy' CRS ('", raster::crs(bathy),
                     "') is different from 'mask' CRS ('", raster::crs(mask), "')."),
              immediate. = TRUE)
    }
    pairs_ud <- raster::mask(pairs_ud, mask)
  }

  #### Visualise KUD
  if(plot) prettyGraphics::pretty_map(add_rasters = list(x = pairs_ud))

  #### Return KUD
  return(pairs_ud)
}
