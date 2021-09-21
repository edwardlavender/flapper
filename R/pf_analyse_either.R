######################################
######################################
#### pf_plot_map()

#' @title Plot `probability of use' from a PF algorithm
#' @description This function creates a \code{\link[raster]{raster}} of the `probability of use' across an area based on (a) particles sampled or (b) paths reconstructed by a particle filtering (PF) algorithm. To implement the function, an (a) \code{\link[flapper]{pf_archive-class}} object that contains connected particles (locations) sampled by \code{\link[flapper]{pf}} and processed by \code{\link[flapper]{pf_simplify}} or (b) \code{\link[flapper]{pf_path-class}} object that contains reconstructed paths must be supplied. The function extracts sampled locations and, for each location, calculates `the probability of use' for that location over the time series (see Details). This is (optionally) plotted and returned (invisibly) as a \code{\link[raster]{raster}}.
#' @param xpf A \code{\link[flapper]{pf_archive-class}} object (from \code{\link[flapper]{pf}} plus \code{\link[flapper]{pf_simplify}} with \code{return = "archive"}) or a \code{\link[flapper]{pf_path-class}} object (from \code{\link[flapper]{pf}} plus \code{\link[flapper]{pf_simplify}} with \code{return = "path"}).
#' @param map A \code{\link[raster]{raster}} that defines a grid across the area of interest.
#' @param transform (optional) A function to transform cell weights (e.g, \code{\link[base]{log}}).
#' @param scale A character that defines how \code{\link[raster]{raster}} values are scaled: \code{"original"} uses the original values; \code{"max"} scales values by the maximum value (so that, if \code{transform = NULL}, they lie between zero and one; and \code{"sum"} scales values by their sum so that they sum to one.
#' @param plot A logical input that defines whether or not to plot the map.
#' @param add_rasters If \code{plot = TRUE}, \code{add_rasters} is a named list, passed to \code{\link[prettyGraphics]{pretty_map}}, to customise the appearance of the plotted surface.
#' @param ... Additional arguments passed to \code{\link[prettyGraphics]{pretty_map}}.
#'
#' @details For particle-based implementations, this function is designed to be implemented for the subset of sampled particles that formed continuous paths from the start to the end of the time series (see \code{\link[flapper]{pf_simplify}}). At each time step, only one record of each location (derived by summarising the probabilities of multiple samples of the same location with the \code{summarise_pf} argument in \code{\link[flapper]{pf_simplify}}) should be passed to \code{\link[flapper]{pf_plot_map}} to ensure that cell scores depend on the number of time steps when the individual could have occupied a given cell, rather than the total number of samples of a location. For each location, the 'probability of use' is calculated as the sum of the number of times (time steps) that the location was sampled, weighted by the associated probabilities of each sample, over the total number of time steps. The benefit of this approach is that all particles that were part of paths from the start to the end of the time series are incorporated in the resultant map. However, this comes at the cost of simplifying cell probabilities for duplicate records and ignoring variation in the overall likelihood of different movement paths.
#'
#' For path-based implementations, the function is designed to be implemented for paths reconstructed by \code{\link[flapper]{pf_simplify}}. For each path, as for particle-based implementations, for each location, the 'probability of use' is calculated as the sum of the number of times (time steps) that the location was sampled, weighted by the associated probabilities of each sample, over the total number of time steps. (This is equivalent to calculating a weighted sum of the paths). Scores are then averaged across paths. This benefit of this approach is that it is possible to account for both location probabilities (in the weighted summation) and path probabilities (in the averaging of cell scores across paths, although this is not yet implemented). However, this approach can usually only be implemented for a subset of all possible paths (see \code{max_n_copies} in \code{\link[flapper]{pf_simplify}}) and these paths may not be independent (they may share substantial sections).
#'
#' For either implementation, raw scores can be transformed or scaled to facilitate comparisons.
#'
#' @return The function invisibly returns a \code{\link[raster]{raster}}, in which each cell contains the `probability of use' score  and (optionally) produces a plot of this surface.
#'
#' @examples
#' #### Prepare data
#' ## Particle-based implementation
#' # The example data 'dat_dcpf_histories' contains all particles sampled
#' # ... by an implementation of the DCPF algorithm. However, not all particles
#' # ... that were sampled at one time step may have been near to particles sampled
#' # ... at the next time step. In addition, some particles may have been sampled
#' # ... multiple times at one time step, but our maps of space use should reflect
#' # ... the number of time steps that the individual could have occupied a location,
#' # ... rather than the total number of samples of a location. Hence, to map
#' # ... space use, we should focus on the subset of particles that were connected
#' # ... between time steps and only retain one record of each particle at each time step
#' # ... using pf_simplify() with return = "archive"
#' dat_dcpf_histories_connected <-
#'   pf_simplify(dat_dcpf_histories,
#'              summarise_pr = TRUE,
#'              return = "archive")
#' ## Path based implementation
#' # The example data 'dat_dcpf_paths' contains a sample of paths reconstructed
#' # ... by the DCPF algorithm and we can also implement the function for these paths.
#'
#' #### Example (1): Implement the function with default options
#' pp <- par(mfrow = c(1, 2))
#' pf_plot_map(dat_dcpf_histories_connected, map = dat_dc$args$bathy)
#' pf_plot_map(dat_dcpf_paths, map = dat_dc$args$bathy)
#' par(pp)
#'
#' #### Example (2): Re-scale the map(s)
#' pp <- par(mfrow = c(2, 2))
#' pf_plot_map(dat_dcpf_histories_connected, map = dat_dc$args$bathy, scale = "max")
#' pf_plot_map(dat_dcpf_histories_connected, map = dat_dc$args$bathy, scale = "sum")
#' pf_plot_map(dat_dcpf_paths, map = dat_dc$args$bathy, scale = "max")
#' pf_plot_map(dat_dcpf_paths, map = dat_dc$args$bathy, scale = "sum")
#' par(pp)
#'
#' #### Example (3): Customise the map(s)
#' pp <- par(mfrow = c(1, 2))
#' pf_plot_map(dat_dcpf_histories_connected, map = dat_dc$args$bathy,
#'             add_rasters = list(col = grDevices::grey.colors(n = 100)),
#'             xlab = "x", ylab = "y")
#' pf_plot_map(dat_dcpf_paths, map = dat_dc$args$bathy,
#'             add_rasters = list(col = grDevices::topo.colors(n = 100)),
#'             xlab = "x", ylab = "y")
#' par(pp)
#'
#' @seealso \code{\link[flapper]{pf}} implements PF. \code{\link[flapper]{pf_simplify}} assembles paths from particle histories. \code{\link[flapper]{pf_plot_history}} visualises particle histories. \code{\link[flapper]{pf_plot_1d}}, \code{\link[flapper]{pf_plot_2d}} and \code{\link[flapper]{pf_plot_3d}} provide plotting routines for paths. \code{\link[flapper]{pf_loglik}} calculates the log-likelihood of each path.
#' @author Edward Lavender
#' @export

pf_plot_map <- function(xpf,
                        map,
                        transform = NULL,
                        scale = c("original", "max", "sum"),
                        plot = TRUE,
                        add_rasters = list(),...){
  #### Check inputs
  check_class(input = xpf, to_class = c("pf_archive", "pf_path"))
  scale <- match.arg(scale)

  #### Extract cell locations and associated probabilities

  ## pf_archive implementation
  # Extract particle histories as a single dataframe (and check for duplicate particles)
  if(inherits(xpf, "pf_archive")){
    dat <- lapply(1:length(xpf$history), function(t) {
      elm <- xpf$history[[t]][, c("id_current", "pr_current"), drop = FALSE]
      if(any(duplicated(elm$id_current))) {
        if(xpf$method != "pf_simplify"){
          warning(paste0("xpf$history[[", t, "]] contains duplicate cells. ",
                         "Implementing pf_simplify() with 'summarise_pr = TRUE' and return = 'archive' specified first is advised."),
                  immediate. = TRUE, call. = FALSE)
        } else {
          warning(paste0("xpf$history[[", t, "]] contains duplicate cells. ",
                         "Did you implement pf_simplify() without 'summarise_pr = TRUE'?"),
                  immediate. = TRUE, call. = FALSE)
        }
      }
      return(elm)
    })
    n <- length(dat)
    dat <- do.call(rbind, dat)
    colnames(dat) <- c("cell_id", "cell_pr")
    dat$path_id <- 1L

    ## pf_path implementation
  } else if(inherits(xpf, "pf_path")){
    dat <- xpf[, c("cell_id", "cell_pr", "path_id")]
    n <- length(which(dat$path_id == dat$path_id[1]))
  }

  #### Calculate cell scores (the frequency with which each cell was sampled, weighted by the probability)
  # For particle-based implementations,
  # ... for each location, we calculate the 'score' as the frequency
  # ... with which the cell was sampled weighted by the probability.
  # For path-based implementations, for each path, for each location
  # ... we calculate the number of times the cell was sampled weighed by the probability.
  # ... We then combine scores across paths by taking the average frequency
  # ... (ideally weighted by the overall likelihood of the path).
  # ... For path based implementations, this approach is equivalent (but faster than)
  # ... to defining a stack of rasters, one for each path, and then calculating the average/a weighed average.
  wt_freq <-
    dat %>%
    dplyr::group_by(.data$path_id, .data$cell_id) %>%
    dplyr::summarise(score = sum(.data$cell_pr)/n) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(.data$cell_id) %>%
    dplyr::summarise(score = mean(.data$score))

  #### Transform and scale cell scores
  # Transform weighted frequencies
  if(!is.null(transform)) {
    wt_freq$score <- transform(wt_freq$score)
    if(any(is.na(wt_freq$score))) stop("'transform' function has created NAs.")
  }
  # Re-scale weighted frequencies e.g. so that the maximum value has a score of one
  if(scale == "max"){
    wt_freq$score <- wt_freq$score/max(wt_freq$score)
  } else if(scale == "sum"){
    wt_freq$score <- wt_freq$score/sum(wt_freq$score)
  }

  #### Assign scores to map
  p <- raster::setValues(map, 0)
  p[wt_freq$cell_id] <- wt_freq$score
  p <- raster::mask(p, map)

  #### Plot map and return outputs
  if(plot){
    if(!is.null(add_rasters)) add_rasters$x <- p
    prettyGraphics::pretty_map(x = p,
                               add_rasters = add_rasters,...)
  }
  return(invisible(p))
}


######################################
######################################
#### pf_kud_*()

#' @title Apply kernel smoothing to particles or paths from a PF algorithm
#' @description These functions are wrappers designed to apply kernel utilisation distribution (KUD) estimation to the outputs of a particle filtering (PF) algorithm. To implement these routines, an (a) \code{\link[flapper]{pf_archive-class}} object from \code{\link[flapper]{pf}} (plus \code{\link[flapper]{pf_simplify}} with the \code{return = "archive"} argument) containing particle histories for connected particles or (b) a \code{\link[flapper]{pf_path-class}} object containing reconstructed paths must be supplied. Depending on the implementation, using a subset, all or an expanded sample of sampled locations, the functions apply KUD smoother(s) via a user-supplied estimation routine (i.e., \code{\link[adehabitatHR]{kernelUD}} or \code{\link[flapper]{kud_around_coastline}}). The functions extract the utilisation distribution(s) as \code{\link[raster]{raster}}(s), combine distribution(s) (if necessary), apply a spatial mask (e.g. coastline), plot the processed distribution (if specified) and return this as a \code{\link[raster]{raster}}.
#'
#' @param xpf A \code{\link[flapper]{pf_archive-class}} object (from \code{\link[flapper]{pf}} plus \code{\link[flapper]{pf_simplify}} with \code{return = "archive"}) or a \code{\link[flapper]{pf_path-class}} object (from \code{\link[flapper]{pf}} plus \code{\link[flapper]{pf_simplify}} with \code{return = "path"}). For particle-based implementations, for \code{\link[flapper]{pf_kud_1}}, \code{\link[flapper]{pf_simplify}} should be implemented with \code{summarise_pr = FALSE}; for \code{\link[flapper]{pf_kud_2}}, \code{\link[flapper]{pf_simplify}} should be implemented with \code{summarise_pr = TRUE} (see Details).
#' @param bathy A \code{\link[raster]{raster}} that defines the grid across which \code{\link[flapper]{pf}} was applied. This used to extract cell coordinates and to express KUD(s).
#' @param sample_size (optional) An integer that defines the number of particles to sample from (a) particle histories or (b) each path in \code{xpf} for KUD estimation. For \code{\link[flapper]{pf_kud_1}}, sampling is used to account for particle uncertainty so \code{sample_size} should be greater than the number of particle samples at each time step (see Details). If specified, for each time step, \code{sample_size} particles are sampled from (a) particle histories or (b) reconstructed paths with replacement in line with their probability. For \code{\link[flapper]{pf_kud_2}}, sampling is used to reduce memory requirements and computation time, so \code{sample_size} should be lower than the total number of particle samples (per path, if applicable). If specified, \code{sample_size} particles are sampled from (a) particle histories or (b) each path with replacement in line with their probability. If \code{sample_size = NULL}, all particles are used.
#' @param estimate_ud A function (either \code{\link[adehabitatHR]{kernelUD}} or \code{\link[flapper]{kud_around_coastline}}) that estimates kernel utilisation distributions. For \code{\link[flapper]{pf_kud_1}}, \code{\link[flapper]{kud_around_coastline_fast}} can also be used for faster estimation.
#' @param grid,... Arguments passed to \code{estimate_ud} (and ultimately \code{\link[adehabitatHR]{kernelUD}}, where they are defined) to estimate the kernel utilisation distribution. If \code{\link[flapper]{kud_around_coastline}} or \code{\link[flapper]{kud_around_coastline_fast}} is supplied to \code{estimate_ud}, then \code{grid} must be a \code{\link[sp]{SpatialPixelsDataFrame}}. However, note that in all cases, KUD(s) are resampled onto \code{bathy}.
#' @param scale For \code{\link[flapper]{pf_kud_1}}, \code{scale} is a logical input that defines whether or not to scale the KUD for each time step such that the most probable locations are assigned a score of one.
#' @param plot_by_time,prompt For \code{\link[flapper]{pf_kud_1}}, \code{plot_by_time} is a logical variable that defines whether or not to plot the cumulative (un-normalised) KUD for each time step. If supplied, \code{prompt} is a logical variable that defines whether or not to pause function execution between sequential plots. These arguments are not implemented in parallel (see \code{cl}, below).
#' @param chunks,cl,varlist For \code{\link[flapper]{pf_kud_1}}, \code{chunks}, \code{cl} and \code{varlist} are chunk-wise implementation controls. \code{chunks} is an integer that defines the number of chunks into which to split particle/path time series. To minimise memory requirements, within each chunk, a blank map is sequentially updated with the KUD for each time step; the cumulative KUD for each chunk is then summed across chunks to create a single KUD. This approach minimises memory use while facilitating improvements in computation time through the parallel processing of each chunk via \code{cl} and \code{varlist}. \code{cl} is cluster object created by \code{\link[parallel]{makeCluster}} to implement the algorithm in parallel. If \code{cl} is supplied, \code{varlist} may also be required. This is a character vector of objects to export. \code{varlist} is passed to the \code{varlist} of \code{\link[parallel]{clusterExport}}. Exported objects must be located in the global environment.
#' @param mask (optional) A spatial mask (see \code{\link[raster]{mask}}).
#' @param plot A logical input that defines whether or not to plot the KUD.
#' @param verbose A logical input that defines whether or not to print messages to the console to monitor function progress.
#'
#' @details
#'
#' \subsection{Methods}{These function create smooth KUD representations of particle or path samples from a PF algorithm (see \code{\link[flapper]{pf_plot_map}}). Two different methods are implemented.
#' }
#'
#' \subsection{Method (1)}{\code{\link[flapper]{pf_kud_1}} implements KUD estimation by fitting a single kernel to sampled locations at each time step; time step-specific kernels are then summed and re-normalised. This method can be implemented for (a) the sampled particles that formed continuous paths from the start to the end of the time series (see \code{\link[flapper]{pf_simplify}}) or (b) reconstructed paths. The method has two main advantages. First, the kernel bandwidth is allowed to vary through time. Second, this implementation permits a Bayesian-style resampling process that can be used to account for particle uncertainty. Specifically, at each time step, large numbers of particles can be re-sampled, in line with their probability and with replacement, from the initial list of sampled particles; more probable locations are sampled more often and consequently have more influence on the KUD for each time step, thus accounting for particle probability.
#'
#' A limitation with this method is that the fitting a KUD to each time step can be memory intensive and computationally intensive. To minimise memory requirements, by default (\code{chunks = 1L}), the function starts with a blank map and iterates over each time step, sequentially adding KUDs to the map at each step. By continuously updating a single map, this option minimises memory requirements but is slow. A faster option is to split the time series into chunks, implement an iterative option within each chunk, and then join maps for each chunk. The advantage of this option is that memory use remains limited while computation time can be improved by the parallel processing of each chunk. This is implemented via \code{chunks}, \code{cl} and \code{varlist}.
#' }
#'
#' \subsection{Method (2)}{\code{\link[flapper]{pf_kud_2}} implements KUD estimation by fitting a single kernel to all sampled locations or to each path (and then aggregating KUDs across paths). The main advantage of this method is speed: unlike \code{\link[flapper]{pf_kud_1}}, KUDs are not fitted to the locations for each time step. The limitations are that kernel bandwidth is constant for all time steps and in most situations re-sampling cannot be used in the same way to account for particle uncertainty (due to memory limitations).
#'
#'  For particle-based implementations, this method is designed to be implemented for the subset of unique, sampled particles that formed continuous paths from the start to the end of the time series (see \code{\link[flapper]{pf_simplify}} and \code{\link[flapper]{pf_plot_map}}). These particles are used for KUD estimation. By default, all particles are used, but \eqn{n =} \code{sample_size} particles can be sampled at random, in line with their probability, if specified, for faster KUD estimation. Selected particles are then used to estimate a KUD by effectively treating samples as `relocations', ultimately via \code{\link[adehabitatHR]{kernelUD}}. This distribution is then processed, plotted and returned.
#'
#' For path-based implementations, this function is designed to be implemented for paths reconstructed by \code{\link[flapper]{pf_simplify}}. As for the particle-based implementation, for each path, all locations, or random sample of \code{sample_size} locations are used to estimate a KUD by treating sampled locations as `relocations'. KUDs are processed and combined across paths into a single, average KUD. The advantage of this approach is that the overall probability of the paths can be incorporated in the estimation procedure via sampling or weights when path-specific KUDs are averaged (although that is not yet implemented).
#' }
#'
#' @return The functions return a \code{\link[raster]{raster}} of the KUD.
#'
#' @examples
#' #### Define particle samples for smoothing
#' # To do this, we will re-implement the pf() for the example dat_dcpf_histories
#' # ... dataset, but with a larger number of particles. This is necessary
#' # ... because kernel smoothing is only appropriate if there are
#' # ... enough locations to permit smoothing.
#' set.seed(1)
#' dcpf_args <- dat_dcpf_histories$args
#' dcpf_args$calc_distance_euclid_fast <- TRUE
#' dcpf_args$n <- 250L
#' out_dcpf_particles <- do.call(pf, dcpf_args)
#'
#' #### Process particles and paths
#' out_dcpf_particles_1 <-
#'   pf_simplify(out_dcpf_particles, return = "archive")
#' out_dcpf_particles_2 <-
#'   pf_simplify(out_dcpf_particles, summarise_pr = TRUE, return = "archive")
#' out_dcpf_paths <-
#'   pf_simplify(out_dcpf_particles, max_n_paths = 100L)
#'
#' #### Define a grid across which to implement estimation
#' # This grid takes values of 0 on land and values of 1 in the sea
#' bathy <- out_dcpf_particles$args$bathy
#' grid <- raster::raster(raster::extent(bathy), nrows = 100, ncols = 100)
#' raster::values(grid) <- 0
#' grid <- raster::mask(grid, dat_coast, updatevalue = 1)
#' grid <- methods::as(grid, "SpatialPixelsDataFrame")
#'
#' #### Example (1): Implement pf_kud_1() using default options
#' ## Implementation based on particles
#' pp <- par(mfrow = c(1, 2), mar = c(3, 3, 3, 3))
#' pf_kud_1(out_dcpf_particles_1,
#'          bathy = bathy,
#'          sample_size = 500,
#'          estimate_ud = kud_around_coastline_fast, grid = grid)
#' ## Implementation based on paths
#' \dontrun{
#' pf_kud_1(out_dcpf_paths,
#'          bathy = bathy,
#'          sample_size = 500,
#'          estimate_ud = kud_around_coastline_fast, grid = grid)
#' prettyGraphics::add_sp_path(x = out_dcpf_paths$cell_x, y = out_dcpf_paths$cell_y,
#'                             length = 0.01)
#' }
#' par(pp)
#'
#' #### Example (2): Implement pf_kud_2() using default options
#' ## Implementation based on particles
#' pp <- par(mfrow = c(1, 2), mar = c(3, 3, 3, 3))
#' pf_kud_2(out_dcpf_particles_2,
#'          bathy = bathy,
#'          estimate_ud = kud_around_coastline, grid = grid)
#' ## Implementation based on paths
#' pf_kud_2(out_dcpf_paths,
#'          bathy = bathy,
#'          estimate_ud = kud_around_coastline, grid = grid)
#' prettyGraphics::add_sp_path(x = out_dcpf_paths$cell_x, y = out_dcpf_paths$cell_y,
#'                             length = 0.01)
#' par(pp)
#'
#' #### Example (3): For improved speed with pf_kud_1(), use parallelisation
#' ## Implementation based on particles
#' pp <- par(mfrow = c(1, 2), mar = c(3, 3, 3, 3))
#' pf_kud_1(out_dcpf_particles_1,
#'          bathy = bathy,
#'          sample_size = 500,
#'          estimate_ud = flapper::kud_around_coastline_fast, grid = grid,
#'          chunks = 2L,
#'          cl = parallel::makeCluster(2L))
#' ## Implementation based on paths
#' \dontrun{
#' cl <- parallel::makeCluster(2L)
#' parallel::clusterEvalQ(cl = cl, library(raster))
#' pf_kud_1(out_dcpf_paths,
#'          bathy = bathy,
#'          sample_size = 500,
#'          estimate_ud = flapper::kud_around_coastline_fast, grid = grid,
#'          chunks = 2L,
#'          cl = cl)
#' }
#' par(pp)
#'
#' #### Example (4): For improved speed with pf_kud_2(), use sample_size
#' ## Implementation based on particles
#' pp <- par(mfrow = c(1, 2), mar = c(3, 3, 3, 3))
#' pf_kud_2(out_dcpf_particles_2,
#'          bathy = bathy,
#'          sample_size = 50, # sample 50 particles overall
#'          estimate_ud = kud_around_coastline, grid = grid)
#' ## Implementation based on paths
#' pf_kud_2(out_dcpf_paths,
#'          bathy = bathy,
#'          sample_size = 50, # sample 50 particles per path
#'          estimate_ud = kud_around_coastline, grid = grid)
#' par(pp)
#'
#' #### Example (5): Compare pf_kud_1() and pf_kud_2()
#' pp <- par(mfrow = c(2, 2), mar = c(3, 3, 3, 3))
#' kud_1a <- pf_kud_1(xpf = out_dcpf_particles_1,
#'                    bathy = out_dcpf_particles$args$bathy,
#'                    sample_size = out_dcpf_particles$args$n,
#'                    estimate_ud = kud_around_coastline_fast, grid = grid,
#'                    plot = TRUE)
#'
#' kud_1b <- pf_kud_1(xpf = out_dcpf_particles_1,
#'                    bathy = out_dcpf_particles$args$bathy,
#'                    sample_size = 500,
#'                    estimate_ud = kud_around_coastline_fast, grid = grid,
#'                    plot = TRUE)
#'
#' kud_1c <- pf_kud_1(xpf = out_dcpf_particles_1,
#'                    bathy = out_dcpf_particles$args$bathy,
#'                    sample_size = 5000,
#'                    estimate_ud = kud_around_coastline_fast, grid = grid,
#'                    plot = TRUE)
#'
#' kud_1d <- pf_kud_2(xpf = out_dcpf_particles_2,
#'                  bathy = out_dcpf_particles$args$bathy,
#'                  estimate_ud = kud_around_coastline, grid = grid,
#'                  plot = TRUE)
#' par(pp)
#'
#' @seealso \code{\link[flapper]{pf}}, \code{\link[flapper]{pf_simplify}}, \code{\link[flapper]{pf_plot_map}}, \code{\link[adehabitatHR]{kernelUD}}, \code{\link[flapper]{kud_around_coastline}}, \code{\link[flapper]{kud_around_coastline_fast}}, \code{\link[flapper]{eval_by_kud}}
#' @author Edward Lavender
#' @name pf_kud
NULL


######################################
#### pf_kud_1()

#' @rdname pf_kud
#' @export

pf_kud_1 <- function(xpf,
                     bathy,
                     sample_size = NULL,
                     estimate_ud = adehabitatHR::kernelUD,
                     grid,..., scale = FALSE,
                     plot_by_time = FALSE, prompt = TRUE,
                     chunks = 1L, cl = NULL, varlist = NULL,
                     mask = NULL, plot = TRUE,
                     verbose = TRUE
){

  #### Checks
  cat_to_console <- function(..., show = verbose){
    if(show) cat(paste(..., "\n"))
  }
  t_onset <- Sys.time()
  cat_to_console(paste0("flapper::pf_kud_1() called (@ ", t_onset, ")..."))
  cat_to_console("... Setting up function...")
  check_class(input = xpf, to_class = c("pf_archive", "pf_path"))
  crs_bathy <- crs_grid <- crs_mask <- NULL
  crs_bathy <- raster::crs(bathy)
  if(!is.null(grid)) if(!inherits(grid, c("numeric", "integer"))) crs_grid <- raster::crs(grid)
  if(!is.null(mask)) crs_mask <- raster::crs(mask)
  crs <- list(crs_bathy, crs_grid, crs_mask)
  crs <- compact(crs)
  if(!length(unique(crs)) == 1L){
    warning("The CRS(s) of spatial layer(s) ('bathy', 'grid' and 'mask', if applicable) are not identical.",
            immediate. = TRUE, call. = FALSE)
  }
  crs_spp <- sapply(crs, function(x) !is.na(x))
  if(any(crs_spp)){
    crs <- crs[[which(crs_spp)[1]]]
  } else {
    crs <- sp::CRS(as.character(NA))
  }
  message("CRS taken as: '", crs, "'.")
  if(plot_by_time & !is.null(cl)) {
    warning("'plot_by_time' ignored for 'cl' != NULL.", immediate. = TRUE, call. = FALSE)
    plot_by_time <- FALSE
  }
  if(chunks <= 1L & !is.null(cl)){
    warning("'cl' supplied but 'chunks <= 1L': ignoring 'cl'.", immediate. = TRUE, call. = FALSE)
    if(inherits(cl, "cluster")) parallel::stopCluster(cl)
    cl <- NULL
  }

  #### Define a dataframe of particle locations through time
  cat_to_console("... Processing sampled locations...")
  if(inherits(xpf, "pf_archive")){
    particles_by_t <-
      lapply(1:length(xpf$history), function(t) {
        elm <- xpf$history[[t]]
        history_for_t <- elm[, c("id_current", "pr_current"), drop = FALSE]
        colnames(history_for_t) <- c("cell_id", "cell_pr")
        history_for_t$timestep <- t
        return(history_for_t)
      })
    particles <- do.call(rbind, particles_by_t)
    particles[, c("cell_x", "cell_y")] <- raster::xyFromCell(bathy, particles$cell_id)
  } else if(inherits(xpf, "pf_path")){
    particles <- data.frame(xpf)
    particles <- particles[, c("timestep", "cell_id", "cell_pr", "cell_x", "cell_y")]
  }
  particles <- particles %>% dplyr::arrange(.data$timestep)

  #### Define chunks
  cat_to_console("... Defining chunk(s)...")
  if(chunks <= 1L){
    particles$chunk <- 1L
  } else{
    particles$chunk <- cut(particles$timestep, chunks)
  }
  particles_by_chunk <- split(particles, particles$chunk)

  #### Implement KUD estimation chunk-wise
  cat_to_console("... Implementing KUD estimation over chunk(s)...")
  # Define blank raster that provides the foundation for KUDs
  blank <- raster::setValues(bathy, 0)
  # Loop over each 'chunk' and calculate cumulative KUD for that chunk
  if(!is.null(cl) & !is.null(varlist)) parallel::clusterExport(cl = cl, varlist = varlist)
  kud_by_chunk <- pbapply::pblapply(particles_by_chunk, cl = cl, function(particles_for_chunk){
    # Isolate particles for that chunk
    particles_by_t <- split(particles_for_chunk, particles_for_chunk$timestep)
    # Define blank UD
    ud <- blank
    # Loop over each time step and add to UD
    for(t in 1:length(particles_by_t)){
      # Get particles for time step
      particles_for_t <- particles_by_t[[t]]
      # If the number of unique locations is fewer than 5, we won't implement kernel smoothing
      if(length(unique(particles_for_t$cell_id)) < 5){
        message(paste0("Fewer than five unique cells for kernel estimation at time ", t, "."))
        denom <- sum(particles_for_t$cell_pr)
        particles_for_t <-
          particles_for_t %>%
          dplyr::group_by(.data$cell_id) %>%
          dplyr::mutate(cell_pr = sum(.data$cell_pr/denom)) %>%
          dplyr::slice(1L) %>%
          data.frame()
        ud_for_t <- blank
        ud_for_t[particles_for_t$cell_id] <- particles_for_t$cell_pr
      # Implement kernel smoothing
      } else {
        # Implement resampling, if specified, within the loop (to minimise memory requirements)
        if(!is.null(sample_size)){
          particles_for_t <-
            particles_for_t %>%
            dplyr::slice_sample(n = sample_size, weight_by = .data$cell_pr, replace = TRUE)
        }
        # Get locations
        xy_for_t <- sp::SpatialPointsDataFrame(particles_for_t[, c("cell_x", "cell_y")],
                                               data = data.frame(ID = factor(rep(1, nrow(particles_for_t)))),
                                               proj4string = crs)
        # Get UD for timestep
        ud_for_t <- estimate_ud(xy_for_t, grid = grid,...)
        if(inherits(ud_for_t, "estUD")){
          ud_for_t <- raster::raster(ud_for_t)
        } else if(inherits(ud_for_t, "estUDm")){
          ud_for_t <- raster::raster(ud_for_t[[1]])
        }
        # Resample UD onto common grid
        ud_for_t <- raster::resample(ud_for_t, blank, method = "bilinear")
      }
      # Add UD to cumulative UD for chunk
      if(scale){
        min_ud <- raster::cellStats(ud_for_t, "min")
        max_ud <- raster::cellStats(ud_for_t, "max")
        ud_for_t <- (ud_for_t - min_ud)/(max_ud - min_ud)
      }
      ud <- sum(ud, ud_for_t, na.rm = TRUE)
      # Visualise cumulative UD (if specified)
      if(plot_by_time){
        raster::plot(ud, main = paste("chunk", particles_for_t$chunk[1], ": time ", t))
        if(prompt) readline(prompt = "Press [Enter] to continue...")
      }
    }
    return(ud)
  })
  if(!is.null(cl)) parallel::stopCluster(cl = cl)

  #### Process KUDs
  cat_to_console("... Processing KUD(s)...")
  ## Summarise across chunks
  if(chunks < 2L){
    kud <- kud_by_chunk[[1]]
  } else {
    kud_stack <- raster::stack(kud_by_chunk)
    kud <- raster::calc(kud_stack, sum, na.rm = TRUE)
  }
  ## Mask KUD
  if(!is.null(mask)) kud <- raster::mask(kud, mask)
  kud <- kud/raster::cellStats(kud, "sum")

  #### Plot KUD
  cat_to_console("... Plotting KUD...")
  ext <- raster::extent(bathy)
  if(plot) prettyGraphics::pretty_map(add_rasters = list(x = kud),
                                      xlim = ext[1:2], ylim = ext[3:4])

  #### Return kud
  t_end <- Sys.time()
  total_duration <- difftime(t_end, t_onset, units = "mins")
  cat_to_console(paste0("... flapper::pf_kud_1() call completed (@ ", t_end, ") after ~", round(total_duration, digits = 2), " minutes."))
  return(kud)

}


######################################
#### pf_kud_2()

#' @rdname pf_kud
#' @export

pf_kud_2 <- function(xpf,
                     bathy,
                     sample_size = NULL,
                     estimate_ud = adehabitatHR::kernelUD,
                     grid, ...,
                     mask = NULL,
                     plot = TRUE,
                     verbose = TRUE){

  #### Checks
  cat_to_console <- function(..., show = verbose){
    if(show) cat(paste(..., "\n"))
  }
  t_onset <- Sys.time()
  cat_to_console(paste0("flapper::pf_kud_2() called (@ ", t_onset, ")..."))
  cat_to_console("... Setting up function...")
  check_class(input = xpf, to_class = c("pf_archive", "pf_path"))
  if(inherits(xpf, "pf_path")){
    if(isTRUE(all.equal(estimate_ud, kud_around_coastline_fast))){
      stop("flapper::kud_around_coastline_fast() cannot be used for 'pf_path' class objects.", call. = FALSE)
    }
  }
  crs_bathy <- crs_grid <- crs_mask <- NULL
  crs_bathy <- raster::crs(bathy)
  if(!is.null(grid)) if(!inherits(grid, c("numeric", "integer"))) crs_grid <- raster::crs(grid)
  if(!is.null(mask)) crs_mask <- raster::crs(mask)
  crs <- list(crs_bathy, crs_grid, crs_mask)
  crs <- compact(crs)
  if(!length(unique(crs)) == 1L){
    warning("The CRS(s) of spatial layer(s) ('bathy', 'grid' and 'mask', if applicable) are not identical.",
            immediate. = TRUE, call. = FALSE)
  }
  crs_spp <- sapply(crs, function(x) !is.na(x))
  if(any(crs_spp)){
    crs <- crs[[which(crs_spp)[1]]]
  } else {
    crs <- sp::CRS(as.character(NA))
  }
  message("CRS taken as: '", crs, "'.")
  blank <- raster::setValues(bathy, 0)


  ######################################
  #### pf_archive approach

  if(inherits(xpf, "pf_archive")){

    #### Get cell probabilities and coordinates
    cat_to_console("... Processing sampled locations...")
    if(xpf$method != "pf_simplify"){
      warning("xpf$method != 'pf_simplify'", immediate. = TRUE, call. = FALSE)
    }
    particles_by_t <-
      lapply(1:length(xpf$history), function(t) {
        elm <- xpf$history[[t]]
        history_for_t <- elm[, c("id_current", "pr_current"), drop = FALSE]
        colnames(history_for_t) <- c("cell_id", "cell_pr")
        history_for_t$timestep <- t
        return(history_for_t)
      })
    particles <- do.call(rbind, particles_by_t)
    particles[, c("cell_x", "cell_y")] <- raster::xyFromCell(bathy, particles$cell_id)

    #### Sample locations according to their probability
    # For consistency with the paths approach below, particles are not sampled within time steps
    if(!is.null(sample_size)){
      n_particles <- length(which(particles$timestep == particles$timestep[1]))
      if(sample_size > n_particles) warning("'sample_size' is greater than the number of particles.",
                                            immediate. = TRUE, call. = FALSE)
      particles <-
        particles %>%
        # dplyr::group_by(.data$timestep) %>%
        dplyr::slice_sample(n = sample_size, weight_by = .data$cell_pr, replace = TRUE)
    }

    #### Implement KUD estimation
    cat_to_console("... Implementing KUD estimation...")
    particles_spdf <- sp::SpatialPointsDataFrame(
      particles[, c("cell_x", "cell_y")],
      data = data.frame(ID = factor(rep(1, nrow(particles)))),
      proj4string = crs)
    ud <- estimate_ud(xy = particles_spdf, grid = grid,...)

    #### Process KUD
    cat_to_console("... Processing KUD(s)...")
    if(inherits(ud, "estUDm")) ud <- ud[[1]]
    if(!inherits(ud, "RasterLayer")){
      ud <- raster::raster(ud)
    }
    ud <- raster::resample(ud, bathy)
    if(!is.null(mask)) ud <- raster::mask(ud, mask)


    ######################################
    #### pf_path approach

  } else if(inherits(xpf, "pf_path")){

    #### Get cell coordinates
    particles <- data.frame(xpf)

    #### Sub sample
    if(!is.null(sample_size)){
      n_particles <- length(which(particles$path_id == particles$path_id[1]))
      if(sample_size > n_particles) warning("'sample_size' is greater than the number of particles.",
                                            immediate. = TRUE, call. = FALSE)
      particles$index <- 1:nrow(particles)
      particles <-
        particles %>%
        dplyr::group_by(.data$path_id, .groups = "drop_last") %>%
        dplyr::slice_sample(n = sample_size, weight_by = .data$cell_pr, replace = TRUE) %>%
        data.frame()
    }

    #### Implement KUD estimation for each path
    # This approach cannot be implemented using kud_around_coastline_fast()
    cat_to_console("... Implementing KUD estimation...")
    particles_spdf <- sp::SpatialPointsDataFrame(
      particles[, c("cell_x", "cell_y")],
      data = data.frame(ID = factor(particles$path_id)),
      proj4string = crs)
    ud_by_path <- estimate_ud(xy = particles_spdf, grid = grid,...)

    #### Combine KUDs across paths
    ## Convert KUDs to rasterStack
    cat_to_console("... Processing KUD(s)...")
    ud_by_path <- pbapply::pblapply(ud_by_path, function(ud){
      ud <- raster::raster(ud)
      ud <- raster::resample(ud, blank)
      if(!is.null(mask)) ud <- raster::mask(ud, mask)
      return(ud)
    })
    ud_by_path <- raster::stack(ud_by_path)
    ud <- raster::calc(ud_by_path, mean)
  }

  #### Renormalise KUDs
  ud <- ud/raster::cellStats(ud, "sum")

  #### Visualise KUD
  cat_to_console("... Plotting KUD...")
  ext <- raster::extent(bathy)
  if(plot) prettyGraphics::pretty_map(add_rasters = list(x = ud),
                                      xlim = ext[1:2], ylim = ext[3:4])

  #### Return KUD
  t_end <- Sys.time()
  total_duration <- difftime(t_end, t_onset, units = "mins")
  cat_to_console(paste0("... flapper::pf_kud_2() call completed (@ ", t_end, ") after ~", round(total_duration, digits = 2), " minutes."))
  return(ud)
}
