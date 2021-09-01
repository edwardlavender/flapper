######################################
######################################
#### pf_plot_map()

#' @title Plot `probability of use' from a PF algorithm
#' @description This function creates a plot of the `probability of use' across an area based on (a) particles sampled or (b) paths reconstructed by a particle filtering (PF) algorithm. To implement the function, an (a) \code{\link[flapper]{pf_archive-class}} object that contains connected particles (locations) sampled by \code{\link[flapper]{pf}} and processed by \code{\link[flapper]{pf_simplify}} or (b) \code{\link[flapper]{pf_path-class}} object that contains reconstructed paths must be supplied. The function extracts sampled locations and, for each location, calculates `the probability of use' for that location over the time series (see Details). This is plotted and returned (invisibly) as a \code{\link[raster]{raster}}.
#' @param xpf A \code{\link[flapper]{pf_archive-class}} object (from \code{\link[flapper]{pf}} plus \code{\link[flapper]{pf_simplify}} with \code{return = "archive"}) or a \code{\link[flapper]{pf_path-class}} object (from \code{\link[flapper]{pf}} plus \code{\link[flapper]{pf_simplify}} with \code{return = "path"}).
#' @param map A \code{\link[raster]{raster}} that defines a grid across the area of interest.
#' @param transform (optional) A function to transform cell weights (e.g, \code{\link[base]{log}}).
#' @param scale A character that defines how \code{\link[raster]{raster}} values are scaled: \code{"original"} uses the original values; \code{"max"} scales values by the maximum value (so that, if \code{transform = NULL}, they lie between zero and one; and \code{"sum"} scales values by their sum so that they sum to one.
#' @param add_rasters A named list, passed to \code{\link[prettyGraphics]{pretty_map}}, to customise the appearance of the plotted surface.
#' @param ... Additional arguments passed to \code{\link[prettyGraphics]{pretty_map}}.
#'
#' @details For particle-based implementations, this function is designed to be implemented for the subset of sampled particles that formed continuous paths from the start to the end of the time series (see \code{\link[flapper]{pf_simplify}}). At each time step, only one record of each location (derived by summarising the probabilities of multiple samples of the same location with the \code{summarise_pf} argument in \code{\link[flapper]{pf_simplify}}) should be passed to \code{\link[flapper]{pf_plot_map}} to ensure that cell scores depend on the number of time steps when the individual could have occupied a given cell, rather than the total number of samples of a location. For each location, the 'probability of use' is calculated as the sum of the number of times (time steps) that the location was sampled, weighted by the associated probabilities of each sample, over the total number of time steps. The benefit of this approach is that all particles that were part of paths from the start to the end of the time series are incorporated in the resultant map. However, this comes at the cost of simplifying cell probabilities for duplicate records and ignoring variation in the overall likelihood of different movement paths.
#'
#' For path-based implementations, the function is designed to be implemented for paths reconstructed by \code{\link[flapper]{pf_simplify}}. For each path, as for particle-based implementations, for each location, the 'probability of use' is calculated as the sum of the number of times (time steps) that the location was sampled, weighted by the associated probabilities of each sample, over the total number of time steps. (This is equivalent to calculating a weighted sum of the paths). Scores are then averaged across paths. This benefit of this approach is that it is possible to account for both location probabilities (in the weighted summation) and path probabilities (in the averaging of cell scores across paths, although this is not yet implemented). However, this approach can usually only be implemented for a subset of all possible paths (see \code{max_n_copies} in \code{\link[flapper]{pf_simplify}}) and these paths may not be independent (they may share substantial sections).
#'
#' For either implementation, raw scores can be transformed or scaled to facilitate comparisons.
#'
#' @return The function invisibly returns a \code{\link[raster]{raster}}, in which each cell contains the `probability of use' score  and produces a plot of this surface.
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
#'              summarise_pr = max,
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
                         "Implementing pf_simplify() with 'summarise_pr' and return = 'archive' specified first is advised."),
                  immediate. = TRUE, call. = FALSE)
        } else {
          warning(paste0("xpf$history[[", t, "]] contains duplicate cells. ",
                         "Did you implement pf_simplify() without specifying 'summarise_pr'?"),
                  immediate. = TRUE, call. = FALSE)
        }
      }
      return(elm)
    })
    dat <- do.call(rbind, dat)
    colnames(dat) <- c("cell_id", "cell_pr")
    dat$path_id <- 1L

    ## pf_path implementation
  } else if(inherits(xpf, "pf_path")){
    dat <- xpf[, c("cell_id", "cell_pr", "path_id")]
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
  n <- length(which(dat$path_id == dat$path_id[1]))
  wt_freq <-
    dat %>%
    dplyr::group_by(.data$path_id, .data$cell_id) %>%
    dplyr::summarise(score = sum(.data$cell_pr)/n) %>%
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
  if(!is.null(add_rasters)) add_rasters$x <- p
  prettyGraphics::pretty_map(x = p,
                             add_rasters = add_rasters,...)
  return(invisible(p))
}


######################################
######################################
#### pf_kud()

#' @title Apply kernel smoothing to particles or paths from a PF algorithm
#' @description This function is a wrapper designed to apply kernel utilisation distribution (KUD) estimation to the outputs of a particle filtering (PF) algorithm. To implement the approach, an (a) \code{\link[flapper]{pf_archive-class}} object from \code{\link[flapper]{pf}} (plus \code{\link[flapper]{pf_simplify}} with the \code{return = "archive"} argument) containing particle histories for connected particles or (b) a \code{\link[flapper]{pf_path-class}} object containing reconstructed paths must be supplied. Using all, or a subset of sampled locations, the function applies KUD smoother(s) via a user-supplied estimation routine (i.e., \code{\link[adehabitatHR]{kernelUD}} or \code{\link[flapper]{kud_around_coastline}}). The function extracts the utilisation distribution(s) as a \code{\link[raster]{raster}}, applies a spatial mask (e.g.  coastline), combines the distributions (if necessary), plots the processed distribution (if specified) and returns this as a a \code{\link[raster]{raster}}.
#'
#' @param xpf A \code{\link[flapper]{pf_archive-class}} object (from \code{\link[flapper]{pf}} plus \code{\link[flapper]{pf_simplify}} with \code{return = "archive"}) or a \code{\link[flapper]{pf_path-class}} object (from \code{\link[flapper]{pf}} plus \code{\link[flapper]{pf_simplify}} with \code{return = "path"}).
#' @param bathy A \code{\link[raster]{raster}} that defines the grid across which \code{\link[flapper]{pf}} was applied. This used to extract cell coordinates and to express KUD(s).
#' @param sample_size (optional) An integer that defines the number of particles to sample from (a) particle histories or (b) each path in \code{xpf} for the estimation. If \code{sample_size = NULL}, all particles are used. If specified, \code{sample_size} particles are sampled from (a) particle histories or (b) each path without replacement in line with their probability. Sampling is implemented to improve estimation speed.
#' @param estimate_ud A function (either \code{\link[adehabitatHR]{kernelUD}} or \code{\link[flapper]{kud_around_coastline}}) that estimates kernel utilisation distributions. The latter option is used when kernels need to be processed to account for barriers to movement that cannot be modelled via \code{\link[adehabitatHR]{kernelUD}}.
#' @param grid,... Arguments passed to \code{estimate_ud} (and ultimately \code{\link[adehabitatHR]{kernelUD}}, where they are defined) to estimate the kernel utilisation distribution. If  \code{\link[flapper]{kud_around_coastline}} is supplied to \code{estimate_ud}, then \code{grid} must be a \code{\link[sp]{SpatialPixelsDataFrame}}. However, note that in all cases, KUD(s) are resampled onto \code{bathy}.
#' @param mask (optional) A spatial mask (see \code{\link[raster]{mask}}).
#' @param plot A logical input that defines whether or not to plot the KUD.
#'
#' @details This function creates smooth KUD representations of particle or path samples from a PF algorithm (see \code{\link[flapper]{pf_plot_map}}).
#'
#' For particle-based implementations, this function is designed to be implemented for the subset of unique, sampled particles that formed continuous paths from the start to the end of the time series (see \code{\link[flapper]{pf_simplify}} and \code{\link[flapper]{pf_plot_map}}). These particles are used for KUD estimation. By default, all particles are used, but \eqn{n =} \code{sample_size} particles can be sampled at random, in line with their probability, if specified, for faster KUD estimation. Selected particles are then used to estimate a KUD by treating samples as `relocations', ultimately via \code{\link[adehabitatHR]{kernelUD}}. This distribution is then processed, plotted and returned.
#'
#' For path-based implementations, this function is designed to be implemented for paths reconstructed by \code{\link[flapper]{pf_simplify}}. As for the particle-based implementation, for each path, all locations, or random sample of \code{sample_size} locations are used to estimate a KUD by treating sampled locations as `relocations'. KUDs are processed and combined across paths into a single, average KUD. The advantage of this approach is that the overall probability of the paths can be incorporated in the estimation procedure as weights when path-specific KUDs are averaged (although that is not yet implemented).
#'
#' For either implementation, unlike maps derived directly from particle or path samples (see \code{\link[flapper]{pf_plot_map}}), KUDs do not account for cell-specific uncertainty.
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
#' ## Implementation based on particles
#' pp <- par(mfrow = c(1, 2))
#' pf_kud(pf_simplify(dat_dcpf_histories, summarise_pr = max, return = "archive"),
#'        bathy = bathy, sample_size = 100,
#'        estimate_ud = kud_around_coastline, grid = grid)
#' ## Implementation based on paths
#' pf_kud(dat_dcpf_paths,
#'        bathy = bathy,
#'        estimate_ud = kud_around_coastline, grid = grid)
#' prettyGraphics::add_sp_path(x = dat_dcpf_paths$cell_x, y = dat_dcpf_paths$cell_y,
#'                             length = 0.01)
#' par(pp)
#'
#' #### Example (2): Implement function using random sampling (for speed)
#' ## Implementation based on particles
#' pp <- par(mfrow = c(1, 2), mar = c(3, 3, 3, 3))
#' pf_kud(pf_simplify(dat_dcpf_histories, return = "archive"),
#'        bathy = bathy,
#'        sample_size = 50, # sample 50 particles overall
#'        estimate_ud = kud_around_coastline, grid = grid)
#' ## Implementation based on paths
#' pf_kud(dat_dcpf_paths,
#'        bathy = bathy,
#'        sample_size = 50, # sample 50 particles per path
#'        estimate_ud = kud_around_coastline, grid = grid)
#' par(pp)
#'
#' @seealso \code{\link[flapper]{pf}}, \code{\link[flapper]{pf_simplify}}, \code{\link[flapper]{pf_plot_map}}, \code{\link[adehabitatHR]{kernelUD}}, \code{\link[flapper]{kud_around_coastline}}, \code{\link[flapper]{eval_by_kud}}
#' @export
#' @author Edward Lavender
#'

pf_kud <- function(xpf,
                   bathy,
                   sample_size = NULL,
                   estimate_ud = adehabitatHR::kernelUD,
                   grid, ...,
                   mask = NULL,
                   plot = TRUE){

  #### Checks
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


  ######################################
  #### pf_archive approach

  if(inherits(xpf, "pf_archive")){

    #### Get cell probabilities and coordinates
    if(xpf$method != "pf_simplify"){
      warning("xpf$method != 'pf_simplify'", immediate. = TRUE, call. = FALSE)
    }
    pairs_df <-
      lapply(xpf$history, function(elm) {
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
                immediate. = TRUE, call. = FALSE)
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
      proj4string = crs)
    pairs_ud <- estimate_ud(xy = pairs_xy_spdf, grid = grid,...)

    #### Process KUD
    pairs_ud <- raster::raster(pairs_ud[[1]])
    if(!is.null(mask)) pairs_ud <- raster::mask(pairs_ud, mask)


    ######################################
    #### pf_path approach

  } else if(inherits(xpf, "pf_path")){

    #### Get cell coordinates
    pairs_df <- data.frame(xpf)
    pairs_xy <- xpf[, c("cell_x", "cell_y")]

    #### Sub sample
    if(!is.null(sample_size)){
      if(sample_size > length(which(pairs_df$path_id == pairs_df$path_id[1]))){
        warning(paste0("'sample_size' (n = ", sample_size,
                       ") exceeds the number of sampled locations (n = ", nrow(pairs_xy),
                       "): implementing estimation using all sampled locations."),
                immediate. = TRUE, call. = FALSE)
        sample_size <- NULL
      }
      if(!is.null(sample_size)){
        pairs_df$index <- 1:nrow(pairs_df)
        pairs_df <-
          pairs_df %>%
          dplyr::group_by(.data$path_id, .groups = "drop_last") %>%
          dplyr::slice_sample(n = sample_size, weight_by = .data$cell_pr) %>%
          data.frame()
        pairs_xy <- pairs_xy[pairs_df$index, ]
      }
    }

    #### Implement KUD estimation for each path
    pairs_xy_spdf <- sp::SpatialPointsDataFrame(
      pairs_xy,
      data = data.frame(ID = factor(pairs_df$path_id)),
      proj4string = crs)
    pairs_ud <- estimate_ud(xy = pairs_xy_spdf, grid = grid,...)

    #### Combine KUDs across paths
    ## Convert KUDs to rasterStack
    blank <- raster::setValues(bathy, 0)
    pairs_ud <- lapply(pairs_ud, function(ud){
      ud <- raster::raster(ud)
      ud <- raster::resample(ud, blank)
      if(!is.null(mask)) ud <- raster::mask(ud, mask)
      return(ud)
    })
    pairs_ud <- raster::stack(pairs_ud)
    ## Define weights associated with each path
    # wts <- pf_loglik(paths = xpf)
    # wts$wt <- (1 - wts$delta/sum(wts$delta)) # non linear
    # wts <- wts[order(wts$path_id), ]
    ## Summarise KUDs across paths
    # pairs_ud <- raster::weighted.mean(pairs_ud, wts$wt)
    pairs_ud <- raster::calc(pairs_ud, mean)
    ## Renormalise KUDs
    pairs_ud <- pairs_ud/raster::cellStats(pairs_ud, "sum")
  }

  #### Visualise KUD
  ext <- raster::extent(pairs_ud)
  if(plot) prettyGraphics::pretty_map(add_rasters = list(x = pairs_ud),
                                      xlim = ext[1:2], ylim = ext[3:4])

  #### Return KUD
  return(pairs_ud)
}
