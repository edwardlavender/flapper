#' @title Smooth POU maps
#' @description This function smooths proportion-of-use (POU) maps (from \code{\link[flapper]{pf_plot_map}}) by applying kernel utilisation distribution (KUD) estimation. Depending on the implementation, following optional initial time trials, using a subset, all or an expanded sample of POU locations, the function applies a KUD smoother via a user-supplied estimation routine (i.e., \code{\link[adehabitatHR]{kernelUD}} or \code{\link[flapper]{kud_around_coastline}}). The function extracts the KUD as a \code{\link[raster]{raster}}, applies a spatial mask (e.g., coastline), plots the processed KUD (if specified) and returns this as a \code{\link[raster]{raster}}.
#'
#' @param xpf A POU \code{\link[raster]{raster}} object (from \code{\link[flapper]{pf_plot_map}}).
#' @param sample_size (optional) An integer expansion factor for the number of locations used for KUD estimation. If supplied, \eqn{n} locations are randomly sampled from \code{xpf} with replacement in line with their probability, where \eqn{n = n_{pou} \times sample_size} and \eqn{n_{pou}} is the number of non-zero POU scores. This resampling approach avoids treating locations as `relocations'.
#' @param estimate_ud A function (either \code{\link[adehabitatHR]{kernelUD}} or \code{\link[flapper]{kud_around_coastline}}) that estimates the KUD.
#' @param grid,... Arguments passed to \code{estimate_ud} (and ultimately \code{\link[adehabitatHR]{kernelUD}}, where they are defined) to estimate the KUD. If \code{\link[flapper]{kud_around_coastline}} is supplied to \code{estimate_ud}, then \code{grid} must be a \code{\link[sp]{SpatialPixelsDataFrame}}. The resultant KUD is resampled onto \code{xpf}.
#' @param trial_cells,trial_grids (Optional) Lists that define the numbers of locations (cells) and the grids used for time trials. If either \code{trial_cells} or \code{trial_grids} is supplied, the function implements \code{estimate_ud} for small (trial) numbers of cells and any grid(s) specified. If \code{trial_cells} is supplied, but not \code{trial_grids}, then \code{grid} is used for estimation. If \code{trial_grids} is supplied but not \code{trial_cells}, then time trials are implemented for each grid with 10, 50 and 100 locations. For each grid, the linear regression of the time required to estimate the KUD is plotted against the number of locations and used to predict the time required to fit the KUD to all locations. The user is then asked whether or not to continue with estimation across all locations.
#' @param mask (optional) A spatial mask (see \code{\link[raster]{mask}}).
#' @param plot A logical input that defines whether or not to plot the KUD.
#' @param verbose A logical input that defines whether or not to print messages to the console to monitor function progress.
#'
#' @details For computational efficiency, it may be necessary to aggregate (and renormalise) POU scores across the \code{grid} used for estimation before implementing this function.
#'
#' @return The function (a) plots the outcomes of time trials (if requested), (b) estimates and plots a KUD (if requested) and (c) returns a \code{\link[raster]{raster}} of the KUD.
#'
#' @examples
#' #### Define POU map for examples
#' out_dcpf_s <- pf_simplify(dat_dcpf_histories,
#'   summarise_pr = TRUE,
#'   return = "archive"
#' )
#' out_dcpf_pou <- pf_plot_map(out_dcpf_s, dat_dcpf_histories$args$bathy)
#'
#' #### Example (1): Implement function using default options
#' pf_kud(xpf = out_dcpf_pou, grid = 10)
#'
#' #### Example (2): Implement function using resampling
#' pf_kud(xpf = out_dcpf_pou, sample_size = 100, grid = 10)
#'
#' #### Example (3): Implement time trials
#' if (interactive()) {
#'   # Implement time trials for specified numbers of cells
#'   pf_kud(
#'     xpf = out_dcpf_pou,
#'     sample_size = 100,
#'     grid = 60,
#'     trial_cells = list(10, 100, 1000, 10000)
#'   )
#'   # Implement time trials for specified grids
#'   pf_kud(
#'     xpf = out_dcpf_pou,
#'     sample_size = 100,
#'     grid = 180,
#'     trial_grids = list(60, 120, 180)
#'   )
#'   # Implement time trials for specified numbers of cells and grids
#'   pf_kud(
#'     xpf = out_dcpf_pou,
#'     sample_size = 100,
#'     grid = 10,
#'     trial_cells = list(10, 100, 1000, 10000),
#'     trial_grids = list(10, 30, 60)
#'   )
#' }
#'
#' #### Example (4): Force alignment between POU scores and grids for speed
#' # Example with numeric grid
#' out_dcpf_pou_agg <-
#'   raster::aggregate(out_dcpf_pou, fact = 2)
#' out_dcpf_pou_agg <-
#'   out_dcpf_pou_agg / raster::cellStats(out_dcpf_pou_agg, "sum")
#' grid <- raster::res(out_dcpf_pou_agg)[1]
#' pf_kud(out_dcpf_pou_agg, grid = grid)
#' # Example with SpatialPixels grid
#' grid <- kud_habitat(out_dcpf_pou_agg)
#' pf_kud(out_dcpf_pou_agg, grid = grid)
#'
#' @author Edward Lavender
#' @export

pf_kud <- function(xpf,
                   sample_size = NULL,
                   estimate_ud = adehabitatHR::kernelUD,
                   grid, ...,
                   trial_cells = list(),
                   trial_grids = list(),
                   mask = NULL,
                   plot = TRUE,
                   verbose = TRUE) {
  #### Checks
  cat_to_console <- function(..., show = verbose) {
    if (show) cat(paste(..., "\n"))
  }
  t_onset <- Sys.time()
  cat_to_console(paste0("flapper::pf_kud called (@ ", t_onset, ")..."))
  cat_to_console("... Setting up function...")
  crs_xpf <- crs_grid <- crs_mask <- NULL
  crs_xpf <- raster::crs(xpf)
  if (!is.null(grid)) if (!inherits(grid, c("numeric", "integer"))) crs_grid <- raster::crs(grid)
  if (!is.null(mask)) crs_mask <- raster::crs(mask)
  crs <- list(crs_xpf, crs_grid, crs_mask)
  crs <- compact(crs)
  if (!length(unique(crs)) == 1L) {
    warning("The CRS(s) of spatial layer(s) ('xpf', 'grid' and 'mask', if applicable) are not identical.",
      immediate. = TRUE, call. = FALSE
    )
  }
  crs_spp <- sapply(crs, function(x) !is.na(x))
  if (any(crs_spp)) {
    crs <- crs[[which(crs_spp)[1]]]
  } else {
    crs <- sp::CRS(as.character(NA))
  }
  message("CRS taken as: '", crs, "'.")
  blank <- raster::setValues(xpf, 0)

  #### Get POU scores
  cat_to_console("... Getting POU scores...")
  pou <- data.frame(cell = raster::Which(xpf > 0, cells = TRUE, na.rm = TRUE))
  pou$score <- raster::extract(xpf, pou$cell)
  pou$cell_x <- raster::xFromCell(xpf, pou$cell)
  pou$cell_y <- raster::yFromCell(xpf, pou$cell)
  cat_to_console("... .... POU scores extracted for", nrow(pou), "locations...")

  #### Sample locations according to their probability
  if (!is.null(sample_size)) {
    cat_to_console("... Sampling cells...")
    pou <-
      pou %>%
      dplyr::slice_sample(
        n = nrow(pou) * sample_size,
        weight_by = pou$score, replace = TRUE
      ) %>%
      data.frame()
    cat_to_console("... ... POU locations expanded to", nrow(pou), "locations...")
  }

  #### Build SPDF for KUD estimation
  cat_to_console("... Building SpatialPointsDataFrame...")
  spdf <- sp::SpatialPointsDataFrame(
    pou[, c("cell_x", "cell_y")],
    data = data.frame(ID = factor(rep(1, nrow(pou)))),
    proj4string = crs
  )

  #### Implement time trials
  if ((length(trial_cells) > 0L) | (length(trial_grids) > 0L)) {
    #### Define default parameters
    cat_to_console("... Implementing time trials...")
    if (length(trial_cells) == 0L) trial_cells <- c(10, 50, 100)
    trial_cells <- unlist(trial_cells)
    trial_cells <- trial_cells[trial_cells > 5 & trial_cells <= nrow(pou)]
    if (length(trial_grids) == 0L) trial_grids <- list(grid)

    #### Estimate times for UD fitting
    trial_time <- lapply(trial_cells, function(trial_cell) {
      trial_time_by_grid <- sapply(trial_grids, function(trial_grid) {
        t1 <- Sys.time()
        invisible(estimate_ud(xy = spdf[1:trial_cell, ], grid = trial_grid, ...))
        t2 <- Sys.time()
        return(as.numeric(difftime(t2, t1, units = "mins")))
      })
      return(trial_time_by_grid)
    }) %>% unlist()
    trial_time <-
      expand.grid(
        cell = trial_cells,
        grid = 1:length(trial_grids)
      ) %>%
      dplyr::arrange(.data$cell) %>%
      dplyr::mutate(time = trial_time) %>%
      data.frame()
    po <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(po), add = TRUE)
    pp <- graphics::par(mfrow = prettyGraphics::par_mf(max(trial_time$grid)))
    lapply(split(trial_time, trial_time$grid), function(tt) {
      mod <- stats::lm(time ~ cell, data = tt)
      pred <- stats::predict(mod, newdata = data.frame(cell = nrow(pou)), se.fit = TRUE)
      ci <- prettyGraphics::list_CIs(pred)
      if (is.na(pred$se.fit)) {
        title <- paste0("Grid (", tt$grid[1], "); Pred (", round(pred$fit, 2), " mins)")
      } else {
        title <- paste0(
          "Grid (", tt$grid[1], "); Pred (", round(pred$fit, 2),
          " [", round(ci$lowerCI, 2), "-", round(ci$upperCI, 2), "] mins)"
        )
      }
      prettyGraphics::pretty_predictions_1d(mod,
        add_main = list(text = title)
      )
    })
    graphics::par(pp)
    readline("Press [Enter] to continue or [Esc] to quit...")
  }

  #### Estimate UD
  cat_to_console("... Implementing KUD estimation based on", nrow(pou), "cells...")
  ud <- estimate_ud(xy = spdf, grid = grid, ...)

  #### Process KUD
  cat_to_console("... Processing KUD(s)...")
  if (inherits(ud, "estUDm")) ud <- ud[[1]]
  if (!inherits(ud, "RasterLayer")) {
    ud <- raster::raster(ud)
  }
  ud <- raster::resample(ud, xpf)
  if (!is.null(mask)) ud <- raster::mask(ud, mask)

  #### Renormalise KUDs
  ud <- ud / raster::cellStats(ud, "sum")

  #### Visualise KUD
  if (plot) {
    cat_to_console("... Plotting KUD...")
    ext <- raster::extent(xpf)
    if (plot) {
      prettyGraphics::pretty_map(
        add_rasters = list(x = ud),
        xlim = ext[1:2], ylim = ext[3:4]
      )
    }
  }

  #### Return KUD
  return(invisible(ud))
}
