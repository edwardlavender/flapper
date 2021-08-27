######################################
######################################
#### kud_around_coastline()

#' @title Process a kernel utilisation distribution around a barrier
#' @description Given an animal movement path over a gridded surface, this function estimates a `raw' kernel utilisation distribution (from \code{\link[adehabitatHR]{kernelUD}}) and then processes the distribution to account for barriers to movement, such as coastline. To implement the function, the movement path(s) should be supplied as a SpatialPointsDataFrame and the grid over which estimation is implemented as a SpatialPixelsDataFrame with values 0 and 1 defining unsuitable and suitable habitat respectively.
#'
#' @param xy A \code{\link[sp]{SpatialPointsDataFrame}} object that defines the movement path(s). This should contain a column that defines the individual ID(s) as a factor.
#' @param grid A \code{\link[sp]{SpatialPixelsDataFrame}} object that defines the grid over which estimation is implemented and binary habitat suitability (0, unsuitable; or 1, suitable).
#' @param ... Additional arguments passed to \code{\link[adehabitatHR]{kernelUD}}.
#'
#' @details Utilisation distributions (UDs) are bivariate probability distributions that describe the probability (density) of locating an individual in any given area at a randomly chosen time. These can be estimated using the \code{\link[adehabitatHR]{kernelUD}} function. The algorithms implemented by this function can incorporate simple barriers, but restrictions on the shapes of barriers mean that in many real-world settings (e.g., in areas with complex coastline) barriers cannot be implemented. As a result, a pragmatic (if somewhat unsatisfactory) approach is to post-process the raw utilisation distribution by removing areas in which movement is impossible and then re-normalise the distribution (so that probabilities sum to one). This function achieves this by implementing the estimation over a grid, which defines whether (1) or not (0) an area is `habitat'. After the estimation of the raw UD across the grid, probability density scores are combined (multiplied) with the habitat suitability score (0, 1) and then renormalised (by dividing by the total score across suitable areas).
#'
#' @return The function returns an object of class `estUDm'. This is a list, with one component per animal, of \code{\link[adehabitatHR]{estUD-class}} objects. The `h' slot of the output (\code{output@h}) has been modified so that the method (`meth') is given as `specified'.
#'
#' @examples
#' #### Set up
#' ## (1) Simulate path for which to compute UD
#' # Focus on a sample of the marine environment off Oban, West Scotland
#' sea <- invert_poly(dat_coast)
#' sea <- raster::crop(sea, update_extent(raster::extent(sea), x_shift = -2500))
#' prettyGraphics::pretty_map(add_polys = list(x = sea, col = "skyblue"))
#' # Simulate path
#' n <- 1000
#' path_ls <- sim_path_sa(n = n,
#'                        area = sea,
#'                        sim_step = function(...) stats::rgamma(1, shape = 20, scale = 20),
#'                        seed = 1,
#'                        plot = FALSE
#'                        )
#' prettyGraphics::add_sp_path(path_ls$xy_mat,
#'                             col = viridis::viridis(n),
#'                             length = 0.02)
#' ## (2) Define path as a SpatialPointsDataFrame (SpatialPoints is not allowed)
#' path <- sp::SpatialPointsDataFrame(
#'   path_ls$xy_mat,
#'   data = data.frame(ID = factor(rep(1, nrow(path_ls$xy_mat)))),
#'                     proj4string = raster::crs(dat_coast))
#' ## (3) Define grid over which to implement estimation
#' # ... The grid needs to be sufficiently small to capture the coastline
#' # ... reasonably while being large enough to enable calculation
#' # ... of the home range.
#' r <- raster::raster(raster::extent(dat_coast), nrows = 100, ncols = 100)
#' raster::values(r) <- 0
#' r <- raster::mask(r, dat_coast, updatevalue = 1)
#' habitat <- methods::as(r, "SpatialPixelsDataFrame")
#' sp::plot(habitat)
#'
#' #### Example (1) Implement estimation and processing
#' ## Estimate raw UD
#' ud_raw <- adehabitatHR::kernelUD(xy = path, grid = habitat)
#' # Object is of class estUDm, which is a list of estUD objects
#' # The outputs for each animal can be accessed by indexing
#' ud_raw[[1]]
#' # Check smoothing parameters
#' ud_raw[[1]]@h
#' ## Estimate raw UD and post-process
#' ud_pro <- kud_around_coastline(xy = path, grid = habitat)
#' # The same type of object is returned
#' ud_pro[[1]]
#' # Smoothing parameters have been modified
#' ud_pro[[1]]@h
#' ## Compare plots
#' # ... Notice that the processed version doesn't 'bleed' onto land
#' # ... and the scale differs due to the re-normalisation
#' ud_raw_r <- raster::raster(ud_raw[[1]])
#' ud_pro_r <- raster::raster(ud_pro[[1]])
#' pp <- graphics::par(mfrow = c(1, 2))
#' prettyGraphics::pretty_map(add_rasters = list(x = ud_raw_r),
#'                            add_polys = list(x = dat_coast),
#'                            add_paths = list(x = path,
#'                                             col = viridis::viridis(n),
#'                                             lwd = 0.25,
#'                                             length = 0.02)
#'                            )
#' prettyGraphics::pretty_map(add_rasters = list(x = ud_pro_r),
#'                            add_polys = list(x = dat_coast),
#'                            add_paths = list(x = path,
#'                                             col = viridis::viridis(n),
#'                                             lwd = 0.25,
#'                                             length = 0.02)
#'                            )
#' graphics::par(pp)
#'
#' #### Further analysis can be implemented as usual
#' # For example we can compute the home range contours from the UD
#' # ... using adehabitatHR::getvolumeUD(). In practice, this converts from the
#' # ... probability density scale to a more intuitive % home range scale.
#' # Get volume
#' vol_raw <- adehabitatHR::getvolumeUD(ud_raw, standardize = TRUE)
#' vol_pro <- adehabitatHR::getvolumeUD(ud_pro, standardize = TRUE)
#' # Get contours
#' ver_raw <- adehabitatHR::getverticeshr(ud_raw[[1]], standardize = TRUE)
#' ver_pro <- adehabitatHR::getverticeshr(ud_pro[[1]], standardize = TRUE)
#' # Rasterise
#' vol_raw_r <- raster::raster(vol_raw[[1]])
#' vol_pro_r <- raster::raster(vol_pro[[1]])
#' # For neatness on the plot, it is convenient to exclude areas beyond 95 %
#' vol_raw_r[vol_raw_r[] > 95] <- NA
#' vol_pro_r[vol_pro_r[] > 95] <- NA
#' # Plot
#' pp <- graphics::par(mfrow = c(1, 2))
#' prettyGraphics::pretty_map(add_rasters = list(x = vol_raw_r),
#'                            add_polys = list(list(x = dat_coast),
#'                                             list(x = ver_raw,
#'                                                  border = "blue",
#'                                                  lwd = 2)),
#'                            add_paths = list(x = path,
#'                                             col = viridis::viridis(n),
#'                                             lwd = 0.25,
#'                                             length = 0.02)
#'                            )
#' prettyGraphics::pretty_map(add_rasters = list(x = vol_pro_r),
#'                            add_polys = list(list(x = dat_coast),
#'                                             list(x = ver_pro,
#'                                                  border = "blue",
#'                                                  lwd = 2)),
#'                            add_paths = list(x = path,
#'                                             col = viridis::viridis(n),
#'                                             lwd = 0.25,
#'                                             length = 0.02)
#'                            )
#' graphics::par(pp)
#'
#' @source This forum is a useful resource: http://r-sig-geo.2731867.n2.nabble.com/Walruses-and-adehabitatHR-class-estUDm-exclusion-of-non-habitat-pixels-and-summary-over-all-animals-td6497315.html.
#' @author Edward Lavender
#' @export

kud_around_coastline <- function(xy, grid,...) {

  ## Define raw UD
  ud <- adehabitatHR::kernelUD(xy, grid = grid,...)
  names_ud <- names(ud)

  ## Checks
  if(!inherits(xy, "SpatialPointsDataFrame")){
    message("'xy' is not a SpatialPointsDataFrame: returning raw UD as post-processing cannot be implemented.")
    return(ud)
  }
  if(!inherits(grid, "SpatialPointsDataFrame")){
    message("'grid' is not a SpatialPixelsDataFrame: returning raw UD as post-processing cannot be implemented.")
    return(ud)
  }

  ## Convert it to SpatialPixelsDataFrame:
  ud <- adehabitatHR::estUDm2spixdf(ud) # ud as SpatialPointsDataFrame required for estUDm2spixdf function.
  ## Convert to SpatialGridDataFrames
  sp::fullgrid(ud) <- TRUE
  sp::fullgrid(grid) <- TRUE

  ## Exclude areas of 'unsuitable habitat' and renormalise
  # ... by multiplying the value of the UD by the habitat (0, 1) and
  # ... dividing by the total value of the UD
  ud_vals_renorm <- lapply(1:ncol(ud), function(i) {
    ud[[i]] * grid[[1]] / sum(ud[[i]] * grid[[1]])
  })
  ud_vals_renorm <- as.data.frame(ud_vals_renorm)
  names(ud_vals_renorm) <- "ud"

  ## Add to UD as data slot
  ud@data <- ud_vals_renorm
  # sum(ud@data) == 1

  ## Coerce object back to estUD-class
  # This is simplify a SpatialPixelsDataFrame with a slot for h and a slot for vol
  sp::fullgrid(ud) <- FALSE
  ud_ade <- lapply(1:ncol(ud), function(i) {
    so <- methods::new("estUD", ud[, i])
    so@h <- list(h = 0, meth = "specified")
    so@vol <- FALSE
    return(so)
  })
  names(ud_ade) <- names_ud
  class(ud_ade) <- "estUDm"

  ## Return
  return(ud_ade)
}
