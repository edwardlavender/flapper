#############################################
#############################################
#### eval_by_kud()

#' @title Evaluate movement path estimates using kernel utilisation distributions
#' @description This function provides a qualitative evaluation of a movement-path estimation procedure when the true movement path is known by graphically comparing the two paths and the resultant patterns of space use that emerge from these using kernel utilisation distributions. To implement the function, the simulated and estimated paths need to be provided as SpatialPointsDataFrames. For each path, the function estimates the resultant home range under a kernel utilisation distribution model and then produces up to four plots that show the two paths and the two estimated home ranges. When the estimation procedure performs well, the patterns of space use inferred from the simulated and estimated paths should be similar. In contrast, differences between these patterns can help to reveal the circumstances under which the movement-path estimation procedure performs less well.
#'
#' @param path_sim,path_est \code{\link[sp]{SpatialPointsDataFrame}} objects that the coordinates of a simulated movement path and estimated movement path respectively. These should contain a column that defines the individual's ID as a factor. The coordinate reference system is taken from \code{path_sim} and should be defined, if applicable.
#' @param estimate_ud A function (either \code{\link[adehabitatHR]{kernelUD}} or \code{\link[flapper]{kud_around_coastline}}) that estimates kernel utilisation distributions. The latter option is used when kernels need to be processed to account for barriers to movement that cannot be modelled via \code{\link[adehabitatHR]{kernelUD}}.
#' @param grid,h,hlim,kern,extent,boundary Arguments passed to \code{estimate_ud} (and ultimately \code{\link[adehabitatHR]{kernelUD}}, where they are defined) to estimate kernel utilisation distributions. If  \code{\link[flapper]{kud_around_coastline}} is supplied to \code{estimate_ud}, then \code{grid} must be a \code{\link[sp]{SpatialPixelsDataFrame}}.
#' @param process (optional) A function, applied to a \code{\link[raster]{raster}}, to process the home range raster (extracted from from \code{\link[adehabitatHR]{getvolumeUD}}). For example, \code{function(r) { r[r > 95] <- NA; return(r) }} can help to produce prettier plots by masking areas beyond the 95 percent contour.
#' @param plot An integer vector (\code{1:4L}) that defines which plots to produce. \code{1} and \code{2} plot the receiver array and the simulated or estimated path respectively. In contrast, \code{3} and \code{4} plot the receiver array and the kernel utilisation distribution based on the simulated or estimated path respectively.
#' @param array A named list that defines the properties of an array (e.g., the `array' element of \code{\link[flapper]{sim_array}}). At a minimum, this must contain an `area' element that defines the area as a \code{\link[sp]{SpatialPolygons-class}} object. To add receivers, land and sea to the plot(s), `xy', `land' and `sea' elements are also required (see \code{\link[flapper]{sim_array}}).
#' @param add_land,add_sea,add_receivers,add_path_sim,add_path_est,add_vol_sim,add_vol_est (optional) Named lists of arguments that customise the appearance of the land and sea, receivers, the simulated and estimated paths and the kernel volumes for the simulated and estimated paths respectively. An empty list (\code{list()}) will implement default graphical options. If \code{add_receivers}, \code{add_land} and \code{add_sea} are specified, then \code{array} must contain `xy', `land' and `sea' elements as described above. Alternatively, \code{add_receivers}, \code{add_land} and \code{add_sea} can be \code{NULL} to suppress their addition to plots, but the other arguments cannot be \code{NULL}.
#' @param one_page A logical variable that defines whether or not to produce all plots on a single page.
#' @param verbose A logical variable that defines whether or not to print messages to the console to relay function progress.
#' @param ... Additional arguments, passed to \code{\link[prettyGraphics]{pretty_map}}, that affect all plots.
#'
#' @details This function can be combined with \code{\link[flapper]{sim_array}}, \code{\link[flapper]{sim_path_*}} functions, \code{\link[flapper]{sim_detections}} and space use algorithms (e.g., \code{\link[flapper]{coa}}) to evaluate the relative performance of different approaches for the inference of patterns of space use under different array designs, movement models and/or detection models.
#'
#' @return The function returns up to four plots that provide a qualitative evaluation of an estimation procedure that infers patterns of space use from observations.
#'
#' @examples
#' #### Steps
#' # We need to simulate an array, a movement path through this arra
#' # ... and detections arising at receivers given this path. Then, we need to
#' # ... implement a space use algorithm (e.g., coa()) to estimate the movement
#' # ... 'path' from the observed data (detections). With the simulated and
#' # ... estimated movement path, we can then compare how well the estimated
#' # ... path recovers the 'true' pattern of space use by passing these
#' # ... paths to the eval_by_kud() function.
#'
#' #### Step (1) Simulate an array in an area
#' # We will pack the area with receivers to generate lots of detections
#' array_ls <- sim_array(boundaries = raster::extent(dat_coast),
#'                       coastline = dat_coast,
#'                       n_receivers = 1000,
#'                       arrangement = "regular",
#'                       seed = 1)
#' raster::lines(dat_coast)
#' array <- array_ls$array
#'
#' #### Step (2) Simulate a movement path in this area
#' n <- 500
#' path_ls <- sim_path_sa(n = n,
#'                        sim_step = function(...) stats::rgamma(1, shape = 25, scale = 25),
#'                        area = array$sea,
#'                        seed = 1,
#'                        plot = FALSE
#'                        )
#' prettyGraphics::add_sp_path(path_ls$xy_mat, col = viridis::viridis(n), length = 0.02)
#'
#' #### Step (3) Simulate and aggregate detections over some delta_t interval
#' # E.g., we can imagine simulating a movement path at 60 s resolution
#' # ... and determining whether or not the individual is detected in that interval
#' # ... and then aggregating detections over some interval e.g., X mins
#' # ... (In reality, we would need to consider this choice carefully).
#' det_ls <- sim_detections(path = path_ls$xy_mat,
#'                          xy = sp::coordinates(array$xy),
#'                          calc_detection_pr = function(dist) ifelse(dist < 700, 1, 0),
#'                          delta_t = 10)
#'
#' str(det_ls)
#'
#' #### Step (4) Implement estimation proceedure e.g., calculate COAs
#' coas <- coa(mat = det_ls$agg$det_mat, xy = sp::coordinates(array$xy))
#'
#' #### Step (5) Define simulated and estimated paths as SpatialPointsDataFrames
#' path_true <- path_ls$xy_mat
#' path_true <-
#'   sp::SpatialPointsDataFrame(path_true,
#'                              data = data.frame(ID = factor(rep(1, nrow(path_true)))),
#'                              proj4string = raster::crs(dat_coast))
#' path_coa <- coas
#' path_coa <-
#'   sp::SpatialPointsDataFrame(path_coa,
#'                              data = data.frame(ID = factor(rep(1, nrow(path_coa)))),
#'                              proj4string = raster::crs(dat_coast))
#'
#' #### Example (1): Implement algorithm using default options
#' eval_by_kud(path_sim = path_true,
#'             path_est = path_coa,
#'             grid = 60,
#'             array = array)
#'
#' #### Example (2): Account for coastline via kud_around_coastline()
#' # Define grid of habitat/non habitat
#' r <- raster::raster(raster::extent(dat_coast), nrows = 100, ncols = 100)
#' raster::values(r) <- 0
#' r <- raster::mask(r, dat_coast, updatevalue = 1)
#' habitat <- methods::as(r, "SpatialPixelsDataFrame")
#' # Implement algorithm
#' eval_by_kud(path_sim = path_true,
#'             path_est = path_coa,
#'             estimate_ud = kud_around_coastline,
#'             grid = habitat,
#'             array = array)
#'
#' #### Example (3): Plot customisation options
#' # Use add_* lists to customise main plotting features
#' # ... such as receivers and land (but also the movement paths etc.)
#' eval_by_kud(path_sim = path_true,
#'             path_est = path_coa,
#'             estimate_ud = kud_around_coastline,
#'             grid = habitat,
#'             array = array,
#'             add_receivers = list(pch = "."),
#'             add_land = list(col = "darkgreen"))
#' # Use 'process' to focus on contours of the home range of interest
#' process <- function(r) { r[r > 95] <- NA; return(r) }
#' eval_by_kud(path_sim = path_true,
#'             path_est = path_coa,
#'             estimate_ud = kud_around_coastline,
#'             grid = habitat,
#'             process = process,
#'             array = array,
#'             add_receivers = list(pch = "."),
#'             add_land = list(col = "darkgreen"))
#' # Pass other arguments to via ... to pretty map
#' eval_by_kud(path_sim = path_true,
#'             path_est = path_coa,
#'             estimate_ud = kud_around_coastline,
#'             grid = habitat,
#'             process = process,
#'             array = array,
#'             add_receivers = list(pch = "."),
#'             add_land = list(col = "darkgreen"),
#'             xlab = "x", ylab = "y")
#'
#' #### Results
#' # In this case, with a regularly spaced and dense array, the COA metric
#' # ... provides good inferences regarding patterns of space use. However,
#' # ... this will not always be the case.
#'
#' @author Edward Lavender
#' @export

eval_by_kud <-
  function(path_sim,
           path_est,
           estimate_ud = adehabitatHR::kernelUD,
           grid, h = "href",
           hlim = c(0.1, 1.5), kern = c("bivnorm", "epa"), extent = 1, boundary = NULL,
           process = NULL,
           plot = 1:4L,
           array,
           add_receivers = list(),
           add_land = NULL,
           add_sea = NULL,
           add_path_sim = list(col = viridis::viridis(nrow(path_sim)), length = 0.02),
           add_path_est = add_path_sim,
           add_vol_sim = list(), add_vol_est = list(),
           one_page = TRUE,
           verbose = TRUE,...){

  #### Get CRS
  crs_area <- raster::crs(path_sim)

  #### Estimate UDs and volumes

  ## Estimate KUDs for path_sim
  # derive the raw UD:
  path_sim_ud <- estimate_ud(xy = path_sim, h = h, grid = grid, hlim = hlim, kern = kern, extent = extent, boundary = boundary)
  if(inherits(path_sim_ud, "estUDm")) path_sim_ud <- path_sim_ud[[1]]
  raster::crs(path_sim_ud) <- crs_area
  # define the home range %:
  path_sim_vol <- adehabitatHR::getvolumeUD(path_sim_ud, standardize = TRUE)
  path_sim_vol <- raster::raster(path_sim_vol)
  if(!is.null(process)) path_sim_vol <- process(path_sim_vol)
  raster::crs(path_sim_vol) <- crs_area

  ## Estimate KUDs for path_est
  path_est_ud <- adehabitatHR::kernelUD(xy = path_est, h = h, grid = grid, hlim = hlim, kern = kern, extent = extent, boundary = boundary)
  if(inherits(path_est_ud, "estUDm")) path_est_ud <- path_est_ud[[1]]
  raster::crs(path_est_ud) <- crs_area
  # define the home range %:
  path_est_vol <- adehabitatHR::getvolumeUD(path_est_ud, standardize = TRUE)
  path_est_vol <- raster::raster(path_est_vol)
  if(!is.null(process)) path_est_vol <- process(path_est_vol)
  raster::crs(path_est_vol) <- crs_area

  #### Set up plotting param
  if(one_page) pp <- graphics::par(mfrow = prettyGraphics::par_mf(length(plot)))
  if(!is.null(add_land)) add_land$x <- array$land
  if(!is.null(add_sea)) add_sea$x <- array$sea
  add_polys <- list(add_land, add_sea)
  add_polys <- compact(add_polys)
  if(length(add_polys) == 0L) add_polys <- NULL
  if(!is.null(add_receivers)) add_receivers$x <- array$xy

  #### Plot 1: receiver array and true trajectory
  if(1L %in% plot) {

    add_path_sim$x <- path_sim
    prettyGraphics::pretty_map(array$area,
                               add_polys = add_polys,
                               add_paths = add_path_sim,
                               add_points = add_receivers,
                               verbose = FALSE,...)

  }

  #### Plot 2: receiver array and path_est (e.g., COAs)
  if(2L %in% plot) {
    add_path_est$x <- path_est
    prettyGraphics::pretty_map(array$area,
                               add_polys = add_polys,
                               add_paths = add_path_est,
                               add_points = add_receivers,
                               verbose = FALSE,...)
  }

  #### Plot 3: receiver array and KUD based on true trajectory
  if(3L %in% plot) {
    add_vol_sim$x <- path_sim_vol
    prettyGraphics::pretty_map(array$area,
                               add_rasters = add_vol_sim,
                               add_polys = add_polys,
                               verbose = FALSE,...)
  }

  #### Plot 4: receiver array and KUD based on path_est (e.g., COAs)
  if(4L %in% plot) {
    add_vol_est$x <- path_est_vol
    prettyGraphics::pretty_map(array$area,
                               add_rasters = add_vol_est,
                               add_polys = add_polys,
                               verbose = FALSE,...)
  }

  #### Close plotting window
  graphics::par(pp)

}
