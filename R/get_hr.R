#' @title Get animal `home ranges'
#' @description These functions extract `home range' estimates from \code{SpatRaster} (or \code{\link[raster]{raster}}) objects that describe the intensity of movements within an area (from example from \code{\link[flapper]{pf_kud}}).
#'
#' @param x A \code{SpatRaster} (or a \code{\link[raster]{raster}}) of the utilisation distribution (UD).
#' @param prop For \code{\link[flapper]{get_hr_prop}}, \code{prop} is a number that defines the home range proportion.
#' @param plot A logical variable that defines whether or not to plot the home range.
#' @param add_raster,add_contour,... Plot customisation options. \code{add_raster} is a named list of arguments to customise the home range surface that is passed to \code{\link[prettyGraphics]{pretty_map}}. \code{...} are additional arguments passed to \code{\link[prettyGraphics]{pretty_map}}. \code{add_contour} is a named list of arguments passed to \code{\link[raster]{contour}} that is called afterwards.
#'
#' @details Animal home ranges are widely quantified as `the smallest subregion that accounts for a specified proportion, \emph{p}, of [the animal's] total distribution' (Jennrich and Turner 1969, page 232). In line with this approach, \code{\link[flapper]{get_hr_prop}} extracts the region within a frequency distribution of space use (i.e., UD) that is enclosed by a specified proportion (\code{prop}) contour. Following the most widely used adopted conventions, three additional wrapper functions facilitate the extraction of core, home and full ranges:
#' \itemize{
#'   \item \code{\link[flapper]{get_hr_core}} extracts the `core range' as the region enclosed by the 50 percent contour of the UD (\code{prop = 0.50});
#'   \item \code{\link[flapper]{get_hr_home}} extracts the `home range' as the 95 percent contour of the UD (\code{prop = 0.95});
#'   \item \code{\link[flapper]{get_hr_full}} extracts the `full' range as the boundaries of the UD (\code{prop = 1.00});
#' }
#'
#' These functions are simple wrappers for \code{\link[spatialEco]{raster.vol}}. They differ from functions in the \code{adehabitatHR} package (namely \code{\link[adehabitatHR]{getverticeshr}}) in that they are designed to input and output \code{\link[raster]{raster}} objects.
#'
#' @return The functions return a \code{\link[raster]{raster}}. Cells with a value of one are inside the specified range boundaries; cells with a value of zero are beyond range boundaries.
#'
#' @examples
#' #### Define an example UD
#' # We will use particles sampled by a particle filtering algorithm
#' # ... to create a UD:
#' particles <- pf_simplify(dat_dcpf_histories,
#'   summarise_pr = max,
#'   return = "archive"
#' )
#' # Define grids for UD estimation
#' map <- dat_dcpf_histories$args$bathy
#' habitat <- kud_habitat(map, plot = FALSE)
#' # Define UD as a raster
#' ud <- pf_kud_2(particles,
#'   bathy = map, grid = habitat,
#'   estimate_ud = kud_around_coastline,
#'   plot = FALSE
#' )
#'
#' #### Plot UD and home range estimators
#' pp <- par(mfrow = c(2, 2))
#' prettyGraphics::pretty_map(add_rasters = list(x = ud), main = "UD")
#' get_hr_full(ud, main = "Full range")
#' get_hr_home(ud, main = "Home range")
#' get_hr_core(ud, main = "Core range")
#' par(pp)
#'
#' #### Extract custom ranges with get_hr_prop()
#' get_hr_prop(ud, prop = 0.25)
#' get_hr_prop(ud, prop = 0.10)
#' get_hr_prop(ud, prop = 0.05)
#'
#' @references Jennrich, R. I. and Turner, F. B. (1969). Measurement of non-circular home range. Journal of Theoretical Biology, 22, 227--237.
#'
#' @author Edward Lavender
#' @name get_hr
NULL


#### get_hr_prop()
#' @name get_hr
#' @export

get_hr_prop <- function(x, prop = 0.5, plot = TRUE, add_raster = list(), add_contour = list(), ...) {
  if (!requireNamespace("spatialEco", quietly = TRUE)) {
    stop("This function requires the 'spatialEco' package.", call. = FALSE)
  }
  check_class(input = x, to_class = c("SpatRaster", "RasterLayer"), type = "stop")
  if (length(prop) != 1L) {
    stop("'prop' should be a single number (proportion).", call. = FALSE)
  }
  if (inherits(x, "RasterLayer")) {
    if (!requireNamespace("terra", quietly = TRUE)) {
      stop("This function requires the 'terra' package.", call. = FALSE)
    }
    x <- terra::rast(x)
  }
  x <- spatialEco::raster.vol(x, p = prop, sample = FALSE)
  if (inherits(x, "SpatRaster")) x <- raster::raster(x)
  if (plot) {
    if (!is.null(add_raster)) add_raster$x <- x
    prettyGraphics::pretty_map(add_rasters = add_raster, ...)
    if (!is.null(add_contour)) {
      add_contour$x <- x
      if (is.null(add_contour$add)) add_contour$add <- TRUE
      if (is.null(add_contour$nlevels)) add_contour$nlevels <- 1
      if (is.null(add_contour$drawlabels)) add_contour$drawlabels <- FALSE
      do.call(raster::contour, add_contour)
    }
  }
  return(invisible(x))
}

#### get_hr_core()
#' @name get_hr
#' @export

get_hr_core <- function(x, plot = TRUE, add_raster = list(), add_contour = list(), ...) {
  get_hr_prop(x = x, prop = 0.5, plot = plot, add_raster = add_raster, add_contour = add_contour, ...)
}

#### get_hr_home()
#' @name get_hr
#' @export

get_hr_home <- function(x, plot = TRUE, add_raster = list(), add_contour = list(), ...) {
  get_hr_prop(x = x, prop = 0.95, plot = plot, add_raster = add_raster, add_contour = add_contour, ...)
}

#### get_hr_full()
#' @name get_hr
#' @export

get_hr_full <- function(x, plot = TRUE, add_raster = list(), add_contour = list(), ...) {
  get_hr_prop(x = x, prop = 1, plot = plot, add_raster = add_raster, add_contour = add_contour, ...)
}
