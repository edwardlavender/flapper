######################################
######################################
#### update_extent()

#' @title Shrink or expand an \code{\link[raster]{extent}} object
#' @description This function updates a \code{\link[raster]{raster}}'s \code{\link[raster]{extent}} by shrinking or expanding the x and y limits.
#'
#' @param x A \code{\link[raster]{raster}} or an \code{\link[raster]{extent}} object.
#' @param x_shift A number that defines the change in the x limits. Positive numbers increase the extent (i.e., the lower limit is reduced and the upper limit is increased) while negative numbers reduce extent (i.e., the lower limit is increased and the upper limit is decreased).
#' @param y_shift A number that defines the change in the y limits. By default, this equals \code{x_shift} (i.e., the \code{\link[raster]{extent}} is changed uniformly in all directions).
#'
#' @return The function returns an updated \code{\link[raster]{extent}} object.
#'
#' @examples
#' # Define example raster
#' r <- raster::raster(nrows = 3, ncols = 3,
#'                     resolution = c(5, 5),
#'                     ext = raster::extent(0, 15, 0, 15))
#' # Increase raster extent by 5 units in every direction
#' update_extent(r, 5)
#' # Decrease raster extent by 5 units in every direction
#' update_extent(r, -5)
#' # Increase x and y extent parameters differently
#' update_extent(r, 5, 10)
#'
#' @author Edward Lavender
#' @export
#'

update_extent <- function(x, x_shift, y_shift = x_shift){
  if(inherits(x, "RasterLayer")) x <- raster::extent(x)
  x[1] <- x[1] - x_shift
  x[2] <- x[2] + x_shift
  x[3] <- x[3] - y_shift
  x[4] <- x[4] + y_shift
  return(x)
}


######################################
######################################
#### buffer_and_crop()

#' @title Buffer and crop a spatial object
#' @description This function creates a buffer around a spatial object and then crops another spatial object to lie within the extent of the buffered object.
#' @param to_buffer A spatial object to be buffered (see \code{\link[rgeos]{gBuffer}}).
#' @param to_crop A spatial object to be cropped by the buffered object (see \code{\link[raster]{crop}}).
#' @param buffer A named list of arguments, passed to \code{\link[rgeos]{gBuffer}} to buffer the \code{to_buffer} object (e.g., \code{buffer = list(width = 10)}).
#' @param ... Additional arguments passed to \code{\link[raster]{crop}}.
#' @details  This is a simple wrapper for \code{\link[rgeos]{gBuffer}} and \code{\link[raster]{crop}}. If \code{buffer = NULL}, the function simply implements \code{\link[raster]{crop}}.
#' @return The function returns the \code{to_crop} object, cropped to the \code{\link[raster]{extent}} of the buffered \code{to_buffer} object.
#' @examples
#' # Define an example raster
#' nrw <- ncl <- 50
#' r <- raster::raster(nrow = nrw, ncol = ncl)
#' r[] <- stats::runif(nrw*ncl, 0, 1)
#' # Buffer and crop the raster around an example location
#' xy <- sp::SpatialPoints(matrix(c(0, 0), ncol = 2))
#' r2  <- buffer_and_crop(to_buffer = xy,
#'                           to_crop = r,
#'                           buffer = list(width = 10))
#' # Visualise outputs
#' pp <- par(mfrow = c(1, 2))
#' raster::plot(r)
#' raster::plot(r2)
#' par(pp)
#'
#' @author Edward Lavender
#' @export
#'

buffer_and_crop <- function(to_buffer,
                            to_crop,
                            buffer = NULL,...){
  if(!is.null(buffer)){
    check_class(input = buffer, to_class = "list", type = "stop")
    if(!rlang::has_name(to_buffer, "spgeom")) buffer$spgeom <- to_buffer
    to_buffer <- do.call(rgeos::gBuffer, buffer)
  }
  out <- raster::crop(to_crop, raster::extent(to_buffer),...)
  return(out)
}


######################################
######################################
#### xy_from_clicks()

#' @title Get location coordinates from mouse click(s)
#' @description This function defines a two-column matrix of x, y coordinates from clicked locations on a map.
#' @return The function returns a two-column matrix with coordinates.
#' @examples
#' \dontrun{
#' raster::plot(dat_gebco)
#' xy <- xy_from_click()
#' graphics::points(xy, col = "red")
#' }
#' @author Edward Lavender
#' @export

xy_from_click <- function(){
  cat("Please click locations on the map and press [Esc] when you are done...\n")
  xy <- graphics::locator()
  xy <- matrix(c(xy$x, xy$y), ncol = 2, byrow = FALSE)
  return(xy)
}


######################################
######################################
#### crop_from_click()

#' @title Interactively crop a \code{\link[raster]{raster}}
#' @description Interactively crop a \code{\link[raster]{raster}} to a boundary box specified by mouse clicks.
#' @param x A \code{\link[raster]{raster}}.
#' @param plot A logical input that defines whether or not to plot \code{x} before and after cropping.
#' @param ... Additional arguments passed to \code{\link[prettyGraphics]{pretty_map}} to customise plots.
#'
#' @details The function is implemented as follows:
#' \enumerate{
#'   \item The function plots the supplied raster (\code{x}), if \code{plot = TRUE};
#'   \item Boundary coordinates delineating the area to which \code{x} should be cropped are defined by mouse clicks on the plot;
#'   \item The plotted raster is cropped using the minimum and maximum x and y coordinates of mouse clicks, via \code{\link[raster]{crop}};
#'   \item The cropped raster is plotted, if \code{plot = TRUE};
#'   \item The cropped raster is returned.
#' }
#' @return The function returns a \code{\link[raster]{raster}} and, if \code{plot = TRUE}, a plot of the area before/after cropping.
#' @examples
#' if(interactive()) crop_from_click(dat_gebco)
#' @author Edward Lavender
#' @export
#'

crop_from_click <- function(x, plot = TRUE,...){
  # Plot the raster
  if(plot){
    prettyGraphics::pretty_map(x,
                               add_rasters = list(x = x),...)
  }
  # Capture interactively defined locations
  cat("Please click four boundary locations on the map and press [Esc] when you are done...")
  dat <- graphics::locator()
  # Get extent of area to be cropped
  dat <- data.frame(x = dat$x, y = dat$y)
  xlim <- range(dat$x)
  ylim <- range(dat$y)
  ext <- raster::extent(xlim[1], xlim[2], ylim[1], ylim[2])
  # Crop raster
  x_crop <- raster::crop(x, ext)
  # Plot cropped raster and return
  if(plot){
    prettyGraphics::pretty_map(x_crop,
                               add_rasters = list(x = x_crop),...)
  }
  return(x_crop)
}


######################################
######################################
#### invert_poly()

#' @title Invert a (spatial) polygon
#' @description This function inverts a (spatial) polygon so that the `inside' of the original polygon becomes the `outside' and vice-versa. The function was motivated by marine applications in which polygons that define the coastline `contain' land and need to be inverted to define the sea.
#'
#' @param x An \code{\link[sp]{SpatialPolygons-class}} or \code{\link[sp]{SpatialPolygonsDataFrame-class}} object.
#' @param boundaries A \code{\link[raster]{extent}} object that defines the boundaries of the area under consideration. By default, this is defined by the extent of \code{x}.
#' @param ... Additional arguments passed to \code{\link[rgeos]{gDifference}}.
#'
#' @return The function returns a \code{\link[sp]{SpatialPolygons-class}} object.
#'
#' @examples
#' #### Example (1): Compare original and inverted polygon
#' # In this example, we have a polygon that defines the coastline
#' # ... with the polygon enclosing the land. We can invert
#' # ... the polygon to return a polygon that defines the sea.
#' pp <- graphics::par(mfrow = c(1, 2))
#' raster::plot(dat_coast, col = "darkgreen")
#' dat_sea <- invert_poly(dat_coast)
#' raster::plot(dat_sea, col = "skyblue")
#' graphics::par(pp)
#' # The CRS of the two objects is identical
#' raster::crs(dat_coast); raster::crs(dat_sea)
#' # Compare the classes of the two objects
#' class(dat_coast); class(dat_sea)
#' @author Edward Lavender
#' @export
#'

invert_poly <- function(x, boundaries = raster::extent(x),...){
  boundary_xy <- raster::coordinates(boundaries)
  boundary_poly <- sp::Polygon(boundary_xy)
  boundary_sp_poly <- sp::SpatialPolygons(list(sp::Polygons(list(boundary_poly), ID = 1)))
  raster::crs(boundary_sp_poly) <- raster::crs(x)
  x <- rgeos::gDifference(boundary_sp_poly, x,...)
  return(x)
}


######################################
######################################
#### cells_from_val()

#' @title Obtain a RasterLayer or the cells of RasterLayer that are equal to or lie within a range of specified values
#' @description This function obtains a RasterLayer or the cells of a RasterLayer that are equal to a specified value or lie within a specified range of values. To implement this function, a \code{\link[raster]{raster}} (\code{x}) and a value or range of values (\code{y}) for which a RasterLayer or the numbers of cells corresponding to those values are desired must be supplied. If a range of values is supplied, an additional argument (\code{interval}) controls whether or not cells within and equal to or simply within the specified range are returned.
#' @param x A \code{\link[raster]{raster}} object.
#' @param y A number or a vector with two numbers representing the lower and upper boundaries of the interval of values within which cells are identified.
#' @param interval If y is a vector of length two, \code{interval} is an integer that controls whether or not to query cells within and equal to (\code{interval = 1L}) or simply within (\code{interval = 2L}) the range specified by \code{y}.
#' @param cells A logical variable that defines whether or not to return a vector of cell numbers (\code{TRUE}) or a RasterLayer of the cells corresponding to \code{y}.
#' @param na.rm A logical variable that defines whether or not to ignore NAs.
#'
#' @return The function returns a RasterLayer (if \code{cells = FALSE}) or an integer vector of numbers (if \code{cells = TRUE}) that defines the cells that are equal to, or lie within, specified value(s) \code{y}.
#'
#' @examples
#' # Define an example RasterLayer
#' ncl <- 10
#' nrw <- 10
#' n   <- ncl*nrw
#' mat <- matrix(1:n, ncol = ncl, nrow = nrw, byrow = TRUE)
#' r <- raster::raster(mat)
#' # Visualise example RasterLayer
#' raster::plot(r)
#' raster::text(r)
#' # Obtain the number(s) of cells corresponding to a particular value
#' cells_from_val(r, 1)
#' # Obtain a RasterLayer of the cells corresponding to a particular value
#' r1 <- cells_from_val(r, 1, cells = FALSE)
#' raster::plot(r1)
#' # Obtain the number(s) of cells within or equal to a range of values
#' cells_from_val(r, c(1, 10))
#' # Obtain the number(s) of cells within a range of values
#' cells_from_val(r, c(1, 10), interval = 2L)
#' # As above but returning a raster layer
#' cells_from_val(r, c(1, 10), interval = 2L, cells = FALSE)
#' @author Edward Lavender
#' @export

cells_from_val <- function(x, y, interval = 1L, cells = TRUE, na.rm = TRUE){
  if(length(y) == 1){
    cells <- raster::Which(x == y, cells = cells, na.rm = na.rm)
  } else if(length(y) == 2){
    interval <- check_value(input = interval, supp = 1:2, warn = TRUE, default = 1)
    if(y[2] <= y[1]) stop("Nonsensical y range supplied: y[2] <= y[1].")
    if(interval == 1){
      cells <- raster::Which(x >= y[1] & x <= y[2], cells = cells, na.rm = na.rm)
    } else if(interval == 2){
      cells <- raster::Which(x > y[1] & x < y[2], cells = cells, na.rm = na.rm)
    }
  } else{
    stop("length(y) does not equal 1 or 2.")
  }
  return(cells)
}


######################################
######################################
#### mask_io()

#' @title Implement \code{\link[raster]{mask}} using the inside or outside of the mask
#' @description This function implements \code{\link[raster]{mask}} using the inside or the outside of the mask. The function is implemented in the same way as \code{\link[raster]{mask}} but with an additional logical argument, \code{mask_inside}. If \code{mask_inside = FALSE} (the default), values beyond the mask are masked and the function simply implements \code{\link[raster]{mask}}. In contrast, if \code{mask_inside = TRUE}, values within the mask are masked.
#'
#' @param x A Raster* (see \code{\link[raster]{mask}}).
#' @param mask A mask of class \code{\link[sp]{SpatialPolygonsDataFrame-class}}.
#' @param mask_inside A logical variable that defines whether or not to mask values in \code{x} outside (\code{FALSE}) or inside (\code{TRUE}) of the \code{mask}.
#' @param ... Additional arguments passed to \code{\link[raster]{mask}}.
#'
#' @details This function was motivated by animal movement datasets from coastal environments. For example, consider the cost of movement over a surface between two points for an exclusively benthic marine animal versus an exclusively terrestrial animal. With a \code{mask} that represents the coastline, for the exclusively marine animal, it is necessary to mask values inside the coastline (i.e., on land), so \code{mask_inside = TRUE}, to prevent movement on land; while for a terrestrial animal, it is necessary to mask values beyond the coastline (i.e., in the sea), so \code{mask_inside = FALSE}, to prevent movement into the sea.
#' @return The function returns a Raster* object in which values outside or inside of a supplied \code{\link[sp]{SpatialPolygonsDataFrame-class}} object have been masked.
#'
#' @examples
#' #### Define an example raster
#' # We will use some bathymetry data from around Oban
#' # ... which we will mask using some coastline data. All values on land
#' # ... are currently NA so we'll set these to 10 m for demonstration purposes.
#' dat_gebco[is.na(dat_gebco[])] <- 10
#' dat_gebco <- raster::crop(dat_gebco, raster::extent(dat_coast))
#' raster::plot(dat_gebco)
#'
#' #### Standard implementation with mask_inside = FALSE simply implements raster::mask()
#' m <- mask_io(dat_gebco, dat_coast)
#' raster::plot(dat_coast)
#' raster::plot(m, add = TRUE)
#'
#' #### Implementation with mask_inside = TRUE implements the mask within the coastline
#' m <- mask_io(dat_gebco, dat_coast, mask_inside = TRUE)
#' raster::plot(dat_coast)
#' raster::plot(m, add = TRUE)
#'
#' #### Additional arguments to raster::mask() can be passed via ... as usual
#' m <- mask_io(dat_gebco, dat_coast, mask_inside = TRUE, updatevalue = 1e5)
#' raster::plot(dat_coast)
#' raster::plot(m, add = TRUE)
#'
#' @author Edward Lavender
#' @export
#'

mask_io <- function(x, mask, mask_inside = FALSE,...){
  # Re-define mask if we are masking points inside the mask:
  if(mask_inside){
    area <- raster::extent(x)
    area <- sp::Polygon(raster::coordinates(area))
    area <- sp::SpatialPolygons(list(sp::Polygons(list(area), ID = 1)))
    raster::crs(area) <- raster::crs(x)
    mask <- rgeos::gDifference(area, mask)
  }
  # Mask raster
  x_masked <- raster::mask(x = x, mask = mask,...)
  return(x_masked)
}


######################################
######################################
#### split_raster_equally()

#' @title Split a raster into equal-area parts
#' @description This function splits a raster object into parts with approximately equal area.
#' @param r A \code{\link[raster]{raster}}.
#' @param n An integer that defines the number of parts into which to split the raster.
#' @details The raster (\code{r}) should not only contain NAs.
#' @note This function requires the `plyr' package.
#' @return The function returns a list containing the split raster components.
#'
#' @examples
#' l <- split_raster_equally(dat_gebco, 2)
#' l <- split_raster_equally(dat_gebco, 3)
#' pp <- graphics::par(mfrow = c(1, 3))
#' lapply(l, function(r) prettyGraphics::pretty_map(add_rasters = list(x = r)))
#' graphics::par(pp)
#'
#' @source The function taken and slightly modified from the `greenbrown' package (see https://rdrr.io/rforge/greenbrown/src/R/SplitRasterEqually.R). The function is defined separately in \code{\link[flapper]{flapper}} to reduce reliance on non-default packages.
#' @references Forkel M, Wutzler T (2015) greenbrown -- land surface phenology and trend analysis. A package for the R software. Version 2.2, 2015-04-15, http://greenbrown.r-forge.r-project.org/.
#' @export

split_raster_equally <- function(r, n) {
  # Check for plyr
  if(!requireNamespace("plyr", quietly = TRUE)) stop("This function requires the 'plyr' package.")
  # get total number of non NA grid cells
  mean.r <- raster::calc(r, function(x) {
    s <- sum(!is.na(x))
    s[s == 0] <- NA
    return(s)
  })
  ntotal <- sum(stats::na.omit(raster::values(mean.r)))
  if (n > ncol(mean.r)/2) {
    n <- ncol(mean.r)/2
    warning(paste("SplitRasterEqually: n changed to", n))
  }

  # compute optimal splitting based on weighted quantile
  xy <- raster::coordinates(mean.r)
  df <- data.frame(lon = xy[,1], val = raster::extract(mean.r, xy))
  df <- stats::na.omit(df)
  x.best <- stats::quantile(df$lon, seq(0, 1, length = n + 1), na.rm = TRUE)

  # split raster according to optimal splitting
  tiles.l <- plyr::llply(as.list(1:(length(x.best)-1)), function(i) {
    xmin <- x.best[i]
    xmax <- x.best[i +1 ]
    tile.r <- raster::crop(r, raster::extent(xmin,
                                             xmax,
                                             raster::extent(mean.r)@ymin,
                                             raster::extent(mean.r)@ymax))
    return(tile.r)
  })

  # Return output
  return(tiles.l)
}


######################################
######################################
#### sim_surface()

#' @title Populate a raster with simulated values
#' @description This function is designed to populate a raster with simulated values. To implement the function, a (blank) raster should be supplied. A user-defined function, or list of functions, is evaluated across this raster, or across sub-regions of this raster, to generate a new raster with simulated values.
#'
#' @param blank A \code{\link[raster]{raster}}.
#' @param n An integer that defines the number of (approximately equal area) pieces into which to split \code{blank}.
#' @param sim_values A function or, if \code{n > 1L}, a list of functions, that, for a given number of cells, simulate new values for those cells.
#' @param mask,mask_inside Arguments required to implement a spatial mask via \code{\link[flapper]{mask_io}}.
#' @param plot An integer that defines whether or not to plot a histogram of simulated values (\code{1L}), a heat map of the simulated raster (\code{2L}) or both (\code{1:2L}).
#'
#' @return The function returns a \code{\link[raster]{raster}}, with the same properties as \code{blank}, with values generated from the \code{sim_values} function(s).
#'
#' @examples
#' #### Example (1): Simulate values across the whole raster
#' sim_surface(dat_gebco,
#'             sim_values = function(n) stats::runif(n = n, 0, 1))
#' sim_surface(dat_gebco,
#'             sim_values = function(n) stats::rnorm(n = n, 0, 1))
#'
#' #### Example (2): Simulate values differently across different areas
#' # .. by defining the number of areas into which to split the raster
#' # .. and a list of function(s)
#' sim_surface(dat_gebco,
#'             n = 2, sim_values = list(function(n) stats::runif(n = n, 0, 1),
#'                                      function(n) stats::runif(n = n, 10, 11))
#'             )
#'
#' #### Example (3): Include a spatial mask
#' sim_surface(dat_gebco,
#'             n = 2, sim_values = list(function(n) stats::runif(n = n, 9, 10),
#'                                     function(n) stats::runif(n = n, 10, 11)),
#'             mask = dat_coast, mask_inside = TRUE
#'             )
#'
#' @author Edward Lavender
#' @export

sim_surface <- function(blank,
                        n = 1L,
                        sim_values,
                        mask = NULL, mask_inside = FALSE,
                        plot = 1:2L){

  #### Simulate values across a single raster
  if(n == 1L) {
    ncells <- raster::ncell(blank)
    surface <- blank
    raster::values(surface) <- sim_values(ncells)

    #### Simulate values across pieces of a raster
  } else {

    # Split the raster equally
    if(!inherits(sim_values, "list")) stop("If n > 1L, 'sim_values' should be a list of function(s).")
    blank_ls <- split_raster_equally(r = blank, n = n)

    # Loop over the list of blank rasters and list of functions
    # ... to apply to each chunk
    surface_ls <- mapply(blank_ls, sim_values, FUN = function(r, sv){
      ncells <- raster::ncell(r)
      raster::values(r) <- sv(ncells)
      return(r)
    }, SIMPLIFY = FALSE)

    # Rejoin rasters
    surface <- do.call(raster::merge, surface_ls)
  }

  #### Mask surface
  if(!is.null(mask)) surface <- mask_io(x = surface, mask = mask, mask_inside = mask_inside, updatevalue = NA)

  #### Plot surface
  if(1L %in% plot & 2L %in% plot) pp <- graphics::par(mfrow = c(1, 2)) else pp <- graphics::par()
  if(1L %in% plot) {
    raster::hist(surface)
  }
  if(2L %in% plot) {
    prettyGraphics::pretty_map(add_rasters = list(x = surface), verbose = FALSE)
  }
  graphics::par(pp)

  #### Return surface
  return(surface)
}


######################################
######################################
#### segments_cross_barrier()

#' @title Determine if Euclidean path segments cross a barrier
#' @description Given a sequence of `starting' and `ending' locations (\code{start} and \code{end}), this function determines whether or not the Euclidean path (`segment') between each location pair crosses a \code{barrier}.
#'
#' @param start A two-column matrix of coordinates that defines the `start' location of each segment.
#' @param end A two-column matrix of coordinates that defines the `end' location of each segment.
#' @param barrier A simple feature geometry that defines the barrier (see \code{\link[sf]{st_intersects}}).
#' @param distance (optional) A \code{\link[raster]{raster}} that defines distances from the \code{barrier}. If supplied, \code{mobility} is required (see below).
#' @param mobility (optional) If \code{distance} is supplied, \code{mobility} is a number that defines the distance threshold. Location pairs for which both locations are further than \code{mobility} from the \code{barrier} are assumed not to overlap with the \code{barrier}. If \code{start} and \code{end} are possible locations for an animal at a given pair of time steps, \code{mobility} is the distance that the individual could move between time steps (see \code{\link[flapper]{pf}}).
#'
#' @details
#'
#' This function was motivated by the need to support internal routines in \code{\link[flapper]{pf_simplify}}. Specifically, the function is used to minimise the number of shortest-distance calculations that are required by restricting calculations (if applicable) to particle pairs that require movement around a barrier, such as the coastline. (In these cases, Euclidean distances may be poor approximations of shortest distances that an aquatic animal must travel.)
#'
#' The function implements a three-step approach to derive barrier overlaps:
#' \enumerate{
#'   \item The function determines whether or not the minimum convex polygon (i.e., boundary box) around \code{start}/\code{end} intersects with \code{barrier}. If it does not, then no location segments can overlap with the barrier. This step is fast.
#'   \item If the locations' minimum convex polygon intersects with the \code{barrier}, and if \code{distance} and \code{mobility} have been supplied, the function extracts the distance of each location in \code{start} and \code{end} from the \code{barrier}. Location pairs for which both locations are further than \code{mobility} from the \code{barrier} are dropped. This step is also fast.
#'   \item For any remaining location pairs, the function links each \code{start} and \code{end} location and determines whether or not each linkage (`segment') intersects with the \code{barrier} using \code{\link[sf]{sf}} routines. This step can be slow for large numbers of locations (hence the first two filtering steps).
#' }
#'
#' The following criteria apply to applications of this function:
#' \enumerate{
#'   \item The number of observations in \code{start} and \code{end} must match.
#'   \item The coordinate reference system for \code{start}, \code{end} and \code{barrier} must match.
#'   \item If \code{distance} is supplied, \code{mobility} must also be supplied.
#'   \item The function requires the \code{\link[sfheaders]{sf_linestring}} and \code{\link[sf]{st_intersects}} functions.
#' }
#'
#' For speed in iterative applications, the function does not check whether or not these criteria are met.
#'
#' @return The function returns a one-column matrix, with each row corresponding to a row in \code{start}/\code{end}, with a logical value that defines whether or not the Euclidean path segment connecting those two locations crosses the \code{barrier} (\code{TRUE}) or not (\code{FALSE}).
#'
#' @examples
#' #### Plot example area and define barrier
#' raster::plot(dat_gebco)
#' raster::lines(dat_coast)
#' barrier <- sf::st_as_sf(dat_coast)
#'
#' #### Example (1): Implement function using barrier only
#'
#' ## Define example starting and ending locations
#' start <- matrix(c(701854.9, 6260399,
#'                   709202.5, 6258892), ncol = 2, byrow = TRUE)
#' end <- matrix(c(706753.3, 6264261,
#'                 709673.5, 6257102), ncol = 2, byrow = TRUE)
#'
#' ## Visualise segments
#' # ... The first segment crosses the coastline (our barrier)
#' # ... The second segment does not cross the coastline (our barrier)
#' graphics::arrows(x0 = start[1, 1], y0 = start[1, 2],
#'                  x1 = end[1, 1], y1 = end[1, 2])
#' graphics::arrows(x0 = start[2, 1], y0 = start[2, 2],
#'                  x1 = end[2, 1], y1 = end[2, 2])
#'
#' ## Implement function
#' segments_cross_barrier(start, end, barrier = barrier)
#'
#' #### Example (2): Implement function using barrier with distance and mobility
#'
#' ## Define distances from barrier
#' dat_dist <- raster::rasterize(dat_coast, dat_gebco)
#' dat_dist <- raster::distance(dat_dist)
#'
#' ## Implement function for a specified mobility parameter
#' segments_cross_barrier(start, end, barrier = barrier,
#'                        distance = dat_dist, mobility = 500)
#'
#' #### Example (3): With many locations, supplying distance improves speed
#'
#' ## Sample a large number of random starting/ending locations
#' start <- raster::sampleRandom(dat_gebco, size = 3000, xy = TRUE)[, 1:2]
#' end   <- raster::sampleRandom(dat_gebco, size = 3000, xy = TRUE)[, 1:2]
#'
#' ## Compare the duration of the function without/with distance included
#' # The first method without distance is much slower than the second method
#' # (~0.714 s versus 0.131 s for 3000 locations)
#' system.time(
#'   int_1 <- segments_cross_barrier(start, end, barrier = barrier)
#'   )
#' system.time(
#'   int_2 <- segments_cross_barrier(start, end, barrier = barrier,
#'                                   distance = dat_dist, mobility = 500)
#'   )
#'
#' ## The two approaches return identical solutions:
#' identical(int_1, int_2)
#'
#' @return Edward Lavender
#' @export

segments_cross_barrier <- function(start, end, barrier, distance = NULL, mobility = NULL){

  ## Define point boundaries
  xlim <- range(c(start[, 1], end[, 1]))
  ylim <- range(c(start[, 2], end[, 2]))
  boundaries <- matrix(c(xlim[1], ylim[1],
                         xlim[2], ylim[1],
                         xlim[2], ylim[2],
                         xlim[1], ylim[2]), ncol = 2, byrow = TRUE)
  boundaries <- rbind(boundaries, boundaries[1, , drop = FALSE])
  boundaries <- sf::st_polygon(list(boundaries))
  # plot(boundaries)

  ## Option (1): If the points' boundaries contain barrier(s),
  # ... we will work out for each line whether or not it crosses the barrier
  if(sf::st_intersects(boundaries, barrier, sparse = FALSE)){
    # Define point matrices as dataframes
    start <- data.frame(start[, 1:2])
    end   <- data.frame(end[, 1:2])
    colnames(start) <- colnames(end) <- c("x", "y")
    # If 'distance' has been supplied, we will focus on the subset of lines
    # ... for which at least one of the points is within 'mobility' of the barrier.
    # ... This should reduce the number of spatial intersections that are
    # ... required in many situations.
    if(!is.null(distance)){
      dat <- data.frame(index = seq_len(nrow(start)),
                        dist_1 = raster::extract(distance, start),
                        dist_2 = raster::extract(distance, end),
                        int    = FALSE)
      dat$bool <- (dat$dist_1 < mobility) | (dat$dist_2 < mobility)
      dat_sbt  <- dat %>% dplyr::filter(.data$bool)
      if(nrow(dat_sbt) > 0L){
        start <- start[dat_sbt$index, , drop = FALSE]
        end   <- end[dat_sbt$index, , drop = FALSE]
      } else start <- data.frame()
    }
    if(nrow(start) > 0L){
      # Assign line IDs
      start$linestring_id <- end$linestring_id <- seq_len(nrow(start))
      # Define lines
      lines <-
        dplyr::bind_rows(start, end) %>%
        dplyr::arrange(.data$linestring_id)
      lines <- sfheaders::sf_linestring(lines,
                                        x = "x", y = "y",
                                        linestring_id = "linestring_id")
      sf::st_crs(lines) <- sf::st_crs(barrier)
      int <- sf::st_intersects(lines, barrier, sparse = FALSE)
      if(!is.null(distance)){
        dat$int[dat_sbt$index] <- int
        return(matrix(dat$int, ncol = 1))
      } else return(int)

    } else return(matrix(FALSE, nrow = nrow(end), ncol = 1))

  ## Option (2): If the points' boundaries do not enclose any barriers
  # ... then none of the Euclidean paths can cross the barrier
  } else return(matrix(FALSE, nrow = nrow(end), ncol = 1))

}
