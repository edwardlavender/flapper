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
#### cells_from_val()

#' @title Obtain a RasterLayer or the cells of RasterLayer that are equal to or lie within a range of specified values
#' @description This function obtains a RasterLayer or the cells of a RasterLayer that are equal to a specified value or lie within a specified range of values. To implement this function, a \code{\link[raster]{raster}} (\code{x}) and a value or range of values (\code{y}) for which a RasterLayer or the numbers of cells corresponding to those values are desired must be supplied. If a range of values is supplied, an additional argument (\code{interval}) controls whether or not cells within and equal to or simply within the specified range are returned.
#' @param x A \code{\link[raster]{raster}} object.
#' @param y A number or a vector with two numbers representing the lower and upper boundaries of the interval of values within which cells are identified.
#' @param interval If y is a vector of length two, \code{interval} is an integer that controls whether or not to query cells within and equal to (\code{interval = 1L} or simply within (\code{interval = 2L}) the range specified by \code{y}.
#' @param cells A logical variable that defines whether or not to return a vector of cell numbers (\code{TRUE}) or a RasterLayer of the cells corresponding to \code{y}.
#' @param na.rm A logical variable that defines whether or not to ignore NAs.
#' @param ... Additional arguments (none implemented).
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

cells_from_val <- function(x, y, interval = 1L, cells = TRUE, na.rm = TRUE,...){
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
#### pythagoras_3d()

#' @title Calculate the Euclidean distance between points in three-dimensional space.
#' @description This function returns the Euclidean distance between points in three-dimensional space.
#' @param x1 A number that defines the x-coordinate of the first point.
#' @param x2 A number that defines the x-coordinate of the second point.
#' @param y1 A number that defines the y-coordinate of the first point.
#' @param y2 A number that defines the y-coordinate of the second point.
#' @param z1 A number that defines the z-coordinate of the first point.
#' @param z2 A number that defines the z-coordinate of the second point.
#' @details The distance between two points in three dimensional space is given by Pythagoras' Theorem: \eqn{\sqrt{(x_2 - x_1)^2 + (y_2 - y_1)^2 + (z_2 - z_1)^2}}.
#' @return A number that is equal to the Euclidean distance between points in three-dimensional space.
#' @examples
#' pythagoras_3d(1, 2, 1, 2, 1, 2)
#' @author Edward Lavender
#' @export
#'

pythagoras_3d <- function(x1, x2, y1, y2, z1, z2){
  td_dist <- sqrt((x2 - x1)^2 + (y2 - y1)^2 + (z2 - z1)^2)
  return(td_dist)
  }


######################################
######################################
#### dist_over_surface()

#' @title Calculate the total distance of a path over a three-dimensional surface
#' @description This function calculates the total distance of a path, whose horizontal coordinates are known, over a three-dimensional surface. To implement this function, the \code{path} should be supplied as a matrix or data.frame of coordinates or a \code{\link[sp]{SpatialLines}} object and the \code{surface} should be supplied as a \code{\link[raster]{raster}}. The function takes the horizontal coordinates of the \code{path} and extracts the values of the surface at these points, and then calculates the total distance of the path as the sum of the paired distances between each pair of points.
#' @param path A matrix or data.frame of horizontal coordinates (x, y) or a \code{\link[sp]{SpatialLines}} object which defines the path over a surface. The coordinate reference system (projection) of \code{path} should be the same as that for the \code{surface}.
#' @param surface A \code{\link[raster]{raster}} over which the movement that generated the path occurred. The coordinate reference system (projection) of \code{path} should be the same as that for the \code{surface} and the values of the \code{surface} should also be expressed in the same units (e.g., metres).
#' @details The total distance of a path over a three-dimensional surface is equal to the sum of the pairwise distances between each point (\eqn{i}) and its successor (\eqn{i + 1}) according to the equation: \deqn{\Sigma_{i = 1}^n [\sqrt{(x_{i+1} - x_i)^2 + (y_{i + 1} - y_i)^2 + (z_{i + 1} - z_i)^2)}]} where \eqn{x}, \eqn{y} and \eqn{z} are the x, y and z coordinates of each point in three-dimensional space. Pairwise distances are calculated via \code{\link[flapper]{pythagoras_3d}}.
#' @return The function returns a number equal to the total distance along the path.
#' @examples
#' #### Simulate a hypothetical landscape
#' # Define a miniature, blank landscape with known dimensions
#' proj_utm <- sp::CRS("+proj=utm +zone=29 ellps=WGS84")
#' r <- raster::raster(nrows = 3, ncols = 3,
#'                     crs = proj_utm,
#'                     resolution = c(5, 5),
#'                     ext = raster::extent(0, 15, 0, 15))
#' # Define a matrix of hypothetical values for the landscape
#' mat <- matrix(c(5, 10, 3,
#'                 2, 1, 4,
#'                 5, 6, 6), ncol = 3, nrow = 3, byrow = TRUE)
#' r[] <- mat
#' # Visualise simulated landscape
#' raster::plot(r)
#' raster::text(r)
#'
#' #### Example (1) Total distance between two example adjacent points
#' path_cells <- c(1, 2)
#' path_matrix <- sp::coordinates(r)[path_cells, ]
#' dist_over_surface(path_matrix, r)
#' sqrt(5^2 + (10-5)^2)
#'
#' #### Example (2) Total distance between two example diagonal points
#' path_cells <- c(1, 5)
#' path_matrix <- sp::coordinates(r)[path_cells, ]
#' dist_over_surface(path_matrix, r)
#' sqrt(sqrt(5^2 + 5^2)^2 + (5 - 1)^2)
#'
#' #### Example (3) Total distance along a longer path
#' path_cells <- c(1, 2, 3)
#' path_matrix <- sp::coordinates(r)[path_cells, ]
#' dist_over_surface(path_matrix, r)
#' sqrt(5^2 + (10-5)^2) + sqrt(5^2 + (10-3)^2)
#'
#' #### Example (4) Total distance along an even longer path
#' path_cells <- c(1, 2, 3, 6, 8)
#' path_matrix <- sp::coordinates(r)[path_cells, ]
#' dist_over_surface(path_matrix, r)
#'
#' #### Example (5) A SpatialLines object can be used for the path
#' path_line <- Orcs::coords2Lines(path_matrix, ID = 1, proj4string = proj_utm)
#' raster::lines(path_line)
#' dist_over_surface(path_line, r)
#'
#' @seealso If the coordinates of the points are known in three dimensions already, \code{\link[flapper]{pythagoras_3d}} can be used directly.
#' @author Edward Lavender
#' @export
#'

dist_over_surface <- function(path, surface){
  # extract the coordinates from the path, if necessary
  check_class(input = path, to_class = c("matrix", "data.frame", "SpatialLines"), type = "stop")
  if(inherits(path, "SpatialLines")) xy <- raster::geom(path)[, c("x", "y")] else xy <- path
  # calculate the number of coordinate pairs:
  n  <- nrow(xy)
  # extract the value of the surface for these coordinates:
  z <- raster::extract(surface, xy)
  # create a dataframe
  xyz <- data.frame(x1 = xy[, 1], y1 = xy[, 2], z1 = z)
  # remove the final row (we can't calculate a distance from the last value)
  xyz <- xyz[1:(n-1), ]
  # add x2, y2, and z2 columns that are one shifted (i.e. the next one in the sequence)
  xyz$x2 <- xy[2:n, 1]
  xyz$y2 <- xy[2:n, 2]
  xyz$z2 <- z[2:n]
  # calculate the distances between each pair of points in 3d space
  xyz$td_dist <- pythagoras_3d(xyz$x1, xyz$x2, xyz$y1, xyz$y2, xyz$z1, xyz$z2)
  # calculate the total distance implied by that spatial line:
  td_dist_tot <- sum(xyz$td_dist)
  # return the total distance:
  return(td_dist_tot)
}



#### End of code.
######################################
######################################
