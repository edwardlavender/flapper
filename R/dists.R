########################################
########################################
#### raster:::..planedist2()

#' @title Calculate a matrix of Euclidean distances between points
#' @param p1 A coordinate matrix.
#' @param p2 A coordinate matrix.
#' @details This function assumes that coordinates are in planar space.
#' @source This function is an internal function in the \code{\link[raster]{raster}} package (https://rdrr.io/rforge/raster/src/R/pointdistance.R). It is defined separately in \code{\link[flapper]{flapper}} for stability.
#' @keywords internal

.planedist2 <- function(p1, p2) {
  z0 <- complex(, p1[, 1], p1[, 2])
  z1 <- complex(, p2[, 1], p2[, 2])
  outer(z0, z1, function(z0, z1) Mod(z0 - z1))
}


########################################
########################################
#### dist_btw_receivers()

#' @title Compute Euclidean distances between receivers
#' @description This function computes Euclidean distances (km) between all combinations of receivers.
#'
#' @param moorings A dataframe that defines each unique receiver deployment. This should contain the columns: `receiver_id', a unique identifier of each receiver; `receiver_long' or `receiver_easting', the longitude or easting of that receiver; and `receiver_lat' or `receiver_northing', the latitude or northing that receiver (see \code{\link[flapper]{dat_moorings}}).
#' @param f (optional) A function to process distances. For example, round distances to the nearest km with \code{f = function(x) round(x, digits = 0)} or simply \code{round}.
#' @param return A character that defines the class of the returned object (\code{"data.frame"} or \code{"matrix"}).
#'
#' @details Distances are calculated via \code{\link[raster]{pointDistance}}. If \code{moorings} contains `receiver_long' and `receiver_lat', \code{\link[raster]{pointDistance}} is implemented with \code{lonlat = TRUE} and distances are in km; otherwise \code{lonlat = FALSE} and distances are in map units over 1000 (i.e., km if map units are metres).
#'
#' To calculate distances between specific receiver pairs, call \code{\link[raster]{pointDistance}} directly.
#'
#' @return The function returns a dataframe, with columns `r1', `r2' and `dist', or a matrix. Dataframe columns define the IDs of each combination of receivers and the associated distance between them, in km (or map units/1000). Note that the dataframe contains duplicate combinations of receivers (e.g., both r1 = 1 and r2 = 2 and r1 = 2 and r2 = 1). Alternatively, matrix rows and columns define receiver IDs and cell values give distances between each combination.
#'
#' @examples
#' #### Example (1): Implementation using lat/long coordinates
#' dat <- data.frame(
#'   receiver_id = dat_moorings$receiver_id,
#'   receiver_long = dat_moorings$receiver_long,
#'   receiver_lat = dat_moorings$receiver_lat
#' )
#' # Compute distances
#' dist_btw_receivers_km <- dist_btw_receivers(dat)
#' head(dist_btw_receivers_km)
#'
#' #### Example (2): Implementation using planar coordinates
#' proj_wgs84 <- sp::CRS(SRS_string = "EPSG:4326")
#' proj_utm <- sp::CRS(SRS_string = "EPSG:32629")
#' xy <- sp::SpatialPoints(
#'   dat[, c("receiver_long", "receiver_lat")],
#'   proj_wgs84
#' )
#' xy <- sp::spTransform(xy, proj_utm)
#' xy <- sp::coordinates(xy)
#' dat <- data.frame(
#'   receiver_id = dat_moorings$receiver_id,
#'   receiver_easting = xy[, 1],
#'   receiver_northing = xy[, 2]
#' )
#' head(dist_btw_receivers(dat))
#'
#' #### Example (3): Post-process distances via the f argument
#' dist_btw_receivers_km_round <- dist_btw_receivers(dat, f = round)
#' head(dist_btw_receivers_km_round)
#' # convert distances to m
#' dist_btw_receivers_m <- dist_btw_receivers(dat, f = function(x) x * 1000)
#' head(dist_btw_receivers_m)
#'
#' #### Example (4) Return distances in a matrix
#' # Get distances
#' dist_btw_receivers(dat, return = "matrix")
#' # Compare sample to dataframe
#' dist_btw_receivers(dat)[1:5, ]
#' dist_btw_receivers(dat, return = "matrix")[1:5, 1:5]
#'
#' @author Edward Lavender
#' @export
#'

dist_btw_receivers <-
  function(moorings, f = NULL, return = c("data.frame", "matrix")) {
    #### Checks
    return <- match.arg(return)
    check_names(input = moorings, req = "receiver_id", extract_names = colnames, type = all)
    check_names(input = moorings, req = c("receiver_long", "receiver_easting"), extract_names = colnames, type = any)
    check_names(input = moorings, req = c("receiver_lat", "receiver_northing"), extract_names = colnames)
    if (any(c("receiver_long", "receiver_lat") %in% colnames(moorings))) {
      check_names(input = moorings, req = c("receiver_long", "receiver_lat"), extract_names = colnames, type = all)
    } else if (any(c("receiver_easting", "receiver_northing") %in% colnames(moorings))) {
      check_names(input = moorings, req = c("receiver_easting", "receiver_northing"), extract_names = colnames, type = all)
    }

    #### Define x and y coordinates
    if ("receiver_long" %in% colnames(moorings)) {
      lonlat <- TRUE
      moorings$x <- moorings$receiver_long
      moorings$y <- moorings$receiver_lat
    } else {
      lonlat <- FALSE
      moorings$x <- moorings$receiver_easting
      moorings$y <- moorings$receiver_northing
    }

    #### Compute the distances between each receiver combination in km (or map units/1000)
    dists <- moorings[, c("x", "y")]
    out <- raster::pointDistance(dists[, c("x", "y")], allpairs = TRUE, lonlat = lonlat)
    rownames(out) <- colnames(out) <- moorings$receiver_id
    out <- out / 1000
    if (!is.null(f)) out <- f(out)
    out[upper.tri(out)] <- t(out)[upper.tri(out)]
    if (return == "data.frame") {
      out <- data.frame(
        r1 = moorings$receiver_id[row(out)],
        r2 = moorings$receiver_id[col(out)],
        dist = as.vector(out)
      )
    }

    #### Return outputs
    return(out)
  }


######################################
######################################
#### dist_btw_clicks()

#' @title Calculate the distance between sequential mouse clicks on a map
#' @description This function calculates the distance between sequential mouse clicks on a plotted map, by combining \code{\link[graphics]{locator}} with a distance calculator, such as \code{\link[raster]{pointDistance}}.
#'
#' @param calc_distance A function that calculates distances between two sets of points, such as \code{\link[raster]{pointDistance}}. The first two arguments of this function must accept a dataframe comprising the x and y coordinates of the first and second set of points respectively. The function must return a numeric vector of distances between these.
#' @param ... Additional arguments passed to \code{calc_distance}, such as \code{lonlat} for \code{\link[raster]{pointDistance}}.
#' @param add_paths (optional) A named list of arguments, passed to \code{\link[prettyGraphics]{add_sp_path}}, to customise the paths added to the plot. \code{add_paths = NULL} suppresses this option.
#'
#' @return The function returns a dataframe with an integer ID for each path segment (`segment'), the first and second x and y coordinates (`x', `x2', `y', `y2`) and the distances between these points (`dist'). If \code{add_paths} is not \code{NULL}, the segments are drawn on the map.
#'
#' @examples
#' \dontrun{
#' raster::plot(dat_gebco)
#' dist_btw_clicks(lonlat = FALSE)
#' }
#'
#' @author Edward Lavender
#' @export
#'

dist_btw_clicks <- function(calc_distance = raster::pointDistance, ..., add_paths = list(length = 0.025)) {
  cat("Please click locations on the map and press [Esc] when you are done...\n")
  dat <- graphics::locator()
  cat("Getting distances...\n")
  dat <- data.frame(
    segment = 1:length(dat$x),
    x = dat$x,
    x2 = dplyr::lead(dat$x),
    y = dat$y,
    y2 = dplyr::lead(dat$y)
  )
  if (!is.null(add_paths)) {
    add_paths$x <- dat$x
    add_paths$y <- dat$y
    do.call(prettyGraphics::add_sp_path, add_paths)
  }
  dat <- dat[stats::complete.cases(dat), ]
  dat$dist <- raster::pointDistance(dat[, c("x", "y")], dat[, c("x2", "y2")], ...)
  return(dat)
}


######################################
######################################
#### dist_btw_points_3d()

#' @title Calculate the Euclidean distance(s) between points in three-dimensional space.
#' @description This function calculates the Euclidean distance(s) between points in three-dimensional space.
#' @param x1 A numeric vector that defines the x-coordinate(s) of the origin (for each point pair).
#' @param x2 A numeric vector that defines the x-coordinate(s) of the destination (for each point pair).
#' @param y1 A numeric vector that defines the y-coordinate(s) of the origin (for each point pair).
#' @param y2 A numeric vector that defines the y-coordinate(s) of the destination (for each point pair).
#' @param z1 A numeric vector that defines the z-coordinate(s) of the origin (for each point pair).
#' @param z2 A numeric vector that defines the z-coordinate(s) of the destination (for each point pair).
#' @details The distance between any two points in three dimensional space is given by Pythagoras' Theorem: \eqn{\sqrt{(x_2 - x_1)^2 + (y_2 - y_1)^2 + (z_2 - z_1)^2}}.
#' @details This function does not currently support non-planar coordinates.
#' @return The function returns a numeric vector that is equal to the Euclidean distance(s) between each pair of points in three-dimensional space.
#' @examples
#' #### Example (1): A single pair of points
#' dist_btw_points_3d(1, 2, 1, 2, 1, 2)
#' sqrt((2 - 1)^2 + (2 - 1)^2 + (2 - 1)^2)
#'
#' #### Example (2): Multiple pairs of points (e.g., an animal movement path)
#' xy <- data.frame(
#'   x = c(1, 2, 3, 4),
#'   y = c(1, 2, 3, 4),
#'   z = c(1, 2, 3, 4)
#' )
#' dist_btw_points_3d(
#'   xy$x, dplyr::lead(xy$x),
#'   xy$y, dplyr::lead(xy$y),
#'   xy$z, dplyr::lead(xy$z)
#' )
#' @author Edward Lavender
#' @export
#'

dist_btw_points_3d <- function(x1, x2, y1, y2, z1, z2) {
  sqrt((x2 - x1)^2 + (y2 - y1)^2 + (z2 - z1)^2)
}


######################################
######################################
#### dist_over_surface()

#' @title Calculate the total distance of a path over a three-dimensional surface
#' @description This function calculates the total distance of a path, whose horizontal coordinates are known, over a three-dimensional surface. To implement the function, the \code{path} should be supplied as a matrix or dataframe of coordinates or a \code{\link[sp]{SpatialLines}} object and the \code{surface} should be supplied as a \code{\link[raster]{raster}}. The function takes the horizontal coordinates of the \code{path} and extracts the values of the surface at these points, and then calculates the total distance of the path as the sum of the paired distances between each pair of points.
#' @param path A matrix or dataframe of horizontal coordinates (x, y) or a \code{\link[sp]{SpatialLines}} object which defines the path over a surface. The coordinate reference system (projection) of \code{path} should be the same as that for the \code{surface}.
#' @param surface A \code{\link[raster]{raster}} over which the movement that generated the path occurred. The coordinate reference system (projection) of \code{path} should be the same as that for the \code{surface} and the values of the \code{surface} should also be expressed in the same units (e.g., metres).
#' @details The total distance of a path over a three-dimensional surface is equal to the sum of the pairwise distances between each point (\eqn{i}) and its successor (\eqn{i + 1}) according to the equation: \deqn{\Sigma_{i = 1}^n \sqrt{(x_{i+1} - x_i)^2 + (y_{i + 1} - y_i)^2 + (z_{i + 1} - z_i)^2)}} where \eqn{x}, \eqn{y} and \eqn{z} are the x, y and z coordinates of each point in three-dimensional space and \eqn{n} is the total number of points minus 1. Pairwise distances are calculated via \code{\link[flapper]{dist_btw_points_3d}}. Note that for realistic distances, some interpolation (e.g., via least-cost paths) between points may be required to generate localisations at sufficiently high resolution to effectively capture the shape of the landscape.
#' @details This function does not currently support non-planar coordinates (\code{path} or \code{surface}).
#' @return The function returns a number equal to the total distance along the path.
#' @examples
#' #### Simulate a hypothetical landscape
#' # Define a miniature, blank landscape with known dimensions
#' proj_utm <- sp::CRS(SRS_string = "EPSG:32629")
#' r <- raster::raster(
#'   nrows = 3, ncols = 3,
#'   crs = proj_utm,
#'   resolution = c(5, 5),
#'   ext = raster::extent(0, 15, 0, 15)
#' )
#' # Define a matrix of hypothetical values for the landscape
#' mat <- matrix(c(
#'   5, 10, 3,
#'   2, 1, 4,
#'   5, 6, 6
#' ), ncol = 3, nrow = 3, byrow = TRUE)
#' r[] <- mat
#' # Visualise simulated landscape
#' raster::plot(r)
#' raster::text(r)
#'
#' #### Example (1): Total distance between two example adjacent points
#' path_cells <- c(1, 2)
#' path_matrix <- sp::coordinates(r)[path_cells, ]
#' dist_over_surface(path_matrix, r)
#' sqrt(5^2 + (10 - 5)^2)
#'
#' #### Example (2): Total distance between two example diagonal points
#' path_cells <- c(1, 5)
#' path_matrix <- sp::coordinates(r)[path_cells, ]
#' dist_over_surface(path_matrix, r)
#' sqrt(sqrt(5^2 + 5^2)^2 + (5 - 1)^2)
#'
#' #### Example (3): Total distance along a longer path
#' path_cells <- c(1, 2, 3)
#' path_matrix <- sp::coordinates(r)[path_cells, ]
#' dist_over_surface(path_matrix, r)
#' sqrt(5^2 + (10 - 5)^2) + sqrt(5^2 + (10 - 3)^2)
#'
#' #### Example (4): Total distance along an even longer path
#' path_cells <- c(1, 2, 3, 6, 8)
#' path_matrix <- sp::coordinates(r)[path_cells, ]
#' dist_over_surface(path_matrix, r)
#'
#' #### Example (5): A SpatialLines object can be used for the path
#' path_line <- Orcs::coords2Lines(path_matrix, ID = 1, proj4string = proj_utm)
#' raster::lines(path_line)
#' dist_over_surface(path_line, r)
#'
#' @seealso If the coordinates of the points are known in three dimensions already, \code{\link[flapper]{dist_btw_points_3d}} can be used directly.
#' @author Edward Lavender
#' @export
#'

dist_over_surface <- function(path, surface) {
  # extract the coordinates from the path, if necessary
  check_class(input = path, to_class = c("matrix", "data.frame", "SpatialLines"), type = "stop")
  if (!is.na(raster::crs(surface))) {
    message("The CRS of 'surface' is not NA. Note that this function only supports planar surfaces.")
  }
  if (inherits(path, "SpatialLines")) {
    if (!identical(raster::crs(surface), raster::crs(path))) {
      warning("The CRS of 'surface' and 'path' are not identical.",
        immediate. = TRUE, call. = FALSE
      )
    }
    xy <- raster::geom(path)[, c("x", "y")]
  } else {
    xy <- path
  }
  # calculate the number of coordinate pairs:
  n <- nrow(xy)
  # extract the value of the surface for these coordinates:
  z <- raster::extract(surface, xy)
  # create a dataframe
  xyz <- data.frame(x1 = xy[, 1], y1 = xy[, 2], z1 = z)
  # remove the final row (we can't calculate a distance from the last value)
  xyz <- xyz[1:(n - 1), ]
  # add x2, y2, and z2 columns that are one shifted (i.e. the next one in the sequence)
  xyz$x2 <- xy[2:n, 1]
  xyz$y2 <- xy[2:n, 2]
  xyz$z2 <- z[2:n]
  # calculate the distances between each pair of points in 3d space
  xyz$td_dist <- dist_btw_points_3d(xyz$x1, xyz$x2, xyz$y1, xyz$y2, xyz$z1, xyz$z2)
  # calculate the total distance implied by that spatial line:
  td_dist_tot <- sum(xyz$td_dist)
  # return the total distance:
  return(td_dist_tot)
}
