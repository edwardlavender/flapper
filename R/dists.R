########################################
########################################
#### dist_btw_receivers()

#' @title Compute Euclidean distances between receivers
#' @description This function computes linear distances (km) between all combinations of receivers.
#'
#' @param moorings A dataframe which defines the location of each unique receiver combination. This should contain the columns: 'receiver_id', a unique identifier of each receiver, 'receiver_lat', the latitude of that receiver in decimal degrees; and 'receiver_long', the longitude of that receiver in decimal degrees (see \code{\link[flapper]{dat_moorings}}).
#' @param f (optional) A function which is used to process distances before these are returned. For example, it may be useful to round distances to nearest km with \code{f = function(x) round(x, digits = 0)}.
#'
#' @return The function returns a dataframe with columns 'r1', 'r2' and 'dist'. These define the IDs of each combination of receivers and the associated distance between them, in km. Note that the dataframe contains duplicate combinations of receivers (e.g., both r1 = 1 and r2 = 2 and r1 = 2 and r2 = 1).
#'
#' @examples
#' #### Example (1): Compute distances between all combinations of receivers in km
#' # Define dataframe with required columns
#' dat <- data.frame(receiver_id = dat_moorings$receiver_id,
#'                   receiver_lat = dat_moorings$receiver_lat,
#'                   receiver_long = dat_moorings$receiver_long)
#' # Compute distances
#' dist_btw_receivers_km <- dist_btw_receivers(dat)
#' head(dist_btw_receivers_km)
#'
#' #### Example (2): Post-process distances via the f argument
# round distances
#' dist_btw_receivers_km_round <- dist_btw_receivers(dat, f = round)
#' head(dist_btw_receivers_km_round)
#' # convert distances to m
#' dist_btw_receivers_m <- dist_btw_receivers(dat, f = function(x) x*1000)
#' head(dist_btw_receivers_m)
#'
#' @author Edward Lavender
#' @export
#'

dist_btw_receivers <-
  function(moorings,
           f = NULL){

    #### Checks
    stopifnot(all(c("receiver_id", "receiver_long", "receiver_lat") %in% colnames(moorings)))

    #### Define all combinations of receivers
    dists <- expand.grid(r1 = moorings$receiver_id, r2 = moorings$receiver_id)

    #### Define lat and long
    dists$lat1 <- moorings$receiver_lat[match(dists$r1, moorings$receiver_id)]
    dists$long1 <- moorings$receiver_long[match(dists$r1, moorings$receiver_id)]
    dists$lat2 <- moorings$receiver_lat[match(dists$r2, moorings$receiver_id)]
    dists$long2 <- moorings$receiver_long[match(dists$r2, moorings$receiver_id)]

    #### Compute the distances between each receiver combination in km
    dists$dist <- NA
    for(i in 1:nrow(dists)){
      # calculate distances
      dists$dist[i] <- geosphere::distGeo(c(dists$long1[i], dists$lat1[i]), c(dists$long2[i], dists$lat2[i]))
      # convert output from m to km
      dists$dist[i] <- dists$dist[i]/1000
    }

    #### Process dataframe
    if(!is.null(f)){
      dists$dist <- f(dists$dist)
    }
    dists <- dists[, c("r1", "r2", "dist")]

    #### Return dataframe
    return(dists)

}

######################################
######################################
#### dist_btw_points_3d()

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
#' dist_btw_points_3d(1, 2, 1, 2, 1, 2)
#' @author Edward Lavender
#' @export
#'

dist_btw_points_3d <- function(x1, x2, y1, y2, z1, z2){
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
#' @details The total distance of a path over a three-dimensional surface is equal to the sum of the pairwise distances between each point (\eqn{i}) and its successor (\eqn{i + 1}) according to the equation: \deqn{\Sigma_{i = 1}^n [\sqrt{(x_{i+1} - x_i)^2 + (y_{i + 1} - y_i)^2 + (z_{i + 1} - z_i)^2)}]} where \eqn{x}, \eqn{y} and \eqn{z} are the x, y and z coordinates of each point in three-dimensional space. Pairwise distances are calculated via \code{\link[flapper]{dist_btw_points_3d}}.
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
#' @seealso If the coordinates of the points are known in three dimensions already, \code{\link[flapper]{dist_btw_points_3d}} can be used directly.
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
  xyz$td_dist <- dist_btw_points_3d(xyz$x1, xyz$x2, xyz$y1, xyz$y2, xyz$z1, xyz$z2)
  # calculate the total distance implied by that spatial line:
  td_dist_tot <- sum(xyz$td_dist)
  # return the total distance:
  return(td_dist_tot)
}


######################################
######################################
#### lcp_over_surface()

#' @title Calculate shortest pathway(s) and/or distance(s) over a surface between origin and destination coordinates
#' @description This function computes the shortest pathway(s) and/or distance(s) over a \code{surface} between \code{origin} and \code{destination} coordinates. To implement this function, \code{origin} and \code{destination} coordinates need to be specified as matrices and the surface over which movement occurs should be supplied as a \code{\link[raster]{raster}}. Since determining shortest pathways can be computationally and memory-intensive, the \code{surface} can be reduced in size and/or resolution before these are computed, by (a) cropping the surface within user-defined extents; (b) focusing on a buffer zone along a Euclidean transect connecting \code{origin} and \code{destination} coordinates; (c) aggregating the surface to reduce the resolution; and/or (d) masking out areas over which movement is impossible (e.g., land for marine animals). Then, the function computes distances between connected cells, given (a) the planar distances between connected cells and (b) their difference in elevation. These distances are taken as a measure of 'cost'. For each pair of \code{origin} and \code{destination} coordinates, or for all combinations of coordinates, these distances are used to compute the least-cost pathway (i.e., the shortest pathway) and/or the distance of this pathway, using functions in the \code{\link[cppRouting]{cppRouting}} or \code{\link[gdistance]{gdistance}} package. The function returns the shortest pathway(s) and/or their distance(s) (m) along with a plot and a list of objects involved in the calculations.
#' @param origin A matrix which defines the coordinates (x, y) of the starting location(s). Coordinates should lie on a plane (i.e., Universal Transverse Mercator projection).
#' @param destination A matrix which defines the coordinates (x, y) of the finishing location(s). Coordinates should lie on a plane (i.e., Universal Transverse Mercator projection).
#' @param surface A \code{\link[raster]{raster}} over which the object (e.g., individual) must move from \code{origin} to \code{destination}. The \code{surface} must be planar (i.e., Universal Transverse Mercator projection) with units of metres in x, y and z directions (m). The \code{surface}'s \code{\link[raster]{resolution}} is taken to define the distance between horizontally and vertically connected cells and must be the same in both x and y directions (for \code{surface}'s with unequal horizontal resolution, \code{\link[raster]{resample}} can be used to equalise resolution: see Examples). Any cells with NA values (e.g., due to missing data) are treated as 'impossible' to move though by the algorithm. In this case, the \code{surface} might need to be pre-processed so that NAs are replaced/removed before implementing the function, depending on their source.
#' @param crop (optional) An \code{\link[raster]{extent}} object that is used to \code{\link[raster]{crop}} the extent of the \code{surface}, before the least-cost algorithms are implemented. This may be useful for large rasters to reduce memory requirements/computation time.
#' @param buffer (optional) A named list of arguments, passed to \code{\link[rgeos]{gBuffer}} (e.g. \code{buffer = list(width = 1000)}) (m) that is used to define a buffer around a Euclidean transect connecting the \code{origin} and \code{destination}. (This option can only be implemented for a single \code{origin} and \code{destination} pair.) The \code{surface} is then cropped to the extent of this buffer before the least-cost algorithms are implemented. This may be useful for large rasters to reduce memory requirements and/or computation time.
#' @param aggregate (optional) A named list of arguments, passed to \code{\link[raster]{aggregate}}, to aggregate raster cells before the least-cost algorithms are implemented. This may be useful for large rasters to reduce memory requirements and/or computation time.
#' @param mask (optional) A Raster or Spatial \code{\link[raster]{mask}} that is used to prevent movement over 'impossible' areas on the \code{surface}. This must also lie on a planar surface (i.e., Universal Transverse Mercator projection). For example, for marine animals, \code{mask} might be a \code{\link[sp]{SpatialPolygonsDataFrame}} which defines the coastline. The effect of the \code{mask} depends on \code{mask_inside} (see below).
#' @param mask_inside A logical input that defines whether or not to mask the \code{surface} inside (\code{TRUE}) or outside (\code{FALSE}) of the \code{mask} (see \code{\link[flapper]{mask_io}}).
#' @param plot A logical input that defines whether or not to plot the inputted and processed surfaces. If \code{TRUE}, the inputted and processed plots are produced side-by-side. For the inputted surface, the \code{mask} and the region selected (via \code{crop} and/or \code{buffer}) are shown along with the \code{origin} and \code{destination}. For the processed surface, the surface and the \code{origin} and \code{destination} are shown, along with the shortest pathway(s) (if and once computed: see \code{goal}). This is useful for checking that any \code{surface} processing steps have been applied correctly and the \code{origin} and \code{destination} are positioned correctly on the \code{surface}.
#' @param goal An integer that defines the output of the function: \code{goal = 1} computes shortest distances, \code{goal = 2} computes shortest pathways and \code{goal = 3} computes both shortest pathways and the corresponding distances. Note that \code{goal = 3} results in least-cost algorithms being implemented twice, which will be inefficient for large problems; in this case, use \code{goal = 2} to compute shortest pathways and then calculate their distance using outputs returned by the function (see Value).
#' @param combination A character string (\code{"pair"} or \code{"matrix"}) that defines whether or not to compute shortest distances/pathways for (a) each sequential \code{origin} and \code{destination} pair of coordinates (\code{combination = "pair"}) or (b) all combinations of \code{origin} and \code{destination} coordinates (\code{combination = "matrix"}). This argument is only applicable if there is more than one \code{origin} and \code{destination}. For \code{combination = "pair"}, the number of \code{origin} and \code{destination} coordinates needs to be the same, since each \code{origin} is matched with each \code{destination}.
#' @param method A character string (\code{"cppRouting"} or \code{"gdistance"}) that defines the method used to compute the shortest distances between the \code{origin} and the \code{destination}. \code{"cppRouting"} is the default \code{method}. Under this option, functions in the \code{\link[cppRouting]{cppRouting}} package are used to compute the shortest pathways (\code{\link[cppRouting]{get_path_pair}} or \code{\link[cppRouting]{get_multi_paths}} for each pair of coordinates or for all combinations of coordinates, respectively) and/or distances (\code{\link[cppRouting]{get_distance_pair}} or \code{\link[cppRouting]{get_distance_matrix}}). This package implements functions written in C++ massively outperforms the other \code{method = "gdistance"} for large problems. Otherwise, if \code{method = "gdistance"}, functions in the \code{\link[gdistance]{gdistance}} are called iteratively to compute shortest pathways (via \code{\link[gdistance]{shortestPath}}) or distances (via \code{\link[gdistance]{costDistance}}).
#' @param cppRouting_algorithm A character string that defines the algorithm used to compute shortest pathways or distances. This is only applicable if \code{method = "cppRouting"}: \code{method = "gdistance"} implements Dijkstra's algorithm only. For shortest pathways or their distances between pairs of coordinates, the options are \code{"Dijkstra"}, \code{"bi"}, \code{"A*"} or \code{"NBA"} for the uni-directional Dijkstra, bi-directional Dijkstra, A star unidirectional search or new bi-directional A star algorithms respectively (see \code{\link[cppRouting]{get_path_pair}} or \code{\link[cppRouting]{get_distance_pair}}). For shortest pathways between all combinations of coordinates, \code{cppRouting_algorithm} is ignored and the Dijkstra algorithm is implemented recursively. For shortest distances between all combinations of coordinates, the options are \code{"phast"} or \code{"mch"} (see \code{\link[cppRouting]{get_distance_matrix}}).
#' @param use_all_cores,cl,varlist Parallelisation arguments for \code{method = "cppRouting"} (\code{use_all_cores}) or \code{method = "gdistance"} (\code{cl} and \code{varlist}) respectively. If \code{method = "cppRouting"}, parallelisation is implemented via \code{use_all_cores} for computing shortest distances only (not computing shortest pathways). \code{use_all_cores} is a logical input that defines whether or not to use all cores for computing shortest distance(s).  If \code{method = "gdistance"}, parallelisation is implemented via \code{cl} and \code{varlist} for both shortest pathways and distances function calls. \code{cl} is a cluster object created by \code{\link[parallel]{makeCluster}}. If supplied, the connection to the cluster is stopped within the function. \code{varlist} is a character vector of containing the names of exported objects. This may be required if \code{cl} is supplied. This is passed to the \code{varlist} argument of \code{\link[parallel]{clusterExport}}. Exported objects must be located in the global environment.
#' @param check A logical input that defines whether or not to check function inputs. If \code{TRUE}, internal checks are implemented to check user-inputs and whether or not inputted coordinates are in appropriate places on the processed \code{surface} (for instance, to ensure inputted coordinates do not lie over masked areas). This helps to prevent intractable error messages. If \code{FALSE}, these checks are not implemented, so function progress may be faster initially (especially for large \code{origin}/\code{destination} coordinate matrices).
#' @param verbose A logical input that defines whether or not to print messages to the console to monitor function progress. This is especially useful with a large \code{surface} since the algorithms are computationally intensive.
#' @param ... Additional arguments (none implemented).
#'
#' @details
#' \subsection{Methods}{
#' This function was motivated by the need to determine the shortest pathways and their distances between points for benthic animals, which must move over the seabed to navigate from A to B. For these animals, especially in areas with heterogeneous bathymetric landscapes and/or coastline, the shortest pathway that an individual must travel to move from A and B may differ substantially from the Euclidean pathway that is often used as a proxy for distance in biological studies. However, this function can still be used in situations where the surface over which an individual must move is irrelevant (e.g., for a pelagic animal), by supplying a flat surface; then shortest pathways/distances simply depend on the planar distances between locations and any barriers (e.g., the coastline). (However, this process will be somewhat inefficient.)
#'
#' The function conceptualises a object moving across a landscape as a queen on a chessboard which can move, in eight directions around its current position, across this surface. Given the potentially large number of possible pathways between an \code{origin} and \code{destination}, the surface may be reduced in extent or size before the game begins. To determine shortest pathway/distance over the surface between each \code{origin} and \code{destination} pair/combination, the function first considers the distance that an object must travel between pairs of connected cells. This depends on the planar distances between cells and their differences in elevation. Planar distances (\eqn{d_p}, m) depend on the movement type: under a rook's movement (i.e., horizontally or vertically), the distance (\eqn{d_{p,r}}) between connected cells is extracted from the raster's resolution (which is assumed to be identical in the x and y directions); under a bishop's movement (i.e., diagonally), the distance between connected cells \eqn{d_{p,b}} is given by Pythagoras' Theorem: \eqn{d_{p,b} = \sqrt{(d_{p, r}^2 + d_{p, r}^2)}}. Vertical distances (\eqn{d_v}, m) are simply the differences in height between cells. The total distance (\eqn{d_t}) between any two connected cells is a combination of these distances given by Pythagoras' Theorem: \eqn{d_t = \sqrt{(d_p^2 + d_v^2)}}. These distances are taken to define the 'cost' of movement between connected cells. Thus, 'costs' are symmetric (i.e., the cost of moving from A to B equals the cost of moving from B to A).
#'
#' This cost surface is then used to compute the shortest pathway and/or distance of the shortest path between each \code{origin} and \code{destination} pair/combination using functions in the \code{\link[cppRouting]{cppRouting}} or \code{\link[gdistance]{gdistance}} package. The functions implemented depend on the \code{goal} (i.e., whether the aim is to compute shortest pathways, shortest distances or both) and, if there is more than one \code{origin}/\code{destination}, the \code{combination} type (i.e., whether to compute shortest pathways/distances for each sequential pair of coordinates or all possible combinations of coordinates).
#' }
#'
#' \subsection{Warnings}{
#' The function returns a warning produced by \code{\link[gdistance]{transition}} which is implemented to facilitate the definition of the cost surface, before shortest pathways/distances are computed by either method: 'In .TfromR(x, transitionFunction, directions, symm) : transition function gives negative values'. This warning arises because the height differences between connecting cells can be negative. It can be safely ignored.
#' }
#'
#' @return
#' \subsection{A named list}{The function returns a named list. The most important element(s) of this list are 'path_lcp' and/or 'dist_lcp', the shortest pathway(s) and/or distance(s) (m) between \code{origin} and \code{destination} coordinate pairs/combinations. 'path_lcp' is returned if \code{goal = 2} or \code{goal = 3} and 'dist_lcp' is returned if \code{goal = 1} or \code{goal = 3}. 'path_lcp' contains (a) a dataframe with the cells comprising each path ('cells'), (b) a named list containing a \code{\link[sp]{SpatialLines}} object for each path ('SpatialLines') and (c) a named list of matrices of the coordinates of each path ('coordinates'). 'dist_lcp' is a (a) numeric vector or (b) matrix with the distances (m) between each pair or combination of coordinates respectively. If 'dist_lcp' is computed, dist_euclid', the Euclidean distances (m) between the \code{origin} and \code{destination}, is also returned for comparison.
#' }
#'
#' \subsection{Common elements}{Other elements of the list record important outputs at sequential stages of the algorithm's progression. These include the following elements: 'args', a named list of user inputs; 'time', a dataframe that defines the times of sequential stages in the algorithm's progression; ; 'surface', the surface over which shortest distances are computed (this may differ from the inputted surface if any of the processing options, such as \code{crop}, have been implemented); 'surface_param', a named list that defines the cell IDs, the number of rows, the number of columns, the coordinates of the implemented surface and the cell IDs of the \code{origin} and \code{destination} nodes; 'cost', a named list of arguments that defines the distances (m) between connected cells under a rook's or bishop's movement ('dist_rook' and 'dist_bishop'), the planar and vertical distances between connected cells ('dist_planar' and 'dist_vertical') and the total distance between connected cells ('dist_total'); and 'cppRouting_param' or 'gdistance_param', a named list of arguments used to compute shortest pathways/distances via \code{\link[cppRouting]{cppRouting}} or \code{\link[gdistance]{gdistance}} (see below).
#' }
#'
#' \subsection{Method-specific elements}{If \code{method = "cppRouting"}, the 'cppRouting_param' list contains a named list of arguments passed to \code{\link[cppRouting]{makegraph}} ('makegraph_param') as well as \code{\link[cppRouting]{get_path_pair}} ('get_path_pair_param') or \code{\link[cppRouting]{get_multi_paths}} ('get_multi_paths_param') and/or \code{\link[cppRouting]{get_distance_pair}} ('get_distance_pair_param') or \code{\link[cppRouting]{get_distance_matrix}} ('get_distance_matrix_param'), depending on whether or not shortest pathways and/or distances have been computed (see \code{goal}) and whether or not shortest pathways/distances have been computed for each pair of coordinates or all combinations of coordinates. If \code{method = "gdistance"}, this list contains a named list of arguments passed iteratively, for each pair/combination of coordinates, to \code{\link[gdistance]{shortestPath}} ('shortestPath_param') or \code{\link[gdistance]{costDistance}} ('costDistance_param'). This includes an object of class TransitionLayer (see \code{\link[gdistance]{Transition-classes}}), in which the \code{transitionMatrix} slot contains a (sparse) matrix that defines the ease of moving between connected cells (the reciprocal of the 'dist_total' matrix).
#' }
#'
#' \subsection{Plot}{If \code{plot = TRUE}, a plot is also produced of the inputted and processed surfaces that are used in the calculations, along with the shortest pathway(s) (if and once computed).
#' }
#'
#' @examples
#' #### Example types
#' # Shortest distances between a single origin and a single destination
#' # Shortest pathways between a single origin and a single destination
#' # Shortest distances/pathways between origin/destination pairs
#' # Shortest distances/pathways between all origin/destination combinations
#'
#' #### Simulate a hypothetical landscape
#' # Define a miniature, blank landscape with appropriate dimensions
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
#' # Extract coordinates of cells
#' rxy <- raster::coordinates(r)
#'
#'
#' ############################################################################
#' #### Shortest distances between a single origin and a single destination
#'
#' #### Example (1): Find the distance between a single origin and destination
#' # ... using the "cppRouting" method:
#' out1 <- lcp_over_surface(origin = rxy[1, , drop = FALSE],
#'                          destination = rxy[2, , drop = FALSE],
#'                          surface = r)
#' # Extract shortest distance
#' out1$dist_lcp
#'
#' #### Example (1) continued: An explanation of function outputs
#' # The function returns a list:
#' # The 'args' element simply contains all inputted arguments
#' out1$args
#' # The 'surface' element contains the surface used to compute shortest distances
#' # ... This may differ from $args$surface if cropped, buffered etc.
#' out1$surface
#' # The 'surface_param' element contains the cell IDs, number of rows, cells and coordinates
#' # ... of this surface
#' out1$surface_param
#' # The 'cost' element is a list of objects that define the cost matrix:
#' # ... 'dist_rook' and 'dist_bishop' are matrices which define the distance of planar
#' # ... ... movement from one cell to any other cell under a rook's or bishop's movement.
#' # ... ... For example, the planar distance of moving from cell 1 to cell 2 to is 5 m:
#' out1$cost$dist_rook
#' out1$cost$dist_bishop
#' # ... 'dist_planar' gives the planar distance between connected cell combinations
#' # ... ... under a queen's movements:
#' out1$cost$dist_planar
#' # ... 'dist_vertical' gives the vertical distance between connected cells
#' out1$cost$dist_vertical
#' # ... and 'dist_total' gives the total distance between connected cells
#' out1$cost$dist_total
#' # 'dist_euclid' is the Euclidean distance between the origin and destination
#' out1$dist_euclid
#' # 'cppRouting_param' contains lists of parameters passed to (a) makegraph() and (b)
#' # ... get_distance_matrix() to compute shortest distances using cppRouting
#' utils::str(out1$cppRouting_param)
#' # 'time' records the time of each stage
#' out1$time
#'
#' #### Example (2): Find the distance between a single origin and destination
#' # ... using the "gdistance" method
#' out2 <- lcp_over_surface(origin = rxy[1, , drop = FALSE],
#'                          destination = rxy[2, , drop = FALSE],
#'                          surface = r,
#'                          method = "gdistance"
#'                          )
#' # Extract distance
#' out2$dist_lcp
#' # Elements of the returned list are the same apart from 'gdistance_param'
#' # ... which contains a list of arguments passed to costDistance() to compute
#' # ... shortest distances
#' utils::str(out2$gdistance_param)
#'
#' #### Example (3): Find the distances between other origins and destinations
#' ## Implement function to determine shortest distances:
#' # shortest distance between cell 1 and 6 via gdistance
#' lcp_over_surface(origin = rxy[1, , drop = FALSE],
#'                  destination = rxy[6, , drop = FALSE],
#'                  surface = r,
#'                  method = "gdistance",
#'                  verbose = FALSE)$dist_lcp
#' # shortest distance between cell 1 and 6 via cppRouting
#' lcp_over_surface(origin = rxy[1, , drop = FALSE],
#'                  destination = rxy[6, , drop = FALSE],
#'                  surface = r,
#'                  method = "cppRouting",
#'                  verbose = FALSE)$dist_lcp
#' ## Compare to manually computed distances
#' # The shortest distance from cell 1 to cell 6 is to move diagonally
#' # ... from cell 1 to 5 and then 6.
#' # Define planar distances of moving like a rook (d_pr) or bishop (d_pb)
#' d_pr <- 5
#' d_pb <- sqrt(5^2 + 5^2)
#' # Define total distance travelled along shortest path, as computed by the
#' # ... algorithm to demonstrate we obtain the same value:
#' sqrt((5 - 1)^2 + d_pb^2) + sqrt((4 - 1)^2 + d_pr^2)
#'
#' #### Example (4): Find the shortest distances around NAs
#' ## Force the 5th cell to be NA
#' rtmp <- r
#' rtmp[5] <- NA
#' raster::plot(rtmp); raster::text(rtmp)
#' ## Compute shortest distances via algorithm:
#' # Now compute shortest distances, which we can see have increased the distance
#' # ... since movement through an NA cell is not allowed. Therefore, if this NA
#' # ... reflects missing data, it may be appropriate to interpolate NAs using
#' # ... surrounding cells (e.g., see raster::approxNA()) so that movement
#' # ... is possible through these cells
#' out1 <- lcp_over_surface(origin = rxy[1, , drop = FALSE],
#'                          destination = rxy[6, , drop = FALSE],
#'                          surface = rtmp,
#'                          method = "gdistance",
#'                          verbose = FALSE)
#' out2 <- lcp_over_surface(origin = rxy[1, , drop = FALSE],
#'                          destination = rxy[6, , drop = FALSE],
#'                          surface = rtmp,
#'                          method = "cppRouting",
#'                          verbose = FALSE)
#' out1$dist_lcp; out2$dist_lcp
#' ## Compare to manual calculations of shortest route:
#' # Route option one: cell 1 to 2, 2 to 6
#' sqrt((5 - 10)^2 + d_pr^2) + sqrt((10 - 4)^2 + d_pb^2)
#' # Or, using the numbers computed in the dist_total object:
#' out1$cost$dist_total[1, 2] + out1$cost$dist_total[2, 6]
#' ## Compare to effect of making a value in the landscape extremely large
#' # In the same way, we can force the shortest pathway away from particular areas
#' # ... by making the height of the landscape in those areas very large or Inf:
#' rtmp[5] <- 1e20
#' raster::plot(rtmp); raster::text(rtmp)
#' lcp_over_surface(origin = rxy[1, , drop = FALSE],
#'                  destination = rxy[6, , drop = FALSE],
#'                  surface = rtmp,
#'                  method = "cppRouting",
#'                  verbose = FALSE)$dist_lcp
#' lcp_over_surface(origin = rxy[1, , drop = FALSE],
#'                  destination = rxy[6, , drop = FALSE],
#'                  surface = rtmp,
#'                  method = "gdistance",
#'                  verbose = FALSE)$dist_lcp
#'
#' #### Example (4): Find the distances between points on a real landscape
#' ## We will use some example bathymetry data:
#' dat_gebco_oban <- prettyGraphics::dat_gebco
#' raster::plot(dat_gebco_oban)
#' ## Process bathymetry data before function implementation
#' # (a) Define utm coordinates:
#' dat_gebco_utm <- raster::projectRaster(dat_gebco_oban, crs = proj_utm)
#' raster::res(dat_gebco_utm)
#' # (b) Resample so that the resolution in the x and y directions is identical
#' dat_gebco_utm_planar <- raster::raster(crs = proj_utm,
#'                                        ext = raster::extent(dat_gebco_utm),
#'                                        resolution = 250)
#' dat_gebco_utm_planar <- raster::resample(dat_gebco_utm, dat_gebco_utm_planar, method = "bilinear")
#' # Examine processed raster
#' pp <- par(mfrow = c(1, 2))
#' raster::plot(dat_gebco_utm, main = "UTM raster")
#' raster::plot(dat_gebco_utm_planar, main = "UTM raster with equal res")
#' par(pp)
#' ## Define example origin and destination
#' set.seed(1)
#' dat_gebco_utm_planar_xy <- raster::coordinates(dat_gebco_utm_planar)
#' index       <- sample(1:nrow(dat_gebco_utm_planar_xy), 2)
#' origin      <- dat_gebco_utm_planar_xy[index[1], , drop = FALSE]
#' destination <- dat_gebco_utm_planar_xy[index[2], , drop = FALSE]
#' ## Implement function to compute shortest distances
#' out_gebco1 <- lcp_over_surface(origin = origin,
#'                                destination = destination,
#'                                surface = dat_gebco_utm_planar,
#'                                method = "gdistance")
#' out_gebco2 <- lcp_over_surface(origin = origin,
#'                                destination = destination,
#'                                surface = dat_gebco_utm_planar,
#'                               method = "cppRouting")
#' # Compare Euclidean and shortest distances
#' out_gebco1$dist_euclid; out_gebco1$dist_lcp
#' out_gebco2$dist_euclid; out_gebco2$dist_lcp
#'
#' #### Example (5A): Reduce the complexity of the landscape by cropping
#' ext <- raster::extent(dat_gebco_utm_planar)
#' ext[3] <- 6255000
#' ext[4] <- 6265000
#' out_gebco1 <- lcp_over_surface(origin = origin,
#'                                destination = destination,
#'                                surface = dat_gebco_utm_planar,
#'                                method = "gdistance",
#'                                crop = ext)
#' out_gebco1$dist_lcp
#'
#' #### Example (5B): Reduce the complexity of the landscape around a buffer
#' out_gebco1 <- lcp_over_surface(origin = origin,
#'                                destination = destination,
#'                                surface = dat_gebco_utm_planar,
#'                                method = "gdistance",
#'                                buffer = list(width = 1000))
#' out_gebco1$dist_lcp
#'
#' #### Example (5C):  Reduce the complexity of the landscape via aggregation
#' out_gebco1 <- lcp_over_surface(origin = origin,
#'                                destination = destination,
#'                                surface = dat_gebco_utm_planar,
#'                                method = "gdistance",
#'                                buffer = list(width = 1000),
#'                                aggregate = list(fact = 5, fun = mean, na.rm = TRUE))
#' out_gebco1$dist_lcp
#'
#' #### Example (6C): Implement a mask
#' # Define coastline
#' dat_coast_around_oban <- prettyGraphics::dat_coast_around_oban
#' coastline <- sp::spTransform(dat_coast_around_oban, proj_utm)
#' # Visualise bathymetry and coastline
#' raster::plot(dat_gebco_utm_planar)
#' raster::lines(coastline)
#' # Define example origin and destination within the sea
#' origin      <- matrix(c(714000, 6260000), ncol = 2)
#' destination <- matrix(c(721000, 6265000), ncol = 2)
#' # Implement algorithm
#' out_gebco1 <- lcp_over_surface(origin = origin,
#'                                destination = destination,
#'                                surface = dat_gebco_utm_planar,
#'                                method = "gdistance",
#'                                crop = raster::extent(coastline),
#'                                mask = coastline,
#'                                mask_inside = TRUE)
#' out_gebco2 <- lcp_over_surface(origin = origin,
#'                                destination = destination,
#'                                surface = dat_gebco_utm_planar,
#'                                method = "cppRouting",
#'                                crop = raster::extent(coastline),
#'                                mask = coastline,
#'                                mask_inside = TRUE)
#' # Compare Euclidean and least-cost distances
#' out_gebco1$dist_euclid; out_gebco1$dist_lcp
#' out_gebco2$dist_euclid; out_gebco2$dist_lcp
#'
#' #### Example (7) Implement shortest distance algorithms in parallel:
#' # With the default method ("cppRouting"), use use_all_cores = TRUE
#' out_gebco1 <- lcp_over_surface(origin = origin,
#'                                destination = destination,
#'                                surface = dat_gebco_utm_planar,
#'                                method = "cppRouting",
#'                                crop = raster::extent(coastline),
#'                                mask = coastline,
#'                                mask_inside = TRUE,
#'                                use_all_cores = TRUE)
#' # With method = "gdistance" use cl argument
#' out_gebco2 <- lcp_over_surface(origin = origin,
#'                                destination = destination,
#'                                surface = dat_gebco_utm_planar,
#'                                method = "cppRouting",
#'                                crop = raster::extent(coastline),
#'                                mask = coastline,
#'                                mask_inside = TRUE,
#'                                cl = parallel::makeCluster(2L))
#' out_gebco1$dist_lcp; out_gebco2$dist_lcp
#'
#'
#' ############################################################################
#' #### Shortest pathways between a single origin and a single destination
#'
#' #### Example (8) Shortest pathways (goal = 2) only using default method
#' # Implement function
#' out1 <- lcp_over_surface(origin = rxy[1, , drop = FALSE],
#'                          destination = rxy[6, , drop = FALSE],
#'                          surface = r,
#'                          goal = 2)
#' # The path is stored in path_lcp, This includes:
#' # ... (a) a dataframe with the cells comprising each path:
#' # ... (b) a SpatialLines object of the path
#' # ... (c) a matrix of the coordinates of the path
#' out1$path_lcp
#' # For method = "cppRouting", paths between pairs of coordinates are computed
#' # ... by cppRouting::get_path_pair(), the arguments of which are retained in
#' # ... this list:
#' out1$cppRouting_param$get_path_pair_param
#' # Note the pathway is also added to the plot produced, if plot = TRUE.
#'
#' #### Example (9) Shortest distances and pathways (goal = 3) using default method
#' out1 <- lcp_over_surface(origin = rxy[1, , drop = FALSE],
#'                          destination = rxy[6, , drop = FALSE],
#'                          surface = r,
#'                          goal = 3)
#' out1$dist_lcp; out1$path_lcp
#'
#' #### Example (10) Shortest distances and pathways (goal = 3) via gdistance
#' out1 <- lcp_over_surface(origin = rxy[1, , drop = FALSE],
#'                          destination = rxy[6, , drop = FALSE],
#'                          surface = r,
#'                          goal = 3,
#'                          method = "gdistance")
#' # As above, paths are stored in path_lcp and distances in dist_lcp
#' out1$path_lcp
#' # For method = "gdistance", paths between pairs of coordinates are computed
#' # ... by repeated calls to gdistance::shortestPath(), for each pair of coordinates
#' # ... in the following list of arguments:
#' out1$gdistance_param$shortestPath_param
#'
#' #### Example (11): Parallelisation proceeds as described above via
#' # ... use_all_cores or cl arguments. For cppRouting, parallelisation
#' # ... is only implemented for distance calculations (so not if goal = 2),
#' # ... while parallelisation is implemented for both distance and shortest
#' # ... pathways for method = "gdistance"
#' out1 <- lcp_over_surface(origin = rxy[1, , drop = FALSE],
#'                          destination = rxy[6, , drop = FALSE],
#'                          surface = r,
#'                          goal = 3,
#'                          use_all_cores = TRUE)
#' out1 <- lcp_over_surface(origin = rxy[1, , drop = FALSE],
#'                          destination = rxy[6, , drop = FALSE],
#'                          surface = r,
#'                          goal = 2,
#'                          method = "gdistance",
#'                          cl = parallel::makeCluster(2L))
#' out1 <- lcp_over_surface(origin = rxy[1, , drop = FALSE],
#'                          destination = rxy[6, , drop = FALSE],
#'                          surface = r,
#'                          goal = 3,
#'                          method = "gdistance",
#'                          cl = parallel::makeCluster(2L))
#'
#'
#' ############################################################################
#' #### Shortest distances/pathways between origin/destination pairs
#'
#' #### Example (12): Shortest distances/pathways computed in sequence:
#' # cppRouting method
#' out1 <- lcp_over_surface(origin = rxy[1:2, , drop = FALSE],
#'                          destination = rxy[5:6, , drop = FALSE],
#'                          surface = r,
#'                          goal = 3)
#' # gdistance method
#' out2 <- lcp_over_surface(origin = rxy[1:2, , drop = FALSE],
#'                          destination = rxy[5:6, , drop = FALSE],
#'                          surface = r,
#'                          goal = 3,
#'                          method = "gdistance")
#' out1$dist_lcp; out1$path_lcp
#' out2$dist_lcp; out2$path_lcp
#'
#' #### Example (13): Shortest distances/pathways computed in parallel:
#' # cppRouting method for goal 3
#' out1 <- lcp_over_surface(origin = rxy[1:2, , drop = FALSE],
#'                          destination = rxy[5:6, , drop = FALSE],
#'                          surface = r,
#'                          goal = 3,
#'                          use_all_cores = TRUE)
#' # gdistance method for goal 2
#' out1 <- lcp_over_surface(origin = rxy[1:2, , drop = FALSE],
#'                          destination = rxy[5:6, , drop = FALSE],
#'                          surface = r,
#'                          goal = 2,
#'                          method = "gdistance",
#'                          cl = parallel::makeCluster(2L))
#' # gdistance method for goal 3
#' out1 <- lcp_over_surface(origin = rxy[1:2, , drop = FALSE],
#'                          destination = rxy[5:6, , drop = FALSE],
#'                          surface = r,
#'                          goal = 3,
#'                          method = "gdistance",
#'                          cl = parallel::makeCluster(2L))
#'
#'
#' ############################################################################
#' #### Shortest distances/pathways between all origin/destination combinations
#'
#' #### Example (14) Compute all combinations via combination = "matrix"
#' # cppRouting goal 3
#' out1 <- lcp_over_surface(origin = rxy[1:2, , drop = FALSE],
#'                          destination = rxy[5:6, , drop = FALSE],
#'                          surface = r,
#'                          goal = 3,
#'                          cppRouting_algorithm = "phast",
#'                          combination = "matrix")
#' out1$dist_euclid; out1$dist_lcp; out1$path_lcp
#' # cppRouting goal 3 parallelised
#' out1 <- lcp_over_surface(origin = rxy[1:2, , drop = FALSE],
#'                          destination = rxy[5:6, , drop = FALSE],
#'                          surface = r,
#'                          goal = 3,
#'                          combination = "matrix",
#'                          cppRouting_algorithm = "phast",
#'                          use_all_cores = TRUE)
#' # gdistance goal 3
#' out1 <- lcp_over_surface(origin = rxy[1:2, , drop = FALSE],
#'                          destination = rxy[5:6, , drop = FALSE],
#'                          surface = r,
#'                          goal = 3,
#'                          method = "gdistance",
#'                          combination = "matrix",
#'                          cl = parallel::makeCluster(2L))
#' out1$dist_euclid; out1$dist_lcp; out1$path_lcp
#'
#' #### Example (15): Real world example with multiple origins/destinations
#' ## Zoom in on an area of interest
#' ext <- raster::extent(715000, 720000, 6250000, 6260000)
#' dat_gebco_utm_planar_zoom <- raster::crop(dat_gebco_utm_planar, ext)
#' ## Define example origins/destinations
#' # Define available coordinates
#' dat_gebco_utm_planar_zoom_xy <- raster::coordinates(dat_gebco_utm_planar_zoom)
#' # Sample random origins
#' set.seed(2019)
#' index       <- sample(1:nrow(dat_gebco_utm_planar_zoom_xy), 2)
#' origin      <- dat_gebco_utm_planar_zoom_xy[index, ]
#' # Sample random destinations
#' set.seed(2020)
#' index       <- sample(1:nrow(dat_gebco_utm_planar_zoom_xy), 3)
#' destination <- dat_gebco_utm_planar_zoom_xy[index, ]
#' # Implement algorithm
#' out_gebco1 <- lcp_over_surface(origin = origin,
#'                                destination = destination,
#'                                surface = dat_gebco_utm_planar,
#'                                method = "gdistance",
#'                                goal = 3,
#'                                combination = "matrix")
#' out_gebco2 <- lcp_over_surface(origin = origin,
#'                                destination = destination,
#'                                surface = dat_gebco_utm_planar,
#'                                method = "cppRouting",
#'                                goal = 3,
#'                                combination = "matrix",
#'                                cppRouting_algorithm = "phast")
#' # The function returns distances between all combinations
#' out_gebco1$dist_lcp; out_gebco2$dist_lcp
#' # The two outputs are the same (though the ordering of factor levels differs)
#' str(out_gebco1$path_lcp$cells); str(out_gebco2$path_lcp$cells)
#'
#' @author Edward Lavender
#' @export
#'

lcp_over_surface <-
  function(origin,
           destination,
           surface,
           crop = NULL,
           buffer = NULL,
           aggregate = NULL,
           mask = NULL,
           mask_inside = TRUE,
           plot = TRUE,
           goal = 1,
           combination = "pair",
           method = "cppRouting",
           cppRouting_algorithm = "bi",
           cl = NULL, varlist = NULL, use_all_cores = FALSE,
           check = TRUE,
           verbose = TRUE,
           ...
  ){

    ########################################
    #### Set up raster

    #### Define a list of outputs
    out      <- list()
    out$time <- data.frame(event = "onset", time = Sys.time())
    out$args <- list(origin = origin,
                     destination = destination,
                     surface = surface,
                     crop = crop,
                     buffer = buffer,
                     aggregate = aggregate,
                     mask = mask,
                     mask_inside = mask_inside,
                     plot = plot,
                     goal = goal,
                     combination = combination,
                     method = method,
                     cppRouting_algorithm = cppRouting_algorithm,
                     cl = cl,
                     varlist = varlist,
                     use_all_cores = use_all_cores,
                     check = check,
                     verbose = TRUE,
                     dots = ...)

    #### Utils
    # Print messages to the console if verbose
    cat_to_console <- function(..., show = verbose){
      if(show) cat(paste(..., "\n"))
    }

    #### Checks
    if(check){
      cat_to_console("Checking user inputs...")
      # Check origin/destination coordinates are inputted as matrices
      check_class(input = origin, to_class = "matrix", type = "stop")
      check_class(input = destination, to_class = "matrix", type = "stop")
      # Check raster resolution is equal in the x and y directions
      surface_res <- raster::res(surface)
      if(!isTRUE(all.equal(surface_res[1], surface_res[2]))){
        stop("surface resolution must be equal in the x and y directions. Use raster::resample() to equalise resolution before implementing function.")
      }
      # Check input options for method, goal and combination
      method <- check_value(input = method, supp = c("cppRouting", "gdistance"))
      goal <- check_value(input = goal, supp = 1:3)
      combination <- check_value(input = combination, supp = c("pair", "matrix"))
      # For method = "pair", the number of origin and destination coordinates should be the same
      if(combination == "pair"){
        n_origin <- nrow(origin)
        n_destination <- nrow(destination)
        if(n_origin != n_destination){
          stop("Different numbers of origin and destination coordinates: combination = 'pair' cannot be implemented.")
        }
      }
      # Check algorithm specification for method = "cppRouting" because this depends on goal and combination input
      if(method == "cppRouting"){
        if(combination == "pair"){
          cppRouting_algorithm <- check_value(input = cppRouting_algorithm,
                                              supp = c("Dijkstra", "bi", "A*", "NBA"))
        } else if(combination == "matrix"){
          if(goal != 2){
            cppRouting_algorithm <- check_value(input = cppRouting_algorithm,
                                                supp = c("phast", "mch"))
          } else{
            message("'cppRouting_algorithm' input ignored for goal = 2 and combination = 'matrix' (cppRouting::get_multi_paths() implements Dijkstra's algorithm recursively.)")
          }
        }
      }
      # Check algorithm parallelisation options
      if(method == "cppRouting"){
        if(!is.null(cl) | !is.null(varlist)) message("'cl' or 'varlist' arguments ignored for method = 'cppRouting': use use_all_cores = TRUE instead.")
        if(goal == 2 & use_all_cores) message("'use_all_cores = TRUE' only implemented for shortest distances (goal = 1 or goal = 3).")
      }
    }

    #### Crop raster
    if(!is.null(crop)){
      cat_to_console("Cropping raster to inputted extent...")
      surface <- raster::crop(surface, crop)
    }

    #### Crop raster around Euclidean transect(s) for speed
    if(!is.null(buffer)){
      cat_to_console("Cropping raster to buffer zone around a Euclidean transect between the origin and destination...")
      # Check only one origin/destination node
      if(nrow(origin) > 1 | nrow(destination) > 1){
        stop("More than one set of origin or destination coordinates inputted; 'buffer' cannot be implemented.")
      } else{
        # Define Euclidean line meeting two points
        connector <- Orcs::coords2Lines(rbind(origin, destination), ID = 1)
        # Define buffer around connector
        buffer$spgeom <- connector
        connector_buf <- do.call(rgeos::gBuffer, buffer)
        # Crop surface around buffer
        surface <- raster::crop(surface, raster::extent(connector_buf))
      }
    }

    #### Aggregate raster
    if(!is.null(aggregate)){
      cat_to_console("Aggregating raster...")
      aggregate$x <- surface
      surface <- do.call(raster::aggregate, aggregate)
    }

    #### Mask raster
    if(!is.null(mask)){
      cat_to_console("Masking raster...")
      surface <- mask_io(x = surface, mask = mask, mask_inside = mask_inside, updatevalue = Inf)
    }

    #### Visualise surface
    if(plot){
      ## Define plotting window to fit two plots
      pp <- graphics::par(mfrow = c(1, 2))
      ## Visualise inputted raster
      # Inputted raster
      raster::plot(out$args$surface, xlab = "x", ylab = "y", main = "Surface (inputted)")
      # Highlight mask, if inputted
      if(!is.null(mask)){
        raster::plot(mask, add = TRUE)
      }
      # Highlight cropped region, if inputted
      if(!is.null(crop)){
        raster::lines(crop, col = "red", lwd = 2)
      }
      # Highlight buffered region, if inputted
      if(!is.null(buffer)){
        raster::lines(connector, col = "red", lwd = 2)
        raster::lines(connector_buf, col = "red", lwd = 2, lty = 3)
      }
      # Add origin/destination
      graphics::points(rbind(origin, destination), pch = 21, bg = "blue", col = "blue")
      ## Visualise processed raster with origin and destination added
      raster::plot(surface, xlab = "x", ylab = "y", main = "Surface (processed)")
      graphics::points(rbind(origin, destination), pch = 21, bg = "blue", col = "blue")
      if(all(goal %in% 1)) graphics::par(pp)
    }

    #### Check coordinates lie within raster extent and are not NA or Inf
    if(check){
      coords <- rbind(origin, destination)
      ext <- raster::extent(surface)
      if(min(coords[, 1]) < ext[1] | max(coords[, 1]) > ext[2]) stop("Coordinates exceed x-limits of processed surface.")
      if(min(coords[, 2]) < ext[3] | max(coords[, 2]) > ext[4]) stop("Coordinates exceed y-limits of processed surface.")
      coord_surface_vals <- raster::extract(surface, coords)
      if(any(is.na(coord_surface_vals)) | any(is.infinite(coord_surface_vals))) stop("Coordinate(s) over NA or Inf cell on processed surface.")
    }

    #### Add processed surface to list of objects to be returned
    out$surface <- surface
    out$time <- rbind(out$time, data.frame(event = "surface_processed", time = Sys.time()))


    ########################################
    #### Set up transition matrix

    #### Define cell connections and distances between adjacent cells
    cat_to_console("Defining cost matrix...")

    #### Define parameters
    # Basis param
    nrw <- raster::nrow(surface)
    ncl <- raster::ncol(surface)
    cells <- 1:raster::ncell(surface)
    rxy   <- raster::coordinates(surface)
    origin_cell <- raster::cellFromXY(surface, origin)
    destination_cell <- raster::cellFromXY(surface, destination)
    # Identify any duplicate pairs
    if(check){
      if(combination == "pair"){
        odcells        <- data.frame(origin_cell = origin_cell, destination_cell = destination_cell)
        odcells$couple <- paste0(odcells$origin_cell, "-", odcells$destination_cell)
        odcells$dup    <- duplicated(odcells$couple)
        if(any(odcells$dup)){
          pos <- which(odcells$dup)
          message(paste0("Duplicated origin/destination pair(s) (i.e., those with identical cell IDs) identified: ", paste(pos, collapse = ", "), ". Removing these may improve algorithm speed."))
        }
      }
    }

    # Add param to out list
    out$surface_param <- list(cells = cells,
                              nrow = nrw,
                              ncol = ncl,
                              coordinates = rxy,
                              origin_cell = origin_cell,
                              destination_cell = destination_cell)

    #### Define functions to compute distances between connected of cells
    calc_height <- function(i) i[2] - i[1]
    lrook   <- raster::res(surface)[1]
    calc_rook <- function(i) lrook
    lbishop <- sqrt(lrook^2 + lrook^2)
    calc_bishop <- function(i) lbishop

    #### Define transition matrices with the distances between connected cells
    tr <- gdistance::transition(surface, calc_height, directions = 8, symm = TRUE)
    tr_rook <- gdistance::transition(surface, calc_rook, directions = 4, symm = TRUE)
    tr_bishop <- gdistance::transition(surface, calc_bishop, directions = "bishop", symm = TRUE)

    ## Extract transitionMatrix with distances
    # vertical distances between connected cells
    vdist_mat  <- gdistance::transitionMatrix(tr)
    # planar rook distances between connected cells
    rook_mat   <- gdistance::transitionMatrix(tr_rook)
    # planar bishop distances between connected cells
    bishop_mat <- gdistance::transitionMatrix(tr_bishop)
    # total planar distances between connected cells
    hdist_mat  <- rook_mat + bishop_mat
    # total distances between connected cells
    tdist_mat <- sqrt(hdist_mat^2 + vdist_mat^2)

    ## Define cost components/surfaces in output matrix
    out$cost <- list(dist_rook = rook_mat,
                     dist_bishop = bishop_mat,
                     dist_planar = hdist_mat,
                     dist_vertical = vdist_mat,
                     dist_total = tdist_mat
    )
    out$time <- rbind(out$time, data.frame(event = "transition_matrix_defined", time = Sys.time()))


    ########################################
    #### Implement LCP algorithm

    #### Calculate Euclidean distances (for comparison to LCP distances)
    if(goal %in% c(1, 3)){
      cat_to_console("Calculating Euclidean distance(s)...")
      if(combination == "pair") allpairs <- FALSE else if(combination == "matrix") allpairs <- TRUE
      dist_euclid <- raster::pointDistance(origin, destination, lonlat = FALSE, allpairs = allpairs)
      if(nrow(origin) == 1 & nrow(destination) == 1) dist_euclid <- as.numeric(dist_euclid)
      if(allpairs){
        rownames(dist_euclid) <- origin_cell
        colnames(dist_euclid) <- destination_cell
      }
      out$dist_euclid <- dist_euclid
    }

    #### LCP distances/pathways
    if(method == "gdistance"){
      cat_to_console("Using method = 'gdistance'...")

      #### Calculate the reciprocal of the costs
      # ... (i.e. the total ease of movement between connected cells):
      cat_to_console("... Defining 'ease' matrix...")
      cells_non_zero <- Matrix::which(tdist_mat > 0, arr.ind = TRUE)
      ease_mat <- tdist_mat
      ease_mat[cells_non_zero] <- 1/ease_mat[cells_non_zero]

      #### Finalise transition object:
      # replace the transitionMatrix slot with the properly computed
      # ... values for the ease_mat of movement between adjacent cells:
      tr@transitionMatrix <- ease_mat

      #### Define a list of origin and destination coordinates to loop over
      # ... if combination == "pair"
      if(combination == "pair"){
        # Join origin and destination coordinates in a single matrix
        coords <- cbind(origin, destination)
        # Define a list of coordinates, with one element for each pair
        coords_by_pair <- lapply(seq_len(nrow(coords)), function(i) coords[i,, drop = FALSE])
      }

      #### Calculate least-cost distance between points
      if(goal %in% c(1, 3)){

        ## Define parameters to compute least-cost distances
        cat_to_console("... Implementing Dijkstra's algorithm to compute least-cost distance(s)...")
        costDistance_param <- list(x = tr, fromCoords = origin, toCoords = destination)
        out$gdistance_param <- list(costDistance_param = costDistance_param)

        ## Least cost distances between pairs
        if(combination == "pair"){

          # Set up cluster to loop over each pair of coordinates
          if(!is.null(cl) & !is.null(varlist)) {
            parallel::clusterExport(cl = cl, varlist = varlist)
          }
          # Loop over each pair of coordinates and compute least-cost distances
          dist_lcp_by_pair <- pbapply::pblapply(coords_by_pair, cl = cl, FUN = function(xy){
            costDistance_param <- list(x = tr, fromCoords = xy[, 1:2, drop = FALSE], toCoords = xy[, 3:4, drop = FALSE])
            dist_lcp <- do.call(gdistance::costDistance, costDistance_param)
            return(dist_lcp)
          })
          if(all(goal %in% 1)) if(!is.null(cl)) parallel::stopCluster(cl = cl)
          # Process least-cost distances and close cluster
          dist_lcp <- as.numeric(unlist(dist_lcp_by_pair))


          ## Least cost differences among all combinations of points
          # ... (The standard implementation of costDistance)
        } else if(combination == "matrix"){
          dist_lcp <- do.call(gdistance::costDistance, costDistance_param)
          rownames(dist_lcp) <- origin_cell
          colnames(dist_lcp) <- destination_cell
        }
        out$time <- rbind(out$time, data.frame(event = "lcp_dist_defined", time = Sys.time()))
      }

      #### Calculate least-cost pathways between points
      if(goal %in% c(2, 3)){
        cat_to_console("... Implementing Dijkstra's algorithm to compute least-cost pathway(s)...")
        shortestPath_param <- list(x = tr, fromCoords = origin, toCoords = destination)
        out$gdistance_param$shortestPath_param <- shortestPath_param

        ## Set up coordinates combinations to loop over (if necessary)
        if(combination == "matrix"){
          # Define all, non duplicate, combinations of coordinates in a dataframe
          origin_index         <- 1:nrow(origin)
          destination_index    <- 1:nrow(destination)
          coords               <- expand.grid(origin_index = origin_index, destination_index = destination_index)
          coords$origin_x      <- origin[coords$origin_index, 1]
          coords$origin_y      <- origin[coords$origin_index, 2]
          coords$destination_x <- destination[coords$destination_index, 1]
          coords$destination_y <- destination[coords$destination_index, 2]
          coords$origin        <- paste0(coords$origin_x, ",", coords$origin_y)
          coords$destination   <- paste0(coords$destination_x, ",", coords$destination_y)
          coords$pair          <- NA
          for(i in 1:nrow(coords)){
            coords$pair[i]       <- paste(sort(c(coords$origin[i], coords$destination[i])), collapse = "-")
          }
          coords <- coords[coords$origin != coords$destination, ]
          coords <- coords[!duplicated(coords$pair), ]
          coords <- cbind(coords$origin_x, coords$origin_y, coords$destination_x, coords$destination_y)
          coords_by_pair <- lapply(seq_len(nrow(coords)), function(i) coords[i,, drop = FALSE])
        }

        ## Set up cluster to loop over each pair of coordinates
        if(!is.null(cl) & !is.null(varlist)) {
          parallel::clusterExport(cl = cl, varlist = varlist)
        }

        ## Loop over each pair of coordinates and compute least-cost pathways
        path_lcp_SpatialLines <- pbapply::pblapply(coords_by_pair, cl = cl, FUN = function(xy){
          param <- list(x = tr, origin = xy[, 1:2, drop = FALSE], goal = xy[, 3:4, drop = FALSE], output = "SpatialLines")
          path_lcp <- do.call(gdistance::shortestPath, param)
          return(path_lcp)
        })
        if(!is.null(cl)) parallel::stopCluster(cl = cl)

        ## Process least-cost pathways
        # Define dataframe with paths and cells
        path_lcp_cells <- mapply(path_lcp_SpatialLines, coords_by_pair, FUN = function(l, dcoords){
          xy <- raster::geom(l)[, c("x", "y")]
          cells <- raster::cellFromXY(surface, xy)
          origin_cells <- raster::cellFromXY(surface, dcoords[, 1:2, drop = FALSE])
          destination_cells <- raster::cellFromXY(surface, dcoords[, 3:4, drop = FALSE])
          dcells <- data.frame(origin = origin_cells, destination = destination_cells, cell = cells)
          return(dcells)
        }, SIMPLIFY = FALSE) %>% dplyr::bind_rows()
        path_lcp_cells$origin      <- as.integer(path_lcp_cells$origin)
        path_lcp_cells$destination <- as.integer(path_lcp_cells$destination)
        path_lcp_cells$cell        <- as.integer(path_lcp_cells$cell)
        path_lcp_cells$path <- paste0(path_lcp_cells$origin, "-", path_lcp_cells$destination)
        path_lcp_cells$path <- factor(path_lcp_cells$path, levels = unique(path_lcp_cells$path))
        # Define a unique identifier for each path (that distinguishes even duplicated paths)
        pos_start <- which(path_lcp_cells$origin == path_lcp_cells$cell)
        if(length(pos_start) >= 2){
          path_lcp_cells$path_id     <- NA
          for(i in 1:(length(pos_start) - 1)){
            path_lcp_cells$path_id[pos_start[i]:pos_start[i+1]] <- i
          }
          path_lcp_cells$path_id[pos_start[length(pos_start)]:nrow(path_lcp_cells)] <- i+1
        } else{
          path_lcp_cells$path_id <- 1
        }
        path_lcp_cells <- path_lcp_cells[, c("path_id", "path", "origin", "destination", "cell")]
        rownames(path_lcp_cells) <- 1:nrow(path_lcp_cells)
        # Define list of coordinates for each path
        path_lcp_coordinates <- lapply(path_lcp_SpatialLines, function(l){
          xy <- raster::geom(l)[, c("x", "y")]
          return(xy)
        })
        out$time <- rbind(out$time, data.frame(event = "path_lcp_defined", time = Sys.time()))

      }


    } else if(method == "cppRouting"){
      cat_to_console("Using method = 'cppRouting'...")

      #### Define a dataframe of edges with associated costs
      cat_to_console("... Defining nodes, edges and costs to make graph...")
      # from, to, cost
      # edges <- expand.grid(from = cells, to = cells)
      # edges$cost <- as.vector(tdist_mat)
      # edges <- edges[edges$from != edges$to, ]
      # edges <- edges[which(edges$cost != 0), ]

      non_zero_indices <- Matrix::which(tdist_mat != 0, arr.ind = TRUE)
      edges <- data.frame(from = cells[non_zero_indices[, 1]],
                          to = cells[non_zero_indices[, 2]],
                          cost = as.vector(tdist_mat[non_zero_indices]))

      #### Specify a dataframe of coordinates for each node in appropriate format
      # ... redefine origin/destination nodes in terms of cell numbers
      coord <- data.frame(node = cells, X = rxy[, 1], Y = rxy[, 2])

      #### Make graph
      cat_to_console("... Constructing graph object...")
      makegraph_param <- list(df = edges, directed = TRUE, coords = coord)
      out$cppRouting_param <- list(makegraph_param = makegraph_param)
      graph <- do.call(cppRouting::makegraph, makegraph_param)
      out$time <- rbind(out$time, data.frame(event = "graph_defined", time = Sys.time()))

      #### Compute all shortest distance between origin and destination nodes
      if(goal %in% c(1, 3)){
        # Define parameters to compute distances
        cat_to_console(paste("... Implementing",  cppRouting_algorithm, "algorithm to compute least-cost distance(s)..."))
        get_distance_param <- list(Graph = graph,
                                   from = origin_cell,
                                   to = destination_cell,
                                   algorithm = cppRouting_algorithm,
                                   allcores = use_all_cores)
        # Define function to compute distances and add parameters to output list
        if(combination == "pair"){
          get_distance <- cppRouting::get_distance_pair
          out$cppRouting_param$get_distance_pair_param <- get_distance_param
        } else if(combination == "matrix"){
          get_distance <- cppRouting::get_distance_matrix
          out$cppRouting_param$get_distance_matrix_param <- get_distance_param
        }
        # Compute distances
        dist_lcp <- do.call(get_distance, get_distance_param)
        out$time <- rbind(out$time, data.frame(event = "dist_lcp_defined", time = Sys.time()))
      }

      #### Compute shortest pathways between origin and destination nodes
      if(goal %in% c(2, 3)){

        ## Set up to compute shortest pathways
        # Define function to compute paths and add parameters to output list
        if(combination == "pair"){
          cat_to_console(paste("... Implementing",  cppRouting_algorithm, "algorithm to compute least-cost pathways(s)..."))
          get_path_param <- list(Graph = graph,
                                 from = origin_cell,
                                 to = destination_cell,
                                 algorithm = cppRouting_algorithm,
                                 constant = 1,
                                 keep = NULL,
                                 long = TRUE)
          get_path <- cppRouting::get_path_pair
          out$cppRouting_param$get_path_pair_param <- get_path_param
        } else if(combination == "matrix"){
          cat_to_console("... Implementing Dijkstra's algorithm recursively to compute least-cost pathways(s)...")
          get_path_param <- list(Graph = graph,
                                 from = origin_cell,
                                 to = destination_cell,
                                 keep = NULL,
                                 long = TRUE)
          get_path <- cppRouting::get_multi_paths
          out$cppRouting_param$get_multi_paths_param <- get_path_param
        }

        ## Compute pathways (origin, destination and cells)
        path_lcp_cells <- do.call(get_path, get_path_param)
        out$time <- rbind(out$time, data.frame(event = "path_lcp_defined", time = Sys.time()))

        ## Process pathways dataframe
        # Define columns
        colnames(path_lcp_cells)   <- c("origin", "destination", "cell")
        path_lcp_cells$origin      <- as.integer(path_lcp_cells$origin)
        path_lcp_cells$destination <- as.integer(path_lcp_cells$destination)
        path_lcp_cells$cell        <- as.integer(path_lcp_cells$cell)
        path_lcp_cells$path <- paste0(path_lcp_cells$origin, "-", path_lcp_cells$destination)
        path_lcp_cells$path <- factor(path_lcp_cells$path, levels = unique(path_lcp_cells$path))
        # Define a unique identifier for each path (that distinguishes even duplicated paths)
        pos_start <- which(path_lcp_cells$destination == path_lcp_cells$cell)
        if(length(pos_start) >= 2){
          path_lcp_cells$path_id     <- NA
          for(i in 1:(length(pos_start) - 1)){
            path_lcp_cells$path_id[pos_start[i]:pos_start[i+1]] <- i
          }
          path_lcp_cells$path_id[pos_start[length(pos_start)]:nrow(path_lcp_cells)] <- i+1
        } else{
          path_lcp_cells$path_id <- 1
        }
        path_lcp_cells <- path_lcp_cells[, c("path_id", "path", "origin", "destination", "cell")]
        # For cppRouting, nodes are reported from the end to the start, so reverse the order
        path_lcp_cells <- lapply(split(path_lcp_cells, path_lcp_cells$path_id), function(d) d[nrow(d):1, ]) %>% dplyr::bind_rows()
        rownames(path_lcp_cells) <- 1:nrow(path_lcp_cells)

        # Define pathway coordinates and SpatialLines
        path_lcp_spatial <- lapply(split(path_lcp_cells, path_lcp_cells$path_id), function(d){
          xy  <- raster::xyFromCell(surface, cell = d$cell)
          spl <- Orcs::coords2Lines(xy, ID = d$path_id[1], proj4string = raster::crs(surface))
          return(list(coordinates = xy, SpatialLines = spl))
        })
        path_lcp_coordinates  <- lapply(path_lcp_spatial, function(elm) elm$coordinates)
        path_lcp_SpatialLines <- lapply(path_lcp_spatial, function(elm) elm$SpatialLines)
      }
    }

    #### Add least-cost distances to output object
    if(goal %in% c(1, 3)){
      if(nrow(origin) == 1 & nrow(destination) == 1) dist_lcp <- as.numeric(dist_lcp)
      out$dist_lcp <- dist_lcp
    }

    #### Add least-cost pathways to output object
    if(goal %in% c(2, 3)){
      # pathways
      names(path_lcp_SpatialLines) <- unique(path_lcp_cells$path_id)
      names(path_lcp_coordinates)  <- unique(path_lcp_cells$path_id)
      out$path_lcp <- list(cells = path_lcp_cells,
                           SpatialLines = path_lcp_SpatialLines,
                           coordinates = path_lcp_coordinates
      )
      # Add pathways to plot, if requested
      if(plot){
        lapply(path_lcp_SpatialLines, function(l) raster::lines(l, col = "royalblue", lwd = 2))
        graphics::par(pp)
      }
    }

    #### Return outputs
    out$time <- rbind(out$time, data.frame(event = "finish", time = Sys.time()))
    out$time$serial_duration <- Tools4ETS::serial_difference(out$time$time)
    cat_to_console("Done.")
    return(out)

  }
