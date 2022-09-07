######################################
######################################
#### lcp_costs()

#' @title Calculate the distances between connected cells in a Raster*
#' @description This function calculates distances between connected cells in a \code{\link[raster]{raster}}, given (a) the planar distances between connected cells and (b) differences in elevation.
#'
#' @param surface A \code{\link[raster]{raster}} for which to calculate distances. The \code{surface} must be planar (i.e., Universal Transverse Mercator projection) with units of metres in the x, y and z directions. The \code{surface}'s \code{\link[raster]{resolution}} is taken to define the distance between connected cells in the x and y directions and must be the same in both cases (for \code{surface}'s with unequal horizontal resolution, \code{\link[raster]{resample}} can be used to equalise resolution).
#' @param verbose A logical function that defines whether or not to print messages to the console to relay function progress.
#'
#' @details This function was motivated by the need to determine the shortest paths between locations over the seabed for benthic animals (see \code{\link[flapper]{lcp_over_surface}}). An animal's movement can be conceptualised as that of a queen on a chessboard, which can move, in eight directions around its current position, across a surface. Movements in the x and y direction are termed `rook's movement' and movements in the diagonal direction are termed `bishop's movement'.
#'
#' Under this framework, the distance that an entity must travel between connected cells depends on the planar distances between cells and their differences in elevation. Planar distances (\eqn{d_p}, m) depend on the movement type: under a rook's movement (i.e., horizontally or vertically), the distance (\eqn{d_{p,r}}) between connected cells is extracted from the raster's resolution (which is assumed to be identical in the x and y directions); under a bishop's movement (i.e., diagonally), the distance between connected cells \eqn{d_{p,b}} is given by Pythagoras' Theorem: \eqn{d_{p,b} = \sqrt{(d_{p, r}^2 + d_{p, r}^2)}}. Vertical distances (\eqn{d_v}, m) are simply the differences in height between cells. The total distance (\eqn{d_t}) between any two connected cells is a combination of these distances given by Pythagoras' Theorem: \eqn{d_t = \sqrt{(d_p^2 + d_v^2)}}.
#'
#' The function returns a warning produced by \code{\link[gdistance]{transition}} which is implemented to facilitate the definition of distances before shortest paths/distances are computed by either method: `In .TfromR(x, transitionFunction, directions, symm) : transition function gives negative values'. This warning arises because the height differences between connecting cells can be negative. It can be safely ignored.
#'
#' @return The function returns a named list of sparse \code{\link[Matrix]{dsCMatrix-class}} matrices that define the distances (m) between connected cells under a rook's or bishop's movement (`dist_rook' and `dist_bishop'), the planar and vertical distances between connected cells (`dist_planar' and `dist_vertical') and the total distance between connected cells (`dist_total').
#'
#' @examples
#' # In this example, consider the distances between connected cells in the example
#' # ... 'dat_gebco' raster. For speed, we will focus on a subset of the area.
#' # ... Within this area, we need to regularise the resolution:
#' boundaries <- raster::extent(707884.6, 709884.6, 6253404, 6255404)
#' blank      <- raster::raster(boundaries, res = c(5, 5))
#' r          <- raster::resample(dat_gebco, blank)
#' # Implement algorithm
#' costs      <- lcp_costs(r)
#' # Examine outputs
#' utils::str(costs)
#'
#' @seealso This routine is implemented by the \code{\link[flapper]{lcp_over_surface}} function to calculate the shortest path(s) and/or the distance(s) of the shortest paths(s) between origin and destination coordinates, and by the \code{\link[flapper]{lcp_from_point}} function to calculate the shortest distances from a point on a \code{\link[raster]{raster}} to surrounding cells.
#' @author Edward Lavender
#' @export
#'

lcp_costs <- function(surface, verbose = TRUE){

  #### Define functions to compute distances between connected of cells
  cat_to_console <- function(..., show = verbose){
    if(show) cat(paste(..., "\n"))
  }
  t_onset <- Sys.time()
  cat_to_console(paste0("flapper::lcp_costs() called (@ ", t_onset, ")..."))
  if(!all.equal(raster::res(surface)[1], raster::res(surface)[2])){
    stop("The resolution of the surface should be equal in the x and y directions. Consider calling raster::resample() before implementing this function.")
  }
  calc_height <- function(i) i[2] - i[1]
  lrook   <- raster::res(surface)[1]
  calc_rook <- function(i) lrook
  lbishop <- sqrt(lrook^2 + lrook^2)
  calc_bishop <- function(i) lbishop

  #### Define transition matrices with the distances between connected cells
  cat_to_console(paste0("... Defining transition matrices..."))
  tr <- gdistance::transition(surface, calc_height, directions = 8, symm = TRUE)
  tr_rook <- gdistance::transition(surface, calc_rook, directions = 4, symm = TRUE)
  tr_bishop <- gdistance::transition(surface, calc_bishop, directions = "bishop", symm = TRUE)

  ## Extract transitionMatrix with distances
  cat_to_console(paste0("... Calculating distance matrices..."))
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
  cat_to_console(paste0("... Assembling LCP costs..."))
  out <- list(dist_rook = rook_mat,
              dist_bishop = bishop_mat,
              dist_planar = hdist_mat,
              dist_vertical = vdist_mat,
              dist_total = tdist_mat
  )

  #### Return outputs
  t_end <- Sys.time()
  total_duration <- difftime(t_end, t_onset, units = "mins")
  cat_to_console(paste0("... flapper::lcp_costs() call completed (@ ", t_end, ") after ~", round(total_duration, digits = 2), " minutes."))
  return(out)
}


######################################
######################################
#### lcp_graph_surface()

#' @title Construct a graph for least-cost paths analysis of movement from a point on a Raster*
#' @description This function constructs a graph for least-cost paths analysis. This function is specifically designed for situations in which the aim is to calculate least-cost distances between points on a \code{\link[raster]{raster}}.
#' @param surface A \code{\link[raster]{raster}} for which to construct the graph. There are some constraints on the form of this raster if the costs of movement between cells are derived from \code{\link[flapper]{lcp_costs}} (see below).
#' @param cost A sparse \code{\link[Matrix]{dsCMatrix-class}} matrix that defines the cost of movement between connected cells on the \code{surface} (see \code{\link[flapper]{lcp_costs}}).
#' @param verbose A logical function that defines whether or not to print messages to the console to relay function progress.
#' @details This is a wrapper for the \code{\link[cppRouting]{makegraph}} function.
#' @return The function returns a named list that defines the graph (see \code{\link[cppRouting]{makegraph}})
#' @examples
#' #### Step (1): Define cost surface
#' # We will consider the distances between connected cells in the example
#' # ... 'dat_gebco' raster as a measure of cost. For this, we will focus on a
#' # ... specific area (for speed), within which we need to
#' # ... we need to regularise the resolution:
#' boundaries <- raster::extent(707884.6, 709884.6, 6253404, 6255404)
#' blank      <- raster::raster(boundaries, res = c(5, 5))
#' r          <- raster::resample(dat_gebco, blank)
#' # Define costs
#' costs      <- lcp_costs(r)
#'
#' #### Step (2): Make graph for LCP analysis
#' graph <- lcp_graph_surface(surface = r, cost = costs$dist_total)
#'
#' #### Step (3): Implement LCP analysis around point on Raster*
#' # ... e.g., via lcp_from_point()
#'
#' @author Edward Lavender
#' @export

lcp_graph_surface <- function(surface, cost, verbose = TRUE){

  #### Set up function
  cat_to_console <- function(..., show = verbose){
    if(show) cat(paste(..., "\n"))
  }
  t_onset <- Sys.time()
  cat_to_console(paste0("flapper::lcp_graph_surface() called (@ ", t_onset, ")..."))
  op <- options()
  options(scipen = 999)
  on.exit(options(op), add = TRUE)
  cells <- 1:raster::ncell(surface)
  surface_df <- data.frame(cell = cells, value = as.vector(surface))

  #### Define a dataframe of edges with associated costs
  cat_to_console("... Defining nodes, edges and costs to make graph...")
  non_zero_indices <- Matrix::which(cost != 0, arr.ind = TRUE)
  edges <- data.frame(from = cells[non_zero_indices[, 1]],
                      to = cells[non_zero_indices[, 2]],
                      cost = as.vector(cost[non_zero_indices]))

  #### Specify a dataframe of coordinates for each node in appropriate format
  # ... redefine origin/destination nodes in terms of cell numbers
  rxy   <- raster::coordinates(surface)
  coord <- data.frame(node = cells, X = rxy[, 1], Y = rxy[, 2])

  #### Make graph
  cat_to_console("... Constructing graph object...")
  makegraph_param <- list(df = edges, directed = FALSE, coords = coord)
  graph <- do.call(cppRouting::makegraph, makegraph_param)
  t_end <- Sys.time()
  total_duration <- difftime(t_end, t_onset, units = "mins")
  cat_to_console(paste0("... flapper::lcp_graph_surface() call completed (@ ", t_end, ") after ~", round(total_duration, digits = 2), " minutes."))
  return(graph)
}


######################################
######################################
#### lcp_from_point()

#' @title Create a Raster* of the least-cost distances around a point
#' @description This function calculates the least-cost distances from a point to all (or some) of the cells in a surrounding \code{\link[raster]{raster}} object, returning a \code{\link[raster]{raster}}. This is the least-cost distance equivalent of \code{\link[raster]{distanceFromPoints}}.
#'
#' @param origin A matrix which defines the coordinates (x, y) of the point from which to calculate least-cost distances. Unlike \code{\link[raster]{distanceFromPoints}}, only a single point is expected.
#' @param surface A \code{\link[raster]{raster}} across which to implement least-cost distance calculations. If the \code{cost} matrix is derived from \code{\link[flapper]{lcp_costs}} (see below), there are some constraints on the form of this \code{surface}; namely, equal resolution in x and y directions and a Universal Transverse Mercator coordinate reference system with units of metres. The \code{surface} defines the properties of the returned \code{\link[raster]{raster}} (see Value).
#' @param destination (optional) An matrix of destination coordinates; an integer vector of cell IDs; or function that defines a subset of destination cells, given their \code{surface} value, for which to implement calculations. For example \code{destination = function(x) x > 0} would restrict least-cost distance calculations to cells of the \code{surface} that have a value of more than zero. Other cells are set to NA. This can improve computational efficiency.
#' @param cost (optional) A sparse \code{\link[Matrix]{dsCMatrix-class}} matrix that defines the cost of movement between connected cells (see \code{\link[flapper]{lcp_costs}}). If unsupplied, if the \code{graph} is also unsupplied (see below), a matrix of distances from \code{\link[flapper]{lcp_costs}} is computed internally and taken to define the cost surface. For this to be appropriate, the \code{surface} should have a Universal Transverse Mercator projection, with equal resolution in the x and y directions and units of metres (see \code{surface}, above). If a \code{graph} is supplied, \code{cost} is unnecessary.
#' @param graph (optional) A graph object that defines cell nodes and edge costs for connected cells within the \code{surface} (see \code{\link[flapper]{lcp_graph_surface}}). If supplied, the calculation of the cost surface and the construction of the graph stages in the computation of least-cost distances are skipped (see Details), which is desirable in iterative applications.
#' @param verbose A logical input that defines whether or not to print messages to the console to relay function progress.
#' @param use_all_cores A logical input that defines whether or not to parallelise least-cost distance calculations across all cores. This is passed to \code{\link[cppRouting]{get_distance_matrix}} which implements calculations.
#'
#' @details This function implements routines provided via \code{\link[flapper]{flapper}} and the \code{\link[cppRouting]{cppRouting}} package to calculate least-cost distances. The main steps are:
#' \enumerate{
#'   \item The calculation of distances between adjacent cells (i.e., \code{cost}, if not supplied, via \code{\link[flapper]{lcp_costs}});
#'   \item The construction of a graph that defines cell connections from the \code{origin} to surrounding cells on the \code{surface} as a network (via \code{\link[cppRouting]{makegraph}});
#'   \item The calculation of shortest distances between the \code{origin} and surrounding cells from the graph (via \code{\link[cppRouting]{get_distance_matrix}});
#'   \item The expression of shortest distances as a \code{\link[raster]{raster}} which is returned.
#' }
#'
#' @return The function returns a \code{\link[raster]{raster}} in which each cell represents the least-cost distance from a specified \code{origin} to that cell. The \code{origin} is assigned a value of zero. Any cells excluded by the \code{destination} filter have a value of NA.
#'
#' @examples
#' #### Step (1): Define example origin
#' proj     <- sp::CRS(SRS_string = "EPSG:4326")
#' proj_utm <- sp::CRS(SRS_string = "EPSG:32629")
#' origin   <- matrix(c(-5.616, 56.388), ncol = 2)
#' origin   <- sp::SpatialPoints(origin, proj)
#' origin   <- sp::spTransform(origin, proj_utm)
#'
#' #### Step (2): Select and process surface
#' # We will focus on an area within the dat_gebco bathymetry raster
#' boundaries <- raster::extent(707884.6, 709884.6, 6253404, 6255404)
#' blank      <- raster::raster(boundaries, res = c(5, 5))
#' r          <- raster::resample(dat_gebco, blank)
#'
#' #### Example (1): Implement function using default options
#' lcp_dist <- lcp_from_point(origin = origin, surface = r)
#' ## Visualise outputs
#' pp <- par(mfrow = c(2, 2))
#' # Plot surface
#' raster::plot(r)
#' # Examine Euclidean distances from point
#' raster::plot(raster::distanceFromPoints(r, origin))
#' # Compare to shortest distances
#' raster::plot(lcp_dist)
#' par(pp)
#'
#' if(flapper_run_slow){
#'
#' #### Example (2): Implement function across specific destinations
#' ## Supply destination cell coordinates/IDs directly
#' # E.g., consider distances to cells shallower than 125 m
#' destination_cells <- raster::Which(r < 125, cells = TRUE, na.rm = TRUE)
#' lcp_dist <- lcp_from_point(origin = origin,
#'                            surface = r,
#'                            destination = destination_cells)
#' raster::plot(lcp_dist)
#' ## Use a function instead to consider distances to cells shallower than 125 m
#' filter_destination_cells <- function(x) x < 125
#' lcp_dist <- lcp_from_point(origin = origin,
#'                            surface = r,
#'                            destination = filter_destination_cells)
#' raster::plot(lcp_dist)
#'
#' #### Example (3): Define cost surfaces for LCP calculations outside of function
#' # This can be implemented internally, but we compute it here via lcp_costs().
#' # Note this imposes restrictions on the nature of the surface, such as equal
#' # ... resolution, which we have forced above.
#' costs <- lcp_costs(r)
#' cost  <- costs$dist_total
#' lcp_dist <- lcp_from_point(origin = origin,
#'                            surface = r,
#'                            destination = filter_destination_cells,
#'                            cost = cost)
#'
#' #### Example (4): Supply a graph object
#' graph    <- lcp_graph_surface(surface = r, cost = cost)
#' lcp_dist <- lcp_from_point(origin = origin, surface = r, graph = graph)
#'
#' #### Example (5): Implement algorithm in parallel via use_all_cores
#' lcp_dist <- lcp_from_point(origin = origin, surface = r, use_all_cores = TRUE)
#'
#' }
#'
#' @seealso This function is similar to \code{\link[raster]{distanceFromPoints}}, which returns a Raster* of Euclidean distances. For iterative applications across the same surface, \code{\link[flapper]{lcp_costs}} and \code{\link[flapper]{lcp_graph_surface}} can be implemented to define the cost matrix and the graph object outside of this function. These can be passed to \code{\link[flapper]{lcp_from_point}}, skipping the need to recompute these objects. For shortest-distances and/or paths between specific origin and destination coordinates, the \code{\link[flapper]{lcp_over_surface}} can be used. The particle filtering movement algorithms in \code{\link[flapper]{flapper}} (i.e., \code{\link[flapper]{pf}}) can implement this approach to ensure that movement paths are biologically realistic.
#'
#' @author Edward Lavender
#' @export
#'

lcp_from_point <- function(origin,
                           surface,
                           destination = NULL,
                           cost = NULL, graph = NULL,
                           use_all_cores = FALSE,
                           verbose = TRUE){

  #### Get surface info
  cat_to_console <- function(..., show = verbose){
    if(show) cat(paste(..., "\n"))
  }
  t_onset <- Sys.time()
  cat_to_console(paste0("flapper::lcp_from_point() called (@ ", t_onset, ")..."))
  op <- options()
  options(scipen = 999)
  on.exit(options(op), add = TRUE)
  cells <- 1:raster::ncell(surface)
  surface_df <- data.frame(cell = cells, value = as.vector(surface))
  origin_cell <- raster::cellFromXY(surface, origin)
  destination_cell <- surface_df$cell
  if(!is.null(destination)){
    if(inherits(destination, "matrix")){
      destination_cell <- raster::cellFromXY(surface, destination)
    } else if(inherits(destination, c("numeric", "integer"))) {
      destination_cell <- destination
    } else if(inherits(destination, "function")){
      destination_cell <- destination_cell[destination(surface_df$value)]
    } else stop("class(destination) is un-supported. Only matrix (coordinates), integer (cell IDs) or functions are supported as inputs.")
  }

  #### Get cost surface (once), if unsupplied
  if(is.null(cost) & is.null(graph)) {
    cat_to_console("... Calculating cost surface...")
    costs <- lcp_costs(surface, verbose = verbose)
    tdist_mat <- costs$dist_total
  } else tdist_mat <- cost

  #### Define a dataframe of edges with associated costs
  if(is.null(graph)){
    cat_to_console("... Defining nodes, edges and costs to make graph...")
    non_zero_indices <- Matrix::which(tdist_mat != 0, arr.ind = TRUE)
    edges <- data.frame(from = cells[non_zero_indices[, 1]],
                        to = cells[non_zero_indices[, 2]],
                        cost = as.vector(tdist_mat[non_zero_indices]))

    #### Specify a dataframe of coordinates for each node in appropriate format
    # ... redefine origin/destination nodes in terms of cell numbers
    rxy   <- raster::coordinates(surface)
    coord <- data.frame(node = cells, X = rxy[, 1], Y = rxy[, 2])

    #### Make graph
    cat_to_console("... Constructing graph object...")
    makegraph_param <- list(df = edges, directed = FALSE, coords = coord)
    graph <- do.call(cppRouting::makegraph, makegraph_param)
  }

  #### Get LCP distances
  cat_to_console("... Calling cppRouting::get_distance_matrix() to get least-cost distances...")
  get_distance_param <- list(Graph = graph,
                             from = origin_cell,
                             to = destination_cell,
                             algorithm = "phast",
                             allcores = use_all_cores)
  dist_lcp <- do.call(cppRouting::get_distance_matrix, get_distance_param)
  dist_lcp <- as.numeric(dist_lcp)

  #### Make a raster
  cat_to_console("... Making raster of least-cost distances from origin ...")
  r <- surface
  r <- raster::setValues(r, NA)
  r[origin_cell] <- 0
  r[destination_cell] <- dist_lcp

  #### Return outputs
  t_end <- Sys.time()
  total_duration <- difftime(t_end, t_onset, units = "mins")
  cat_to_console(paste0("... flapper::lcp_from_point() call completed (@ ", t_end, ") after ~", round(total_duration, digits = 2), " minutes."))
  return(r)

}


######################################
######################################
#### lcp_over_surface()

#' @title Calculate shortest path(s) and/or distance(s) over a surface between origin and destination coordinates
#' @description This function computes the shortest path(s) and/or distance(s) over a \code{surface} between \code{origin} and \code{destination} coordinates. To implement this function, \code{origin} and \code{destination} coordinates need to be specified as matrices and the surface over which movement occurs should be supplied as a \code{\link[raster]{raster}}. Since determining shortest paths can be computationally and memory-intensive, the \code{surface} can be reduced in size and/or resolution before these are computed, by (a) cropping the surface within user-defined extents; (b) focusing on a buffer zone along a Euclidean transect connecting \code{origin} and \code{destination} coordinates; (c) aggregating the surface to reduce the resolution; and/or (d) masking out areas over which movement is impossible (e.g., land for marine animals). Then, the function computes distances between connected cells, given (a) the planar distances between connected cells and (b) their difference in elevation. These distances are taken as a measure of `cost'. For each pair of \code{origin} and \code{destination} coordinates, or for all combinations of coordinates, these distances are used to compute the least-cost (i.e., shortest) path and/or the distance of this path, using functions in the \code{\link[cppRouting]{cppRouting}} or \code{\link[gdistance]{gdistance}} package. The function returns the shortest path(s) and/or their distance(s) (m) along with a plot and a list of objects involved in the calculations.
#' @param origin A matrix which defines the coordinates (x, y) of the starting location(s). Coordinates should lie on a plane (i.e., Universal Transverse Mercator projection).
#' @param destination A matrix which defines the coordinates (x, y) of the finishing location(s). Coordinates should lie on a plane (i.e., Universal Transverse Mercator projection).
#' @param surface A \code{\link[raster]{raster}} over which the object (e.g., individual) must move from \code{origin} to \code{destination}. The \code{surface} must be planar (i.e., Universal Transverse Mercator projection) with units of metres in x, y and z directions (m). The \code{surface}'s \code{\link[raster]{resolution}} is taken to define the distance between horizontally and vertically connected cells and must be the same in both x and y directions (for \code{surface}'s with unequal horizontal resolution, \code{\link[raster]{resample}} can be used to equalise resolution: see Examples). Any cells with NA values (e.g., due to missing data) are treated as `impossible' to move though by the algorithm. In this case, the \code{surface} might need to be pre-processed so that NAs are replaced/removed before implementing the function, depending on their source.
#' @param crop (optional) An \code{\link[raster]{extent}} object that is used to \code{\link[raster]{crop}} the extent of the \code{surface}, before the least-cost algorithms are implemented. This may be useful for large rasters to reduce memory requirements/computation time.
#' @param buffer (optional) A named list of arguments, passed to \code{\link[rgeos]{gBuffer}} (e.g. \code{buffer = list(width = 1000)}) (m) that is used to define a buffer around a Euclidean transect connecting the \code{origin} and \code{destination}. (This option can only be implemented for a single \code{origin} and \code{destination} pair.) The \code{surface} is then cropped to the extent of this buffer before the least-cost algorithms are implemented. This may be useful for large rasters to reduce memory requirements and/or computation time.
#' @param aggregate (optional) A named list of arguments, passed to \code{\link[raster]{aggregate}}, to aggregate raster cells before the least-cost algorithms are implemented. This may be useful for large rasters to reduce memory requirements and/or computation time.
#' @param mask (optional) A Raster or Spatial \code{\link[raster]{mask}} that is used to prevent movement over `impossible' areas on the \code{surface}. This must also lie on a planar surface (i.e., Universal Transverse Mercator projection). For example, for marine animals, \code{mask} might be a \code{\link[sp]{SpatialPolygonsDataFrame}} which defines the coastline. The effect of the \code{mask} depends on \code{mask_inside} (see below).
#' @param mask_inside A logical input that defines whether or not to mask the \code{surface} inside (\code{TRUE}) or outside (\code{FALSE}) of the \code{mask} (see \code{\link[flapper]{mask_io}}).
#' @param plot A logical input that defines whether or not to plot the inputted and processed surfaces. If \code{TRUE}, the inputted and processed plots are produced side-by-side. For the inputted surface, the \code{mask} and the region selected (via \code{crop} and/or \code{buffer}) are shown along with the \code{origin} and \code{destination}. For the processed surface, the surface and the \code{origin} and \code{destination} are shown, along with the shortest path(s) (if and once computed: see \code{goal}). This is useful for checking that any \code{surface} processing steps have been applied correctly and the \code{origin} and \code{destination} are positioned correctly on the \code{surface}.
#' @param goal An integer that defines the output of the function: \code{goal = 1} computes shortest distances, \code{goal = 2} computes shortest paths and \code{goal = 3} computes both shortest paths and the corresponding distances. Note that \code{goal = 3} results in least-cost algorithms being implemented twice, which will be inefficient for large problems; in this case, use \code{goal = 2} to compute shortest paths and then calculate their distance using outputs returned by the function (see Value).
#' @param combination A character string (\code{"pair"} or \code{"matrix"}) that defines whether or not to compute shortest distances/paths for (a) each sequential \code{origin} and \code{destination} pair of coordinates (\code{combination = "pair"}) or (b) all combinations of \code{origin} and \code{destination} coordinates (\code{combination = "matrix"}). This argument is only applicable if there is more than one \code{origin} and \code{destination}. For \code{combination = "pair"}, the number of \code{origin} and \code{destination} coordinates needs to be the same, since each \code{origin} is matched with each \code{destination}.
#' @param method A character string (\code{"cppRouting"} or \code{"gdistance"}) that defines the method used to compute the shortest distances between the \code{origin} and the \code{destination}. \code{"cppRouting"} is the default \code{method}. Under this option, functions in the \code{\link[cppRouting]{cppRouting}} package are used to compute the shortest paths (\code{\link[cppRouting]{get_path_pair}} or \code{\link[cppRouting]{get_multi_paths}} for each pair of coordinates or for all combinations of coordinates, respectively) and/or distances (\code{\link[cppRouting]{get_distance_pair}} or \code{\link[cppRouting]{get_distance_matrix}}). This package implements functions written in C++ and massively outperforms the other \code{method = "gdistance"} for large problems. Otherwise, if \code{method = "gdistance"}, functions in the \code{\link[gdistance]{gdistance}} are called iteratively to compute shortest paths (via \code{\link[gdistance]{shortestPath}}) or distances (via \code{\link[gdistance]{costDistance}}).
#' @param cppRouting_algorithm A character string that defines the algorithm used to compute shortest paths or distances. This is only applicable if \code{method = "cppRouting"}: \code{method = "gdistance"} implements Dijkstra's algorithm only. For shortest paths or their distances between pairs of coordinates, the options are \code{"Dijkstra"}, \code{"bi"}, \code{"A*"} or \code{"NBA"} for the uni-directional Dijkstra, bi-directional Dijkstra, A star unidirectional search or new bi-directional A star algorithms respectively (see \code{\link[cppRouting]{get_path_pair}} or \code{\link[cppRouting]{get_distance_pair}}). For shortest paths between all combinations of coordinates, \code{cppRouting_algorithm} is ignored and the Dijkstra algorithm is implemented recursively. For shortest distances between all combinations of coordinates, the options are \code{"phast"} or \code{"mch"} (see \code{\link[cppRouting]{get_distance_matrix}}).
#' @param use_all_cores,cl,varlist Parallelisation options for \code{method = "cppRouting"} (\code{use_all_cores}) or \code{method = "gdistance"} (\code{cl} and \code{varlist}) respectively. If \code{method = "cppRouting"}, parallelisation is implemented via \code{use_all_cores} for computing shortest distances only (not computing shortest paths). \code{use_all_cores} is a logical input that defines whether or not to use all cores for computing shortest distance(s). If \code{method = "gdistance"}, parallelisation is implemented via \code{cl} and \code{varlist} for both shortest paths and distances function calls. \code{cl} is (a) a cluster object from \code{\link[parallel]{makeCluster}} or (b) an integer that defines the number of child processes. \code{varlist} is a character vector of variables for export (see \code{\link[flapper]{cl_export}}). Exported variables must be located in the global environment. If a cluster is supplied, the connection to the cluster is closed within the function (see \code{\link[flapper]{cl_stop}}). For further information, see \code{\link[flapper]{cl_lapply}} and \code{\link[flapper]{flapper-tips-parallel}}.
#' @param check A logical input that defines whether or not to check function inputs. If \code{TRUE}, internal checks are implemented to check user-inputs and whether or not inputted coordinates are in appropriate places on the processed \code{surface} (for instance, to ensure inputted coordinates do not lie over masked areas). This helps to prevent intractable error messages. If \code{FALSE}, these checks are not implemented, so function progress may be faster initially (especially for large \code{origin}/\code{destination} coordinate matrices).
#' @param verbose A logical input that defines whether or not to print messages to the console to monitor function progress. This is especially useful with a large \code{surface} since the algorithms are computationally intensive.
#'
#' @details
#' \subsection{Methods}{
#' This function was motivated by the need to determine the shortest paths and their distances between points for benthic animals, which must move over the seabed to navigate from A to B. For these animals, especially in areas with heterogeneous bathymetric landscapes and/or coastline, the shortest path that an individual must travel to move from A and B may differ substantially from the Euclidean path that is often used as a proxy for distance in biological studies. However, this function can still be used in situations where the surface over which an individual must move is irrelevant (e.g., for a pelagic animal), by supplying a flat surface; then shortest paths/distances simply depend on the planar distances between locations and any barriers (e.g., the coastline). (However, this process will be somewhat inefficient.)
#'
#' The function conceptualises an object moving across a landscape as a queen on a chessboard which can move, in eight directions around its current position, across this surface. Given the potentially large number of possible paths between an \code{origin} and \code{destination}, the surface may be reduced in extent or size before the game begins. To determine shortest path/distance over the surface between each \code{origin} and \code{destination} pair/combination, the function first considers the distance that an object must travel between pairs of connected cells. This depends on the planar distances between cells and their differences in elevation. Planar distances (\eqn{d_p}, m) depend on the movement type: under a rook's movement (i.e., horizontally or vertically), the distance (\eqn{d_{p,r}}) between connected cells is extracted from the raster's resolution (which is assumed to be identical in the x and y directions); under a bishop's movement (i.e., diagonally), the distance between connected cells \eqn{d_{p,b}} is given by Pythagoras' Theorem: \eqn{d_{p,b} = \sqrt{(d_{p, r}^2 + d_{p, r}^2)}}. Vertical distances (\eqn{d_v}, m) are simply the differences in height between cells. The total distance (\eqn{d_t}) between any two connected cells is a combination of these distances given by Pythagoras' Theorem: \eqn{d_t = \sqrt{(d_p^2 + d_v^2)}}. These distances are taken to define the `cost' of movement between connected cells. Thus, `costs' are symmetric (i.e., the cost of moving from A to B equals the cost of moving from B to A).
#'
#' This cost surface is then used to compute the shortest path and/or distance of the shortest path between each \code{origin} and \code{destination} pair/combination using functions in the \code{\link[cppRouting]{cppRouting}} or \code{\link[gdistance]{gdistance}} package. The functions implemented depend on the \code{goal} (i.e., whether the aim is to compute shortest paths, shortest distances or both) and, if there is more than one \code{origin}/\code{destination}, the \code{combination} type (i.e., whether to compute shortest paths/distances for each sequential pair of coordinates or all possible combinations of coordinates).
#' }
#'
#' \subsection{Warnings}{
#' The function returns a warning produced by \code{\link[gdistance]{transition}} which is implemented to facilitate the definition of the cost surface, before shortest paths/distances are computed by either method: `In .TfromR(x, transitionFunction, directions, symm) : transition function gives negative values'. This warning arises because the height differences between connecting cells can be negative. It can be safely ignored.
#' }
#'
#' @return
#' \subsection{A named list}{The function returns a named list. The most important element(s) of this list are `path_lcp' and/or `dist_lcp', the shortest path(s) and/or distance(s) (m) between \code{origin} and \code{destination} coordinate pairs/combinations. `path_lcp' is returned if \code{goal = 2} or \code{goal = 3} and `dist_lcp' is returned if \code{goal = 1} or \code{goal = 3}. `path_lcp' contains (a) a dataframe with the cells comprising each path (`cells'), (b) a named list containing a \code{\link[sp]{SpatialLines}} object for each path (`SpatialLines') and (c) a named list of matrices of the coordinates of each path (`coordinates'). `dist_lcp' is a (a) numeric vector or (b) matrix with the distances (m) between each pair or combination of coordinates respectively. If `dist_lcp' is computed, `dist_euclid', the Euclidean distances (m) between the \code{origin} and \code{destination}, is also returned for comparison.
#' }
#'
#' \subsection{Common elements}{Other elements of the list record important outputs at sequential stages of the algorithm's progression. These include the following elements: `args', a named list of user inputs; `time', a dataframe that defines the times of sequential stages in the algorithm's progression; `surface', the surface over which shortest distances are computed (this may differ from the inputted surface if any of the processing options, such as \code{crop}, have been implemented); `surface_param', a named list that defines the cell IDs, the number of rows, the number of columns, the coordinates of the implemented surface and the cell IDs of the \code{origin} and \code{destination} nodes; `cost', a named list of arguments that defines the distances (m) between connected cells under a rook's or bishop's movement (`dist_rook' and `dist_bishop'), the planar and vertical distances between connected cells (`dist_planar' and `dist_vertical') and the total distance between connected cells (`dist_total'); and `cppRouting_param' or `gdistance_param', a named list of arguments used to compute shortest paths/distances via \code{\link[cppRouting]{cppRouting}} or \code{\link[gdistance]{gdistance}} (see below).
#' }
#'
#' \subsection{Method-specific elements}{If \code{method = "cppRouting"}, the `cppRouting_param' list contains a named list of arguments passed to \code{\link[cppRouting]{makegraph}} (`makegraph_param') as well as \code{\link[cppRouting]{get_path_pair}} (`get_path_pair_param') or \code{\link[cppRouting]{get_multi_paths}} (`get_multi_paths_param') and/or \code{\link[cppRouting]{get_distance_pair}} (`get_distance_pair_param') or \code{\link[cppRouting]{get_distance_matrix}} (`get_distance_matrix_param'), depending on whether or not shortest paths and/or distances have been computed (see \code{goal}) and whether or not shortest paths/distances have been computed for each pair of coordinates or all combinations of coordinates. If \code{method = "gdistance"}, this list contains a named list of arguments passed iteratively, for each pair/combination of coordinates, to \code{\link[gdistance]{shortestPath}} (`shortestPath_param') or \code{\link[gdistance]{costDistance}} (`costDistance_param'). This includes an object of class TransitionLayer (see \code{\link[gdistance]{Transition-classes}}), in which the \code{transitionMatrix} slot contains a (sparse) matrix that defines the ease of moving between connected cells (the reciprocal of the `dist_total' matrix).
#' }
#'
#' \subsection{Plot}{If \code{plot = TRUE}, a plot is also produced of the inputted and processed surfaces that are used in the calculations, along with the shortest path(s) (if and once computed).
#' }
#'
#' @examples
#' #### Example types
#' # Shortest distances between a single origin and a single destination
#' # Shortest paths between a single origin and a single destination
#' # Shortest distances/paths between origin/destination pairs
#' # Shortest distances/paths between all origin/destination combinations
#'
#' #### Simulate a hypothetical landscape
#' # Define a miniature, blank landscape with appropriate dimensions
#' proj_utm <- sp::CRS(SRS_string = "EPSG:32629")
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
#' # In the same way, we can force the shortest path away from particular areas
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
#' #### Example (5C): Reduce the complexity of the landscape via aggregation
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
#' #### Shortest paths between a single origin and a single destination
#'
#' #### Example (8) Shortest paths (goal = 2) only using default method
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
#' # Note the path is also added to the plot produced, if plot = TRUE.
#'
#' #### Example (9) Shortest distances and paths (goal = 3) using default method
#' out1 <- lcp_over_surface(origin = rxy[1, , drop = FALSE],
#'                          destination = rxy[6, , drop = FALSE],
#'                          surface = r,
#'                          goal = 3)
#' out1$dist_lcp; out1$path_lcp
#'
#' #### Example (10) Shortest distances and paths (goal = 3) via gdistance
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
#' # ... paths for method = "gdistance"
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
#' #### Shortest distances/paths between origin/destination pairs
#'
#' #### Example (12): Shortest distances/paths computed in sequence:
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
#' #### Example (13): Shortest distances/paths computed in parallel:
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
#' #### Shortest distances/paths between all origin/destination combinations
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
           verbose = TRUE
  ){

    ########################################
    #### Set up raster

    #### Define a list of outputs
    cat_to_console <- function(..., show = verbose){
      if(show) cat(paste(..., "\n"))
    }
    t_onset <- Sys.time()
    cat_to_console(paste0("flapper::lcp_over_surface() called (@ ", t_onset, ")..."))
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
                     verbose = TRUE)

    #### Checks
    op <- options()
    options(scipen = 999)
    on.exit(options(op), add = TRUE)
    if(check){
      cat_to_console("... Checking user inputs...")
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
      cat_to_console("... Cropping raster to inputted extent...")
      surface <- raster::crop(surface, crop)
    }

    #### Crop raster around Euclidean transect(s) for speed
    if(!is.null(buffer)){
      cat_to_console("... Cropping raster to buffer zone around a Euclidean transect between the origin and destination...")
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
      cat_to_console("... Aggregating raster...")
      aggregate$x <- surface
      surface <- do.call(raster::aggregate, aggregate)
    }

    #### Mask raster
    if(!is.null(mask)){
      cat_to_console("... Masking raster...")
      surface <- mask_io(x = surface, mask = mask, mask_inside = mask_inside)
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
    cat_to_console("... Defining cost matrix...")

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
      cat_to_console("... Calculating Euclidean distance(s)...")
      if(combination == "pair") allpairs <- FALSE else if(combination == "matrix") allpairs <- TRUE
      dist_euclid <- raster::pointDistance(origin, destination, lonlat = FALSE, allpairs = allpairs)
      if(nrow(origin) == 1 & nrow(destination) == 1) dist_euclid <- as.numeric(dist_euclid)
      if(allpairs){
        rownames(dist_euclid) <- origin_cell
        colnames(dist_euclid) <- destination_cell
      }
      out$dist_euclid <- dist_euclid
    }

    #### LCP distances/paths
    if(method == "gdistance"){
      cat_to_console("... Using method = 'gdistance'...")

      #### Calculate the reciprocal of the costs
      # ... (i.e. the total ease of movement between connected cells):
      cat_to_console("... ... Defining 'ease' matrix...")
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
        cat_to_console("... ... Implementing Dijkstra's algorithm to compute least-cost distance(s)...")
        costDistance_param <- list(x = tr, fromCoords = origin, toCoords = destination)
        out$gdistance_param <- list(costDistance_param = costDistance_param)

        ## Least cost distances between pairs
        if(combination == "pair"){

          # Loop over each pair of coordinates and compute least-cost distances
          cl_export(cl, varlist)
          dist_lcp_by_pair <- pbapply::pblapply(coords_by_pair, cl = cl, FUN = function(xy){
            costDistance_param <- list(x = tr, fromCoords = xy[, 1:2, drop = FALSE], toCoords = xy[, 3:4, drop = FALSE])
            dist_lcp <- do.call(gdistance::costDistance, costDistance_param)
            return(dist_lcp)
          })
          if(all(goal %in% 1)) cl_stop(cl)
          # Process least-cost distances
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

      #### Calculate least-cost paths between points
      if(goal %in% c(2, 3)){
        cat_to_console("... ... Implementing Dijkstra's algorithm to compute least-cost path(s)...")
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

        ## Loop over each pair of coordinates and compute least-cost paths
        if(all(goal %in% 2)) cl_export(cl, varlist)
        path_lcp_SpatialLines <- pbapply::pblapply(coords_by_pair, cl = cl, FUN = function(xy){
          param <- list(x = tr, origin = xy[, 1:2, drop = FALSE], goal = xy[, 3:4, drop = FALSE], output = "SpatialLines")
          path_lcp <- do.call(gdistance::shortestPath, param)
          return(path_lcp)
        })
        cl_stop(cl)

        ## Process least-cost paths
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
      cat_to_console("... Using method = 'cppRouting'...")

      #### Define a dataframe of edges with associated costs
      cat_to_console("... ... Defining nodes, edges and costs to make graph...")
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
      cat_to_console("... ... Constructing graph object...")
      makegraph_param <- list(df = edges, directed = FALSE, coords = coord)
      out$cppRouting_param <- list(makegraph_param = makegraph_param)
      graph <- do.call(cppRouting::makegraph, makegraph_param)
      out$time <- rbind(out$time, data.frame(event = "graph_defined", time = Sys.time()))

      #### Compute all shortest distance between origin and destination nodes
      if(goal %in% c(1, 3)){
        # Define parameters to compute distances
        cat_to_console(paste("... ... Implementing",  cppRouting_algorithm, "algorithm to compute least-cost distance(s)..."))
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

      #### Compute shortest paths between origin and destination nodes
      if(goal %in% c(2, 3)){

        ## Set up to compute shortest paths
        # Define function to compute paths and add parameters to output list
        if(combination == "pair"){
          cat_to_console(paste("... ... Implementing",  cppRouting_algorithm, "algorithm to compute least-cost paths(s)..."))
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
          cat_to_console("... .... Implementing Dijkstra's algorithm recursively to compute least-cost paths(s)...")
          get_path_param <- list(Graph = graph,
                                 from = origin_cell,
                                 to = destination_cell,
                                 keep = NULL,
                                 long = TRUE)
          get_path <- cppRouting::get_multi_paths
          out$cppRouting_param$get_multi_paths_param <- get_path_param
        }

        ## Compute paths (origin, destination and cells)
        path_lcp_cells <- do.call(get_path, get_path_param)
        out$time <- rbind(out$time, data.frame(event = "path_lcp_defined", time = Sys.time()))

        ## Process paths dataframe
        # Define columns
        connections <- factor(paste0(origin_cell, "-", destination_cell))
        colnames(path_lcp_cells)   <- c("origin", "destination", "cell")
        path_lcp_cells$origin      <- as.integer(path_lcp_cells$origin)
        path_lcp_cells$destination <- as.integer(path_lcp_cells$destination)
        path_lcp_cells$cell        <- as.integer(path_lcp_cells$cell)
        path_lcp_cells$path <- paste0(path_lcp_cells$origin, "-", path_lcp_cells$destination)
        path_lcp_cells$path <- factor(path_lcp_cells$path, levels = unique(path_lcp_cells$path))
        if(any(!(connections %in% path_lcp_cells$path))){
          problematic_connections <- as.character(connections[!(connections %in% path_lcp_cells$path)])
          warning(paste0("LCPs could not be resolved for ",
                         length(problematic_connections),
                         " unique origin-destination pairs: ",
                         paste0(problematic_connections, collapse = ", "),
                         ". These connections are not included in function outputs."),
                  immediate. = TRUE, call. = FALSE)
        }
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

        # Define path coordinates and SpatialLines
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

    #### Add least-cost paths to output object
    if(goal %in% c(2, 3)){
      # paths
      names(path_lcp_SpatialLines) <- unique(path_lcp_cells$path_id)
      names(path_lcp_coordinates)  <- unique(path_lcp_cells$path_id)
      out$path_lcp <- list(cells = path_lcp_cells,
                           SpatialLines = path_lcp_SpatialLines,
                           coordinates = path_lcp_coordinates
      )
      # Add paths to plot, if requested
      if(plot){
        lapply(path_lcp_SpatialLines, function(l) raster::lines(l, col = "royalblue", lwd = 2))
        graphics::par(pp)
      }
    }

    #### Return outputs
    t_end <- Sys.time()
    total_duration <- difftime(t_end, t_onset, units = "mins")
    out$time <- rbind(out$time, data.frame(event = "finish", time = t_end))
    out$time$serial_duration <- Tools4ETS::serial_difference(out$time$time)
    cat_to_console(paste0("... flapper::lcp_over_surface() call completed (@ ", t_end, ") after ~", round(total_duration, digits = 2), " minutes."))
    return(out)

  }


######################################
######################################
#### lcp_interp()

#' @title Interpolate shortest (least-cost) paths between locations along a movement path
#' @description This function is a wrapper for \code{\link[flapper]{lcp_over_surface}} designed to interpolate shortest (least-cost) paths between sequential locations along an animal movement path. The function is specifically motivated for the interpolation of paths between locations for benthic animals, whose movement is restricted by the seabed (see \code{\link[flapper]{lcp_over_surface}}).
#' @param paths A dataframe that defines movement paths. This must contain a unique, integer identifier for each path from 1 to the number of paths (`path_id'); an integer time step index (`timestep') and the coordinates of sequential locations (`cell_x' and `cell_y'). The example dataset comprising movement paths reconstructed over the seabed by the depth-contour particle filtering algorithm (\code{\link[flapper]{dat_dcpf_paths}}) exemplifies of the appropriate format.
#' @param surface A \code{\link[raster]{raster}} of the surface over which the object (i.e., individual) moved between sequential locations. \code{\link[flapper]{lcp_over_surface}} forces some restrictions on the form of this raster. The \code{surface} must be planar (i.e., Universal Transverse Mercator projection) with units of metres and equal resolution in the x, y and z directions (see \code{\link[flapper]{lcp_over_surface}}).
#' @param ... Additional arguments, passed to \code{\link[flapper]{lcp_over_surface}}, to control the interpolation (excluding \code{origin}, \code{destination}, \code{combination} and \code{goal}).
#' @param keep_cols A logical input that defines whether or not to retain all columns in \code{paths} in the returned dataframe (see Value).
#' @param calc_distance A logical input that defines whether or not to calculate distances between sequential positions along the shortest paths (see Value).
#' @param verbose A logical input that defines whether or not to print messages to the console to monitor function progress.
#'
#' @details If \code{calc_distances = TRUE}, distances are calculated with movement from the previous location to the `current' location (see Value).
#'
#' A useful application of this function in \code{\link[flapper]{flapper}} is the post-hoc evaluation of particle filtering movement algorithms (see \code{\link[flapper]{pf}}). These can be implemented using movement models based on Euclidean or shortest distances. Since the former is typically much faster, a useful starting point is to implement the chosen algorithm using Euclidean distances and then, for the sample of paths reconstructed by the algorithm, use \code{\link[flapper]{lcp_interp}} to examine the similarity between Euclidean and shortest distances and the effects of updated distance values on movement probabilities. If a Euclidean distances implementation of an algorithm is acceptable (i.e., minimum swimming distances are not too large under the movement model), shortest distances from \code{\link[flapper]{lcp_interp}} can be used to adjust the movement probabilities to more realistic values. Alternatively, a shortest distances implementation of the movement path algorithm may be necessary.
#'
#' This function can also be useful for visualising movement paths (e.g., via \code{\link[flapper]{pf_plot_3d}}).
#'
#' @return The function returns a dataframe or a list depending on the input to \code{calc_distance}. If \code{calc_distance = FALSE}, the function returns a dataframe that defines, for each path (`path_id'), for each time step (`timestep'), a vector of the cell coordinates (`cell_x', `cell_y', `cell_z') on the \code{surface} that define the shortest path from the location at the previous time step to the location at the current time step. Thus, the first row contains the individual's initial location and subsequent rows for that time step (\code{timestep = 2} if the first time step is defined as \code{timestep =  1}, since movement is considered from the location at the `previous' time step to the location at the `current' time step) include the sequential locations of an individual, were it to have moved along the shortest path, to the location at the next time step.  If \code{keep_cols = TRUE}, coordinate columns are suffixed by `.x' and the coordinates for each time step (as inputted) are included (with the `.y' suffix) along with any other columns in the inputted dataframe (\code{paths}).
#'
#' If \code{calc_distance = TRUE}, a named list is returned with two elements. The `path_lcp' element contains the dataframe of interpolated path coordinates, as described above, with an extra `dist' column that defines the distance over the surface between sequential positions along the path. The `dist_lcp' element contains a dataframe, exactly as inputted via \code{paths}, but with the total distance along each shortest path (from the `previous' location to the `current' location) included as a `dist' column (this is simply the sum of the distances provided in \code{path_lcp$dist}) for all movements within each timestep for each path.
#'
#' @examples
#' #### Define movement paths
#' # We will interpolate LCPs between sequential locations
#' # ... of a skate on the seabed
#' # ... reconstructed by the DCPF algorithm (dc() & pf()), using the
#' # ... example dat_dcpf_* datasets. We extract the paths
#' # ... and the surface over which movement occurred from this object:
#' paths   <- dat_dcpf_paths
#' surface <- dat_dcpf_histories$args$bathy
#'
#' #### Example (1): Implement lcp_interp() for an example path
#' # ... with calc_distance = FALSE
#' # Implement approach
#' paths_1 <- paths[paths$path_id == 1, ]
#' paths_interp_1 <- lcp_interp(paths = paths_1,
#'                              surface = surface,
#'                              calc_distance = FALSE)
#' # With calc_distance = FALSE, we get a dataframe with sequential locations
#' utils::head(paths_interp_1)
#'
#' #### Example (2): Keep the original columns in 'paths'
#' # ... via keep_cols = TRUE
#' paths_interp_2 <- lcp_interp(paths = paths_1,
#'                              surface = surface,
#'                              calc_distance = FALSE,
#'                              keep_cols = TRUE)
#' utils::head(paths_interp_2)
#'
#' #### Example (3): Calculate shortest distances along each path
#' # ... via calc_distance = TRUE (the default)
#' paths_interp_3 <- lcp_interp(paths = paths_1,
#'                              surface = surface)
#' # A named list is returned
#' names(paths_interp_3)
#' # The 'dist_lcp' element is as inputted but with a 'dist' column that defines
#' # ... the total distance along the shortest path between sequential locations
#' utils::head(paths_interp_3$dist_lcp)
#' # The 'path_lcp' element contains the paths, as described above, with an
#' # ... additional 'dist' column for the distances between sequential locations
#' # ... along each shortest path (i.e., time step)
#' utils::head(paths_interp_3$path_lcp)
#'
#' @author Edward Lavender
#' @export

lcp_interp <- function(paths, surface, ..., keep_cols = FALSE, calc_distance = TRUE, verbose = TRUE){

  #### Set up
  cat_to_console <- function(..., show = verbose){
    if(show) cat(paste(..., "\n"))
  }
  t_onset <- Sys.time()
  cat_to_console(paste0("flapper::lcp_interp() called (@ ", t_onset, ")..."))
  cat_to_console("... Setting up function...")
  out <- list(paths = paths)
  check_names(input = paths, req = c("path_id", "timestep", "cell_x", "cell_y"))
  check...(c("origin", "destination", "combination", "goal"),...)

  #### Reformat dataframe with the 'previous' and 'current' step for each path
  cat_to_console("... Processing paths...")
  paths_xy <- lapply(split(paths, paths$path_id), function(d){
    d <- data.frame(d)
    d$cell_id_previous <- dplyr::lag(d$cell_id)
    d$cell_x_previous  <- dplyr::lag(d$cell_x)
    d$cell_y_previous  <- dplyr::lag(d$cell_y)
    d$cell_id_current  <- d$cell_id
    d$cell_x_current   <- d$cell_x
    d$cell_y_current   <- d$cell_y
    return(d)
  }) %>% dplyr::bind_rows()
  paths_xy <- paths_xy[which(!is.na(paths_xy$cell_x_previous)), ]

  #### Assign origin/destination pair IDs
  # Assign IDs based on cell IDs (rather than origin-destination coordinates) to avoid rounding issues
  # ... when deriving the corresponding IDs for reconstructed paths (see below)
  paths_xy$pair <- paste0(paths_xy$cell_id_previous, "-", paths_xy$cell_id_current)

  #### Get non-duplicate pairs for speed
  paths_xy$dup <- duplicated(paths_xy$pair)
  paths_xy_sbt <- paths_xy[which(!paths_xy$dup), ]

  #### Calculate least-cost paths between coordinates
  cat_to_console("... Calculating least-cost paths via flapper::lcp_over_surface()...")
  origin <- as.matrix(paths_xy_sbt[, c("cell_x_previous", "cell_y_previous")])
  destination <- as.matrix(paths_xy_sbt[, c("cell_x_current", "cell_y_current")])
  paths_lc <- lcp_over_surface(origin = origin,
                               destination = destination,
                               surface = surface,
                               combination = "pair",
                               goal = 2L,
                               verbose = verbose,...)
  # Extract xy coordinates for LCP between each origin and destination pairs
  paths_lc_xy <- paths_lc$path_lcp$coordinates
  names(paths_lc_xy) <- sapply(paths_lc_xy, function(xy) paste0(raster::cellFromXY(surface, xy[1, ]), "-",
                                                                raster::cellFromXY(surface, xy[nrow(xy), ])))
  # Reformat to add back duplicate paths
  paths_xy <- paths_xy[paths_xy$pair %in% names(paths_lc_xy), ]
  paths_lc_xy <- lapply(split(paths_xy, 1:nrow(paths_xy)), function(d){
    paths_lc_xy[d$pair]
  })

  #### Define a dataframe of LCP coordinates
  # Define skeleton dataframe with path IDs, time steps and cell coordinates
  paths_lc_xy <- lapply(1:length(paths_lc_xy), function(i){
    xy <- paths_lc_xy[[i]]
    xy <- data.frame(xy)
    colnames(xy) <- c("cell_x", "cell_y")
    xy$path_id <- paths_xy$path_id[i]
    xy$timestep <- paths_xy$timestep[i]
    return(xy)
  })
  paths_lc_xyz <- do.call(rbind, paths_lc_xy)
  # Add z coordinates
  paths_lc_xyz$cell_z <- raster::extract(surface, paths_lc_xyz[, c("cell_x", "cell_y")])
  # Calculate distances
  if(calc_distance){
    paths_lc_xyz <-
      paths_lc_xyz %>%
      dplyr::arrange(.data$path_id, .data$timestep) %>%
      dplyr::group_by(.data$path_id) %>%
      dplyr::mutate(dist = dist_btw_points_3d(dplyr::lag(.data$cell_x), .data$cell_x,
                                              dplyr::lag(.data$cell_y), .data$cell_y,
                                              dplyr::lag(.data$cell_z), .data$cell_z)
      ) %>%
      as.data.frame()
  }
  # Organise columns
  cols <- c("path_id", "timestep", "cell_x", "cell_y", "cell_z", "dist")
  if(!calc_distance) cols <- cols[!("dist" == cols)]
  paths_lc_xyz <- paths_lc_xyz[, cols]
  if(keep_cols){
    paths_lc_xyz <- dplyr::left_join(paths_lc_xyz, out$paths, by = c("path_id", "timestep"))
  }
  # Add to outputs
  out$path_lcp <- paths_lc_xyz

  #### Summarise distances
  if(calc_distance){
    cat_to_console("... Summarising distances for each least-cost path...")
    dists <-
      paths_lc_xyz %>%
      dplyr::group_by(.data$path_id, .data$timestep) %>%
      dplyr::summarise(dist = sum(.data$dist, na.rm = TRUE)) %>%
      as.data.frame()
    dists$key      <- paste0(dists$path_id, "-", dists$timestep)
    out$paths$key  <- paste0(out$paths$path_id, "-", out$paths$timestep)
    out$paths$dist <- dists$dist[match(out$paths$key, dists$key)]
    out$paths$key  <- NULL
    out$dist_lcp   <- out$paths
    out$paths <- NULL
  }

  #### Return outputs
  if(!calc_distance) out <- out$path_lcp
  t_end <- Sys.time()
  total_duration <- difftime(t_end, t_onset, units = "mins")
  cat_to_console(paste0("... flapper::lcp_interp() call completed (@ ", t_end, ") after ~", round(total_duration, digits = 2), " minutes."))
  return(out)

}


######################################
######################################
#### lcp_comp()

#' @title Compare Euclidean and shortest distances
#' @description This function compares Euclidean distances to shortest distances for randomly sampled pairs of points across a \code{\link[raster]{raster}}.
#'
#' @param surface A \code{\link[raster]{raster}} from which to sample locations and calculate distances. The \code{surface} must be planar with units of metres in x, y and z directions. The \code{\link[raster]{resolution}}  must be the same in both x and y directions (see \code{\link[flapper]{lcp_over_surface}}). If applicable, the \code{surface} should be masked by \code{barrier} (see \code{\link[flapper]{mask_io}}) prior to implementation of this function.
#' @param barrier (optional) A simple features geometry that defines a barrier (see \code{\link[flapper]{segments_cross_barrier}}). If supplied, the function determines whether or not each sampled pair of coordinates crosses the barrier.
#' @param distance A numeric vector of distances. If supplied, for each distance, \code{n} location pairs that are that approximately that distance apart (see \code{interval}) are sampled and, for each pair, Euclidean and shortest distances are calculated. If supplied, \code{mobility} should not be supplied (see below).
#' @param interval If \code{distance} is supplied, \code{interval} is a number that defines the range of distances around each distance that are admissible. For example, if \code{distance = 50} (m) and \code{interval = 5}, \code{n} location pairs that are within 45--55 m of each other are sampled.
#' @param mobility (optional) A number that defines the maximum distance between sampled locations. If supplied, \code{n} locations are sampled; the Euclidean distances between all combinations of locations are calculated; and for the location pairs for which the Euclidean distance is less than \code{mobility} (up to a maximum of \code{n_max} pairs), shortest distances are calculated. This method is potentially fast the \code{distance}-vector-based method that is implemented is \code{distance} is supplied, but the final number of locations cannot be predetermined. If supplied, \code{distance} should not be supplied (see above).
#' @param n An integer that defines (a) the number of location pairs for each \code{distance} or (b) the initial number of sampled locations if \code{mobility} is supplied.
#' @param n_max If \code{mobility} is supplied, \code{n_max} is an integer that defines the maximum number of location pairs for which shortest distances are calculated.
#' @param graph (optional) A graph object that defines cell nodes and edge costs for connected cells within the surface (see \code{\link[flapper]{lcp_graph_surface}}). This is used for shortest-distance calculations.
#' @param ... Additional arguments passed to \code{\link[cppRouting]{get_distance_pair}} for shortest-distance calculations (\code{algorithm}, \code{constant} and/or \code{allcores}).
#' @param cl,varlist (optional) Parallelisation options for distance-based location sampling (i.e., if \code{distance != NULL}). \code{cl} is (a) a cluster object from \code{\link[parallel]{makeCluster}} or (b) an integer that defines the number of child processes. \code{varlist} is a character vector of variables for export (see \code{\link[flapper]{cl_export}}). Exported variables must be located in the global environment. If a cluster is supplied, the connection to the cluster is closed within the function (see \code{\link[flapper]{cl_stop}}). For further information, see \code{\link[flapper]{cl_lapply}} and \code{\link[flapper]{flapper-tips-parallel}}.
#' @param verbose A logical input that defines whether or not to print messages to the console to monitor function progress.
#'
#' @details This function was motivated by the need to determine the extent to which Euclidean distances are a suitable approximation of shortest distances in movement models for benthic animals.
#'
#' To address this issue, this function samples multiple pairs of points on a \code{\link[raster]{raster}} (\code{surface}) and compares the Euclidean and shortest distances between sampled locations.
#'
#' Two sampling methods are implemented. The first method (implemented if \code{distance} is supplied), samples \code{n} pairs of points for each distance that are approximately that distance apart. (Sampling across a grid is approximate because grid resolution is finite. The approximation is controlled by the \code{interval} parameter.) The advantage of this method is that the number of sampled location pairs for which Euclidean and shortest distances are compared is predetermined (by \code{n}), but the method can be slow if \code{n} is large. The second method (implemented if \code{mobility} is supplied), simply samples \code{n} locations and calculates the Euclidean distances between all combinations of sampled locations. Any location pairs that are more than \code{mobility} apart are then dropped and for the remaining location pairs (up to a maximum of \code{n_max} randomly sampled pairs) shortest distances are calculated. A potential advantage of this method is speed, but the number of location pairs below \code{mobility} for which Euclidean and shortest distances can be compared is not predetermined.
#'
#' Regardless of sampling method, the expectation for the analysis is that at short Euclidean distances, Euclidean distances are likely to approximate shortest distances well (unless there is a barrier, such as the coastline, in the way). At longer Euclidean distances, shortest distances are likely to become much longer than Euclidean distances. For some movement frameworks, especially the particle filtering (\code{\link[flapper]{pf}}) and processing (\code{\link[flapper]{pf_simplify}}) routines in \code{\link[flapper]{flapper}}, the Euclidean distance at which the shortest distances start to exceed the distance that an animal could move (over some time interval of interest) is an important parameter (termed \code{calc_distance_limit} in \code{\link[flapper]{pf_simplify}}). For sampled locations that are a Euclidean distance less than this threshold apart, it is reasonable to assume that there exists a `valid' shortest path over the \code{surface}, given the animal's mobility, without having to calculate that path; for sampled locations that are a Euclidean distance more than this threshold apart, this assumption is not valid. For areas with a clear-cut threshold, for Euclidean-based sampling frameworks (see \code{\link[flapper]{pf}}), it may be reasonable to reduce the animal's mobility to the threshold value to account for this difference; in other cases, it may be necessary to compute shortest distances for potentially problematic locations. Either way, an understanding of the extent to which Euclidean distances effectively approximate shortest distances at spatial scales relevant to an animal can help to minimise the number of shortest-distance calculations that are required within a movement modelling framework (which are much more computationally demanding than Euclidean-based calculations).
#'
#' @return The function returns a dataframe with sampled location pairs and the Euclidean and shortest distances between them. This contains the following columns:
#' \itemize{
#'   \item index -- A number that defines the row.
#'   \item id -- If \code{distance} is specified, \code{id} is number that defines the sample (from \code{1:n} for each \code{distance}).
#'   \item cell_x0, cell_y0, cell_x1, cell_y1 -- Numbers that define the coordinates of sampled locations.
#'   \item cell_id0, cell_id1 -- Integers that define the cell IDs for sampled locations.
#'   \item dist_sim -- If \code{distance} is specified, \code{dist_sim} is a number that defines the specified distance.
#'   \item dist_euclid -- A number that defines the Euclidean distance from \code{(cell_x0, cell_y0)} to \code{(cell_x1, cell_y1)}.
#'   \item dist_lcp -- A number that defines the shortest distance from \code{(cell_x0, cell_y0)} to \code{(cell_x1, cell_y1)}.
#'   \item barrier -- If applicable, \code{barrier} is factor that defines whether (\code{"1"}) or not (\code{"0"}) the Euclidean line connecting \code{(cell_x0, cell_y0)} and \code{(cell_x1, cell_y1)} crosses the \code{barrier}.
#' }
#'
#' @examples
#' #### Load packages
#' require(prettyGraphics)
#'
#' #### Define surface for examples
#' # We will focus on a relatively small area for speed
#' # The raster resolution should be equal in x and y directions
#' bathy      <- flapper::dat_gebco
#' boundaries <- raster::extent(700652.2, 708401.2, 6262905, 6270179)
#' blank      <- raster::raster(boundaries, res = c(5, 5))
#' bathy      <- raster::resample(bathy, blank)
#' bathy      <- mask_io(bathy, flapper::dat_coast, mask_inside = TRUE)
#' # raster::plot(bathy)
#' # raster::lines(dat_coast)
#'
#' #### Define barrier to movement
#' coastline <- raster::crop(dat_coast, raster::extent(bathy))
#' coastline <- sf::st_as_sf(coastline)
#'
#' #### Define graph for shortest-distance calculations
#' bathy_costs <- lcp_costs(bathy)                                 # ~0.38 mins
#' bathy_graph <- lcp_graph_surface(bathy, bathy_costs$dist_total) # ~0.26 mins
#'
#' #### Example (1): Distance-vector-based implementation
#'
#' ## Implement function
#' out_1 <- lcp_comp(surface = bathy,
#'                   distance = c(50, 100),
#'                   graph = bathy_graph,
#'                   barrier = coastline)
#'
#' ## Examine outputs
#' # The function returns a dataframe
#' utils::str(out_1)
#'
#' ## Visualise sampled locations on map
#' raster::plot(bathy)
#' raster::lines(dat_coast)
#' invisible(lapply(split(out_1, seq_len(nrow(out_1))), function(d){
#'   arrows(x0 = d$cell_x0, y0 = d$cell_y0,
#'          x1 = d$cell_x1, y1 = d$cell_y1, length = 0.01, lwd = 2)
#' }))
#'
#' ## Compare simulated, Euclidean and shortest distances
#' pp <- graphics::par(mfrow = c(1, 2))
#' pretty_plot(out_1$dist_sim, out_1$dist_euclid,
#'             xlab = "Distance (simulated) [m]",
#'             ylab = "Distance (Euclidean)[m]")
#' graphics::abline(0, 1)
#' pretty_plot(out_1$dist_euclid, out_1$dist_lcp,
#'             xlab = "Distance (Euclidean) [m]",
#'             ylab = "Distance (shortest) [m]")
#' graphics::abline(0, 1)
#' graphics::par(pp)
#'
#' #### Example (2): Distance-vector-based implementation in parallel
#' out_2 <- lcp_comp(surface = bathy,
#'                   distance = c(50, 100),
#'                   graph = bathy_graph,
#'                   barrier = coastline,
#'                   cl = parallel::makeCluster(2L)
#'                   )
#'
#' #### Example (3): Mobility-based implementation
#'
#' ## Implement function for mobility = 500
#' mob <- 500
#' out_3 <- lcp_comp(surface = bathy,
#'                   mobility = mob,
#'                   n = 1000,
#'                   graph = bathy_graph,
#'                   barrier = coastline
#'                   )
#'
#' ## Visualise sampled locations on map
#' raster::plot(bathy)
#' raster::lines(dat_coast)
#' invisible(lapply(split(out_3, seq_len(nrow(out_3))), function(d){
#'   arrows(x0 = d$cell_x0, y0 = d$cell_y0,
#'          x1 = d$cell_x1, y1 = d$cell_y1, length = 0.01, lwd = 0.5)
#' }))
#'
#' ## Compare Euclidean and shortest distances
#' # Extract the data for the paths that do not versus do cross a barrier
#' out_3_barrier0 <- out_3[out_3$barrier == 0, ]
#' out_3_barrier1 <- out_3[out_3$barrier == 1, ]
#' # Set up plotting window
#' pp <- graphics::par(mfrow = c(1, 2), oma = c(3, 3, 3, 3), mar = c(4, 4, 4, 4))
#' # Results for paths that do not cross a barrier
#' # ... Visualisation
#' pretty_plot(out_3_barrier0$dist_euclid, out_3_barrier0$dist_lcp,
#'             xlab = "Distance (Euclidean) [m]",
#'             ylab = "Distance (shortest) [m]",
#'             pch = ".")
#' graphics::abline(0, 1, col = "red")
#' graphics::abline(h = mob, col = "royalblue", lty = 3)
#' # ... Euclidean distance parameter at which mobility is exceeded
#' limit0 <- min(out_3_barrier0$dist_euclid[out_3_barrier0$dist_lcp > mob]); limit0
#' graphics::abline(v = limit0, col = "royalblue", lty = 3)
#' # Results for paths that cross a barrier
#' # ... Visualisation
#' pretty_plot(out_3_barrier1$dist_euclid, out_3_barrier1$dist_lcp,
#'            xlab = "Distance (Euclidean) [m]",
#'            ylab = "Distance (shortest) [m]",
#'            pch = ".")
#' graphics::abline(0, 1, col = "red")
#' graphics::abline(h = mob, col = "royalblue", lty = 3)
#' # ... Euclidean distance parameter at which mobility is exceeded
#' limit1 <-
#'   min(out_3_barrier1$dist_euclid[out_3_barrier1$dist_lcp > mob]); limit1
#' graphics::abline(v = limit1, col = "royalblue", lty = 3)
#' graphics::par(pp)
#'
#' ## Plot the difference between Euclidean and shortest distances with distance
#' pp <- graphics::par(mfrow = c(1, 2), oma = c(3, 3, 3, 3), mar = c(4, 4, 4, 4))
#' pretty_plot(out_3_barrier0$dist_euclid,
#'             out_3_barrier0$dist_lcp - out_3_barrier0$dist_euclid,
#'             xlab = "Distance (Euclidean) [m]",
#'             ylab = "Distance (shortest - Euclidean) [m]",
#'             pch = ".")
#' pretty_plot(out_3_barrier1$dist_euclid,
#'             out_3_barrier1$dist_lcp - out_3_barrier1$dist_euclid,
#'             xlab = "Distance (Euclidean) [m]",
#'             ylab = "Distance (shortest - Euclidean) [m]",
#'             pch = ".")
#' graphics::par(pp)
#'
#' ## Approximate shortest distances using a model based on Euclidean distances
#' # Fit model
#' mod <- lm(dist_lcp ~ 0 + dist_euclid:barrier, data = out_3)
#' # Examine model predictions
#' x <- seq(0, 500)
#' nd <- data.frame(barrier = factor(0, levels = c(0, 1)), dist_euclid = x)
#' ci <- list_CIs(predict(mod, nd, se.fit = TRUE))
#' pretty_plot(out_3$dist_euclid, out_3$dist_lcp,
#'                             col = c("grey", "darkred")[out_3$barrier],
#'                             cex = c(0.1, 1.5)[out_3$barrier])
#' add_error_envelope(x, ci)
#' nd$barrier <- factor(1, levels = c(0, 1))
#' ci <- list_CIs(predict(mod, nd, se.fit = TRUE))
#' add_error_envelope(x, ci,
#'                    add_fit = list(col = "red"),
#'                    add_ci = list(col = scales::alpha("red", 0.5),
#'                                  border = FALSE))
#' # Define predictive model based on this analysis
#' lcp_predict <- function(x, m, barrier){
#'   stopifnot(length(m) == 2L)
#'   stopifnot(all(barrier %in% c(0, 1)))
#'   return(m[1] * x * (barrier == 0) + m[2] * x * (barrier == 1))
#'   }
#'
#' @author Edward Lavender
#' @export

lcp_comp <- function(surface,
                     barrier = NULL,
                     distance = NULL,
                     interval = raster::res(surface)[1],
                     mobility = NULL,
                     n = 10L, n_max = n,
                     graph = NULL,...,
                     cl = NULL, varlist = NULL, verbose = TRUE){

  #### Initiate function
  t_onset <- Sys.time()
  cat_to_console <- function(..., show = verbose){
    if(show) cat(paste(..., "\n"))
  }
  cat_to_console(paste0("flapper::lcp_comp() called (@ ", t_onset, ")..."))

  #### Define parameters
  cat_to_console("... Defining parameters...")
  ## Check surface and barrier properties
  if(!all.equal(raster::res(surface)[1], raster::res(surface)[2])){
    stop("The resolution of the surface should be equal in the x and y directions. Consider calling raster::resample() before implementing this function.")
  }
  if(!is.null(barrier)){
    crs_surface <- raster::crs(surface)
    crs_barrier <- raster::crs(barrier)
    crs_check   <- all.equal(crs_surface, crs_barrier)
    if(isFALSE(crs_check)){
      warning("The CRS of 'surface' and 'barrier' is not identical.",
              immediate. = TRUE, call. = FALSE)
      message("... details: ", crs_check, ".")
      message("... surface CRS: '", crs_surface, "'.")
      message("... barrier CRS: '", crs_barrier, "'.")
    }
  }
  ## Check distances or mobility
  if(is.null(distance) & is.null(mobility))
    stop("Either 'distance' or 'mobility' must be supplied.", call. = FALSE)
  if(!is.null(distance) & !is.null(mobility))
    stop("Only one of 'distance' or 'mobility' must be supplied.", call. = FALSE)
  ## Check cluster
  if(!is.null(cl) & !is.null(mobility)){
    warning("'cl' is not implemented if 'mobility' is supplied.",
            immediate. = TRUE, call. = FALSE)
    cl <- varlist <- NULL
  }
  ## Check dots
  dots <- list(...)
  if(length(dots) > 0L) {
    dots_bool <- names(dots) %in% c("algorithm", "constant", "allcores")
    if(!(all(dots_bool))){
      issues <- names(dots)[!dots_bool]
      stop(paste0("The following argument(s) passed via ... are not supported: ",
                  paste(issues, collapse = ", "),
                  ". The only supported arguments are 'algorithm', 'constant' and 'allcores' (see cppRouting::get_distance_pair())."),
           call. = FALSE)
    }
  }

  #### Define a list of location (start and end) samples, for each distance
  if(!is.null(distance)){

    cat_to_console("... Sampling n location pairs for each distance value...")
    points_by_dist <-
      cl_lapply(distance, cl = cl, varlist = varlist, fun = function(step){

        #### Sample n random locations within area
        cat_to_console(paste0("... ... Sampling for distance = ", step, "..."))
        cat_to_console("... ... ... Sampling starting locations...")
        # step <- distance[1]
        start_xy <- raster::sampleRandom(surface, n, xy = TRUE, na.rm = TRUE)[, 1:2]
        colnames(start_xy) <- c("cell_x0", "cell_y0")

        #### Define n ending locations
        cat_to_console("... ... ... Sampling ending locations...")
        points_for_dist_by_row <-
          lapply(split(start_xy, seq_len(nrow(start_xy))), function(start_xyi){
            start_xyi <- matrix(start_xyi, ncol = 2)
            colnames(start_xyi) <- c("cell_x0", "cell_y0")
            around <- raster::distanceFromPoints(surface, start_xyi)
            around <- raster::mask(around, surface)
            around_mask <- around <= (step - interval) | around >= (step + interval)
            around  <- raster::mask(around, around_mask, maskvalue = 1)
            end_xyi <- raster::sampleRandom(around, 1, xy = TRUE, na.rm = TRUE)[, 1:2, drop = FALSE]
            colnames(end_xyi) <- c("cell_x1", "cell_y1")
            if(nrow(end_xyi) == 0L) {
              warning("There are no non-NA cells within ",
                      step - interval, " and ", step + interval,
                      " m of the current sampled location (",
                      paste0(start_xyi, collapse = ", "), ").",
                      immediate. = TRUE, call. = FALSE)
              points_for_dist_for_row <- data.frame(dist_sim = numeric(),
                                                    x0 = numeric(), y0 = numeric(),
                                                    x1 = numeric(), y1 = numeric())
            } else {
              points_for_dist_for_row <- data.frame(dist_sim = step, start_xyi, end_xyi)
            }
            return(points_for_dist_for_row)
          })

        #### Define single dataframe with point samples for specified distance
        points_for_dist <- dplyr::bind_rows(points_for_dist_by_row)
        rownames(points_for_dist) <- NULL
        return(points_for_dist)
      })

    #### Join dataframes
    point_samples <- dplyr::bind_rows(points_by_dist)

  } else {
    cat_to_console("... Sampling n locations...")
    xy  <- raster::sampleRandom(surface, n, xy = TRUE, na.rm = TRUE)[, 1:2]
    cat_to_console("... Calculating Euclidean distances between sampled locations...")
    mat <- .planedist2(xy, xy)
    ind <- which(upper.tri(mat, diag = FALSE), arr.ind = TRUE)
    cat_to_console("... Filtering by mobility...")
    point_samples <-
      data.frame(cell_x0 = xy[ind[, 1], 1],
                 cell_y0 = xy[ind[, 1], 2],
                 cell_x1 = xy[ind[, 2], 1],
                 cell_y1 = xy[ind[, 2], 2],
                 dist_euclid = mat[ind]) %>%
      dplyr::filter(.data$dist_euclid <= mobility)
    if(nrow(point_samples) == 0L)
      stop("'No samples within 'mobility' drawn.", call. = FALSE)
    if(nrow(point_samples) > n)
      point_samples <- point_samples %>% dplyr::slice_sample(n = n_max)
  }
  #### Define a dataframe for Euclidean/LCP distances
  cat_to_console("... Processing samples...")
  # Get raster cells
  point_samples$cell_id0  <-
    raster::cellFromXY(surface, point_samples[, c("cell_x0", "cell_y0")])
  point_samples$cell_id1  <-
    raster::cellFromXY(surface, point_samples[, c("cell_x1", "cell_y1")])
  # Get Euclidean distances on the grid
  if(!rlang::has_name(point_samples, "dist_euclid")){
    point_samples$dist_euclid <- raster::pointDistance(point_samples[, c("cell_x0", "cell_y0")],
                                                       point_samples[, c("cell_x1", "cell_y1")],
                                                       lonlat = FALSE)
  }

  #### Get segment overlaps with the coastline
  if(!is.null(barrier)){
    cat_to_console("... Determining barrier overlaps...")
    point_samples$barrier <-
      segments_cross_barrier(start = point_samples[, c("cell_x0", "cell_y0"), drop = FALSE],
                             end = point_samples[, c("cell_x1", "cell_y1"), drop = FALSE],
                             barrier = barrier)
    point_samples$barrier <- factor(as.character(point_samples$barrier + 0), levels = c("0", "1"))
  }

  #### Get LCPs between coordinate pairs
  cat_to_console("... Implementing shortest-distance calculations...")
  # Define graph if unsupplied
  if(is.null(graph)){
    cat_to_console("... ... Defining surface graph...")
    cat_to_console("... ... ... Defining cost surface via flapper::lcp_costs()...")
    costs <- lcp_costs(surface)
    cat_to_console("... ... ... Assembling graph via flapper::lcp_graph_surface()...")
    graph <- lcp_graph_surface(surface = surface, cost = costs$dist_total)
  }
  # Define LCPs
  cat_to_console(paste0("... ... Calculating shortest-distances via cppRouting::get_distance_pair() for ",
                        nrow(point_samples), " location pair(s) ..."))
  dist_lcp <-
    tryCatch(
      cppRouting::get_distance_pair(Graph = graph,
                                    from = point_samples$cell_id0,
                                    to = point_samples$cell_id1,...),
      error = function(e) return(e)
      )
  if(inherits(dist_lcp, "error")){
    warning("Error during shortest-distance calculations:",
            immediate. = TRUE, call. = FALSE)
    message(dist_lcp)
    message("Returning the dataframe for diagnostics and closing the function...")
    return(point_samples)
  } else point_samples$dist_lcp <- dist_lcp

  #### Tidy dataframe
  cat_to_console("... Cleaning results...")
  point_samples$index <- seq_len(nrow(point_samples))
  if(!is.null(distance))
    point_samples <-
    point_samples %>%
    dplyr::group_by(.data$dist_sim) %>%
    dplyr::mutate(id = 1:dplyr::n()) %>%
    data.frame()
  cols <- c("index", "id",
            "cell_x0", "cell_y0", "cell_x1", "cell_y1",
            "cell_id0", "cell_id1",
            "dist_sim",
            "dist_euclid",
            "dist_lcp",
            "barrier")
  point_samples <- point_samples[, cols[cols %in% colnames(point_samples)]]
  rownames(point_samples) <- NULL

  #### Return outputs
  t_end <- Sys.time()
  total_duration <- difftime(t_end, t_onset, units = "mins")
  cat_to_console(paste0("... flapper::lcp_comp() call completed (@ ", t_end, ") after ~", round(total_duration, digits = 2), " minutes."))
  return(point_samples)

}
