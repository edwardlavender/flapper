% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spatial_tools.R
\name{get_intersection}
\alias{get_intersection}
\title{Spatial intersections with one or more geometries}
\usage{
get_intersection(x, ...)
}
\arguments{
\item{x}{A list of Spatial* geometries supported by \code{\link[rgeos]{gIntersection}}. If \code{x} is a list with a single Spatial* element, that element is simply returned. If \code{x} is a list with two Spatial* elements, \code{\link[rgeos]{gIntersection}} is implemented directly to obtain the intersection. If \code{x} is a list with more than two Spatial* elements, \code{\link[rgeos]{gIntersection}} is implemented iteratively, effectively intersecting all geometries. \code{NULL} elements are ignored.}

\item{...}{Arguments passed to \code{\link[rgeos]{gIntersection}}.}
}
\value{
The function returns a Spatial* object of the intersection between all inputted geometries.
}
\description{
This function is a wrapper for \code{\link[rgeos]{gIntersection}} that extracts the intersection between inputted geometries. Unlike \code{\link[rgeos]{gIntersection}}, the function requires geometries to be specified as a list. The intersections between all (one or more) specified geometries, excluding \code{NULL} elements, are returned.
}
\examples{
## Define raster parameters (for sampling/plotting)
set.seed(1)
n <- raster::ncell(dat_gebco)
xy <- raster::coordinates(dat_gebco)
utm <- raster::crs(dat_gebco)
cols <- c("black", "red", "orange", "green", lwd = 2)
col_int <- scales::alpha("dimgrey", 0.5)

## Define example polygons
n_poly <- 4
polys <- lapply(1:n_poly, function(i) {
  xy_i <- xy[sample(1:n, 1), , drop = FALSE]
  xy_i <- sp::SpatialPoints(xy_i, utm)
  xy_i_buf <- rgeos::gBuffer(xy_i, width = 10000)
  return(xy_i_buf)
})

#### Example (1): Define intersections between multiple (here, four) polygons
# Plot area with polygons
raster::plot(dat_gebco)
lapply(1:length(polys), function(i) raster::lines(polys[[i]], col = cols[i]))
# Get intersection and add
int <- get_intersection(polys)
raster::plot(int, add = TRUE, col = col_int)

#### Example (2): Define intersections between fewer polygons
raster::plot(dat_gebco)
lapply(1:length(polys), function(i) raster::lines(polys[[i]], col = cols[i]))
raster::plot(get_intersection(polys[1:3]), add = TRUE, col = col_int)
raster::plot(get_intersection(list(polys[[1]], polys[[3]])), add = TRUE, col = col_int)

#### Example (3): The function can handle a single element
raster::plot(dat_gebco)
lapply(1:length(polys), function(i) raster::lines(polys[[i]], col = cols[i]))
raster::plot(get_intersection(polys[1]), add = TRUE, col = col_int)

#### Example (4): The function can handle NULL elements
raster::plot(dat_gebco)
lapply(1:length(polys), function(i) raster::lines(polys[[i]], col = cols[i]))
raster::plot(get_intersection(list(polys[[1]], NULL, NULL, polys[[4]])),
  add = TRUE, col = col_int
)

}
\seealso{
\code{\link[rgeos]{gIntersection}}
}
\author{
Edward Lavender
}
