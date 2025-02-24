% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spatial_tools.R
\name{cells_from_val}
\alias{cells_from_val}
\title{Obtain a RasterLayer or the cells of RasterLayer that are equal to or lie within a range of specified values}
\usage{
cells_from_val(x, y, interval = 1L, cells = TRUE, na.rm = TRUE)
}
\arguments{
\item{x}{A \code{\link[raster]{raster}} object.}

\item{y}{A number or a vector with two numbers representing the lower and upper boundaries of the interval of values within which cells are identified.}

\item{interval}{If y is a vector of length two, \code{interval} is an integer that controls whether or not to query cells within and equal to (\code{interval = 1L}) or simply within (\code{interval = 2L}) the range specified by \code{y}.}

\item{cells}{A logical variable that defines whether or not to return a vector of cell numbers (\code{TRUE}) or a RasterLayer of the cells corresponding to \code{y}.}

\item{na.rm}{A logical variable that defines whether or not to ignore NAs.}
}
\value{
The function returns a RasterLayer (if \code{cells = FALSE}) or an integer vector of numbers (if \code{cells = TRUE}) that defines the cells that are equal to, or lie within, specified value(s) \code{y}.
}
\description{
This function obtains a RasterLayer or the cells of a RasterLayer that are equal to a specified value or lie within a specified range of values. To implement this function, a \code{\link[raster]{raster}} (\code{x}) and a value or range of values (\code{y}) for which a RasterLayer or the numbers of cells corresponding to those values are desired must be supplied. If a range of values is supplied, an additional argument (\code{interval}) controls whether or not cells within and equal to or simply within the specified range are returned.
}
\examples{
# Define an example RasterLayer
ncl <- 10
nrw <- 10
n <- ncl * nrw
mat <- matrix(1:n, ncol = ncl, nrow = nrw, byrow = TRUE)
r <- raster::raster(mat)
# Visualise example RasterLayer
raster::plot(r)
raster::text(r)
# Obtain the number(s) of cells corresponding to a particular value
cells_from_val(r, 1)
# Obtain a RasterLayer of the cells corresponding to a particular value
r1 <- cells_from_val(r, 1, cells = FALSE)
raster::plot(r1)
# Obtain the number(s) of cells within or equal to a range of values
cells_from_val(r, c(1, 10))
# Obtain the number(s) of cells within a range of values
cells_from_val(r, c(1, 10), interval = 2L)
# As above but returning a raster layer
cells_from_val(r, c(1, 10), interval = 2L, cells = FALSE)
}
\author{
Edward Lavender
}
