% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spatial_tools.R
\name{xy_from_click}
\alias{xy_from_click}
\title{Get location coordinates from mouse click(s)}
\usage{
xy_from_click()
}
\value{
The function returns a two-column matrix with coordinates.
}
\description{
This function defines a two-column matrix of x, y coordinates from clicked locations on a map.
}
\examples{
\dontrun{
raster::plot(dat_gebco)
xy <- xy_from_click()
graphics::points(xy, col = "red")
}
}
\author{
Edward Lavender
}
