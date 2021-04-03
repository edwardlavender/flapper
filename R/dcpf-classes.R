######################################
######################################
#### .dcpf-class

#' @title ".dcpf" class
#' @description An S3 class that defines the object returned by \code{\link[flapper]{dcpf}}.
#'
#' @return This object is a named list, with two elements, which records the parameters used in the call to \code{\link[flapper]{dcpf}} and the particles sampled at each time step:
#'
#' \describe{
#'   \item{args}{A named list that records the function arguments used to generate outputs. This is as inputted, but with the `cells_by_time' and `calc_distance_graph' elements added if unsupplied and applicable.}
#'   \item{history}{A list of dataframes, one for each time step, that record particle samples. Each dataframe comprises \code{n} rows (one for each particle) and the following three columns:}
#'   \itemize{
#'     \item{id_current}{ An integer that uniquely defines each location on the bathymetry \code{\link[raster]{raster}} (\code{bathy}).}
#'     \item{pr_current} A double that defines the probability of movement into each cell.
#'     \item{timestep}{ An integer that defines each time step.}
#'   }
#' }
#'
#' @seealso \code{\link[flapper]{dat_dcpf_histories}} provides an example of a \code{\link[flapper]{.dcpf-class}} object. \code{\link[flapper]{dcpf_simplify}} converts \code{\link[flapper]{.dcpf-class}} objects into \code{\link[flapper]{dcpf-class}} objects that comprise the reconstructed movement paths.
#'
#' @author Edward Lavender
#' @docType package
#' @name .dcpf-class
NULL


######################################
######################################
#### dcpf-class

#' @title "dcpf" class
#' @description An S3 class that defines the object returned by \code{\link[flapper]{dcpf_simplify}}, following \code{\link[flapper]{dcpf}}.
#'
#' @return A dataframe with that records the reconstructed paths. This includes a unique identifier for each path, the time step, the location (cell ID and three-dimensional coordinates) on \code{bathy} and the probability associated with that cell, given movement from the previous cell, in the following columns:
#' \describe{
#'     \item{path_id}{ An integer that uniquely defines each path.}
#'     \item{timestep}{ An integer that defines each time step.}
#'     \item{cell_id}{ An integer that defines the cells ID of the surface raster (over which paths were reconstructed) retained by the algorithm.}
#'     \item{cell_x}{ A double that defines the cell x coordinate.}
#'     \item{cell_y}{ A double that defines the cell y coordinate.}
#'     \item{cell_z}{ A double that defines the value of the surface raster in each cell.}
#'     \item{cell_pr}{ A double that defines the probability of movement into each cell, given immediately preceding sampled locations (or the origin).}
#' Rows are ordered by path and then time step.
#'   }
#'
#' @seealso \code{\link[flapper]{dat_dcpf_paths}} are derived from \code{\link[flapper]{.dcpf-class}} objects via \code{\link[flapper]{dcpf_simplify}}. \code{\link[flapper]{dcpf_simplify}}provides an example.
#'
#' @author Edward Lavender
#' @docType package
#' @name dcpf-class
NULL
