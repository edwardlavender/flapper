######################################
######################################
#### pf_archive-class

#' @title "pf_archive" class
#' @description An S3 class that defines the object returned by \code{\link[flapper]{pf}}.
#'
#' @return This object is a named list, with two elements, which records the parameters used in the call to \code{\link[flapper]{pf}} and the particles sampled at each time step:
#'
#' \describe{
#'   \item{args}{A named list that records the function arguments used to generate outputs. This is as inputted, but with the `calc_distance_graph' elements added if unsupplied and applicable.}
#'   \item{history}{A list of dataframes, one for each time step, that record particle samples. Each dataframe comprises \code{n} rows (one for each particle) and the following five columns:}
#'   \itemize{
#'     \item{id_previous}{ An integer that uniquely defines a previous location on the \code{record} \code{\link[raster]{raster}}s. This is absent for the first time step. For subsequent time steps, this is \code{NA} if the fast Euclidean distances method has been used, which does not `remember' previous locations, but specified otherwise.)}
#'     \item{pr_current}{ A double that defines the movement probabilities associated with previous locations. This is absent for the first time step. For subsequent time steps, this is \code{NA} if the fast Euclidean distances method has been used, which does not `remember' previous locations, but specified otherwise.}
#'     \item{id_current}{ An integer that uniquely defines each location on the \code{record} \code{\link[raster]{raster}}s.}
#'     \item{pr_current} A double that defines the probability of movement into each cell.
#'     \item{timestep}{ An integer that defines each time step.}
#'   }
#' }
#'
#' @seealso \code{\link[flapper]{dat_dcpf_histories}} provides an example of a \code{\link[flapper]{pf_archive-class}} object. \code{\link[flapper]{pf_simplify}} converts \code{\link[flapper]{pf_archive-class}} objects into \code{\link[flapper]{pf_path-class}} objects that comprise the reconstructed movement paths.
#'
#' @author Edward Lavender
#' @docType package
#' @name pf_archive-class
NULL


######################################
######################################
#### pf_path-class

#' @title "pf" class
#' @description An S3 class that defines the object returned by \code{\link[flapper]{pf_simplify}}, following \code{\link[flapper]{pf}}.
#'
#' @return A dataframe with that records the reconstructed paths. This includes a unique identifier for each path, the time step, the location (cell ID and three-dimensional coordinates) on the specified surface (e.g., \code{bathy}) and the probability associated with that cell, given movement from the previous cell, in the following columns:
#' \describe{
#'     \item{path_id}{ An integer that uniquely defines each path.}
#'     \item{timestep}{ An integer that defines each time step. If an origin is supplied to \code{\link[flapper]{pf}} and \code{add_origin = TRUE} in \code{\link[flapper]{pf_simplify}}, then time step 0 refers to the origin. Later time steps refer to sequential depth observations.}
#'     \item{cell_id}{ An integer that defines the cells ID of the surface \code{\link[raster]{raster}} (over which paths were reconstructed) retained by the algorithm.}
#'     \item{cell_x}{ A double that defines the cell x coordinate.}
#'     \item{cell_y}{ A double that defines the cell y coordinate.}
#'     \item{cell_z}{ A double that defines the value of the surface \code{\link[raster]{raster}} in each cell. If \code{bathy} is unsupplied to \code{\link[flapper]{pf_simplify}} and unavailable in the inputted object, this is NA.}
#'     \item{cell_pr}{ A double that defines the probability of localisation in each cell, given the `intrinsic' probabilities associated with each cell and the immediately preceding sampled locations (or the origin)'.}
#' Rows are ordered by path and then time step.
#'   }
#'
#' @seealso \code{\link[flapper]{pf_path-class}} objects are derived from \code{\link[flapper]{pf_archive-class}} objects via \code{\link[flapper]{pf_simplify}}. \code{\link[flapper]{dat_dcpf_paths}} provides an example.
#'
#' @author Edward Lavender
#' @docType package
#' @name pf_path-class
NULL
