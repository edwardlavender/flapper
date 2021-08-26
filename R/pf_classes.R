######################################
######################################
#### pf_archive-class

#' @title "pf_archive" class
#' @description An S3 class that defines the object returned by \code{\link[flapper]{pf}} or \code{\link[flapper]{pf}} plus \code{\link[flapper]{pf_simplify}} with the \code{return = "archive"} argument.
#'
#' @return This object is a named list, with three elements, which records the particles sampled at each time step, the method and the parameters used in the call to \code{\link[flapper]{pf}}:
#'
#' \describe{
#'   \item{history}{A list of dataframes, one for each time step, that record particle samples. Each dataframe comprises \code{n} rows (one for each particle) and the following five columns:
#'   \itemize{
#'     \item{id_previous}{ An integer that uniquely defines a previous location on the \code{record} \code{\link[raster]{raster}}s. This is absent for the first time step. For subsequent time steps, this is \code{NA} if the fast Euclidean distances method has been used, which does not `remember' previous locations, but specified otherwise.)}
#'     \item{pr_current}{ A double that defines the movement probabilities associated with previous locations. This is absent for the first time step. For subsequent time steps, this is \code{NA} if the fast Euclidean distances method has been used, which does not `remember' previous locations, but specified otherwise.}
#'     \item{id_current}{ An integer that uniquely defines each location on the \code{record} \code{\link[raster]{raster}}s.}
#'     \item{pr_current} A double that defines the probability of movement into each cell.
#'     \item{timestep}{ An integer that defines each time step.}
#'   }
#' }
#'
#' \item{method}{A character that defines whether or not \code{history} was derived directly from \code{\link[flapper]{pf}} (\code{method = "pf"}), in which case \code{history} contains all of the particles sampled at each time step, or via \code{\link[flapper]{pf}} plus \code{\link[flapper]{pf_simplify}} with \code{return = "archive"} (\code{method = "pf_simplify"}), in which case \code{history} contains the subset of particles at each time step that were re-sampled at the next time step; for particles that were sampled multiple times on a given time step, this only contains the most probable sample (see \code{\link[flapper]{pf_simplify}}).}
#'
#' \item{args}{A named list that records the function arguments passed to \code{\link[flapper]{pf}}. This is as inputted to \code{\link[flapper]{pf}}, but with the `calc_distance_graph' elements added if unsupplied and applicable.}
#'
#' }
#'
#' @seealso \code{\link[flapper]{dat_dcpf_histories}} provides an example of a \code{\link[flapper]{pf_archive-class}} object. \code{\link[flapper]{pf}} and \code{\link[flapper]{pf_simplify}} return \code{\link[flapper]{pf_archive-class}} objects. \code{\link[flapper]{pf_simplify}} can also convert \code{\link[flapper]{pf_archive-class}} objects into \code{\link[flapper]{pf_path-class}} objects that comprise the reconstructed movement paths.
#'
#' @author Edward Lavender
#' @docType package
#' @name pf_archive-class
NULL


######################################
######################################
#### pf_path-class

#' @title "pf" class
#' @description An S3 class that defines the object returned by \code{\link[flapper]{pf_simplify}} with \code{return = "path"}, following \code{\link[flapper]{pf}}.
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
