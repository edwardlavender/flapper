#' @title flapper: An R Package to Explore Animal Space Use Within Passive Acoustic Telemetry Arrays
#' @author Edward Lavender
#'
#' @description flapper is an R package designed to facilitate the use of passive acoustic telemetry (PAT) data for ecological inferences, especially those pertaining to animal space use. This includes functions for processing PAT data, spatial tools, new algorithms for inferring space use and simulations designed to evaluate the efficacy of existing and new algorithms for inferring space use. Package development has been motivated by the collection of PAT data for a Critically Endangered benthopelagic elasmobranch off the West Coast of Scotland.
#'
#' @section Data processing:
#' A set of functions is used for processing PAT data, including exploring false detections, assembling detection -- transmission (i.e., range testing) datasets and implementing quality checks.
#'
#' @section Spatial tools:
#' A set of functions is used to implement standard spatial operations, distance calculations and the (rapid) computation of least-cost pathways over a surface between origin and destination coordinates.
#'
#' @section Space use algorithms:
#' A set of functions is used to implement existing and new algorithms for inferring patterns of space use from PAT and/or archival data.
#'
#' @section Simulations:
#' A set of functions is used to simulate PAT arrays, movement and detections to evaluate data processing and space use algorithms.
#'
#' @seealso https://github.com/edwardlavender/flapper
#'
#' @docType package
#' @name flapper
NULL
