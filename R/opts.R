#####################################
#####################################
#### flapper_opts

#' @title Global example controls
#' @description These variables provide convenient switches for turning examples on/off during package testing. They can be overridden by the user.
#'
#' @details A set of logical variables:
#' \itemize{
#'   \item \code{\link[flapper]{flapper_run_parallel}} is a logical variable that defines whether or not to run selected parallel examples.
#'   \item \code{\link[flapper]{flapper_run_slow}} is logical variable that defines whether or not to run selected `slow' examples.
#' }
#'
#' @examples
#' #### Example (1): flapper_run_parallel
#' ## Illustration of intended usage for a parallelised function
#' if (flapper_run_parallel) {
#'   cl_lapply(1:10,
#'     function(x) x + 1,
#'     cl = parallel::makeCluster(2L)
#'   )
#' }
#'
#' #### Example (2): flapper_run_slow
#' ## Illustration of intended usage for a slow function
#' if (flapper_run_slow) {
#'   pf_args <- dat_dcpf_histories$args
#'   pf_args$n <- 100
#'   pf_args$calc_distance <- "lcp"
#'   pf_args$seed <- 1
#'   out_pf <- do.call(pf, pf_args)
#' }
#'
#' @author Edward Lavender
#' @name flapper_opts
NULL

#### flapper_run_parallel
#' @name flapper_opts
"flapper_run_parallel"

#### flapper_run_slow
#' @name flapper_opts
"flapper_run_slow"
