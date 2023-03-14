######################################
######################################
#### flapper-tips-parallel documentation

#' @title Parallelisation in \code{\link[flapper]{flapper}}
#' @description
#'
#' Parallelisation in the base R installation (since v. 2.14.0) can be implemented via the \code{\link[parallel]{parallel}} package. This is an amalgamation of two earlier packages (\code{snow} and \code{multicore}) that extends the \code{apply} family of functions with parallelisation. (Beyond the \code{apply} family, other packages further develop parallelisation, such as \code{foreach}, which integrates parallelisation into the familiar \code{for} loop.)
#'
#' Of particular interest in \code{\link[flapper]{flapper}} are the two parallelised versions of \code{\link[base]{lapply}} provided by \code{\link[parallel]{parallel}}: \code{\link[parallel]{parLapply}} and \code{\link[parallel]{mclapply}}. These functions implement parallelisation in two distinct ways.
#'
#' \code{\link[parallel]{parLapply}} implements parallelisation via socket clusters*, which essentially involves launching an additional R process on each core on which code is executed. This approach is initiated by the definition of a cluster object with \code{\link[parallel]{makeCluster}}. Since each R process is unique, objects defined in the `main' workspace that are required by the parallelised code have to be exported to each core (see \code{\link[parallel]{clusterEvalQ}} and \code{\link[parallel]{clusterExport}}) before that code is executed. During code execution, each R process occupies memory. At the end of the code, the cluster needs to be closed (see \code{\link[parallel]{stopCluster}}).
#'
#' In contrast, \code{\link[parallel]{mclapply}} implements parallelisation using fork clusters. Rather than launching an additional R process on each core, this approach clones the existing `parent' workspace on each core. Both `parent' and `child' processes refer to the same address in the system's memory. This simplifies implementation and reduces overhead costs: in \code{\link[parallel]{mclapply}}, the number of cores (or `child processes') is simply specified and everything else is taken care of. On POSIX systems, such as MacOS, this approach can be faster than the socket approach and is typically preferable. However, error messages may be less informative (see Examples). Forking is not available on Windows.
#'
#' In the \code{\link[flapper]{flapper}} package, parallelisation is implemented using \code{\link[pbapply]{pblapply}}, which is a convenient wrapper for the \code{\link[parallel]{parLapply}} and \code{\link[parallel]{mclapply}} functions, via a \code{cl} argument. If \code{cl} is a cluster object, a character vector of objects required for export can also be supplied via \code{varlist}. If supplied, \code{\link[flapper]{flapper}} functions will call \code{\link[parallel]{clusterExport}} before \code{\link[pbapply]{pblapply}}. The latter is then called and parallelisation is implemented via \code{\link[parallel]{parLapply}}. Following code execution, \code{\link[flapper]{flapper}} functions close the cluster. In contrast, if \code{cl} is an integer, \code{\link[flapper]{flapper}} functions pass this directly to \code{\link[pbapply]{pblapply}} and parallelisation is implemented via \code{\link[parallel]{mclapply}}.
#'
#' Despite potential benefits, the overhead costs associated with parallelisation can be substantial. To minimise overheads, some \code{\link[flapper]{flapper}} functions can implement parallelisation over `chunks', within which serial iteration via \code{\link[base]{lapply}} is implemented. For calls to \code{\link[pbapply]{pblapply}}, there is an additional overhead associated with the progress bar, which can be suppressed via \code{pbapply::pboptions(type = "none")}. Given overhead costs, in some situations, parallelisation may be undesirable. For computationally intensive functions in \code{\link[flapper]{flapper}} that support parallelisation, such as \code{\link[flapper]{ac}}, \code{\link[flapper]{dc}}, \code{\link[flapper]{acdc}} and \code{\link[flapper]{pf_simplify}}, time trials of alternative approaches using subsets of the data are strongly recommenced.
#'
#' In the future, parallelisation routines in \code{\link[flapper]{flapper}} may be improved to reduce memory requirements and expanded beyond `perfectly parallel' problems (in which each iteration is entirely independent) to other situations.
#'
#' *Forking can be implemented via \code{\link[parallel]{makeCluster}} for  \code{\link[parallel]{parLapply}} but it is usually preferable to call \code{\link[parallel]{mclapply}} directly (e.g., as in \code{\link[pbapply]{pblapply}}).
#'
#' @author Edward Lavender
#'
#' @examples
#' if (flapper_run_slow) {
#'   #### Background
#'   # In this example, we explore, for a small time series
#'   # ... of particles created via pf(), the speed of alternative
#'   # ... particle processing routines without/with parallelisation
#'   # ... and with sockets versus forking.
#'
#'   #### Implement pf() for example data using a large number of particles
#'   pf_args <- dat_dcpf_histories$args
#'   pf_args$n <- 200L
#'   pf_args$calc_distance_euclid_fast <- TRUE
#'   out_pf <- do.call(pf, pf_args)
#'
#'   #### Compare the speed of alternative parallelisation approaches
#'   # ... for particle processing
#'   if (requireNamespace("rbenchmark", quietly = TRUE)) {
#'     pb_op <- pbapply::pboptions(type = "none")
#'     rbenchmark::benchmark(
#'       baseline = pf_simplify(out_pf,
#'         cl = NULL,
#'         calc_distance = "lcp",
#'         verbose = FALSE
#'       ),
#'       fork = pf_simplify(out_pf,
#'         cl = 2L,
#'         calc_distance = "lcp",
#'         verbose = FALSE
#'       ),
#'       socket = pf_simplify(out_pf,
#'         calc_distance = "lcp",
#'         cl = parallel::makeCluster(2L),
#'         verbose = FALSE
#'       ),
#'       replications = 1L
#'     )
#'     pbapply::pboptions(pb_op)
#'   }
#' }
#'
#' #### Note a difference in error handling between socket clusters and forking:
#' \dontrun{
#' # Implement pf() with a small number of particles
#' pf_args <- dat_dcpf_histories$args
#' pf_args$n <- 10L
#' pf_args$calc_distance_euclid_fast <- TRUE
#' out_pf <- do.call(pf, pf_args)
#' # With a socket cluster, we get an informative error message that there are no
#' # ... particles that obey shortest-distance constraints
#' # ... (which is because we have implemented pf() using Euclidean sampling
#' # ... with only 10 particles).
#' pf_simplify(out_pf,
#'   calc_distance = "lcp",
#'   cl = parallel::makeCluster(2L)
#' )
#' # With forking, we get an apparently less helpful error message:
#' pf_simplify(out_pf,
#'   calc_distance = "lcp",
#'   cl = 2L
#' )
#' # This results from a difference in error handling
#' # ... between parLapply() and mclapply(), as shown here:
#' out <- pbapply::pblapply(1:10, cl = parallel::makeCluster(2L), function(i) {
#'   if (i == 5L) stop("error") else return(i)
#' })
#' out <- pbapply::pblapply(1:10, cl = 2L, function(i) {
#'   if (i == 5L) stop("error") else return(i)
#' })
#' # In the example above, connected cells are being identified at each time step
#' # ... based on previous cells. The function hits a point at which there are no
#' # ... allowed cells, and returns an error, but mclapply() handles this with
#' # ... tryCatch and attempts to continue the filtration process on the next time step.
#' # ... The filter operation then fails because it expects a dataframe, not a character
#' # ... (error message).
#' }
#'
#' @references
#'
#' Errickson, J. Parallel Processing in R. https://dept.stat.lsa.umich.edu/~jerrick/courses/stat701/notes/parallel.html. Accessed 29 November 2021.
#'
#' Hallquist, M. (2018). Parallel computing in R. https://psu-psychology.github.io/r-bootcamp-2018/talks/parallel_r.html. Accessed 29 November 2021.
#'
#' @name flapper-tips-parallel
NULL
