#' @title Parallelisation helpers
#'
#' @description A set of wrappers for \code{\link[parallel]{parallel}} functions that facilitate the implementation of parallel routines in functions via \code{\link[pbapply]{pblapply}}.
#'
#' @param x A \code{\link[base]{list}} over which to iterate.
#' @param fun,... A function that is applied to elements of \code{x} alongside any optional arguments to \code{fun}.
#' @param cl (optional) A cluster from \code{\link[parallel]{makeCluster}} or an integer that defines the number of child processes (see \code{\link[pbapply]{pblapply}}).
#' @param varlist (optional) A character vector of objects for export (see \code{\link[parallel]{clusterExport}}). If \code{cl} is a cluster, this may be required. Exported objects must be located in the global environment.
#' @param use_chunks A logical vector that defines whether to parallelise over `chunks' (\code{TRUE}) or over the elements of \code{x} (\code{FALSE}). If \code{use_chunks = TRUE}, \code{x} is split into \emph{n} chunks (one per core) that are processed in parallel; within each chunk \code{x} is updated iteratively.
#' @param length An integer that defines the number of elements in the iteration.
#'
#' @details
#'
#' \code{\link[flapper]{cl_lapply}} is a wrapper for \code{\link[pbapply]{pblapply}} that handles cluster checking, set up and closure, using the following functions:
#'
#' \itemize{
#'   \item \code{\link[flapper]{cl_check}} checks \code{cl} and \code{varlist} arguments, as inputted to a parent function. For example, if \code{cl = NULL}, \code{varlist} should also be \code{NULL}.
#'   \item  \code{\link[flapper]{cl_chunks}} defines a list, with one element for core specified, that contains an integer vector of the positions of an object over which to iterate serially in each chunk.
#'   \item \code{\link[flapper]{cl_export}} implements \code{\link[parallel]{clusterExport}} if both \code{cl} and \code{varlist} are specified.
#'   \item \code{\link[flapper]{cl_stop}} implements \code{\link[parallel]{stopCluster}} if \code{cl} is a cluster object from \code{\link[parallel]{makeCluster}}.
#' }
#'
#' @examples
#' #### Examples of cl_lapply()
#' # Implement cl_lapply() without cluster
#' z <- cl_lapply(1:10, function(x) x + 1)
#' # Implement cl_lapply() with forking (not on Windows)
#' z <- cl_lapply(1:10, function(x) x + 1, cl = 2L)
#' # Implement cl_lapply() with socket cluster
#' z <- cl_lapply(1:10, function(x) x + 1, cl = parallel::makeCluster(2L))
#'
#' #### Catch mistakes
#' z <- cl_lapply(1:10, function(x) x + 1, cl = 2L, varlist =  list())
#' z <- cl_lapply(1:10, function(x) x + 1, varlist = list())
#'
#' #### Compare time trials for chunk-wise versus element-wise parallelisation
#'
#' if(flapper_run_parallel){
#'
#' ## Background
#' # In this simple example, we will sample 'size' cells n times from a raster
#' # The output is a list of cell samples. We compare the time taken to complete
#' # sampling using different approaches.
#'
#' ## Define a dataframe of time trial scenarios
#' require(dplyr)
#' dat    <- expand.grid(n = 1e4,
#'                       method = c("socket", "fork"),
#'                       cores = 2L,
#'                       use_chunks = c(FALSE, TRUE),
#'                       time = NA)
#'
#' ## Estimate the duration of each scenario
#' dat_by_trial <-
#'   lapply(split(dat, seq_len(nrow(dat))), function(d){
#'     if(d$method == "socket"){
#'       t1 <- Sys.time()
#'       z  <- cl_lapply(x = 1:d$n,
#'                       fun = function(i)
#'                         raster::sampleRandom(flapper::dat_gebco, size = 5),
#'                       cl = parallel::makeCluster(d$cores),
#'                       use_chunks = d$use_chunks)
#'       t2 <- Sys.time()
#'     } else if(d$method == "fork"){
#'       t1 <- Sys.time()
#'       z  <- cl_lapply(x = 1:d$n,
#'                       fun = function(i)
#'                         raster::sampleRandom(flapper::dat_gebco, size = 5),
#'                       cl = d$cores,
#'                       use_chunks = d$use_chunks)
#'       t2 <- Sys.time()
#'     }
#'     d$time <- as.numeric(difftime(t2, t1, "secs"))
#'     return(d)
#'   })
#'
#' ## Examine the results
#' dat_for_trials <-
#'   dat_by_trial %>%
#'   dplyr::bind_rows() %>%
#'   dplyr::arrange(.data$n, .data$time) %>%
#'   print()
#'
#' }
#'
#' @return
#' \itemize{
#'   \item \code{\link[flapper]{cl_lapply}} returns a list.
#'   \item \code{\link[flapper]{cl_chunks}} returns a list of integers.
#'   \item \code{\link[flapper]{cl_check}}, \code{\link[flapper]{cl_export}} and \code{\link[flapper]{cl_stop}} return \code{invisible()}.
#' }
#'
#' @author Edward Lavender
#' @name cl
NULL


#### cl_lapply()
#' @rdname cl
#' @export

cl_lapply <- function(x, fun, ..., cl = NULL, varlist = NULL, use_chunks = FALSE){
  # Check cluster
  cl_check(cl = cl, varlist = varlist)
  if(use_chunks){
    # Define list of indices by chunk
    index_by_chunk <- cl_chunks(cl = cl, length = length(x))
    # Loop over chunks in parallel
    cl_export(cl = cl, varlist = varlist)
    y_by_chunks <- pbapply::pblapply(index_by_chunk, cl = cl, function(index_for_chunk){
      # Get indices for chunk
      x_for_chunk <- x[index_for_chunk]
      # Loop over chunk in serial
      y_for_chunk <- lapply(x_for_chunk, function(xi) return(fun(xi)))
      return(y_for_chunk)
    })
    # Close cluster
    cl_stop(cl = cl)
    # Flatten list-by-chunk into a single level list
    y <- purrr::flatten(y_by_chunks)
  } else {
    # Loop over x elements in parallel
    cl_export(cl = cl, varlist = varlist)
    y <- pbapply::pblapply(x, cl = cl, function(xi) return(fun(xi)))
    cl_stop(cl = cl)
  }
  return(y)
}

#### cl_check()
#' @rdname cl
#' @export

cl_check <- function(cl = NULL, varlist = NULL){
  if(is.null(cl)){
    if(!is.null(varlist)){
      warning("'cl' is NULL: input to 'varlist' ignored.",
               immediate. = TRUE, call. = FALSE)
    }
  } else {
    if(!inherits(cl, "cluster")){
      if(.Platform$OS.type == "windows"){
        warning("Integer specifications for 'cl' (i.e., forking) on Windows are not supported.",
                immediate. = TRUE, call. = FALSE)
      }
      if(!is.null(varlist)){
        warning("'cl' is an integer: input to 'varlist' ignored.",
                immediate. = TRUE, call. = FALSE)
      }
    }
  }
  return(invisible())
}


#### cl_chunks()
#' @rdname cl
#' @export

cl_chunks <- function(cl = NULL, length){
  if(is.null(cl)){
    chunks <- 1L
  } else {
    if(inherits(cl, "cluster")) chunks <- length(cl) else chunks <- cl
  }
  index <- parallel::splitIndices(length, chunks)
  return(invisible(index))
}


#### cl_export()
#' @rdname cl
#' @export

cl_export <- function(cl = NULL, varlist = NULL){
  if(!is.null(cl) && inherits(cl, "cluster") && !is.null(varlist))
    parallel::clusterExport(cl = cl, varlist = varlist)
  return(invisible())
}


#### cl_stop()
#' @rdname cl
#' @export

cl_stop <- function(cl = NULL){
  if(!is.null(cl) && inherits(cl, "cluster"))
    parallel::stopCluster(cl = cl)
  return(invisible())
}

