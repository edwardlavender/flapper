######################################
######################################
#### make_matrix_ids()

#' @title Matricise individual deployment time series
#' @importFrom lubridate `%within%`
#'
#' @description This function creates a matrix that, for each time step (matrix row) in a sequence of user-defined times, defines whether or not each individual (matrix column) was at liberty. To implement the function, a dataframe with individual IDs and deployment start and end times must be supplied (via \code{ids}). The times for which to express whether or not each individual was at liberty are provided by optionally defining a \code{start} and \code{end} date (these can be taken from the range of deployment times in \code{ids} if unspecified) and the interval (\code{delta_t}) between time steps.
#'
#' @param ids A dataframe that defines individual IDs and deployment times. This must contain the following columns: an identifier for individuals (named 'individual_id'), the start time of individuals' deployment periods ('tag_start_date') and the end time of individuals' deployment periods ('tag_end_date'). Deployment times can be recorded as Date or POSIXct objects but, if the former is provided, they need to be coerced to POSIXct objects. This can be done automatically (see \code{as_POSIXct}).
#' @param start,end Date or POSIXct objects that define the start and end time. If unspecified, these are taken from the range of deployment times in \code{ids}.
#' @param delta_t A number or character that defines the time interval between successive time steps. This is passed to the 'by' argument of \code{\link[base]{seq.POSIXt}}.
#' @param as_POSIXct A function that coerces supplied any supplied times (\code{ids$tag_start_date}, \code{ids$tag_end_date}, \code{start} and \code{end}) that are not POSIXct objects to POSIXct objects.
#' @param set_names A logical variable that defines whether or not to set the row and column names of the matrix to the time steps and the individual IDs respectively.
#' @param ... Additional arguments (none implemented).
#'
#' @return The function returns a matrix with one row for each time step and one column for each individual. Each cell defines whether (1) or not (0) each individual was at liberty during that time step.
#'
#' @examples
#' dat_ids$tag_end_date <- as.Date("2017-06-02")
#' mat_hours <- make_matrix_ids(dat_ids, delat_t = "hours")
#' mat_days  <- make_matrix_ids(dat_ids, delta_t = "days")
#' utils::str(mat_hours)
#' utils::str(mat_days)
#' @author Edward Lavender
#' @export


make_matrix_ids <- function(ids,
                            start = NULL,
                            end = NULL,
                            delta_t = "120 mins",
                            as_POSIXct = as.POSIXct,
                            set_names = TRUE,...) {
  check_names(input = ids, req = c("individual_id", "tag_start_date", "tag_end_date"),
              extract_names = colnames, type = all)
  if(is.null(start)) start <- min(ids$tag_start_date, na.rm = TRUE)
  if(is.null(end)) end <- max(ids$tag_end_date, na.rm = TRUE)
  if(!is.null(as_POSIXct)) {
    if(!inherits(ids$tag_start_date, "POSIXct")) ids$tag_start_date <- as_POSIXct(ids$tag_start_date)
    if(!inherits(ids$tag_end_date, "POSIXct")) ids$tag_end_date <- as_POSIXct(ids$tag_end_date)
    if(!inherits(start, "POSIXct")) start <- as_POSIXct(start)
    if(!inherits(end, "POSIXct")) end <- as_POSIXct(end)
  }
  bin <- seq.POSIXt(start, end, delta_t)
  ids$interval <- lubridate::interval(ids$tag_start_date, ids$tag_end_date)
  mat <- matrix(NA, nrow = length(bin), ncol = length(unique(ids$individual_id)))
  for(j in 1:ncol(mat)) {
    mat[, j] <- (bin %within% ids$interval[j]) + 0
  }
  if(set_names) {
    rownames(mat) <- as.character(bin)
    colnames(mat) <- as.character(ids$individual_id)
  }
  return(mat)
}

