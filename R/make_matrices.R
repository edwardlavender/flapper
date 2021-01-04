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


######################################
######################################
#### make_matrix_receivers()

#' @title Matricise receiver deployment time series
#' @importFrom lubridate `%within%`
#'
#' @description This function creates a matrix that, for each time step (matrix row) in a sequence of user-defined times, defines whether or not each receiver (matrix column) was active. To implement the function, a dataframe with receiver IDs and deployment start and end times must be supplied (via \code{moorings}). The times for which to express whether or not each receiver was active are provided by optionally defining a \code{start} and \code{end} date (these can be taken from the range of deployment times in \code{moorings} if unspecified) and the interval (\code{delta_t}) between time steps.
#'
#' @param moorings A dataframe that defines receiver IDs and deployment times. This must contain the following columns: an identifier for receivers (named 'receiver_id'), the start time of receiver' deployment periods ('receiver_start_date') and the end time of receivers' deployment periods ('receiver_end_date'). Deployment times can be recorded as Date or POSIXct objects but, if the former is provided, they need to be coerced to POSIXct objects. This can be done automatically (see \code{as_POSIXct}).
#' @param start,end Date or POSIXct objects that define the start and end time. If unspecified, these are taken from the range of deployment times in \code{moorings}.
#' @param delta_t A number or character that defines the time interval between successive time steps. This is passed to the 'by' argument of \code{\link[base]{seq.POSIXt}}.
#' @param as_POSIXct A function that coerces supplied any supplied times (\code{moorings$receiver_start_date}, \code{moorings$receiver_end_date}, \code{start} and \code{end}) that are not POSIXct objects to POSIXct objects.
#' @param set_names A logical variable that defines whether or not to set the row and column names of the matrix to the time steps and the receiver IDs respectively.
#' @param ... Additional arguments (none implemented).
#'
#' @return The function returns a matrix with one row for each time step and one column for each receiver. Each cell defines whether (1) or not (0) each receiver was at active during that time step.
#'
#' @examples
#' mat_hours <- make_matrix_receivers(dat_moorings, delat_t = "hours")
#' mat_days  <- make_matrix_receivers(dat_moorings, delta_t = "days")
#' utils::str(mat_hours)
#' utils::str(mat_days)
#' @author Edward Lavender
#' @export

make_matrix_receivers <- function(moorings,
                                  start = NULL,
                                  end = NULL,
                                  delta_t = "120 mins",
                                  as_POSIXct = as.POSIXct,
                                  set_names = TRUE,...){
  check_names(input = moorings, req = c("receiver_id", "receiver_start_date", "receiver_end_date"),
              extract_names = colnames, type = all)
  if(is.null(start)) start <- min(moorings$receiver_start_date, na.rm = TRUE)
  if(is.null(end)) end <- max(moorings$receiver_end_date, na.rm = TRUE)
  if(!is.null(as_POSIXct)) {
    if(!inherits(moorings$receiver_start_date, "POSIXct")) moorings$receiver_start_date <- as_POSIXct(moorings$receiver_start_date)
    if(!inherits(moorings$receiver_end_date, "POSIXct")) moorings$receiver_end_date <- as_POSIXct(moorings$receiver_end_date)
    if(!inherits(start, "POSIXct")) start <- as_POSIXct(start)
    if(!inherits(end, "POSIXct")) end <- as_POSIXct(end)
  }
  bin <- seq.POSIXt(start, end, delta_t)
  moorings$interval <- lubridate::interval(moorings$receiver_start_date, moorings$receiver_end_date)
  mat <- matrix(NA, nrow = length(bin), ncol = length(unique(moorings$receiver_id)))
  for(j in 1:ncol(mat)) {
    mat[, j] <- (bin %within% moorings$interval[j]) + 0
  }
  if(set_names) {
    rownames(mat) <- as.character(bin)
    colnames(mat) <- as.character(moorings$receiver_id)
  }
  return(mat)
}


######################################
######################################
#### make_matrix_detections()

#' @title Matricise detection time series
#' @description This function creates a list of matrices that, for each individual (list element), defines the number of detections of that individual in each time interval (matrix row) at each receiver (matrix column). To implement the function, a dataframe with acoustic detection time series must be provided via \code{acoustics}. The time intervals over which to count detections are provided by optionally defining a \code{start} and \code{end} date (these can be taken from the range of times in \code{acoustics} if unspecified) and the interval (\code{delta_t}) between time steps. By default, matrix elements that are 'outside' individual or receiver deployment periods are defined as 0 (not detected) but can be changed to another value (e.g., NA) via \code{set_outside}. In this case, the \code{acoustics} dataframe also needs to include the deployment times for each individual and an additional dataframe must be supplied with the same information for receivers via \code{moorings}.
#'
#' @param acoustics A dataframe that defines passive acoustic telemetry detection time series. This should contain the following columns: a vector of individual IDs, named 'individual_id'; a (factor) vector of receiver IDs, named 'receiver_id'. If \code{set_outside} is specified, this should also contain POSIXct vectors of the start and end time of each individual's time at liberty, named 'tag_start_date' and 'tag_end_date' respectively.
#' @param moorings (optional)  If \code{set_outside} is specified, \code{moorings} is dataframe that defines passive acoustic telemetry receiver metadata. This should contain the following columns: a vector of receiver IDs, named 'receiver_id (as in \code{acoustics}); and POSIXct vectors of the start and end times of each receiver's deployment time, named 'receiver_start_date' and 'receiver_end_date' respectively.
#' @param start,end POSIXct objects that define the start and end time. If unspecified, these are taken from the range of detection times in \code{acoustics}.
#' @param delta_t A number or character that defines the time interval between successive time steps. This is passed to the 'by' argument of \code{\link[base]{seq.POSIXt}}.
#' @param simplify (optional) A function that simplifies detection matrices, such as \code{function(x) (x > 0) + 0} to convert counts into a boolean outcomes that define whether or a detection was made.
#' @param set_outside (optional) A value (e.g., \code{NA}) that is assigned to any matrix element that is outside of an individual's or receiver's deployment period when detection was not possible.
#' @param as_POSIXct A function that coerces supplied any supplied times that are not POSIXct objects to POSIXct objects.
#' @param set_names A logical variable that defines whether or not to set the row and column names of the matrix to the time steps and the receiver IDs respectively.
#' @param verbose A logical variable that defines whether or not to print messages to the console to relay function progress.
#' @param ... Additional arguments (none implemented).
#'
#' @return A matrix, or a list of matrices (one for each individual in \code{acoustics} is there is more than one individual) with one column for each time step and one column for each receiver. Each cell defines whether (1) or not (0) the individual was detected (or, if \code{set_outside} is supplied, it could not have been detected) during that time step. All matrices are expressed across the same sequence of time steps and receivers.
#'
#' @examples
#' #### Example (1) Construct matrix from detected individuals using detection time series
#' dat_acoustics$receiver_id <- factor(dat_acoustics$receiver_id)
#' mat_by_id <- make_matrix_detections(dat_acoustics)
#' summary(mat_by_id)
#' range(mat_by_id[[1]])
#'
#' #### Example (2) Construct matrix across all receivers and use set_outside
#' dat_moorings$receiver_id <- factor(dat_moorings$receiver_id)
#' dat_acoustics$receiver_id <- factor(dat_acoustics$receiver_id,
#'                                     levels = levels(dat_moorings$receiver_id))
#'
#' match_index <- match(dat_acoustics$individual_id, dat_ids$individual_id)
#' dat_acoustics$tag_start_date <- dat_ids$tag_start_date[match_index]
#' dat_acoustics$tag_end_date   <- as.Date("2017-06-02")
#' mat_by_id <- make_matrix_detections(dat_acoustics,
#'                                     moorings = dat_moorings,
#'                                     set_outside = NA)
#' summary(mat_by_id)
#'
#' @author Edward Lavender
#' @export

make_matrix_detections <- function(acoustics,
                                   moorings = NULL,
                                   start = NULL, end = NULL,
                                   delta_t = "120 mins",
                                   simplify = NULL,
                                   set_outside = NULL,
                                   as_POSIXct = as.POSIXct,
                                   set_names = TRUE,
                                   verbose = TRUE,...) {

  #### Initiate function
  t_onset <- Sys.time()
  cat_to_console <- function(..., show = verbose) if(show) cat(paste(..., "\n"))
  cat_to_console(paste0("flapper::make_matrix_detections() called (@ ", t_onset, ")..."))

  #### Checks
  cat_to_console("... Checking user inputs...")
  by_id <- "individual_id"
  check_names(input = acoustics, req = c(by_id, "timestamp", "receiver_id"))
  check_class(input = acoustics$timestamp, to_class = "POSIXct", type = "stop")
  acoustics$receiver_id <- check_class(input = acoustics$receiver_id, to_class = "factor",
                                       type = "warning", coerce_input = factor)

  if(is.null(start)) {
    start <- min(acoustics$timestamp, na.rm = TRUE)
  } else {
    start <- check_class(input = start, to_class = "POSIXct",
                         type = "warning", coerce_input = function(x) as_POSIXct(x))
  }
  if(is.null(end)) {
    end <- max(acoustics$timestamp, na.rm = TRUE)
  } else {
    end <- check_class(input = end, to_class = "POSIXct",
                         type = "warning", coerce_input = function(x) as_POSIXct(x))
  }
  if(!is.null(set_outside)) {
    check_names(input = acoustics, req = c("tag_start_date", "tag_end_date"))
    acoustics$tag_start_date <- check_class(input = acoustics$tag_start_date, to_class = "POSIXct",
                                            type = "warning", coerce_input = function(x) as_POSIXct(x))
    acoustics$tag_end_date <- check_class(input = acoustics$tag_end_date, to_class = "POSIXct",
                                          type = "warning", coerce_input = function(x) as_POSIXct(x))
    check_names(input = moorings, req = c("receiver_id", "receiver_start_date", "receiver_end_date"))
    moorings$receiver_id <- check_class(input = moorings$receiver_id, to_class = "factor",
                                        type = "warning", coerce_input = factor) # needed
    moorings <- moorings[order(moorings$receiver_id), ]
    moorings$receiver_start_date <- check_class(input = moorings$receiver_start_date, to_class = "POSIXct",
                                                type = "warning", coerce_input = function(x) as_POSIXct(x))
    moorings$receiver_end_date <- check_class(input = moorings$receiver_end_date, to_class = "POSIXct",
                                              type = "warning", coerce_input = function(x) as_POSIXct(x))
  }

  #### Define bins for which to determine detections
  cat_to_console("... Defining time bins given 'delta_t'...")

  bin <- seq.POSIXt(start, end, delta_t)
  acoustics$bin <- cut(acoustics$timestamp, bin)
  acoustics$bin <- factor(acoustics$bin, levels = levels(factor(bin)))

  #### Determine receiver activity in each bin
  if(!is.null(set_outside)) {
    cat_to_console("... Making receiver matrix...")
    if(!all(moorings$receiver_id %in% levels(acoustics$receiver_id))) {
      warning("Not all receivers in moorings$receiver_id are included within levels(acoustics$receiver_id).")
      moorings <- moorings[which(moorings$receiver_id %in% levels(acoustics$receiver_id)), ]
    }
    if(!all(levels(acoustics$receiver_id %in% moorings$receiver_id))) {
      stop("Not all receivers in levels(acoustics$receiver_id) are in moorings$receiver_id")
    }
    if(any(duplicated(moorings$receiver_id))) {
      warning("Duplicated receiver IDs in moorings will be dropped.")
      moorings <- moorings$receiver_id[!duplicated(moorings$receiver_id), ]
    }
    moorings$receiver_id <- factor(moorings$receiver_id, levels = levels(acoustics$receiver_id))
    moorings <- moorings[moorings$receiver_id, ]
    moorings$interval <- lubridate::interval(moorings$receiver_start_date, moorings$receiver_end_date)
    receiver_mat <- matrix(NA, nrow = length(bin), ncol = length(moorings$receiver_id))
    for(j in 1:ncol(receiver_mat)) {
      receiver_mat[, j] <- !(bin %within% moorings$interval[j])
    }
    outside_receiver_int <- is.na(receiver_mat)
  }

  #### Compute matrices
  cat_to_console("... Making detection matrix...")
  # Define a list of outcomes for each individual
  acoustics_ls <- split(acoustics, acoustics[, by_id])
  det_ls <- pbapply::pblapply(acoustics_ls, function(d){
    # Define detection matrix (timestamps x receivers)
    detection <- table(d$bin, d$receiver_id)
    if(!is.null(simplify)) detection <- simplify(detection)
    # Process detection matrix
    if(!is.null(set_outside)) {
      # Determine any intervals outside of the individual's deployment time and force NA
      id_int <- lubridate::interval(d$tag_start_date[1], d$tag_end_date[1])
      outside_id_int <- !(bin %within% id_int)
      if(any(outside_id_int)) detection[which(outside_id_int), ] <- set_outside
      # Set bins outside of receivers' deployment time to NA
      if(any(outside_receiver_int)) detection[outside_receiver_int] <- set_outside
    }
    if(set_names) {
      rownames(detection) <- as.character(bin)
      colnames(detection) <- levels(acoustics$receiver_id)
    }
    return(detection)
  })

  #### Return outputs
  out <- det_ls
  if(length(out) == 1) out <- out[[1]]
  t_end <- Sys.time()
  duration <- difftime(t_end, t_onset, units = "mins")
  cat_to_console(paste0("... flapper::make_matrix_detections() call completed (@ ", t_end, ") after ~", round(duration, digits = 2), " minutes."))
  return(out)
}
