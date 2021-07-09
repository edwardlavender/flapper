######################################
######################################
#### make_matrix_ids()

#' @title Matricise individual deployment time series
#' @importFrom lubridate `%within%`
#'
#' @description This function creates a matrix that, for each time step (matrix row) in a sequence of user-defined times, defines whether or not each individual (matrix column) was at liberty. To implement the function, a dataframe with individual IDs and deployment start and end times must be supplied (via \code{ids}). The times for which to express whether or not each individual was at liberty are provided by optionally defining a \code{start} and \code{end} date (these can be taken from the range of deployment times in \code{ids} if unspecified) and the interval (\code{delta_t}) between time steps.
#'
#' @param ids A dataframe that defines individual IDs and deployment times. This must contain the following columns: an identifier for individuals (named `individual_id'), the start time of individuals' deployment periods (`tag_start_date') and the end time of individuals' deployment periods (`tag_end_date'). Deployment times can be recorded as Date or POSIXct objects but, if the former is provided, they need to be coerced to POSIXct objects. This can be done automatically (see \code{as_POSIXct}).
#' @param start,end Date or POSIXct objects that define the start and end time. If unspecified, these are taken from the range of deployment times in \code{ids}.
#' @param delta_t A number or character that defines the time interval between successive time steps. This is passed to the `by' argument of \code{\link[base]{seq.POSIXt}}.
#' @param as_POSIXct A function that coerces any supplied times (\code{ids$tag_start_date}, \code{ids$tag_end_date}, \code{start} and \code{end}) that are not POSIXct objects to POSIXct objects.
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
#' @description This function creates a matrix that, for each time step (matrix row) in a sequence of user-defined times, defines whether or not each receiver (matrix column) was active. To implement the function, a dataframe with receiver IDs and deployment start and end times must be supplied (via \code{moorings}). Servicing dates can also be accounted for via a dataframe with receiver IDs and servicing times (\code{services}). The times for which to express whether or not each receiver was active are provided by optionally defining a \code{start} and \code{end} date (these can be taken from the range of deployment times in \code{moorings} if unspecified) and the interval (\code{delta_t}) between time steps.
#'
#' @param moorings A dataframe that defines receiver IDs and deployment times. This must contain the following columns: an identifier for receivers (named `receiver_id'), the start time of receiver' deployment periods (`receiver_start_date') and the end time of receivers' deployment periods (`receiver_end_date') (see \code{\link[flapper]{dat_moorings}} for an example). Deployment times can be recorded as Date or POSIXct objects.
#' @param services (optional) A dataframe that defines receiver IDs and servicing dates (times during the deployment period of a receiver when it was not active due to servicing). If provided, this must contain the following columns: an identifier for serviced receivers (named `receiver_id') and two columns that define the time of the service(s) (`service_start_date' and `service_end_date'). Times can be recorded as Date or POSIXct objects. Before/after service events, receivers are assumed to have been deployed in the same locations; receiver deployments in different locations before/after servicing should be treated as distinct deployments in \code{moorings}.
#' @param start,end Date or POSIXct objects that define the start and end time. If unspecified, these are taken from the range of deployment times in \code{moorings}.
#' @param delta_t A number or character that defines the time interval between successive time steps. This is passed to the `by' argument of \code{\link[base]{seq.POSIXt}} or \code{\link[base]{seq.Date}} (depending on \code{as_POSIXct}, below).
#' @param as_POSIXct (optional) A function that coerces supplied any supplied times (\code{moorings$receiver_start_date}, \code{moorings$receiver_end_date}, \code{services$service_start_date}, \code{services$service_end_date}, \code{start} and \code{end}) that are not POSIXct objects to POSIXct objects. This can be suppressed via \code{as_POSIXct = NULL} if supplied times are Date objects and \code{delta_t} is not less than one day.
#' @param set_names A logical variable that defines whether or not to set the row and column names of the matrix to the time steps and the receiver IDs respectively.
#' @param ... Additional arguments (none implemented).
#'
#' @return The function returns a matrix with one row for each time step and one column for each receiver. Each cell defines whether (1) or not (0) each receiver was at active during that time step. A `bins' attribute is included, which defines the time steps as a Date or POSIXct vector.
#'
#' @examples
#' #### Example (1): Illustration using fake data
#'
#' ## Define some example 'moorings' data
#' # ... with receiver IDs and deployment times
#' moorings <- data.frame(receiver_id = c(1, 2, 3, 4, 5),
#'                        receiver_start_date = as.Date(c("2016-01-01",
#'                                                        "2016-01-02",
#'                                                        "2016-01-03",
#'                                                        "2016-01-04",
#'                                                        "2016-01-05")),
#'                        receiver_end_date = as.Date(c("2016-01-06",
#'                                                      "2016-01-07",
#'                                                      "2016-01-08",
#'                                                      "2016-01-09",
#'                                                      "2016-01-09"))
#'                        )
#'
#' ## Define some example 'servicing' data
#' # ... with receiver IDs and servicing times
#' # ... Here, receiver 1 was serviced twice
#' # ... ... from 2016-01-02--3 and 2016-01-04--5
#' # ... and receiver 5 was serviced
#' # ... ... on 2016-01-08.
#' services <- data.frame(receiver_id = c(1, 1, 5),
#'                        service_start_date = as.Date(c("2016-01-02",
#'                                                       "2016-01-04",
#'                                                       "2016-01-08")),
#'                        service_end_date = as.Date(c("2016-01-03",
#'                                                     "2016-01-05",
#'                                                     "2016-01-08"))
#'                        )
#'
#' ## Get daily receiver status (0, 1) matrix
#' make_matrix_receivers(moorings, delta_t = "days", as_POSIXct = NULL)
#'
#' ## Get daily receiver status (0, 1) matrix
#' # ... accounting for servicing dates
#' make_matrix_receivers(moorings, services, delta_t = "days", as_POSIXct = NULL)
#'
#' #### Example (2): Illustration using actual data
#' # ... for different time windows
#' mat_days  <- make_matrix_receivers(dat_moorings, delta_t = "days", as_POSIXct = NULL)
#' mat_hours <- make_matrix_receivers(dat_moorings, delta_t = "hours")
#' utils::str(mat_days)
#' utils::str(mat_hours)
#'
#' @author Edward Lavender
#' @export

make_matrix_receivers <- function(moorings,
                                  services = NULL,
                                  start = NULL,
                                  end = NULL,
                                  delta_t = "120 mins",
                                  as_POSIXct = as.POSIXct,
                                  set_names = TRUE,...){

  #### Check dates
  check_names(input = moorings, req = c("receiver_id", "receiver_start_date", "receiver_end_date"),
              extract_names = colnames, type = all)
  if(!is.null(services)){
    check_names(input = services, req = c("receiver_id", "service_start_date", "service_end_date"),
                extract_names = colnames, type = all)
    if(!all(unique(services$receiver_id) %in% unique(moorings$receiver_id))){
      message("Not all receivers in services$receiver_id are in moorings$receiver_id.")
    }
  }

  #### Process dates, if necessary
  if(is.null(start)) start <- min(moorings$receiver_start_date, na.rm = TRUE)
  if(is.null(end)) end <- max(moorings$receiver_end_date, na.rm = TRUE)
  if(!is.null(as_POSIXct)) {
    if(!inherits(moorings$receiver_start_date, "POSIXct")) moorings$receiver_start_date <- as_POSIXct(moorings$receiver_start_date)
    if(!inherits(moorings$receiver_end_date, "POSIXct")) moorings$receiver_end_date <- as_POSIXct(moorings$receiver_end_date)
    if(!inherits(start, "POSIXct")) start <- as_POSIXct(start)
    if(!inherits(end, "POSIXct")) end <- as_POSIXct(end)
    if(!is.null(services)){
      if(!inherits(services$service_start_date, "POSIXct")) services$service_start_date <- as_POSIXct(services$service_start_date)
      if(!inherits(services$service_end_date, "POSIXct")) services$service_end_date <- as_POSIXct(services$service_end_date)
    }
  }

  #### Define time steps over which to consider operations and define operational intervals
  # Define time steps over which to check operational status
  bin <- seq(start, end, delta_t)
  # Define intervals of moorings deployment
  moorings$interval <- lubridate::interval(moorings$receiver_start_date, moorings$receiver_end_date)
  # Define a corresponding list, in the same order as receivers in moorings, with servicing intervals
  if(!is.null(services)) {
    services$interval <- lubridate::interval(services$service_start_date, services$service_end_date)
    services_ls <- lapply(split(moorings, 1:nrow(moorings)), function(d){
      out <- NULL
      if(d$receiver_id %in% services$receiver_id){
        out <- services[which(services$receiver_id %in% moorings$receiver_id), ]
      }
      return(out)
    })
  }

  #### Define matrix of activity status (0, 1)
  # Blank matrix
  mat <- matrix(NA, nrow = length(bin), ncol = length(unique(moorings$receiver_id)))
  attr(mat, "bins") <- bin
  # Fill matrix
  for(j in 1:ncol(mat)) {
    # Fill matrix with 0,1 according to the overlap between bins and the deployment interval
    mat[, j] <- (bin %within% moorings$interval[j]) + 0
    # Post-hoc adjustment to suppress any time steps (bins) during this interval when receivers were being serviced)
    if(!is.null(services)){
      if(!is.null(services_ls[[j]])){
        services_for_j <- services_ls[[j]]
        for(k in 1:nrow(services_for_j)){
          mat[(bin %within% services_for_j$interval[k]), j] <- 0
        }
      }
    }
  }
  # Define matrix names
  if(set_names) {
    rownames(mat) <- as.character(bin)
    colnames(mat) <- as.character(moorings$receiver_id)
  }

  #### Return outputs
  return(mat)
}


######################################
######################################
#### make_matrix_detections()

#' @title Matricise detection time series
#' @description This function creates a list of matrices that, for each individual (list element), defines the number of detections of that individual in each time interval (matrix row) at each receiver (matrix column). To implement the function, a dataframe with acoustic detection time series must be provided via \code{acoustics}. The time intervals over which to count detections are provided by optionally defining a \code{start} and \code{end} date (these can be taken from the range of times in \code{acoustics} if unspecified) and the interval (\code{delta_t}) between time steps. By default, matrix elements that are `outside' individual or receiver deployment periods are defined as 0 (not detected) but can be changed to another value (e.g., NA) via \code{set_outside}. In this case, the \code{acoustics} dataframe also needs to include the deployment times for each individual and an additional dataframe must be supplied with the same information for receivers via \code{moorings}.
#'
#' @param acoustics A dataframe that defines passive acoustic telemetry detection time series. This should contain the following columns: a vector of individual IDs, named `individual_id'; a (factor) vector of receiver IDs, named `receiver_id'. If \code{set_outside} is specified, this should also contain POSIXct vectors of the start and end time of each individual's time at liberty, named `tag_start_date' and `tag_end_date' respectively.
#' @param moorings (optional)  If \code{set_outside} is specified, \code{moorings} is dataframe that defines passive acoustic telemetry receiver metadata. This should contain the following columns: a vector of receiver IDs, named `receiver_id (as in \code{acoustics}); and POSIXct vectors of the start and end times of each receiver's deployment time, named `receiver_start_date' and `receiver_end_date' respectively.
#' @param start,end POSIXct objects that define the start and end time. If unspecified, these are taken from the range of detection times in \code{acoustics}.
#' @param delta_t A number or character that defines the time interval between successive time steps. This is passed to the `by' argument of \code{\link[base]{seq.POSIXt}}.
#' @param simplify (optional) A function that simplifies detection matrices, such as \code{function(x) (x > 0) + 0} to convert counts into a boolean outcomes that define whether or a detection was made.
#' @param set_outside (optional) A value (e.g., \code{NA}) that is assigned to any matrix element that is outside of an individual's or receiver's deployment period when detection was not possible.
#' @param as_POSIXct A function that coerces any supplied times that are not POSIXct objects to POSIXct objects.
#' @param set_names A logical variable that defines whether or not to set the row and column names of the matrix to the time steps and the receiver IDs respectively.
#' @param verbose A logical variable that defines whether or not to print messages to the console to relay function progress.
#' @param ... Additional arguments (none implemented).
#'
#' @return A matrix, or a list of matrices (one for each individual in \code{acoustics} if there is more than one individual) with one column for each time step and one column for each receiver. Each cell defines whether (1) or not (0) the individual was detected (or, if \code{set_outside} is supplied, it could not have been detected) during that time step. All matrices are expressed across the same sequence of time steps and receivers.
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


#######################################
#######################################
#### make_matrix_cooccurence()

#' @title Compute a detection history similarity matrix
#' @description The function computes a detection history similarity matrix. For all combinations of individuals, this shows the total number (or percentage) of detections `nearby' in space and time, which can help to elucidate possible interactions among individuals that affect space use (see Details). To compute this matrix, the function pairs detections for each individual with the detections nearest in time for each other individual. The function computes the time (minutes) between paired detection time series, and the distance (m) between the receiver(s) at which paired detections occurred, dropping any detection pairs that are further apart in time or space than user-defined thresholds (which depend on the mobility of the species under investigation). For each combination of individuals, the function returns total number (or percentage) of detections that are closely associated in time and space. For very large combinations of individuals, especially those with long, overlapping time series, the function may take some time to run; therefore, testing the function on a small subset of individuals first is advisable. Parallelisation can be used to improve computation time. Similarity matrices can be visualised with \code{\link[prettyGraphics]{pretty_mat}}.
#'
#' @param acoustics_ls A list of dataframes, with one element for each individual, which contains each individual's detection time series. Each dataframe must include the following columns: `individual_id', a factor which specifies unique individuals; `timestamp', a POSIXct object which specifies the time of each detection; `receiver_long', the longitude (decimal degrees) of the receiver(s) at the individual was detected; and `receiver_lat', the latitude (decimal degrees) of the receiver(s) at which individual was detected. Each dataframe should be ordered by `individual_id' and then by `timestamp'. Careful ordering of `individual_id' factor levels (e.g. perhaps by population group, then by the number of detections of each individual) can aid visualisation of similarity matrices, in which the order or rows/columns corresponds directly to the order of individuals in \code{acoustics_ls}. Sequential elements in \code{acoustics_ls} should correspond to sequential factor levels for `individual_id', which should be the same across all dataframes.
#' @param thresh_time A number which specifies the time, in minutes, after which detections at nearby receivers are excluded.
#' @param thresh_dist A number which specifies the (Euclidean) distance between receivers, in metres, beyond which detections are excluded (see Details).
#' @param cl (optional) A cluster object created by \code{\link[parallel]{makeCluster}}. This is required if you want to run the algorithm in parallel. If supplied, the connection to the cluster is stopped within the function.
#' @param varlist (optional) A character vector of names of objects to export, to be passed to the \code{varlist} argument of \code{\link[parallel]{clusterExport}}. This may be required if \code{cl} is supplied. Exported objects must be located in the global environment.
#' @param output A number which specifies the output type. If \code{output = 1}, the function returns a (usually) symmetric similarity matrix in which each number represents the number of detections nearby in space and time for each pair of individuals. Row and column names are assigned from the `individual_id' column in \code{acoustics_ls} dataframes. This matrix is usually symmetric, but this is not necessarily the case for data collected from tags which transmit at random intervals around a nominal delay: under this scenario, the tag for a given individual (i) may transmit multiple signals in the space of time that the tag for another individual (j) only releases a single transmission. In this case, the pairing i,j will comprise all unique transmissions for individual i, paired with the nearest observations for individual j, some of which will be duplicate observations. Therefore, this pairing will contain more `shared observations' than the pairing j,i. However, even under random transmission, the matrix will usually be either symmetric or very nearly symmetric, with only small differences between identical pairs. If \code{output = 2}, the function returns a list with the following elements: (1) `mat_sim', the symmetric similarity matrix (see above); `mat_nobs', a matrix with the same dimensions as `mat_sim' which specifies the number of observations for each individual (by row, used to calculate `mat_pc', see later); `mat_pc', a non-symmetric matrix in which each cell represents the percent of observations of the individual in the i'th row that are shared with the individual in the j'th column; and `dat', a nested list, with one element for each individual which comprises a list of dataframes, one for each other individual, each one of which contains the subset of observations that are shared between the two individuals. Each dataframe contains the same columns as in the \code{acoustics_ls} dataframes with the following columns added: `pos_in_acc2', `timestamp_acc2', `receiver_lat_acc2' and `receiver_long_acc2', which represent the positions, time stamps and locations of corresponding observations in the second individual's dataframe to the first individual's dataframe, and `difftime_abs' and `dist_btw_receivers' which represent the duration (minutes) and distances (m) between corresponding observations. When there are no shared observations between a pair of individuals, the element simply contains \code{NULL}. Note that `mat_pc' is computed by (mat_sim/mat_nobs)*100. The matrix is therefore non-symmetric (if individuals have differing numbers of observations); i.e., mat_pc[i, j] is the percent of individual i's observations that are shared with individual j; while mat_pc[j, i] is the percent of individual j's observations that are shared with individual i. NaN elements are possible in `mat_pc' for levels of the factor `individual_id' without observations.
#' @param verbose A logical input which specifies whether or not to print messages to the console which relay function progress. This is ignored if \code{cl} is supplied.
#'
#' @details
#' \subsection{Background}{
#' Passive acoustic telemetry is widely used to study animal space use, and the possible drivers of spatiotemporal patterns in space use, in aquatic environments. Patterns in space use are widely related to environmental conditions, such as temperature, but the role of interactions among individuals is often more challenging to investigate due to a paucity of data, despite their likely importance. However, discrete detections also contain information on interactions among individuals that may influence space use through their similarities and differences among individuals over time and space.}
#'
#' \subsection{Implications}{
#' Similarities and differences can take different forms with differing ecological implications. For example, for individuals that are frequently detected in similar areas, detections may indicate (a) prolonged associations among individuals, if detections are usually closely associated in time and space (for example, due to parent-offspring relationships, group-living and/or mating); or (b) avoidance and/or territorial-like behaviour if detections, while close in space, are usually at different receivers and/or disjointed in time. Likewise, detection similarities among individuals that are rarely detected, or usually detected at disparate receivers, may reflect important interactions among those individuals at particular times (e.g. mating).}
#'
#'\subsection{Methods}{
#' To explore similarities and differences in patterns of space use, visualisation of detection histories with abacus plots and maps is beneficial. However, with many individuals and large receiver arrays, quantification of the similarities in detections over time and space is challenging. To this end, \code{\link[flapper]{make_matrix_cooccurence}} computes a similarity matrix across all individuals, defining the number (or percentage) of detections for each individual that are nearby in time, or space, to detections for each other individual.}
#'
#' \subsection{Assumptions}{
#' The distances beyond which detections of different individuals nearby in time are considered to demonstrate that those individuals are not closely associated are Euclidean. This may be problematic (e.g. when receivers hug complex coastlines).}
#'
#' @return The function returns a matrix or a list, depending on the input to \code{output} (see above).
#'
#' @examples
#' #### Prepare data
#' # acoustics_ls requires a dataframe with certain columns
#' # individual_id should be a factor
#' dat_acoustics$individual_id <- factor(dat_acoustics$individual_id)
#' # ensure dataframe ordered by individual, then time stamp
#' dat_acoustics <- dat_acoustics[order(dat_acoustics$individual_id, dat_acoustics$timestamp), ]
#' # define list of dataframes
#' acoustics_ls <- split(dat_acoustics, factor(dat_acoustics$individual_id))
#'
#' #### Example (1): Compute detection similarity matrix using default options
#' # mat_sim contains the number of observations shared among individuals
#' mat_sim <- make_matrix_cooccurence(acoustics_ls = acoustics_ls,
#'                            thresh_time = 90,
#'                            thresh_dist = 0)
#' prettyGraphics::pretty_mat(mat_sim, col_diag = "dimgrey")
#'
#' #### Example (2): Return list of outputs
#' out_ls <- make_matrix_cooccurence(acoustics_ls = acoustics_ls,
#'                           thresh_time = 90,
#'                           thresh_dist = 0,
#'                           output = 2)
#' names(out_ls)
#' # Examine number of observations for each individual
#' prettyGraphics::pretty_mat(out_ls$mat_nobs)
#' # Examine % shared detections between individuals
#' prettyGraphics::pretty_mat(out_ls$mat_pc, col_diag = "dimgrey")
#'
#' #### Example (3): Turn off messages with verbose = FALSE
#' out_ls_non_verb <- make_matrix_cooccurence(acoustics_ls = acoustics_ls,
#'                                    thresh_time = 90,
#'                                    thresh_dist = 0,
#'                                    verbose = FALSE)
#'
#' #### Example (4): Implement algorithm in parallel
#' out_ls_pl <- make_matrix_cooccurence(acoustics_ls = acoustics_ls,
#'                              thresh_time = 90,
#'                              thresh_dist = 0,
#'                              cl = parallel::makeCluster(2L),
#'                              output = 2)
#' names(out_ls_pl)
#'
#' @author Edward Lavender
#' @export

make_matrix_cooccurence <-
  function(acoustics_ls,
           thresh_time,
           thresh_dist,
           cl = NULL,
           varlist = NULL,
           output = 1,
           verbose = TRUE
  ){

    #### Initial checks
    # Check that acoustics dataframes contains required columns
    mapply(acoustics_ls, 1:length(acoustics_ls), FUN = function(acoustics, i){
      if(!is.null(acoustics)){
        if(any(!(c("individual_id", "timestamp", "receiver_long", "receiver_lat") %in% colnames(acoustics)))){
          stop(paste0("acoustic_ls[[", i, "]] does not contain all required column names."))
        }
      }
    })

    #### Set up
    # Check that acoustics_ls is a factor
    first_non_NULL <- min(which(sapply(acoustics_ls, function(acoustics) return(!is.null(acoustics)))))
    if(!inherits(acoustics_ls[[first_non_NULL]]$individual_id, "factor")) stop("acoustics_ls dataframe must contain individual_id column that is a factor.")
    # Define factor levels (i.e. the names of individuals )
    id_names <- levels(acoustics_ls[[first_non_NULL]]$individual_id)
    # id_names <- as.character(sapply(acoustics_ls, function(acoustics) return(acoustics$individual_id[1])))
    nid <- length(id_names)
    # Check that there is one element in acoustics_ls for every individual
    if(length(acoustics_ls) != nid) stop("'acoustics_ls' needs one element for every individual individual_id level. Add NULL elements to 'acoustics_ls' for remaining individuals.")
    # Create a blank similarity matrix which we'll fill in
    mat_sim <- matrix(NA, nrow = nid, ncol = nid, dimnames = list(id_names, id_names))

    #### Loop over all combinations of individuals, pair time series and identify
    # ... nearby observations in time and space
    if(all(!is.null(varlist), !is.null(cl))) parallel::clusterExport(cl = cl, varlist = varlist)
    lout <-
      pbapply::pblapply(acoustics_ls, cl = cl, function(acc1){

        #### Testing:
        # acc1 = acoustics_ls[[1]]; acc2 = acoustics_ls[[2]];
        # acc1 = acoustics_ls[[2]]; acc2 = acoustics_ls[[1]];

        #### For each individual, loop over each other individual...
        lint <-
          lapply(acoustics_ls, function(acc2){

            #### Ignore NULL/empty elements
            if(any(is.null(nrow(acc1)), is.null(nrow(acc2)), nrow(acc1) == 0, nrow(acc2) == 0)) return(NULL)

            if(acc1$individual_id[1] != acc2$individual_id[1]){

              #### Print individual
              if(verbose){
                cat("\n===================================================================================\n")
                cat(paste("Individual (", as.character(acc1$individual_id[1]),
                          ") and individual (", as.character(acc2$individual_id[1]), ").\n"))
              }

              #### Match time series based on closest observations in time using pair_ts()
              if(verbose) cat("Matching detection time series...\n")
              # Remove any observations from the second individual more than some limit outside of the time series of the first individual
              # ... and vice versa (for speed when matching).
              acc2 <- acc2[acc2$timestamp >= (min(acc1$timestamp) - thresh_time*60*2) &
                             acc2$timestamp <= (max(acc1$timestamp) + thresh_time*60*2), ]
              if(nrow(acc2) == 0) return(NULL)
              acc1 <- acc1[acc1$timestamp >= (min(acc2$timestamp) - thresh_time*60*2) &
                             acc1$timestamp <= (max(acc2$timestamp) + thresh_time*60*2), ]
              if(nrow(acc1) == 0) return(NULL)

              # Check for duplicated time stamps in each individual's dataframe, and adjust these
              # ... by a small fraction prior to matching, so that all are included. Unless there are 100,000s
              # ... of duplicate time stamps, this approach does not produce any duplicated observations.
              dup1 <- duplicated(acc1$timestamp)
              dup2 <- duplicated(acc2$timestamp)
              adj <- (thresh_time*60)/4
              if(any(dup1)){
                pos_dups1 <- which(dup1)
                lpd1 <- length(pos_dups1)
                adj_dups1 <- stats::runif(lpd1, -adj, adj)
                acc1$timestamp[pos_dups1] <- acc1$timestamp[pos_dups1] + adj_dups1
                if(any(duplicated(acc1$timestamp))) warning(paste("Duplicate time stamps in, ", as.character(acc1$individual_id[1]), "element in acoustic_ls."))
              }
              if(any(dup2)){
                pos_dups2 <- which(dup2)
                lpd2 <- length(pos_dups2)
                adj_dups2 <- stats::runif(lpd2, -adj, adj)
                acc2$timestamp[pos_dups2] <- acc2$timestamp[pos_dups2] + adj_dups2
                if(any(duplicated(acc2$timestamp))) warning(paste("Duplicate time stamps in, ", as.character(acc2$individual_id[1]), "element in acoustic_ls."))
              }
              # Match time series, readjusting any adjusted time stamps back to their original values
              # ... before these are added to the dataframe.
              acc1$pos_in_acc2 <- Tools4ETS::match_ts_nearest(acc1$timestamp, acc2$timestamp)
              if(any(dup1)) acc1$timestamp[pos_dups1] <- acc1$timestamp[pos_dups1] - adj_dups1
              if(any(dup2)) acc2$timestamp[pos_dups2] <- acc2$timestamp[pos_dups2] - adj_dups2
              acc1$timestamp_acc2 <- acc2$timestamp[acc1$pos_in_acc2]

              #### Exclude any time stamps more than time stamp beyond each other (could be 0 mins)
              # Implement this now, before a threshold based on distance, below, for speed.
              if(verbose) cat("Processing time series by theshold time difference...\n")
              acc1$difftime_abs <- abs(difftime(acc1$timestamp, acc1$timestamp_acc2, units = "mins"))
              acc1 <- acc1[acc1$difftime_abs <= thresh_time, ]
              if(nrow(acc1) == 0) return(NULL)

              #### Distances between pairs of receivers
              if(verbose) cat("Computing differences between pairs of receivers...\n")
              # Add receivers
              acc1$receiver_lat_acc2 <- acc2$receiver_lat[acc1$pos_in_acc2]
              acc1$receiver_long_acc2 <- acc2$receiver_long[acc1$pos_in_acc2]
              # Compute distances
              acc1$dist_btw_rec <- geosphere::distGeo(acc1[, c("receiver_long", "receiver_lat")],
                                                      acc1[, c("receiver_long_acc2", "receiver_lat_acc2")])

              #### Exclude any receivers more than some threshold distance beyond each other (could be 0 m):
              if(verbose) cat("Processing time series by threshold distance...\n")
              acc1 <- acc1[acc1$dist_btw_rec <= thresh_dist, ]
              if(nrow(acc1) == 0) return(NULL) else return(acc1)

            }
          })

        return(lint)

      })
    if(!is.null(cl)) parallel::stopCluster(cl = cl)

    #### Populate similarity matrix
    for(i in 1:nid){
      for(j in 1:nid){
        if(i != j){
          d <- lout[[i]][[j]]
          if(!is.null(d)) mat_sim[i, j] <- nrow(d) else mat_sim[i, j] <- 0
        }
      }
    }

    #### Matrix based on % similarity
    mat_nobs <- mat_sim[]
    for(i in 1:nrow(mat_nobs)){
      if(!is.null(acoustics_ls[[i]])) mat_nobs[i, ] <- nrow(acoustics_ls[[i]]) else mat_nobs[i, ] <- 0
    }
    mat_pc <- (mat_sim/mat_nobs)*100

    #### Return outputs
    if(!(output %in% 1:2)){
      warning(paste("'output ", output, " not supported; defaulting to output = 1."))
      output <- 1
    }
    if(output == 1) {
      out <- mat_sim
    } else if(output == 2) {
      out <- list(mat_sim = mat_sim, mat_nobs = mat_nobs, mat_pc = mat_pc, dat = lout)
    }
    return(out)

  }
