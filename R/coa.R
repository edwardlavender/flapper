#######################################
#######################################
#### coa_setup_delta_t()

#' @title Suggest time intervals over which to calculate centres of activity
#' @description This function implements the two approaches introduced by Simpfendorfer et al. (2002) to suggest suitable time interval(s) over which to estimate centres of activity (COAs). Both approaches rely on prior specification of a candidate time interval (\code{delta_t}) and provide a qualitative assessment of its suitability. The manual iteration of this scheme over multiple intervals provides a measure of their relative suitability.
#'
#' The first method (\code{method = 1L}) is based on the premise that estimated COAs are more accurate when they are estimated from detections at more receivers and thus examines the frequency distribution (over the whole time series) of the number of receivers at which detections were made within a given time interval (e.g., one hour). The second method (\code{method = 2L}) is based on the premise that estimated COAs are more accurate when they are estimated from a larger number of detections and thus examines the frequency distribution (over the whole time series) of the number of detections in a given time interval. In a comparison of multiple candidate intervals, the most appropriate interval is the one which optimises a trade off between sufficiently meeting these criteria while remaining sufficiently small in duration such that COAs remain a meaningful representation of the individual's location or short-term centre of activity.
#'
#' To implement these approaches, it is necessary to supply a dataframe with detections for a particular individual (\code{acoustics}), the interval to be evaluated (\code{delta_t}) and the \code{method}. For each method, two implementations are possible. The first implementation (\code{implementation = 1L}) follows Simpfendorfer et al.'s (2002) approach and examines the frequency distribution(s) for the total number of receivers and/or detections, which is appropriate if the total number of receivers is constant over the duration of detections. The second implementation (\code{implementation = 2L}) examines the distribution(s) of the percent of receivers with detections or the number of detections per receiver and can be more appropriate if the number of receivers changes substantially through time. Using these inputs, the function plots the specified distributions and returns a named list with the data plotted. This approach can be applied iteratively over multiple \code{delta_t} values to evaluate their relative suitability.
#'
#' @param acoustics A dataframe with passive acoustic telemetry detection time series (see \code{\link[flapper]{dat_acoustics}} for an example) for a specific individual. This should contain the following columns: a vector of receiver IDs, named 'receiver_id'; and a POSIXct vector of time stamps, named 'timestamp'.
#' @param delta_t A number or character, passed to the 'breaks' argument of \code{\link[base]{cut.POSIXt}}, that defines the time interval to be evaluated.
#' @param method An integer (\code{1L}, \code{2L} or \code{1:2L}) that defines the method(s) used to evaluate the time intervals (see Description).
#' @param implementation An integer (\code{1L} or \code{2L}) that defines the metric(s) used in the evaluation: \code{1}, the total number of receivers/detections; or \code{2}, the percent of receivers/number of detections per receiver (see Description). If \code{implementation = 2L} and the total number of receivers is constant through time, the function reverts to \code{implementation = 1L}, the results of which are more interpretable.
#' @param moorings (optional) If \code{implementation = 2L}, \code{moorings} is a dataframe that defines, for each receiver, the deployment period. This is required to calculate the number of operational receivers over time (via \code{\link[flapper]{get_n_operational_ts}}). This must contain the following columns: a vector of receiver IDs that includes all receivers in \code{acoustics}, named 'receiver_id'; and the times of receiver deployment and retrieval, named 'receiver_start_date' and 'receiver_end_date' respectively. These are coerced to POSIXct vectors, if required, to match \code{acoustics$timestamp}.
#' @param xlim,ylim X and y axis limits. If \code{method = 1:2L}, a single vector will affect both plots identically, whereas a list with one element for each method will affect plots differently.
#' @param xlab,ylab,main Character strings that define the x and y axis labels and the plot title. If \code{method = 1:2L}, a single input will label both plots identically, whereas a list with one element for each method will label the plots differently.
#' @param add_additional (optional) A stand-alone function, to be executed after a plot has been made, to customise the result.
#' @param type,... Additional arguments, passed to \code{\link[prettyGraphics]{pretty_plot}}, that affect all plots.
#' @param verbose A logical variable that defines whether or not to print messages to the console to relay function progress.
#'
#' @details The selection of the most appropriate time interval over which to calculate COAs reflects a trade off between (a) larger intervals, over which time more detections at more receivers are made, which enable more accurate estimates, versus (b) smaller intervals, over which time movement is more restricted and centres of activity are representative of the individual's location or short-term centre of activity (Simpfendorfer et al., 2002). As a starting point, Simpfendorfer et al. (2002) suggest intervals from 5 - 60 minutes may often be suitable, although longer intervals may be required.
#'
#' @return The function returns a named list and a plot for each \code{method}. The list contains a 'data' element that is a named list of dataframes for each method ('m1', 'm2') and implementation ('i1', 'i2'). For 'm1_i1', the dataframe contains the number of receivers with detections ('n_receiver_with_detections') and the percent of the time series for which that number of receivers made detections ('pc_of_ts'). For 'm1_i2', the percentage of receivers with detections ('pc_receiver_with_detections') is given instead. For 'm2_i1', the number of detections ('n_detections') is given; for 'm2_i2', it is the number of detections per receiver ('n_detections_per_receiver'). The list also contains an 'args' element that records the inputs to \code{acoustics} and \code{delta_t} for reference.
#'
#' @examples
#' #### Example (1): For a specified delta_t, use a specific method and implementation
#' pp <- graphics::par(mfrow = c(2, 2))
#' dat <- coa_setup_delta_t(acoustics = dat_acoustics[dat_acoustics$individual_id == 25, ],
#'                          delta_t = "6 hours",
#'                          method = 1L,
#'                          implementation = 1L,
#'                          main = "1(1)")
#' utils::str(dat)
#' dat <- coa_setup_delta_t(acoustics = dat_acoustics[dat_acoustics$individual_id == 25, ],
#'                          delta_t = "6 hours",
#'                          method = 1L,
#'                          implementation = 2L,
#'                          moorings = dat_moorings,
#'                          main = "1(2)")
#' utils::str(dat)
#' dat <- coa_setup_delta_t(acoustics = dat_acoustics[dat_acoustics$individual_id == 25, ],
#'                          delta_t = "6 hours",
#'                          method = 2L,
#'                          implementation = 1L,
#'                          main = "2(1)")
#' utils::str(dat)
#' dat <- coa_setup_delta_t(acoustics = dat_acoustics[dat_acoustics$individual_id == 25, ],
#'                          delta_t = "6 hours",
#'                          method = 2L,
#'                          implementation = 2L,
#'                          moorings = dat_moorings,
#'                          main = "2(2)")
#' utils::str(dat)
#' graphics::par(pp)
#'
#' #### Example (2) For a specified delta_t, use both methods
#' pp <- graphics::par(mfrow = c(1, 2))
#' dat <- coa_setup_delta_t(acoustics = dat_acoustics[dat_acoustics$individual_id == 25, ],
#'                          delta_t = "6 hours",
#'                          method = 1:2L,
#'                          implementation = 1L)
#' graphics::par(pp)
#'
#' #### Example (3) Plot customisation options
#' # xlim,ylim,xlab,ylab,main accept a vector that affects all plots
#' # ... or a list that affects each plot (see 'ylim' versus 'xlim' below).
#' # ... other arguments can be passed to prettyGraphics::pretty_plot() via...
#' pp <- graphics::par(mfrow = c(1, 2))
#' dat <- coa_setup_delta_t(acoustics = dat_acoustics[dat_acoustics$individual_id == 25, ],
#'                          delta_t = "6 hours",
#'                          xlim = list(c(0, 50), c(0, 450)), ylim = c(0, 100),
#'                          method = 1:2L,
#'                          implementation = 1L,
#'                          main = list("A", "B"),
#'                          col = "royalblue", lwd = 2)
#' graphics::par(pp)
#'
#' #### Example (4) Compare multiple delta_t values
#' pp <- graphics::par(mfrow = c(3, 2))
#' delta_t_opts <- c("6 hours", "12 hours", "24 hours")
#' lapply(delta_t_opts, function(delta_t){
#'   coa_setup_delta_t(acoustics = dat_acoustics[dat_acoustics$individual_id == 25, ],
#'                     delta_t = delta_t,
#'                     method = 1:2L,
#'                     implementation = 2L,
#'                     moorings = dat_moorings,
#'                     main = delta_t, col = "royalblue")
#' })
#' graphics::par(pp)
#'
#' @seealso \code{\link[flapper]{coa}} calculates centres of activity.
#' @references Simpfendorfer, C. A., M. R. Heupel, and R. E. Hueter. 2002. Estimation of short-term centers of activity from an array of omnidirectional hydrophones and its use in studying animal movements. Canadian Journal of Fisheries and Aquatic Sciences 59:23-32.
#'
#' @author Edward Lavender
#' @export
#'

coa_setup_delta_t <- function(acoustics,
                              delta_t,
                              method = 1L,
                              implementation = 1L,
                              moorings = NULL,
                              xlab = NULL, ylab = NULL, main = NULL,
                              xlim = NULL, ylim = NULL,
                              add_additional = NULL,
                              type = "b",...,
                              verbose = TRUE){

  #### Initiate function
  ## Function call
  t_onset <- Sys.time()
  cat_to_console <- function(..., show = verbose) if(show) cat(paste(..., "\n"))
  cat_to_console(paste0("flapper::coa_setup_delta_t() called (@ ", t_onset, ")..."))
  ## List to store outputs
  out <- list(data = list("m1_i1" = NULL, "m1_i2" = NULL, "m2_i1" = NULL, "m2_i2" = NULL),
              args = list(acoustics = acoustics, delta_t = delta_t)
  )
  ## Checks
  check_names(input = acoustics, req = c("timestamp", "receiver_id"))
  acoustics$timestamp <- check_tz(input = acoustics$timestamp)
  tz <- lubridate::tz(acoustics$timestamp)

  #### Define limits
  if(!inherits(xlim, "list")) xlim <- list(xlim, xlim)
  if(!inherits(ylim, "list")) ylim <- list(ylim, ylim)

  #### Define intuitive axis labels
  if(!inherits(main, "list")) main <- list(main, main)
  if(!is.null(xlab)){
    if(length(xlab) == 1L) xlab <- list(xlab, xlab)
  } else {
    xlab_ls <- list(method_1 = list(implementation_1 = "N (receivers)",
                                    implementation_2 = "% (receivers)"),
                    method_2 = list(implementation_1 = "N (detections)",
                                    implementation_2 = "N (detections per receiver)")
    )
    xlab <- list("", "")
    for(i in method) xlab[[i]] <- xlab_ls[[i]][[implementation]]
  }
  if(is.null(ylab)) ylab <- "Frequency (%)"
  if(length(ylab) == 1L) ylab <- list(ylab, ylab)

  #### Define a time series of receiver deployment periods
  if(implementation == 2L) {
    cat_to_console("... Calculating the number of operational receivers through time for implementation = 2L...")
    if(is.null(moorings)) stop("'moorings' is required for implementation 2L.")
    check_names(input = moorings, req = c("receiver_start_date", "receiver_end_date"),
                extract_names = colnames, type = all)
    # Ensure acoustics timestamps and moorings deployment periods are in the same format
    if(!inherits(moorings$receiver_start_date, "POSIXct")) {
      moorings$receiver_start_date <- as.POSIXct(moorings$receiver_start_date, tz = tz)
    }
    if(!inherits(moorings$receiver_end_date, "POSIXct")) {
      moorings$receiver_end_date   <- as.POSIXct(moorings$receiver_end_date, tz = tz)
    }
    # Get the number of receivers that were operational on each time step
    dat <- get_n_operational_ts(moorings,
                                start = "receiver_start_date", stop = "receiver_end_date",
                                times = seq(min(acoustics$timestamp), max(acoustics$timestamp), by = delta_t),
                                plot = FALSE
    )
    colnames(dat) <- c("interval", "n_receiver")
    # Revert to implementation = 1L if the number of receivers is constant (more interpretable)
    if(length(unique(dat$n)) == 1L) {
      message("The number of receivers deployed in each time interval is identical: reverting to implementation = 1L.")
      implementation <- 1L
    }
  } else {
    dat <- data.frame(interval = seq(min(acoustics$timestamp), max(acoustics$timestamp), by = delta_t))
  }

  #### Receiver statistics
  cat_to_console("... Calculating detection statistics across time intervals....")
  # Calculate the total number of detections at each receiver in each interval across the whole study
  acoustics$intervals <- cut(acoustics$timestamp, breaks = delta_t)
  n_intervals <- length(levels(acoustics$intervals))
  n_detections_per_receiver <- table(acoustics$intervals, acoustics$receiver_id)

  # Calculate the total number of receivers which received a detection in each X minute interval
  n_receiver_with_detections     <- rowSums(n_detections_per_receiver > 0)
  dat$n_receiver_with_detections <- as.numeric(n_receiver_with_detections)

  #### Method (1): The frequency distribution of the number (or %) of receivers with detections
  if(1L %in% method) {
    cat_to_console("... Implementing method 1...")

    #### Implementation (1)
    # We assume the number of receivers is constant through time
    # We calculate the frequency (%) distribution of the number of receivers with detections in each interval
    # ... (out of the total number of intervals)
    if(implementation == 1L) {
      cat_to_console("... ... Using implementation = 1L...")

      # Calculate the frequency distribution (%) for number of receivers which received a detection in each X minute interval
      # ... by taking the distribution of the number of receivers with detections divided by the number of time intervals (*100)
      dat_1 <- data.frame(table(n_receiver_with_detections)/n_intervals*100)
      colnames(dat_1) <- c("n_receiver_with_detections", "pc_of_ts")
      dat_1$n_receiver_with_detections <- as.numeric(as.character(dat_1$n_receiver_with_detections))

      # Define dataframe for plotting
      # Blank dataframe
      tmp <- data.frame(n_receiver_with_detections = 0:length(unique(acoustics$receiver_id)),
                        pc_of_ts = 0)
      # add calculated frequencies to this temporary dataframe
      tmp$pc_of_ts <- dat_1$pc_of_ts[match(tmp$n_receiver_with_detections, dat_1$n_receiver_with_detections)]
      # set NAs to 0
      tmp$pc_of_ts[is.na(tmp$pc_of_ts)] <- 0

      # Redefine dat_1 and add to list
      dat_1 <- tmp
      out$data$m1_i1 <- dat_1

      # Make plot
      prettyGraphics::pretty_plot(dat_1$n_receiver_with_detections, dat_1$pc_of_ts,
                                  xlim = xlim[[1]], ylim = ylim[[1]], main = main[[1]],
                                  xlab = xlab[[1]], ylab = ylab[[1]],
                                  type = type,...)

      #### Implementation (2)
      # ... The number of receivers varies through time
      # ... We calculate the frequency distribution of the percent of receivers that made a detection in each interval
    } else if (implementation == 2L){
      cat_to_console("... ... Using implementation = 2L...")

      # Calculate the percent of receivers which made a detection in each interval
      # ... (out of the total number of receivers that could have made a detection in that signal)
      dat$pc_receiver_with_detections <- (dat$n_receiver_with_detections/dat$n_receiver)*100

      # Get the overall frequency (%) distribution for the percent of receivers that made a detection by in each X minute interval
      # ... by taking the distribution of the percent of receivers with detections divided by the total number of time intervals (*100)
      dat_1 <- data.frame(table(dat$pc_receiver_with_detections)/n_intervals*100)
      colnames(dat_1) <- c("pc_receiver_with_detections", "pc_of_ts")
      dat_1$pc_receiver_with_detections <- as.numeric(as.character(dat_1$pc_receiver_with_detections))
      # Add to output list
      out$data$m1_i2 <- dat_1

      # Make plot
      prettyGraphics::pretty_plot(dat_1$pc_receiver_with_detections, dat_1$pc_of_ts,
                                  xlim = xlim[[1]], ylim = ylim[[1]], main = main[[1]],
                                  xlab = xlab[[1]], ylab = ylab[[1]],
                                  type = type,...)
    }

    #### Customise plot
    if(!is.null(add_additional)) add_additional()
  }

  #### Method (2): The frequency distribution of the number (or fraction per-receiver) of detections over the time series
  if(2L %in% method) {
    cat_to_console("... Implementing method 2...")

    # Calculate the total number of detections in each interval
    dat$n_detections <- rowSums(n_detections_per_receiver)

    #### Implementation (1)
    # ... The number of receivers is constant and we simply express the distribution for the total number of detections
    if(implementation == 1L) {
      cat_to_console("... ... Using implementation = 1L...")

      # Calculate the frequency (%) distribution of the total number of detections in each interval
      # ... (of the specific individual) over the whole time series
      dat_2 <- data.frame(table(dat$n_detections)/n_intervals * 100)
      colnames(dat_2) <- c("n_detections", "pc_of_ts")
      dat_2$n_detections <- as.numeric(as.character(dat_2$n_detections))
      # Add to list
      out$data$m2_i1 <- dat_2

      # Make plot
      prettyGraphics::pretty_plot(dat_2$n_detections, dat_2$pc_of_ts,
                                  xlim = xlim[[2]], ylim = ylim[[2]], main = main[[2]],
                                  xlab = xlab[[2]], ylab = ylab[[2]],
                                  type = type,...)

      #### Implementation (2)
      # ... The number of receivers varies and we express the distribution of the fraction of detections per receiver
    } else if(implementation == 2L){
      cat_to_console("... ... Using implementation = 2L...")

      # Calculate the fraction of detections in each interval per receiver
      dat$pc_detections_per_receiver <- dat$n_detections/dat$n_receiver

      # Calculate the frequency distribution (%) for the fraction of detections per receiver in each interval
      dat_2 <- data.frame(table(dat$pc_detections_per_receiver)/n_intervals * 100)
      colnames(dat_2) <- c("pc_detections_per_receiver", "pc_of_ts")
      dat_2$pc_detections_per_receiver <- as.numeric(as.character(dat_2$pc_detections_per_receiver))
      # Add to list
      out$data$m2_i2 <- dat_2

      # Make plot
      prettyGraphics::pretty_plot(dat_2$pc_detections_per_receiver, dat_2$pc_of_ts,
                                  xlim = xlim[[2]], ylim = ylim[[2]], main = main[[2]],
                                  xlab = xlab[[2]], ylab = ylab[[2]],
                                  type = type,...)
    }

    #### Customise plot
    if(!is.null(add_additional)) add_additional()
  }

  #### Return outputs
  t_end <- Sys.time()
  duration <- difftime(t_end, t_onset, units = "mins")
  cat_to_console(paste0("... flapper::coa_setup_delta_t() call completed (@ ", t_end, ") after ~", round(duration, digits = 2), " minutes."))
  return(out)

}


#######################################
#######################################
#### coa()

#' @title The centres of activity (COA) algorithm
#' @description This function implements the arithmetic mean-position algorithm to calculate an individual's centres of activity (COAs) though time from detections at passive acoustic telemetry receivers, as described by Simpfendorfer et al (2002). Under this approach, each COA is calculated as the arithmetic mean of the locations of receivers at which an individual was detected over a specified time interval, weighted by the frequency of detections at each of those receivers. To implement the function, a detection matrix that defines the number of detections of an individual along sequence of time steps at each receiver (e.g., from \code{\link[flapper]{make_matrix_detections}}) needs to be supplied, along a matrix of receiver locations using a planar coordinate system. For each time interval, the function calculates the centre of activity and returns a matrix or dataframe with this information.
#'
#' @param mat A detection matrix, with one row for each time step and one column for each receiver, in which each cell defines the number of detections at each time step/receiver (for a particular individual) (see \code{\link[flapper]{make_matrix_detections}}). It is advisable that the rows and columns of this matrix are labelled by time stamp and receiver respectively, especially if there are any rows without detections (see \code{output} below).
#' @param xy A matrix that defines receiver locations (x, y). This should contain one row for each receiver (column) in \code{mat} (in the same order as in \code{mat}) and two columns for the coordinates. Planar coordinates (i.e., Universal Transverse Mercator) projection are required for the averaging process.
#' @param cl A cluster object created by \code{\link[parallel]{makeCluster}} to implement the algorithm in parallel. If supplied, the connection to the cluster is closed within the function.
#' @param na_omit,as_POSIXct Processing options. \code{na_omit} is a logical variable that defines whether or not to omit NAs (i.e., rows in \code{mat} for which no detections were made and thus for which COAs cannot be calculated) from returned coordinates. If \code{output = "data.frame"} (see below), \code{as_POSIXct} is a function to convert timestamps, taken from the row names of \code{mat}, to a POSIXct vector. \code{as_POSIXct = NULL} suppresses this conversion.
#' @param output A character that defines the output format. Currently supported options are: \code{"matrix"}, which returns a matrix of the coordinates of COAs; and \code{"data.frame"}, which returns a dataframe with timestamps (taken from the row names of \code{mat}) and COA coordinates.
#' @param ... Additional arguments (none implemented).
#'
#' @details Centres of activity (COA) are a widely used metric for the reconstruction of patterns of space use from passive acoustic telemetry detections. Several methods have been developed to calculate COAs, but the mean-position algorithm is the commonest. Under this approach, COAs are estimated as an average of the locations of receivers at which an individual is detected over a specified time interval, weighted by the frequency of detections at each of those receivers. Within \code{\link[flapper]{flapper}}, COAs are calculated in three stages by first exploring possible time intervals over which to calculate COAs with \code{\link[flapper]{coa_setup_delta_t}}; then summarising detections over those intervals with \code{\link[flapper]{make_matrix_detections}}; and finally passing the resultant detection matrix to the \code{\link[flapper]{coa}} function to calculate COAs. This implements the arithmetic version of the mean-position algorithm, calculating the arithmetic mean of the receiver locations, weighted by the frequency of detections at each receiver.
#'
#' To generate estimates of space use, COAs are usually taken as point estimates from which utilisation distributions (typically kernel utilisation distributions, KUDs) are estimated. Thus, in the case of a coupled COA-KUD approach, usually the estimate of space use is a (kernel) utilisation distribution which describes the probability of relocating an individual in any given area at a randomly chosen time. Alternative methods of home range analysis, including those which incorporate time, such as dynamic Brownian bridge movement models, can be used to estimate the utilisation distribution. Generally, the COA approach is most suitable when detections are relatively frequent, and receivers are regularly distributed across an area. Under other conditions, its performance as a method for estimating space use has been subject to relatively limited evaluation but it can be problematic (e.g., in clustered arrays).
#'
#' @return The function returns a matrix or a dataframe, depending on the \code{output} argument, that represents a time series of COAs. If \code{na_omit = TRUE}, the time series may have 'gaps' over which COAs could not be calculated due to the absence of detections.
#'
#' @examples
#' #### Define data for the calculation of COAs
#'
#' ## (1) Define the period over which to calculate COAs
#' # ... by focusing on the range of time over which IDs were at liberty.
#' dat_ids$tag_start_date <- as.POSIXct(dat_ids$tag_start_date)
#' dat_ids$tag_end_date   <- as.POSIXct("2017-06-02")
#'
#' ## (2) Define the receivers over which to calculate detections
#' # ... Here we use factors to ensure that the order of receiver coordinates
#' # ... (see below) and the order or receivers in the detection matrix
#' # ... matches.
#' dat_moorings$receiver_id <- factor(dat_moorings$receiver_id)
#' dat_acoustics$receiver_id <- factor(dat_acoustics$receiver_id,
#'                                     levels = levels(dat_moorings$receiver_id))
#'
#' ## (3) Define the detection matrix
#' # ... Here we simply create a detection matrix across all IDs and using
#' # ... a convenient (but unjustified) delta_t value. In reality, we would need
#' # ... to consider the time series over which it is appropriate to calculate
#' # ... COAs and the appropriate delta t value(s) more carefully.
#' detection_matrix_by_id <- make_matrix_detections(dat_acoustics,
#'                                                  delta_t = "days",
#'                                                  start = min(dat_ids$tag_start_date),
#'                                                  end = max(dat_ids$tag_end_date))
#'
#' ## (4) Define receiver coordinates with UTM projection
#' proj_wgs84 <- sp::CRS("+init=epsg:4326")
#' proj_utm <- sp::CRS(paste("+proj=utm +zone=29 +datum=WGS84",
#'                           "+units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))
#' xy <- sp::SpatialPoints(dat_moorings[, c("receiver_long", "receiver_lat")],
#'                         proj_wgs84)
#' xy <- sp::spTransform(xy, proj_utm)
#' xy <- sp::coordinates(xy)
#'
#' #### Example (1): Implement the COA algorithm for an example individual
#' coa_mat <- coa(mat = detection_matrix_by_id[[1]], xy = xy)
#' utils::str(coa_mat)
#'
#' #### Example (2): Change the output format and coerce timestamps to output format
#' coa_dat <- coa(mat = detection_matrix_by_id[[1]], xy = xy, output = "data.frame")
#' utils::str(coa_dat)
#'
#' #### Example (3): Implement the algorithm on a cluster
#' # This will only be faster for very large detection time series.
#' \dontrun{
#' coa_mat <- coa(mat = detection_matrix_by_id[[1]],
#'                xy = xy,
#'                cl = parallel::makeCluster(2L))
#' }
#'
#' @seealso \code{\link[flapper]{coa_setup_delta_t}} suggests suitable time intervals over which to calculate COAs. \code{\link[flapper]{make_matrix_detections}} makes the detection matrices from detection time series data required by this function. For data in the VEMCO Vue export format, the 'COA' function in the VTrack package (https://github.com/RossDwyer/VTrack) can also be used to calculate centres of activity.
#'
#' @references Simpfendorfer, C. A., M. R. Heupel, and R. E. Hueter. 2002. Estimation of short-term centers of activity from an array of omnidirectional hydrophones and its use in studying animal movements. Canadian Journal of Fisheries and Aquatic Sciences 59:23-32.
#'
#' @author Edward Lavender
#' @export
#'

coa <- function(mat, xy, cl = NULL, na_omit = TRUE, as_POSIXct = as.POSIXct, output = "matrix",...){

  #### Define objects
  if(all(mat == 0)) stop("No detection(s) identified in mat: unable to calculate COA(s).")
  output <- check_value(input = output, supp = c("matrix", "data.frame"))
  if(is.null(rownames(mat))) {
    rownames(mat) <- 1:nrow(mat)
    if(output == "data.frame" & !is.null(as_POSIXct)) {
      message("'output = 'data.frame' and 'as_POSIXct' has been defined but rownames(mat) is NULL: ignoring input to 'as_POSIXct' argument.")
      as_POSIXct <- NULL
    }
  }
  if(is.null(colnames(mat))) colnames(mat) <- 1:ncol(mat)
  det_mat <- mat
  tmp <- matrix(c(NA, NA), nrow = 1, ncol = 2)
  rx <- xy[, 1]
  ry <- xy[, 2]

  #### Calculate COAs
  coa_ls <- pbapply::pblapply(1:nrow(det_mat), cl = cl, function(i){
    # Extract parameters for speed
    det_mat_i     <- det_mat[i, ]
    det_mat_i_sum <- sum(det_mat_i)
    if(det_mat_i_sum != 0) {
      # Calculate the average x location
      # ... as an average of the receiver locations, weighted by the number of detections at each receiver
      tmp[1, 1] <- sum(rx * det_mat_i)/det_mat_i_sum
      # Calculate the average y location
      tmp[1, 2] <- sum(ry * det_mat_i)/det_mat_i_sum
    } else {
      # If there are no detections in that interval, define NA (otherwise it will appear as 0, 0):
      tmp[ , ] <- NA
    }
    return(tmp)
  })
  if(!is.null(cl)) parallel::stopCluster(cl = cl)
  coa_mat <- do.call(rbind, coa_ls)
  rownames(coa_mat) <- rownames(det_mat)
  colnames(coa_mat) <- c("x", "y")

  #### Process COAs
  out <- coa_mat
  if(na_omit) out <- stats::na.omit(out)
  if(output == "data.frame") {
    out <- data.frame(timestamp = rownames(out), x = out[, 1], y = out[, 2])
    if(!is.null(as_POSIXct)) out$timestamp <- as_POSIXct(out$timestamp)
  }

  #### Return COAs
  return(out)
}
