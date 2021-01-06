######################################
######################################
#### get_detection_pr()

#' @title A detection probability function based on distance
#' @description This function calculates detection probability (e.g., of an acoustic detection) at specified distances from the sampling device (e.g., a passive acoustic telemetry receiver) using user-defined parameters (a model intercept, an coefficient for the effect of distance and an inverse link function). The function returns a plot of detection probability with distance and/or a vector of detection probabilities.
#'
#' @param distance A numeric vector of distances at which to calculate detection probability.
#' @param beta_0,beta_1 Single numbers that define the model coefficients (i.e., the intercept and gradient on the scale of the link function).
#' @param inv_link A function that defines the inverse link function.The default function is the logistic (inverse logit) function.
#' @param output An integer (\code{1L}, \code{2L} or \code{3L}) that defines the output type. \code{1L} returns a plot of detection probability against distance; \code{2L} returns a numeric vector of detection probabilities; and \code{3L} returns both of the above.
#' @param ... Additional arguments passed to \code{\link[prettyGraphics]{pretty_plot}} to customise the plot. These are only implemented if \code{output = 1L} or \code{output = 3L}.
#'
#' @return The function calculates detection probability at each specified distance and returns a plot, a vector of detection probabilities, or both, depending on the value of the \code{output} argument. If a vector of detection probabilities is returned, this contains the following attributes: X', the model matrix; 'beta', the regression coefficients; and 'inv_link', the inverse link function.
#'
#' @examples
#' #### Example (1): Implement the function using the default parameters
#' # The function returns a graph and a vector of detection probabilities
#' det_pr <- get_detection_pr()
#' utils::head(det_pr)
#' # The vector has attributes:
#' # ... 'X' (the model matrix),
#' # ... 'beta' (the regression coefficient)
#' # ... 'inv_link' (the inverse link function)
#' utils::str(det_pr)
#'
#' #### Example (2): Adjust model parameters
#' # Change regression coefficients
#' det_pr <- get_detection_pr(beta_0 = 2.5, beta_1 = -0.006)
#' # Use inverse probit
#' det_pr <- get_detection_pr(beta_0 = 2.5, beta_1 = -0.006, inv_link = stats::pnorm)
#'
#' #### Example (3): Modify graphical properties
#' det_pr <- get_detection_pr(beta_0 = 2.5,
#'                           beta_1 = -0.006,
#'                           type = "l",
#'                           xlab = "Distance (m)",
#'                           ylab = "Detection Probability")
#'
#' #### Example (4): Modify return options
#' # Only graph
#' get_detection_pr(output = 1L)
#' # Only values
#' get_detection_pr(output = 2L)
#' # Both graph and values (the default)
#' get_detection_pr(output = 3L)
#'
#' @author Edward Lavender
#' @export
#'

get_detection_pr <- function(distance = 1:1000,
                             beta_0 = 2.5,
                             beta_1 = -0.01,
                             inv_link = stats::plogis,
                             output = 3L,...){
  #### Checks
  stopifnot(length(beta_0) == 1 & length(beta_1) == 1)
  output <- check_value(input = output, supp = 1:3, warn = TRUE, default = 3L)

  #### Calculate detection probabilities
  X <- matrix(c(rep(1, length(distance)), distance), ncol = 2, byrow = FALSE)
  beta <- matrix(c(beta_0, beta_1), nrow = 2)
  rownames(beta) <- c("beta_0", "beta_1")
  Y <- inv_link(X %*% beta)
  Y <- as.numeric(Y)
  # Add attributes
  attributes(Y)$X        <- X
  attributes(Y)$beta     <- beta
  attributes(Y)$inv_link <- inv_link

  #### Visualise detection probabilities
  if(output %in% c(1, 3)){
    prettyGraphics::pretty_plot(X[, 2], Y,...)
  }

  #### Return detection probabilities
  if(output %in% c(2, 3)){
    return(Y)
  } else return(invisible())
}


######################################
######################################
#### get_detection_centroids()

#' @title Define detection centroids around receivers
#' @description This function defines the areas surveyed by receivers (termed 'detection centroids') as a spatial object, based on an estimate of the detection range (m) and any barriers to detection. To implement the function, receiver locations must be supplied as a SpatialPoints or SpatialPointsDataFrame object with the Universe Transverse Mercator coordinate reference system. The function defines a spatial buffer around each receiver according to the estimated detection range, cuts out any barriers to detection, such as the coastline, and returns a SpatialPolygons object that defines the combined detection centroid across all receivers or receiver-specific detection centroids.
#'
#' @param xy A \code{\link[sp]{SpatialPoints-class}} or \code{\link[sp]{SpatialPointsDataFrame-class}} object that defines receiver locations. The coordinate reference system should be the Universe Transverse Mercator coordinate reference system.
#' @param detection_range A number that defines the detection range (m) of receivers.
#' @param coastline (optional) A \code{\link[sp]{SpatialPolygonsDataFrame-class}} object that defines barriers (such as the coastline) that block receivers from surveying areas within their detection range.
#' @param plot A logical input that defines whether or not to plot receivers, their centroids, and the buffer (if specified).
#' @param ... Additional arguments passed to \code{\link[rgeos]{gBuffer}}, such as \code{byid} and/or \code{quadsegs}.
#'
#' @return The function returns a \code{\link[sp]{SpatialPolygons-class}} object of the detection centroids around receivers that represents the area they survey under the assumption of a constant detection range, accounting for any barriers to detection. By default, this will contain a single feature, which is suitable for the calculation of the total area surveyed by receivers (see \code{\link[flapper]{get_detection_area_sum}}) because it accounts for the overlap in the detection ranges of receivers. However, if \code{byid = TRUE} is passed via \code{...} to \code{\link[rgeos]{gBuffer}}, the returned object will have a feature for each pair of coordinates in (\code{xy}) (i.e., receiver). This is less appropriate for calculating the area surveyed by receivers, since areas surveyed by multiple receivers will be over-counted, but it is suitable when the centroids for particular receivers are required (e.g., to extract environmental conditions within a specific receiver's detection ranges) (see \code{\link[flapper]{get_detection_centroids_envir}}.
#'
#' @examples
#' #### Define receiver locations as a SpatialPoints object with a UTM CRS
#' proj_wgs84 <- sp::CRS("+init=epsg:4326")
#' proj_utm <- sp::CRS(paste("+proj=utm +zone=29 +datum=WGS84",
#'                           "+units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))
#' xy <- sp::SpatialPoints(dat_moorings[, c("receiver_long", "receiver_lat")],
#'                         proj_wgs84)
#' xy <- sp::spTransform(xy, proj_utm)
#'
#' #### Example (1): Calculate the total area sampled by receivers
#' get_detection_centroids(xy)
#'
#' #### Example (2): Account for barriers in the study area
#' get_detection_centroids(xy, coastline = dat_coast)
#'
#' #### Example (3): Adjust the detection range
#' get_detection_centroids(xy, detection_range = 400, coastline = dat_coast)
#' get_detection_centroids(xy, detection_range = 500, coastline = dat_coast)
#'
#' #### Example (4): Suppress the plot
#' get_detection_centroids(xy, coastline = dat_coast, plot = FALSE)
#'
#' #### Example (5): Output characteristics are controlled via byid
#' # A SpatialPolygons with one feature is the implicit output
#' sp_1 <- get_detection_centroids(xy, coastline = dat_coast, byid = FALSE)
#' sp_1
#' # An SpatialPolygons with one feature for each element in xy
#' # ... can be returned via byid = TRUE
#' sp_2 <- get_detection_centroids(xy, coastline = dat_coast, byid = TRUE)
#' sp_2
#' # The total area of the former will be smaller, since areas covered
#' # ... by multiple receivers are merged
#' rgeos::gArea(sp_1); rgeos::gArea(sp_2)
#' # But it can be more convenient to use the latter format in some cases
#' # ... because it is easy to isolate specific centroids
#' raster::plot(dat_coast)
#' raster::plot(sp_1[1], add = TRUE, col = "red")  # single feature
#' raster::plot(sp_2[1], add = TRUE, col = "blue") # isolate specific features
#'
#' @author Edward Lavender
#' @export

get_detection_centroids <- function(xy,
                                detection_range = 425,
                                coastline = NULL,
                                plot = TRUE,...){

  #### Checks
  # Check xy is a SpatialPoints object or similar
  check_class(input = xy, to_class = c("SpatialPoints", "SpatialPointsDataFrame"), type = "stop")

  #### Define buffers around receivers equal to detection radius
  xy_buf <- rgeos::gBuffer(xy, width = detection_range,...)

  #### Clip around coastline (if applicable)
  if(!is.null(coastline)) {
    if(length(xy_buf) == 1) {
      xy_buf <- rgeos::gDifference(xy_buf, coastline)
    } else {
      xy_buf <- lapply(1:length(xy_buf), function(i) rgeos::gDifference(xy_buf[i], coastline))
      xy_buf <- do.call(raster::bind, xy_buf)
    }
  }

  #### Plot [update to use pretty_map()]
  if(plot){
    if(!is.null(coastline)) {
      raster::plot(coastline)
      graphics::points(xy, pch = 4)
    } else {
      raster::plot(xy, pch = 4)
    }
    raster::lines(xy_buf)
  }

  #### Return outputs
  return(xy_buf)

}


######################################
######################################
#### get_detection_area_sum()

#' @title Calculate the total area sampled by acoustic receivers
#' @description This function calculates the total area sampled by receivers, under the assumption of a constant detection range. To implement the function, receiver locations must be supplied as a SpatialPoints or SpatialPointsDataFrame object with the Universe Transverse Mercator coordinate reference system. The \code{\link[flapper]{get_detection_centroids}} is used to calculate the detection centroids around receivers, given a specified detection range (m) and any barriers to detection, such as coastline, and then the total area covered by receivers is calculated, accounting for overlapping centroids.
#'
#' @param xy,detection_range,coastline,plot,... Arguments required to calculate and visualise detection centroids via \code{\link[flapper]{get_detection_centroids}}; namely, receiver locations (\code{xy}), the detection range (\code{range}), barriers to detection (\code{coastline}), and whether or not to plot the centroids (\code{plot}).
#' @param scale A number that scales the total area (m). The default (\code{1/(1000^2)}) converts the units of \eqn{m^2} to \eqn{km^2}.
#'
#' @details This is a simple metric of the overall receiver sampling effort. This may be a poor metric if the assumption of a single detection range across all receivers is substantially incorrect or if there are substantial changes in the receiver array over the course of a study.
#'
#' @return The function returns a number that represents the total area surveyed by receivers (by default in \eqn{km^2}) and, if \code{plot = TRUE}, a plot of the area with receivers and their detection ranges.
#'
#' @seealso \code{\link[flapper]{get_detection_centroids}} defines detection centroids, across which the detection area is calculated. \code{\link[flapper]{get_detection_area_ts}} calculates the area sampled by receivers through time.
#'
#' @examples
#' #### Define receiver locations as a SpatialPoints object with a UTM CRS
#' proj_wgs84 <- sp::CRS("+init=epsg:4326")
#' proj_utm <- sp::CRS(paste("+proj=utm +zone=29 +datum=WGS84",
#'                           "+units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))
#' xy <- sp::SpatialPoints(dat_moorings[, c("receiver_long", "receiver_lat")],
#'                         proj_wgs84)
#' xy <- sp::spTransform(xy, proj_utm)
#'
#' #### Example (1): Calculate the total area sampled by receivers
#' get_detection_area_sum(xy)
#'
#' #### Example (2): Account for barriers in the study area
#' get_detection_area_sum(xy, coastline = dat_coast)
#'
#' #### Example (3): Adjust the detection range
#' get_detection_area_sum(xy, detection_range = 400, coastline = dat_coast)
#' get_detection_area_sum(xy, detection_range = 500, coastline = dat_coast)
#'
#' #### Example (4): Adjust the units
#' get_detection_area_sum(xy, coastline = dat_coast, scale = 1) # m2
#'
#' #### Example (5): Suppress the plot
#' get_detection_area_sum(xy, coastline = dat_coast, plot = FALSE)
#'
#' @author Edward Lavender
#' @export
#'

get_detection_area_sum <- function(xy,
                                   detection_range = 425,
                                   coastline = NULL,
                                   scale = 1/(1000^2),
                                   plot = TRUE,...){

  #### Checks
  # If xy is empty, return area = 0
  if(length(xy) == 0) return(0)

  #### Define detection centroids
  xy_buf <- get_detection_centroids(xy = xy,
                                    detection_range = detection_range,
                                    coastline = coastline,
                                    plot = plot,...)

  #### Calculate area (m2) covered by buffers overall
  xy_area <- rgeos::gArea(xy_buf) * scale
  return(xy_area)
}


######################################
######################################
#### get_detection_area_ts()

#' @title Calculate the area sampled by receivers through time
#' @description This function extends \code{\link[flapper]{get_detection_area_sum}} to calculate how the total area sampled by receivers changes through time.
#'
#' @param xy,detection_range,coastline,scale Arguments required to calculate the total area surveyed by receivers (at each time point) via \code{\link[flapper]{get_detection_area_sum}}.
#' @param plot A logical input that defines whether or not to plot a time series of the total area sampled by receivers.
#' @param verbose A logical input that defines whether or not to print messages to the console to relay function progress.
#' @param cl (optional) A cluster object created by \code{\link[parallel]{makeCluster}} to implement the algorithm in parallel. The connection to the cluster is closed within the function.
#' @param varlist (optional) A character vector of names of objects to export, to be passed to the \code{varlist} argument of \code{\link[parallel]{clusterExport}}. This may be required if \code{cl} is supplied. Exported objects must be located in the global environment.
#' @param ... Additional arguments, passed to \code{\link[prettyGraphics]{pretty_plot}}, to customise the plot produced.
#'
#' @return The function returns a dataframe with, for each date ('date') from the time of the first receiver's deployment to the time of the last receiver's retrieval, the number of receivers operational on that date ('n') and the total area sampled ('receiver_area'). If \code{plot = TRUE}, the function also returns a plot of the area sampled by receivers through time.
#'
#' @examples
#' #### Define SpatialPointsDataFrame with receiver locations and deployment dates
#' proj_wgs84 <- sp::CRS("+init=epsg:4326")
#' proj_utm <- sp::CRS(paste("+proj=utm +zone=29 +datum=WGS84",
#'                           "+units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))
#' xy <- sp::SpatialPoints(dat_moorings[, c("receiver_long", "receiver_lat")],
#'                         proj_wgs84)
#' xy <- sp::spTransform(xy, proj_utm)
#' xy <- sp::SpatialPointsDataFrame(xy, data = dat_moorings)
#'
#' #### Example (1): Implement function with default arguments
#' dat <- get_detection_area_ts(xy)
#'
#' #### Example (2): Adjust detection range, include coastline and use parallel processing
#' # For areas with complex coastline, this will reduce the speed of the algorithm
#' # So we will also supply a cluster to improve the computation time.
#' \dontrun{
#' dat <- get_detection_area_ts(xy,
#'                        detection_range = 500,
#'                        coastline = dat_coast,
#'                        cl = parallel::makeCluster(2L),
#'                        varlist = "dat_coast"
#'                        )
#' }
#'
#' #### Example (3) Hide or customise the plot
#' dat <- get_detection_area_ts(xy, plot = FALSE)
#' dat <- get_detection_area_ts(xy,
#'                        pretty_axis_args = list(axis = list(list(format = "%b-%y"),
#'                                                            list())),
#'                        xlab = "Time (month-year)",
#'                        ylab = expression(paste("Area (", m^2, ")")),
#'                        type = "l")
#'
#' @author Edward Lavender
#' @export
#'

get_detection_area_ts <- function(xy,
                              detection_range = 425,
                              coastline = NULL,
                              scale = 1/(1000^2),
                              plot = TRUE,
                              verbose = TRUE,
                              cl = NULL,
                              varlist = NULL,...){

  #### Checks
  t_onset <- Sys.time()
  cat_to_console <- function(..., show = verbose) if(show) cat(paste(..., "\n"))
  cat_to_console(paste0("flapper::get_detection_area_ts() called (@ ", t_onset, ")..."))
  cat_to_console("... Implementing function checks...")
  check_class(input = xy, to_class = "SpatialPointsDataFrame", type = "stop")
  check_names(input = xy, req = c("receiver_start_date", "receiver_end_date"))
  if(is.null(cl) & !is.null(varlist)) message("cl = NULL: input to 'varlist' ignored.")

  #### Implement algorithm
  cat_to_console("... Implementing algorithm...")
  rdate_seq <- seq.Date(min(xy$receiver_start_date), max(xy$receiver_end_date), 1)
  if(!is.null(cl) & !is.null(varlist)) parallel::clusterExport(cl = cl, varlist = varlist)
  rcov <- pbapply::pblapply(rdate_seq, cl = cl, function(rdate){
    pos <- which(xy$receiver_start_date <= rdate & xy$receiver_end_date >= rdate)
    receiver_area <- get_detection_area_sum(xy = xy[pos, ],
                                        detection_range = detection_range,
                                        coastline = coastline,
                                        scale = scale,
                                        plot = FALSE)
    d <- data.frame(date = rdate, n = length(pos), receiver_area = receiver_area)
    return(d)
  })
  if(!is.null(cl)) parallel::stopCluster(cl)
  rcov <- dplyr::bind_rows(rcov)


  #### Visualise time series
  cat_to_console("... Visualising time series")
  if(plot) prettyGraphics::pretty_plot(rcov$date, rcov$receiver_area,...)

  #### Return outputs
  t_end <- Sys.time()
  duration <- difftime(t_end, t_onset, units = "mins")
  cat_to_console(paste0("... flapper::get_detection_area_ts() call completed (@ ", t_end, ") after ~", round(duration, digits = 2), " minutes."))
  return(rcov)

}


######################################
######################################
#### get_n_operational_ts()

#' @title Calculate the number of operational units through time
#' @importFrom lubridate `%within%`
#'
#' @description This function calculates the number of operational units through time (e.g., the number of individuals at liberty or the number of active acoustic receivers over the course of a study). To implement the function, a dataframe (\code{data}) must be supplied that defines the start and end time of each unit's operational period. Then, for each time step in a user-specified sequence of times, or an automatically generated sequence of times from the earliest start time to the latest end time in \code{data}, the function counts the number of units that were operational on each time step and returns a dataframe (and, if specified, a plot) with this information.
#'
#' @param data A dataframe of observations. At a minimum, this should contain two columns that define the start and end of each unit's operational period, identified by \code{start} and \code{stop} below.
#' @param start A character string that defines the name of the column in \code{data} that defines the start time of each unit's operational period.
#' @param stop A character string that defines the name of the column in \code{data} that defines the end time of each unit's operational period.
#' @param times A vector of times for which to calculate the number of units that were operational at each time. If \code{times = NULL}, a regular sequence of dates from the first to the last date in \code{data} is used.
#' @param plot A logical variable that defines whether or not to plot the time series of the number of operational units.
#' @param ... Additional arguments, passed to \code{\link[prettyGraphics]{pretty_plot}}, to customise the plot.
#'
#' @details This is a simple metric of sampling effort. For acoustic receivers, \code{\link[flapper]{get_detection_area_ts}} provides another metric of sampling effort.
#'
#' @return The function returns a dataframe that, for each time step ('time'), defines the number of operational units at that time ('n'). If \code{plot = TRUE}, the function also plots a time series of the number of operational units.
#'
#' @examples
#' #### Example (1): Number of operational receivers over an acoustic telemetry study
#' dat_n <- get_n_operational_ts(data = dat_moorings,
#'                           start = "receiver_start_date",
#'                           stop = "receiver_end_date")
#' utils::head(dat_n)
#'
#' #### Example (2): Number of individuals at liberty over a tagging study
#' # Define 'tag_end_date' as hypothetical end date of a study
#' # ...  and assume that all individuals remained tagged until this time
#' dat_ids$tag_end_date <- as.Date("2017-06-05")
#' dat_n <- get_n_operational_ts(data = dat_ids,
#'                           start = "tag_start_date",
#'                           stop = "tag_end_date")
#'
#' #### Example (3): Specify the time period under consideration
#' dat_n <- get_n_operational_ts(data = dat_ids,
#'                           start = "tag_start_date",
#'                           stop = "tag_end_date",
#'                           times = seq(min(dat_moorings$receiver_start_date),
#'                                       max(dat_moorings$receiver_end_date), 1)
#'                           )
#'
#' #### Example (4): Suppress or customise the plot
#' dat_n <- get_n_operational_ts(data = dat_ids,
#'                           start = "tag_start_date",
#'                           stop = "tag_end_date",
#'                           plot = FALSE)
#' dat_n <- get_n_operational_ts(data = dat_ids,
#'                           start = "tag_start_date",
#'                           stop = "tag_end_date",
#'                           xlab = "Time", ylab = "N (individuals)",
#'                           type = "l")
#'
#' @author Edward Lavender
#' @export
#'

get_n_operational_ts <- function(data, start, stop, times = NULL, plot = TRUE,...){

  #### Define a sequence of times that span the range of the data, if required.
  if(is.null(times)){
    times <- seq(min(data[, start], na.rm = TRUE), max(data[, stop], na.rm = TRUE), by = "days")
  }

  #### Count the number of operational units at each time step
  count <- data.frame(time = times)
  count$n <- 0
  data$interval <- lubridate::interval(data[, start], data[, stop])
  for(i in 1:nrow(count)){
    count$n[i] <- sum(count$time[i] %within% data$interval)
  }

  #### Plot the results
  if(plot) prettyGraphics::pretty_plot(count$time, count$n,...)

  #### Return the dataframe
  return(count)
}


######################################
######################################
#### get_id_rec_overlap()

#' @title Calculate the overlap between individuals' time at liberty and receivers' operational periods
#' @description This function calculates the duration of the overlap (in days) between individuals' time at liberty and receivers' operational periods. To implement this function, a dataframe with individual deployment periods and another with receiver deployment periods must be specified. The duration of the overlap between these intervals can be calculated for all combinations of individuals and receivers within these two dataframes, for all combinations of specified individuals and receivers, or for specific individual/receiver pairs. The function returns a dataframe of the overlap duration for these individual/receiver combinations or a vector of values that is matched against another dataframe.
#' @param ids A dataframe that defines individual deployment periods. This must contain a column that defines individual IDs (named 'individual_id') and the time of tagging (named 'tag_start_date') and time of tag retrieval ('tag_end_date') (see \code{\link[flapper]{dat_ids}} for an example).
#' @param moorings A dataframe that defines receiver deployment periods. This must contain a column that defines receiver IDs (named 'receiver_id') and the time of receiver deployment (named 'receiver_start_date') and retrieval (named 'receiver_end_date') (see \code{\link[flapper]{dat_moorings}} for an example).
#' @param individual_id (optional) A vector of individuals for which to calculate overlap duration.
#' @param receiver_id (optional) A vector of receivers for which to calculate overlap duration.
#' @param type If both \code{individual_id} and \code{receiver_id} are specified, then \code{type} is an integer that defines whether or not to calculate overlap duration for (a) each individual/receiver pair (\code{type = 1L}) or (b) all combinations of specified individuals/receivers (\code{type = 2L}).
#' @param match_to (optional) A dataframe against which to match the calculated overlap duration(s). This must contain an 'individual_id' and 'receiver_id' column, as in \code{ids} and \code{moorings} respectively. If supplied, an integer vector of overlap durations for individual/receiver combinations, matched against the individuals/receivers in this dataframe, is returned (see also Value).
#' @param ... Additional arguments (none implemented).
#'
#' @return The function returns a dataframe with the deployment overlap duration for specific or all combinations of individuals and receivers, with the 'individual_id', 'receiver_id', 'tag_start_date', 'tag_end_date', 'receiver_start_date' and 'receiver_end_date' columns retained. The 'get_id_rec_overlap' column defines the temporal overlap (days). Alternatively, if \code{match_to} is supplied, a vector of overlap durations that matches each individual/receiver observation in that dataframe is returned.
#'
#' @examples
#' #### Prepare data to include required columns
#' # moorings requires 'receiver_id', 'receiver_start_date' and 'receiver_end_date'
#' # ids requires 'individual_id', 'tag_start_date' and 'tag_end_date'
#' # These columns are already supplied in the example datasets
#' # ... except tag_end_date:
#' dat_ids$tag_end_date <- as.Date("2017-06-05")
#'
#' #### Example (1): Temporal between all combinations
#' # ... of individuals and receivers
#' dat <- get_id_rec_overlap(dat_ids, dat_moorings)
#'
#' #### Example (2) Temporal overlap between all combinations of specified
#' #... individuals/receivers
#' dat <- get_id_rec_overlap(dat_ids,
#'                       dat_moorings,
#'                       individual_id = c(25, 26),
#'                       receiver_id = c(3, 4),
#'                       type = 2L)
#'
#' #### Example (3) Temporal overlap between specified individual/receiver pairs
#' dat <- get_id_rec_overlap(dat_ids,
#'                       dat_moorings,
#'                       individual_id = c(25, 26),
#'                       receiver_id = c(3, 4),
#'                       type = 1L)
#'
#' #### Example (4) Match temporal overlap to another dataframe
#' dat_acoustics$get_id_rec_overlap <- get_id_rec_overlap(dat_ids,
#'                                                dat_moorings,
#'                                                match_to = dat_acoustics,
#'                                                type = 1L)
#'
#' @author Edward Lavender
#' @export
#'

get_id_rec_overlap <- function(ids,
                           moorings,
                           individual_id = NULL,
                           receiver_id = NULL,
                           type = 1L,
                           match_to = NULL,...){

  #### Checks
  # Dataframes must contains required names
  check_names(input = ids, req = c("individual_id", "tag_start_date", "tag_end_date"), extract_names = colnames, type = all)
  check_names(input = moorings, req = c("receiver_id", "receiver_start_date", "receiver_end_date"), extract_names = colnames, type = all)
  if(!is.null(match_to)) check_names(input = match_to, req = c("individual_id", "receiver_id"), extract_names = colnames, type = all)
  # Check input to type
  type <- check_value(input = type, supp = 1:2L)

  #### Define dataframe with individuals and receivers
  ## Option (1) Both individual_id and receiver_id have been supplied
  # ... in which case we will define a dataframe for these specific individuals based on type
  if(!is.null(individual_id) & !is.null(receiver_id)) {

    # If type == 1, then we will consider each pair of individuals and receivers
    if(type == 1L) {
      if(length(individual_id) != length(receiver_id)) {
        stop("Both 'individual_id' and 'receiver_id' have been specified and type = 1L but length(individual_id) != length(receiver_id).")
      }
      dat <- data.frame(individual_id = individual_id, receiver_id = receiver_id)

      # Otherwise, we will focus on all combinations of specified individuals and receivers
    } else if (type == 2L) {
      dat <- expand.grid(individual_id = individual_id, receiver_id = receiver_id)
    }

    ## Option (2) Use all combinations of individuals/receivers
  } else {

    # Filter out any unwanted individuals or receivers
    if(!is.null(individual_id)) ids <- ids[which(ids$individual_id %in% individual_id), ]
    if(!is.null(receiver_id)) moorings <- moorings[which(moorings$receiver_id %in% receiver_id), ]
    # Define dataframe
    # This will only include receivers that recorded detections
    dat <- expand.grid(individual_id = unique(ids$individual_id),
                       receiver_id = unique(moorings$receiver_id))

  }

  #### Define dates
  # Define start/end dates for individuals' time at liberty
  dat$tag_start_date      <- ids$tag_start_date[match(dat$individual_id, ids$individual_id)]
  dat$tag_end_date        <- ids$tag_end_date[match(dat$individual_id, ids$individual_id)]
  # Define start/end dates for receivers' deployment time
  dat$receiver_start_date <- moorings$receiver_start_date[match(dat$receiver_id, moorings$receiver_id)]
  dat$receiver_end_date   <- moorings$receiver_end_date[match(dat$receiver_id, moorings$receiver_id)]
  # Define intervals
  dat$tag_interval       <- lubridate::interval(dat$tag_start_date, dat$tag_end_date)
  dat$receiver_interval  <- lubridate::interval(dat$receiver_start_date, dat$receiver_end_date)

  #### Calculate overlap
  # Define the overlap in days, including the first day of overlap (+1)
  dat$id_rec_overlap <- lubridate::day(lubridate::as.period(lubridate::intersect(dat$tag_interval, dat$receiver_interval), "days")) + 1

  #### Match detection days to another dataframe, if requested
  if(!is.null(match_to)) {
    dat$key <- paste0(dat$individual_id, "-", dat$receiver_id)
    match_to$key  <- paste0(match_to$individual_id, "-", match_to$receiver_id)
    match_to$id_rec_overlap <- dat$id_rec_overlap[match(match_to$key, dat$key)]
    out <- match_to$id_rec_overlap
    if(any(is.na(out))){
      message(sum(is.na(out)), "NAs identified in matched vector of id_rec_overlap.")
    }
  } else{
    out <- dat
  }

  #### Return outputs
  return(out)

}


######################################
######################################
#### get_detection_centroids_envir()

#' @title Sample environmental conditions around receivers
#' @description This function is used to sample environmental conditions from within the detection centroids of receivers. To implement the function, a SpatialPoints object that defines receiver locations (\code{xy}) must be provided, along with the detection range (\code{detection_range}) of receivers. This information is used to define detection centroids, via \code{\link[flapper]{get_detection_centroids}}. Within each receiver's centroid, all values of an environmental variable, or a random sample of values, are extracted from a user-defined \code{\link[raster]{raster}} (\code{envir}). Under random sampling, values can be sampled according to a detection probability function (\code{sample_probs}). The function returns a list of dataframes, one for each receiver, that include the sampled values.

#' @param xy,detection_range,coastline,plot,... Arguments required to calculate and visualise detection centroids via \code{\link[flapper]{get_detection_centroids}}; namely, receiver locations (\code{xy}), the detection range (\code{detection_range}), barriers to detection (\code{coastline}) and whether or not to plot the centroids (\code{plot}). Additional arguments can be passed via \code{...} but note that \code{byid} is necessarily \code{TRUE} and should not be provided.
#' @param envir A \code{\link[raster]{raster}} that defines the values of an environmental variable across the study area. The coordinate reference system should be the Universal Transverse Mercator system.
#' @param sample_size (optional) An integer that defines the number of samples of the environmental variable to draw from the area around each receiver (see the 'size' argument of \code{\link[base]{sample}}). If this is provided, \code{sample_size} samples are taken from this area; otherwise, all values are extracted.
#' @param sample_replace (optional) If \code{sample_size} is specified, \code{sample_replace} is a logical input that defines whether to implement sampling with (\code{sample_replace = TRUE}, the default) or without (\code{sample_replace = FALSE}) replacement (see the 'replace' argument of \code{\link[base]{sample}}).
#' @param sample_probs (optional) If \code{sample_size} is specified, \code{sample_probs} is a function that calculates the detection probability given the distance (m) between an cell and a receiver.
#' @param cl A cluster object created by \code{\link[parallel]{makeCluster}}. If supplied, the connection to the cluster is closed within the function.
#' @param varlist A character vector of names of objects to export, to be passed to the \code{varlist} argument of \code{\link[parallel]{clusterExport}}. This may be required if \code{cl} is supplied. Exported objects must be located in the global environment.
#' @param verbose A logical variable that defines whether or not relay messages to the console to monitor function progress.
#'
#' @return The function returns a list of dataframes (one for each element in \code{xy} i.e., each receiver), each of which includes the cell IDs of \code{envir} from which values were extracted ('cell'), the value of the environmental variable in that cell ('envir') and, if applicable, the distance between that cell and the receiver ('dist', m) and the detection probability in that cell ('prob').
#'
#' @examples
#' #### Define receiver locations as a SpatialPoints object with a UTM CRS
#' proj_wgs84 <- sp::CRS("+init=epsg:4326")
#' proj_utm <- sp::CRS(paste("+proj=utm +zone=29 +datum=WGS84",
#'                           "+units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))
#' xy <- sp::SpatialPoints(dat_moorings[, c("receiver_long", "receiver_lat")],
#'                         proj_wgs84)
#' xy <- sp::spTransform(xy, proj_utm)
#'
#' #### Example (1): Extract all depth values within each receiver's centroid
#' depths_by_centroid <- get_detection_centroids_envir(xy = xy,
#'                                                 detection_range = 425,
#'                                                 coastline = dat_coast,
#'                                                 envir = dat_gebco
#'                                                 )
#' # The function returns a list of dataframes, one for each receiver
#' # ... with the cell IDs and the value of the environmental variable
#' utils::str(depths_by_centroid)
#' # Collapse the list and compare conditions across receivers
#' depths_by_centroid <-
#'   lapply(1:length(depths_by_centroid), function(i){
#'     d <- depths_by_centroid[[i]]
#'     d$receiver_id <- dat_moorings$receiver_id[i]
#'     return(d)
#'   })
#' depths_by_centroid <- dplyr::bind_rows(depths_by_centroid)
#' prettyGraphics::pretty_boxplot(depths_by_centroid$receiver_id,
#'                                depths_by_centroid$envir)
#'
#' #### Example (2): Extract a random sample of values
#' # (We'll keep the values small for speed)
#' depths_by_centroid <- get_detection_centroids_envir(xy = xy,
#'                                                 detection_range = 425,
#'                                                 coastline = dat_coast,
#'                                                 envir = dat_gebco,
#'                                                 sample_size = 2
#'                                                  )
#' utils::str(depths_by_centroid)
#'
#' #### Example (3) Extract a random sample of values with weighted probabilities
#' # Define detection probability function based only on distance
#' calc_detection_pr <-
#'   function(dist){
#'     dpr <- get_detection_pr(distance = dist,
#'                         beta_0 = 2.5,
#'                         beta_1 = -0.01,
#'                         inv_link = stats::plogis,
#'                         output = 2L)
#'     return(dpr)
#'   }
#' # Implement sampling with replacement according to detection probability
#' depths_by_centroid <- get_detection_centroids_envir(xy = xy,
#'                                                 detection_range = 425,
#'                                                 coastline = dat_coast,
#'                                                 envir = dat_gebco,
#'                                                 sample_size = 2,
#'                                                 sample_probs = calc_detection_pr,
#'                                                 )
#' # Each element of the outputted list includes the 'cell' and 'envir' column
#' # ... as well as 'dist' and 'prob' that define the distance of that cell
#' # ... from the location in xy and the corresponding detection probability
#' # ... at that distance respectively
#' utils::str(depths_by_centroid)
#'
#' #### Example (4) Sampling without replacement via sample_replace = FALSE
#' depths_by_centroid <- get_detection_centroids_envir(xy = xy,
#'                                                 detection_range = 425,
#'                                                 coastline = dat_coast,
#'                                                 envir = dat_gebco,
#'                                                 sample_size = 2,
#'                                                 sample_probs = calc_detection_pr,
#'                                                 sample_replace = FALSE
#'                                                 )
#' utils::str(depths_by_centroid)
#'
#' #### Example (5) Parallelise the algorithm via cl and varlist arguments
#' \dontrun{
#' depths_by_centroid <- get_detection_centroids_envir(xy = xy,
#'                                                 detection_range = 425,
#'                                                 coastline = dat_coast,
#'                                                 envir = dat_gebco,
#'                                                 sample_size = 2,
#'                                                 sample_probs = calc_detection_pr,
#'                                                 sample_replace = FALSE,
#'                                                 cl = parallel::makeCluster(2L),
#'                                                 varlist = c("dat_gebco",
#'                                                             "calc_detection_pr")
#'                                                 )
#' utils::str(depths_by_centroid)
#' }
#'
#' @author Edward Lavender
#' @export
#'

get_detection_centroids_envir <- function(xy,
                                      detection_range,
                                      coastline,
                                      plot = FALSE,
                                      envir,
                                      sample_size = NULL,
                                      sample_replace = TRUE,
                                      sample_probs = NULL,
                                      cl = NULL,
                                      varlist = NULL,
                                      verbose = TRUE,...){

  #### Checks
  t_onset <- Sys.time()
  cat_to_console <- function(..., show = verbose) if(show) cat(paste(..., "\n"))
  cat_to_console(paste0("flapper::get_detection_centroids_envir() called (@ ", t_onset, ")..."))
  cat_to_console("... Implementing function checks...")
  if(is.null(sample_size)){
    if(!is.null(sample_probs)) message("sample_size = NULL: input to 'sample_probs' ignored.")
    sample_probs <- NULL
  }
  if(is.null(cl) & !is.null(varlist)) message("cl = NULL: input to 'varlist' ignored.")
  check...("byid",...)

  #### Define detection centroids
  cat_to_console("... Defining detection centroid(s)...")
  xy_buf <- get_detection_centroids(xy = xy,
                                detection_range = detection_range,
                                coastline = coastline,
                                plot = plot,
                                byid = TRUE,...)
  xy_buf_ls <- lapply(1:length(xy_buf), function(i) xy_buf[i])

  #### Extract conditions
  cat_to_console("... Extracting environmental conditions from detection area(s)...")
  if(!is.null(cl) & !is.null(varlist)) {
    parallel::clusterExport(cl = cl, varlist = varlist)
  }
  ls_envir <-
    pbapply::pblapply(xy_buf_ls, cl = cl, FUN = function(centroid){

      ## Extract conditions
      # Create list of conditions sampled by each receiver
      envir_sampled <- raster::extract(envir, centroid, cellnumbers = TRUE)
      envir_sampled <- envir_sampled[[1]]
      dat <- data.frame(envir_sampled)
      colnames(dat) <- c("cell", "envir")

      ## Define distances
      if(!is.null(sample_probs)){
        rdist <- raster::distanceFromPoints(envir, xy)
        dist_sampled <- raster::extract(rdist, centroid, cellnumbers = TRUE)
        dist_sampled <- dist_sampled[[1]]
        dist_sampled <- data.frame(dist_sampled)
        colnames(dist_sampled) <- c("cell", "dist")
        dat$dist <- dist_sampled$dist[match(dat$cell, dist_sampled$cell)]
      }

      ## Return outputs
      return(dat)

    })
  if(!is.null(cl)) parallel::stopCluster(cl)

  #### Sample values according to their probability
  if(!is.null(sample_size)) {
    ls_envir_sample <- lapply(ls_envir, function(d){
      if(!is.null(sample_probs)){
        d$prob <- sample_probs(d$dist)
      } else d$prob <- NULL
      envir_sampled <- d[sample(1:nrow(d), size = sample_size, replace = sample_replace, prob = d$prob), ]
      return(envir_sampled)
    })
  } else{
    ls_envir_sample <- ls_envir
  }

  #### Return outputs
  return(ls_envir_sample)

}


######################################
######################################
#### get_detection_days()

#' @title Calculate detection days
#' @description The function calculates the total number of days (termed 'detection days') during which individuals were detected at passive acoustic telemetry receivers. To implement the function, a dataframe with passive acoustic telemetry detections of individuals at receivers must be supplied. Detection days can be calculated for all combinations of individuals/receivers in this dataframe, for all combinations of specified individuals and receivers, or for specific individual/receiver pairs. The function returns a dataframe of detection days for these individual/receiver combinations or a vector of detection days that is matched against another dataframe.
#'
#' @param acoustics A dataframe that contains passive acoustic telemetry detection time series (see \code{\link[flapper]{dat_acoustics}} for an example). This should contain the following columns: a vector of individual IDs, named 'individual_id'; a vector of receiver IDs, named 'receiver_id'; and a POSIXct vector of time stamps when detections were made, named 'timestamp'.
#' @param individual_id (optional) A vector of individuals for which to calculate detection days.
#' @param receiver_id (optional) A vector of receivers for which to calculate detection days.
#' @param type If both \code{individual_id} and \code{receiver_id} are specified, then \code{type} is an integer that defines whether or not to calculate detection days for (a) each individual/receiver pair (\code{type = 1L}) or (b) all combinations of individuals/receivers.
#' @param match_to (optional) A dataframe against which to match detection days. This must contain the 'individual_id' and 'receiver_id' column, as in \code{acoustics}. If supplied, an integer vector of detection days for individual/receiver combinations, matched against the individuals/receivers in \code{match_to}, is returned (see also Value).
#' @param ... Additional arguments passed to \code{\link[pbapply]{pblapply}}.
#'
#' @return The function returns a dataframe with the detection days for all, or specified, combinations of individuals and receivers. Note that if \code{acoustics} only contains individuals/receivers that made detections, then this will only contain individuals/receivers for/at which detections were made. Alternatively, if \code{match_to} is supplied, a vector of detection days, matched against each individual/receiver observation in that dataframe, is returned.
#'
#' @examples
#' #### Example (1): Detection days between all combinations
#' # ... of detected individuals and receivers with detections
#' dat <- get_detection_days(dat_acoustics)
#'
#' #### Example (2) Detection days between all combinations of specified
#' #... individuals/receivers
#' dat <- get_detection_days(dat_acoustics,
#'                       individual_id = 25,
#'                       receiver_id = c(25, 21),
#'                       type = 2L)
#'
#' #### Example (3) Detection days between specified individual/receiver pairs
#' dat <- get_detection_days(dat_acoustics,
#'                       individual_id = c(25, 26),
#'                       receiver_id = c(25, 21),
#'                       type = 1L)
#'
#' #### Example (4) Match detection days to another dataframe
#' dat_acoustics$detection_days <- get_detection_days(dat_acoustics,
#'                                                match_to = dat_acoustics,
#'                                                type = 1L)
#'
#' utils::head(dat_acoustics)
#'
#' @author Edward Lavender
#' @export
#'

get_detection_days <- function(acoustics,
                           individual_id = NULL,
                           receiver_id = NULL,
                           type = 1L,
                           match_to = NULL,...){
  #### Checks
  # Dataframes must contains required names
  check_names(input = acoustics, req = c("individual_id", "receiver_id", "timestamp"), extract_names = colnames, type = all)
  if(!is.null(match_to)) check_names(input = match_to, req = c("individual_id", "receiver_id"), extract_names = colnames, type = all)
  # Check input to type
  type <- check_value(input = type, supp = 1:2L)

  #### Define dataframe with individuals and receivers
  ## Option (1) Both individual_id and receiver_id have been supplied
  # ... in which case we will define a dataframe for these specific individuals based on type
  if(!is.null(individual_id) & !is.null(receiver_id)) {

    # If type == 1, then we will consider each pair of individuals and receivers
    if(type == 1L) {
      if(length(individual_id) != length(receiver_id)) {
        stop("Both 'individual_id' and 'receiver_id' have been specified and type = 1L but length(individual_id) != length(receiver_id).")
      }
      dat <- data.frame(individual_id = individual_id, receiver_id = receiver_id)

      # Otherwise, we will focus on all combinations of specified individuals and receivers
    } else if (type == 2L) {
      dat <- expand.grid(individual_id = individual_id, receiver_id = receiver_id)
    }

    ## Option (2) Use all combinations of individuals/receivers
  } else {

    # Filter out any unwanted individuals or receivers
    if(!is.null(individual_id)) acoustics <- acoustics[which(acoustics$individual_id %in% individual_id), ]
    if(!is.null(receiver_id)) acoustics <- acoustics[which(acoustics$receiver_id %in% receiver_id), ]
    # Define dataframe
    # This will only include receivers that recorded detections
    dat <- expand.grid(individual_id = unique(acoustics$individual_id),
                       receiver_id = unique(acoustics$receiver_id))

  }

  #### Calculate detections days
  acoustics$date <- as.Date(acoustics$timestamp)
  detection_days_ls <-
    pbapply::pblapply(split(dat, 1:nrow(dat)), FUN = function(d){
      # Subset acoustics to consider correct individual and receiver
      aco_tmp <- acoustics[acoustics$individual_id == d$individual_id &
                             acoustics$receiver_id == d$receiver_id, ]
      # Count the number of detections and return
      if(nrow(aco_tmp) > 0) n <- length(unique(aco_tmp$date)) else n <- 0
      return(n)
    },...)
  dat$detection_days <- unlist(detection_days_ls)

  #### Match detection days to another dataframe, if requested
  if(!is.null(match_to)) {
    dat$key <- paste0(dat$individual_id, "-", dat$receiver_id)
    match_to$key  <- paste0(match_to$individual_id, "-", match_to$receiver_id)
    match_to$detection_days <- dat$detection_days[match(match_to$key, dat$key)]
    out <- match_to$detection_days
    if(any(is.na(out))){
      message(sum(is.na(out)), "NAs identified in matched vector of detection days.")
    }
  } else{
    out <- dat
  }

  #### Return outputs
  return(out)

}
