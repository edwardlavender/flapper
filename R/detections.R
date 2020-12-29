######################################
######################################
#### detection_pr()

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
#' det_pr <- detection_pr()
#' utils::head(det_pr)
#' # The vector has attributes:
#' # ... 'X' (the model matrix),
#' # ... 'beta' (the regression coefficient)
#' # ... 'inv_link' (the inverse link function)
#' utils::str(det_pr)
#'
#' #### Example (2): Adjust model parameters
#' # Change regression coefficients
#' det_pr <- detection_pr(beta_0 = 2.5, beta_1 = -0.006)
#' # Use inverse probit
#' det_pr <- detection_pr(beta_0 = 2.5, beta_1 = -0.006, inv_link = stats::pnorm)
#'
#' #### Example (3): Modify graphical properties
#' det_pr <- detection_pr(beta_0 = 2.5,
#'                           beta_1 = -0.006,
#'                           type = "l",
#'                           xlab = "Distance (m)",
#'                           ylab = "Detection Probability")
#'
#' #### Example (4): Modify return options
#' # Only graph
#' detection_pr(output = 1L)
#' # Only values
#' detection_pr(output = 2L)
#' # Both graph and values (the default)
#' detection_pr(output = 3L)
#'
#' @author Edward Lavender
#' @export
#'

detection_pr <- function(distance = 1:1000,
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
#### detection_centroids()

#' @title Define detection centroids around receivers
#' @description This function defines the areas surveyed by receivers (termed 'detection centroids') as a spatial object, based on an estimate of the detection range (m) and any barriers to detection. To implement the function, receiver locations must be supplied as a SpatialPoints or SpatialPointsDataFrame object with the Universe Transverse Mercator coordinate reference system. The function defines a spatial buffer around each receiver according to the estimated detection range, cuts out any barriers to detection, such as the coastline, and returns a SpatialPolygons object that defines the combined detection centroid across all receivers or receiver-specific detection centroids.
#'
#' @param xy A \code{\link[sp]{SpatialPoints-class}} or \code{\link[sp]{SpatialPointsDataFrame-class}} object that defines receiver locations. The coordinate reference system should be the Universe Transverse Mercator coordinate reference system.
#' @param detection_range A number that defines the detection range (m) of receivers.
#' @param coastline (optional) A \code{\link[sp]{SpatialPolygonsDataFrame-class}} object that defines barriers (such as the coastline) that block receivers from surveying areas within their detection range.
#' @param plot A logical input that defines whether or not to plot receivers, their centroids, and the buffer (if specified).
#' @param ... Additional arguments passed to \code{\link[rgeos]{gBuffer}}, such as \code{byid} and/or \code{quadsegs}.
#'
#' @return The function returns a \code{\link[sp]{SpatialPolygons-class}} object of the detection centroids around receivers that represents the area they survey under the assumption of a constant detection range, accounting for any barriers to detection. By default, this will contain a single feature, which is suitable for the calculation of the total area surveyed by receivers (see \code{\link[flapper]{detection_area_sum}}) because it accounts for the overlap in the detection ranges of receivers. However, if \code{byid = TRUE} is passed via \code{...} to \code{\link[rgeos]{gBuffer}}, the returned object will have a feature for each pair of coordinates in (\code{xy}) (i.e., receiver). This is less appropriate for calculating the area surveyed by receivers, since areas surveyed by multiple receivers will be over-counted, but it is suitable when the centroids for particular receivers are required (e.g., to extract environmental conditions within a specific receiver's detection ranges) (see \code{\link[flapper]{detection_centroids_envir}}.
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
#' detection_centroids(xy)
#'
#' #### Example (2): Account for barriers in the study area
#' detection_centroids(xy, coastline = dat_coast)
#'
#' #### Example (3): Adjust the detection range
#' detection_centroids(xy, detection_range = 400, coastline = dat_coast)
#' detection_centroids(xy, detection_range = 500, coastline = dat_coast)
#'
#' #### Example (4): Suppress the plot
#' detection_centroids(xy, coastline = dat_coast, plot = FALSE)
#'
#' #### Example (5): Output characteristics are controlled via byid
#' # A SpatialPolygons with one feature is the implicit output
#' sp_1 <- detection_centroids(xy, coastline = dat_coast, byid = FALSE)
#' sp_1
#' # An SpatialPolygons with one feature for each element in xy
#' # ... can be returned via byid = TRUE
#' sp_2 <- detection_centroids(xy, coastline = dat_coast, byid = TRUE)
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

detection_centroids <- function(xy,
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
#### detection_area_sum()

#' @title Calculate the total area sampled by acoustic receivers
#' @description This function calculates the total area sampled by receivers, under the assumption of a constant detection range. To implement the function, receiver locations must be supplied as a SpatialPoints or SpatialPointsDataFrame object with the Universe Transverse Mercator coordinate reference system. The \code{\link[flapper]{detection_centroids}} is used to calculate the detection centroids around receivers, given a specified detection range (m) and any barriers to detection, such as coastline, and then the total area covered by receivers is calculated, accounting for overlapping centroids.
#'
#' @param xy,detection_range,coastline,plot,... Arguments required to calculate and visualise detection centroids, using \code{\link[flapper]{detection_centroids}}; namely, receiver locations (\code{xy}), the detection range (\code{range}), barriers to detection (\code{coastline}), and whether or not to plot the centroids (\code{plot}).
#' @param scale A number that scales the total area (m). The default (\code{1/(1000^2)}) converts the units of \eqn{m^2} to \eqn{km^2}.
#'
#' @details This is a simple metric of the overall receiver sampling effort. This may be a poor metric if the assumption of a single detection range across all receivers is substantially incorrect or if there are substantial changes in the receiver array over the course of a study.
#'
#' @return The function returns a number that represents the total area surveyed by receivers (by default in \eqn{km^2}) and, if \code{plot = TRUE}, a plot of the area with receivers and their detection ranges.
#'
#' @seealso \code{\link[flapper]{detection_centroids}} defines detection centroids, across which the detection area is calculated. \code{\link[flapper]{detection_area_ts}} calculates the area sampled by receivers through time.
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
#' detection_area_sum(xy)
#'
#' #### Example (2): Account for barriers in the study area
#' detection_area_sum(xy, coastline = dat_coast)
#'
#' #### Example (3): Adjust the detection range
#' detection_area_sum(xy, detection_range = 400, coastline = dat_coast)
#' detection_area_sum(xy, detection_range = 500, coastline = dat_coast)
#'
#' #### Example (4): Adjust the units
#' detection_area_sum(xy, coastline = dat_coast, scale = 1) # m2
#'
#' #### Example (5): Suppress the plot
#' detection_area_sum(xy, coastline = dat_coast, plot = FALSE)
#'
#' @author Edward Lavender
#' @export
#'

detection_area_sum <- function(xy,
                               detection_range = 425,
                               coastline = NULL,
                               scale = 1/(1000^2),
                               plot = TRUE,...){

  #### Checks
  # If xy is empty, return area = 0
  if(length(xy) == 0) return(0)

  #### Define detection centroids
  xy_buf <- detection_centroids(xy = xy,
                                detection_range = detection_range,
                                coastline = coastline,
                                plot = plot,...)

  #### Calculate area (m2) covered by buffers overall
  xy_area <- rgeos::gArea(xy_buf) * scale
  return(xy_area)
}


######################################
######################################
#### detection_area_ts()

#' @title Calculate the area sampled by receivers through time
#' @description This function extends \code{\link[flapper]{detection_area_sum}} to calculate how the total area sampled by receivers changes through time.
#'
#' @param xy,detection_range,coastline,scale Arguments required to calculate the total area surveyed by receivers (at each time point) via \code{\link[flapper]{detection_area_sum}}.
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
#' dat <- detection_area_ts(xy)
#'
#' #### Example (2): Adjust detection range, include coastline and use parallel processing
#' # For areas with complex coastline, this will reduce the speed of the algorithm
#' # So we will also supply a cluster to improve the computation time.
#' \dontrun{
#' dat <- detection_area_ts(xy,
#'                        detection_range = 500,
#'                        coastline = dat_coast,
#'                        cl = parallel::makeCluster(2L),
#'                        varlist = "dat_coast"
#'                        )
#' }
#'
#' #### Example (3) Hide or customise the plot
#' dat <- detection_area_ts(xy, plot = FALSE)
#' dat <- detection_area_ts(xy,
#'                        pretty_axis_args = list(axis = list(list(format = "%b-%y"),
#'                                                            list())),
#'                        xlab = "Time (month-year)",
#'                        ylab = expression(paste("Area (", m^2, ")")),
#'                        type = "l")
#'
#' @author Edward Lavender
#' @export
#'

detection_area_ts <- function(xy,
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
  cat_to_console(paste0("flapper::detection_area_ts() called (@ ", t_onset, ")..."))
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
    receiver_area <- detection_area_sum(xy = xy[pos, ],
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
  cat_to_console(paste0("... flapper::detection_area_ts() call completed (@ ", t_end, ") after ~", round(duration, digits = 2), " minutes."))
  return(rcov)

}


######################################
######################################
#### n_operational_ts()

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
#' @details This is a simple metric of sampling effort. For acoustic receivers, \code{\link[flapper]{detection_area_ts}} provides another metric of sampling effort.
#'
#' @return The function returns a dataframe that, for each time step ('time'), defines the number of operational units at that time ('n'). If \code{plot = TRUE}, the function also plots a time series of the number of operational units.
#'
#' @examples
#' #### Example (1): Number of operational receivers over an acoustic telemetry study
#' dat_n <- n_operational_ts(data = dat_moorings,
#'                           start = "receiver_start_date",
#'                           stop = "receiver_end_date")
#' utils::head(dat_n)
#'
#' #### Example (2): Number of individuals at liberty over a tagging study
#' # Define 'tag_end_date' as hypothetical end date of a study
#' # ...  and assume that all individuals remained tagged until this time
#' dat_ids$tag_end_date <- as.Date("2017-06-05")
#' dat_n <- n_operational_ts(data = dat_ids,
#'                           start = "tag_start_date",
#'                           stop = "tag_end_date")
#'
#' #### Example (3): Specify the time period under consideration
#' dat_n <- n_operational_ts(data = dat_ids,
#'                           start = "tag_start_date",
#'                           stop = "tag_end_date",
#'                           times = seq(min(dat_moorings$receiver_start_date),
#'                                       max(dat_moorings$receiver_end_date), 1)
#'                           )
#'
#' #### Example (4): Suppress or customise the plot
#' dat_n <- n_operational_ts(data = dat_ids,
#'                           start = "tag_start_date",
#'                           stop = "tag_end_date",
#'                           plot = FALSE)
#' dat_n <- n_operational_ts(data = dat_ids,
#'                           start = "tag_start_date",
#'                           stop = "tag_end_date",
#'                           xlab = "Time", ylab = "N (individuals)",
#'                           type = "l")
#'
#' @author Edward Lavender
#' @export
#'

n_operational_ts <- function(data, start, stop, times = NULL, plot = TRUE,...){

  #### Define a sequence of times that span the range of the data, if required.
  if(is.null(times)){
    times <- seq(as.Date(min(data[, start]), na.rm = TRUE), as.Date(max(data[, stop]), na.rm = TRUE), by = "days")
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


