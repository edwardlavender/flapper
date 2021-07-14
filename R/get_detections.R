######################################
######################################
#### get_detection_pr()

#' @title A detection probability function based on distance
#' @description This function calculates detection probability (e.g., of an acoustic detection) at specified distances from the sampling device (e.g., a passive acoustic telemetry receiver) using user-defined parameters (i.e., a model intercept, a coefficient for the effect of distance and an inverse link function). The function returns a plot of detection probability with distance and/or a vector of detection probabilities.
#'
#' @param distance A numeric vector of distances at which to calculate detection probability.
#' @param beta_0,beta_1 Single numbers that define the model coefficients (i.e., the intercept and gradient on the scale of the link function).
#' @param inv_link A function that defines the inverse link function.The default function is the logistic (inverse logit) function.
#' @param output An integer (\code{1L}, \code{2L} or \code{3L}) that defines the output type. \code{1L} returns a plot of detection probability against distance; \code{2L} returns a numeric vector of detection probabilities; and \code{3L} returns both of the above.
#' @param ... Additional arguments, passed to \code{\link[prettyGraphics]{pretty_plot}}, to customise the plot. These are only implemented if \code{output = 1L} or \code{output = 3L}.
#'
#' @return The function calculates detection probability at each specified distance and returns a plot, a vector of detection probabilities, or both, depending on the value of the \code{output} argument. If a vector of detection probabilities is returned, this contains the following attributes: `X', the model matrix; `beta', the regression coefficients; and `inv_link', the inverse link function.
#'
#' @examples
#' #### Example (1): Implement the function using the default parameters
#' # The function returns a graph and a vector of detection probabilities
#' det_pr <- get_detection_pr()
#' utils::head(det_pr)
#' # The vector has attributes:
#' # ... 'X' (the model matrix)
#' # ... 'beta' (the regression coefficient)
#' # ... 'inv_link' (the inverse link function)
#' utils::str(det_pr)
#'
#' #### Example (2): Adjust model parameters
#' # Change regression coefficients
#' det_pr <- get_detection_pr(beta_0 = 2.5, beta_1 = -0.006)
#' # Use inverse probit link function
#' det_pr <- get_detection_pr(beta_0 = 2.5, beta_1 = -0.006, inv_link = stats::pnorm)
#'
#' #### Example (3): Modify graphical properties
#' det_pr <- get_detection_pr(beta_0 = 2.5,
#'                            beta_1 = -0.006,
#'                            type = "l",
#'                            xlab = "Distance (m)",
#'                            ylab = "Detection Probability")
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
#' @description This function defines the areas surveyed by receivers (termed `detection centroids') as a spatial object, based on an estimate of the detection range (m) and any barriers to detection. To implement the function, receiver locations must be supplied as a SpatialPoints or SpatialPointsDataFrame object with the Universe Transverse Mercator coordinate reference system. The function defines a spatial buffer around each receiver according to the estimated detection range, cuts out any barriers to detection, such as the coastline, and returns a SpatialPolygons object that defines the combined detection centroid across all receivers or receiver-specific detection centroids.
#'
#' @param xy A \code{\link[sp]{SpatialPoints-class}} or \code{\link[sp]{SpatialPointsDataFrame-class}} object that defines receiver locations. The coordinate reference system should be the Universe Transverse Mercator coordinate reference system.
#' @param detection_range A number that defines the detection range (m) of receivers.
#' @param coastline (optional) A \code{\link[sp]{SpatialPolygonsDataFrame-class}} object that defines barriers (such as the coastline) that block receivers from surveying areas within their detection range.
#' @param plot A logical input that defines whether or not to plot receivers, their centroids, and the buffer (if specified).
#' @param ... Additional arguments passed to \code{\link[rgeos]{gBuffer}}, such as \code{byid} and/or \code{quadsegs}.
#'
#' @return The function returns a \code{\link[sp]{SpatialPolygons-class}} object of the detection centroids around receivers that represents the area they survey under the assumption of a constant detection range, accounting for any barriers to detection. By default, this will contain a single feature, which is suitable for the calculation of the total area surveyed by receivers (see \code{\link[flapper]{get_detection_area_sum}}) because it accounts for the overlap in the detection ranges of receivers. However, if \code{byid = TRUE} is passed via \code{...} to \code{\link[rgeos]{gBuffer}}, the returned object will have a feature for each pair of coordinates in \code{xy} (i.e., receiver). This is less appropriate for calculating the area surveyed by receivers, since areas surveyed by multiple receivers will be over-counted, but it is suitable when the centroids for particular receivers are required (e.g., to extract environmental conditions within a specific receiver's detection range) (see \code{\link[flapper]{get_detection_centroids_envir}}).
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
#' #### Example (1): Get the simplest centroids around receivers
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
#' # A SpatialPolygons object with one feature is the implicit output
#' sp_1 <- get_detection_centroids(xy, coastline = dat_coast, byid = FALSE)
#' sp_1
#' # A SpatialPolygons object with one feature for each element in xy
#' # ... can be returned via byid = TRUE
#' sp_2 <- get_detection_centroids(xy, coastline = dat_coast, byid = TRUE)
#' sp_2
#' # The total area of the former will be smaller, since areas covered
#' # ... by multiple receivers are merged
#' rgeos::gArea(sp_1); rgeos::gArea(sp_2)
#' # But it can be more convenient to use the latter format in some cases
#' # ... because it is easy to isolate specific centroids:
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


########################################
########################################
#### get_detection_centroids_overlap()

#' @title Get detection centroid overlaps
#' @importFrom lubridate `%within%`
#' @description This functions identifies receivers with overlapping detection centroids in space and time.
#'
#' @param centroids A \code{\link[sp]{SpatialPolygonsDataFrame}} that defines detection centroids (see \code{\link[flapper]{get_detection_centroids}}). The \code{data} slot must include a dataframe with the following columns: an unique, integer identifier for each receiver (`receiver_id') and receiver deployment \code{\link[base]{Dates}} (`receiver_start_date' and `receiver_end_date').
#' @param services (optional) A dataframe that defines receiver IDs and servicing \code{\link[base]{Dates}} (times during the deployment period of a receiver when it was not active due to servicing). If provided, this must contain the following columns: an integer identifier for serviced receivers (named ‘receiver_id’) and two columns that define the time of the service(s) (‘service_start_date’ and ‘service_end_date’) (see \code{\link[flapper]{make_matrix_receivers}}).
#' @param ... Additional arguments (none implemented).
#'
#' @details This function requires the \code{\link[tidyr]{tidyr-package}} (specifically \code{\link[tidyr]{pivot_longer}}).
#'
#' @return The function returns a list with two elements:
#' \itemize{
#'   \item \strong{overlap_by_receiver} is list, with one element for all integers from \code{1:max(centroids$receiver_id)}. Any elements that do not correspond to receivers contain a NULL element. List elements that correspond to receivers contain a dataframe that defines, for each day over the deployment period (defined in `timestamp') of that receiver (defined in `receiver_id'), whether (1) or not (0) that receiver overlapped in space with every other receiver (defined in the remaining columns by their receiver IDs).
#'   \item \strong{overlap_by_date} is a named list, with one element for each date from the start until the end of the study (\code{min(centroids$receiver_start_date):max(centroids$receiver_end_date)}), that records an integer vector of all receivers with overlapping centroids on that date. In this vector, each receiver overlaps with at least one other receiver (but not every receiver will necessarily overlap with every other receiver).
#' }
#'
#' @examples
#' #### Define receiver centroids
#' ## Define receiver locations as a SpatialPoints object with a UTM CRS
#' proj_wgs84 <- sp::CRS("+init=epsg:4326")
#' proj_utm <- sp::CRS(paste("+proj=utm +zone=29 +datum=WGS84",
#'                           "+units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))
#' rownames(dat_moorings) <- dat_moorings$receiver_id
#' xy <- sp::SpatialPoints(dat_moorings[, c("receiver_long", "receiver_lat")],
#'                         proj_wgs84)
#' xy <- sp::spTransform(xy, proj_utm)
#' ## Get receiver-specific detection centroids
#' # ... via get_detection_centroids with byid = TRUE
#' centroids <- get_detection_centroids(xy, byid = TRUE)
#' ## Link detection centroids with receiver IDs and deployment dates
#' # ... in a SpatialPointsDataFrame, as required for this function.
#' centroids_df <- dat_moorings[, c("receiver_id",
#'                                  "receiver_start_date",
#'                                  "receiver_end_date")]
#' row.names(centroids_df) <- names(centroids)
#' centroids <- sp::SpatialPolygonsDataFrame(centroids, centroids_df)
#'
#' ## Simulate some receiver 'servicing' dates for demonstration purposes
#' set.seed(1)
#' # Loop over each receiver...
#' services_by_receiver <- lapply(split(dat_moorings, 1:nrow(dat_moorings)), function(din){
#'   # For the receiver, simulate the number of servicing events
#'   n <- sample(0:3, 1)
#'   dout <- NULL
#'   if(n > 0){
#'     # simulate the timing of servicing events
#'     dates <- sample(seq(min(din$receiver_start_date), max(din$receiver_end_date), "days"), n)
#'     dout <- data.frame(receiver_id = rep(din$receiver_id, length(dates)),
#'                        service_start_date = dates,
#'                        service_end_date = dates)
#'   }
#'   return(dout)
#' })
#' services <- do.call(rbind, services_by_receiver)
#' rownames(services) <- NULL
#' if(nrow(services) == 0) services <- NULL
#'
#' #### Example (1): Implement function using centroids alone
#' overlaps_1 <- get_detection_centroids_overlap(centroids = centroids)
#' summary(overlaps_1)
#'
#' #### Example (2): Account for servicing dates
#' overlaps_2 <- get_detection_centroids_overlap(centroids = centroids,
#'                                               services = services)
#' # Examine the first few simulated servicing events
#' services[1:3, ]
#' # Show that the list_by_date element for the first servicing event
#' # ... includes the first receiver in services$receiver_id
#' # ... for overlaps_1 but not overlaps_2,
#' # ... which accounts for that fact that the receiver was serviced then:
#' overlaps_1$list_by_date[[as.character(services$service_start_date[1])]]
#' overlaps_2$list_by_date[[as.character(services$service_start_date[1])]]
#' # Likewise, show that the list_by_receiver element for that receiver
#' # ... includes overlapping receivers in overlaps_1 but not overlaps_2:
#' r_id <- services$receiver_id[1]
#' overlaps_1$list_by_receiver[[r_id]][overlaps_1$list_by_receiver[[r_id]]$timestamp %in%
#'                                       services$service_start_date[services$receiver_id == r_id], ]
#' overlaps_2$list_by_receiver[[r_id]][overlaps_2$list_by_receiver[[r_id]]$timestamp %in%
#'                                       services$service_start_date[services$receiver_id == r_id], ]
#'
#' @seealso \code{\link[flapper]{get_detection_centroids}} creates detection centroids.
#' @author Edward Lavender
#' @export


get_detection_centroids_overlap <- function(centroids, services = NULL,...){

  #### Checks
  ## packages
  if(!requireNamespace("tidyr", quietly = TRUE)) stop("Please install 'tidyr': this function requires the tidyr::pivot_longer() routine.")
  ## centroids
  if(!inherits(centroids, "SpatialPolygonsDataFrame")) stop("'centroids' must be a SpatialPolygonsDataFrame.")
  ## moorings
  moorings <- data.frame(centroids)
  check_names(input = moorings, req = c("receiver_id", "receiver_start_date", "receiver_end_date"),
              extract_names = colnames, type = all)
  if(is.numeric(moorings$receiver_id)) moorings$receiver_id <- as.integer(moorings$receiver_id)
  if(!is.integer(moorings$receiver_id))
    stop(paste("Argument 'xy$receiver_id' must be of class 'integer', not class(es):"), class(moorings$receiver_id))
  if(any(moorings$receiver_id <= 0))
    stop("Argument 'xy$receiver_id' cannot contain receiver IDs <= 0.")
  if(any(duplicated(moorings$receiver_id )))
    stop("Argument 'xy$receiver_id' contains duplicate elements.")
  ## services
  if(!is.null(services)){
    check_names(input = services, req = c("receiver_id", "service_start_date", "service_end_date"),
                extract_names = colnames, type = all)
    if(is.numeric(services$receiver_id)) services$receiver_id <- as.integer(services$receiver_id)
    if(!is.integer(services$receiver_id))
      stop(paste("Argument 'services$receiver_id' must be of class 'integer', not class(es):"), class(services$receiver_id))
    if(!all(unique(services$receiver_id) %in% unique(moorings$receiver_id))){
      message("Not all receivers in services$receiver_id are in moorings$receiver_id.")
    }
    services$interval <- lubridate::interval(services$service_start_date,
                                             services$service_end_date)
  }
  ## dots
  if(!is.null(names(list(...)))){
    warning(paste0("The following argument(s) passed via ... are not supported: ",
                   paste(names(list(...)), collapse = ", "), "."),
            call. = FALSE, immediate. = TRUE)
  }

  #### Define receiver activity status matrix
  # This defines whether or not each receiver was active on each date (0, 1)
  # We'll start from this point because it accounts for receiver activity status (including servicing).
  # We'll then update this, for each receiver, to define whether or not, if that receiver was active on a given date
  # ... which other receivers (if any) it overlapped in space (and time) with.
  rs_active_mat <- make_matrix_receivers(moorings = moorings,
                                         services = services,
                                         delta_t = "days",
                                         as_POSIXct = NULL)

  #### Define a list, with one dataframe element per receiver, that defines, for each time step, the overlapping receivers (0, 1)
  # Loop over each centroid...
  centroids_ls <- lapply(1:length(centroids), function(i) centroids[i, ])
  list_by_receiver <- pbapply::pblapply(centroids_ls, function(centroid) {

    #### Collect centroid and receiver status information
    # centroid <- centroids_ls[[2]]
    # Copy receiver activity status matrix
    info <- rs_active_mat
    # Focus on receiver's deployment window
    info <- info[as.Date(rownames(info)) %within% lubridate::interval(centroid$receiver_start_date, centroid$receiver_end_date), ]

    #### Convert receiver 'active' index (0, 1) to 'overlapping' index
    # ... A) Check for overlapping receivers
    # ... B) If there are overlapping receivers,
    # ... ... then we force all dates when the receiver of interest was not active to take 0 (no receivers could overlap with it then)
    # ... ... and we force all receivers that didn't overlap in space to 0
    # ... C) If there are no overlapping receivers, the whole matrix just gets forced to 0

    ## (A) Get an index of the receivers that intersected with the current receiver
    centroids_sbt <- centroids[!(centroids$receiver_id %in% centroid$receiver_id), ]
    int_1 <- rgeos::gIntersects(centroid, centroids_sbt, byid = TRUE)

    ## (B) If there are any overlapping receivers,
    if(any(int_1)){

      ## Process 'overlap' when the receiver was not active
      # ... Any time there is a '0' for activity status of the current receiver (e.g., due to servicing),
      # ... there can be no overlap with that receiver
      # ... so we will set a '0' to all other receivers
      # ... some of which may have been active on that date
      # ... Note the implementation of this step before the step below, when all rows for the receiver
      # ... of interest are (inadvertently) set to 0.
      info[which(info[, as.character(centroid$receiver_id)] == 0), ] <- 0

      ## Process 'overlap' for overlapping/non-overlapping receivers
      # For overlapping receivers, we'll leave these as defined in the activity matrix
      # ... (if active, then they overlap;
      # ... if not active, e.g., due to a servicing event for that receiver, then they can't overlap).
      # ... For the non-overlapping receivers, we'll force '0' for the overlap (even if they were active).
      # Get receiver IDs
      centroids_that_overlapped <- data.frame(centroids_sbt[which(int_1), ])
      # For all non-overlapping receivers, set '0' for overlap
      # ... Note that this will include the receiver of interest
      # ... But that doesn't matter because we'll drop that column anyway
      info[, !(colnames(info) %in% centroids_that_overlapped$receiver_id)] <- 0

      ## (C) If there aren't any spatially overlapping receivers, then the whole matrix just takes on 0
    } else  info[] <- 0

    #### Process dataframe
    rnms <- rownames(info)
    info <- data.frame(info)
    colnames(info) <- colnames(rs_active_mat)
    info[, as.character(centroid$receiver_id)] <- NULL
    cnms <- colnames(info)
    info$timestamp <- as.Date(rnms)
    info$receiver_id <- centroid$receiver_id
    info <- info[, c("timestamp", "receiver_id", cnms)]
    return(info)
  })
  names(list_by_receiver) <- as.character(centroids$receiver_id)

  #### On each date, get the vector of overlapping receivers
  # Note that not every receiver in this list will necessarily overlap with every other receiver though.
  lbd <- lapply(list_by_receiver, function(d){
    tidyr::pivot_longer(data = d,
                        cols = 3:ncol(d),
                        names_to = "receiver_id_2",
                        names_transform = list(receiver_id_2 = as.integer))
  })
  lbd <- dplyr::bind_rows(lbd) %>% dplyr::filter(.data$value == 1)
  lbd <- lapply(split(lbd, lbd$timestamp), function(d) unique(c(d$receiver_id[1], d$receiver_id_2)))

  ##### Process outputs
  # For the list_by_receiver, we will have one element for each receiver from 1:max(moorings$receiver_id)
  # ... (for each indexing)
  list_by_receiver <- lapply(as.integer(1:max(moorings$receiver_id)), function(i){
    if(i %in% moorings$receiver_id) return(list_by_receiver[[as.character(i)]]) else return(NULL)
  })
  # For the list_by_date (lbd), we will have one element for each date from the start to the end of the array
  list_by_date <- list()
  for(day in as.character(seq(min(moorings$receiver_start_date), max(moorings$receiver_end_date), "days"))){
    if(is.null(lbd[[day]])) list_by_date[[day]] <- NULL else list_by_date[[day]] <- lbd[[day]]
  }

  #### Return outputs
  out <- list()
  out$list_by_receiver <- list_by_receiver
  out$list_by_date     <- list_by_date
  return(out)
}


######################################
######################################
#### get_detection_area_sum()

#' @title Calculate the total area sampled by acoustic receivers
#' @description This function calculates the total area sampled by receivers, under the assumption of a constant detection range. To implement the function, receiver locations must be supplied as a SpatialPoints or SpatialPointsDataFrame object with the Universe Transverse Mercator coordinate reference system. The \code{\link[flapper]{get_detection_centroids}} is used to calculate the detection centroids around receivers, given a specified detection range (m) and any barriers to detection, such as coastline, and then the total area covered by receivers is calculated, accounting for overlapping centroids.
#'
#' @param xy,detection_range,coastline,plot,... Arguments required to calculate and visualise detection centroids via \code{\link[flapper]{get_detection_centroids}}; namely, receiver locations (\code{xy}), the detection range (\code{detection_range}), barriers to detection (\code{coastline}), and whether or not to plot the centroids (\code{plot}).
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
#' @param xy,detection_range,coastline,scale Arguments required to calculate the total area surveyed by receivers (at each time point) via \code{\link[flapper]{get_detection_area_sum}}. For this function, \code{xy} should be a \code{\link[sp]{SpatialPolygonsDataFrame-class}} object that includes both receiver locations and corresponding deployment times (in columns named `receiver_start_date' and `receiver_end_date' respectively).
#' @param plot A logical input that defines whether or not to plot a time series of the total area sampled by receivers.
#' @param verbose A logical input that defines whether or not to print messages to the console to relay function progress.
#' @param cl (optional) A cluster object created by \code{\link[parallel]{makeCluster}} to implement the algorithm in parallel. The connection to the cluster is closed within the function.
#' @param varlist (optional) A character vector of names of objects to export, to be passed to the \code{varlist} argument of \code{\link[parallel]{clusterExport}}. This may be required if \code{cl} is supplied. Exported objects must be located in the global environment.
#' @param ... Additional arguments, passed to \code{\link[prettyGraphics]{pretty_plot}}, to customise the plot produced.
#'
#' @return The function returns a dataframe with, for each date (`date') from the time of the first receiver's deployment to the time of the last receiver's retrieval, the number of receivers operational on that date (`n') and the total area sampled (`receiver_area'). If \code{plot = TRUE}, the function also returns a plot of the area sampled by receivers through time.
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
#'                              detection_range = 500,
#'                              coastline = dat_coast,
#'                              cl = parallel::makeCluster(2L),
#'                              varlist = "dat_coast"
#'                              )
#' }
#'
#' #### Example (3) Hide or customise the plot
#' dat <- get_detection_area_ts(xy, plot = FALSE)
#' dat <-
#'   get_detection_area_ts(xy,
#'                         pretty_axis_args =
#'                           list(
#'                             axis = list(list(format = "%b-%y"),
#'                                         list()
#'                                         )),
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
#' @return The function returns a dataframe that, for each time step (`time'), defines the number of operational units at that time (`n'). If \code{plot = TRUE}, the function also plots a time series of the number of operational units.
#'
#' @examples
#' #### Example (1): Number of operational receivers over an acoustic telemetry study
#' dat_n <- get_n_operational_ts(data = dat_moorings,
#'                               start = "receiver_start_date",
#'                               stop = "receiver_end_date")
#' utils::head(dat_n)
#'
#' #### Example (2): Number of individuals at liberty over a tagging study
#' # Define 'tag_end_date' as hypothetical end date of a study
#' # ...  and assume that all individuals remained tagged until this time
#' dat_ids$tag_end_date <- as.Date("2017-06-05")
#' dat_n <- get_n_operational_ts(data = dat_ids,
#'                               start = "tag_start_date",
#'                               stop = "tag_end_date")
#'
#' #### Example (3): Specify the time period under consideration
#' dat_n <- get_n_operational_ts(data = dat_ids,
#'                               start = "tag_start_date",
#'                               stop = "tag_end_date",
#'                               times = seq(min(dat_moorings$receiver_start_date),
#'                                           max(dat_moorings$receiver_end_date), 1)
#'                               )
#'
#' #### Example (4): Suppress or customise the plot
#' dat_n <- get_n_operational_ts(data = dat_ids,
#'                               start = "tag_start_date",
#'                               stop = "tag_end_date",
#'                               plot = FALSE)
#' dat_n <- get_n_operational_ts(data = dat_ids,
#'                               start = "tag_start_date",
#'                               stop = "tag_end_date",
#'                               xlab = "Time", ylab = "N (individuals)",
#'                               type = "l")
#'
#' #### Example (5): Additional examples with simulated data
#' # Example with one unit deployed on each day
#' tmp <- data.frame(id = 1:3L,
#'                   start = as.Date(c("2016-01-01", "2016-01-02", "2016-01-03")),
#'                   stop = as.Date(c("2016-01-01", "2016-01-02", "2016-01-03"))
#'                   )
#' get_n_operational_ts(data = tmp, start = "start", stop = "stop")
#' # Example with one unit deployed over a longer period
#' tmp <- data.frame(id = 1:3L,
#'                   start = as.Date(c("2016-01-01", "2016-01-02", "2016-01-03")),
#'                   stop = as.Date(c("2016-01-10", "2016-01-02", "2016-01-03"))
#'                   )
#' get_n_operational_ts(data = tmp, start = "start", stop = "stop")
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
#' @param ids A dataframe that defines individual deployment periods. This must contain a column that defines individual IDs (named `individual_id') and the time of tagging (named `tag_start_date') and time of tag retrieval (`tag_end_date') (see \code{\link[flapper]{dat_ids}} for an example).
#' @param moorings A dataframe that defines receiver deployment periods. This must contain a column that defines receiver IDs (named `receiver_id') and the time of receiver deployment (named `receiver_start_date') and retrieval (named `receiver_end_date') (see \code{\link[flapper]{dat_moorings}} for an example).
#' @param individual_id (optional) A vector of individuals for which to calculate overlap duration.
#' @param receiver_id (optional) A vector of receivers for which to calculate overlap duration.
#' @param type If both \code{individual_id} and \code{receiver_id} are specified, then \code{type} is an integer that defines whether or not to calculate overlap duration for (a) each individual/receiver pair (\code{type = 1L}) or (b) all combinations of specified individuals/receivers (\code{type = 2L}).
#' @param match_to (optional) A dataframe against which to match the calculated overlap duration(s). This must contain an `individual_id' and `receiver_id' column, as in \code{ids} and \code{moorings} respectively. If supplied, an integer vector of overlap durations for individual/receiver combinations, matched against the individuals/receivers in this dataframe, is returned (see also Value).
#' @param ... Additional arguments (none implemented).
#'
#' @return The function returns a dataframe with the deployment overlap duration for specific or all combinations of individuals and receivers, with the `individual_id', `receiver_id', `tag_start_date', `tag_end_date', `receiver_start_date' and `receiver_end_date' columns retained, plus `tag_interval' and `receiver_interval' columns that define individual and receiver deployment periods as \code{\link[lubridate]{Interval-class}} objects. The `id_rec_overlap' column defines the temporal overlap (days). Alternatively, if \code{match_to} is supplied, a vector of overlap durations that matches each individual/receiver observation in that dataframe is returned.
#'
#' @examples
#' #### Prepare data to include required columns
#' # moorings requires 'receiver_id', 'receiver_start_date' and 'receiver_end_date'
#' # ids requires 'individual_id', 'tag_start_date' and 'tag_end_date'
#' # These columns are already supplied in the example datasets
#' # ... except tag_end_date:
#' dat_ids$tag_end_date <- as.Date("2017-06-05")
#'
#' #### Example (1): Temporal overlap between all combinations
#' # ... of individuals and receivers
#' dat <- get_id_rec_overlap(dat_ids, dat_moorings)
#'
#' #### Example (2): Temporal overlap between all combinations of specified
#' #... individuals/receivers
#' dat <- get_id_rec_overlap(dat_ids,
#'                           dat_moorings,
#'                           individual_id = c(25, 26),
#'                           receiver_id = c(3, 4),
#'                           type = 2L)
#'
#' #### Example (3): Temporal overlap between specified individual/receiver pairs
#' dat <- get_id_rec_overlap(dat_ids,
#'                           dat_moorings,
#'                           individual_id = c(25, 26),
#'                           receiver_id = c(3, 4),
#'                           type = 1L)
#'
#' #### Example (4): Match temporal overlap to another dataframe
#' dat_acoustics$get_id_rec_overlap <-
#'   get_id_rec_overlap(dat_ids,
#'                      dat_moorings,
#'                      match_to = dat_acoustics,
#'                      type = 1L)
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
#' @param sample_size (optional) An integer that defines the number of samples of the environmental variable to draw from the area around each receiver (see the `size' argument of \code{\link[base]{sample}}). If this is provided, \code{sample_size} samples are taken from this area; otherwise, all values are extracted.
#' @param sample_replace (optional) If \code{sample_size} is specified, \code{sample_replace} is a logical input that defines whether to implement sampling with (\code{sample_replace = TRUE}, the default) or without (\code{sample_replace = FALSE}) replacement (see the `replace' argument of \code{\link[base]{sample}}).
#' @param sample_probs (optional) If \code{sample_size} is specified, \code{sample_probs} is a function that calculates the detection probability given the distance (m) between a cell and a receiver.
#' @param cl A cluster object created by \code{\link[parallel]{makeCluster}}. If supplied, the connection to the cluster is closed within the function.
#' @param varlist A character vector of names of objects to export, to be passed to the \code{varlist} argument of \code{\link[parallel]{clusterExport}}. This may be required if \code{cl} is supplied. Exported objects must be located in the global environment.
#' @param verbose A logical variable that defines whether or not relay messages to the console to monitor function progress.
#'
#' @return The function returns a list of dataframes (one for each element in \code{xy}; i.e., each receiver), each of which includes the cell IDs of \code{envir} from which values were extracted (`cell'), the value of the environmental variable in that cell (`envir') and, if applicable, the distance between that cell and the receiver (`dist', m) and the detection probability in that cell (`prob').
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
#'                                                     detection_range = 425,
#'                                                     coastline = dat_coast,
#'                                                     envir = dat_gebco
#'                                                     )
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
#'                                                    detection_range = 425,
#'                                                    coastline = dat_coast,
#'                                                    envir = dat_gebco,
#'                                                    sample_size = 2
#'                                                    )
#' utils::str(depths_by_centroid)
#'
#' #### Example (3) Extract a random sample of values with weighted probabilities
#' # Define detection probability function based only on distance
#' calc_detection_pr <-
#'   function(dist){
#'     dpr <- get_detection_pr(distance = dist,
#'                             beta_0 = 2.5,
#'                             beta_1 = -0.01,
#'                             inv_link = stats::plogis,
#'                             output = 2L)
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
#' @description The function calculates the total number of days (termed `detection days') during which individuals were detected at passive acoustic telemetry receivers. To implement the function, a dataframe with passive acoustic telemetry detections of individuals at receivers must be supplied. Detection days can be calculated for all combinations of individuals/receivers in this dataframe, for all combinations of specified individuals and receivers, or for specific individual/receiver pairs. The function returns a dataframe of detection days for these individual/receiver combinations or a vector of detection days that is matched against another dataframe.
#'
#' @param acoustics A dataframe that contains passive acoustic telemetry detection time series (see \code{\link[flapper]{dat_acoustics}} for an example). This should contain the following columns: a vector of individual IDs, named `individual_id'; a vector of receiver IDs, named `receiver_id'; and a POSIXct vector of time stamps when detections were made, named `timestamp'.
#' @param individual_id (optional) A vector of individuals for which to calculate detection days.
#' @param receiver_id (optional) A vector of receivers for which to calculate detection days.
#' @param type If both \code{individual_id} and \code{receiver_id} are specified, then \code{type} is an integer that defines whether or not to calculate detection days for (a) each individual/receiver pair (\code{type = 1L}) or (b) all combinations of individuals/receivers.
#' @param match_to (optional) A dataframe against which to match detection days. This must contain the `individual_id' and `receiver_id' column, as in \code{acoustics}. If supplied, an integer vector of detection days for individual/receiver combinations, matched against the individuals/receivers in \code{match_to}, is returned (see also Value).
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
#'                           individual_id = 25,
#'                           receiver_id = c(25, 21),
#'                           type = 2L)
#'
#' #### Example (3) Detection days between specified individual/receiver pairs
#' dat <- get_detection_days(dat_acoustics,
#'                           individual_id = c(25, 26),
#'                           receiver_id = c(25, 21),
#'                           type = 1L)
#'
#' #### Example (4) Match detection days to another dataframe
#' dat_acoustics$detection_days <- get_detection_days(dat_acoustics,
#'                                                    match_to = dat_acoustics,
#'                                                    type = 1L)
#'
#' utils::head(dat_acoustics)
#'
#' @seealso \code{\link[flapper]{get_detection_clump_lengths}}, \code{\link[flapper]{get_residents}}
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


######################################
######################################
#### get_detection_clump_lengths()

#' @title Get detection clump lengths
#' @importFrom rlang .data
#' @description This function gets the frequency distribution of the duration of unique, uninterrupted series of detections (`clumps') over a given time interval--for example, the number of occasions when detections clumped into periods of one, two, three, ..., n days in duration without gaps longer than one day in length.
#' @param acoustics A dataframe that contains passive acoustic telemetry detection time series (see \code{\link[flapper]{dat_acoustics}} for an example). At a minimum, this should contain a POSIXct vector of time stamps when detections were made named `timestamp'.
#' @param fct (optional) A character that defines the name of a column in \code{acoustics} that distinguishes acoustic time series for different individuals.
#' @param interval A character that defines the time interval. This must be supported by \code{\link[lubridate]{round_date}}. The default is \code{"days"}, in which case the function gets the frequency distribution of the duration of clump lengths in units of days, with detections within the same clump occurring within one day of each other.
#' @param ... Additional arguments (none implemented).
#'
#' @details Detection clumps provide a means to assess residency within an array by determining how long individuals tended to spend around receivers.
#'
#' @return The function returns a dataframe that defines the frequency distribution of the duration of detection clumps (i.e., the number of occasions when detection clumps lasted different lengths of time). This includes the following columns: (1) `n_interval', the duration of the clump, in units of \code{interval}; (2) `n_occasions', the number of occasions when detections clumped into windows of that duration; and, if \code{fct} is supplied, a column that distinguishes each individual.
#'
#' @examples
#' #### Define a hypothetical series of detections
#' # with one occasion when the 'chunk length' is 1 day
#' # ... two separate occasions when the chunk length is two days
#' # ... ... (i.e., the individual was detected on two consecutive day
#' # ... .... on two different occasions)
#' # ... one occasion when the chunk length is 5 days
#' # ... one occasion when the chunk lenghth is 7 days
#' eg <-
#'   data.frame(
#'     timestamp =
#'       as.Date(
#'         c("2016-01-01", # one week of continuous detections
#'           "2016-01-02",
#'           "2016-01-03",
#'           "2016-01-04",
#'           "2016-01-05",
#'           "2016-01-06",
#'           "2016-01-07",
#'           "2016-02-01", # one day with an isolated detection
#'           "2016-02-03", # two days with detections
#'           "2016-02-04",
#'           "2016-02-15", # another two days with detections
#'           "2016-02-16",
#'           "2016-03-01", # five days of continuous detections
#'           "2016-03-02",
#'           "2016-03-03",
#'           "2016-03-04",
#'           "2016-03-05")))
#'
#' #### Example (1): Implement function with default options
#' # ... (for one individual, with a daily time interval)
#' get_detection_clump_lengths(eg)
#'
#' #### Example (2): Implement function for multiple individuals
#' # Use a factor to distinguish IDs. As an example:
#' eg$individual_id <- 1L
#' get_detection_clump_lengths(eg, fct = "individual_id")
#'
#' #### Example (3): Change the time interval
#' ## E.g. Use an hourly interval:
#' # There are 17 unique clumps in this dataset, each comprising a single hour
#' eg$timestamp <- as.POSIXct(eg$timestamp)
#' get_detection_clump_lengths(eg, interval = "hours")
#' ## E.g. Use a monthly interval
#' # There is a single 'three-month' clump of detections for this individual
#' # ... when viewed at a monthly timescale:
#' get_detection_clump_lengths(eg, interval = "months")
#'
#' @seealso \code{\link[flapper]{get_detection_days}}, \code{\link[flapper]{get_residents}}
#' @author Edward Lavender
#' @export


get_detection_clump_lengths <-
  function(acoustics,
           fct = NULL,
           interval = "days",...){

    #### Checks
    check_names(input = acoustics, req = c("timestamp", fct), type = all)
    splitter <- fct
    if(is.null(fct)) {
      splitter <- "individual_id"
      acoustics[, splitter] <- 1L
    }

    #### Define a sequence of times at which to identify whether or not detections were recorded
    times <- seq(min(lubridate::floor_date(acoustics$timestamp, interval)),
                 max(lubridate::ceiling_date(acoustics$timestamp, interval)),
                 interval)
    acoustics$timestamp <- lubridate::round_date(acoustics$timestamp, interval)
    if(!identical(class(acoustics$timestamp), class(times)))
      stop("class(acoustics$timestamp) and class(times) are not identical.")

    #### Define a dataframe with .....
    counts_by_id <-
      # Loop over each individual
      lapply(split(acoustics, acoustics[, splitter]), function(d){
        # Identify whether detections were recorded at each timestamp
        bool <- times %in% d$timestamp
        # For each group of consecutive times when detections were made
        # ... count the number of consecutive times in that group. This
        # ... gives a dataframe with a unique identifier for each group of
        # ... consecutive detections and the number of times in that group
        # ... with detections (e.g., 1 day, 2 days etc.)
        tmp <- data.frame(grp = data.table::rleid(bool), bool = bool)
        tmp <-
          tmp %>%
          dplyr::filter(.data$bool) %>%
          dplyr::group_by(.data$grp) %>%
          dplyr::summarise(n = dplyr::n(), .groups = "drop_last")
        # Get the frequency distribution of the lengths of strings of consecutive detections
        # This shows the number of occasions when strings of consecutive detections
        # ... of a certain length were made
        # ... e.g., 15 occasions when the number of intervals in a string of consecutive detections
        # ... was only one, 5 occasions when consecutive detections lasted for two days (etc.)
        out <- table(tmp$n)
        out <- data.frame(n_intervals = as.integer(names(out)), n_occasions = as.numeric(out))
        out <- out %>% dplyr::arrange(.data$n_intervals)
        if(!is.null(fct)) out[, fct] <- d[1, fct]
        return(out)
      })
    counts <- do.call(rbind, counts_by_id)
    rownames(counts) <- NULL
    return(counts)
  }


########################################
########################################
#### get_detection_overlaps()

#' @title Get `overlapping' detections
#' @description This function isolates detections that occurred at `effectively the same time' at different receivers with overlapping or non-overlapping detection centroids. To implement the function, a dataframe of acoustic detections for a specific individual (\code{acoustics}) is required. Within this dataframe, the function isolates any detections that occurred within a user-specified time interval (\code{clock_drift}) at different receivers. If a list of the receivers with overlapping detection centroids in space and time is supplied (\code{overlaps}), the function also flags the subset of these detections that occurred at receivers with overlapping or non-overlapping detection centroids. This information is important for identifying potential false detections and the implementation of the AC* algorithm(s) (see Details). The function returns this information via summary message(s) along with a dataframe of the overlapping detections.
#'
#' @param acoustics A dataframe of passive acoustic telemetry detection time series (see \code{\link[flapper]{dat_acoustics}} for an example) for a single individual. This must contain an integer vector of receiver IDs, named `receiver_id' and a POSIXct vector of time stamps when detections were made, named `timestamp'.
#' @param overlaps (optional) A named list, from \code{\link[flapper]{get_detection_centroids_overlap}}, that defines, for each receiver, for each day over its deployment period, whether or not its detection centroid overlapped with those of other receivers.
#' @param clock_drift A number that defines the time (s) between sequential detections at which they are considered to have occurred at `effectively the same time'.
#' @param ... Additional arguments (none implemented).
#'
#' @details Detections at different receivers that occur at effectively the same time have important implications for inferences of animal movement patterns, especially via the AC* algorithm(s) in \code{\link[flapper]{flapper}} (e.g. \code{\link[flapper]{acdc}}). Within the AC* algorithm(s), when an individual is detected at the same time at two different receivers, detection probability kernels dictate that these receivers must have overlapping detection centroids. If the detection centroids do not overlap, this suggests that either (a) receiver clocks are not well-aligned; (b) the definition of `effectively the same time' is be overly large (such that the individual could move from within the detection centroid of one receiver into the detection centroid of another); (c) detection centroids are too small and detection probability is higher than realised; and/or (d) one or more of the detections are false. The most likely cause may be guided by knowledge of the array design, detection probability and false detection algorithms. For example, if it is plausible that detection centroids are too small, repeating the implementation of this function with larger centroids may indicate whether or not this is likely to have been the case: if so, all detections that occurred at effectively the same time at receivers with non-overlapping detection centroids should be captured by the large centroids (though this does not rule out other explanations). The purpose of this function is to flag any such detections so that they can be investigated and addressed prior to the implementation of the AC* algorithm(s).
#'
#' @return The function returns a message that defines the number of detections at different receivers that occurred at effectively the same time (within \code{clock_drift}) and, if \code{overlaps} is supplied, the subset of these that occurred at non-overlapping receivers. A dataframe is also invisibly returned that records the details of overlapping detections. This includes the following columns:
#' \itemize{
#'    \item \code{timestamp_1}, a POSIXct vector of time stamps at which detections were made;
#'    \item \code{receiver_id_1}, an integer identifier of the receivers at which detections were made;
#'    \item \code{timestamp_2}, a POSIXct vector of the time stamps at which immediately subsequent detections were made;
#'    \item \code{receiver_id_2}, an integer identifier of the receivers at which the immediately subsequent detections were made;
#'    \item \code{diff_time}, a numeric vector that defines the time (s) between consecutive detections (\code{timestamp_1} and \code{timestamp_2});
#'    \item \code{detection_in_overlapping_centroid}, a binary vector that defines whether (1) or not (0) the detection centroids of \code{receiver_id_1} and \code{receiver_id_2} overlapped at the time of the detection (this is only included if \code{overlaps} is provided);
#' }
#'
#' @examples
#' #### Example (1): Implement function using detection data for an example individual
#' dat_acoustics_25 <- dat_acoustics[dat_acoustics$individual_id == 25, ]
#' dat <- get_detection_overlaps(acoustics = dat_acoustics_25)
#' utils::head(dat)
#'
#' #### Example (2): Implement function, including information on detection centroids
#'
#' ## Get detection centroid overlaps to include in function
#' ## (see ?flapper::get_detection_centroid_overlaps)
#' # Define receiver locations
#' proj_wgs84 <- sp::CRS("+init=epsg:4326")
#' proj_utm <- sp::CRS(paste("+proj=utm +zone=29 +datum=WGS84",
#'                           "+units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))
#' rownames(dat_moorings) <- dat_moorings$receiver_id
#' xy <- sp::SpatialPoints(dat_moorings[, c("receiver_long", "receiver_lat")],
#'                         proj_wgs84)
#' xy <- sp::spTransform(xy, proj_utm)
#' # Get receiver-specific detection centroids
#' centroids <- get_detection_centroids(xy, byid = TRUE)
#' centroids_df <- dat_moorings[, c("receiver_id",
#'                                  "receiver_start_date",
#'                                  "receiver_end_date")]
#' row.names(centroids_df) <- names(centroids)
#' centroids <- sp::SpatialPolygonsDataFrame(centroids, centroids_df)
#' # Define detection centroid overlaps
#' overlaps <- get_detection_centroids_overlap(centroids = centroids)
#'
#' ## Implement function with detection centroid overlaps included
#' dat <- get_detection_overlaps(acoustics = dat_acoustics_25, overlaps = overlaps)
#' utils::head(dat)
#'
#' #### Implement function across all individuals
#' # For some individuals, there are simultaneous detections at receivers with
#' # ... non overlapping detection centroids, suggesting these are probably too small.
#' dat_by_id <-
#'   lapply(split(dat_acoustics, dat_acoustics$individual_id),
#'          function(acc_for_id){
#'            print(paste("individual_id", acc_for_id$individual_id[1], "-------"))
#'            get_detection_overlaps(acoustics = acc_for_id, overlaps = overlaps)
#'          })
#' ## Test this hypothesis by re-implementing approach with larger centroids
#' # Re-define centroids and receiver overlap
#' centroids <- get_detection_centroids(xy, detection_range = 750, byid = TRUE)
#' centroids <- sp::SpatialPolygonsDataFrame(centroids, centroids_df)
#' overlaps <- get_detection_centroids_overlap(centroids = centroids)
#' # Re-implement algorithm
#' dat_by_id <-
#'   lapply(split(dat_acoustics, dat_acoustics$individual_id),
#'          function(acc_for_id){
#'            print(paste("individual_id", acc_for_id$individual_id[1], "-------"))
#'            get_detection_overlaps(acoustics = acc_for_id, overlaps = overlaps)
#'          })
#' # There are now no observations within clock_drift at receivers with
#' # ... non-overlapping centroids.
#'
#' @author Edward Lavender
#' @export
#'

get_detection_overlaps <- function(acoustics, overlaps = NULL, clock_drift = 5,...){

  #### Checks
  check_names(input = acoustics, req = c("timestamp", "receiver_id"), type = all)
  if(!is.null(overlaps)) {
    check_names(input = overlaps, req = "list_by_receiver", type = all)
    overlaps <- overlaps$list_by_receiver
  }

  #### Process acoustics
  # Select columns
  acoustics$timestamp_1 <- acoustics$timestamp
  acoustics$timestamp <- NULL
  acoustics$receiver_id_1 <- acoustics$receiver_id
  acoustics <- acoustics[, c("timestamp_1", "receiver_id_1")]
  # Get original order
  acoustics$index <- 1:nrow(acoustics)
  if(is.unsorted(acoustics$timestamp_1)) {
    message("'acoustics' is not sorted by timestamp: are there detections for multiple individuals (there shouldn't be)? Sorting 'acoustics' by timestamp...")
    acoustics <- acoustics[order(acoustics$timestamp_1), ]
  }
  # Calculate the duration between sequential detections
  acoustics$receiver_id_2 <- dplyr::lead(acoustics$receiver_id_1)
  acoustics$timestamp_2 <- dplyr::lead(acoustics$timestamp_1)
  acoustics <- acoustics[1:(nrow(acoustics) - 1), ]
  acoustics$diff_time <- as.numeric(difftime(acoustics$timestamp_2, acoustics$timestamp_1, units = "secs"))
  # Filter observations that occurred within 'clock_drift'
  acoustics <- acoustics[acoustics$receiver_id_1 != acoustics$receiver_id_2 &
                           acoustics$diff_time <= clock_drift, ]

  #### Determine the number of observations that occurred within the clock drift
  n <- nrow(acoustics)
  message(n, " observation(s) identified at another receiver within ", clock_drift, " secs.")

  #### Examine which of these occurred at receivers with overlapping detection centroids
  # For each receiver, for each date, we will identify whether detections within the clock drift
  # ... occurred at a spatially overlapping receiver.
  if(n > 0 & !is.null(overlaps)) {
    # receiver_id_2 as a character
    acoustics$receiver_id_2_char <- as.character(acoustics$receiver_id_2)
    # Define a column to distinguish detection dates
    acoustics$timestamp_date <- as.Date(acoustics$timestamp_1)
    # Define blank column to store whether or not detections occurred at an overlapping receiver
    acoustics$detection_in_overlapping_centroid <- 0
    # Update acoustics with information on whether or not detections occurred at an overlapping receiver
    acc_by_receiver <-
      # For each receiver...
      lapply(split(acoustics, acoustics$receiver_id_1), function(acc_by_receiver){
        # Get  receiver-specific overlaps matrix
        # acc_by_receiver <- split(acoustics, acoustics$receiver_id_1)[[1]]
        overlap_for_receiver <- overlaps[[acc_by_receiver$receiver_id_1[1]]]
        stopifnot(acc_by_receiver$receiver_id[1] == overlap_for_receiver$receiver_id[1])
        # For each date...
        acc_by_receiver_by_date <-
          lapply(split(acc_by_receiver, acc_by_receiver$timestamp_date), function(acc_by_receiver_on_date){
            # Get date-specific overlaps matrix
            # acc_by_receiver_on_date <- split(acc_by_receiver, acc_by_receiver$timestamp_date)[[1]]
            overlap_for_receiver_on_date <-
              overlap_for_receiver[which(overlap_for_receiver$timestamp %in% acc_by_receiver_on_date$timestamp_date[1]),
                                   , drop = FALSE]
            # For each detection...
            for(i in 1:nrow(acc_by_receiver_on_date)){
              # Work out whether or not that detection occurred at an overlapping receiver
              acc_by_receiver_on_date$detection_in_overlapping_centroid[i] <-
                overlap_for_receiver_on_date[, acc_by_receiver_on_date$receiver_id_2_char[i]] == 1
            }
            return(acc_by_receiver_on_date)
          })
        acc_by_receiver <- do.call(rbind, acc_by_receiver_by_date)
        return(acc_by_receiver)
      })
    acoustics <- do.call(rbind, acc_by_receiver)
    acoustics$timestamp_date <- NULL

    #### Determine the number of detections at non-overlapping receivers
    n_at_non_overlapping_receivers <- length(which(acoustics$detection_in_overlapping_centroid == 0))
    message("Of these, there are ", n_at_non_overlapping_receivers, " observation(s) within ", clock_drift, " secs that are not in overlapping centroids.")
    # Examine the receiver combinations at which simultaneous detections at non-overlapping receivers occurred:
    if(n_at_non_overlapping_receivers > 0){
      acoustics$key <- paste0(acoustics$receiver_id_1, "-", acoustics$receiver_id_2)
      keys <- do.call(rbind, strsplit(unique(acoustics$key), split = "-", fixed = TRUE))
      keys <- data.frame(keys)
      keys$key <- NA
      for(i in 1:nrow(keys)){
        keys$key[i] <- paste0(sort(c(keys[i, 1], keys[i, 2])), collapse = "-")
      }
      keys <- keys[!duplicated(keys$key), ]
      message("These simulataneous detections at non-overlapping receivers occurred at ", nrow(keys), " unique receiver pair(s): ", paste0(keys$key, collapse = ", "), ".")
      acoustics$key <- NULL
    }
  }

  #### Return acoustics
  acoustics <- acoustics[order(acoustics$index), ]
  acoustics$index <- NULL
  acoustics$receiver_id_2_char <- NULL
  rownames(acoustics) <- NULL
  return(invisible(acoustics))
}


######################################
######################################
#### get_residents()

#' @title Identify `resident' individuals in a passive acoustic telemetry array
#' @importFrom rlang .data
#' @description This function identifies `resident' individuals in a passive acoustic telemetry array from acoustic detection time series.
#'
#' @param acoustics A dataframe that contains passive acoustic telemetry detection time series (see \code{\link[flapper]{dat_acoustics}} for an example). At a minimum, this should contain a POSIXct vector of time stamps when detections were made named `timestamp'.
#' @param fct (optional) A character that defines the name of a column in \code{acoustics} that distinguishes acoustic time series for different individuals.
#' @param resident_threshold_gap A number that specifies the maximum gap (days) between detections for an individual to be considered `resident' over the period of detections.
#' @param resident_threshold_duration A sorted, numeric vector that specifies the duration (days) over which `regular' detections (spaced less than \code{resident_threshold_gap} apart) must occur for an individual to be considered as `resident'. A vector can be supplied to distinguish multiple residency categories (e.g., short-term residents and long-term residents).
#' @param resident_labels A character vector of labels, one for each \code{resident_threshold_duration}, for residency categories (e.g., short-term and long-term residents).
#' @param non_resident_label A character that defines the label given to non-residents.
#' @param keep A character vector of column names in \code{acoustics} to retain in the returned dataframe (e.g., `sex', `maturity' etc.).
#'
#' @details This function assigns residency categories (e.g., non-resident, short-term resident, long-term resident) to individuals in \code{acoustics}, returning a dataframe with one row for each individual and a column for assigned residency categories. `Resident' individuals are defined as individuals with (a) `regular' detections (i.e. detections spaced less than some threshold time (\code{resident_threshold_duration}) apart) that (b) last for a requisite time interval (\code{resident_threshold_duration}). Multiple residency categories can be specified. Under the default options, non-residents are labelled \code{"N"}. Short-term residents (labelled \code{"S"}) are defined as individuals for which sequential detections were less than 31 days apart for more than 91 days (three months). Long-term residents (labelled \code{"L"}) are defined as individuals for which sequential detections were less than 31 days apart for more than 365 days (12 months).
#'
#' A \code{method} argument may be added in due course to implement other methods to define `resident' individuals (e.g., based on detection days, see \code{\link[flapper]{get_detection_days}}) .
#'
#' @return The function returns a dataframe that defines, for each individual (\code{fct}), the duration (days) between the first detection and the last detection in a series of detections spaced less than \code{resident_threshold_gap} days apart (`time'), and the residency category (as defined by \code{resident_threshold_duration} and \code{resident_labels}).
#'
#' @examples
#' #### Example (1): Implement default options
#' get_residents(dat_acoustics, fct = "individual_id")
#' # Compare results to abacus plot
#' prettyGraphics::pretty_plot(dat_acoustics$timestamp, dat_acoustics$individual_id)
#'
#' @seealso \code{\link[flapper]{get_detection_days}}, \code{\link[flapper]{get_detection_clump_lengths}}
#' @references Lavender, et al. (in review). Movement patterns of a Critically Endangered elasmobranch (\emph{Dipturus intermedius}) in a Marine Protected Area. Aquatic Conservation: Marine and Freshwater Ecosystems.
#'
#' @author Edward Lavender
#' @export

get_residents <- function(acoustics,
                          fct = NULL,
                          resident_threshold_gap = 31,
                          resident_threshold_duration = c(91.2501, 365.25),
                          resident_labels = c("S", "L"),
                          non_resident_label = "N",
                          keep = NULL){
  #### Checks
  check_names(input = acoustics, req = c("timestamp", fct), type = all)
  if(is.unsorted(resident_threshold_duration))
    stop("'resident_threshold_duration' should be a sorted, numeric vector.")
  if(length(resident_threshold_duration) != length(resident_labels))
    stop("The number of resident_threshold_duration(s) and resident_label(s) should be the same.")

  #### Split based on individual ID
  splitter <- fct
  if(is.null(fct)) {
    splitter <- "individual_id"
    acoustics[, splitter] <- 1L
  }

  #### Define a dataframe with residence times
  residents <-
    acoustics %>%
    dplyr::group_by(.data[[splitter]]) %>%
    dplyr::arrange(.data$timestamp) %>%
    dplyr::mutate(gap = Tools4ETS::serial_difference(.data$timestamp, units = "days"),
                  flag = Tools4ETS::flag_ts(.data$timestamp, duration_threshold = 24 * 60 * resident_threshold_gap, flag = 3)$flag3,
                  key = paste0(.data[[splitter]], "-", .data$flag)
    ) %>%
    dplyr::group_by(.data$key) %>%
    dplyr::mutate(time = as.numeric(difftime(max(.data$timestamp), min(.data$timestamp), units = "days"))) %>%
    dplyr::slice(1L) %>%
    dplyr::group_by(.data[[splitter]]) %>%
    dplyr::mutate(time = max(.data$time)) %>%
    dplyr::slice(1L) %>%
    dplyr::mutate(resident = factor(non_resident_label, levels = sort(c(non_resident_label, resident_labels)))) %>%
    data.frame()
  # Drop excess columns
  residents <- residents[,  c(splitter, keep, "time", "resident")]
  if(is.null(fct)) residents[, splitter] <- NULL

  #### Assign residency categories
  for(i in 1:length(resident_threshold_duration)){
    residents$resident[residents$time >= resident_threshold_duration[i]] <- resident_labels[i]
  }

  #### Return dataframe
  return(residents)
}

