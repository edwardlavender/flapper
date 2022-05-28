#######################################
#######################################
#### get_mvt_mobility_*()

#' @title Estimate individual swimming speeds from acoustic and archival data
#' @description These functions are designed to provide `ballpark' estimates of individual swimming speeds from (a) acoustic detections at passive acoustic telemetry receivers (\code{\link[flapper]{get_mvt_mobility_from_acoustics}}) and archival (depth) time series (\code{\link[flapper]{get_mvt_mobility_from_archival}}).
#'
#' @param data A dataframe of animal movement time series used to estimate swimming speeds. For \code{\link[flapper]{get_mvt_mobility_from_acoustics}}, \code{data} should contain passive acoustic telemetry detection time series (see \code{\link[flapper]{dat_acoustics}} for an example). This must contain a vector of receiver IDs, named `receiver_id', and a POSIXct vector of time stamps when detections were made, named `timestamp'. For \code{\link[flapper]{get_mvt_mobility_from_archival}}, \code{data} is a dataframe that contains depth time series (see \code{\link[flapper]{dat_archival}} for an example). This must contain a numeric vector of observed depths (m), named `depth', and a POSIXct vector of regular time stamps when observations were made, named `timestamp'. In either case, an additional column can be included to distinguish time series for different individuals (see \code{fct}).
#' @param fct (optional) A character variable that defines the name of a column in \code{data} that distinguishes levels of a grouping variable (e.g., individuals).
#' @param moorings In \code{\link[flapper]{get_mvt_mobility_from_acoustics}}, \code{moorings} is \code{\link[sp]{SpatialPointsDataFrame}} that defines receiver locations (in the Universe Transverse Mercator coordinate reference system) and receiver IDs (in a column named `receiver_id').
#' @param detection_range In \code{\link[flapper]{get_mvt_mobility_from_acoustics}}, \code{detection_range} is a number that defines the detection range (m) of receivers.
#' @param calc_distance In \code{\link[flapper]{get_mvt_mobility_from_acoustics}}, \code{calc_distance} is a character that defines the method used to calculate distances between receivers. Currently supported options are Euclidean distances (\code{"euclid"}) or least-cost (shortest) path distances ("lcp") over a surface (see \code{bathy}).
#' @param bathy In \code{\link[flapper]{get_mvt_mobility_from_acoustics}}, if \code{calc_distance = "lcp"}, \code{bathy} is \code{\link[raster]{raster}} that defines the surface over which individual(s) moved. \code{bathy} must be planar (i.e., Universal Transverse Mercator projection) with units of metres in x, y and z directions (m). The surface's resolution is taken to define the distance between horizontally and vertically connected cells and must be the same in both x and y directions. Any cells with NA values (e.g., due to missing data) are treated as `impossible' to move though by the algorithm (see \code{\link[flapper]{lcp_over_surface}}).
#' @param transmission_interval (optional) In \code{\link[flapper]{get_mvt_mobility_from_acoustics}}, \code{transmission_interval} is the (maximum) time delay between sequential acoustic transmissions.
#' @param step (optional) In \code{\link[flapper]{get_mvt_mobility_from_acoustics}}, \code{step} is a number that defines the time (s) interval over which to represent speeds. By default, speeds are calculated in m/s. If \code{step} is supplied, speeds are also calculated in m/step s, which can be easier to interpret. For \code{\link[flapper]{get_mvt_mobility_from_archival}}, `step' is set automatically as the duration (s) between sequential depth observations, which is assumed to remain constant through time.
#'
#' @details
#'
#' \subsection{Acoustic estimates}{Speed estimates from passive acoustic telemetry time series are derived from examination of detection patterns via \code{\link[flapper]{get_mvt_mobility_from_acoustics}}. For each \code{fct} group (e.g., individual), the function identifies sequential detections, which exceeded the maximum transmission interval, at receivers with non-overlapping detection containers (areas within the detection range of each receiver). Assuming that all detections are true, these can only result from movement between receivers. For these movements, transition distances are calculated as Euclidean distances, via \code{\link[raster]{pointDistance}}, or shortest swimming distances, assuming movement over a surface, via \code{\link[flapper]{lcp_over_surface}}. Given that detections can arise from movement anywhere within the detection range of a receiver, three distances are calculated: an average distance, from receiver to receiver; a lower-bound distance, between the closest edges of receiver detection containers (i.e., the average distance minus two times the detection range); and an upper-bound distance from the furthest edges of receiver detection ranges (i.e., the average distance plus two times the detection range). These estimates assume a uniform detection ranges over space. For each transition, distances are converted into an average, lower and upper speed (m/s) estimate (termed `speed_avg_ms`, `speed_min_ms` and `speed_max_ms` respectively). If \code{step} is supplied, speeds are also expressed per step. On many occasions, individuals will take indirect routes between receivers, resulting in inappropriately low speed estimates. However, if receivers are sufficiently close together such that individuals sometimes effectively transition directly between receivers, the faster speed estimates derived via this method may be quite informative about actual swimming speeds.}
#'
#' \subsection{Archival estimates}{Speed estimates from \code{archival} time series are derived from examination of changes in depth through time (vertical activity) via \code{\link[flapper]{get_mvt_mobility_from_archival}}. For each individual, speed is calculated from the vertical distances (m) between sequential, regular depth (m) observations over time (s). Speeds are also expressed per step, where `step' is the duration (s) between sequential observations.}
#'
#' @return Both functions print a statistical summary of the speed estimates to the console and return a dataframe with observations and the corresponding speed estimates invisibly.
#'
#' For \code{\link[flapper]{get_mvt_mobility_from_acoustics}}, for each speed estimate, the minimum, mean and maximum speeds are printed, in units of m/s and, if specified, m/step. The returned dataframe defines all transitions between receivers and includes the following columns:
#' \itemize{
#'   \item \code{fct} is the grouping factor (if specified);
#'   \item \code{receiver_id_1} and \code{receiver_id_2} are the receiver IDs for each transition between receivers;
#'   \item \code{timestamp_1} and \code{timestamp_2} are the time stamps of detections for each transition between receivers;
#'   \item \code{time} is the duration (s) between \code{timestamp_1} and \code{timestamp_2};
#'   \item \code{dist_min}, \code{dist_avg} and \code{dist_max} are the estimates for the distance (m) travelled to transition between receivers;
#'   \item \code{speed_min_ms}, \code{speed_avg_ms} and \code{speed_min_ms} are the corresponding speed (m/s) estimates;
#'   \item \code{speed_min_mstep}, \code{speed_avg_mstep} and \code{speed_max_mstep} are the corresponding speed estimates, in m per step (if specified);
#' }
#'
#' For \code{\link[flapper]{get_mvt_mobility_from_archival}}, the function prints a statistical summary of estimated speeds in m/s and m per `step', where `step' is the duration (s) between sequential depth observations. The function also invisibly returns dataframe of the estimates. This is as inputted, without any \code{fct} levels with fewer than two observations and with the following additional columns:
#' \itemize{
#'   \item \code{dist} is the distance (m) between sequential depth observations;
#'   \item \code{speed_ms} is the speed (m/s) of vertical movement between sequential depth observations;
#'   \item \code{speed_m_per_step} is the speed of vertical movement between sequential depth observations in units of m per step, where `step' is the duration (s) between sequential depth observations;
#' }
#'
#' @examples
#' #### Estimate mobility from acoustic data using Euclidean distances
#' ## (A) Define receiver coordinates as SPDF in UTM CRS
#' proj     <- sp::CRS(SRS_string = "EPSG:4326")
#' proj_utm <- sp::CRS(SRS_string = "EPSG:32629")
#' moorings <- sp::SpatialPoints(dat_moorings[, c("receiver_long", "receiver_lat")], proj)
#' moorings <- sp::spTransform(moorings, proj_utm)
#' moorings <- sp::SpatialPointsDataFrame(moorings, data.frame(receiver_id = dat_moorings$receiver_id))
#' ## (B) Implement algorithm using Euclidean distances
#' mob_1 <- get_mvt_mobility_from_acoustics(data = dat_acoustics,
#'                                          fct = "individual_id",
#'                                          moorings = moorings,
#'                                          detection_range = 750,
#'                                          transmission_interval = 120,
#'                                          step = 120)
#'
#' #### Estimate mobility from acoustic data using LCP distances
#' ## (A) Define receiver coordinates
#' # ... (as above)
#' ## (B) Define bathymetry surface for LCP calculations
#' # Requirements: mask land; planar UTM projection; equal resolution;
#' bathy <- raster::raster(ext = raster::extent(dat_gebco), resolution = 50)
#' bathy <- raster::resample(dat_gebco, bathy, method = "bilinear")
#' ## (C) Implement algorithm using LCP distances
#' mob_2 <- get_mvt_mobility_from_acoustics(data = dat_acoustics,
#'                                          fct = "individual_id",
#'                                          moorings = moorings,
#'                                          detection_range = 750,
#'                                          calc_distance = "lcp",
#'                                          bathy = bathy,
#'                                          transmission_interval = 120,
#'                                          step = 120)
#'
#' #### Estimate mobility from archival data
#' # Note the use of 'individual_id' for 'fct' here is only appropriate
#' # ... when there are no breaks in the time series.
#' mob_3 <- get_mvt_mobility_from_archival(dat_archival, fct = "individual_id")
#'
#' @author Edward Lavender
#' @name get_mvt_mobility
NULL


#### get_mvt_mobility_from_acoustics()
#' @rdname get_mvt_mobility
#' @export

get_mvt_mobility_from_acoustics <- function(data,
                                            fct = NULL,
                                            moorings,
                                            detection_range,
                                            calc_distance = c("euclid", "lcp"),
                                            bathy = NULL,
                                            transmission_interval,
                                            step = NULL){

  #### Checks
  check_names(input = data, req = c("receiver_id", "timestamp", fct), type = all)
  if(!is.null(fct)) data$fct <- data[, fct] else data$fct <- 1L
  moorings_data <- data.frame(moorings)
  check_names(input = moorings_data, req = "receiver_id", type = all)
  if(!all(unique(data$receiver_id) %in% moorings_data$receiver_id)) stop("Not all unique data$receiver_id are in moorings$receiver_id.")
  calc_distance <- match.arg(calc_distance)
  if(calc_distance == "lcp" & is.null(bathy)) stop("'bathy' is required for calc_distance = 'lcp'.")
  msg_no_suitable_data <- "No 'data' following filtering with which to estimate mobility."

  #### Add receiver coordinates to data
  coords <- sp::coordinates(moorings)
  colnames(coords)   <- c("easting", "northing")
  match_index        <- match(data$receiver_id, moorings$receiver_id)
  data$easting  <- coords[match_index, "easting"]
  data$northing <- coords[match_index, "northing"]

  #### Define a dataframe of transitions between receivers
  # These are the observations we will use to calculate speeds
  # Note that we will only use detections separated by more than transmission_interval
  # ... because some detections very close together in time tend to occur
  # ... at receivers with overlapping detection ranges
  # ... which causes issues for distance calculations, given we do not know
  # ... the exact speed of the animals.
  transitions <-
    data %>%
    dplyr::select(.data$fct, .data$timestamp, .data$receiver_id, .data$easting, .data$northing) %>%
    dplyr::group_by(.data$fct) %>%
    dplyr::arrange(.data$timestamp) %>%
    dplyr::mutate(r_1 = .data$receiver_id,
                  r_2 = dplyr::lead(.data$receiver_id),
                  r_key = paste0(.data$r_1, "-", .data$r_2),
                  t_1 = .data$timestamp,
                  t_2 = dplyr::lead(.data$timestamp),
                  time = as.numeric(difftime(.data$t_2, .data$t_1, units = "s")),
                  easting_1 = .data$easting,
                  northing_1 = .data$northing,
                  easting_2 = dplyr::lead(.data$easting),
                  northing_2 = dplyr::lead(.data$northing)) %>%
    dplyr::select(-c(.data$timestamp, .data$receiver_id, .data$easting, .data$northing)) %>%
    dplyr::filter(.data$r_1 != .data$r_2) %>%
    dplyr::filter(.data$time >= transmission_interval) %>%
    dplyr::filter(!is.na(.data$r_2))
  if(nrow(transitions) == 0L){
    message(msg_no_suitable_data)
    return(invisible())
  }

  #### Define unique receiver transitions
  # We will use these to calculate distances
  receiver_pairwise_dists <- transitions %>%
    dplyr::ungroup() %>%
    dplyr::group_by(.data$r_key) %>%
    dplyr::slice(1L) %>%
    dplyr::select(.data$r_1, .data$r_2, .data$r_key, .data$easting_1, .data$northing_1, .data$easting_2, .data$northing_2)

  #### Get distances between receivers
  # Calculate distances between receivers using Euclidean/LCP distances
  # Then process distances to account for detection radii
  # ... Minus detection_range * 2 m from the distance between receivers to get a adjusted distance
  # ... Then drop any detections at receivers that are within each other's distance
  ## Define coordinates
  origins      <- receiver_pairwise_dists[, c("easting_1", "northing_1")] %>% as.matrix()
  destinations <- receiver_pairwise_dists[, c("easting_2", "northing_2")] %>% as.matrix()
  ## Get distances
  if(calc_distance == "euclid"){
    receiver_pairwise_dists$dist_avg <- raster::pointDistance(origins, destinations, lonlat = FALSE)
  } else {
    receiver_pairwise_dists_lcps <- lcp_over_surface(origin = origins,
                                                     destination = destinations,
                                                     surface = bathy,
                                                     goal = 3,
                                                     combination = "pair",
                                                     verbose = FALSE)
    receiver_pairwise_dists$dist_avg <- receiver_pairwise_dists_lcps$dist_lcp
  }

  ## Adjust distances and drop spatially overlapping receivers
  receiver_pairwise_dists <-
    receiver_pairwise_dists %>%
    dplyr::mutate(dist_min = .data$dist_avg - detection_range * 2,
                  dist_max = .data$dist_avg + detection_range * 2) %>%
    dplyr::filter(.data$dist_min > 0)
  if(nrow(receiver_pairwise_dists) == 0L){
    message(msg_no_suitable_data)
    return(invisible())
  }

  #### Update transitions
  ## Filter transitions to include transitions between appropriate receivers
  # ... (whose adjusted distance is more than 0 m!)
  transitions <-
    transitions %>% dplyr::filter(.data$r_key %in% receiver_pairwise_dists$r_key)
  if(nrow(transitions) == 0L){
    message(msg_no_suitable_data)
    return(invisible())
  }
  ## Add distances
  match_index <- match(transitions$r_key, receiver_pairwise_dists$r_key)
  transitions$dist_min <- receiver_pairwise_dists$dist_min[match_index]
  transitions$dist_avg <- receiver_pairwise_dists$dist_avg[match_index]
  transitions$dist_max <- receiver_pairwise_dists$dist_max[match_index]

  #### Tidy transitions
  transitions <- transitions %>%
    dplyr::rename(receiver_id_1 = .data$r_1,
                  receiver_id_2 = .data$r_2,
                  timestamp_1 = .data$t_1,
                  timestamp_2 = .data$t_2) %>%
    dplyr::select(.data$fct,
                  .data$receiver_id_1, .data$receiver_id_2,
                  .data$timestamp_1, .data$timestamp_2,
                  .data$time, .data$dist_min, .data$dist_avg, .data$dist_max
                  )
  if(is.null(fct)){
    transitions$fct <- NULL
  } else {
    cnms <- colnames(transitions)
    cnms[which(colnames(transitions) == "fct")] <- fct
    colnames(transitions) <- cnms
  }

  #### Calculate speeds
  message("Mobility estimates from n = ", nrow(transitions), " observation(s).")
  ## Speeds (m/s)
  transitions$speed_min_ms <- transitions$dist_min/transitions$time
  transitions$speed_avg_ms <- transitions$dist_avg/transitions$time
  transitions$speed_max_ms <- transitions$dist_max/transitions$time
  cat("--------------------------------------\n")
  cat("Estimates (m/s)-----------------------\n")
  stats <- data.frame(variable = c("speed_min_ms", "speed_avg_ms", "speed_max_ms"),
                      min = c(min(transitions$speed_min_ms), min(transitions$speed_avg_ms), min(transitions$speed_max_ms)),
                      mean = c(mean(transitions$speed_min_ms), mean(transitions$speed_avg_ms), mean(transitions$speed_max_ms)),
                      max = c(max(transitions$speed_min_ms), max(transitions$speed_avg_ms), max(transitions$speed_max_ms))
                      )
  stats$min  <- round(stats$min, 2)
  stats$mean <- round(stats$mean, 2)
  stats$max <- round(stats$max, 2)
  print(stats)
  cat("--------------------------------------\n")
  ## Speeds (m per step)
  if(!is.null(step)) {
    transitions$speed_min_mstep <- transitions$speed_min_ms * step
    transitions$speed_avg_mstep <- transitions$speed_avg_ms * step
    transitions$speed_max_mstep <- transitions$speed_max_ms * step
    cat("Estimates (m/step)--------------------\n")
    stats <- data.frame(variable = c("speed_min_mstep", "speed_avg_mstep", "speed_max_mstep"),
                        min = c(min(transitions$speed_min_mstep), min(transitions$speed_avg_mstep), min(transitions$speed_max_mstep)),
                        mean = c(mean(transitions$speed_min_mstep), mean(transitions$speed_avg_mstep), mean(transitions$speed_max_mstep)),
                        max = c(max(transitions$speed_min_mstep), max(transitions$speed_avg_mstep), max(transitions$speed_max_mstep))
                        )
    stats$min  <- round(stats$min, 2)
    stats$mean <- round(stats$mean, 2)
    stats$max <- round(stats$max, 2)
    print(stats)
    cat("--------------------------------------\n")
  }

  #### Return outputs
  return(invisible(transitions))
}


#### get_mvt_mobility_from_archival()
#' @rdname get_mvt_mobility
#' @export

get_mvt_mobility_from_archival <- function(data, fct = NULL){
  check_names(input = data, req = c("timestamp", "depth", fct), type = all)
  if(!is.null(fct)) data$fct <- data[, fct] else data$fct <- 1L
  id_n_obs <-
    data %>%
    dplyr::group_by(.data$fct) %>%
    dplyr::summarise(n = dplyr::n(), .groups = "drop_last")
  if(any(id_n_obs$n <= 2)){
    id_to_drop <- id_n_obs$fct[which(id_n_obs$n <= 2)]
    warning(paste0("Too few observations for the following individual(s) for analysis: ", paste0("'", id_to_drop, collapse = "', "), "'."))
    data <- data %>% dplyr::filter(.data$fct != id_to_drop)
    if(nrow(data) <= 2) stop("Insufficient data data remaining for analysis.")
  }
  data_1 <- data %>% dplyr::filter(fct == unique(data$fct)[1])
  step   <- as.numeric(difftime(data_1$timestamp[2], data_1$timestamp[1]), units = "secs")
  data <-
    data %>%
    dplyr::group_by(.data$fct) %>%
    dplyr::mutate(dist = abs(Tools4ETS::serial_difference(.data$depth)),
                  speed_ms = .data$dist/step,
                  speed_mstep = .data$dist) %>%
    dplyr::ungroup() %>%
    dplyr::filter(!is.na(.data$dist)) %>%
    dplyr::select(-.data$fct)
  message("Mobility estimates from n = ", nrow(data), " observation(s).")
  cat("--------------------------------------\n")
  cat("Estimates (m/s)-----------------------\n")
  stats <- data.frame(variable = "speed", min = min(data$speed_ms), mean = mean(data$speed_ms), max = max(data$speed_ms))
  stats[, 2:ncol(stats)] <- round(stats[, 2:ncol(stats)], 2)
  print(stats)
  cat("Estimates (m/step)--------------------\n")
  stats <- data.frame(variable = "speed", min = min(data$speed_mstep), mean = mean(data$speed_mstep), max = max(data$speed_mstep))
  stats[, 2:ncol(stats)] <- round(stats[, 2:ncol(stats)], 2)
  print(stats)
  cat("--------------------------------------\n")
  return(invisible(data))
}


#########################################
#########################################
#### get_mvt_resting()

#' @title Identify resting behaviour within depth time series
#' @importFrom rlang .data
#' @description This function implements a simple threshold-based approach to identify periods of resting behaviour within depth time series collected from benthic animals. In these data, resting is usually inferred from low vertical activity and sinusoidal tidal traces that reflect the rise-and-fall of the tides on an individual's depth when it is stationary on the seabed. To implement the approach, a dataframe comprising depth time series (\code{archival}) is required. The algorithm identifies any time steps in which the change in depth is less than a specified value (\code{th_depth}) as potentially exhibiting resting behaviour. Candidate states (resting = 0, not resting = 1) assigned from this depth threshold can be smoothed by supplying a function (\code{weights}) that is applied to states over a time window (\code{th_time}). For example (\code{weights = mean}) calculates the mean state for each observation from all observations within \code{th_time}. Average state scores can be discretised into `resting' (0) or `not resting' (1) states by specifying the maximum score of an observation for it to be defined as `resting' via \code{discrete}. For example, if \code{weights = mean} and \code{discrete = 0}, then an observation is only defined as `resting' if all of the observations in the surrounding window (controlled by \code{th_time}) are also equal to zero. Higher scores dampen the effects of occasional movements that exceed the depth threshold. The function returns a numeric vector of average or discrete states, with one element for each observation in \code{archival}.
#'
#' @param archival A dataframe that defines depth time series. This should contain a numeric column of depth observations named `depth'.
#' @param fct (optional) A character that defines the name of a column in \code{archival} that distinguishes independent time series.
#' @param th_depth A double that defines the depth threshold. Candidate resting behaviour is identified when the absolute change in depth is \eqn{\leq} \code{th_depth}.
#' @param th_time (optional) An integer that defines the number of time steps (i.e., rows in \code{archival}) around each depth observation that should influence state assignment, via the \code{weights} function and related arguments (\code{align}, \code{...} and \code{discrete}).
#' @param weights (optional) If \code{th_time} is provided, \code{weights} is a function applied to the candidate binary states assigned from the depth threshold (\code{th_depth}) within each time interval (\code{th_time}). The default is \code{\link[base]{mean}}, which applies a rolling mean to the states within each window.
#' @param align,... (optional) Additional arguments, passed to \code{\link[data.table]{frollapply}}, that control the implementation of \code{weights} to each window (\code{th_time}). By default, \code{align = "centre"}, which applies \code{weights} to the \code{th_time/2} observation(s) centred around each observation. (The states for observations that are too close to the start or end of the time series for averaging are retained as initially assigned by \code{th_depth} but not averaged.)
#' @param discrete (optional) If \code{th_time} is provided, \code{discrete} is double between 0 and 1 that defines the maximum score of an observation for it to be defined as `resting'. With the default \code{weights = mean} function, this is equivalent to the maximum proportion of `not active' observations in the window around each depth observation (\code{th_time}) for an observation to be defined as `resting'. For instance, \code{discrete = 0} means that a `resting' is only defined if all of the observations around in the window around each observation are less than or equal to the depth threshold (\code{th_depth}), while \code{discrete = 0.05} means `resting' is defined if \eqn{\leq} 5 percent of observations in the time window are below the depth threshold.
#'
#' @details This approach was motivated by the need to identify efficiently `resting' periods in  depth time series data, sampled at a resolution of two-minutes, collected from flapper skate (\emph{Dipturus intermedius}). While there are other approaches for the identification of resting behaviour, such as Hidden Markov models, for the flapper skate depth time series these are much more computationally demanding and appear to be much less effective at correctly assigning states.
#'
#' @return The function returns a numeric vector, of the same length as \code{archival$depth}, that defines, for each depth observation, either (a) a discrete behavioural state, if \code{th_time = NULL} or \code{discrete} is supplied (i.e., resting = 0, not resting = 1) or (b) an averaged behavioural score (0 - 1) over \code{th_time}.
#'
#' @examples
#' #### Example (1): Assign 'resting' based on a simple depth threshold
#' dat_archival$state_1 <- get_mvt_resting(archival = dat_archival,
#'                                         fct = "individual_id",
#'                                         th_depth = 0.25)
#'
#' #### Example (2): Assign 'resting' based on depth and time thresholds
#' # ... Under the default settings, all of observations in the time threshold
#' # ... must be below the depth threshold to qualify as 'resting'
#' dat_archival$state_2 <- get_mvt_resting(archival = dat_archival,
#'                                         fct = "individual_id",
#'                                         th_depth = 0.25,
#'                                         th_time = 30)
#'
#' #### Example (3): Dampen the effects of occasionally exceeding the depth threshold
#' # ... by increasing the proportion of observations that are allowed to
#' # ... exceed the depth threshold in a each time window
#' dat_archival$state_3 <- get_mvt_resting(archival = dat_archival,
#'                                         fct = "individual_id",
#'                                         th_depth = 0.25,
#'                                         th_time = 30,
#'                                         discrete = 0.05)
#'
#' #### Example (4): Return average state scores via discrete = NULL
#' dat_archival$state_4 <- get_mvt_resting(archival = dat_archival,
#'                                         fct = "individual_id",
#'                                         th_depth = 0.25,
#'                                         th_time = 30,
#'                                         discrete = NULL)
#'
#' #### Compare the frequency distribution of states among methods
#' # In the first example, a large number of 'resting' states are assigned
#' # In the second example, there are far fewer 'resting' states because all of the
#' # ... observations in the specified time window around each depth observation need
#' # ... to meet the depth threshold for the observation to be defined as 'resting',
#' # ... not just the observation itself.
#' # In the third example, there are sightly more 'resting' states because the
#' # ... criterion for all observations to meet the depth threshold has been weakened.
#' # In the final example, average scores have been returned.
#' pp <- graphics::par(mfrow = c(2, 2))
#' prettyGraphics::pretty_hist(dat_archival$state_1, main = "method_1")
#' prettyGraphics::pretty_hist(dat_archival$state_2, main = "method_2")
#' prettyGraphics::pretty_hist(dat_archival$state_3, main = "method_3")
#' prettyGraphics::pretty_hist(dat_archival$state_4, main = "method_4")
#' graphics::par(pp)
#'
#' #### Compare the time series for an example individual for each method
#' ## Filter results for specific individual
#' dat_archival_1 <- dat_archival[dat_archival$individual_id == 25, ]
#' ## Define helper functions
#' # Define helper function to plot blank depth time series
#' plot_blank <- function(...){
#'   prettyGraphics::pretty_plot(dat_archival_1$timestamp, dat_archival_1$depth*-1,
#'                               pretty_axis_args = list(side = 3:2),
#'                               type = "n",...)
#' }
#' # Define helper function to add depth time series, coloured by state, to the plot
#' add_lines_for_state <- function(state,...){
#'   prettyGraphics::add_lines(x = dat_archival_1$timestamp,
#'                             y1 = dat_archival_1$depth*-1,
#'                             y2 = dat_archival_1[, state],...)
#' }
#' ## Make plots
#' pp <- graphics::par(mfrow = c(2, 2))
#' # state 1
#' plot_blank(main = "method_1")
#' add_lines_for_state("state_1")
#' # state 2
#' plot_blank(main = "method_2")
#' add_lines_for_state("state_2")
#' # state 3
#' plot_blank(main = "method_3")
#' add_lines_for_state("state_3")
#' # state 4
#' plot_blank(main = "method_4")
#' add_lines_for_state("state_4")
#' graphics::par(pp)
#'
#' @author James Thorburn conceptualised and developed this method. Edward Lavender extended the initial approach and wrote the implementation for the \code{\link[flapper]{flapper}} package.
#' @export

get_mvt_resting <- function(archival,
                            fct = NULL,
                            th_depth = 0.25,
                            th_time = NULL,
                            weights = mean, align = "center",...,
                            discrete = 0
                            ){

  #### Implement function checks
  check_names(input = archival, req = c("depth", fct))
  if(!is.null(fct)){
    check_class(input = fct, to_class = "character", coerce_input = as.character)
  }
  if(!is.null(discrete)){
    if(!(discrete >= 0 & discrete <= 1)) stop("'discrete' should be between 0 and 1 inclusive.")
  }

  #### Define states based on depth only
  if(is.null(fct)) archival$fct <- 1L else archival$fct <- archival[, fct]
  archival <-
    archival %>%
    dplyr::group_by(.data$fct) %>%
    dplyr::mutate(va_abs = abs(Tools4ETS::serial_difference(.data$depth))) %>%
    dplyr::mutate(state_1 = dplyr::if_else(.data$va_abs <= th_depth, 0, 1)) %>%
    dplyr::mutate(state_2 = .data$state_1)

  #### Average states over time interval
  if(!is.null(th_time)){
    # Average states over time window via rollapply()
    archival <-
      archival %>%
      dplyr::mutate(state_2 = data.table::frollapply(.data$state_1,
                                                     n = th_time,
                                                     FUN = weights, align = align,...)) %>%
      dplyr::mutate(state_2 = dplyr::if_else(is.na(.data$state_2), .data$state_1, .data$state_2))

    # Assign discrete states based on threshold proportion parameter
    if(!is.null(discrete)){
      archival <- archival %>% dplyr::mutate(state_2 = dplyr::if_else(.data$state_2 <= discrete, 0, 1))
    }
  }

  #### Return outputs
  return(archival$state_2)
}
