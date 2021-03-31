######################################
######################################
#### process_receiver_id()

#' @title Add unique receiver IDs to passive acoustic telemetry time series
#' @description This function extracts the unique receiver deployment IDs (e.g., 1,...n) from a dataframe containing receiver attributes (e.g. unique receiver deployment IDs, receiver codes and locations), termed \code{moorings}, that correspond to passive acoustic telemetry (PAT) receiver codes in detection time series, termed \code{acoustics}. This is especially useful if the same receiver has been deployed multiple times (e.g. in different locations) over the course of the study. In this scenario, the receiver code in the PAT detection time series does not uniquely identify the unique receiver deployment, which means that receiver codes in PAT data and metadata cannot simply be matched: instead, both receiver code(s) and the time of detection(s) need to be included in the matching procedure. In this case, \code{\link[flapper]{process_receiver_id}} examines each receiver with detections in \code{acoustics}. If that receiver was redeployed, the function examines the time of each detection and returns the ID of the receiver that recorded that detection (given the timing of receiver deployment). If each receiver was only deployed once, this function simply matches receiver codes in \code{moorings} and receiver codes in \code{acoustics} and returns the receiver IDs (which, in this case, may be the same as the receiver codes) in \code{moorings} that correspond to each receiver code in \code{acoustics}.
#'
#' @param acoustics A dataframe which comprises passive acoustic telemetry detection time series. This must contain two named columns: `timestamp', a POSIXct vector which defines the time of each detection; and `receiver', a vector which defines the receiver at which each detection was made. The column `receiver' should also be found in \code{moorings} (see below).
#' @param moorings A dataframe which contains passive acoustic telemetry receiver metadata. This must contain four named columns: `receiver', as above for \code{acoustics}; `receiver_id', a unique identifier for each receiver deployment; `receiver_start_date', a POSIXct vector which defines the time of each receiver deployment; and `receiver_end_date', a POSIXct vector which defines the end of each receiver deployment. If objects of class Date are provided for `receiver_start_date' and `receiver_end_date', these are coerced to POSIXct objects with a warning so that the timing of receiver detections and deployment is comparable (i.e., of the same object type).
#'
#' @details The function implements \code{\link[data.table]{foverlaps}} to overlap the timing of detections with the timing of receiver deployments to account for receiver deployment in the assignment of receiver IDs.
#'
#' @return The function returns a vector of receiver IDs, as defined in the \code{moorings$receiver_id} column, which correspond to each detection in the \code{acoustics} dataframe.
#'
#' @examples
#' #### Define example data
#' # In this example, we have two receivers, but one has been re-deployed:
#' moorings <- data.frame(receiver = c(1, 1, 2),
#'                        receiver_id = c(1, 2, 3),
#'                        start_date = as.POSIXct(c("2016-01-01",
#'                                                  "2016-01-02",
#'                                                  "2016-01-01"), tz = "UTC"),
#'                        end_date = as.POSIXct(c("2016-01-02",
#'                                                "2016-01-03",
#'                                                "2016-01-02"), tz = "UTC")
#'                        )
#' # Our observational dataframe contains receivers but not unique receiver IDs:
#' acoustics <- data.frame(receiver = c(1, 1, 2),
#'                        timestamp = as.POSIXct(c("2016-01-01 00:30:00",
#'                                                  "2016-01-02 00:30:00",
#'                                                  "2016-01-01 00:30:00"), tz = "UTC")
#'                                               )
#'
#' #### Example (1): Add unique receiver IDs to the observational dataframe
#' # The first observation corresponds to receiver 1;
#' # The second observation corresponds to the same receiver
#' # ... but a different deployment, and has receiver_id = 2
#' # The third observation corresponds to receiver id 3;
#' acoustics$receiver_id <- process_receiver_id(acoustics, moorings)
#' acoustics
#'
#' @seealso \code{\link[Tools4ETS]{add_unit_id}} is a slightly more general function.
#'
#' @author Edward Lavender
#' @export

process_receiver_id <-
  function(acoustics, moorings){

    #### Checks
    # Check that moorings and acoustics contain required columns
    if(!all(c("receiver", "receiver_id", "start_date", "end_date") %in% colnames(moorings))) stop("moorings does not contain all required column names.")
    if(!all(c("receiver", "timestamp") %in% colnames(acoustics))) stop("acoustics does not contain all required column names.")

    #### If some receivers were deployed, we need to
    # ... account for the receiver code and the time of deployment when we add receiver IDs to the acoustics df.
    if(any(duplicated(moorings$receiver))){

      #### Check acoustics/moorings dataframes contain required information in correct format
      if(!inherits(acoustics$timestamp, "POSIXct")) stop("acoustics$timestamp must an POSIXct object")
      if(lubridate::tz(acoustics$timestamp) == ""){
        warning("acoustics$timestamp lacking timezone (see lubridate::tz()). tz = 'UTC' forced.")
        lubridate::tz(acoustics$timestamp) <- "UTC"
      }
      tz <- lubridate::tz(acoustics$timestamp)
      if(!inherits(moorings$start_date, "POSIXct")){
        warning("moorings$start_date must be a POSIXct object; attempting to coerce moorings$start_date into a POSIXct object.")
        moorings$start_date <- lubridate::round_date(as.POSIXct(moorings$start_date, tz = tz), unit = "day")
        lubridate::tz(moorings$start_date) <- tz
      }
      if(!inherits(moorings$end_date, "POSIXct")){
        warning("moorings$end_date must be a POSIXct object; attempting to coerce moorings$end_date into a POSIXct object.")
        moorings$end_date <- lubridate::round_date(as.POSIXct(moorings$end_date, tz = tz), unit = "day")
        lubridate::tz(moorings$end_date) <- tz
      }
      if(lubridate::tz(moorings$start_date) == ""){
        warning("moorings$start_date lacking timezone (see lubridate::tz()). tz = lubridate::tz(acoustics$timestamp) forced.")
        lubridate::tz(moorings$start_date) <- tz
      }
      if(lubridate::tz(moorings$end_date) == ""){
        warning("moorings$timestamp lacking timezone (see lubridate::tz()). tz = lubridate::tz(acoustics$timestamp) forced.")
        lubridate::tz(moorings$end_date) <- tz
      }
      if(class(moorings$receiver)[1] != class(acoustics$receiver)[1]){
        warning("class(moorings$receiver)[1] != class(acoustics$receiver)[1]; both coerced to character vectors.")
        moorings$receiver <- as.character(moorings$receiver)
        acoustics$receiver <- as.character(acoustics$receiver)
      }

      #### Define data.tables
      acoustics$start_date <- acoustics$timestamp
      acoustics$end_date <- acoustics$timestamp
      acoustics <- data.table::data.table(acoustics)
      moorings <- data.table::data.table(moorings)
      receiver <- NULL; start_date <- NULL; end_date <- NULL
      data.table::setkey(moorings, receiver, start_date, end_date)

      #### Implement data.table::foverlaps()
      order <- data.table::foverlaps(acoustics, moorings, type = "within", nomatch = NA, which = TRUE)

      #### Define acoustics$receiver_id
      acoustics$receiver_id <- moorings$receiver_id[order$yid]

    #### If there are no redeployed receivers, we can simply use match()
    # ... to obtain the receiver IDs for the acoustics dataframe.
    } else{
      acoustics$receiver_id <- moorings$receiver_id[match(acoustics$receiver, moorings$receiver)]
    }

    #### Return the receiver IDs
    lna <- length(which(is.na(acoustics$receiver_id)))
    if(lna > 0) warning(paste("receiver IDs returned contain,", lna, "NAs (e.g. possibly due to incorrect start/end dates)."))
    return(acoustics$receiver_id)
  }


#########################################
#########################################
#### process_false_detections_sf()

#' @title Pass putative false detections through a spatial filter
#' @description The identification of false detections in acoustic telemetry data is an important aspect of processing and/or modelling these data. False detections can be identified using the short interval criterion, whereby any detection of an individual at a receiver which is not accompanied by other detections at the same receiver in a specified time window (depending on the nominal acoustic transmission delay) are flagged. This approach can be implemented using the \code{\link[glatos]{false_detections}} function in the \code{\link[glatos]{glatos}} package. Flag detections can then be examined, or passed through other filters, to examine their plausibility.
#'
#' This function passes false detections flagged by \code{\link[glatos]{false_detections}} through a spatial filter. The key idea is that detections at nearby receivers within a defined time window may, in fact, be plausible. To implement this approach, the user must define a dataframe comprising detections, a temporal threshold, a spatial threshold and a dataframe of distances between receivers. The function examines whether any putative false detections are accompanied by additional detections at other receivers within a user-defined time window and Euclidean distance of that receiver. If so, these could be explained by an individual that dips in-and-out of the detection ranges of receivers (e.g. in a sparse acoustic array) and may not, in fact, be false.
#'
#' @param det A dataframe containing detection time series. Following \code{\link[glatos]{false_detections}}, this should contain the following columns: `detection_timestamp_utc', `transmitter_codespace', `transmitter_id' and `receiver_sn', as well as `passed_filter' (see \code{\link[glatos]{false_detections}}).
#' @param tf A number that defines the time threshold (s) used to flag false detections (see \code{\link[glatos]{false_detections}}).
#' @param sf A number that defines the threshold Euclidean distance between receivers beyond which, even if a false detection is accompanied by detections at other receivers, it is likely to be a true false detection, because the individual could not have moved between receivers separated by more than this threshold over the specified time interval. \code{sf} should be defined in the same units as the distances provided in \code{dist_btw_receivers} (see below).
#' @param dist_btw_receivers A dataframe that defines the distances between receiver pairs. This should contain the columns: `r1', `r2' and `dist', whereby `r1' and `r2' contain the unique receiver serial number for all combinations of receivers and `dist' contains the distance between them. This dataframe should include duplicate combinations (e.g., both r1 = 1 and r2 = 2 and r1 = 2 and r2 = 1). This can be created with \code{\link[flapper]{dist_btw_receivers}}.
#'
#' @details There are limitations with the application of this spatial filter to false detections. First, the spatial threshold beyond which false detections are likely to be false is based on Euclidean distances at present. These may be problematic (e.g. when receivers hug complex coastlines). Second, for small arrays, fast-swimming organisms and/or a large nominal transmission delay (i.e., time threshold), the spatial filter is a poor filter because individuals can access the whole area over the whole time interval.
#'
#' @return The function returns a vector, of the same length as \code{det} with three possible values: NA, which identifies detections which have not been flagged as false detections (i.e., \code{passed_filter = 0}, see \code{\link[glatos]{false_detections}}) and are therefore not passed through the spatial filter; 1, which identifies detections which `passed' the spatial filter (i.e., false detections which are accompanied by detections at nearby receivers within the defined spatial and temporal thresholds); or 0, which defines detections which `failed' the spatial filter (i.e., false detections which are not accompanied by detections at nearby receivers within the defined spatial and temporal thresholds). This vector has one attribute, `details', a dataframe with the same number of rows as \code{det} with the following columns: `passed_filter_sf', `n_wn_sf', `detection_timestamp_utc_sf' and `receiver_sn_sf'. `passed_filter_sf' takes a value of NA, 0 or 1 depending on whether or not the detection was flagged as a false detection (if not, NA) and whether each false detection passed the spatial filter (no, 0; yes, 1). If the detection did pass the spatial filter, `n_wn_sf' provides the number of detections at nearby receivers within \code{tf} and \code{sf}; and `detection_timestamp_utc_sf', `receiver_sn_sf' and `dist_sf' define the timestamp of the detection at the nearest receiver, the receiver at which the detection was made and the distance between the two receivers respectively.
#'
#' @examples
#' #### Add some false detections for demonstration purposes
#' # Add three rows to dat_acoustics which, below, we'll make 'false detections'
#' dat_acoustics_with_false_det <-
#'   rbind(dat_acoustics, dat_acoustics[rep(nrow(dat_acoustics), 3), ])
#' pos_false <- (nrow(dat_acoustics_with_false_det) - 2):nrow(dat_acoustics_with_false_det)
#' # Add an isolated detection accompanied by a detection at a nearby receiver
#' dat_acoustics_with_false_det$timestamp[pos_false[1:2]] <-
#'   dat_acoustics_with_false_det$timestamp[pos_false[1:2]] + 60*60*60
#' dat_acoustics_with_false_det$receiver_id[pos_false[2]] <- 33
#' # Add an isolated detection not accompanied by a detection at a nearby receiver
#' dat_acoustics_with_false_det$timestamp[pos_false[3]] <-
#'   dat_acoustics_with_false_det$timestamp[pos_false[3]] + 60*60*60*2
#'
#' #### Define necessary columns to compute false detections using glatos::false_detections()
#' dat_acoustics_with_false_det$detection_timestamp_utc <-
#'   dat_acoustics_with_false_det$timestamp
#' dat_acoustics_with_false_det$transmitter_codespace <-
#'   substr(dat_acoustics_with_false_det$transmitter_id, 1, 8)
#' dat_acoustics_with_false_det$transmitter_id <-
#'   substr(dat_acoustics_with_false_det$transmitter_id, 10, 13)
#' dat_acoustics_with_false_det$receiver_sn <- dat_acoustics_with_false_det$receiver_id
#' det <- dat_acoustics_with_false_det[, c("detection_timestamp_utc",
#'                                         "transmitter_codespace",
#'                                         "transmitter_id",
#'                                         "receiver_sn")]
#'
#' #### Compute false detections
#' # 3 false detections returned, as expected:
#' det <- glatos::false_detections(det, tf = 3600)
#' tail(det$passed_filter)
#'
#' #### Pass false detections through a spatial filter
#' # distances between receivers are required
#' dist_btw_receivers_km <-
#'   dist_btw_receivers(dat_moorings[, c("receiver_id", "receiver_long", "receiver_lat")])
#' # Implement spatial filter.
#' # Note the function returns a vector, unlike glatos::false_detections():
#' det$passed_filter_sf <- process_false_detections_sf(det,
#'                                             tf = 3600,
#'                                             sf = 0.5,
#'                                             dist_btw_receivers = dist_btw_receivers_km)
#' # Only the last observation failed the spatial filter, as expected:
#' tail(det$passed_filter_sf)
#' # Additional information is available from the attributes dataframe:
#' tail(attr(det$passed_filter_sf, "details"))
#'
#'
#' @author Edward Lavender
#' @export

process_false_detections_sf <-
  function(det,
           tf,
           sf,
           dist_btw_receivers
  ){

    #### Initial checks
    # Check that there are at least some detections which failed glatos' false detection filter
    if(!any(det$passed_filter == 0)){
      warning("No detections failed glatos' false detection filter; det returned unchanged.")
      return(det)
    }
    # Sort the dataframe by transmitter then timestamp
    # This is necessary for later.
    rownames(det) <- 1:nrow(det)
    det$original_order <- 1:nrow(det)
    pos_fail_original_order <- which(det$passed_filter == 0)
    det <- det[order(det$transmitter_id, det$detection_timestamp_utc), ]

    #### Define the lower and upper time around each detection, defining a window
    # ... within which we'll search for other detections at other receivers.
    det$t1 <- det$detection_timestamp_utc - tf
    det$t2 <- det$detection_timestamp_utc + tf

    #### Extract the data which failed the checks.
    pos_fail <- which(det$passed_filter == 0)
    min_lag_fail <- det[pos_fail, ]

    #### Define blank dataframe
    dout <- data.frame(passed_filter = NA, n_wn_sf = NA, detection_timestamp_utc = as.POSIXct(NA), receiver_sn = NA, dist_btw_receivers = NA)
    dout_ls <- lapply(pos_fail, function(i) return(dout))

    #### Loop over all data pertaining to failed detections
    dout_ls <-
      pbapply::pbmapply(
        min_lag_fail$t1,
        min_lag_fail$t2,
        min_lag_fail$transmitter_id,
        min_lag_fail$receiver_sn,
        dout_ls,
        FUN = function(t1, t2, trans_id, rec_id, dout){

          #### Testing
          # testing <- FALSE
          # if(testing){
          #  p <- 4
          #  pos <- pos_fail[p]
          #  t1 = min_lag_fail$t1[p]
          #  t2 = min_lag_fail$t2[p]
          #  trans_id = min_lag_fail$transmitter_id[p]
          #  rec_id = min_lag_fail$receiver_sn[p]
          # }

          #### Define a subsetted detection dataframe which contains
          # ... all detections of a transmitted within the necessary time interval
          # ... i.e., include possible detections at other receivers.
          det_tmp <- det[det$transmitter_id == trans_id &
                           det$detection_timestamp_utc >= t1 &
                           det$detection_timestamp_utc <= t2, ]

          #### If det_tmp only contains one detection (i.e., at the receiver in question)
          # ... then no detections were made at other receivers, so we'll return a 0
          # ... i.e., it still looks like this is a false detection.
          if(nrow(det_tmp) == 1){
            dout$passed_filter <- 0

            #### If there are detections at multiple receivers,
            # ... then we'll extract the distances between the receiver with the worrying detection and
            # ... all of those other receivers. If the distance is below a user-defined threshold,
            # ... perhaps this isn't a false detection.
          } else{

            #### Remove the detection receiver in question
            # Otherwise this will cause issues with the calculation of distances between receivers
            # ... and the identification of the smallest distances.
            det_tmp <- det_tmp[det_tmp$receiver_sn != rec_id, ]

            #### Extract distances between rec_id and all other receiver in det_tmp with a detection
            # (Note calculations for all receivers, not just unique receivers, which makes the
            # ... code below simpler e.g. reporting the number of detections within the threshold distance.)
            dist_sbt <- dist_btw_receivers[which(dist_btw_receivers$r1 == rec_id & dist_btw_receivers$r2 %in% det_tmp$receiver_sn), ]
            det_tmp$dist <- dist_sbt$dist[match(det_tmp$receiver_sn, dist_sbt$r2)]
            # |dist$r1 %in% det_tmp$receiver_sn & dist$r2 == rec_id), ]

            #### If any of the distances are less than a threshold, this suggests they
            # ... may not be false after all. Hence, return a 1. Otherwise, return a 0
            # ... i.e., it still looks like a false detection.
            if(any(det_tmp$dist <= sf)){
              # Define passed filter and the number of detections at other receivers within the necessary distance
              dout$passed_filter <- 1
              dout$n_wn_sf <- length(which(det_tmp$dist <= sf))
              # Define the receiver, distance and detection timestamp of the (first) nearest detection in space
              # (There may be others at slightly further receivers sooner or later in time, but, because the dataframe
              # ... has been sorted, only the first is identified).
              pos_min_dist <- which.min(det_tmp$dist)
              dout$dist <- min(det_tmp$dist)
              dout$receiver_sn <- det_tmp$receiver_sn[pos_min_dist]
              dout$detection_timestamp_utc <- det_tmp$detection_timestamp_utc[pos_min_dist]

            } else{
              dout$passed_filter <- 0
            }

          }

          # Return the updated dataframe
          return(dout)

        }, SIMPLIFY = FALSE)

    #### Join the list into a single dataframe
    dout <- dplyr::bind_rows(dout_ls)

    #### Message
    n_false_detections <- nrow(det) - sum(det$passed_filter)
    n_false_detections_passed_spatial <- sum(dout$passed_filter)
    n_false_detections_failed_spatial <- nrow(dout) - sum(dout$passed_filter)

    message(paste0("The spatial filter retained ",
                   n_false_detections_failed_spatial,
                   " detections, out of ", n_false_detections, " previously identified false detections ",
                   " (", round(n_false_detections_failed_spatial/n_false_detections * 100, 2),
                   " %) as 'true' false detections."
    )
    )

    #### Return a vector of whether or not each
    # ... observation passed the spatial filter,
    # ...with other details remaining as an attribute.
    det$passed_filter_sf <- NA
    det$n_wn_sf <- NA
    det$detection_timestamp_utc_sf <- as.POSIXct(NA)
    det$receiver_sn_sf <- NA
    det$dist_sf <- NA
    det$passed_filter_sf[pos_fail] <- dout$passed_filter
    det$n_wn_sf[pos_fail] <- dout$n_wn_sf
    det$detection_timestamp_utc_sf[pos_fail] <- dout$detection_timestamp_utc
    det$receiver_sn_sf[pos_fail] <- dout$receiver_sn
    det$dist_sf[pos_fail] <- dout$dist
    det <- det[order(det$original_order), ]
    attr(det$passed_filter_sf, "details") <- det[, c("passed_filter_sf",
                                                     "n_wn_sf",
                                                     "receiver_sn_sf",
                                                     "detection_timestamp_utc_sf",
                                                     "dist_sf")]
    return(det$passed_filter_sf)

  }


#########################################
#########################################
#### process_quality_check()

#' @title Basic quality checks of passive acoustic telemetry datasets
#' @description This function passes through passive acoustic telemetry datasets through some basic quality checks (see Details). Following data processing, these provide a useful `final check' prior to analysis.
#'
#' @param acoustics A dataframe which comprises passive acoustic telemetry detection time series. This must contain the following named columns: `timestamp', a POSIXct vector which defines the time of each detection; `receiver_id', a unique identifier of each receiver; and `individual_id', a unique identifier of each individual (see \code{\link[flapper]{dat_acoustics}} for an example).
#' @param moorings A dataframe which contains passive acoustic telemetry receiver metadata. This must contain the following named columns: `receiver_id', a unique identifier for each receiver deployment; `receiver_start_date', a POSIXct vector which defines the time of each receiver's deployment; and `receiver_end_date', a POSIXct vector which defines the end of each receiver's deployment (see \code{\link[flapper]{dat_moorings}} for an example). If objects of class Date are provided for `receiver_start_date' and `receiver_end_date', these are coerced to POSIXct objects with a warning.
#' @param ids A dataframe which contains the passive acoustic telemetry individual metadata. This must contain the following named columns: `individual_id', a unique identifier for each individual; `tag_start_date', a POSIXct vector which defines the time of each tag's deployment; and `tag_end_date', a POSIXct vector which defines the end of each tag's deployment (see \code{\link[flapper]{dat_ids}} for an example). If objects of class Date are provided for `tag_start_date' and `tag_end_date', these are coerced to POSIXct objects with a warning.
#'
#' @details The function implements the following checks:
#' \enumerate{
#' \item Valid receivers. \code{acoustics} should only contain receivers recorded in \code{moorings}; other receivers may be included in centralised databases (e.g., from other projects) and often need to be removed.
#' \item Valid detections (at receivers). Observations at receivers should occur during their deployment windows; other observations may be included in centralised databases due to receiver checks, range testing or re-deployment elsewhere.
#' \item Valid tags. \code{acoustics} is checked for any unknown tag IDs. These may be due unrecorded use of tags for receiver checking or range testing, other tagging programmes or type A false detections.
#' \item Valid detections (of tags). As for detections at receivers, all detections of tags should occur during their deployment windows.
#' \item False detections. False detections should be flagged.
#' }
#'
#' \code{\link[flapper]{process_quality_check}} is mainly designed to be implemented after data-processing has already taken place as a basic `final check' for common issues in passive acoustic telemetry datasets. For each check, the function returns a message or warning depending on the outcome; subsequently, the most appropriate course of action (e.g., retention versus removal of flagged observations in \code{acoustics} will depend on the context). Other important checks -- such as checking for receivers which were lost and later recovered, excluding observations during receiver servicing dates, excluding observations during tag capture events and further investigation of false detections -- may be required.
#'
#' @return For each check, the function returns a message or a warning with relevant details.
#'
#' @examples
#' #### Prepare data
#' ## All data have previously passed false detection filters (see glatos::false_detections())
#' dat_acoustics$passed_filter <- 1
#' ## Times should be in POSIXct format
#' dat_moorings$receiver_start_date <-  as.POSIXct(dat_moorings$receiver_start_date)
#' lubridate::tz(dat_moorings$receiver_start_date) <- "UTC"
#' dat_moorings$receiver_end_date   <-  as.POSIXct(dat_moorings$receiver_end_date)
#' lubridate::tz(dat_moorings$receiver_end_date) <- "UTC"
#' dat_ids$tag_start_date <- as.POSIXct(dat_ids$tag_start_date)
#' lubridate::tz(dat_ids$tag_start_date) <- "UTC"
#' ## tag_end_date column needed in dat_ids
#' dat_ids$tag_end_date   <- as.POSIXct("2020-01-01", tz = "UTC")
#'
#' #### Implement process_quality_check() on processed data as a final check for any issues
#' process_quality_check(dat_acoustics, dat_moorings, dat_ids)
#'
#' #### Add erroneous data to acoustics for demonstrating process_quality_check()
#' ## Define a convenience function to add erroneous data to
#' # ... acoustics to demonstrate process_quality_check()
#' add_erroneous_row <- function(acoustics, row = nrow(acoustics), col, val){
#'   tmp_ls <- lapply(val, function(v){
#'     tmp <- acoustics[row, ]
#'     tmp[1, col] <- v
#'     return(tmp)
#'   })
#'   tmp <- dplyr::bind_rows(tmp_ls)
#'   acoustics <- rbind(acoustics, tmp)
#'   return(acoustics)
#' }
#' ## Add erroneous receiver ids
#' nrw <- nrow(dat_acoustics)
#' acoustics_wth_errors <- add_erroneous_row(dat_acoustics,
#'                                           row = nrw,
#'                                           col = "receiver_id",
#'                                           val = c(100, 200, 300))
#' ## Add erroneous time stamps (outside receiver/individual id deployment periods )
#' acoustics_wth_errors <- add_erroneous_row(acoustics_wth_errors,
#'                                          row = nrw,
#'                                           col = "timestamp",
#'                                           val = as.POSIXct(c("2019-01-01", "2019-03-01"),
#'                                                            tz = "UTC"))
#' ## Add erroneous individual ids
#' acoustics_wth_errors <- add_erroneous_row(acoustics_wth_errors,
#'                                           row = nrw,
#'                                           col = "individual_id",
#'                                           val = c(100, 200, 300))
#' ## Examine erroneous data:
#' utils::tail(acoustics_wth_errors, 10)
#'
#' #### Implement process_quality_check()
#' process_quality_check(acoustics_wth_errors, dat_moorings, dat_ids)
#'
#' @author Edward Lavender
#' @export

process_quality_check <-
  function(acoustics,
           moorings,
           ids){

    #### Checks
    ## acoustics colnames
    check_names(arg = "acoustics", input = acoustics,
                req = c("timestamp", "receiver_id", "individual_id"),
                extract_names = colnames, type = all)
    ## moorings colnames
    check_names(arg = "moorings",
                input = moorings,
                req = c("receiver_id", "receiver_start_date", "receiver_end_date"),
                extract_names = names, type = all)
    ## ids colnames
    check_names(arg = "ids",
                input = ids,
                req = c("individual_id", "tag_start_date", "tag_end_date"),
                extract_names = names, type = all)
    ## format of times and timezones
    times_ls <-
      mapply(list("acoustics$timestamp",
                  "moorings$receiver_start_date", "moorings$receiver_end_date",
                  "ids$tag_start_date", "ids$tag_end_date"),
             list(acoustics$timestamp,
                  moorings$receiver_start_date, moorings$receiver_end_date,
                  ids$tag_start_date, ids$tag_end_date),
             FUN = function(arg, elm){
               out <- check_class(arg = arg,
                                  input = elm,
                                  if_class = NULL,
                                  to_class = "POSIXct",
                                  type = "warning",
                                  coerce_input = function(x) as.POSIXct(x, tz = "UTC"))
               out <- check_tz(arg, elm)
               return(out)
             }, SIMPLIFY = FALSE)
    acoustics$timestamp          <- times_ls[[1]]
    moorings$receiver_start_date <- times_ls[[2]]
    moorings$receiver_end_date   <- times_ls[[3]]
    ids$tag_start_date           <- times_ls[[4]]
    ids$tag_end_date             <- times_ls[[5]]

    #### Receiver identity
    # All receivers should have been deployed in the study in question.
    runknown <- unique(acoustics$receiver_id)[!(unique(acoustics$receiver_id) %in% moorings$receiver_id)]
    lrunknown <- length(runknown)
    if(lrunknown > 0){
      pos <- which(acoustics$receiver_id %in% runknown)
      lpos <- length(pos)
      warn <- paste0("Check 1 (receiver identity): failed. ",
                     lrunknown, " receiver identities unknown (", paste(runknown, collapse = ", "), "), ",
                     "corresponding to ", lpos, " observations in acoustics. \n")
      warning(warn, immediate. = TRUE)
      acoustics <- acoustics[-pos, ]
    } else{
      message("Check 1 (receiver identity): passed. \n")
    }

    #### Receiver loss
    # Not currently implemented.

    #### Receiver operation window
    match_receiver <- match(acoustics$receiver_id, moorings$receiver_id)
    acoustics$receiver_start_date  <- moorings$receiver_start_date[match_receiver]
    acoustics$receiver_end_date    <- moorings$receiver_end_date[match_receiver]
    acoustics$interval   <- lubridate::interval(acoustics$receiver_start_date, acoustics$receiver_end_date)
    acoustics$not_within <- !(acoustics$timestamp %within% acoustics$interval)
    if(length(which(acoustics$not_within)) > 0){
      pos  <- which(acoustics$not_within)
      lpos <- length(pos)
      warn <- paste0("Check 3 (receiver deployment windows): failed. ",
                     lpos, " observations in acoustics outside of receiver deployment windows. \n")
      warning(warn, immediate. = TRUE)
      acoustics <- acoustics[-pos, ]
    } else{
      message("Check 2 (receiver deployment windows): passed. \n")
    }

    #### Receiver servicing dates.
    # Not currently implemented.

    #### Tag identity
    tunknown <- unique(acoustics$individual_id)[!(unique(acoustics$individual_id) %in% ids$individual_id)]
    ltunknown <- length(tunknown)
    if(ltunknown > 0){
      pos <- which(acoustics$individual_id %in% tunknown)
      lpos <- length(pos)
      warn <- paste0("Check 3 (tag identity): failed. ",
                     ltunknown, " tag identities unknown (", paste(tunknown, collapse = ", "), "), ",
                     "corresponding to ", lpos, " observations in acoustics. \n")
      warning(warn, immediate. = TRUE)
      acoustics <- acoustics[-pos, ]
    } else{
      message("Check 3 (tag identity): passed. \n")
    }

    #### Tag loss
    # Not currently implemented.

    #### Tag operation window
    match_tag <- match(acoustics$individual_id, ids$individual_id)
    acoustics$tag_start_date <- ids$tag_start_date[match_tag]
    acoustics$tag_end_date   <- ids$tag_end_date[match_tag]
    acoustics$interval   <- lubridate::interval(acoustics$tag_start_date, acoustics$tag_end_date)
    acoustics$not_within <- !(acoustics$timestamp %within% acoustics$interval)
    if(length(which(acoustics$not_within)) > 0){
      pos  <- which(acoustics$not_within)
      lpos <- length(pos)
      warn <- paste0("Check 3 (tag deployment windows): failed. ",
                     lpos, " observations in acoustics outside of tag deployment windows. \n")
      warning(warn, immediate. = TRUE)
      acoustics <- acoustics[-pos, ]
    } else{
      message("Check 4 (tag deployment windows): passed. \n")
    }

    #### Tag recapture
    # Not currently implemented.

    #### False detections
    if(!rlang::has_name(acoustics, "passed_filter")){
      warn <- "Check 5 (false detections): 'passed_filter' column not found in acoustics. See glatos::false_detections() to analyse false detections."
      warning(warn, immediate. = TRUE)
    } else{
      lpos <- length(which(acoustics$passed_filter == 0))
      if(lpos > 1){
        warn <- paste0("There are ", lpos, " false detections in the acoustics$passed_filter column.")
        warning(warn, immediate. = TRUE)
      } else{
        message("Check 5 (false detections) passed: acoustics$passed_filter does not contain false detections. \n")
      }
    }

  }


#########################################
#########################################
#### process_behav_rest()

#' @title Identify resting behaviour within depth time series
#' @importFrom rlang .data
#' @description This function implements a simple threshold-based approach to identify periods of resting behaviour within  depth time series collected from benthic animals. In these data, resting is usually inferred from low vertical activity and sinusoidal tidal traces that reflect the rise-and-fall of the tides on an individual's depth when it is stationary on the seabed. To implement the approach, a dataframe comprising depth time series (\code{archival}) is required. The algorithm identifies any time steps in which the change in depth is less than a specified value (\code{th_depth}) as potentially exhibiting resting behaviour. Candidate states (resting = 0, not resting = 1) assigned from this depth threshold can be smoothed by supplying a function (\code{weights}) that is applied to states over a time window (\code{th_time}). For example (\code{weights = mean}) calculates the mean state for each observation from all observations within \code{th_time}. Average state scores can be discretised into `resting' (0) or `not resting' (1) states by specifying the maximum score of an observation for it to be defined as `resting' via \code{discrete}. For example, if \code{weights = mean} and \code{discrete = 0}, then an observation is only defined as `resting' if all of the observations in the surrounding window (controlled by \code{th_time}) are also equal to zero. Higher scores dampen the effects of occasional movements that exceed the depth threshold. The function returns a numeric vector of average or discrete states, with one element for each observation in \code{archival}.
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
#' dat_archival$state_1 <- process_behav_rest(archival = dat_archival,
#'                                            fct = "individual_id",
#'                                            th_depth = 0.25)
#'
#' #### Example (2): Assign 'resting' based on depth and time thresholds
#' # ... Under the default settings, all of observations in the time threshold
#' # ... must be below the depth threshold to qualify as 'resting'
#' dat_archival$state_2 <- process_behav_rest(archival = dat_archival,
#'                                            fct = "individual_id",
#'                                            th_depth = 0.25,
#'                                            th_time = 30)
#'
#' #### Example (3): Dampen the effects of occasionally exceeding the depth threshold
#' # ... by increasing the proportion of observations that are allowed to
#' # ... exceed the depth threshold in a each time window
#' dat_archival$state_3 <- process_behav_rest(archival = dat_archival,
#'                                            fct = "individual_id",
#'                                            th_depth = 0.25,
#'                                            th_time = 30,
#'                                            discrete = 0.05)
#'
#' #### Example (4): Return average state scores via discrete = NULL
#' dat_archival$state_4 <- process_behav_rest(archival = dat_archival,
#'                                            fct = "individual_id",
#'                                            th_depth = 0.25,
#'                                            th_time = 30,
#'                                            discrete = NULL)
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

process_behav_rest <- function(archival,
                               fct = NULL,
                               th_depth = 0.25,
                               th_time = NULL,
                               weights = mean, align = "center",...,
                               discrete = 0
){
  #### Implement function checks
  check_names(input = archival, req = c("depth", fct))
  check_class(input = fct, to_class = "character", coerce_input = as.character)
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


#########################################
#########################################
#### process_surface()

#' @title Process a Raster* by aggregation and quantify the error induced by this process
#' @description This function reduces the resolution of a \code{\link[raster]{raster}} by multiple aggregation methods and then quantifies the relative error induced by each method from the differences between the original values and the aggregated values. To implement the function, a \code{\link[raster]{raster}} (\code{x}) must be supplied as well as an aggregation factor (\code{fact}) and a named list of functions (\code{stat}) used to aggregate the \code{\link[raster]{raster}}. The \code{\link[raster]{raster}} is aggregated using each method (function) and mapped back onto the original resolution for calculation of the differences between the original \code{\link[raster]{raster}} and the aggregated \code{\link[raster]{raster}}(s). The function returns a visual statistical summary of the differences (if \code{plot = TRUE}) and a named list comprising the aggregated \code{\link[raster]{raster}}(s) and the re-sampled version(s) of those mapped back onto the original resolution.
#'
#' @param x A \code{\link[raster]{raster}} to be processed. For implementations preceding a call to of one of \code{\link[flapper]{flapper}}'s particle filtering algorithms, \code{x} should be planar (i.e., Universal Transverse Mercator projection) with equal resolution in the x, y directions and identical units in the x, y and z directions (e.g., see \code{\link[flapper]{dcpf}}).
#' @param fact A positive integer that defines by how much \code{x} should be aggregated (see \code{\link[raster]{aggregate}}).
#' @param stat A named list of functions used to aggregate \code{x} (see the \code{fun} argument of \code{\link[raster]{aggregate}}).
#' @param ... Additional arguments passed to \code{\link[raster]{aggregate}} to control aggregation.
#' @param plot A logical input that defines whether or not to plot a summary of the differences between the original \code{\link[raster]{raster}} (\code{x}) and the aggregated \code{\link[raster]{raster}}(s). If specified, the minimum, median and maximum difference are shown for each statistic (\code{stat}).
#' @param cl (optional) A cluster object created by \code{\link[parallel]{makeCluster}}. If supplied, the connection to the cluster is stopped within the function.
#' @param varlist (optional) A character vector of names of objects to export, to be passed to the \code{varlist} argument of \code{\link[parallel]{clusterExport}}. This may be required if \code{cl} is supplied. Exported objects must be located in the global environment.
#' @param verbose A logical input that defines whether or not to print messages to the console to relay function progress.
#'
#' @details This function was motivated by the particle filtering algorithms in \code{\link[flapper]{flapper}} (e.g., \code{\link[flapper]{dcpf}}). For these algorithms, it is computationally beneficial to reduce \code{\link[raster]{raster}} resolution, where possible, by aggregation. To facilitate this process, this function quantifies the relative error induced by different aggregation functions. If appropriate, the particle filtering algorithm(s) can then be implemented using the aggregated \code{\link[raster]{raster}} that minimises the error, with the magnitude of that error propagated via the \code{depth_error} parameter.
#'
#' @return The function returns a plot of the differences between the original and aggregated \code{\link[raster]{raster}}(s), if \code{plot = TRUE}, and a named list of (a) the aggregated \code{\link[raster]{raster}}(s) (`agg_by_stat'), (b) the aggregated, resampled \code{\link[raster]{raster}}(s) (`agg_by_stat_rs') and (c) the summary statistics plotted.
#'
#' @examples
#' # Define the raster for which to implement the function
#' x <- dat_gebco
#' blank <- raster::raster(raster::extent(x), crs = raster::crs(x), resolution = 250)
#' x <- raster::resample(x, blank, method = "bilinear")
#' # Implement function using a list of statistics
#' out <- process_surface(x, fact = 2, stat = list(min = min, mean = mean, median = median, max = max))
#' summary(out)
#'
#' @seealso \code{\link[raster]{aggregate}}, \code{\link[raster]{resample}}
#' @author Edward Lavender
#' @export

process_surface <- function(x,
                            fact = 2L,
                            stat = list(mean = mean),...,
                            plot = TRUE,
                            cl = NULL, varlist = NULL,
                            verbose = TRUE){

  # Set up function
  t_onset <- Sys.time()
  cat_to_console <- function(..., show = verbose) if(show) cat(paste(..., "\n"))
  cat_to_console(paste0("flapper::process_surface() called (@ ", t_onset, ")..."))
  check_named_list(input = stat)

  # Define blank raster with same extent
  x_blank <- x

  # Aggregate raster by each statistic
  cat_to_console("... Aggregating raster...")
  if(!is.null(cl) & !is.null(varlist)) parallel::clusterExport(cl = cl, varlist = varlist)
  x_agg_by_stat <- pbapply::pblapply(stat, cl = cl, function(foo){
    x_agg <- raster::aggregate(x, fact = fact, fun = foo,...)
    return(x_agg)
  })

  # Re-sample aggregated rasters to original resolution
  cat_to_console("... Resampling aggregated raster(s) back onto the original resolution...")
  x_agg_by_stat_rs <- pbapply::pblapply(x_agg_by_stat, cl = cl, function(x_agg){
    x_agg_rs <- raster::resample(x_agg, x_blank, method = "ngb")
    return(x_agg_rs)
  })
  if(!is.null(cl)) parallel::stopCluster(cl)

  # Get differences between original raster and aggregated (resampled) rasters for each statistic
  cat_to_console("... Computing differences between the original and aggregated raster(s)...")
  x_agg_by_stat_rs_diff <- pbapply::pblapply(x_agg_by_stat_rs, function(x_agg_rs){
    x_agg_rs_diff <- x - x_agg_rs
    return(x_agg_rs_diff)
  })

  # Summarise differences
  if(plot){
    cat_to_console("... Summarising the differences between rasters across statistic(s)...")
    mins <- sapply(x_agg_by_stat_rs_diff, raster::cellStats, stat = "min")
    meds <- sapply(x_agg_by_stat_rs_diff, raster::cellStats, stat = "mean")
    maxs <- sapply(x_agg_by_stat_rs_diff, raster::cellStats, stat = "max")
    xp <- factor(names(stat), levels = names(stat))
    prettyGraphics::pretty_plot(xp, meds,
                                ylim = range(c(mins, maxs)),
                                type = "n", xlab= "Statistic", ylab = "Difference [x - x_agg]")
    prettyGraphics::add_error_bars(x = xp, fit = meds, lwr = mins, upr = maxs)
    x_summary_stats <- data.frame(stat = names(stat),
                                min = mins,
                                median = meds,
                                max = maxs)
  } else x_summary_stats <- NULL


  # Return outputs
  out <- list(agg_by_stat = x_agg_by_stat,
              agg_by_stat_rs_diff = x_agg_by_stat_rs_diff,
              summary_stats = x_summary_stats)
  t_end <- Sys.time()
  duration <- difftime(t_end, t_onset, units = "mins")
  cat_to_console(paste0("... flapper::process_surface() call completed (@ ", t_end, ") after ~", round(duration, digits = 2), " minutes."))
  return(out)

}
