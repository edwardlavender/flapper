#' @title Pass putative false detections through a spatial filter
#' @description The identification of false detections in acoustic telemetry data is an important aspect of processing and/or modelling these data. False detections can be identified using the short interval criterion, whereby any detection of an individual at a receiver which is not accompanied by other detections at the same receiver in a specified time window (depending on the nominal acoustic transmission delay) are flagged. This approach can be implemented using the \code{\link[glatos]{false_detections}} function in the \code{\link[glatos]{glatos}} package. Flag detections can then be examined, or passed through other filters, to examine their plausibility. This function passes false detections flagged by \code{\link[glatos]{false_detections}} through a spatial filter. The key idea is that detections at nearby receivers within the define time window can be used to suggest 'false' detections which are, in fact, plausible. To implement this approach, the user must define a dataframe comprising detections, a temporal threshold, a spatial threshold and a dataframe of distances between receivers. The function examines whether any putative false detection are accompanied by additional detections at other receivers within a user-defined euclidean distance of that receiver. If so, these could be explained by an individual which dips in-and-out of the detection ranges of receivers (e.g. in a sparse acoustic array) and may not, in fact, be false.
#'
#' @param det A dataframe containing detection data. This should contain the following columns: 'detection_timestamp_utc', 'transmitter_codespace', 'transmitter_id' and 'receiver_sn', as well as 'passed_filter' (see \code{\link[glatos]{false_detections}}).
#' @param tf A number which defines the time threshold (s) used to flag false detections (see \code{\link[glatos]{false_detections}}).
#' @param sf A number which defines the threshold Euclidean distance between receivers beyond which, even if a false detection is accompanied by detections at other receivers, it is likely to be a true false detection, because the the individual could not have moved between receivers separated by more than this threshold over the specified time interval. \code{sf} should be defined in the same units as the distances provided in \code{dist} (see below).
#' @param dist_btw_receivers A dataframe which defines the distances between receiver pairs. This should contain the columns: 'r1', 'r2' and 'dist', whereby 'r1' and 'r2' contain the unique receiver serial number for all combinations of receivers and 'dist' contains the distance between them. This dataframe should include duplicate combinations (e.g., both r1 = 1 and r2 = 2 and r1 = 2 and r2 = 1). This can be created with \code{\link[flapper]{dist_btw_receivers}}.
#'
#' @details There are limitations with the application of this spatial filter to false detections. First, the spatial threshold beyond which false detections are likely to be false is based on Euclidean distances at present. These may be problematic (e.g. when receivers hug complex coastlines). Second, for small arrays, fast-swimming organisms and/or a large nominal transmission delay (i.e., time threshold), the spatial filter is a poor filter because individuals can access the whole area over the whole time interval.
#'
#' @return The function returns a vector, of the same length as \code{det} with three possible values: NA, which identifies detections which have not been flagged as false detections (i.e., \code{passed_filter = 0}, see \code{\link[glatos]{false_detections}}) and are therefore not passed through the spatial filter; 1, which identifies detections which 'passed' the spatial filter (i.e., false detections which are accompanied by detections at nearby receivers within the defined spatial and temporal thresholds); or 0, which defines detections which 'failed' the spatial filter (i.e., false detections which are not accompanied by detections at nearby receivers within the defined spatial and temporal thresholds). This vector has one attribute, 'details', a dataframe with the same number of rows as \code{det} with the following columns: 'passed_filter_sf', 'n_wn_sf', 'detection_timestamp_utc_sf' and 'receiver_sn_sf'. 'passed_filter_sf' takes a value of NA, 0 or 1 depending on whether or not the detection was flagged as a false detection (if not, NA) and whether each false detection passed the spatial filter (no, 0; yes, 1). If the detection did pass the spatial filter, 'n_wn_sf' provides the number of detections at nearby receivers within \code{tf} and \code{sf}; and 'detection_timestamp_utc_sf', 'receiver_sn_sf' and 'dist_sf' define the timestamp of the detection at the nearest receiver, the receiver at which the detection was made and the distance between the two receivers respectively.
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
#' det$passed_filter_sf <- false_detections_sf(det,
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
#'


#######################################
#######################################
#### false_detections_sf (false detections spatial filter)

false_detections_sf <-
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


#### End of code.
#######################################
#######################################
