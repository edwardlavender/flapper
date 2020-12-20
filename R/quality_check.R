#' @title Basic quality checks of passive acoustic telemetry datasets
#' @description This function passes through passive acoustic telemetry datasets through some basic quality checks (see Details). Following data processing, these provide a useful 'final check' prior to analysis.
#' @param acoustics A dataframe which comprises passive acoustic telemetry detection time-series. This must contain the following named columns: 'timestamp', a POSIXct vector which defines the time of each detection; 'receiver_id', a unique identifier of each receiver; and 'individual_id', a unique identifier of each individual (see \code{\link[flapper]{dat_acoustics}} for an example).
#' @param moorings A dataframe which contains passive acoustic telemetry receiver metadata. This must contain the following named columns: 'receiver_id', a unique identifier for each receiver deployment; 'receiver_start_date', a POSIXct vector which defines the time of each receiver deployment; and 'receiver_end_date', a POSIXct vector which defines the end of each receiver deployment (see \code{\link[flapper]{dat_moorings}} for an example). If objects of class Date are provided for 'receiver_start_date' and 'receiver_end_date', these are coerced to POSIXct objects with a warning.
#' @param ids A dataframe which contains the passive acoustic telemetry individual metadata. This must contain the following named columns: 'individual_id', a unique identifier for each individual; 'tag_start_date', a POSIXct vector which defines the time of each tag deployment; and 'tag_end_date', a POSIXct vector which defines the end of each tag deployment (see \code{\link[flapper]{dat_ids}} for an example). If objects of class Date are provided for 'tag_start_date' and 'tag_end_date', these are coerced to POSIXct objects with a warning.
#' @details The function implements the following checks. (1) \code{acoustics} should only contain receivers recorded in \code{moorings}; other receivers may be included in centralised databases (e.g., from other projects) and often need to be removed. (2) Observations at receivers should occur during their deployment windows; other observations may be included in centralised databases due to receiver checks, range testing or re-deployment elsewhere. (3) \code{acoustics} is checked for any unknown tag IDs. These may be due unrecorded use of tags for receiver checking or range testing, other tagging programmes or type A false detections. (4) As for detections at receivers, all detections of tags should occur during their deployment windows. (5) False detections should be flagged. Other important checks - such as checking for receivers which were lost and later recovered, excluding observations during receiver servicing dates, excluding observations during tag recapture events and further investigation of false detections - may be required. \code{\link[flapper]{quality_check}} is mainly designed to be implemented after data-processing has already taken place as a basic 'final check' for common issues in passive acoustic telemetry datasets. For each check, the function returns a message or warning depending on the outcome; subsequently, the most appropriate course of action (e.g., retention versus removal of flagged observations in \code{acoustics} will depend on the context).
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
#' #### Implement quality_check() on processed data as a final check for any issues
#' quality_check(dat_acoustics, dat_moorings, dat_ids)
#'
#' #### Add erroneous data to acoustics for demonstrating quality_check()
#' ## Define a convenience function to add erroneous data to
#' # ... acoustics to demonstrate quality_check()
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
#' ## Add erroneous timestamps (outside receiver/individual id deployment periods )
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
#' #### Implement quality_check()
#' quality_check(acoustics_wth_errors, dat_moorings, dat_ids)
#'
#' @author Edward Lavender
#' @export
#'

#########################################
#########################################
#### quality_check()


quality_check <-
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


#### End of code.
#########################################
#########################################
