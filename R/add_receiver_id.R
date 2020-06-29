#' @title Add unique receiver IDs to passive acoustic telemetry timeseries
#' @description This function extracts the unique receiver deployment IDs (e.g., 1,...n) from a dataframe containing receiver attributes (e.g. unique receiver deployment IDs, receiver codes and locations), termed \code{moorings}, that correspond to passive acoustic telemetry (PAT) receiver codes in detection timeseries, termed \code{acoustics}. This is especially useful if the same receiver has been deployed multiple times (e.g. in different locations) over the course of the study. In this scenario, the receiver code in the PAT detection timeseries does not uniquely identify the unique receiver deployment, which means that receiver codes in PAT data and metadata cannot simply be matched: instead, both receiver code(s) and the time of detection(s) need to be included in the matching procedure. In this case, \code{add_receiver_ids()} examines each receiver with detections in \code{acoustics}. If that receiver was redeployed, the function examines the time of each detection and returns the ID of the receiver that recorded that detection (given the timing of receiver deployment). If each receiver was only deployed once, this function simply matches receiver codes in \code{moorings} and receiver codes in \code{acoustics} and returns the receiver IDs (which, in this case, may be the same as the receiver codes) in \code{moorings} that correspond to each receiver code in \code{acoustics}.
#'
#' @param acoustics A dataframe which comprises passive acoustic telemetry detection timeseries. This must contain two named columns: 'timestamp', a POSIXct vector which defines the time of each detection; and 'receiver', a vector which defines the receiver at which the individual was detected. The column 'receiver' should also be found in \code{moorings} (see below).
#' @param moorings A dataframe which contains passive acoustic telemetry receiver metadata. This must contain four named columns: 'receiver', as above for \code{acoustics}; 'receiver_id', a unique identifier for each receiver deployment; 'receiver_start_date', a POSIXct vector which defines the time of each receiver deployment; and 'receiver_end_date', a POSIXct vector which defines the end of each receiver deployment. If objects of class Date are provided for 'receiver_start_date' and 'receiver_end_date', these are coerced to POSIXct objects with a warning so that the timing of receiver detections and deployment is comparable (i.e., of the same object type).
#'
#' @details The function implements \code{\link[data.table]{foverlaps}} to overlap the timing of detections with the timing of receiver deployments to account for receiver deployment in the assignment of receiver IDs.
#'
#' @return The function returns a vector of receiver IDs, as defined in the \code{moorings$receiver_id} column, which correspond to each detection in the \code{acoustics} dataframe.
#'
#' @examples
#' #### Define example data
#' # In this example, we have two receivers, but one has been re-deployed:
#' moorings <- data.frame(receiver = c(1, 1, 2),
#'                         receiver_id = c(1, 2, 3),
#'                         start_date = as.POSIXct(c("2016-01-01",
#'                                                   "2016-01-02",
#'                                                   "2016-01-01"), tz = "UTC"),
#'                         end_date = as.POSIXct(c("2016-01-02",
#'                                                 "2016-01-03",
#'                                                 "2016-01-02"), tz = "UTC")
#'                         )
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
#' acoustics$receiver_id <- add_receiver_id(acoustics, moorings)
#' acoustics
#'
#' @seealso see \code{\link[Tools4ETS]{add_unit_id}} for a slightly more general function.
#'
#' @author Edward Lavender
#' @export
#'


#############################################
#############################################
#### add_receiver_id()

add_receiver_id <-
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



#### End of code.
#############################################
#############################################
