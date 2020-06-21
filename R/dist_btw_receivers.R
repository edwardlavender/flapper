#' @title Compute euclidean distances between receivers
#' @description This function computes linear distances (km) between all combinations of receivers.
#'
#' @param moorings A dataframe which defines the location of each unique receiver combination. This should contain the columns: 'receiver_id', a unique identifier of each receiver, 'receiver_lat', the latitude of that receiver in decimal degrees; and 'receiver_long', the longitude of that receiver in decimal degrees (see \code{\link[flapper]{dat_moorings}}.
#' @param f (optional) A function which is used to process distances before these are returned. For example, it may be useful to round distances to nearest km with \code{f = function(x) round(x, digits = 0)}.
#'
#' @return The function returns a dataframe with columns 'r1', 'r2' and 'dist'. These define the IDs of each combination of receivers and the associated distance between them, in km. Note that the dataframe contains duplicate combinations of receivers (e.g., both r1 = 1 and r2 = 2 and r1 = 2 and r2 = 1).
#'
#' @examples
#' #### Example (1): Compute distances between all combinations of receivers in km
#' # Define dataframe with required columns
#' dat <- data.frame(receiver_id = dat_moorings$receiver_id,
#'                   receiver_lat = dat_moorings$receiver_lat,
#'                   receiver_long = dat_moorings$receiver_long)
#' # Compute distances
#' dist_btw_receivers_km <- dist_btw_receivers(dat)
#' head(dist_btw_receivers_km)
#'
#' #### Example (2): Post-process distances via the f argument
# round distances
#' dist_btw_receivers_km_round <- dist_btw_receivers(dat, f = round)
#' head(dist_btw_receivers_km_round)
#' # convert distances to m
#' dist_btw_receivers_m <- dist_btw_receivers(dat, f = function(x) x*1000)
#' head(dist_btw_receivers_m)
#'
#' @author Edward Lavender
#' @export
#'

########################################
########################################
#### dist_btw_receivers

dist_btw_receivers <-
  function(moorings,
           f = NULL){

    #### Checks
    stopifnot(all(c("receiver_id", "receiver_long", "receiver_lat") %in% colnames(moorings)))

    #### Define all combinations of receivers
    dists <- expand.grid(r1 = moorings$receiver_id, r2 = moorings$receiver_id)

    #### Define lat and long
    dists$lat1 <- moorings$receiver_lat[match(dists$r1, moorings$receiver_id)]
    dists$long1 <- moorings$receiver_long[match(dists$r1, moorings$receiver_id)]
    dists$lat2 <- moorings$receiver_lat[match(dists$r2, moorings$receiver_id)]
    dists$long2 <- moorings$receiver_long[match(dists$r2, moorings$receiver_id)]

    #### Compute the distances between each receiver combination in km
    dists$dist <- NA
    for(i in 1:nrow(dists)){
      # calculate distances
      dists$dist[i] <- geosphere::distGeo(c(dists$long1[i], dists$lat1[i]), c(dists$long2[i], dists$lat2[i]))
      # convert output from m to km
      dists$dist[i] <- dists$dist[i]/1000
    }

    #### Process dataframe
    if(!is.null(f)){
      dists$dist <- f(dists$dist)
    }
    dists <- dists[, c("r1", "r2", "dist")]

    #### Return dataframe
    return(dists)

}

#### End of code.
########################################
########################################
