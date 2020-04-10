#' @title Compute euclidean distances between receivers
#' @description This function computes linear distances (km) between all combinations of receivers. 
#' 
#' @param moorings A dataframe which defines the location of each unique receiver combination. This should contain the columns: 'ID', a unique identifier of each receiver, 'lat', the latitude of that receiver in decimal degrees; and 'long', the longitude of that receiver in decimal degrees. 
#' @param f (optional) A function which is used to process distances before these are returned. For example, it may be useful to round distances to nearest km with \code{f = function(x) round(x, digits = 0)}. 
#' 
#' @return The function returns a dataframe with columns 'r1', 'r2' and 'dist'. These define the IDs of each combination of receivers and the associated distance between them, in km. Note that the dataframe contains duplicate combinations of receivers (e.g., both r1 = 1 and r2 = 2 and r1 = 2 and r2 = 1). 
#' 
#' @examples 
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
    
    #### Define all combinations of receivers 
    dists <- expand.grid(r1 = moorings$ID, r2 = moorings$ID)
    
    #### Define lat and long 
    dists$lat1 <- moorings$lat[match(dists$r1, moorings$ID)]
    dists$long1 <- moorings$long[match(dists$r1, moorings$ID)]
    dists$lat2 <- moorings$lat[match(dists$r2, moorings$ID)]
    dists$long2 <- moorings$long[match(dists$r2, moorings$ID)]
    
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
