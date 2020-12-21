######################################
######################################
#### setup_acdc()

#' @title Setup the ACDC algorithm
#' @description This function produces the acoustic contours required by the acoustic-centroid depth-contour (ACDC) algorithm.
#' @param rs A integer vector of receiver IDs.
#' @param xy A \code{\link[sp]{SpatialPoints}} object that defines the locations of each receiver. The order of points in this object should match the order of receivers defined in \code{rs}. The coordinate reference system should be the Universal Transverse Mercator system with distances in metres (to match \code{detection_range}, see below).
#' @param detection_range A number that defines the maximum detection range at which an individual could be detected from a receiver.
#' @param mobility A number that defines the distance that an individual could move in the time period between archival observations.
#' @param n_timesteps An integer that defines the the number of timesteps after a hypothetical detection for which centroids will be created, where the duration of each timestep is given by the duration between archival observations.
#' @param coastline (optional) A \code{\link[sp]{SpatialPolygonsDataFrame-class}} object that defines the coastline in an area. If provided, acoustic centroids are processed to remove any areas on land. Algorithm speed declines with the complexity of the coastline.
#' @param boundaries (optional) A \code{\link[raster]{extent}} object that defines the boundaries of an area within which individuals are assumed to have remained. If provided, acoustic centroids are processed to remain within this area.
#' @param plot A logical input that defines whether or not to produce a plot of the area, including receivers, the coastline and the area boundaries (if provided), and acoustic centroids. This is useful for checking purposes but it can reduce algorithm speed.
#' @param cl (optional) A cluster object created by \code{\link[parallel]{makeCluster}} to implement the algorithm in parallel. The connection to the cluster is closed within the function.
#' @param verbose A logical input that defines whether or not to print messages to the console to relay function progress.
#'
#' @return The function returns a list of \code{\link[sp]{SpatialPolygonsDataFrame-class}} objects, with one element for all numbers from 1 to the maximum receiver number (\code{rx}). Any list elements that do not correspond to receivers contain a \code{NULL} element. List elements that correspond to receivers contain a \code{\link[sp]{SpatialPolygonsDataFrame-class}} object containing all the centroids for that receiver.
#'
#' @examples
#' #### Define data for setup_acdc()
#' ## Define coordinates of receivers as SpatialPoints with UTM CRS
#' # CRS of receiver locations as recorded in dat_moorings
#' proj_wgs84 <- sp::CRS("+init=epsg:4326")
#' # CRS of receiver locations required
#' proj_utm <- sp::CRS(paste("+proj=utm +zone=29 +datum=WGS84",
#'                           "+units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))
#' # Define SpatialPoints object
#' xy_wgs84 <- sp::SpatialPoints(dat_moorings[, c("receiver_long", "receiver_lat")], proj_wgs84)
#' xy_utm <- sp::spTransform(xy_wgs84, proj_utm)
#'
#' #### Example (1): Define a list of centroids with specified parameters
#' # ... (Argument values are small to reduce computation time for examples)
#' centroids <- setup_acdc(rs = dat_moorings$receiver_id,
#'                         xy = xy_utm,
#'                         detection_range = 500,
#'                         mobility = 250,
#'                         n_timesteps = 3
#'                         )
#' # A list of SpatialPolygonsDataFrames is returned with elements from 1:max(rs)
#' # NULL elements correspond to numbers in this sequence that do not refer to receivers
#' # Otherwise a SpatialPolygonsDataFrame is returned with all the centroids for that receiver
#' centroids
#'
#' #### Example (2): Visualise the acoustic centroids produced via plot = TRUE
#' centroids <- setup_acdc(rs = dat_moorings$receiver_id,
#'                         xy = xy_utm,
#'                         detection_range = 500,
#'                         mobility = 250,
#'                         n_timesteps = 3,
#'                         plot = TRUE
#'                         )
#'
#' #### Example (3): Remove areas of the centroids that overlap with coastline
#' centroids <- setup_acdc(rs = dat_moorings$receiver_id,
#'                         xy = xy_utm,
#'                         detection_range = 500,
#'                         mobility = 250,
#'                         n_timesteps = 3,
#'                         plot = TRUE,
#'                         coastline = dat_coast
#'                         )
#'
#' #### Example (4): Remove areas of the centroids beyond a boundary
#' xy_utm_coords <- sp::coordinates(xy_utm)
#' boundaries <- raster::extent(min(xy_utm_coords[, 1]),
#'                              max(xy_utm_coords[, 1]),
#'                              min(xy_utm_coords[, 2]),
#'                              max(xy_utm_coords[, 2])
#'                         )
#' centroids <- setup_acdc(rs = dat_moorings$receiver_id,
#'                         xy = xy_utm,
#'                         detection_range = 500,
#'                         mobility = 250,
#'                         n_timesteps = 3,
#'                         plot = TRUE,
#'                         coastline = dat_coast,
#'                         boundaries = boundaries
#'                         )
#'
#' #### Example (5): Implement the algorithm in parallel
#' centroids <- setup_acdc(rs = dat_moorings$receiver_id,
#'                         xy = xy_utm,
#'                         detection_range = 500,
#'                         mobility = 250,
#'                         n_timesteps = 3,
#'                         plot = TRUE,
#'                         coastline = dat_coast,
#'                         boundaries = boundaries,
#'                         cl = parallel::makeCluster(2L)
#'                         )
#'
#' #### Example (6): Acoustic centroids can be saved to file using rlist::list.save()
#' # rlist::list.save(centroids, paste0(tempdir(), "/centroids.RData"))

#' @author Edward Lavender
#' @export
#'

setup_acdc <- function(
  rs,
  xy,
  detection_range,
  mobility,
  n_timesteps = 250,
  coastline = NULL,
  boundaries = NULL,
  plot = FALSE,
  cl = NULL,
  verbose = TRUE
){

  #### Initiate function
  cat_to_console <- function(..., show = verbose) if(show) cat(paste(..., "\n"))
  cat_to_console("flapper::setup_acdc() called...")
  if(is.numeric(rs)) rs <- as.integer(rs)
  if(!is.integer(rs)) stop(paste("Argument 'rs' must be of class 'integer', not class(es):"), class(rs))
  if(any(rs <= 0)) stop("Argument 'rs' cannot contain receiver IDs <= 0.")
  if(any(duplicated(rs))){
    message("Argument 'rs' contains duplicate elements. rs has been simplified to unique(rs).")
    rs <- unique(rs)
  }
  if(!is.null(coastline) & !is.null(boundaries)) {
    coastline <- raster::crop(coastline, boundaries)
    if(is.null(coastline)) message("No coastline within defined boundaries. \n")
  }
  if(plot){
    cat_to_console("... Plotting background map of area...")
    if(!is.null(coastline)) {
      raster::plot(coastline, col = "lightgreen", border = "darkgreen", lwd = 1.5)
      graphics::points(xy, pch = 21, col = "royalblue", bg = "royalblue")
    } else {
      raster::plot(xy, pch = 21, col = "royalblue", bg = "royalblue")
    }
    if(!is.null(boundaries)) raster::lines(boundaries, col = "red", lty = 3)
  }


  #### Define a list of acoustic centroids for each receiver
  cat_to_console("... Building a nested list of acoustic centroids. This is the slow step...")

  ## Define a sequence of centroid sizes
  # Around each receiver, we'll create a polygon of this size
  size_seq <- seq(detection_range, length.out = n_timesteps, by = mobility)

  ## Define a list of receivers, with a list of centroids for each receiver
  bathy_ls <- pbapply::pblapply(1:length(rs), cl = cl, function(i){

      centroids_ls <- lapply(size_seq, function(size){

        # Define a buffer around the current receiver of appropriate size
        bathy_poly <- rgeos::gBuffer(xy[i], width = size)

        # Reduce the size of the polygon by overlapping to remove areas on land
        # This keeps the polygons as small as possible which is important for ACDC/MP algorithm computation efficiency.
        if(!is.null(coastline)) {
          bathy_poly <- rgeos::gDifference(bathy_poly, coastline, byid = FALSE)
        }

        # Remove any areas beyond specified boundaries
        # Again this keeps polygon size to a minimum
        if(!is.null(boundaries)) {
          bathy_poly <- raster::crop(bathy_poly, boundaries)
        }

        # Return acoustic centroid
        return(bathy_poly)

      })

    # Define names of the rasters forn receiver, i, based on size
    names(centroids_ls) <- paste0("s_", size_seq)
    return(centroids_ls)

    })
  if(!is.null(cl)) parallel::stopCluster(cl)
  names(bathy_ls) <- rs

  #### Add NULL elements to the list for any receivers in the range 1:max(rs) that are not in rs
  # This means we can use receiver numbers to go straight to the correct element in the list in ACDC/MP algorithms.
  bathy_ls <- lapply(as.integer(1:max(rs)), function(i){
    if(i %in% rs){
      return(bathy_ls[[as.character(i)]])
    } else{
      return(NULL)
    }
  })

  #### Convert nested list of polygons to a SpatialPolygonsDataFrame
  # ... with one element for each receiver and each dataframe containing all the polygons for that receiver
  cat_to_console("... Converting the nested list of acoustic centroids to a SpatialPolygonsDataFrame...")
  if(is.null(bathy_ls)) {
    stop("There are no acoustic centroids within defined spatial boundaries.")
  }
  spdf_ls <- pbapply::pblapply(bathy_ls, cl = NULL, function(element){
    if(!is.null(element)){
      # bind all the sub-elements in each element together into a single spatial polygon
      sp <- raster::bind(element)
      # convert this into a spatial polygons dataframe
      spdf <- sp::SpatialPolygonsDataFrame(sp, data.frame(id = 1:length(sp)))
      return(spdf)
    }})

  #### Visualise centroids
  if(plot){
    cat_to_console("... Plotting centroids on map...")
    pbapply::pblapply(spdf_ls, function(spdf) if(!is.null(spdf)) raster::lines(spdf, col = "dimgrey", lwd = 0.75))
  }

  #### Return list of SpatialPolygonsDataFrame
  return(spdf_ls)

}



######################################
######################################
#### .acdc()




######################################
######################################
#### acdc()




######################################
######################################
#### animate_acdc()


