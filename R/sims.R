######################################
######################################
#### sim_array()

#' @title Simulate (marine) monitoring arrays
#' @description This function is designed to simulate different kinds of array designs for monitoring stations. The function has been particularly inspired by the need to simulate passive acoustic telemetry array designs, which comprise networks of acoustic hydrophones that listen for acoustic transmissions from tagged marine animals. To implement the function, it is necessary to define the boundaries of the area (\code{boundaries}), within which barriers to movement, such as coastline, can be simulated or included from real datasets. Within this area, a specified number of receivers (\code{n_receivers}) can be simulated under different array designs (namely, uniform, regular or random arrangements) or incorporated from real data. The function returns a list of spatial objects that define the array and, if requested, a plot of the area.
#'
#' @param boundaries A \code{\link[raster]{extent}} object that defines the boundaries of simulated study site.
#' @param coastline (optional) This argument is used to incorporate the presence of barriers, such as coastline, within an array. There are three options. If \code{coastline = NULL}, no barriers are incorporated. If \code{coastline = "simple_random"}, then some coastline is simulated in the study area. Alternatively, a spatial object, such as a SpatialPolygonsDataFrame that defines the coastline in an area can be incorporated into the array design by passing this \code{coastline}.
#' @param land_inside_coastline A logical variable that defines whether or not the land is 'inside' the polygon defined by \code{coastline} (\code{land_inside_coastline = TRUE}) or the sea is 'inside' (\code{land_inside_coastline = FALSE}).
#' @param n_receivers An integer that defines the number of receivers in the array. This is ignored if receiver locations are specified via \code{arrangement}.
#' @param arrangement,... A character string or a SpatialPoints object that defines the arrangement of receivers. Supported character strings options for simulated arrays are \code{"regular"}, \code{"random"} and \code{"stratified"}, \code{"nonaligned"}, \code{"hexagonal"} and \code{"clustered"} (see \code{\link[sp]{spsample}}, which is used to simulate receiver locations). Additional arguments can be passed to this function via \code{...} for further control. Otherwise, a SpatialPoints object that defines the coordinates of receivers (in the same coordinate reference system as \code{boundaries} and, if applicable, \code{coastline}) is assumed to have been provided.
#' @param seed An integer that is used to set the seed to enable reproducible simulations (see \code{\link[base]{set.seed}}).
#' @param plot A logical variable that defines whether or not plot the array.
#' @param xlim,ylim (optional) Axis limits for the plot. These can be specified in any way supported by \code{\link[prettyGraphics]{pretty_axis}}.
#' @param add_sea (optional) If \code{plot = TRUE}, \code{add_sea} is a named list of arguments, passed to \code{\link[raster]{plot}}, to customise the appearance of the sea on the plot. \code{add_sea = NULL} suppresses the addition of the sea to the plot. To use the default graphical parameters, simply specify \code{add_sea = list()}.
#' @param add_land (optional) If \code{plot = TRUE}, \code{add_land} is a named list of arguments, passed to \code{\link[raster]{plot}}, to customise the appearance of the land on the plot. \code{add_sea = NULL} suppresses the addition of the land to the plot. To use the default graphical parameters, simply specify \code{add_sea = list()}.
#' @param add_receivers (optional) If \code{plot = TRUE}, \code{add_receivers} is a named list of arguments, passed to \code{\link[graphics]{points}}, to customise the appearance of receivers on the plot. \code{add_receivers = NULL} suppresses the addition of the receivers to the plot. To use the default graphical parameters, simply specify \code{add_receivers = list()}.
#' @param verbose A logical variable that defines whether or not to print messages to the console to relay function progress.
#'
#' @return The function returns a named list of (a) spatial objects that define the simulated array ('array') and (b) the arguments used to generate this array ('args'). The 'array' element if this list contains the following elements: 'boundaries', a \code{\link[raster]{Extent-class}} object that defines the boundaries of the area (as inputted); 'area', a \code{\link[sp]{SpatialPolygons-class}} object that defines the boundaries of the area; 'land' and 'sea' are \code{\link[sp]{SpatialPolygons-class}} or \code{\link[sp]{SpatialPolygonsDataFrame-class}} objects that define the land and sea respectively; and 'xy' is a \code{\link[sp]{SpatialPoints-class}} object that defines receiver locations. If \code{plot = TRUE}, the function also returns a plot.
#'
#' @examples
#' #### Example (1): Simulate an array using default parameters
#' # ... And force reproducible simulations by defining seed
#' seed <- 1
#' array <- sim_array(boundaries = raster::extent(-10, 10, -10, 10),
#'                    seed = 1)
#'
#' #### Example (2): Simulate coastline and customise plot
#' # ... via add_land and add_sea
#' array <- sim_array(boundaries = raster::extent(-10, 10, -10, 10),
#'                    coastline = "simple_random",
#'                    add_land = list(col = "darkgreen"),
#'                    add_sea = list(col = scales::alpha("skyblue", 0.2)),
#'                    seed = 1
#'                    )
#'
#' #### Example (3) Add custom coastline
#' array <- sim_array(boundaries = raster::extent(dat_coast),
#'                    coastline = dat_coast,
#'                    add_land = list(col = "darkgreen"),
#'                    add_sea = list(col = scales::alpha("skyblue", 0.2)),
#'                    seed = 1
#'                    )
#'
#' #### Example (4) Change the number of receivers
#' array <- sim_array(n_receivers = 5)
#' array <- sim_array(n_receivers = 25)
#'
#' #### Example (5) Change the arrangement of receivers
#' ## Explore different arrangements
#' array <- sim_array(n_receivers = 25, arrangement = "random")
#' array <- sim_array(n_receivers = 25, arrangement = "regular")
#' array <- sim_array(n_receivers = 25, arrangement = "clustered", nclusters = 5)
#' array <- sim_array(n_receivers = 25, arrangement = "stratified")
#' array <- sim_array(n_receivers = 25, arrangement = "nonaligned")
#' array <- sim_array(n_receivers = 25, arrangement = "hexagonal")
#' ## Force arrangements around coastline
#' # Simulated island
#' array <- sim_array(n_receivers = 25,
#'                    coastline = "simple_random",
#'                    arrangement = "regular",
#'                    add_land = list())
#' # Real coastline
#' array <- sim_array(boundaries = raster::extent(dat_coast),
#'                    n_receivers = 25,
#'                    coastline = dat_coast,
#'                    arrangement = "regular",
#'                    add_land = list())
#' ## Incorporate custom arrangements
#' # Define receiver locations as a SpatialPoints object with a UTM CRS
#' # ... to match other spatial datasets
#' proj_wgs84 <- sp::CRS("+init=epsg:4326")
#' proj_utm <- sp::CRS(paste("+proj=utm +zone=29 +datum=WGS84",
#'                           "+units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))
#' xy <- sp::SpatialPoints(dat_moorings[, c("receiver_long", "receiver_lat")],
#'                         proj_wgs84)
#' xy <- sp::spTransform(xy, proj_utm)
#' # Make array
#' array <- sim_array(boundaries = raster::extent(dat_coast),
#'                    coastline = dat_coast,
#'                    arrangement = xy,
#'                    add_land = list()
#'                    )
#'
#' @author Edward Lavender
#' @export
#'

sim_array <- function(boundaries = raster::extent(-10, 10, -10, 10),
                      coastline = NULL,
                      land_inside_coastline = TRUE,
                      n_receivers = 10,
                      arrangement = "random",
                      seed = NULL,
                      plot = TRUE,
                      xlim = NULL, ylim = NULL,
                      add_sea = NULL,
                      add_land = NULL,
                      add_receivers = list(),
                      verbose = TRUE,...
){

  #### Initiate function
  t_onset <- Sys.time()
  cat_to_console <- function(..., show = verbose) if(show) cat(paste(..., "\n"))
  cat_to_console(paste0("flapper::sim_array() called (@ ", t_onset, ")..."))
  if(!is.null(seed)) set.seed(1)

  #### Define area
  cat_to_console("... Defining area...")
  area <- methods::as(boundaries, "SpatialPolygons")
  if(!is.null(coastline)) {
    if(!is.character(coastline)) {
      area_crs <- raster::crs(coastline)
    } else if(!is.character(arrangement)) {
      area_crs <- raster::crs(arrangement)
    } else area_crs <- NA
  } else{
    area_crs <- NA
  }
  if(is.na(area_crs)) message("CRS of area is NA.")
  raster::crs(area) <- area_crs

  #### Define the land and sea
  ## If coastline has been specified, we will simulate or incorporate coastline
  if(!is.null(coastline)) {
    cat_to_console("... Incorporating coastline...")
    ## Simulate coastline
    if(is.character(coastline)) {
      if(coastline == "simple_random") {
        cat_to_console("... ... Simulating coastline...")
        # randomly sample a three points in area
        points <- sp::spsample(area, n = 3, type = "random")
        # convert to a spatial polygon
        land <- sp::Polygon(points)
        land <- sp::SpatialPolygons(list(sp::Polygons(list(land), ID = 1)))
        # cut the land out of the area to make the sea
        sea <- rgeos::gDifference(area, land)
      } else {
        stop("Input to 'coastline' is not supported.")
      }

      ## Or incorporate coastline
    } else {
      if(land_inside_coastline) {
        land <- coastline
        sea <- rgeos::gDifference(area, land)
      } else{
        land <- rgeos::gDifference(area, coastline)
        sea <- coastline
      }
    }

    ## Otherwise, the whole area is effectively sea
  } else{
    land <- NULL
    sea <- area
  }

  #### Simulate receiver arrangement
  cat_to_console("... Incorporating receivers...")
  if(is.character(arrangement)) {
    cat_to_console("... ... Simulating receivers...")
    rxy <- sp::spsample(sea, n = n_receivers, type = arrangement,...)
  } else {
    rxy <- arrangement
  }

  #### Plot array
  if(plot){
    cat_to_console("... Plotting array...")
    # Get pretty axes
    x <- lapply(list(area, land, sea, rxy), function(x) if(!is.null(x)) raster::extent(x)[1:2])
    x <- unlist(x)
    y <- lapply(list(area, land, sea, rxy), function(x) if(!is.null(x)) raster::extent(x)[3:4])
    y <- unlist(y)
    pretty_axis_args <- list(side = 1:4,
                             axis = list(list(),
                                         list(),
                                         list(labels = FALSE),
                                         list(labels = FALSE)),
                             control_sci_notation = list(magnitude = 16L, digits = 0)
                             )
    axis_param <- prettyGraphics::implement_pretty_axis_args(x = list(x, y),
                                                             pretty_axis_args = pretty_axis_args,
                                                             xlim = xlim,
                                                             ylim = ylim)
    # Draw background plot
    raster::plot(area,
                 xlim = axis_param[[1]]$lim, ylim = axis_param[[2]]$lim,
                 axes = FALSE, border = NA)
    # Add spatial layers
    if(!is.null(add_land)){
      add_land$x <- land
      add_land$add <- TRUE
      do.call(raster::plot, add_land)
    }
    if(!is.null(add_sea)){
      add_sea$x <- sea
      add_sea$add <- TRUE
      do.call(raster::plot, add_sea)
    }
    if(!is.null(add_receivers)){
      add_receivers$x <- rxy
      do.call(graphics::points, add_receivers)
    }
    # Add pretty axes
    prettyGraphics::pretty_axis(axis_ls = axis_param, add = TRUE)

  }

  #### Return outputs
  ## Define outputs
  cat_to_console("... Defining outputs...")
  set.seed(NULL)
  out <- list()
  out$array = list(boundaries = boundaries,
                   area = area,
                   land = land,
                   sea = sea,
                   xy = rxy)
  out$args = list(boundaries = boundaries,
                  coastline = coastline,
                  land_inside_coastline = land_inside_coastline,
                  n_receivers = n_receivers,
                  arrangement = "arrangement",
                  plot = plot,
                  add_sea = add_sea,
                  add_land = add_land,
                  add_receivers = add_receivers,
                  verbose = verbose,
                  dots = list(...))
  ## Return outputs
  t_end <- Sys.time()
  duration <- difftime(t_end, t_onset, units = "mins")
  cat_to_console(paste0("... flapper::sim_array() call completed (@ ", t_end, ") after ~", round(duration, digits = 2), " minutes."))
  return(out)

}

