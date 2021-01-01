######################################
######################################
#### sim_array()

#' @title Simulate (marine) monitoring arrays
#' @description This function is designed to simulate different kinds of array designs for monitoring stations. The function has been particularly inspired by the need to simulate passive acoustic telemetry array designs, which comprise networks of acoustic hydrophones that listen for acoustic transmissions from tagged marine animals. To implement the function, it is necessary to define the boundaries of the area (\code{boundaries}). Barriers to movement, such as coastline, within this area can be simulated or included from real datasets. The area is populated with a specified number of receivers (\code{n_receivers}) that are simulated under different array designs (namely, uniform, regular or random arrangements) or incorporated from real data. The function returns a list of spatial objects that define the array and, if requested, a plot of the area.
#'
#' @param boundaries An \code{\link[raster]{extent}} object that defines the boundaries of (simulated) study area.
#' @param coastline (optional) This argument is used to incorporate the presence of barriers, such as coastline, in the area. There are three options. If \code{coastline = NULL}, no barriers are incorporated. If \code{coastline = "simple_random"}, then a triangular island is simulated in the study area. Alternatively, a spatial object, such as a SpatialPolygonsDataFrame, that defines the coastline in an area can be incorporated into the array design by passing this \code{coastline}.
#' @param land_inside_coastline A logical variable that defines whether or not the land is 'inside' the polygon(s) defined by \code{coastline} (\code{land_inside_coastline = TRUE}) or the sea is 'inside' the polygon(s) (\code{land_inside_coastline = FALSE}).
#' @param n_receivers An integer that defines the number of receivers in the array. This is ignored if receiver locations are specified via \code{arrangement}.
#' @param arrangement,... A character string or a SpatialPoints object that defines the arrangement of receivers. Supported character string options for simulated arrays are \code{"regular"}, \code{"random"} and \code{"stratified"}, \code{"nonaligned"}, \code{"hexagonal"} and \code{"clustered"} (see \code{\link[sp]{spsample}}, which is used to simulate receiver locations). Additional arguments can be passed to this function via \code{...} for further control. Otherwise, a SpatialPoints object that defines the coordinates of receivers (in the same coordinate reference system as \code{boundaries} and, if applicable, \code{coastline}) is assumed to have been provided.
#' @param seed An integer that is used to set the seed to enable reproducible simulations (see \code{\link[base]{set.seed}}).
#' @param plot A logical variable that defines whether or not plot the array.
#' @param xlim,ylim (optional) Axis limits for the plot. These can be specified in any way supported by \code{\link[prettyGraphics]{pretty_axis}}.
#' @param add_sea (optional) If \code{plot = TRUE}, \code{add_sea} is a named list of arguments, passed to \code{\link[raster]{plot}}, to customise the appearance of the sea on the plot. \code{add_sea = NULL} suppresses the addition of the sea to the plot. To use the default graphical parameters, simply specify \code{add_sea = list()}.
#' @param add_land (optional) If \code{plot = TRUE}, \code{add_land} is a named list of arguments, passed to \code{\link[raster]{plot}}, to customise the appearance of the land on the plot. \code{add_land = NULL} suppresses the addition of the land to the plot. To use the default graphical parameters, simply specify \code{add_land = list()}.
#' @param add_receivers (optional) If \code{plot = TRUE}, \code{add_receivers} is a named list of arguments, passed to \code{\link[graphics]{points}}, to customise the appearance of receivers on the plot. \code{add_receivers = NULL} suppresses the addition of the receivers to the plot. To use the default graphical parameters, simply specify \code{add_receivers = list()}.
#' @param verbose A logical variable that defines whether or not to print messages to the console to relay function progress.
#'
#' @return The function returns a named list of (a) the spatial objects that define the simulated array ('array') and (b) the arguments used to generate this array ('args'). The 'array' element is a named list contains the following elements: 'boundaries', an \code{\link[raster]{Extent-class}} object that defines the boundaries of the area (as inputted); 'area', a \code{\link[sp]{SpatialPolygons-class}} object that defines the boundaries of the area; 'land' and 'sea', \code{\link[sp]{SpatialPolygons-class}} or \code{\link[sp]{SpatialPolygonsDataFrame-class}} objects that define the land and sea respectively; and 'xy', a \code{\link[sp]{SpatialPoints-class}} object that defines receiver locations. If \code{plot = TRUE}, the function also returns a plot of the array.
#'
#' @examples
#' #### Example (1): Simulate an array using default parameters
#' # ... And force reproducible simulations by defining the seed
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
                      n_receivers = 10L,
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
  if(!is.null(seed)) set.seed(seed)

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


######################################
######################################
#### sim_path_sa()

#' @title Simulate discrete-time movement paths from step lengths and turning angles
#' @description This function simulates movement paths from step lengths and turning angles. To implement the function, the number of time steps (\code{n}) needs to be specified and, if applicable, the area within which movement should occur. For example, in marine environments, the inclusion of the sea as a spatial layer would restrict movement within the sea. The starting location (\code{p_1}) can be provided or simulated. At each time step, user-defined functions are used to simulate step lengths and turning angles, which can depend previous values of those variables via a \code{lag} parameter, from which the next position is calculated. This implementation enables movement paths to be simulated under a variety of movement models, including random walks and correlated random walks, providing that they are conceptualised in terms of step lengths and turning angles. The function returns a list of outputs that includes the simulated path and, if requested, produces a plot of the simulated path.
#'
#' @param n An integer that defines the number of time steps in the simulation.
#' @param area (optional) A \code{\link[sp]{SpatialPolygons-class}} or \code{\link[sp]{SpatialPolygonsDataFrame-class}} object that defines the area(s) within which movement is allowed.
#' @param p_1 (optional) A matrix with one row and two columns that defines the starting location (x, y). If \code{p_1 = NULL}, then a random location is sampled from \code{area}, if applicable, or simulated from a uniform distribution with a minimum and maximum value of 0 and 1 respectively.
#' @param sim_angle A function that is used to simulate turning angles. This must accept a single number that represents some previous turning angle (degrees), even if this is simply ignored (see \code{lag}, below). For example, \code{sim_angle = function() 1} will break but \code{sim_angle = function(...) 1} is fine. For convenience, a default function is included that simulates angles from a wrapped normal circular distribution with a mean and standard deviation of 1 (see \code{\link[circular]{rwrappednormal}}). Functions that actually depend on some previous angle also need to be able to generate initial angles before enough previous angles have been generated for the function to depend on those (see \code{lag}, below). All functions should return a single number that defines the turning angle in degrees.
#' @param sim_step A function that is used to simulate step lengths. This follows the same rules as for \code{sim_angle}. For convenience, a default function is included that simulates angles from a gamma distribution with shape and scale parameters of 15 (see \code{\link[stats]{rgamma}}).
#' @param lag If \code{sim_angle} and/or \code{sim_step} have been defined such that they depend on some previous angle/step length, then \code{lag} is an integer that defines the number of time steps between the current time step and some previous time step that affects the current turning angle and/or step length.
#' @param plot A logical variable that defines whether or not to produce a plot of the area (if provided) and the simulated movement path.
#' @param add_points,add_path (optional) Named lists of arguments that are used to customise the appearance of points (the starting location) and the path on the map.
#' @param seed (optional) An integer that defines the seed (for reproducible simulations: see \code{\link[base]{set.seed}}).
#' @param verbose A logical variable that defines whether or not to print messages to the console that relay function progress.
#' @param ... Additional arguments. For \code{\link[flapper]{sim_path_sa}}, these are passed to \code{\link[prettyGraphics]{pretty_map}} to customise the map. For the default \code{\link[flapper]{sim_angles}} and \code{\link[flapper]{sim_steps}} functions, \code{...} is required but additional parameters are ignored.
#'
#' @details This function requires the \code{\link[circular]{circular}} package.
#'
#' @return The function returns a named list of arguments that defines the simulated path ('xy_mat', 'angle_mat', 'step_mat' and 'path') and a named list of arguments that were used to generate the path ('args'). 'xy_mat' is an n-row, two-column matrix that defines the simulated position (x, y) at each time step; 'angle_mat' and 'step_mat' are n-row, one-column matrices that define the simulated turning angle (degrees) and step length (in map units) at each time step; and 'path' is a \code{\link[sp]{SpatialLines}} representation of the movement path.
#'
#' @seealso \code{\link[flapper]{sim_path_ou_1}} simulates a movement path based on past locations according to an Ornstein-Uhlenbeck process (which is not based on step lengths and turning angles).
#'
#' @examples
#' #### Example (1): Simulate movement path under default parameters
#' # Simulate path
#' path <- sim_path_sa()
#' # The function returns a list of parameters that define the array and a plot
#' summary(path)
#'
#' #### Example (2): Change the number of time steps
#' path <- sim_path_sa(n = 100)
#'
#' #### Example (3): Change the characteristics of the study area
#' # .. and define the starting location of the individual
#' sea  <- invert_poly(dat_coast)
#' path <- sim_path_sa(n = 100,
#'                     area = sea,
#'                     p_1 = matrix(c(706529.1, 6262293), ncol = 2),
#'                     add_polys = list(x = sea, col = "skyblue"))
#'
#' #### Example (4): Change the movement model(s) to use alternative distributions/parameters
#'
#' ## Step lengths
#' # Define new function to simulate step lengths
#' sim_step_lengths <- function(...) stats::rgamma(1, shape = 10, scale = 1)
#' # Check outputs suitable values
#' prettyGraphics::pretty_hist(replicate(n = 1000, expr = sim_step_lengths()))
#' # Implement simulation
#' path <- sim_path_sa(n = 100, sim_step = sim_step_lengths)
#' prettyGraphics::pretty_hist(as.numeric(path$step_mat))
#'
#' ## Turning angles
#' # E.g., Random walk: draw turning angle from von Mises distribution
#' sim_angles_vmd <- function(...){
#'   angle <- circular::rvonmises(n = 1,
#'                                mu = circular::circular(0),
#'                                kappa = 0,
#'                                control.circular = list(units = "degrees"))
#'   return(as.numeric(angle))
#' }
#' path <- sim_path_sa(n = 100, sim_angle = sim_angles_vmd)
#'
#' # E.g., Correlated random walk: draw turning angle from wrapped normal distribution
#' sim_angles_rwn <- function(...){
#'   angle <- circular::rwrappednormal(n = 1,
#'                                     mu = circular::circular(0),
#'                                     rho = 0.999,
#'                                     sd = 0,
#'                                     control.circular = list(units = "degrees"))
#'   return(as.numeric(angle))
#' }
#' path <- sim_path_sa(n = 100, sim_angle = sim_angles_rwn)
#'
#' #### Example (5) Change the movement models to depend on some lagged value
#' # ... of the variable in question
#' # Define a sim_angle function that depends on some previous angle
#' # While the time step is less than the lag, the function needs to be
#' # ... able to handle missing angles and return sensible values in these
#' # ... cases e.g., via is.null structure:
#' sim_angles_wrn_with_lag <- function(angle = NULL,...){
#'   if(is.null(angle)) {
#'     cat("\n... ... method (1) activated...\n") # useful check
#'     angle_out <- circular::circular(0)
#'   } else{
#'     angle_out <- circular::rwrappednormal(n = 1,
#'                                           mu = circular::circular(angle, units = "degrees"),
#'                                           rho = 0.9,
#'                                           sd = 0.1,
#'                                           control.circular = list(units = "degrees"))
#'   }
#'   return(as.numeric(angle_out))
#' }
#' # Check function
#' sim_angles_wrn_with_lag(NULL)
#' sim_angles_wrn_with_lag(1)
#' # Implement algorithm
#' path <- sim_path_sa(sim_angle = sim_angles_wrn_with_lag, lag = 1)
#' path <- sim_path_sa(sim_angle = sim_angles_wrn_with_lag, lag = 2)
#'
#' @author Edward Lavender
#' @name sim_path_sa
NULL

#' @rdname sim_path_sa
#' @export

sim_path_sa <- function(n = 10,
                       area = NULL,
                       p_1 = NULL,
                       sim_angle = sim_angles,
                       sim_step = sim_steps,
                       lag = 0L,
                       plot = TRUE,
                       add_points = list(pch = 21, bg = "darkgreen"),
                       add_path = list(length = 0.05, col = viridis::viridis(n)),
                       seed = NULL,
                       verbose = TRUE,...){

  #### Utils
  t_onset <- Sys.time()
  cat_to_console <- function(..., show = verbose) if(show) cat(paste(..., "\n"))
  cat_to_console(paste0("flapper::sim_path_sa() called (@ ", t_onset, ")..."))
  if(!requireNamespace("circular", quietly = TRUE)) stop("This function requires the 'circular' package. Please install it with `install.packages('circular')` first.")
  if(!is.null(seed)) set.seed(seed)
  out <- list(xy_mat = NULL, angle_mat = NULL, step_mat = NULL, time = NULL, args = NULL)
  out$time <- data.frame(event = "onset", time = Sys.time())
  if(!is.null(seed)) set.seed(seed)

  #### Define matrices to store movement path param
  cat_to_console("... Setting up simulation...")
  # Simulated step lengths
  step_mat <- matrix(NA, ncol = 1, nrow = n - lag)
  # Simulated angles
  angle_mat <- matrix(NA, ncol = 1, nrow = n - lag)
  # Simulated locations
  xy_mat <- matrix(NA, ncol = 2, nrow = n - lag)

  #### Simulate movement
  cat_to_console("... Simulating movement path...")
  pb <- utils::txtProgressBar(min = 0, max = n - lag, style = 3)
  for(t in 1:(n - lag)){

    #### Simulate starting position, if required
    if(t <= (lag + 1)){
      ## Generate starting angles and step lengths (required if lag > 0)
      step <- sim_step()
      angle <- sim_angle()
      ## Define positions
      if(is.null(p_1)) {
        if(!is.null(area)) {
          p_1 <- sp::spsample(area, n = 1, type = "random")
          p_1 <- sp::coordinates(p_1)
        } else {
          p_1 <- matrix(stats::runif(2, 0, 1), ncol = 2)
        }
      }
      px <- p_1[1]
      py <- p_1[2]

      #### Simulate movement from previous positions
    } else if(t > lag) {

      #### Simulation within an unbounded area
      if(is.null(area)) {

        ## Simulate step lengths and turning angles using user-supplied functions
        step <- sim_step(step_mat[t - lag])
        angle <- sim_angle(angle_mat[t - lag])

        ## Calculate new location
        dx <- step * cos(angle) # change in x
        dy <- step * sin(angle) # change in y
        px <- xy_mat[t-1, 1] + dx # new x position is previous x + change in x
        py <- xy_mat[t-1, 2] + dy # new y position is previous y + change in y

      } else {

        #### Simulation within a bounded area
        # ... We will repeatedly simulate the next position until
        # ... one that is inside the domain of interest is generated.

        repeat_count <- 0
        repeat{

          ## Simulate step lengths and turning angles using user-supplied functions
          step <- sim_step(step_mat[t - lag])
          angle <- sim_angle(angle_mat[t - lag])

          ## Calculate new location
          dx <- step * cos(angle) # change in x
          dy <- step * sin(angle) # change in y
          px <- xy_mat[t-1, 1] + dx # new x position is previous x + change in x
          py <- xy_mat[t-1, 2] + dy # new y position is previous y + change in y

          ## Identify whether the point is outside the domain
          psp <- sp::SpatialPoints(matrix(c(px, py), ncol = 2))
          raster::crs(psp) <- raster::crs(area)
          outside <- is.na(sp::over(psp, area))

          ## If the point is not outside of the area, then break the loop
          repeat_count <- repeat_count + 1
          if(!outside) break
        }
      }
    }

    #### Update matrices
    angle_mat[t, ] <- angle
    step_mat[t, ]  <- step
    xy_mat[t, ] <- c(px, py)
    utils::setTxtProgressBar(pb, t)

  }

  #### Define path as a SpatialLines object
  if(!is.null(area)) area_crs <- raster::crs(area) else area_crs <- NA
  path <- Orcs::coords2Lines(xy_mat, ID = 1)
  raster::crs(path) <- area_crs

  #### Plot
  if(plot) {
    cat_to_console("... Plotting simulated path...")
    if(!is.null(add_points)) add_points$x <- xy_mat[1, ]
    if(!is.null(add_path)) add_path$x <- xy_mat
    prettyGraphics::pretty_map(x = area,
                               add_points = add_points,
                               add_path = add_path,...)
  }

  #### Return outputs
  ## Define outputs
  set.seed(NULL)
  out <- list()
  out$xy_mat    <- xy_mat
  out$angle_mat <- angle_mat
  out$step_mat  <- step_mat
  out$path      <- path
  out$args = list(n = n,
                  area = area,
                  p_1 = p_1,
                  sim_angle = sim_angle,
                  sim_step = sim_step,
                  lag = lag,
                  plot = plot, add_points = add_points, add_path = add_path,
                  seed = seed,
                  verbose = verbose,
                  dots = list(...))
  ## Return outputs
  t_end <- Sys.time()
  duration <- difftime(t_end, t_onset, units = "mins")
  cat_to_console(paste0("... flapper::sim_path_sa() call completed (@ ", t_end, ") after ~", round(duration, digits = 2), " minutes."))
  return(out)
}

#' @rdname sim_path_sa
#' @export
sim_steps <- function(...) stats::rgamma(1, shape = 15, scale = 15)

#' @rdname sim_path_sa
#' @export
sim_angles <- function(...) {
  angle <- circular::rwrappednormal(n = 1,
                                    mu = circular::circular(0),
                                    rho = NULL,
                                    sd = 1,
                                    control.circular = list(units = "degrees"))
  return(as.numeric(angle))
}


######################################
######################################
#### sim_path_ou_1()

#' @title Simulate discrete-time movement paths under a Ornstein-Uhlenbeck process (1)
#' @description This function simulates movement paths under a discrete-time Ornstein-Uhlenbeck process in which the parameters of the movement model are assumed to remain constant through time.
#'
#' @param n An integer that defines the number of time steps in the simulation.
#' @param r_1 (optional) A matrix with one row and two columns that defines the starting location (x, y).  If \code{r_1 = NULL}, then a random location is simulated from a uniform distribution with a minimum and maximum value of 0 and 1 respectively.
#' @param r_h  A matrix with one row and two columns that defines the 'home range' centre location (x, y).
#' @param k A number that defines the strength of the 'central harmonic force' that pulls an individual back towards its home range centre.
#' @param delta_t A number that defines the number of time units between each time step.
#' @param eps A number that defines the variance.
#' @param plot A logical variable that defines whether or not to plot the simulated path.
#' @param add_paths,add_points (optional) Named lists of arguments that are used to customise the appearance of points (the home range and starting location, in that order, shown as filled green and black circles by default) and the path on the map.
#' @param verbose A logical variable that defines whether or not to print messages to the console that relay function progress.
#' @param ... Additional arguments, passed to \code{\link[prettyGraphics]{pretty_map}}, to customise the map.
#'
#' @details Ornstein-Uhlenbeck processes are convenient models for animal movement around a 'home range' centre. In the model, a parameter (\code{k}) 'pulls' the location of the individual (\eqn{\vec{r}}) back towards the centre of its home range (\eqn{\vec{r}^H}) as it moves away from this centre. This function implements the discretised form of the Ornstein-Uhlenbeck model in which the parameters of the movement model remain constant through time and in which movement is not constrained by barriers described by Alós et al. (2016) (see equations (8) and (9) in particular). Under this model, the position \eqn{\vec{r}} of the animal at time \eqn{n + 1} is given by:
#'
#' \deqn{\vec{r}_{n + 1} = \vec{r}^H + e^{-k \Delta t} (\vec{r}_n - \vec{r}^H) + \vec{R}_n,}
#'
#' where \eqn{\vec{r}^H} is the location of the home range centre; \eqn{k} is the strength of the central harmonic force; \eqn{\Delta t} is the duration between time steps; and \eqn{\vec{R}_n} is a bi-dimensional normally distributed random variable with mean zero and standard deviation (\eqn{\sigma}) given by:
#'
#' \deqn{\sigma = \sqrt{ \frac{\epsilon (1 - e^{-2 k \Delta t}}) {2 k}}.}
#'
#' Note that the default plotting parameters for this function require the \code{\link[viridis]{viridis}} package for pretty visualisation.
#'
#' @return The function returns an n-row, two-column matrix that defines the simulated location (x, y) at each time step and, if \code{plot = TRUE}, a plot of the path.
#'
#' @seealso \code{\link[flapper]{sim_path_sa}} simulates a movement path based on step lengths and turning angles. This can support movement within restricted areas.
#'
#' @examples
#' #### Example (1): Implement simulation with default options
#' path <- sim_path_ou_1()
#'
#' #### Example (2): Change the number of time steps
#' path <- sim_path_ou_1(n = 10000)
#'
#' #### Example (3): Change model parameters
#' # esp parameter
#' path <- sim_path_ou_1(n = 10000, eps = 10)
#' path <- sim_path_ou_1(n = 10000, eps = 500)
#' # central harmonic parameter
#' path <- sim_path_ou_1(n = 1000, eps = 1, k = 1)
#' path <- sim_path_ou_1(n = 1000, eps = 1, k = 3)
#' path <- sim_path_ou_1(n = 1000, eps = 1, k = 500)
#'
#' #### Example (4): Customise the plot via add_paths, add_points and ...
#' n <- 1000
#' path <- sim_path_ou_1(n = n,
#'                       add_points = list(pch = c(1, 4), lwd = 3),
#'                       add_paths = list(col = viridis::magma(n)),
#'                       pretty_axis_args = list(1:4)
#'                       )
#'
#' @references Alós, J. et al. (2016) Bayesian State-Space Modelling of Conventional Acoustic Tracking Provides Accurate Descriptors of Home Range Behavior in a Small-Bodied Coastal Fish Species. Plos One 11, e0154089
#'
#' @author Edward Lavender
#' @export

sim_path_ou_1 <-
  function(n = 1000,
           r_1 = NULL,
           r_h = matrix(c(0, 0), ncol = 2),
           k = 1,
           delta_t = 1,
           eps = 1,
           plot = TRUE,
           add_paths = list(length = 0.04, col = viridis::viridis(n)),
           add_points = list(pch = 21, cex = 2, lwd = 2, bg = c("darkgreen", "black")),
           verbose = TRUE,...){

    #### Initiate function
    t_onset <- Sys.time()
    cat_to_console <- function(..., show = verbose) if(show) cat(paste(..., "\n"))
    cat_to_console(paste0("flapper::sim_path_ou_1() called (@ ", t_onset, ")..."))

    #### Define matrices to hold results
    cat_to_console("... Setting up simulation...")
    xy_mat <- matrix(NA, nrow = n, ncol = 2)
    Rn_mat <- matrix(NA, nrow = n, ncol = 2)

    #### Simulate starting position
    if(is.null(r_1)) r_1 <- matrix(stats::runif(2, 0, 1), ncol = 2)
    xy_mat[1, ] <- r_1

    #### Simulate Rn parameter (uncorrelated 'errors' that relate to extent of movement)
    sd <- sqrt((eps * (1 - exp(-2 * k * delta_t)))/(2 * k))
    Rn_mat[ , ] <- stats::rnorm(n*2, 0, sd)

    #### Simulate process
    cat_to_console("... Simulating movement path...")
    # Define constant exp term for efficiency
    exp_term <- exp(-k*delta_t)
    # For each time step, simulate the position based on the OU process
    for(t in 1:(n - 1)){
      xy_mat[t+1, ] <- r_h + exp_term * (xy_mat[t, ] - r_h) + Rn_mat[t, ]
    }

    #### Visualise the simulated path
    if(plot) {
      cat_to_console("... Plotting simulated path...")
      add_paths$x <- xy_mat
      add_points$x <- matrix(rbind(r_h, xy_mat[1, ]), ncol = 2)
      prettyGraphics::pretty_map(add_paths = add_paths,
                                 add_points = add_points,
                                 verbose = verbose,...)
    }

    #### Return the path
    t_end <- Sys.time()
    duration <- difftime(t_end, t_onset, units = "mins")
    cat_to_console(paste0("... flapper::sim_path_ou_1() call completed (@ ", t_end, ") after ~", round(duration, digits = 2), " minutes."))
    return(xy_mat)

  }

