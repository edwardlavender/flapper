######################################
######################################
#### check_class_acdc_record()

#' @title Check an object is of class \code{"acdc_record"}
#' @description This function checks whether or not an object is of class \code{\link[flapper]{acdc_record-class}}.
#' @param x An object.
#' @param arg A character that defines the name of a function argument.
#' @return If \code{x} is not an \code{\link[flapper]{acdc_record-class}} object, the function returns a helpful error message.
#' @author Edward Lavender
#' @keywords internal
#'

check_class_acdc_record <- function(x, arg = deparse(substitute(x))){
  if(!inherits(x, "acdc_record")){
    if(inherits(x, "acdc_archive")) {
      stop(paste0("'", arg, "' must be converted from class 'acdc_archive' to class 'acdc_record' via flapper::acdc_simplify() before implementing this function."), call. = FALSE)
    } else{
      stop(paste0("'", arg, "' must be of class 'acdc_record'."), call. = FALSE)
    }
  }
}


######################################
######################################
#### acdc_access_*()

#' @title Short-cuts to elements of an \code{"acdc_record"} object
#' @description These functions provide short-cuts to elements of an \code{\link[flapper]{acdc_record-class}} object.
#'
#' @param record An \code{\link[flapper]{acdc_record-class}} object.
#' @param type For \code{\link[flapper]{acdc_access_maps}}, \code{type} is a character that specifies whether or not to access time step-specific maps (\code{type = "map_timestep"}) or  cumulative maps (\code{type = "map_cumulative"}).
#' @param select (optional) For \code{\link[flapper]{acdc_access_maps}}, \code{select} is an integer vector that defines the cumulative time steps (i.e., accounting for both primary (acoustic) and secondary (archival) time steps) for which maps are required.
#'
#' @details
#' \itemize{
#'   \item \code{\link[flapper]{acdc_access_dat}} accesses the \code{record$dat} elements of an \code{\link[flapper]{acdc_record-class}} object for all time steps.
#'   \item \code{\link[flapper]{acdc_access_timesteps}} accesses the total number of time steps stored in an \code{\link[flapper]{acdc_record-class}} object, accounting for both primary (acoustic) and secondary (archival) time steps.
#'   \item \code{\link[flapper]{acdc_access_maps}} accesses all, or specified, maps from the \code{record$spatial} elements of an \code{\link[flapper]{acdc_record-class}} object.
#' }
#'
#' @return
#' \itemize{
#'   \item \code{\link[flapper]{acdc_access_dat}} returns a single dataframe for all time steps.
#'   \item \code{\link[flapper]{acdc_access_timesteps}} returns an integer that defines the total number of time steps.
#'   \item \code{\link[flapper]{acdc_access_maps}} returns a single list of time-step specific or cumulative maps for specified or all time steps.
#' }
#'
#' @examples
#' #### Example (1): acdc_access_dat()
#' acdc_access_dat(acdc_simplify(dat_acdc))
#'
#' #### Example (2): acdc_access_timesteps()
#' acdc_access_timesteps(acdc_simplify(dat_acdc))
#'
#' #### Example (3): acdc_access_maps()
#' acdc_access_maps(acdc_simplify(dat_acdc))
#'
#' @author Edward Lavender
#' @name acdc_access
NULL

#### acdc_access_dat()
#' @rdname acdc_access
#' @export

acdc_access_dat <- function(record){
  check_class_acdc_record(record)
  dat <- lapply(record$record, function(record_elm) record_elm$dat)
  dat <- do.call(rbind, dat)
  return(dat)
}

#### acdc_access_timesteps()
#' @rdname acdc_access
#' @export

acdc_access_timesteps <- function(record){
  check_class_acdc_record(record)
  max(record$record[[length(record$record)]]$dat$timestep_cumulative)
}

#### acdc_access_maps()
#' @rdname acdc_access
#' @export

acdc_access_maps <- function(record, type = c("map_timestep", "map_cumulative"), select = NULL){
  check_class_acdc_record(record)
  type <- match.arg(type)
  maps <- lapply(record$record, function(record_elm){
    lapply(record_elm$spatial, function(spatial_elm) spatial_elm[[type]])
  })
  maps <- unlist(maps)
  if(!is.null(select)) maps <- maps[select]
  return(maps)
}


######################################
######################################
#### acdc_plot_trace()

#' @title Plot AC* centroid dynamics
#' @description This function visually reconstructs the dynamics of an acoustic-centroid algorithm (i.e., \code{\link[flapper]{ac}} or \code{\link[flapper]{acdc}}).
#'
#' To implement the function, an \code{\link[flapper]{acdc_record-class}} object (\code{record}) from \code{\link[flapper]{ac}} or \code{\link[flapper]{acdc}} plus \code{\link[flapper]{acdc_simplify}} that defines the outputs of the AC* algorithm is required. A \code{\link[sp]{SpatialPointsDataFrame}} that defines receiver locations and a matrix that defines the daily operational status of each receiver are also required.
#'
#' For each time step, the function plots the location-probability surface, the receiver(s) at which the individual was detected and the acoustic centroids, illustrating how the expansion, contraction and intersection of acoustic centroids capture the boundaries of an individual's location through time.
#'
#' @param record A \code{\link[flapper]{acdc_record-class}} object from \code{\link[flapper]{ac}} or \code{\link[flapper]{acdc}} plus \code{\link[flapper]{acdc_simplify}}.
#' @param plot An integer vector that defines the time steps for which to make plots. If \code{plot = NULL}, the function will make a plot for all time steps for which the necessary information is available in \code{record}.
#' @param moorings A \code{\link[sp]{SpatialPointsDataFrame}} that defines receiver locations (in the Universe Transverse Mercator coordinate reference system) and receiver IDs (in a column named `receiver_id').
#' @param moorings_matrix A matrix that defines, for each day of the study (rows) and each receiver (columns), receivers' operational status (see \code{\link[flapper]{make_matrix_receivers}}).
#' @param add_raster A named list of arguments, passed to \code{\link[prettyGraphics]{add_sp_raster}}, to plot the location-probability surface.
#' @param add_receiver_1,add_receiver_2,add_receiver_3,add_receiver_n Named lists of arguments, passed to \code{\link[graphics]{points}}, to customise the appearance of receivers. \code{add_receiver_1} controls the appearance of the `current' receiver (at which the individual was last or has just been detected); \code{add_receiver_2} controls the appearance of the receiver at which the individual was next detected; \code{add_receiver_3} controls the appearance of the third receiver at which the individual was detected; and \code{add_receiver_n} controls the appearance of all remaining active receivers.
#' @param add_centroid_ap,add_centroid_an,add_centroid_b,add_centroid_c Named lists of arguments that control the appearance of acoustic centroids. (\code{centroid_ap} defines the boundaries of the individual's location from the perspective of its previous location; \code{centroid_an} defines the boundaries of the individual's location at the moment of detection from the perspective of the receiver that recorded the detection; \code{centroid_b} defines the boundaries of the individual's location from the perspective of the receiver at which the individual was next detected; and \code{centroid_c} defines the boundaries of the individual's location integrated across all of these perspectives; see \code{\link[flapper]{acdc_record-class}}.) \code{add_centroid_ap,add_centroid_an} and \code{add_centroid_b} are passed to \code{\link[raster]{lines,SpatialPolygons-method}} and \code{add_centroid_c} is passed to \code{\link[raster]{plot}}.
#' @param add_coastline A named list of arguments, passed to \code{\link[raster]{plot}}, to add coastline to each plot.
#' @param add_main A named list of arguments, passed to \code{\link[graphics]{mtext}}, to customise the appearance of the plot title. Default plot titles are structured as follows: `Map for t = cumulative time step (detection time step [intermediate time step] time stamp'. When the intermediate time step is one, the individual is detected. At subsequent intermediate time steps, acoustic centroids expand and contract.
#' @param ... Additional plot customisation options passed to \code{\link[prettyGraphics]{pretty_map}}.
#' @param par_param A named list of arguments, passed to \code{\link[graphics]{par}}, to set the plotting window.
#' @param png_param (optional) A named list of arguments, passed to \code{\link[grDevices]{png}}, to save plots to file. If supplied, the plot for each time step is saved separately. The `filename' argument should be the directory in which plots are saved. Plots are then saved as "1.png", "2.png" and so on. If supplied, \code{prompt} is ignored (see below).
#' @param prompt If \code{png_param} is not specified, \code{prompt} is a logical variable that defines whether or not to pause function execution between plots and between centroids within plots to facilitate interpretation.
#'
#' @return The function returns, for each time step, a plot of the location-probability surface and acoustic centroids.
#'
#' @examples
#' #### Prepare example AC algorithm outputs with spatial files
#'
#' ## Define example time series
#' id <- 25
#' acc <- dat_acoustics[dat_acoustics$individual_id == id, ]
#' acc$timestamp <- lubridate::round_date(acc$timestamp, "2 mins")
#' acc$key       <- paste0(acc$timestamp, "-", acc$receiver_id)
#' acc           <- acc[!duplicated(acc$key), ][1:20, ]
#'
#' ## Define receiver locations ('moorings' SPDF) and activity status matrix
#' # Receiver locations
#' proj_wgs84 <- sp::CRS(SRS_string = "EPSG:4326")
#' proj_utm   <- sp::CRS(SRS_string = "EPSG:32629")
#' xy <- sp::SpatialPoints(dat_moorings[, c("receiver_long", "receiver_lat")],
#'                         proj_wgs84)
#' xy <- sp::spTransform(xy, proj_utm)
#' moorings <- sp::SpatialPointsDataFrame(xy, data = dat_moorings)
#' # Daily activity status matrix
#' as_POSIXct <- function(x) as.POSIXct(paste0(x, "00:00:00"), tz = "UTC")
#' moorings_mat <- make_matrix_receivers(dat_moorings,
#'                                       delta_t = "days",
#'                                       as_POSIXct = as_POSIXct)
#'
#' ## Prepare grid
#' # We will use a regular, relatively high resolution grid,
#' # focused on a small area around the receivers at which the ID was detected
#' grid <- raster::raster(raster::extent(dat_gebco),
#'                        res = c(25, 25),
#'                        crs = raster::crs(dat_gebco))
#' grid <- raster::resample(dat_gebco, grid)
#' ext <-
#'   raster::extent(
#'     rgeos::gBuffer(moorings[moorings$receiver_id %in% acc$receiver_id, ],
#'                    width = 10000)
#'   )
#' grid <- raster::crop(grid, ext)
#' grid <- raster::trim(grid)
#' raster::plot(grid)
#'
#' ## Define detection centroids/probability kernels
#' # Define detection centroids
#' moorings <- raster::crop(moorings, grid)
#' dat_centroid <- acs_setup_centroids(xy = moorings,
#'                                     detection_range = 425,
#'                                     coastline = dat_coast,
#'                                     boundaries = ext,
#'                                     plot = TRUE,
#'                                     resolution = 10,
#'                                     verbose = TRUE)
#' # Define detection centroid overlaps
#' centroids_spdf      <- do.call(raster::bind, plyr::compact(dat_centroids))
#' centroids_spdf@data <- dat_moorings
#' dat_centroids_overlaps <-
#'   get_detection_centroids_overlap(centroids = centroids_spdf,
#'                                   services = NULL)
#' # Define detection probability kernels
#' calc_dpr <-
#'   function(x){
#'     ifelse(x <= 425, stats::plogis(2.5 + -0.02 * x), 0)
#'   }
#' dat_kernels <- acs_setup_detection_kernels(xy = moorings,
#'                                            services = NULL,
#'                                            centroids = dat_centroids,
#'                                            overlaps = dat_centroids_overlaps,
#'                                            calc_detection_pr = calc_dpr,
#'                                            bathy = grid)
#'
#' ## Implement AC algorithm
#' out_ac <- ac(acoustics = acc,
#'              step = 120,
#'              bathy = grid,
#'              detection_centroids = dat_centroids,
#'              detection_kernels = dat_kernels,
#'              detection_kernels_overlap = dat_centroids_overlaps,
#'              mobility = 200,
#'              save_record_spatial = NULL
#'              )
#'
#' ## Simplify outputs
#' record <- acdc_simplify(out_ac)
#'
#' #### Example (1): Implement the function with default arguments
#' if(interactive()){
#'   acdc_plot_trace(record, 1:10, moorings, moorings_mat)
#' }
#'
#' #### Example (2): Customise plot via add_* arguments and ...
#' if(interactive()){
#'   acdc_plot_trace(record, 1:10, moorings, moorings_mat,
#'                   add_raster =
#'                     list(plot_method = raster::plot, legend = FALSE),
#'                   add_coastline =
#'                     list(x = dat_coast, col = scales::alpha("dimgrey", 0.5)),
#'                   xlim = c(699000, 711000),
#'                   ylim = c(6250000, 6269500))
#' }
#'
#' #### Example (3): Save plots to file via the png_param argument
#' con <- paste0(tempdir(), "/acdc_trace/")
#' if(!dir.exists(con)) dir.create(con)
#' acdc_plot_trace(record, 1:10, moorings, moorings_mat,
#'                 add_raster =
#'                   list(plot_method = raster::plot, legend = FALSE),
#'                 add_coastline =
#'                   list(x = dat_coast, col = scales::alpha("dimgrey", 0.5)),
#'                 xlim = c(699000, 711000),
#'                 ylim = c(6250000, 6269500),
#'                 png_param = list(filename = con),
#'                 prompt = FALSE)
#' list.files(con)
#'
#' @seealso \code{\link[flapper]{ac}} and \code{\link[flapper]{acdc}} implement the AC and ACDC algorithms and \code{\link[flapper]{acdc_simplify}} simplifies the results. \code{\link[flapper]{acdc_plot_record}} and \code{\link[flapper]{acdc_animate_record}} provide additional visualisation routines.
#' @author Edward Lavender
#' @export

acdc_plot_trace <- function(record,
                            plot = NULL,
                            moorings, moorings_matrix,
                            add_raster = list(),
                            add_receiver_1 = list(pch = 4, lwd = 4, col = "darkgreen"),
                            add_receiver_2 = list(pch = 4, lwd = 2, col = "darkorange"),
                            add_receiver_3 = list(pch = 4, lwd = 1, col = "darkred"),
                            add_receiver_n = list(pch = 4, lwd = 2),
                            add_centroid_ap = list(col = "darkgreen"),
                            add_centroid_an = list(col = "darkgreen"),
                            add_centroid_b = list(col = "darkorange"),
                            add_centroid_c = list(col = scales::alpha("forestgreen", 0.5), density = 20),
                            add_coastline = list(),
                            add_main = list(),...,
                            par_param = list(),
                            png_param = list(),
                            prompt = TRUE){

  #### Check user inputs
  check_class_acdc_record(record)
  if(is.null(record$args))
    stop("Please add the list of ac() or acdc() arguments to 'record'.", call. = FALSE)
  check_class(input = moorings, to_class = "SpatialPointsDataFrame")
  check_names(input = moorings, req = "receiver_id", type = all)
  if(length(png_param) > 0L && prompt){
    warning("Either 'png_param' or 'prompt' should be supplied: ignoring 'prompt'.",
            immediate. = TRUE, call. = FALSE)
    prompt <- FALSE
  }
  save <- FALSE
  if(length(png_param) > 0L){
    check_named_list(input = png_param)
    check_names(input = png_param, req = "filename")
    png_param$filename <- check_dir(input = png_param$filename, check_slash = TRUE)
    png_filename <- png_param$filename
    save <- TRUE
  }

  #### Define objects
  # Helper functions
  cat_to_console <- function(..., show = TRUE) if(show) cat(paste(..., "\n"))
  continue <- function(.prompt) if(.prompt) readline(prompt = "Press [Enter] to continue or [Esc] to quit...")
  # The acoustic/archival dataframes (in their processed form, as used by the AC* algorithm)
  dat_record <- acdc_access_dat(record)
  if(all(!(c(dat_record$receiver_1_id, dat_record$receiver_2_id) %in% moorings$receiver_id)))
    stop("'moorings' does not contain all of the receivers found in 'record'.", call. = FALSE)
  # Define timestamps for 'intermediate' time steps if undefined
  if(!rlang::has_name(dat_record, "archival_timestamp")){
    step <- record$args$step
    dat_record <-
      dat_record %>%
      dplyr::group_by(.data$timestep_detection) %>%
      dplyr::mutate(archival_timestamp =
                      seq(min(.data$receiver_1_timestamp), max(.data$receiver_2_timestamp),
                          by = step)[.data$timestep_archival]) %>%
      data.frame()
  }
  # Get the spatial record as a contiguous record, with one element for each cumulative time step
  record_spatial <- lapply(record$record, function(record_elm){
    lapply(record_elm$spatial, function(spatial_elm) spatial_elm)
  })
  record_spatial <- purrr::flatten(record_spatial)
  record_spatial <- compact(record_spatial)
  # Define the number of time steps with spatial records (before a NULL element)
  pos_wo_spatial <- which(sapply(record$record, function(elm) length(elm$spatial) > 0L) == FALSE)
  if(length(pos_wo_spatial) > 1L){
    max_t <- min(pos_wo_spatial) - 2L
  } else{
    max_t <- nrow(dat_record) - 1L
  }
  time_steps <- 1:max_t
  if(!is.null(plot)){
    if(any(!(plot %in% time_steps))){
      warning("Not all indices in 'plot' can be plotted.",
              immediate. = TRUE, call. = FALSE)
      time_steps <- time_steps[time_steps %in% plot]
    }
  }

  #### Define global plotting param
  pp <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(pp), add = TRUE)
  if(!save){
    pb <- pbapply::pboptions(type = "none")
    on.exit(pbapply::pboptions(pb), add = TRUE)
  }
  if(is.null(add_main$side)) add_main$side <- 3

  #### Loop over cumulative time steps (with spatial data)
  cl_lapply(x = time_steps, fun = function(t){
    cat_to_console(paste0("On timestep: ", t, "..."), show = prompt)

    #### Define time specific parameters
    # t = 1
    detection_timestep_for_t  <- dat_record$timestep_detection[t]
    detection_timesteps_for_t <- detection_timestep_for_t:(detection_timestep_for_t + 2)
    dat_record_for_t <- dat_record[dat_record$timestep_detection %in% detection_timesteps_for_t, ]
    spatial_for_t    <- record_spatial[[t]]

    #### Define plotting param
    ## Set up plot to save
    if(save) {
      png_param$filename <- paste0(png_filename, t, ".png")
      do.call(grDevices::png, png_param)
    }
    ## Plotting window
    graphics::par(par_param)

    #### Plot map
    ## Base map
    cat_to_console("...Plot base map...", show = prompt)
    add_raster$x <- spatial_for_t$map_timestep
    axis_ls <-
      suppressMessages(
        prettyGraphics::pretty_map(add_rasters = add_raster,...)
      )
    ext <- raster::extent(c(axis_ls[[1]]$lim, axis_ls[[2]]$lim))
    if(length(add_coastline) > 0L){
      add_coastline$add <- TRUE
      add_coastline$x <- raster::crop(add_coastline$x, ext)
      do.call(raster::plot, add_coastline)
    }
    ## Add receivers
    # All active receivers
    active <- moorings_matrix[rownames(moorings_matrix) ==
                                as.character(as.Date(dat_record$archival_timestamp[t])), , drop = FALSE]
    add_receiver_n$x <-
      moorings[moorings$receiver_id %in% as.integer(colnames(active)[which(active == 1)]), ]
    do.call(graphics::points, add_receiver_n)
    # Add the receiver which recorded the detection, the next receiver and the receiver after that
    filter_moorings_by_detection_timestamp <-
      function(index)
        moorings[moorings$receiver_id ==
                   dat_record_for_t$receiver_id[dat_record_for_t$timestep_detection == detection_timesteps_for_t[index]],
        ]
    add_receiver_1$x <- filter_moorings_by_detection_timestamp(1)
    add_receiver_2$x <- filter_moorings_by_detection_timestamp(2)
    add_receiver_3$x <- filter_moorings_by_detection_timestamp(3)
    do.call(graphics::points, add_receiver_1)
    do.call(graphics::points, add_receiver_2)
    do.call(graphics::points, add_receiver_3)
    ## Add title
    # Define generic plot title with
    # ... cumulative time step
    # ... time step of the detection
    # ... time step of the archival step
    # ... the actual time step
    if(is.null(add_main$text)){
      main <- paste0("Map for t = ", t,
                     " (", dat_record$timestep_detection[t],
                     "[", dat_record$timestep_archival[t], "]: ",
                     format(dat_record$archival_timestamp[t], "%y-%m-%d %H:%M:%S"), ")")
      add_main$text <- main
    }
    do.call(graphics::mtext, add_main)

    #### Add centroids
    ## Add centroid (Ap)
    if(!is.null(spatial_for_t$centroid_ap)){
      cat_to_console("...Add centroid (An)...", show = prompt)
      add_centroid_ap$x <- spatial_for_t$centroid_ap
      add_centroid_ap$x <- raster::crop(add_centroid_ap$x, ext)
      do.call(raster::lines, add_centroid_ap)
      continue(prompt)
    }

    ## Add centroid (An)
    if(!is.null(spatial_for_t$centroid_an)){
      cat_to_console("...Add centroid (An)...", show = prompt)
      add_centroid_an$x <- spatial_for_t$centroid_an
      add_centroid_an$x <- raster::crop(add_centroid_an$x, ext)
      do.call(raster::lines, add_centroid_an)
      continue(prompt)
    }

    ## Add centroid (B)
    if(!is.null(spatial_for_t$centroid_b)){
      cat_to_console("...Add centroid (B)...", show = prompt)
      add_centroid_b$x <- spatial_for_t$centroid_b
      add_centroid_b$x <- raster::crop(add_centroid_b$x, ext)
      do.call(raster::lines, add_centroid_b)
      continue(prompt)
    }

    ## Add centroid (C)
    cat_to_console("...Add centroid (C)...", show = prompt)
    add_centroid_c$x <- spatial_for_t$centroid_c
    add_centroid_c$x <- raster::crop(add_centroid_c$x, ext)
    add_centroid_c$add <- TRUE
    do.call(raster::plot, add_centroid_c)
    continue(prompt)

    ## Add back coastline at end if necessary (for tidiness)
    if(save && length(add_coastline) > 0L){
      add_coastline$add <- TRUE
      do.call(raster::plot, add_coastline)
    }

    ## Save plot
    if(save) grDevices::dev.off()
    return(invisible())
  })

  #### Return blank
  return(invisible())

}


######################################
######################################
#### acdc_plot_record()

#' @title Plot time-specific maps from the AC/DC algorithm(s)
#' @description This function is used to plot time-specific maps from the AC/DC algorithm(s). To implement the function, an \code{\link[flapper]{acdc_record-class}} list from \code{\link[flapper]{ac}}, \code{\link[flapper]{dc}} or \code{\link[flapper]{acdc}} plus \code{\link[flapper]{acdc_simplify}} must be supplied, from which the results can be extracted and plotted for specified time steps. For each time step, the function extracts the necessary information; sets up a blank background plot using \code{\link[raster]{plot}} and \code{\link[prettyGraphics]{pretty_axis}} and then adds requested spatial layers to this plot. Depending on user-inputs, this will usually show a cumulative map of where the individual could have spent more or less time, summed from the start of the algorithm to each time point. Coastline, receivers and acoustic centroids (if applicable) can be added and customised and the finalised plots can be returned or saved to file.
#' @param record An \code{\link[flapper]{acdc_record-class}} object.
#' @param type A character that defines the plotted surface(s): \code{"map_cumulative"} plots the cumulative surface and \code{"map_timestep"} plots time step-specific surfaces.
#' @param plot An integer vector that defines the time steps for which to make plots. If \code{plot = NULL}, the function will make a plot for all time steps for which the necessary information is available in \code{record}.
#' @param add_coastline (optional) A named list of arguments, passed to \code{\link[raster]{plot}}, to add a polygon (i.e., of the coastline), to the plot. If provided, this must contain an `x' element that contains the coastline as a spatial object (e.g., a SpatialPolygonsDataFrame: see \code{\link[flapper]{dat_coast}} for an example).
#' @param add_receivers (optional) A named list of arguments, passed to \code{\link[graphics]{points}}, to add points (i.e., receivers) to the plot. If provided, this must contain an `x' element that is a SpatialPoints object that specifies receiver locations (in the same coordinate reference system as other spatial data).
#' @param add_raster (optional) A named list of arguments, passed to \code{\link[fields]{image.plot}}, to plot the RasterLayer of possible locations that is extracted from \code{record}.
#' @param add_centroids (optional) For outputs from the AC* algorithms (\code{\link[flapper]{ac}} or \code{\link[flapper]{acdc}}), \code{centroids} is a named list of arguments, passed to \code{\link[raster]{plot}}, to add the acoustic centroid to the plot.
#' @param add_additional (optional) A stand-alone function, to be executed after the background plot has been made and any specified spatial layers have been added to this, to customise the result (see Examples).
#' @param crop_spatial A logical variable that defines whether or not to crop spatial data to lie within the axis limits.
#' @param xlim,ylim,fix_zlim,pretty_axis_args Axis control arguments. \code{xlim} and \code{ylim} control the axis limits, following the rules of the \code{lim} argument in \code{\link[prettyGraphics]{pretty_axis}}. \code{fix_zlim} is a logical input that defines whether or not to fix z axis limits across all plots (to facilitate comparisons), or a vector of two numbers that define a custom range for the z axis which is fixed across all plots. \code{fix_zlim = FALSE} produces plots in which the z axis is allowed to vary flexibly between time units. Other axis options supported by \code{\link[prettyGraphics]{pretty_axis}} are implemented by passing a named list of arguments to this function via \code{pretty_axis_args}.
#' @param par_param (optional) A named list of arguments, passed to \code{\link[graphics]{par}}, to control the plotting window. This is executed before plotting is initiated and therefore affects all plots.
#' @param png_param (optional) A named list of arguments, passed to \code{\link[grDevices]{png}}, to save plots to file. If supplied, the plot for each time step is saved separately. The `filename' argument should be the directory in which plots are saved. Plots are then saved as "1.png", "2.png" and so on.
#' @param cl,varlist (optional) Parallelisation options. \code{cl} is (a) a cluster object from \code{\link[parallel]{makeCluster}} or (b) an integer that defines the number of child processes. \code{varlist} is a character vector of variables for export (see \code{\link[flapper]{cl_export}}). Exported variables must be located in the global environment. If a cluster is supplied, the connection to the cluster is closed within the function (see \code{\link[flapper]{cl_stop}}). For further information, see \code{\link[flapper]{cl_lapply}} and \code{\link[flapper]{flapper-tips-parallel}}. If supplied, the function loops over specified time steps in parallel to make plots. This is only implemented if plots are saved to file (i.e., \code{png_param} is supplied).
#' @param verbose A logical variable that defines whether or not relay messages to the console to monitor function progress.
#' @param check A logical variable that defines whether or not to check user inputs to the function before its initiation.
#' @param ... Additional arguments, passed to \code{\link[raster]{plot}}, to customise the blank background plot onto which spatial layers are added, such as \code{xlab}, \code{ylab} and \code{main}.
#'
#' @return The function plots the results of the AC* algorithm(s) at specified time steps, with one plot per time step. Plots are saved to file if \code{png_param} is supplied.
#' @examples
#' #### Step (1): Implement AC/DC algorithm(s)
#' # ... see examples via ac(), dc() and acdc()
#'
#' #### Step (2): Simplify outputs of the AC/DC algorithm(s)
#' dat_acdc <- acdc_simplify(dat_acdc)
#'
#' #### Example (1): The default options simply plot the first surface
#' acdc_plot_record(record = dat_acdc)
#'
#' #### Example (2): Define the number of plots to be produced and control the plotting window
#' acdc_plot_record(record = dat_acdc,
#'                  plot = 1:2,
#'                  par_param = list(mfrow = c(1, 2), mar = c(8, 8, 8, 8)))
#'
#' #### Example (3): Add and customise spatial information via add_* args
#' ## Define a SpatialPoints object of receiver locations
#' proj_wgs84 <- sp::CRS(SRS_string = "EPSG:4326")
#' proj_utm   <- sp::CRS(SRS_string = "EPSG:32629")
#' rsp <- sp::SpatialPoints(dat_moorings[, c("receiver_long", "receiver_lat")], proj_wgs84)
#' rsp <- sp::spTransform(rsp, proj_utm)
#' ## Plot with receiver locations and coastline, customise the centroids and the raster
#' acdc_plot_record(record = dat_acdc,
#'                  add_coastline = list(x = dat_coast, col = "darkgreen"),
#'                  add_receivers = list(x = rsp, pch = 4, col = "royalblue"),
#'                  add_centroids = list(col = "red"),
#'                  add_raster = list(col = rev(topo.colors(100)))
#'                  )
#'
#' #### Example (4): Control axis properties
#' # ... via smallplot argument for raster, pretty_axis_args, xlim, ylim and fix_zlim
#' # ... set crop_spatial = TRUE to crop spatial data within adjusted limits
#' acdc_plot_record(record = dat_acdc,
#'                  add_coastline = list(x = dat_coast, col = "darkgreen"),
#'                  add_receivers = list(x = rsp, pch = 4, col = "royalblue"),
#'                  add_centroids = list(col = "red"),
#'                  add_raster = list(smallplot= c(0.85, 0.9, 0.25, 0.75)),
#'                  crop_spatial = TRUE,
#'                  pretty_axis_args = list(side = 1:4,
#'                                          control_sci_notation =
#'                                            list(magnitude = 16L, digits = 0)),
#'                  xlim = raster::extent(dat_coast)[1:2],
#'                  ylim = raster::extent(dat_coast)[3:4],
#'                  fix_zlim = c(0, 1)
#'                  )
#'
#' #### Example (5): Modify each plot after it is produced via add_additional
#' # Specify a function to add titles to a plot
#' add_titles <- function(){
#'   mtext(side = 1, "x (UTM)", line = 2)
#'   mtext(side = 2, "y (UTM)", line = -8)
#' }
#' # Make plots with added titles
#' acdc_plot_record(record = dat_acdc,
#'                  plot = 1:2,
#'                  par_param = list(mfrow = c(1, 2), mar = c(8, 8, 8, 8)),
#'                  add_coastline = list(x = dat_coast, col = "darkgreen"),
#'                  add_receivers = list(x = rsp, pch = 4, col = "royalblue"),
#'                  add_centroids = list(col = "red"),
#'                  add_raster = list(),
#'                  crop_spatial = TRUE,
#'                  xlim = raster::extent(dat_coast)[1:2],
#'                  ylim = raster::extent(dat_coast)[3:4],
#'                  add_additional = add_titles
#'                  )
#'
#' #### Example (6) Save plots via png_param
#' list.files(tempdir())
#' acdc_plot_record(record = dat_acdc,
#'                  plot = 1:2,
#'                  png_param = list(filename = tempdir())
#'                  )
#' list.files(tempdir())
#'
#' #### Example (7) To plot the overall map, you can also just use a
#' # ... a raster plotting function like prettyGraphics::pretty_map()
#' ext <- update_extent(raster::extent(dat_coast), -1000)
#' prettyGraphics::pretty_map(x = ext,
#'                            add_rasters = list(x = dat_acdc$map),
#'                            add_points = list(x = rsp, pch = "*", col = "red"),
#'                            add_polys = list(x = dat_coast, col = "lightgreen"),
#'                            crop_spatial = TRUE,
#'                            xlab = "Easting", ylab = "Northing"
#'                            )
#'
#' @seealso This function is typically used following calls to \code{\link[flapper]{ac}}, \code{\link[flapper]{dc}} or \code{\link[flapper]{acdc}} and \code{\link[flapper]{acdc_simplify}}.
#' @author Edward Lavender
#' @export
#'

acdc_plot_record <- function(record,
                             type = c("map_cumulative", "map_timestep"),
                             plot = 1,
                             add_coastline = NULL,
                             add_receivers = NULL,
                             add_raster = list(col = rev(grDevices::terrain.colors(255))),
                             add_centroids = list(),
                             add_additional = NULL,
                             crop_spatial = FALSE,
                             xlim = NULL, ylim = NULL, fix_zlim = FALSE,
                             pretty_axis_args = list(side = 1:4,
                                                     axis = list(list(),
                                                                 list(),
                                                                 list(labels = FALSE),
                                                                 list(labels = FALSE)),
                                                     control_axis = list(las = TRUE),
                                                     control_sci_notation = list(magnitude = 16L, digits = 0)),
                             par_param = list(),
                             png_param = list(),
                             cl = NULL, varlist = NULL,
                             verbose = TRUE,
                             check = TRUE,...){

  #### Checks
  cat_to_console <- function(..., show = verbose) if(show) cat(paste(..., "\n"))
  cat_to_console("flapper::acdc_plot_record() called...")
  type <- match.arg(type)
  if(check){
    cat_to_console("... Checking function inputs...")
    ## Check object class
    check_class_acdc_record(record)
    ## Check plots to be produced
    if(any(plot <= 0L)) stop("Input to 'plot' must be > 0.")
    ## Check spatial data have been provided correctly
    if(!is.null(add_coastline)) {
      check_named_list(input = add_coastline)
      check_names(input = add_coastline, req = "x")
    }
    if(!is.null(add_receivers)) {
      check_named_list(input = add_receivers)
      check_names(input = add_receivers, req = "x")
      check_class(input = add_receivers$x, to_class = "SpatialPoints", type = "stop")
    }
    if(!is.null(add_raster)) check_named_list(input = add_raster)
    if(!is.null(add_centroids)) check_named_list(input = add_centroids)
    ## Check plotting window param
    check_named_list(input = par_param, ignore_empty = TRUE)
    ## Check png_param, if provided
    if(length(png_param) > 0){
      check_named_list(input = png_param)
      check_names(input = png_param, req = c("filename"), extract_names = names, type = all)
      png_param$filename <- check_dir(input = png_param$filename, check_slash = TRUE)
      save_png <- TRUE
    } else save_png <- FALSE
    ## Check cluster
    if(!is.null(cl) & is.null(png_param$filename)){
      message("Input to 'cl' ignored as png_param$filename is unspecified.")
      cl <- NULL
    }
    ## Check spatial information
    if(all(c(is.null(add_raster), is.null(add_coastline), is.null(add_receivers)))){
      stop("At least one argument out of 'add_raster', 'add_coastline' and 'add_receivers' should be provided.")
    }
    ## Other plotting param
    if(!is.null(pretty_axis_args$side)) {
      if(length(pretty_axis_args$side) == 1) stop("At least two sides in pretty_axis_args$side should be specified for a map.")
    }
    # Check dots
    check...("zlim",...)
  }

  #### Define data for background plot
  cat_to_console("... Defining data for background plot...")
  ## Unpack information for plotting and isolate relevant plots
  acdc_plot <- lapply(record$record, function(elm) elm$spatial)
  acdc_plot <- lapply(acdc_plot, function(x) if(length(x) > 0) return(x))
  acdc_plot <- purrr::flatten(acdc_plot)
  if(!is.null(plot)) acdc_plot <- acdc_plot[plot]
  acdc_plot <- acdc_plot[which(!sapply(acdc_plot, is.null))]
  if(is.null(plot)) plot <- 1:length(acdc_plot)
  if(check) {
    if(length(acdc_plot) <= 0) stop("No plotting data available for selected plot(s).")
    if(any(length(plot) > length(acdc_plot))) stop("'plot' exceeds the number of available plots.")
    if(!is.null(add_centroids)){
      if(!rlang::has_name(acdc_plot[[1]], "centroid")){
        add_centroids <- NULL
        message("add_centroids = NULL implemented: 'record' does not contain centroids.")
      }
    }
  }

  ## Extent of area
  if(!is.null(add_raster)) {
    first_raster <- acdc_plot[[1]][[type]]
    ext_ras <- raster::extent(first_raster)
    x_ras   <- ext_ras[1:2]
    y_ras   <- ext_ras[3:4]
  } else{
    x_ras <- NULL
    y_ras <- NULL
  }
  ## Extent of coastline provided
  if(!is.null(add_coastline)) {
    ext_coastline <- raster::extent(add_coastline$x)
    x_coastline   <- ext_coastline[1:2]
    y_coastline   <- ext_coastline[3:4]
  } else {
    x_coastline <- NULL
    y_coastline <- NULL
  }
  ## Extent of receiver locations
  if(!is.null(add_receivers)) {
    rxy <- add_receivers$x
    rxy <- sp::coordinates(rxy)
    x_rxy <- range(rxy[, 1])
    y_rxy <- range(rxy[, 2])
  } else {
    x_rxy <- NULL
    y_rxy <- NULL
  }
  ## Define x and and y values to the extremes of these possible limits
  x <- range(x_ras, x_coastline, x_rxy)
  y <- range(y_ras, y_coastline, y_rxy)
  ## Define axis limits
  axis_param <- prettyGraphics::implement_pretty_axis_args(x = list(x, y),
                                                           pretty_axis_args = pretty_axis_args,
                                                           xlim = xlim,
                                                           ylim = ylim)
  xlim <- axis_param[[1]]$lim
  ylim <- axis_param[[2]]$lim
  ext <- raster::extent(xlim, ylim)
  ## Crop spatial data within limits
  if(!is.null(add_coastline)) if(crop_spatial) add_coastline$x <- raster::crop(add_coastline$x, ext)

  ## Define zlim, if requested
  if(!is.null(add_raster)) {
    if(is.logical(fix_zlim)) {
      if(fix_zlim){
        range_use <- lapply(acdc_plot, function(map_info){
          min_use <- raster::minValue(map_info[[type]])
          max_use <- raster::maxValue(map_info[[type]])
          return(c(min_use, max_use))
        })
        range_use <- do.call(rbind, range_use)
        min_use <- min(range_use[, 1])
        max_use <- max(range_use[, 2])
        zlim <- c(min_use, max_use)
      }
    }
  }

  #### Define plotting window
  cat_to_console("... Setting plotting window...")
  pp <- do.call(graphics::par, par_param)

  #### Loop over every detection
  cat_to_console("... Making plots for each time step ...")
  cl_lapply(1:length(acdc_plot), cl = cl, varlist = varlist, fun = function(i){

    #### Set up image to save
    if(save_png){
      title <- paste0(i, ".png")
      png_param_tmp <- png_param
      png_param_tmp$filename <- paste0(png_param_tmp$filename, title)
      do.call(grDevices::png, png_param_tmp)
    }

    #### Define background plot
    map_info <- acdc_plot[[i]]
    map      <- map_info[[type]]
    area <- sp::SpatialPoints(raster::coordinates(ext), raster::crs(first_raster))
    raster::plot(x = area, xlim = xlim, ylim = ylim,...)

    #### Define time-specific zlim, if requested
    if(!is.null(add_raster)) {
      if(is.logical(fix_zlim)){
        if(!fix_zlim) {
          min_use <- raster::minValue(map)
          max_use <- raster::maxValue(map)
          zlim <- c(min_use, max_use)
        }
      } else {
        zlim <- fix_zlim
      }
    }

    #### Add spatial objects
    # Add spatial use surface
    if(!is.null(add_raster)) {
      add_raster$x <-  map
      add_raster$zlim <- zlim
      if(crop_spatial) add_raster$x <- raster::crop(add_raster$x, ext)
      add_raster$add <- TRUE
      do.call(fields::image.plot, add_raster)
    }
    # Add acoustic centroid
    if(!is.null(add_centroids)) {
      add_centroids$x   <- map_info$centroid
      if(crop_spatial) add_centroids$x <- raster::crop(add_centroids$x, ext)
      add_centroids$add <- TRUE
      do.call(raster::plot, add_centroids)
    }
    # Add the coastline (note that the coastline has already been cropped, if necessary)
    if(!is.null(add_coastline)) {
      add_coastline$add <- TRUE
      do.call(raster::plot, add_coastline)
    }
    # Add receivers
    if(!is.null(add_receivers)) {
      do.call(graphics::points, add_receivers)
    }
    # Add additional
    if(!is.null(add_additional)) {
      add_additional()
    }

    #### Add axes back at end
    prettyGraphics::pretty_axis(axis_ls = axis_param, add = TRUE)

    #### Save fig
    if(save_png) grDevices::dev.off()
  })

  return(invisible())

}


######################################
######################################
#### acdc_animate_record()

#' @title Create a html animation of the AC/DC algorithm(s)
#' @description This function is a simple wrapper for \code{\link[flapper]{acdc_plot_record}} and \code{\link[animation]{saveHTML}} which creates an animation of the AC* algorithm(s) over time. To implement this function, a named list of arguments for \code{\link[flapper]{acdc_plot_record}}, which creates the plots, must be supplied. This is embedded within \code{\link[animation]{saveHTML}}, which creates a folder in the specified directory named `images' that contains a .png file for each time step and an animation as a .html file.
#' @param expr_param A named list of arguments, passed to \code{\link[flapper]{acdc_plot_record}}, to create plots.
#' @param dir (optional) A string that defines the directory in which to save files. If unsupplied, if available, \code{dir} is taken from \code{html_name} using \code{\link[base]{dirname}}.
#' @param html_name A string that defines the name of the html file (see `htmlfile' argument in \code{\link[animation]{saveHTML}}).
#' @param image_name A string that defines the names of the individual .png files creates (see `img.name' argument in \code{\link[animation]{saveHTML}}).
#' @param html_title,html_description Character strings that provide a title and a description that are displayed within the html animation (see `title' and `description' arguments in \code{\link[animation]{saveHTML}}).
#' @param navigator A logical variable that defines whether or not to add a navigator panel to the animation (see `navigator' argument in \code{\link[animation]{saveHTML}}).
#' @param ani_height,ani_width,ani_res Numbers that define the size and the resolution of the animation (see `ani.height' `ani.width' and `ani.res' arguments in \code{\link[animation]{ani.options}}).
#' @param interval A number that defines the time interval between sequential frames (see `interval' argument in \code{\link[animation]{ani.options}}).
#' @param verbose A logical or character variable that defines whether or not, or what, to write as a footer to the html animation (see `verbose' argument in \code{\link[animation]{ani.options}}).
#' @param ... Additional arguments passed to \code{\link[animation]{ani.options}}.
#'
#' @return The function produces an animation in .html format in the specified directory. A folder named `images' is also produced which contains the images for each time step. The `css' and `js' folders are also produced by \code{\link[animation]{saveHTML}} which creates the animation.
#'
#' @examples
#' dir_current <- getwd()
#' setwd(tempdir())
#' acdc_record <- acdc_simplify(dat_acdc)
#' acdc_animate_record(expr_param =
#'                      list(record = acdc_record,
#'                           add_coastline = list(x = dat_coast, col = "darkgreen"),
#'                           plot = 1:5,
#'                           fix_zlim = FALSE)
#'                     )
#' setwd(dir_current)
#' @details This function requires the \code{\link[animation]{animation}} package.
#' @author Edward Lavender
#' @export
#'

acdc_animate_record <-
  function(expr_param,
           dir = NULL,
           html_name = "ACDC_algorithm_demo.html",
           image_name = "ACDC",
           html_title = "Demonstration of the ACDC Algorithm",
           html_description = "",
           navigator = FALSE,
           ani_height = 800,
           ani_width = 800,
           ani_res = 1200,
           interval = 0.1,
           verbose = FALSE,
           ...){
    #### Checks
    ## animation package
    if (!requireNamespace("animation", quietly = TRUE)) {
      stop("This function requires the 'animation' package. Please install it before continuing with install.packages('animation').")
    }
    #### Set directory
    if(is.null(dir)) dir <- dirname(html_name)
    wd <- getwd()
    check_dir(input = dir)
    setwd(dir)
    on.exit(setwd(wd), add = TRUE)
    html_name <- basename(html_name)
    #### Make plot
    animation::saveHTML({
      do.call(acdc_plot_record, expr_param)
    },
    htmlfile = html_name,
    img.name = image_name,
    title = html_title,
    description = html_description,
    navigator = navigator,
    ani.height = ani_height,
    ani.width = ani_width,
    ani.res = ani_res,
    interval = interval,
    verbose = verbose,...
    )
    return(invisible())
  }

