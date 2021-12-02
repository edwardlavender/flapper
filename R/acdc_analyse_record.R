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
#### acdc_plot_record()

#' @title Plot time-specific maps from the AC/DC algorithm(s)
#' @description This function is used to plot time-specific maps from the AC/DC algorithm(s). To implement the function, an \code{\link[flapper]{acdc_record-class}} list from \code{\link[flapper]{ac}}, \code{\link[flapper]{dc}} or \code{\link[flapper]{acdc}} plus \code{\link[flapper]{acdc_simplify}} must be supplied, from which the results can be extracted and plotted for specified time steps. For each time step, the function extracts the necessary information; sets up a blank background plot using \code{\link[raster]{plot}} and \code{\link[prettyGraphics]{pretty_axis}} and then adds requested spatial layers to this plot. Depending on user-inputs, this will usually show a cumulative map of where the individual could have spent more or less time, summed from the start of the algorithm to each time point. Coastline, receivers and acoustic centroids (if applicable) can be added and customised and the finalised plots can be returned or saved to file.
#' @param record An \code{\link[flapper]{acdc_record-class}} object.
#' @param type A character that defines the plotted surface(s): \code{"map_cumulative"} plots the cumulative surface and \code{"map_timestep"} plots time step-specific surfaces.
#' @param plot An integer vector that defines the time steps for which to make plots. If \code{plot = NULL}, the function will make a plot for all time steps for which the necessary information is available in \code{archive}.
#' @param add_coastline (optional) A named list of arguments, passed to \code{\link[raster]{plot}}, to add a polygon (i.e., of the coastline), to the plot. If provided, this must contain an `x' element that contains the coastline as a spatial object (e.g., a SpatialPolygonsDataFrame: see \code{\link[flapper]{dat_coast}} for an example).
#' @param add_receivers (optional) A named list of arguments, passed to \code{\link[graphics]{points}}, to add points (i.e., receivers) to the plot. If provided, this must contain an `x' element that is a SpatialPoints object that specifies receiver locations (in the same coordinate reference system as other spatial data).
#' @param add_raster (optional) A named list of arguments, passed to \code{\link[fields]{image.plot}}, to plot the RasterLayer of possible locations that is extracted from \code{archive}.
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
#' proj_wgs84 <- sp::CRS("+init=epsg:4326")
#' proj_utm <- sp::CRS(paste("+proj=utm +zone=29 +datum=WGS84",
#'                           "+units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))
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

