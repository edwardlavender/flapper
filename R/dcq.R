#' @title The `quick' depth-contour (DCQ) algorithm
#' @description This function implements the `quick' depth-contour (DCQ) algorithm. As for the DC algorithm (\code{\link[flapper]{dc}}), the DCQ algorithm relates one-dimensional depth time series to a two-dimensional bathymetry surface to determine the extent to which different parts of an area might have (or have not) been used, or effectively represent occupied depths, over time. However, a modified (binning) approach is used that is much faster and normalisation is not implementedd. To implement the function, a list of depth time series, one for each time unit (e.g. month/season) need to be supplied, along with a bathymetry \code{\link[raster]{raster}}. For each time unit, the function counts the number of depth observations in each user-defined depth bin (e.g., 10 m depth bins) and then relates these counts to the local bathymetry to produce a raster in which the value of each cell is given by the number of times in which the depth bin for that cell was used. The function returns a list of rasters, one for each time unit, and a plot of these rasters, if requested.
#'
#' @param archival_ls A list of dataframes, with one element for each time unit (e.g., month), which contain depth time series to be related to the local bathymetry. Each dataframe should contain a column of depths (`depth') and a column that defines the time unit (`time_unit'). Depth should be recorded using absolute values in the same units as the bathymetry (see below).
#' @param bathy A \code{\link[raster]{raster}} of the bathymetry in an area within which the animal is likely to have been located over the study. Bathymetry values should be recorded as absolute values and in the same units as for depths (see \code{archival_ls}).
#' @param bin A number that defines the interval into which depth observations are binned. This should be provided in the same units as depths (see \code{archival_ls}) and the bathymetry (see \code{bathy}). The appropriate value of \code{bin} depends on the measurement error of the \code{bathy} data and the depth time series, the tidal range in an area, computational requirements, and scale of biological research objectives. For large depth time series and/or large, high-resolution bathymetry rasters, it can be useful to test the algorithm's speed using a relatively large bin.
#' @param transform (optional) A function, such as \code{sqrt}, to transform counts. This affects the returned rasters and any plots produced (see Value). Be careful with some functions, such as \code{log}, which can generate problematic outputs (such as z axis limits if these are not defined manually), in some situations (for example, if some cells in the area are not visited).
#' @param plot A logical input that defines whether or not to plot the rasters. If \code{plot = TRUE}, the function produces a plot for each time unit.
#' @param before_plot,after_plot (optional) Stand-alone functions that are executed before and after the plot for each time unit is created, respectively. For example, it may be useful to plot the coast in an area before each raster is plotted, or add custom axes after each plot has been produced.
#' @param fix_zlim,one_page,... (optional) Plot customisation options. \code{fix_zlim} is a logical input that defines whether or not to fix z axis limits across all plots (to facilitate comparisons), or a vector of two numbers that define a custom range for the z axis which is fixed across all plots. \code{fix_zlim = FALSE} produces plots in which the z axis is allowed to vary flexibly between time units. \code{one_page} is a logical input that defines whether or not to produce all plots on one page; this is only implemented if there are fewer than 25 time units, beyond which there are typically to many plots to fit on one page. Additional plot customisation arguments can be passed to \code{\link[fields]{image.plot}} via \code{...}.
#' @param cl,varlist (optional) Parallelisation options. \code{cl} is (a) a cluster object from \code{\link[parallel]{makeCluster}} or (b) an integer that defines the number of child processes. \code{varlist} is a character vector of variables for export (see \code{\link[flapper]{cl_export}}). Exported variables must be located in the global environment. If a cluster is supplied, the connection to the cluster is closed within the function (see \code{\link[flapper]{cl_stop}}). For further information, see \code{\link[flapper]{cl_lapply}} and \code{\link[flapper]{flapper-tips-parallel}}.
#' @param verbose A logical input that defines whether or not to relay messages to the console to monitor function progress.
#'
#' @return The function returns a named list of rasters, one for each time unit, in which the value of each cell is the number of times that that cell was represented by the corresponding depth bin in the depth time series.
#'
#' @seealso This is a modified version of the DC algorithm implemented by \code{\link[flapper]{dc}}. The ACDC algorithm (see \code{\link[flapper]{acdc}}) extends the depth-contour algorithm by integrating information from acoustic detections of individuals at each time step to restrict the locations in which depth contours are identified.
#'
#' @examples
#' #### Define data for examples
#' # Define archival time series with required columns ('depth' and 'time_unit')
#' dat_archival <- dat_archival[order(dat_archival$timestamp), ]
#' dat_archival$time_unit <- cut(dat_archival$timestamp, "weeks")
#' # Define a list of dataframes with one element for each time unit
#' archival_ls <- split(dat_archival, f = dat_archival$time_unit)
#' # Define bathymetry data (and coastline data for plotting)
#' bathy <- prettyGraphics::dat_gebco
#' coastline <- prettyGraphics::dat_coast_around_oban
#'
#' #### Example (1) Implement the dcq() algorithm with 25 m bins
#' dcq_maps <- dcq(archival_ls = archival_ls,
#'               bathy = bathy,
#'               bin = 25,
#'               plot = FALSE)
#' # The function returns a list of rasters, with one raster
#' # ... for each time unit.
#' dcq_maps
#'
#' #### Example (2): Implement the algorithm in parallel:
#' dcq_maps <- dcq(archival_ls = archival_ls,
#'               bathy = bathy,
#'               bin = 25,
#'               plot = FALSE,
#'               cl = parallel::makeCluster(2L))
#'
#' #### Example (3): Visualise the function outputs on one page
#' # ... using standard options.
#' # Examine results with 25 m bin
#' dcq_maps <- dcq(archival_ls = archival_ls,
#'               bathy = bathy,
#'               bin = 25,
#'               plot = TRUE,
#'               one_page = TRUE)
#' # Examine results with a higher resolution bin
#' dcq_maps <- dcq(archival_ls = archival_ls,
#'               bathy = bathy,
#'               bin = 5,
#'               plot = TRUE,
#'               one_page = TRUE)
#'
#' #### Example (4): Plot customisation options
#' # fix zlim to be constant across all plots to enable comparability
#' dcq_maps <- dcq(archival_ls = archival_ls,
#'               bathy = bathy,
#'               bin = 5,
#'               plot = TRUE,
#'               one_page = TRUE,
#'               fix_zlim = TRUE)
#' # fix zlim using custom limits across all plots
#' dcq_maps <- dcq(archival_ls = archival_ls,
#'               bathy = bathy,
#'               bin = 5,
#'               plot = TRUE,
#'               one_page = TRUE,
#'               fix_zlim = c(0, 5000))
#' # Transform the returned and plotted rasters by supplying a function to the
#' # ... transform argument
#' dcq_maps <- dcq(archival_ls = archival_ls,
#'               bathy = bathy,
#'               bin = 5,
#'               plot = TRUE,
#'               one_page = TRUE,
#'               transform = sqrt)
#' # Customise the plot further via before_plot, after_plot functions
#' # ... and other arguments passed via ... E.g., note the need to include
#' # ... add = TRUE because the raster plot is added to the plot of the coastline.
#' dcq_maps <- dcq(archival_ls = archival_ls,
#'               bathy = bathy,
#'               bin = 5,
#'               plot = TRUE,
#'               one_page = TRUE,
#'               transform = sqrt,
#'               fix_zlim = FALSE,
#'               before_plot = function(x) raster::plot(coastline),
#'               after_plot = function(x) raster::lines(coastline),
#'               add = TRUE,
#'               col = topo.colors(100))
#'
#' @export
#' @author Edward Lavender
#'

dcq <- function(archival_ls,
               bathy,
               bin = 10, transform = NULL,
               plot = TRUE, before_plot = NULL, after_plot = NULL, fix_zlim = FALSE, one_page = FALSE,
               cl = NULL, varlist = NULL,
               verbose = TRUE,...){

  #### Checks
  cat_to_console <- function(..., show = verbose) if(show) cat(paste(..., "\n"))
  cat_to_console("flapper::dcq() called...")
  cat_to_console("... Step 1: Checking user inputs...")
  check_class(input = archival_ls, to_class = "list", type = "stop")
  sapply(1:length(archival_ls), function(i){
    check_names(arg = paste0("archival_ls[[", i, "]]"),
                input = archival_ls[[i]],
                req = c("depth", "time_unit")
                )
  })

  #### Define the frequency of use of different depth bins
  # Define a list called 'use_freq_by_time_unit'
  # ... with one element for each element of 'archival_ls' (i.e., each time unit)
  # ... which will contain a dataframe that, for each bin, defines how many times
  # ... cells in that depth bin could have been used.
  cat_to_console("... Step 2: Calculating the number of observations within each depth bin for each time unit...")
  max_depth <- max(sapply(archival_ls, function(d) max(d$depth, na.rm = TRUE)))
  use_freq_by_time_unit <-
    pbapply::pblapply(archival_ls, cl = NULL, function(df){
      ## Define histogram breaks from from 0 to the max depth by the size of the depth bin specified.
      breaks <- seq(0, max_depth, by = bin)
      ## Because we've specified a regular sequence, the maximum value of the break
      # ... might be less than the maximum depth, which will cause errors when we create the histogram,
      # ... so, if this is the case, we'll add an extra interval:
      if(max(breaks) < max_depth) breaks <- c(breaks, breaks[length(breaks)] + bin)
      ## Create a histogram of counts of observations at each depth in a vector of breaks
      # Use right = TRUE so the frequency in a bin defined by lower and upper
      # ... values x and y is over the interval >= x but < y.
      h <- graphics::hist(df$depth, breaks = breaks, right = TRUE, plot = FALSE)
      # Create a dataframe that includes the middle value of each depth bin
      # ... and the corresponding count of depth records in that bin.
      use_freq <- data.frame(time_unit = df$time_unit[1], mids = h$mids, counts = h$counts)
      # Add the lower and upper values associated with each bin
      use_freq$lower <- use_freq$mids - bin/2
      use_freq$upper <- use_freq$mids + bin/2
      # Return the dataframe
      return(use_freq)
    })

  #### Define a function to get the cell numbers of cells whose value lies within a specified range
  .cells_from_val <- function(y){
    return(cells_from_val(x = bathy, y = y, interval = 1L, cells = TRUE, na.rm = TRUE))
  }

  #### Use depth data to create rasters describing possible patterns in space use:
  # Create a list of objects, one for each time unit,
  # ... one element of which will be a raster of the
  # ... potential number of times in which each cell could have been used
  # ... based on the frequency with which that depth was visited:
  cat_to_console("... Step 3: Translating counts of observations within depth bins into maps...")
  area_use <- bathy
  area_use <- raster::setValues(area_use, 0)
  area_use_ls <- cl_lapply(use_freq_by_time_unit,
                           cl = cl,
                           varlist = varlist,
                           fun = function(use_freq){
    # use_freq <- use_freq_by_time_unit[[1]]
    # Determine the cell IDs of cells which lie within the lower and upper values
    # ... of each depth bin:
    cells_by_interval <- lapply(split(use_freq[, c("lower", "upper")], 1:nrow(use_freq)), FUN = function(y){
      cells <- .cells_from_val(as.numeric(y))
      return(cells)
    })
    # Update the area use raster in each of these cells based on the number of times that depth bin was used:
    # Loop over every element in cells_by_interval (i.e. every depth bin...)
    for(i in 1:length(cells_by_interval)) {
      area_use[cells_by_interval[[i]]] <- use_freq$counts[i]
    }
    # Transform the raster, if required
    if(!is.null(transform)) area_use <- transform(area_use)
    # Return a list of objects
    ls <- list(cells_by_interval = cells_by_interval,
               area_use = area_use)
    return(ls)
  })

  #### Visualise map of depth/space use for each time unit
  if(plot){

    ## Define plotting area
    cat_to_console("... Step 4: Mapping the results...")
    if(one_page) {
      if(length(area_use_ls) > 25) {
        message("The number of time units > 25: ignoring one_page = TRUE...")
        one_page <- FALSE
      }
    }
    if(one_page) {
      par_param <- graphics::par(no.readonly = TRUE)
      pp <- graphics::par(mfrow = prettyGraphics::par_mf(length(area_use_ls)))
      on.exit(graphics::par(par_param), add = TRUE)
    }

    ## Define zlim, if requested
    if(is.logical(fix_zlim)) {
      if(fix_zlim){
        range_use <- lapply(area_use_ls, function(time_unit){
          area_use <- time_unit$area_use
          min_use <- raster::cellStats(area_use, stat = "min", na.rm = TRUE)
          max_use <- raster::cellStats(area_use, stat = "max", na.rm = TRUE)
          return(c(min_use, max_use))
        })
        range_use <- do.call(rbind, range_use)
        min_use <- min(range_use[, 1])
        max_use <- max(range_use[, 2])
        zlim <- c(min_use, max_use)
      }
    }

    ## Loop over each time unit and product a plot
    lapply(area_use_ls, function(time_unit) {
      # Isolate raster
      area_use <- time_unit$area_use
      # Initial plot (e.g., plot coastline)
      if(!is.null(before_plot)) before_plot()
      # Define time-specific zlim, if requested
      if(is.logical(fix_zlim)){
        if(!fix_zlim) {
          min_use <- raster::cellStats(area_use, stat = "min", na.rm = TRUE)
          max_use <- raster::cellStats(area_use, stat = "max", na.rm = TRUE)
          zlim <- c(min_use, max_use)
        }
      } else {
        zlim <- fix_zlim
      }
      # Create map
      fields::image.plot(area_use, zlim = zlim,...)
      # Updates (e.g., re-add coastline)
      if(!is.null(after_plot)) after_plot()
    })
  }

  #### Return outputs
  out <- lapply(area_use_ls, function(elm) elm$area_use)
  return(out)

}
