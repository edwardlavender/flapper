#' @title Query the Open Topo Data API for elevation data
#' @description This function queries the Open Topo Data elevation Application Programming Interface (API) to extract elevation data for inputted coordinates, rasters, or areas.
#'
#' @param x A two-column matrix of coordinates (x, y), an \code{\link[raster]{extent}} object or a \code{\link[raster]{raster}} that defines the coordinates/area for which elevation data are desired. If a matrix is supplied, a query is sent for all pairs of coordinates in this matrix. If an \code{\link[raster]{extent}} object is provided, a query is sent for all cells spanning this area, depending on the specified dimensions (see \code{...}). If a \code{\link[raster]{raster}} is supplied, a query is only sent for non NA cells, which can be more efficient (since areas for which data are not required can be masked prior to function implementation). The coordinate reference system must be latitude/longitude.
#' @param db A character string that defines the database to be queried. Any option supported by Open Topo Data can be inputted, including ASTER (\code{"aster30m"}), ETOPO1 (\code{"etopo1"}), EU-DEM (\code{"eudem25m"}), Mapzen (\code{"mapzen"}), NED (\code{"ned10m"}), NZ DEM (\code{"nzdem8m"}), SRTM (\code{"srtm90m"}), EMOD bathymetry (\code{"emod2018"}) and GEBCO bathymetry (\code{"gebco2020"}) (see https://www.opentopodata.org for further details).
#' @param interpolation A character (\code{"nearest"}, \code{"bilinear"} or \code{"cubic"}) that defines the interpolation method that is used to interpolate elevation values to inputted \code{x} locations.
#' @param encoding (optional) The character encoding (e.g., \code{UTF-8}).
#' @param verbose A logical variable that defines whether or not to print messages to the console to relay function progress.
#' @param ... Additional arguments passed to \code{\link[raster]{raster}} if \code{x} is an \code{\link[raster]{Extent}} object, such as the \code{resolution}.
#'
#' @details Open Topo Data is an elevation API. Further information, including on supported datasets, supported numbers of locations (which, at the time of writing, is limited to 100) and other details are provided here: https://www.opentopodata.org/. This function requires the \code{\link[httr]{httr}} and \code{\link[jsonlite]{jsonlite}} packages to query databases.
#'
#' @return The function returns elevation (`z') values from the specified database as a matrix, if code \code{x} is a matrix, or a \code{\link[raster]{raster}}, if code \code{x} is a \code{\link[raster]{raster}} or an \code{\link[raster]{Extent}} object. Coordinates/areas without data are returned as NAs.
#'
#' @examples
#' \dontrun{
#' #### Set up example spatial data with lat/long projection
#' proj_wgs84 <- sp::CRS(SRS_string = "EPSG:4326")
#' dat_gebco_wgs84 <- raster::projectRaster(dat_gebco, crs = proj_wgs84)
#' dat_coast_wgs84 <- sp::spTransform(dat_coast, proj_wgs84)
#'
#' #### Example (1): Queries with a single set of coordinates
#' # Define coordinates
#' x <- dat_gebco_wgs84
#' x <- matrix(c(-5.616532, 56.50279), ncol = 2)
#' # Plot area
#' prettyGraphics::pretty_map(
#'   add_rasters = list(x = dat_gebco_wgs84),
#'   add_polys = list(x = dat_coast_wgs84),
#'   add_points = list(x = x),
#'   verbose = FALSE
#' )
#' # Check depth in area using available data
#' raster::extract(dat_gebco_wgs84, x)
#' # Query database
#' query_open_topo(x = x, db = "gebco2020")
#'
#' #### Example (2): Use alternative options
#' # Alternative databases, such as EMOD bathymetry
#' query_open_topo(x = x, db = "emod2018", verbose = FALSE)
#' # Set interpolation
#' query_open_topo(
#'   x = x, db = "emod2018",
#'   interpolation = "cubic", verbose = FALSE
#' )
#'
#' #### Example (2): Queries with multiple coordinates
#' # Define a random sample of coordinates
#' x <- raster::coordinates(dat_gebco_wgs84)
#' index <- sample(1:nrow(x), 25)
#' x <- x[index, ]
#' # Query database
#' depth <- query_open_topo(x = x, db = "gebco2020")
#' # Compare to manually extracted values
#' # ... (some of these are NA because dat_gebco has been masked over land)
#' depth <- cbind(depth, raster::extract(dat_gebco_wgs84, x))
#'
#' #### Example (3): Queries using an Extent object
#' # Note that only 100 locations can be queried at a time
#' # ... hence the restrictions on the resolution specified here.
#' x <- raster::extent(dat_coast_wgs84)
#' depth <- query_open_topo(x = x, nrows = 10, ncols = 10, db = "gebco2020")
#' prettyGraphics::pretty_map(
#'   add_rasters = list(x = depth),
#'   add_polys = list(x = dat_coast_wgs84),
#'   verbose = FALSE
#' )
#'
#' #### Example (4): Queries from a masked raster
#' # Focus on a small area
#' ext <- raster::extent(c(-5.709508, -5.648977, 56.48656, 56.50267))
#' x <- raster::crop(dat_gebco_wgs84, ext)
#' prettyGraphics::pretty_map(
#'   add_rasters = list(x = x),
#'   add_polys = list(x = dat_coast_wgs84),
#'   verbose = FALSE
#' )
#' # Query database
#' depth <- query_open_topo(x = x, db = "gebco2020")
#' prettyGraphics::pretty_map(
#'   add_rasters = list(x = depth),
#'   add_polys = list(x = dat_coast_wgs84),
#'   verbose = FALSE
#' )
#' }
#'
#' @seealso Open Topo Data (https://www.opentopodata.org/).
#' @author Edward Lavender
#' @export
#'

query_open_topo <- function(x,
                            db = "gebco2020",
                            interpolation = "bilinear",
                            encoding = NULL,
                            verbose = TRUE, ...) {
  #### Set up
  t_onset <- Sys.time()
  cat_to_console <- function(..., show = verbose) if (show) cat(paste(..., "\n"))
  cat_to_console(paste0("flapper::query_open_topo() called (@ ", t_onset, ")..."))

  #### Checks
  cat_to_console("... Processing 'x'...")
  # Check required packages
  if (!requireNamespace("httr", quietly = TRUE)) stop("This function requires the 'httr' package. Please install it with `install.packages('httr')` first.")
  if (!requireNamespace("jsonlite", quietly = TRUE)) stop("This function requires the 'jsonlite' package. Please install it with `install.packages('jsonlite')` first.")
  # Save inputted 'x'
  .x <- x
  # Check class of 'x' is appropriate
  if (!inherits(.x, c("matrix", "RasterLayer", "Extent"))) {
    stop("class(x) is unsupported: only 'RasterLayer', 'Extent' and 'matrix' are supported.")
  }

  #### Process coordinates
  if (!inherits(.x, "matrix")) {
    ## Get coordinates of raster
    if (inherits(.x, "Extent")) x <- raster::raster(x, ...)
    dat <- raster::coordinates(x)
    rownames(dat) <- 1:nrow(dat)

    ## If a raster is supplied, drop any NA cells
    if (inherits(.x, "RasterLayer")) {
      is_na <- is.na(x)
      id_na <- raster::Which(is_na, cells = TRUE)
      if (length(id_na) > 0) dat <- dat[-c(id_na), ]
    }
  } else {
    dat <- x
    rownames(dat) <- 1:nrow(dat)
  }

  ## Define a dataframe of cell IDs and coordinates
  dat <- data.frame(id = as.integer(rownames(dat)), x = dat[, 1], y = dat[, 2])
  dat$lat_long <- paste0(dat$y, ",", dat$x)

  #### Define query
  cat_to_console("... Setting up RESTful API request...")
  if (nrow(dat) > 100) {
    warning("More than 100 locations supplied: this may exceed max URI length.")
  }
  base <- "https://api.opentopodata.org/v1/"
  endpoint <- db
  query <- paste0("?locations=", paste0(dat$lat_long, collapse = "|"))
  query <- paste0(query, "&interpolation=", interpolation, "&nodata_value=-9999")
  request <- paste0(base, endpoint, query)

  #### Send and retrieve request
  cat_to_console("... Sending the query using HTTP...")
  response <- httr::GET(request)
  httr::warn_for_status(response)
  httr::stop_for_status(response)
  cat_to_console("... Decoding the response...")
  if (is.null(encoding)) {
    results <- jsonlite::fromJSON(httr::content(response, "text"))
  } else {
    results <- jsonlite::fromJSON(httr::content(response, "text"), encoding = encoding)
  }
  cat_to_console(paste0("... Getting status of the query: ", results$status, "."))
  results <- results$results
  dat$z <- results$elevation
  no_data <- dat$z %in% -9999
  if (any(no_data)) {
    dat$z[which(no_data)] <- NA
  }

  #### Process results
  cat_to_console("... Processing and returning data...")
  if (inherits(.x, "Extent")) {
    raster::values(x) <- dat$z
    raster::crs(x) <- sp::CRS("+init=epsg:4326")
    out <- x
  } else if (inherits(.x, "RasterLayer")) {
    x[dat$id] <- dat$z
    out <- x
  } else if (inherits(.x, "matrix")) {
    out <- dat[, c("x", "y", "z")]
  }

  #### Return results
  t_end <- Sys.time()
  duration <- difftime(t_end, t_onset, units = "mins")
  cat_to_console(paste0("... flapper::query_open_topo() call completed (@ ", t_end, ") after ~", round(duration, digits = 2), " minutes."))
  return(out)
}
