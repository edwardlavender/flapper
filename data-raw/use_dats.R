#####################################
#####################################
#### use_dats.R

#### This code:
# 1) Adds example datasets (dats) to flapper R package.


#####################################
#####################################
#### Set up

#### Wipe workspace
rm(list = ls())
try(pacman::p_unload("all"), silent = TRUE)

#### Essential packages
library(dplyr)

#### Define helper functions
## Define file paths
here_data_raw <- function(...) here::here("data-raw", ...)
## Drop non-non-ASCII characters in PROJ comments in spatial data
# See here to check for non-ASCII: https://pages.cs.wisc.edu/~markm/ascii.html
# See also issue here: https://github.com/r-spatial/sf/issues/1341
drop_non_ascii <- function(x) {
  x |>
    stringi::stri_trans_general("latin-ascii") |>
    stringr::str_replace_all("°|º", "degrees")
}


#####################################
#####################################
#### Movement data

#### Load data
dat_ids <- readRDS(here_data_raw("dat_ids.rds"))
dat_moorings <- readRDS(here_data_raw("dat_moorings.rds"))
dat_acoustics <- readRDS(here_data_raw("dat_acoustics.rds"))
dat_sentinel <- readRDS(here_data_raw("dat_sentinel.rds"))
dat_archival <- readRDS(here_data_raw("dat_archival.rds"))

#### Processing
## dat_acoustics
# Add individual-specific detection index to dat_acoustics
dat_acoustics <-
  dat_acoustics %>%
  dplyr::group_by(.data$individual_id) %>%
  dplyr::mutate(index = dplyr::row_number()) %>%
  dplyr::select(
    individual_id, transmitter_id,
    index, timestamp,
    receiver_id, receiver,
    receiver_long, receiver_lat, receiver_depth
  ) %>%
  data.frame()
# Visualise detection indices
lapply(
  split(dat_acoustics, dat_acoustics$individual_id),
  function(d) plot(d$index, type = "l")
)
# Round detection time series to the nearest two minutes and drop duplicate observations
# ... This is beneficial in the AC* algorithms.
dat_acoustics <-
  dat_acoustics %>%
  dplyr::group_by(.data$individual_id, .data$receiver_id) %>%
  dplyr::mutate(timestamp = lubridate::round_date(.data$timestamp, unit = "2 mins")) %>%
  dplyr::filter(!duplicated(.data$timestamp)) %>%
  data.frame()
nrow(dat_acoustics)
## moorings
rownames(dat_moorings) <- dat_moorings$receiver_id


#####################################
#####################################
#### Spatial data

#### Obtain spatial data
## Open Source coastline data
dat_coast_fol <- readRDS(url("https://biogeo.ucdavis.edu/data/gadm3.6/Rsp/gadm36_GBR_0_sp.rds"))
## Open source bathymetry data
# Define a buffered region around the receivers for which to obtain bathymetry data
proj_wgs84 <- sp::CRS(SRS_string = "EPSG:4326")
attr(proj_wgs84, "comment") <- drop_non_ascii(attr(proj_wgs84, "comment"))
proj_utm <- sp::CRS(SRS_string = "EPSG:32629")
attr(proj_utm, "comment") <- drop_non_ascii(attr(proj_utm, "comment"))
rxy_wgs84 <- sp::SpatialPoints(dat_moorings[, c("receiver_long", "receiver_lat")], proj_wgs84)
rxy_utm <- sp::spTransform(rxy_wgs84, proj_utm)
rxy_utm_buf <- rgeos::gBuffer(rxy_utm, width = 1000)
bounds_utm <- raster::extent(rxy_utm_buf)
bounds_utm <- sp::SpatialPoints(sp::coordinates(bounds_utm), proj_utm)
bounds_wgs84 <- sp::spTransform(bounds_utm, proj_wgs84)
# Examine coordinates
# -5.786025, -5.562533, 56.34059, 56.53355
# Use these coordinates to download manually Open source GEBCO bathymetry data from https://download.gebco.net
# Data source: GEBCO Compilation Group (2019) GEBCO 2019 Grid (doi:10.5285/836f016a-33be-6ddc-e053-6c86abc0788e)
# Data saved in /data_raw/
dat_gebco_fol <-
  raster::raster(here_data_raw("gebco_2020_n56.53355_s56.34059_w-5.786025_e-5.562533.tif"))

#### Process spatial data
# Crop coastline to boundaries
area <- raster::extent(bounds_wgs84)
dat_coast_fol <- raster::crop(dat_coast_fol, area)
# Use spatial data with UTM coordinates
dat_coast <- sp::spTransform(dat_coast_fol, proj_utm)
dat_gebco <- raster::projectRaster(dat_gebco_fol, crs = proj_utm)
# Check raster validity (no missing slots)
str(dat_gebco)
# Process bathymetry data to remove observations on land
dat_gebco[dat_gebco[] >= 0] <- NA
# Use absolute values
dat_gebco <- abs(dat_gebco)
# Visual checks
raster::plot(dat_coast)
raster::plot(dat_gebco, add = TRUE)
raster::lines(dat_coast)
points(rxy_utm, cex = 2)
axis(side = 1)
axis(side = 2)


#####################################
#####################################
#### Example function outputs

#### Comments
# In this section, some of the outputs of functions in the flapper package
# ... are created and included as examples in the flapper package
# ... to help illustrate examples quickly.
library(flapper)


#####################################
#### dat_dcpf

#### Define data [example 1 in package]

## Sample species
# In this example, we consider flapper skate (Dipturus intermedius)
# ... off the west coast of Scotland.

## Define starting location (optional) in UTM coordinates
xy <- matrix(c(708886.3, 6254404), ncol = 2)
## Define 'observed' depth time series using absolute values
# Imagine these are observations made every two minutes
depth <- c(
  163.06, 159.71, 153.49, 147.04, 139.86, 127.19, 114.75,
  99.44, 87.01, 78.16, 70.03, 60.23, 49.96, 35.39,
  27.75, 20.13, 12.73, 11.32
)
depth <- data.frame(depth = depth)

## Define surface over which movement must occur
# We will use the example dat_gebco bathymetry dataset
# This is relatively coarse in resolution, so we need to re-sample
# ... the raster to generate a finer-resolution raster such that
# ... our animal (e.g., a skate) can transition between cells
# ... in the two-minutes between depth observations. We'll still restrict
# ... raster resolution to minimise package dataset sizes though. For speed,
# ... we will focus on a small area around the origin. We could
# ... also process the raster in other ways (e.g., mask any areas of land)
# ... to improve efficiency.
surface <- dat_gebco
boundaries <- raster::extent(707884.6, 709884.6, 6253404, 6255404)
blank <- raster::raster(boundaries, res = c(25, 25))
surface <- raster::resample(surface, blank)

## Define movement model
# The default movement model is suitable, with skate moving typically
# ... less than 200 m in a two-minute period.

## Visualise movement surface, with starting location overlaid
prettyGraphics::pretty_map(
  add_rasters = list(x = surface),
  add_points = list(x = xy),
  verbose = FALSE
)

#### (A) Implement DC algorithm
# Use save_args = TRUE for package examples
dat_dc <- dc(
  archival = depth,
  bathy = surface,
  calc_depth_error = function(...) matrix(c(-30, 30), nrow = 2),
  normalise = TRUE,
  save_record_spatial = 0L,
  save_args = TRUE
)

#### Example (1): Implement algorithm using default options
# Note that because the bathymetry data is very coarse, we have to
# ... force the depth_error to be high in this example.
dc_1 <- dc(
  archival = depth,
  bathy = surface,
  calc_depth_error = function(...) matrix(c(-30, 30), nrow = 2),
  save_record_spatial = NULL
)
dcpf_1 <- pf(
  record = lapply(dc_1$archive[[1]]$record, function(elm) elm$spatial[[1]]$map_timestep),
  origin = xy,
  calc_distance = "euclid",
  calc_distance_euclid_fast = FALSE,
  bathy = surface,
  calc_movement_pr =
    function(distance, ...) {
      pr <- stats::plogis(5 + distance * -0.05)
      pr[distance > 200] <- 0
      return(pr)
    },
  mobility = 200,
  n = 10L,
  seed = 5
)

# The function returns a list with particle histories and arguments
summary(dcpf_1)
utils::str(dcpf_1)
dat_dcpf_histories <- dcpf_1
dat_dcpf_paths <- pf_simplify(dat_dcpf_histories)
saveRDS(dat_dcpf_histories, paste0(tempdir(), "/dat_dcpf_histories.rds"))
file.size(paste0(tempdir(), "/dat_dcpf_histories.rds")) / 1e6


#####################################
#### dat_acdc_contours:
# ... an example dataset of detection containers required for ACDC algorithm
# ... useful for demonstrating the acdc algorithm

## Define data for setup_acdc()
# Define coordinates of receivers as SpatialPoints with UTM CRS
# CRS of receiver locations as recorded in dat_moorings
proj_wgs84 <- sp::CRS(SRS_string = "EPSG:4326")
# CRS of receiver locations required
proj_utm <- sp::CRS(SRS_string = "EPSG:32629")
# Define SpatialPointsDataFrame object
xy_wgs84 <- sp::SpatialPoints(dat_moorings[, c("receiver_long", "receiver_lat")], proj_wgs84)
xy_utm <- sp::spTransform(xy_wgs84, proj_utm)
xy_utm <-
  sp::SpatialPointsDataFrame(
    xy_utm,
    dat_moorings[, "receiver_id", drop = FALSE]
  )

## Define a list of containers with specified parameters
# Use reduced 'quadsegs' param to minimise file size
dat_containers <- acs_setup_containers(
  xy = xy_utm,
  detection_range = 425,
  coastline = dat_coast,
  boundaries = raster::extent(dat_gebco),
  plot = TRUE,
  resolution = 5,
  verbose = TRUE
)
# Check size (Mb) of file prior to inclusion in package
saveRDS(dat_containers, paste0(tempdir(), "/dat_containers.rds"))
file.size(paste0(tempdir(), "/dat_containers.rds")) / 1e6


#####################################
#### dat_acdc(): an example output from the .acs() function
# ... useful for demonstrating plotting functions

## Prepare movement time series
# Focus on an example individual
id <- 25
acc <- dat_acoustics[dat_acoustics$individual_id == id, ]
arc <- dat_archival[dat_archival$individual_id == id, ]
# Focus on a sample of data to keep package contents small
acc <- acc[acc$timestamp <= as.POSIXct("2016-03-19"), ]
arc <- arc[arc$timestamp <= as.POSIXct("2016-03-19"), ]
# Focus on the subset of data for which we have both acoustic and archival detections
acc <- acc[acc$timestamp >= min(arc$timestamp) - 2 * 60 &
  acc$timestamp <= max(arc$timestamp) + 2 * 60, ]
arc <- arc[arc$timestamp >= min(acc$timestamp) - 2 * 60 &
  arc$timestamp <= max(acc$timestamp) + 2 * 60, ]
range(arc$timestamp)
range(acc$timestamp)
# Examine data in this region
pp <- par(oma = c(3, 3, 3, 3), xpd = NA)
prettyGraphics::pretty_plot(arc$timestamp, arc$depth * -1, type = "l")
prettyGraphics::pretty_line(acc$timestamp, add = TRUE)
par(pp)

#### dat_acdc
# Use save_args = FALSE to minimise space requirements
dat_acdc <- acdc(
  acoustics = acc,
  archival = arc,
  bathy = dat_gebco,
  detection_containers = dat_containers,
  mobility = 200,
  calc_depth_error = function(...) matrix(c(-2.5, 2.5), nrow = 2),
  normalise = TRUE,
  save_record_spatial = 1:3,
  progress = 3L,
  verbose = TRUE,
  con = "",
  split = "2 hours",
  save_args = FALSE
)

# Check size of file prior to inclusion in package
saveRDS(dat_acdc, paste0(tempdir(), "/dat_acdc.rds"))
file.size(paste0(tempdir(), "/dat_acdc.rds")) / 1e6


#####################################
#####################################
#### Global flapper options

flapper_run_parallel <- TRUE # TRUE
flapper_run_slow <- TRUE # TRUE


#####################################
#####################################
#### Add data to pkg

#### Use data
# movement time series
usethis::use_data(dat_ids, overwrite = TRUE)
usethis::use_data(dat_moorings, overwrite = TRUE)
usethis::use_data(dat_acoustics, overwrite = TRUE)
usethis::use_data(dat_sentinel, overwrite = TRUE)
usethis::use_data(dat_archival, overwrite = TRUE)
# spatial data
usethis::use_data(dat_coast, overwrite = TRUE)
usethis::use_data(dat_gebco, overwrite = TRUE)
# function example data
usethis::use_data(dat_containers, overwrite = TRUE)
usethis::use_data(dat_dc, overwrite = TRUE)
usethis::use_data(dat_acdc, overwrite = TRUE)
usethis::use_data(dat_dcpf_histories, overwrite = TRUE)
usethis::use_data(dat_dcpf_paths, overwrite = TRUE)
# global options
usethis::use_data(flapper_run_parallel, overwrite = TRUE)
usethis::use_data(flapper_run_slow, overwrite = TRUE)
# check size of data directory
sum(sapply(list.files(here::here("data"), full.names = TRUE), file.size)) / 1e6


#### End of code.
#####################################
#####################################
