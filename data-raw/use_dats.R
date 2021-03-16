#####################################
#####################################
#### use_dats.R

#### This code:
# 1) Adds example datasets (dats) to flapper R package.
# setwd(paste0(getwd(), "/data-raw/"))


#####################################
#####################################
#### Movement data

#### Load data
dat_ids       <- readRDS("dat_ids.rds")
dat_moorings  <- readRDS("dat_moorings.rds")
dat_acoustics <- readRDS("dat_acoustics.rds")
dat_sentinel  <- readRDS("dat_sentinel.rds")
dat_archival  <- readRDS("dat_archival.rds")


#####################################
#####################################
#### Spatial data

#### Obtain spatial data
## Open Source coastline data
dat_coast_fol <- readRDS(url("https://biogeo.ucdavis.edu/data/gadm3.6/Rsp/gadm36_GBR_0_sp.rds"))
## Open source bathymetry data
# Define a buffered region around the receivers for which to obtain bathymetry data
proj_wgs84  <- sp::CRS("+init=epsg:4326")
proj_utm    <- sp::CRS(paste("+proj=utm +zone=29 +datum=WGS84",
                              "+units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))
rxy_wgs84    <- sp::SpatialPoints(dat_moorings[, c("receiver_long", "receiver_lat")], proj_wgs84)
rxy_utm      <- sp::spTransform(rxy_wgs84, proj_utm)
rxy_utm_buf  <- rgeos::gBuffer(rxy_utm, width = 1000)
bounds_utm   <- raster::extent(rxy_utm_buf)
bounds_utm   <- sp::SpatialPoints(sp::coordinates(bounds_utm), proj_utm)
bounds_wgs84 <- sp::spTransform(bounds_utm, proj_wgs84)
# Examine coordinates
# -5.786025, -5.562533, 56.34059, 56.53355
# Use these coordinates to download manually Open source GEBCO bathymetry data from https://download.gebco.net
# Data source: GEBCO Compilation Group (2019) GEBCO 2019 Grid (doi:10.5285/836f016a-33be-6ddc-e053-6c86abc0788e)
# Data saved in /data_raw/
dat_gebco_fol <- raster::raster("./data-raw/gebco_2020_n56.53355_s56.34059_w-5.786025_e-5.562533.tif")

#### Process spatial data
# Crop coastline to boundaries
area <- raster::extent(bounds_wgs84)
dat_coast_fol <- raster::crop(dat_coast_fol, area)
# Use spatial data with UTM coordinates
dat_coast <- sp::spTransform(dat_coast_fol, proj_utm)
dat_gebco <- raster::projectRaster(dat_gebco_fol, crs = proj_utm)
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
depth <- c(163.06, 159.71, 153.49, 147.04, 139.86, 127.19, 114.75,
           99.44,  87.01,  78.16,  70.03,  60.23,  49.96,  35.39,
           27.75,  20.13,  12.73,  11.32)

## Define surface over which movement must occur
# We will use the example dat_gebco bathymetry dataset
# This is relatively coarse in resolution, so we need to re-sample
# ... the raster to generate a finer-resolution raster such that
# ... our animal (e.g., a skate) can transition between cells
# ... in the two-minutes between depth observations. For speed,
# ... we will focus on a small area around the origin. We could
# ... also process the raster in other ways (e.g., mask any areas of land)
# ... to improve efficiency.
surface    <- dat_gebco
boundaries <- raster::extent(707884.6, 709884.6, 6253404, 6255404)
blank      <- raster::raster(boundaries, res = c(5, 5))
surface    <- raster::resample(surface, blank)

## Define movement model
# The default movement model is suitable, with skate moving typically
# ... less than 200 m in a two-minute period.

## Visualise movement surface, with starting location overlaid
prettyGraphics::pretty_map(add_rasters = list(x = surface),
                           add_points = list(x = xy),
                           verbose = FALSE)

#### Example (1): Implement algorithm using default options
# Note that because the bathymetry data is very coarse, we have to
# ... force the depth_error to be high in this example.
paths <- dcpf(archival = depth,
              bathy = surface,
              origin  = xy,
              depth_error = 30,
              calc_distance = "euclid",
              calc_movement_pr =
                function(distance) {
                  pr <- stats::plogis(10 + distance * -0.05)
                  pr[distance > 500] <- 0
                  return(pr)
                },
              n = 10L,
              seed = 1)
# The function returns a dataframe with information for each path
utils::str(paths)
dat_dcpf <- paths


#####################################
#### dat_acdc_contours:
# ... an example dataset of acoustic centroids required for ACDC algorithm
# ... useful for demonstrating the acdc algorithm

## Define data for setup_acdc()
# Define coordinates of receivers as SpatialPoints with UTM CRS
# CRS of receiver locations as recorded in dat_moorings
proj_wgs84 <- sp::CRS("+init=epsg:4326")
# CRS of receiver locations required
proj_utm <- sp::CRS(paste("+proj=utm +zone=29 +datum=WGS84",
                          "+units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))
# Define SpatialPoints object
xy_wgs84 <- sp::SpatialPoints(dat_moorings[, c("receiver_long", "receiver_lat")], proj_wgs84)
xy_utm <- sp::spTransform(xy_wgs84, proj_utm)

## Define a list of centroids with specified parameters
dat_centroids <- acdc_setup(rs = dat_moorings$receiver_id,
                            xy = xy_utm,
                            detection_range = 425,
                            mobility = 200,
                            n_timesteps = 25,
                            coastline = dat_coast,
                            boundaries = raster::extent(dat_gebco),
                            plot = TRUE,
                            cl = parallel::makeCluster(11L),
                            verbose = TRUE
)
# Check size (Mb) of file prior to inclusion in package
saveRDS(dat_centroids, paste0(tempdir(), "/dat_centroids.rds"))
file.size(paste0(tempdir(), "/dat_centroids.rds"))/1e6


#####################################
#### dat_acdc(): an example output from the .acdc() function
# ... useful for demonstrating plotting functions

## Prepare movement time series
# Add required columns to dataframes:
dat_acoustics$timestamp_num <- as.numeric(dat_acoustics$timestamp)
dat_archival$timestamp_num  <- as.numeric(dat_archival$timestamp)
# Focus on an example individual
id <- 25
acc <- dat_acoustics[dat_acoustics$individual_id == id, ]
arc <- dat_archival[dat_archival$individual_id == id, ]
# Focus on the subset of data for which we have both acoustic and archival detections
acc <- acc[acc$timestamp >= min(arc$timestamp) - 2*60 &
             acc$timestamp <= max(arc$timestamp) + 2*60, ]
arc <- arc[arc$timestamp >= min(acc$timestamp) - 2*60 &
             arc$timestamp <= max(acc$timestamp) + 2*60, ]
# Focus on a sample of data to keep package contents small
acc <- acc[acc$timestamp <= as.POSIXct("2016-03-19"), ]
arc <- arc[arc$timestamp <= as.POSIXct("2016-03-19"), ]
range(arc$timestamp)
range(acc$timestamp)
# Examine data in this region
prettyGraphics::pretty_plot(arc$timestamp, arc$depth*-1, type = "l")
prettyGraphics::pretty_line(acc$timestamp, add = TRUE)

#### dat_acdc
dat_acdc <- acdc(acoustics = acc,
                 archival = arc,
                 bathy = dat_gebco,
                 space_use = NULL,
                 detection_range = 425,
                 mobility = 200,
                 depth_error = 2.5,
                 acc_centroids = dat_centroids,
                 plot = 1:50,
                 png_param = list(),
                 progress = 3L,
                 verbose = TRUE,
                 con = "",
                 split = "2 hours"
)

# Check size of file prior to inclusion in package
saveRDS(dat_acdc, paste0(tempdir(), "/dat_acdc.rds"))
file.size(paste0(tempdir(), "/dat_acdc.rds"))/1e6


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
usethis::use_data(dat_dcpf, overwrite = TRUE)
usethis::use_data(dat_centroids, overwrite = TRUE)
usethis::use_data(dat_acdc, overwrite = TRUE)

#### End of code.
#####################################
#####################################
