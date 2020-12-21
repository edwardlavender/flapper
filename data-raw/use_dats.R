#####################################
#####################################
#### use_dats.R

#### This code:
# 1) Adds example datasets (dats) to flapper R package.
# setwd(paste0(getwd(), "/data-raw/"))

#### Load data
dat_ids       <- readRDS("dat_ids.rds")
dat_moorings  <- readRDS("dat_moorings.rds")
dat_acoustics <- readRDS("dat_acoustics.rds")
dat_sentinel  <- readRDS("dat_sentinel.rds")
dat_archival  <- readRDS("dat_archival.rds")

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

#### Use data
usethis::use_data(dat_ids, overwrite = TRUE)
usethis::use_data(dat_moorings, overwrite = TRUE)
usethis::use_data(dat_acoustics, overwrite = TRUE)
usethis::use_data(dat_sentinel, overwrite = TRUE)
usethis::use_data(dat_archival, overwrite = TRUE)
usethis::use_data(dat_coast, overwrite = TRUE)
usethis::use_data(dat_gebco, overwrite = TRUE)
