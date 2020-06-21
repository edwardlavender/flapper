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

#### Use data
usethis::use_data(dat_ids, overwrite = TRUE)
usethis::use_data(dat_moorings, overwrite = TRUE)
usethis::use_data(dat_acoustics, overwrite = TRUE)
usethis::use_data(dat_sentinel, overwrite = TRUE)
usethis::use_data(dat_archival, overwrite = TRUE)

