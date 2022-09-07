########################################
########################################
#### test-dists.R

# library(testthat)


########################################
########################################
#### raster:::..planedist2()

# .planedist2 <- flapper:::.planedist2
test_that(".planedist2 returns correct answers", {
  expect_equal(
    .planedist2(matrix(c(1, 2), ncol = 2), matrix(c(3, 4), ncol = 2)),
    matrix(sqrt((3 - 1)^2 + (4 - 2)^2))
  )
})


########################################
########################################
#### dist_btw_receivers()

test_that("dist_btw_receivers() fails without required columns.", {
  # Test error with missing names
  expect_error(dist_btw_receivers(data.frame(receiver_id = integer(),
                                             receiver_long = numeric())))
  # Test error with mismatched names
  expect_error(dist_btw_receivers(data.frame(receiver_id = integer(),
                                             receiver_long = numeric(),
                                             receiver_northing = numeric())))
})


test_that("dist_btw_receivers() returns correct distances.", {

  #### Test distance calculations with a planar landscape/coordinates
  # Define planar landscape with known properties
  proj_utm <- sp::CRS(SRS_string = "EPSG:32629")
  r <- raster::raster(nrows = 3, ncols = 3,
                      crs = proj_utm,
                      resolution = c(5, 5),
                      ext = raster::extent(0, 15, 0, 15))
  r[] <- 1:raster::ncell(r)
  # Define receivers
  rxy <- raster::coordinates(r)
  rxy <- rxy[c(1, 8, 9), ]
  raster::plot(r)
  points(rxy)
  # Check distances
  answer <- raster::pointDistance(rxy, allpairs = TRUE, lonlat = FALSE) %>% as.vector()
  expect_equal(
    answer,
    dist_btw_receivers(data.frame(receiver_id = seq_len(nrow(rxy)),
                                  receiver_easting = rxy[, 1],
                                  receiver_northing = rxy[, 2]),
                       f = function(x) x*1000)$dist
  )

  #### Test distance calculations with a lonlat landscape/coordinates
  # Define lon/lat landscape with known properties
  proj_wgs84 <- sp::CRS(SRS_string = "EPSG:4326")
  r <- raster::projectRaster(dat_gebco, crs = proj_wgs84)
  r[] <- 1:raster::ncell(r)
  # Define receivers
  rxy <- raster::coordinates(r)
  rxy <- rxy[sample(seq_len(nrow(rxy)), 5), ]
  raster::plot(r)
  points(rxy)
  # Check distances
  answer <- raster::pointDistance(rxy, allpairs = TRUE, lonlat = TRUE)
  answer[upper.tri(answer)] <- t(answer)[upper.tri(answer)]
  answer <- as.vector(answer)
  expect_equal(
    answer,
    dist_btw_receivers(data.frame(receiver_id = seq_len(nrow(rxy)),
                                  receiver_long = rxy[, 1],
                                  receiver_lat = rxy[, 2]),
                       f = function(x) x*1000)$dist
  )

  #### Test distance calculations with return = "matrix"
  expect_equal(
    answer,
    dist_btw_receivers(data.frame(receiver_id = seq_len(nrow(rxy)),
                                  receiver_long = rxy[, 1],
                                  receiver_lat = rxy[, 2]),
                       f = function(x) x*1000,
                       return = "matrix") %>% as.vector()
  )
})


########################################
########################################
#### dist_btw_clicks()



########################################
########################################
#### dist_btw_points_3d()

test_that("dist_btw_points_3d() returns correct answers", {

  #### Expect distance 0
  expect_equal(dist_btw_points_3d(1, 1, 1, 1, 1, 1), 0)

  #### Expect NA
  expect_equal(dist_btw_points_3d(1, 1, 1, NA, 1, 1), NA_integer_)

  #### Expect correct solution with random points
  coords <- runif(6, 0, 1)
  print(coords)
  expect_equal(
    dist_btw_points_3d(coords[1], coords[2],
                       coords[3], coords[4],
                       coords[5], coords[6]),
    sqrt((coords[2] - coords[1])^2 + (coords[4] - coords[3])^2 + (coords[6] - coords[5])^2)
    )

  #### Expect correct solutions with multiple points
  # Define a dataframe with coordinates (e.g., animal movement path)
  xy <- data.frame(x = c(1, 2, 3, 4),
                   y = c(1, 2, 3, 4),
                   z = c(1, 2, 3, 4))
  # Calculate distances manually
  dist_0 <- rep(NA, 4)
  for(i in seq_len(nrow(xy))){
    j <- i+1
    dist_0[i] <- sqrt((xy$x[j] - xy$x[i])^2 + (xy$y[j] - xy$y[i])^2 + (xy$z[j] - xy$z[i])^2)
  }
  # Calculate distances using function
  dist_1 <- dist_btw_points_3d(xy$x, dplyr::lead(xy$x),
                               xy$y, dplyr::lead(xy$y),
                               xy$z, dplyr::lead(xy$z))
  # Show that distances are equal
  expect_equal(dist_0, dist_1)
})


########################################
########################################
#### dist_over_surface()

test_that("dist_over_surface() returns correct answers", {

  #### Simulate a hypothetical landscape (as in Examples)
  # Define a miniature, blank landscape with known dimensions
  proj_utm <- sp::CRS(SRS_string = "EPSG:32629")
  r <- raster::raster(nrows = 3, ncols = 3,
                      crs = proj_utm,
                      resolution = c(5, 5),
                      ext = raster::extent(0, 15, 0, 15))
  # Define a matrix of hypothetical values for the landscape
  mat <- matrix(c(5, 10, 3,
                  2, 1, 4,
                  5, 6, 6), ncol = 3, nrow = 3, byrow = TRUE)
  r[] <- mat
  # Visualise simulated landscape
  # raster::plot(r)
  # raster::text(r)

  #### Check the distance between two example adjacent points
  # Define points and distance between them
  path_cells   <- c(1, 2)
  path_matrix  <- sp::coordinates(r)[path_cells, ]
  answer       <- sqrt(5^2 + (10-5)^2)
  # Expect error
  # ... Argument 'path' must be of class(es) 'matrix'...
  expect_error(dist_over_surface(c(1, 2), r))
  # Expect message
  # ... The CRS of the surface is not NA. Note that this function only supports planar surfaces.
  expect_message(dist_over_surface(path_matrix, r))
  # Check solutions using point matrix and spatial lines object
  expect_equal(dist_over_surface(path_matrix, r),
               answer)
  expect_equal(dist_over_surface(Orcs::coords2Lines(path_matrix, ID = 1), r),
               answer)
  # Expect warning
  # ... The CRS of the 'surface' and the 'path' are not identical.
  expect_warning(
    dist_over_surface(Orcs::coords2Lines(path_matrix, ID = 1, proj4string = sp::CRS(NA_character_)),
                      r))

  #### Check the distance between two example diagonal points
  # Define points and distance between them
  path_cells  <- c(1, 5)
  path_matrix <- sp::coordinates(r)[path_cells, ]
  answer      <- sqrt(sqrt(5^2 + 5^2)^2 + (5 - 1)^2)
  expect_equal(dist_over_surface(path_matrix, r),
               answer)
  expect_equal(dist_over_surface(Orcs::coords2Lines(path_matrix, ID = 1), r),
               answer)

  #### Check the distance along a longer path
  # Define points and distance between them
  path_cells  <- c(1, 2, 3)
  path_matrix <- sp::coordinates(r)[path_cells, ]
  answer      <- sqrt(5^2 + (10-5)^2) + sqrt(5^2 + (10-3)^2)
  # Check solutions using point matrix and spatial lines object
  expect_equal(dist_over_surface(path_matrix, r),
               answer)
  expect_equal(dist_over_surface(Orcs::coords2Lines(path_matrix, ID = 1), r),
               answer)
})

