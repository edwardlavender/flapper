########################################
########################################
#### test-get_movement.R

# library(testthat)


########################################
########################################
#### get_mvt_mobility_from_acoustics()



########################################
########################################
#### get_mvt_mobility_from_archival()

#### Test errors
test_that("get_mvt_mobility_from_archival() fails as expected", {
  # Argument data does not contain all required names.
  expect_error(get_mvt_mobility_from_archival(data.frame(depth = numeric(), time = numeric())))
  expect_error(get_mvt_mobility_from_archival(data.frame(depth = numeric(), timestamp = numeric(), ID = integer()), fct = "id"))

  # Argument 'data$timestamp' is of the wrong class
  arc <- data.frame(depth = c(10, 100), timestamp = c(10, 100))
  expect_error(get_mvt_mobility_from_archival(arc))

  # Insufficient data for analysis
  arc <- data.frame(
    depth = c(10, 100),
    timestamp = as.POSIXct(c("2016-01-01", "2016-01-02"))
  )
  expect_error(get_mvt_mobility_from_archival(arc))

  # 'data' should comprise regular time step
  arc <- data.frame(
    depth = c(10, 100, 10),
    timestamp = as.POSIXct(c("2016-01-01", "2016-01-02", "2016-01-04"))
  )
  expect_error(get_mvt_mobility_from_archival(arc))
})

#### Test warnings
test_that("get_mvt_mobility_from_archival() throws warnings as expected", {
  # data$depth contains negative values: using abs(data$depth)
  # Too few observations for analysis for some individuals
  arc <- data.frame(
    depth = c(10, 100, 10, 100, 1000) * -1,
    timestamp = as.POSIXct(c("2016-01-01", "2016-01-02", "2016-01-01", "2016-01-02", "2016-01-03")),
    fct = factor(c(1, 1, 2, 2, 2))
  )
  expect_warning(get_mvt_mobility_from_archival(arc, fct = "fct"))
})

#### Test calculations
test_that("get_mvt_mobility_from_archival() calculations are correct", {
  arc <- data.frame(
    depth = c(10, 100, 5000),
    timestamp = as.POSIXct(c("2016-01-01 00:00:00", "2016-01-01 00:02:00", "2016-01-01 00:04:00"))
  )
  est <- get_mvt_mobility_from_archival(arc)
  expect_equal(est$speed_mstep, c(100 - 10, 5000 - 100))
  expect_equal(est$speed_ms, c((100 - 10) / 120, c(5000 - 100) / 120))
  est <- get_mvt_mobility_from_archival(dat_archival, fct = "individual_id")
  # Check printed statistics qualitatively
  min(est$speed_ms)
  mean(est$speed_ms)
  max(est$speed_ms)
  min(est$speed_mstep)
  mean(est$speed_mstep)
  max(est$speed_mstep)
  # Check the dataframe inputted is the same as the one returned (plus the extra columns)
  expect_equal(colnames(dat_archival), colnames(est)[!(colnames(est) %in% c("dist", "speed_ms", "speed_mstep"))])
})


########################################
########################################
#### get_mvt_resting()
