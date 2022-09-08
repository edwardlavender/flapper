########################################
########################################
#### test-get_detections.R

# library(testthat)

########################################
########################################
#### get_detection_pr()

test_that("get_detection_pr() performs properly.", {
  #### Check errors
  # Coefficients should be length one
  expect_error(get_detection_pr(beta_0 = 1:2))

  #### Check warnings
  # The output should be supported
  expect_warning(get_detection_pr(output = 0))

  #### Check calculations
  # Define function to check calculated values are correct at every distance
  all_equal_vec <- Vectorize(all.equal)
  # Define example data
  x  <- 1:100
  b0 <- 3
  b1 <- -0.5
  f <- gtools::inv.logit
  # Check calculations are correct (output = 3L)
  p1 <- gtools::inv.logit(b0 + b1 * x)
  p2 <- get_detection_pr(distance = x, beta_0 = b0, beta_1 = b1, inv_link = f)
  expect_true(all(all_equal_vec(p1, p2)))
  # Check calculations are correct (output = 2L)
  p3 <- get_detection_pr(distance = x, beta_0 = b0, beta_1 = b1, inv_link = f, output = 2L)
  all(all_equal_vec(p1, p3))

  #### Check outputs
  # Check only plot produced with output = 1L
  expect_type(p2, "double")
  expect_equal(names(attributes(p2)), c("X", "beta", "inv_link"))
  expect_invisible(get_detection_pr(output = 1L))
})


########################################
########################################
#### get_detection_containers()


########################################
########################################
#### get_detection_containers_overlap()


########################################
########################################
#### get_detection_area_sum()


########################################
########################################
#### get_detection_area_ts()


########################################
########################################
#### get_n_operational_ts()


########################################
########################################
#### get_id_rec_overlap()


########################################
########################################
#### get_detection_containers_envir()


########################################
########################################
#### get_detection_days()

#### Check errors
test_that("get_detection_days() throws expected errors.", {
  # Error for missing column names
  expect_error(get_detection_days(data.frame(a = 1)))
  expect_error(get_detection_days(dat_acoustics, match_to = data.frame(a = 1)))
  # Errors for no data for specified individuals/receivers
  expect_error(expect_warning(get_detection_days(dat_acoustics, individual_id = "blah")))
  expect_error(expect_warning(get_detection_days(dat_acoustics, receiver_id = "blah")))
  # Error for different numbers of individual_id and receiver_id with type = 1L
  expect_error(
    get_detection_days(dat_acoustics,
                       individual_id = unique(dat_acoustics$individual)[1],
                       receiver_id = unique(dat_acoustics$receiver_id)[1:2],
                       type = 1L)
  )
})

#### Check warnings
test_that("get_detection_days() throws expected warnings.", {
  # Warning that dots are unused
  expect_warning(get_detection_days(dat_acoustics, ... =  "extra_argument"))
})

#### Check calculations
test_that("get_detection_days() calculations are correct.", {

  #### Define data for calculation checks
  # Define dates
  dat_acoustics$date <- as.Date(dat_acoustics$timestamp)
  # Define example receivers and individuals
  rid_1 <- 3
  iid_1 <- 25
  rid_2 <- 24
  iid_2 <- 28
  # Get detection days for rid_1 and iid_1 and rid_2 and iid_2
  ans_1 <- length(unique(dat_acoustics$date[dat_acoustics$receiver_id == rid_1 & dat_acoustics$individual_id == iid_1]))
  ans_2 <- length(unique(dat_acoustics$date[dat_acoustics$receiver_id == rid_2 & dat_acoustics$individual_id == iid_2]))

  #### Check calculations
  # Detection days between all detected individuals and receivers
  dat <- get_detection_days(dat_acoustics)
  dat$detection_days[dat$receiver_id == rid_1 & dat$individual_id == iid_1] == ans_1
  dat$detection_days[dat$receiver_id == rid_2 & dat$individual_id == iid_2] == ans_2
  # Detection days between selected individuals and receiver pairs
  dat <- get_detection_days(dat_acoustics,
                            individual_id = c(iid_1, iid_2),
                            receiver_id = c(rid_1, rid_2),
                            type = 1L)
  expect_true(all(dat$individual_id %in% c(iid_1, iid_2)))
  expect_true(all(dat$receiver_id %in% c(rid_1, rid_2)))
  expect_true(dat$detection_days[dat$receiver_id == rid_1 & dat$individual_id == iid_1] == ans_1)
  expect_true(dat$detection_days[dat$receiver_id == rid_2 & dat$individual_id == iid_2] == ans_2)
  # Detection days between selected individual and receiver combinations
  dat <- get_detection_days(dat_acoustics,
                            individual_id =  c(iid_1, iid_2),
                            receiver_id = c(rid_1, rid_2),
                            type = 2L)
  expect_true(all(dat$individual_id %in% c(iid_1, iid_2)))
  expect_true(all(dat$receiver_id %in% c(rid_1, rid_2)))
  expect_true(dat$detection_days[dat$receiver_id == rid_1 & dat$individual_id == iid_1] == ans_1)
  expect_true(dat$detection_days[dat$receiver_id == rid_2 & dat$individual_id == iid_2] == ans_2)
  # Detection days for receivers/individuals not in acoustics
  expect_warning(expect_warning(
    get_detection_days(dat_acoustics,
                       individual_id = c(iid_1, "blah2", "blah3"),
                       receiver_id = c(rid_1, "blah2")))

  )

  #### Detection days matched to another dataframe
  dat <- dat_acoustics
  dat$detection_days <- get_detection_days(dat_acoustics,
                                           match_to = dat,
                                           type = 1L)
  expect_true(all(dat$detection_days[dat$receiver_id == rid_1 & dat$individual_id == iid_1] == ans_1))
  expect_true(all(dat$detection_days[dat$receiver_id == rid_2 & dat$individual_id == iid_2] == ans_2))

})


########################################
########################################
#### get_detection_clumps()


########################################
########################################
#### get_detection_overlaps()


########################################
########################################
#### get_residents()
