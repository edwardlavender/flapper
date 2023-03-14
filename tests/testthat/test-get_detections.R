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
  x <- 1:100
  b0 <- 3
  b1 <- -0.5
  f <- stats::plogis
  # Check calculations are correct (output = 3L)
  p1 <- f(b0 + b1 * x)
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
      type = 1L
    )
  )
})

#### Check warnings
test_that("get_detection_days() throws expected warnings.", {
  # Warning that dots are unused
  expect_warning(get_detection_days(dat_acoustics, ... = "extra_argument"))
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
    type = 1L
  )
  expect_true(all(dat$individual_id %in% c(iid_1, iid_2)))
  expect_true(all(dat$receiver_id %in% c(rid_1, rid_2)))
  expect_true(dat$detection_days[dat$receiver_id == rid_1 & dat$individual_id == iid_1] == ans_1)
  expect_true(dat$detection_days[dat$receiver_id == rid_2 & dat$individual_id == iid_2] == ans_2)
  # Detection days between selected individual and receiver combinations
  dat <- get_detection_days(dat_acoustics,
    individual_id = c(iid_1, iid_2),
    receiver_id = c(rid_1, rid_2),
    type = 2L
  )
  expect_true(all(dat$individual_id %in% c(iid_1, iid_2)))
  expect_true(all(dat$receiver_id %in% c(rid_1, rid_2)))
  expect_true(dat$detection_days[dat$receiver_id == rid_1 & dat$individual_id == iid_1] == ans_1)
  expect_true(dat$detection_days[dat$receiver_id == rid_2 & dat$individual_id == iid_2] == ans_2)
  # Detection days for receivers/individuals not in acoustics
  expect_warning(expect_warning(
    get_detection_days(dat_acoustics,
      individual_id = c(iid_1, "blah2", "blah3"),
      receiver_id = c(rid_1, "blah2")
    )
  ))

  #### Detection days matched to another dataframe
  dat <- dat_acoustics
  dat$detection_days <- get_detection_days(dat_acoustics,
    match_to = dat,
    type = 1L
  )
  expect_true(all(dat$detection_days[dat$receiver_id == rid_1 & dat$individual_id == iid_1] == ans_1))
  expect_true(all(dat$detection_days[dat$receiver_id == rid_2 & dat$individual_id == iid_2] == ans_2))
})


########################################
########################################
#### get_detection_clumps()

test_that("get_detection_clumps() throws expected errors.", {
  # Error for missing columns
  expect_error(get_detection_clumps(data.frame(date = as.Date("2016-01-01"))))
  expect_error(
    get_detection_clumps(
      data.frame(timestamp = as.POSIXct("2016-01-01"), ID = 1),
      fct = "id"
    )
  )
})

test_that("get_detection_clumps() example calculations are correct.", {
  #### Define a hypothetical series of detections
  # ... following function examples.
  eg <-
    data.frame(
      timestamp =
        as.POSIXct(
          c(
            "2016-01-01", # one week of continuous detections
            "2016-01-02",
            "2016-01-03",
            "2016-01-04",
            "2016-01-05",
            "2016-01-06",
            "2016-01-07",
            "2016-02-01", # one day with an isolated detection
            "2016-02-03", # two days with detections
            "2016-02-04",
            "2016-02-15", # another two days with detections
            "2016-02-16",
            "2016-03-01", # five days of continuous detections
            "2016-03-02",
            "2016-03-03",
            "2016-03-04",
            "2016-03-05"
          )
        )
    )

  #### Example (1): Implement function with default options
  # Check colnames match description
  expect_true(all(colnames(get_detection_clumps(eg)) %in%
    c("n_intervals", "n_occasions", "eg_occasions")))
  # Check calculations
  expect_true(
    dplyr::all_equal(
      get_detection_clumps(eg),
      data.frame(
        n_intervals = as.integer(c(1, 2, 5, 7)),
        n_occasions = c(1, 2, 1, 1),
        eg_occasions = as.POSIXct(c(
          "2016-02-01", "2016-02-03",
          "2016-03-01", "2016-01-01"
        ))
      )
    )
  )

  #### Example (2): Implement function for multiple individuals
  eg$individual_id <- 1L
  expect_true(rlang::has_name(
    get_detection_clumps(eg, fct = "individual_id"),
    "individual_id"
  ))

  #### Example (3): Change the time interval
  ## E.g. Use an hourly interval:
  eg$timestamp <- as.POSIXct(eg$timestamp)
  expect_true(
    dplyr::all_equal(
      data.frame(
        n_intervals = 1L,
        n_occasions = 17,
        eg_occasions = as.POSIXct("2016-01-01")
      ),
      get_detection_clumps(eg, interval = "hours")
    )
  )
  ## E.g. Use a monthly interval
  get_detection_clumps(eg, interval = "months")
  expect_true(
    dplyr::all_equal(
      data.frame(
        n_intervals = 1L,
        n_occasions = 17,
        eg_occasions = as.POSIXct("2016-01-01")
      ),
      get_detection_clumps(eg, interval = "hours")
    )
  )

  #### Example (4): Identify the timing of each clump with summarise = FALSE
  expect_true(
    dplyr::all_equal(
      get_detection_clumps(eg, summarise = FALSE),
      data.frame(
        n_intervals = as.integer(c(7, 1, 2, 2, 5)),
        timestamp = as.POSIXct(c(
          "2016-01-01", "2016-02-01",
          "2016-02-03", "2016-02-15",
          "2016-03-01"
        ))
      )
    )
  )
})


########################################
########################################
#### get_detection_overlaps()


########################################
########################################
#### get_residents()

test_that("get_residents works as expected.", {
  #### Define example dataframe
  acoustics <-
    data.frame(
      timestamp = seq.POSIXt(
        as.POSIXct("2016-01-01"),
        as.POSIXct("2017-01-01"), "days"
      ),
      ID = "A",
      random = "A"
    )

  #### Check errors
  # Column names
  expect_error(get_residents(data.frame(blah = numeric())))
  # Only one value for resident_threshold_gap should be supplied
  expect_error(get_residents(acoustics, resident_threshold_gap = c(10, 1)))
  # 'resident_threshold_duration' should be a sorted, numeric vector.
  expect_error(get_residents(acoustics, resident_threshold_duration = c(10, 1)))
  # The number of resident_threshold_duration(s) and resident_label(s) should be the same.
  expect_error(get_residents(acoustics, resident_threshold_duration = 1:4))

  #### Check correct column names and labels
  # Check column names
  expect_true(all(colnames(get_residents(acoustics)) %in% c("time", "resident")))
  expect_true(rlang::has_name(get_residents(acoustics, fct = "ID"), "ID"))
  expect_true(rlang::has_name(get_residents(acoustics, keep = "random"), "random"))
  # Check labels
  get_residents(acoustics,
    resident_labels = c("STR", "LTR")
  )$resident == "LTR"



  #### Check calculations for example individual
  res <- get_residents(acoustics)
  expect_true(res$time == 366 & res$resident == "L")
  res <- get_residents(acoustics[1:10, , drop = FALSE],
    resident_threshold_duration = 365,
    resident_labels = "resident"
  )
  expect_true(res$time == 9 & res$resident == "N")

  #### Check calculations for multiple individuals
  ## Define dataframes for non residents/short-term residents and long term residents
  non <-
    data.frame(
      timestamp = as.POSIXct(c("2016-01-01", "2016-02-05")),
      id = 1
    )
  short <-
    data.frame(
      timestamp = seq.POSIXt(
        as.POSIXct("2016-01-01"),
        as.POSIXct("2016-05-01"), "hours"
      ),
      id = 2
    )
  long <-
    data.frame(
      timestamp = seq.POSIXt(
        as.POSIXct("2016-01-01"),
        as.POSIXct("2018-01-01"),
        "weeks"
      ),
      id = 3
    )

  ## Check calculations with default options
  # Check residency categories
  res <- get_residents(rbind(non, short, long), fct = "id")
  expect_true(all(res$resident == c("N", "S", "L")))
  # Check residency categories if individual need to be detected every day
  # ... The individual that is detected weekly is now no longer resident
  res <- get_residents(rbind(non, short, long),
    fct = "id",
    resident_threshold_gap = 1
  )
  expect_true(all(res$resident == c("N", "S", "N")))
})
