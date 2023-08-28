######################################
######################################
#### .acs_pl()

#' @title Intermediate wrapper for \code{\link[flapper]{.acs}} that supports parallelisation
#' @description This function implements the acoustic-container (AC) and acoustic-container depth-contour (ACDC) algorithms. This is called via a front-end function (i.e. \code{\link[flapper]{ac}} or \code{\link[flapper]{acdc}}). It checks and processes inputs and implements the selected algorithm via calls to \code{\link[flapper]{.acs}}. Outputs are returned in a named list.
#'
#' @param acoustics A dataframe, or a list of dataframes, that contains passive acoustic telemetry detection time series (see \code{\link[flapper]{dat_acoustics}} for an example) for a single individual. Each dataframe should contain the following columns: an integer vector of receiver IDs, named `receiver_id'; an integer vector of detection indices, named `index'; and a POSIXct vector of time stamps when detections were made, named `timestamp'. If a list of dataframes is supplied, dataframes must be refer to the detections of a single individual and be ordered by time (e.g., in hourly chunks). In addition, sequential list elements must be linked by identical receiver pairs (i.e., the final receiver at which the individual was detected for any given chunk must be the same as the receiver at which the individual was next detected at the start of the next chunk) because it is only in this specific scenario that information does not need to be shared across time steps (see \code{split}). The algorithm will be implemented on each dataframe, termed `chunk', either in sequence or parallel. Any empty or \code{NULL} elements will be removed automatically.
#' @param archival For the ACDC algorithm, \code{archival} is a dataframe that contains depth time series (see \code{\link[flapper]{.acs}}).
#' @param step A number that defines the time step length (s) between consecutive detections (see \code{\link[flapper]{.acs}}).
#' @param plot_ts A logical input that defines whether or not to plot movement time series (see \code{\link[flapper]{.acs}}).
#' @param bathy A \code{\link[raster]{raster}} that defines the area (for the AC algorithm) or bathymetry (for the ACDC algorithm) across the area within which the individual could have moved (see \code{\link[flapper]{.acs}}).
#' @param detection_containers A list of detection containers (see \code{\link[flapper]{.acs}}).
#' @param detection_kernels A named list of detection probability kernels (see \code{\link[flapper]{.acs}}).
#' @param detection_kernels_overlap A named list of detection probability kernel overlaps, directly from \code{\link[flapper]{get_detection_containers_overlap}}. This must contain an element named `list_by_receiver' with the data for each receiver.
#' @param detection_time_window A number that defines the detection time window (see \code{\link[flapper]{.acs}})
#' @param mobility The mobility parameter (see \code{\link[flapper]{.acs}}).
#' @param calc_depth_error The depth error function (see \code{\link[flapper]{.acs}}).
#' @param normalise A logical input that defines whether or not to normalise maps (see \code{\link[flapper]{.acs}}).
#' @param save_record_spatial An integer of the spatial layers to save (see \code{\link[flapper]{.acs}}).
#' @param write_record_spatial_for_pf A named list used to write time step-specific maps to file (see \code{\link[flapper]{.acs}}).
#' @param verbose A logical variable that defines whether or not to print messages to the console or to file to relay function progress. If \code{con = ""}, messages are printed to the console (which is only supported if the algorithm is not implemented in parallel: see below); otherwise, they are written to file (see below).
#' @param con If \code{verbose = TRUE}, \code{con} is character string defines how messages relaying function progress are returned. If \code{con = ""}, messages are printed to the console (unless redirected by \code{\link[base]{sink}}), an approach that is only implemented if the function is not implemented in parallel. Otherwise, \code{con} defines the directory into which to write .txt files, into which messages are written to relay function progress. This approach, rather than printing to the console, is recommended for clarity, speed and debugging. If the algorithm is implemented step-wise, then a single file is written to the specified directory named acdc_log.txt. If the algorithm is implemented chunk-wise, then an additional file is written for each chunk (named dot_acdc_log_1.txt, dot_acdc_log_2.txt and so on), with the details for each chunk.
#' @param progress (optional) If the algorithm is implemented step-wise, \code{progress} is an integer (\code{1}, \code{2} or \code{3}) that defines whether or not to display a progress bar in the console as the algorithm moves over acoustic time steps (\code{1}), the `archival' time steps between each pair of acoustic detections (\code{2}) or both acoustic and archival time steps (\code{3}), in which case the overall acoustic progress bar is punctuated by an archival progress bar for each pair of acoustic detections. This option is useful if there is a large number of archival observations between acoustic detections. Any other input will suppress the progress bar. If the algorithm is implemented for chunks, inputs to \code{progress} are ignored and a single progress bar is shown of the progress across acoustic chunks.
#' @param split A character string that defines the (approximate) time unit used to split acoustic time series into chunks (e.g., \code{"12 hours"}). If provided, this must be supported by \code{\link[base]{cut.POSIXt}} (otherwise, a pre-defined list of acoustic time series can be passed to \code{acoustics}, e.g., specifying seasonal chunks). If \code{split = NULL} and a cluster has been specified (see \code{cl}) (and \code{acoustics} is a dataframe), then the acoustic time series is automatically split into chunks and the algorithm implemented for each chunk in parallel. In all cases, splitting is subject to the constraint that  chunks must join at identical receiver pairs (i.e., the last receiver at which the individual was detected on one chunk must match the first receiver at which the individual was next detected at the start of the next chunk): in these specific scenarios, information does not need to transfer from one time step to the next.
#' @param cl,varlist (optional) Parallelisation options. \code{cl} is (a) a cluster object from \code{\link[parallel]{makeCluster}} or (b) an integer that defines the number of child processes to implement the algorithm in parallel. If supplied, the algorithm is implemented for each chunk in a list of acoustic time series, either (a) as supplied by the user (if \code{acoustics} is a list), (b) as defined by the input to \code{split}, or (c) as defined automatically from the number of nodes in the cluster if \code{split = NULL}. If \code{cl} is supplied, \code{varlist} may also be required. This is a character vector of objects to export (see \code{\link[flapper]{cl_export}}). Exported variables must be located in the global environment. If a cluster is supplied, the connection to the cluster is closed within the function (see \code{\link[flapper]{cl_stop}}). For further information, see \code{\link[flapper]{cl_lapply}} and \code{\link[flapper]{flapper-tips-parallel}}.
#'
#' @return The function returns an \code{\link[flapper]{acdc_archive-class}} object. If a connection to write files has also been specified, an overall log (acdc_log.txt) as well as chunk-specific logs from calls to \code{\link[flapper]{.acs}}, if applicable, are written to file.
#'
#' @seealso The front-end functions \code{\link[flapper]{ac}} and \code{\link[flapper]{acdc}} call this function, which in turn calls \code{\link[flapper]{.acs}}. \code{\link[flapper]{acs_setup_containers}} defines the detection containers required by this function. \code{\link[flapper]{acs_setup_mobility}} is used to examine the assumption of the constant `mobility' parameter. \code{\link[flapper]{acs_setup_detection_kernels}} produces detection probability kernels for incorporation into the function. For calls via \code{\link[flapper]{ac}} and \code{\link[flapper]{acdc}}, \code{\link[flapper]{acdc_simplify}} simplifies the outputs and \code{\link[flapper]{acdc_plot_trace}}, \code{\link[flapper]{acdc_plot_record}} and \code{\link[flapper]{acdc_animate_record}} visualise the results.
#'
#' @examples
#' # For examples, see ?ac and ?acdc which call this function directly.
#'
#' @author Edward Lavender
#' @keywords internal
#'

.acs_pl <- function(
    acoustics,
    archival = NULL,
    step = 120,
    plot_ts = TRUE,
    bathy,
    detection_containers,
    detection_kernels = NULL, detection_kernels_overlap = NULL, detection_time_window = 5,
    mobility,
    calc_depth_error = function(...) matrix(c(-2.5, 2.5), nrow = 2),
    normalise = TRUE,
    save_record_spatial = 1L,
    write_record_spatial_for_pf = NULL,
    verbose = TRUE,
    con = "",
    progress = 1L,
    split = NULL,
    cl = NULL,
    varlist = NULL) {
  ######################################
  #### Set up

  #### A list to store overall outputs
  t_onset <- Sys.time()
  out <- list(archive = NULL, ts_by_chunk = NULL, time = NULL, args = NULL)

  #### Check parallelisation options
  n_cores <- cl_cores(cl)
  if (inherits(acoustics, "list") & !is.null(split)) {
    message("Input to 'split' ignored since inherits(acoustics, 'list') == TRUE.")
  }

  #### Define function for printing messages to file or console
  ## Check the connection for writing files, if applicable
  if (con != "") {
    if (!verbose) {
      message("Input to 'con' ignored since verbose = FALSE.")
    } else {
      # Check directory
      con <- check_dir(input = con, check_slash = TRUE)
      con_dir <- con
      # Define file
      con <- paste0(con_dir, "acdc_log.txt")
      if (!file.exists(con)) {
        message(paste0(con, " does not exist: attempting to write file in specified directory..."))
        file.create(file1 = con)
        message("... Blank file successfully written to file.")
      }
    }
  } else {
    if (verbose & n_cores > 1) stop("con = '' is not implemented in parallel (!is.null(cl)). Please supply a directory.")
  }
  ## Define function
  append_messages <- ifelse(con == "", FALSE, TRUE)
  cat_to_cf <- function(..., message = verbose, file = con, append = append_messages) {
    if (message) cat(paste(..., "\n"), file = con, append = append)
  }

  #### Checks
  ## Formally initiate function and implement remaining checks
  cat_to_cf(paste0("flapper::.acs_pl() called (@ ", t_onset, ")..."))
  out$time <- data.frame(event = "onset", time = t_onset)
  cat_to_cf("... Checking user inputs...")
  # Check acoustics contains required column names and correct variable types
  if (!inherits(acoustics, "list")) acoustics_tmp <- list(acoustics) else acoustics_tmp <- acoustics
  length_acoustics_tmp <- length(acoustics_tmp)
  n_obs_by_chunk <-
    lapply(acoustics_tmp, function(acc) {
      check_names(
        arg = "acoustics",
        input = acc,
        req = c("timestamp", "receiver_id", "index"),
        extract_names = colnames,
        type = all
      )
      check_class(input = acc$timestamp, to_class = "POSIXct", type = "stop")
      check_class(input = acc$receiver_id, to_class = "integer", type = "stop")
      return(nrow(acc))
    })
  if (min(acoustics_tmp[[1]]$index) != 1L) {
    stop("The first 'index' value in 'acoustics' should be 1L.", call. = FALSE)
  }
  # Check the number of observations
  n_obs <- sum(unlist(n_obs_by_chunk))
  if (n_obs <= 1L) stop("At least two observations are required to implement an AC* algorithm.", call. = FALSE)
  # Check the ratio of the number of observations to the number of cores
  # ... The number of cores should not exceed the maximum number of chunks.
  # ... Since a minimum of two observations per chunk is required,
  # ... the number of cores should be less than half of the number of observations.
  # ... If not, we will reset n_cores to a more appropriate value.
  # ... This helps to guide splitting (if necessary), but ultimately the function is still implemented
  # ... using the user-defined cluster or number of child processes.
  if (n_obs / n_cores < 2L) {
    warning("The number of specified cores exceeds the maximum number of chunks.", immediate. = TRUE, call. = FALSE)
    n_cores <- floor(n_obs / n_cores)
  }
  # Check archival contains required column names and correct variable types
  if (!is.null(archival)) {
    check_names(
      input = archival,
      req = c("timestamp", "depth"),
      extract_names = colnames,
      type = all
    )
    check_class(input = archival$timestamp, to_class = "POSIXct", type = "stop")
    archival$timestamp_num <- as.numeric(archival$timestamp)
    check_class(input = archival$depth, to_class = "numeric", type = "stop")
    # Check data volume
    if (nrow(archival) <= 1) stop("'archival' dataframe only contains one or fewer rows.")
    # Check archival step length
    step_est <- as.numeric(difftime(archival$timestamp[2], archival$timestamp[1], units = "s"))
    if (!all.equal(step, step_est)) {
      stop("'step' does not equal difftime(archival$timestamp[2], archival$timestamp[1], units = 's'.)")
    }
    acoustics_tmp_df <- do.call(rbind, acoustics_tmp)
    if (any(c(
      min(acoustics_tmp_df$timestamp) < min(archival$timestamp),
      max(acoustics_tmp_df$timestamp) > max(archival$timestamp)
    ))) {
      stop("Archival observations do not span the time range of acoustic observations.", call. = FALSE)
    }
  }
  # Check detection containers have been supplied as a list
  check_class(input = detection_containers, to_class = "list", type = "stop")
  out$time <- rbind(out$time, data.frame(event = "initial_checks_passed", time = Sys.time()))
  # Check detection_kernels_overlap
  if (!is.null(detection_kernels_overlap)) {
    if (!("list_by_receiver" %in% names(detection_kernels_overlap))) stop("'detection_kernels_overlap' must contain a 'list_by_receiver' element.")
    detection_kernels_overlap <- detection_kernels_overlap$list_by_receiver
    # For the detections that occurred at the same point in time, check they occurred at receivers with overlapping centroids
    # If not, this suggests detection centroids are too small
    # We want to catch this issue early to prevent errors following centroid intersections later.
    multiples <-
      acoustics |>
      dplyr::filter(.data$timestamp %in% .data$timestamp[duplicated(.data$timestamp)])
    if (nrow(multiples) > 0L) {
      check_overlaps <-
        sapply(split(multiples, multiples$timestamp), function(d) {
          issue <- FALSE
          .overlaps <- detection_kernels_overlap[[d$receiver_id[1]]]
          if (any(.overlaps[, colnames(.overlaps)[colnames(.overlaps) %in% d$receiver_id]] == 0L)) {
            issue <- TRUE
            warning(
              paste0(
                "The individual was detected at multiple receivers ('",
                paste0(d$receiver_id, collapse = "', '"),
                "') at the same moment in time ('",
                d$timestamp[1],
                "') but not all detection radii overlap."
              ),
              immediate. = TRUE, call. = FALSE
            )
          }
          issue
        })
      if (any(check_overlaps)) stop("Detection radii do not overlap.", call. = FALSE)
    }
  }
  # Check depth error
  if (!is.null(archival)) {
    de <- calc_depth_error(archival$depth)
    if (inherits(de, "matrix")) {
      if (nrow(de) == 2) {
        if (ncol(de) == 1) {
          message("'calc_depth_error' function taken to be independent of depth.")
        } else {
          message("'calc_depth_error' taken to depend on depth.")
        }
        if (any(de[1, ] > 0) | any(de[2, ] < 0)) stop("'calc_depth_error' should be a function that returns a two-row matrix with lower (negative) adjustment(s) (top row) and upper (positive) adjustment(s) (bottom row).'", call. = FALSE)
      } else {
        stop("'calc_depth_error' should return a two-row matrix.", call. = FALSE)
      }
    } else {
      stop("'calc_depth_error' should return a two-row matrix.", call. = FALSE)
    }
  }
  # Check write opts
  if (!is.null(write_record_spatial_for_pf)) {
    check_named_list(input = write_record_spatial_for_pf)
    check_names(input = write_record_spatial_for_pf, req = "filename")
    write_record_spatial_for_pf$filename <- check_dir(input = write_record_spatial_for_pf$filename, check_slash = TRUE)
    if (length(list.files(write_record_spatial_for_pf$filename)) != 0L) {
      warning("write_record_spatial_for_pf$filename' is not an empty directory.", immediate. = TRUE, call. = FALSE)
    }
  }

  #### Study site rasters
  ## Check for zeros
  if (is.null(calc_depth_error) && length(raster::Which(bathy == 0, cells = TRUE)) > 0L) {
    warning("'bathy' contains zero values.",
      immediate. = TRUE, call. = TRUE
    )
  }
  ## Blank map for space use over the study area
  map <- raster::setValues(bathy, 0)
  map <- raster::mask(map, bathy)


  ######################################
  #### Implement splitting (if necessary)

  #### Round time series to the nearest step
  round_to <- paste0(step / 60, "mins")
  if (inherits(acoustics, "list")) {
    acoustics <- lapply(1:length(acoustics), function(i) {
      acc <- acoustics[[i]]
      acc_timestamps_round <- lubridate::round_date(acc$timestamp, round_to)
      if (!isTRUE(all.equal(acc_timestamps_round, acc$timestamp))) {
        cat_to_cf(paste0("... ... acoustics[[", i, "]]$timestamp rounded to the nearest ", step, " second(s)."))
        message(paste0("acoustics[[", i, "]]$timestamp rounded to the nearest ", step, " second(s)."))
        acc$timestamp <- lubridate::round_date(acc$timestamp, round_to)
      }
      acc$dup <- duplicated(paste0(acc$receiver_id, "-", acc$timestamp))
      if (any(acc$dup)) {
        cat_to_cf(paste0("... ... Duplicate observations in acoustics[[", i, "]] identified and dropped."))
        message(paste0("Duplicate observations in acoustics[[", i, "]] identified and dropped."))
        acc <- acc[!acc$dup, ]
      }
      return(acc)
    })
  } else {
    acoustics_timestamps_round <- lubridate::round_date(acoustics$timestamp, round_to)
    if (!isTRUE(all.equal(acoustics_timestamps_round, acoustics$timestamp))) {
      cat_to_cf(paste("... ... acoustics$timestamp rounded to the nearest", step, "second(s)."))
      message(paste("acoustics$timestamp rounded to the nearest", step, "second(s)."))
      acoustics$timestamp <- lubridate::round_date(acoustics$timestamp, round_to)
    }
    acoustics$dup <- duplicated(paste0(acoustics$receiver_id, "-", acoustics$timestamp))
    if (any(acoustics$dup)) {
      cat_to_cf("... ... Duplicate observations in 'acoustics' identified and dropped.")
      message("Duplicate observations in 'acoustics' identified and dropped.")
      acoustics <- acoustics[!acoustics$dup, ]
    }
  }
  if (!is.null(archival)) {
    arc_timestamp_round <- lubridate::round_date(archival$timestamp, round_to)
    if (!isTRUE(all.equal(arc_timestamp_round, archival$timestamp))) {
      cat_to_cf(paste("... ... archival$timestamp rounded to the nearest", step, "second(s)."))
      message(paste("archival$timestamp rounded to the nearest", step, "second(s)."))
      archival$timestamp <- lubridate::round_date(archival$timestamp, round_to)
    }
    archival$dup <- duplicated(archival$timestamp)
    if (any(archival$dup)) {
      cat_to_cf("... ... Duplicate observations in 'archival' identified and dropped.")
      message(paste0("Duplicate observations in 'archival' identified and dropped."))
      archival <- archival[!archival$dup, ]
    }
  }
  # Re-define acoustics_tmp
  if (!inherits(acoustics, "list")) acoustics_tmp <- list(acoustics) else acoustics_tmp <- acoustics

  #### Define a list of dataframes
  # .. If the algorithm is to be implemented in parallel
  if (inherits(acoustics, "list") | n_cores > 1 | !is.null(split)) {
    #### Implement splitting
    if (length_acoustics_tmp == 1) {
      cat_to_cf("... Splitting 'acoustics' into chunks...")

      ## Define 'permitted' splitting positions
      # ... Splitting is only permitted at positions
      # ... when the 'previous' receiver == the 'current' receiver
      # ... because at these time steps information does not need
      # ... to pass from the previous chunk to the next chunk).
      # ... Splitting is also subject to the constraint that each chunk should comprise
      # ... at least three observations. Currently, this constraint is not enforced here,
      # ... but issues are captured below.
      acoustics$split_permit <- dplyr::lag(acoustics$receiver_id) == acoustics$receiver_id
      acoustics$split_permit[1] <- FALSE
      if (!any(acoustics$split_permit)) {
        stop("The algorithm cannot be implemented chunkwise for this dataset.", call. = FALSE)
      }

      ## Define (approximate) splitting time interval ('split') if unspecified
      # If 'split' hasn't been specified, then the user doesn't care
      # ... about splitting the outputs into biologically interpretable time intervals
      # ... so we will simply select an interval based on the desired number of chunks
      if (is.null(split)) {
        dft <- ceiling((max(acoustics$timestamp) - min(acoustics$timestamp)) / n_cores)
        dft_num <- as.numeric(dft)
        dft_units <- attr(dft, "units")
        split <- paste(dft_num, dft_units)
      }

      ## Split dataframe
      # Define suitable splitting positions
      acoustics$bin <- cut(acoustics$timestamp, split)
      acc_by_bin <- split(acoustics, acoustics$bin)
      split_ind <- lapply(2:length(acc_by_bin), function(i) {
        acc_for_bin <- acc_by_bin[[i]]
        pos <- which(acc_for_bin$split_permit)
        if (length(pos) >= 1L) {
          return(acc_for_bin$index[min(pos)])
        } else {
          return(NULL)
        }
      })
      split_ind <- unlist(compact(split_ind))
      # Split acoustics at specified positions
      acoustics_ls <- split(acoustics, findInterval(seq_len(nrow(acoustics)), split_ind))
      if (nrow(acoustics_ls[[length(acoustics_ls)]]) < 2L) {
        ind <- length(acoustics_ls)
        acoustics_ls[[ind - 1]] <- rbind(acoustics_ls[[ind - 1]], acoustics_ls[[ind]])
        acoustics_ls[[ind]] <- NULL
      }

      ## Report the number and duration of chunks
      if (length(acoustics_ls) > 1L) {
        dft <- difftime(max(acoustics_ls[[1]]$timestamp), min(acoustics_ls[[1]]$timestamp))
        dft_units <- attr(dft, "units")
        dft_by_chunk <-
          sapply(
            acoustics_ls,
            function(elm) {
              as.numeric(difftime(max(elm$timestamp),
                min(elm$timestamp),
                units = dft_units
              ))
            }
          )
        dft_med <- round(stats::median(dft_by_chunk))
        dft_min <- round(min(dft_by_chunk))
        dft_max <- round(max(dft_by_chunk))
        message(paste0(
          "'acoustics' dataframe split into ", length(acoustics_ls), " chunks of ~",
          dft_med, " (", dft_min, "--", dft_max, ") ", dft_units,
          " across ", n_cores, " core(s)."
        ))
      }
    } else {
      acoustics_ls <- acoustics
    }

    #### Process split dataframes
    cat_to_cf("... Processing acoustics chunks...")

    ## Remove NULL/length 0 elements
    cat_to_cf("... ... Checking for NULL/empty chunks...")
    empty_elms <- sapply(acoustics_ls, function(x) is.null(x) | nrow(x) == 0)
    if (any(empty_elms)) {
      msg <- paste0("acoustics_ls[c(", paste0(which(empty_elms), collapse = ","), ")] chunks are empty/NULL and will be removed...")
      message(msg)
      cat_to_cf(paste("... ... ...", msg))
      acoustics_ls <- acoustics_ls[which(!empty_elms)]
    }

    ## Check splitting with respect to receiver IDs
    # Splitting is only supported between identical receiver pairs
    # ... because it is only for these pairs at which information from the previous time step
    # ... is not required for the next time step (and information can't be shared across chunks).
    for (i in 2:length(acoustics_ls)) {
      acc_1 <- acoustics_ls[[i - 1]]
      acc_2 <- acoustics_ls[[i]]
      acc_1_r <- acc_1[nrow(acc_1), "receiver_id"]
      acc_2_r <- acc_2[1, "receiver_id"]
      if (acc_1_r != acc_2_r) {
        stop("'acoustics' must be split at identical receiver pairs.", call. = FALSE)
      }
    }

    ## Force overlapping time series
    # If we naively split the dataframe into number of different windows,
    # ... on every run, we have to stop before the last acoustic reading
    # ... (because we can't identify the next receiver - their isn't one in the split dataframe)
    # ...which means we're not including some information when we estimate space use.
    # To get around this, in the list dataframes, we need to add the first line of every dataframe
    # ... to the previous dataframe. Then, when we split the dataframe, we won't be loosing information
    # ...because we've copied the last line.
    cat_to_cf("... ... Overlapping chunks...")
    if (length(acoustics_ls) > 1) {
      acoustics_ls_wth_overlap <-
        lapply(2:(length(acoustics_ls)), function(i) {
          # define an adjusted dataframe, binds the previous dataframe
          # ... with the first row of the dataframe in question:
          adj <- rbind(acoustics_ls[[i - 1]], acoustics_ls[[i]][1, ])
          return(adj)
        })
      # Add back the final element:
      acoustics_ls_wth_overlap[[length(acoustics_ls)]] <- acoustics_ls[[length(acoustics_ls)]]
      names(acoustics_ls_wth_overlap) <- names(acoustics_ls)
    } else {
      acoustics_ls_wth_overlap <- acoustics_ls
    }

    #### Additional checks
    # Check the number of rows in acoustics_ls_wth_overlap. This cannot be less than two.
    # ... If there are no rows or only one row,
    # ... then we can't calculate where the individual was
    # ... next detected, which will cause problems.
    l <- length(acoustics_ls_wth_overlap)
    lapply(1:l, function(i) {
      nrw <- nrow(acoustics_ls_wth_overlap[[i]])
      if (nrw < 2) {
        stop(paste("acoustics_ls_wth_overlap[[", i, "]] has less than two rows. This is not allowed."))
      }
    })

    out$time <- rbind(out$time, data.frame(event = "acoustics_chunks_defined", time = Sys.time()))
  } else {
    acoustics_ls_wth_overlap <- acoustics_tmp
  }


  ######################################
  #### Visualise time series

  #### Focus on the data for which we have both acoustic and archival observations

  cat_to_cf("... Processing movement time series...")
  movement_ts <- lapply(1:length(acoustics_ls_wth_overlap), function(i) {
    # Isolate acoustics time series
    acc <- acoustics_ls_wth_overlap[[i]]
    acc$timestamp_num <- as.numeric(acc$timestamp)

    ## AC algorithm implementation
    if (is.null(archival)) {
      ls <- list(acoustics = acc, archival = NULL)

      ## ACDC algorithm implementation
    } else {
      nrw_acc_pre <- nrow(acc)
      nrw_arc_pre <- nrow(archival)
      acc <- acc[acc$timestamp >= min(archival$timestamp) - 2 * 60 &
        acc$timestamp <= max(archival$timestamp) + 2 * 60, ]
      arc <- archival[archival$timestamp >= min(acc$timestamp) - 2 * 60, ]
      if (i < length(acoustics_ls_wth_overlap)) {
        arc <- arc[arc$timestamp <= max(acc$timestamp) + 2 * 60, ]
      }
      nrw_acc_post <- nrow(acc)
      nrw_arc_post <- nrow(arc)
      nrw_acc_delta <- nrw_acc_pre - nrw_acc_post
      nrw_arc_delta <- nrw_arc_pre - nrw_arc_post
      if (nrw_acc_post == 0 | nrw_arc_post == 0) stop("No overlapping acoustic/archival observations to implement algorithm.")
      if (nrw_acc_delta != 0) {
        cat_to_cf(paste("... ...  Chunk", i, ":", nrw_acc_delta, "acoustic observation(s) beyond the ranges of archival observations ignored."))
        message(paste("Chunk", i, ":", nrw_acc_delta, "acoustic observation(s) beyond the ranges of archival observations ignored."))
      }
      if (nrw_arc_delta != 0) {
        cat_to_cf(paste("... ... Chunk", i, ":", nrw_arc_delta, "archival observation(s) beyond the ranges of (processed) acoustic detections ignored."))
        message(paste("Chunk", i, ":", nrw_arc_delta, "archival observation(s) beyond the ranges of (processed) acoustic detections ignored."))
      }

      ## Return processed time series
      ls <- list(acoustics = acc, archival = arc)
    }
    return(ls)
  })
  out$ts_by_chunk <- movement_ts
  out$time <- rbind(out$time, data.frame(event = "movement_time_series_processed", time = Sys.time()))

  #### Visualise processed time series
  if (plot_ts) {
    cat_to_cf("... Plotting movement time series (for each chunk)...")
    if (length(movement_ts) < 25) pp <- graphics::par(mfrow = prettyGraphics::par_mf(length(movement_ts)))
    lapply(movement_ts, function(move) {
      ## Get acoustics time series
      acoustics <- move$acoustics
      ## AC algorithm implementation
      if (is.null(archival)) {
        prettyGraphics::pretty_line(acoustics$timestamp,
          pch = 21, col = "royalblue", bg = "royalblue"
        )
        ## ACDC algorithm implementation
      } else {
        archival <- move$archival
        axis_ls <- prettyGraphics::pretty_plot(archival$timestamp, abs(archival$depth) * -1,
          pretty_axis_args = list(
            side = 3:2,
            axis = list(
              list(format = "%H:%M:%S %d-%m-%y"),
              list()
            )
          ),
          xlab = "Time stamp", ylab = "Depth (m)",
          type = "l"
        )
        prettyGraphics::pretty_line(acoustics$timestamp,
          pretty_axis_args = list(axis_ls = axis_ls),
          inherit = TRUE,
          replace_axis = list(side = 1, pos = axis_ls[[2]]$lim[1]),
          add = TRUE,
          pch = 21, col = "royalblue", bg = "royalblue"
        )
      }
    })

    if (length(movement_ts) < 25) graphics::par(pp)
    out$time <- rbind(out$time, data.frame(event = "time_series_plotted", time = Sys.time()))
  }


  ######################################
  #### Implement ACDC algorithm

  #### Checks
  n_chunks <- length(acoustics_ls_wth_overlap)
  # Define a list of files, one for each chunk
  if (verbose & con != "") {
    con_ls <- lapply(1:n_chunks, function(i) {
      file <- paste0(con_dir, "dot_acdc_log_", i, ".txt")
      return(file)
    })
  } else {
    con_ls <- lapply(1:n_chunks, function(i) {
      return("")
    })
  }

  # Write blank files to directory if required
  if (verbose & con != "" & n_chunks > 1) {
    cat_to_cf("... Defining chunk-specific log files as dot_acdc_log_1.txt, dot_acdc_log_2.txt etc...")
    lapply(con_ls, function(file) {
      if (!file.exists(file)) {
        msg1 <- paste(file, "does not exist: attempting to write file in specified directory...")
        cat_to_cf(paste("... ...", msg1))
        message(msg1)
        file.create(file1 = file)
        cat_to_cf("... ... ... Blank file successfully written to file.")
        message("... Blank file successfully written to file.")
      }
    })
  }

  #### Implement ACDC algorithm directly via .acs back-end
  if (length(acoustics_ls_wth_overlap) == 1) {
    #### Implement algorithm
    cat_to_cf("... Calling .acs() to implement ACDC algorithm on one chunk...")
    out$time <- rbind(out$time, data.frame(event = "calling_.acs()", time = Sys.time()))
    .out <- .acs(
      acoustics = movement_ts[[1]]$acoustics,
      archival = movement_ts[[1]]$archival,
      step = step,
      plot_ts = FALSE,
      round_ts = FALSE,
      bathy = bathy,
      map = map,
      detection_containers = detection_containers,
      detection_kernels = detection_kernels,
      detection_kernels_overlap = detection_kernels_overlap,
      detection_time_window = detection_time_window,
      mobility = mobility,
      calc_depth_error = calc_depth_error,
      normalise = normalise,
      save_record_spatial = save_record_spatial,
      write_record_spatial_for_pf = write_record_spatial_for_pf,
      save_args = FALSE,
      verbose = verbose,
      con = con,
      progress = progress,
      check = FALSE
    )
    .out <- list(.out)
  } else {
    #### Implement algorithm in parallel
    cat_to_cf(paste("... Calling .acs() to implement ACDC algorithm on", length(acoustics_ls_wth_overlap), "chunks, using", n_cores, "cores..."))
    out$time <- rbind(out$time, data.frame(event = "calling_.acs()", time = Sys.time()))
    .out <- cl_lapply(1:length(acoustics_ls_wth_overlap), cl = cl, varlist = varlist, function(i) {
      #### Implement algorithm
      .out <- .acs(
        acoustics = movement_ts[[i]]$acoustics,
        archival = movement_ts[[i]]$archival,
        plot_ts = FALSE,
        round_ts = FALSE,
        step = step,
        bathy = bathy,
        map = map,
        detection_containers = detection_containers,
        detection_kernels = detection_kernels,
        detection_kernels_overlap = detection_kernels_overlap,
        detection_time_window = detection_time_window,
        mobility = mobility,
        calc_depth_error = calc_depth_error,
        normalise = normalise,
        save_record_spatial = save_record_spatial,
        write_record_spatial_for_pf = write_record_spatial_for_pf,
        save_args = FALSE,
        chunk = i,
        verbose = verbose,
        con = con_ls[[i]],
        progress = 0L,
        check = FALSE
      )
      return(.out)
    })
  }

  #### Return outputs
  names(.out) <- 1:length(.out)
  out$archive <- .out
  t_end <- Sys.time()
  out$time <- rbind(out$time, data.frame(event = "algorithm_competion", time = t_end))
  out$time$serial_duration <- Tools4ETS::serial_difference(out$time$time, units = "mins")
  out$time$total_duration <- NA
  total_duration <- sum(as.numeric(out$time$serial_duration), na.rm = TRUE)
  out$time$total_duration[nrow(out$time)] <- total_duration
  cat_to_cf(paste0("... flapper::.acs_pl() call completed (@ ", t_end, ") after ~", round(total_duration, digits = 2), " minutes."))
  class(out) <- c(class(out), "acdc_archive")
  return(out)
}
