########################################
########################################
#### pf_setup_movement_pr()

#' @title A simple movement model dependent on distance
#' @description This function provides a simple movement model that calculates the probability of movement between two locations according to the distance between them, using an logistic equation with pre-defined parameters.
#'
#' @param distance A numeric vector of distances (m).
#' @param ... Additional arguments (none implemented).
#'
#' @details Under this model, for distance(s) \eqn{ \leq 500 } m, \eqn{Pr(distance) = logistic(10 + distance  -0.05)}; otherwise, \eqn{Pr(distance) = 0}. This particular model is designed for flapper skate (\emph{Dipturus intermedius}) and represents a reasonable model for the probability of moving a given distance in a two-minute period (in the absence of additional information).
#'
#' @return The function returns a numeric vector of probabilities that represent the probability of movement between two or more areas given the distances between them.
#'
#' @examples
#' pr <- pf_setup_movement_pr(1:1000)
#' plot(pr, type = "l", xlab = "Distance (m)", ylab = "Pr(distance)")
#' @seealso This function is used as the default movement model in \code{\link[flapper]{pf}}.
#' @author Edward Lavender
#' @export

pf_setup_movement_pr <- function(distance,...) {
  pr <- stats::plogis(10 + distance * -0.05)
  pr[distance > 500] <- 0
  return(pr)
}


######################################
######################################
#### pf_setup_record()

#' @title List `record' files from an AC/DC algorithm for PF
#' @description This function creates an ordered list of `record' files derived from the AC/DC/ACDC algorithm (\code{\link[flapper]{ac}}, \code{\link[flapper]{dc}} or \code{\link[flapper]{acdc}}) for particle filtering (PF) via \code{\link[flapper]{pf}}.
#'
#' @param root A string that defines the directory in which files are loaded.
#' @param type A character that defines the source of the files (\code{type = "acs"} refers to an AC* algorithm and \code{type = "dc"} refers to the DC algorithm).
#' @param use_absolute_paths A logical variable that defines whether to return relative paths (\code{FALSE}) or absolute paths (\code{TRUE}) (see \code{\link[tools]{file_path_as_absolute}}).
#' @param ... Additional arguments passed to \code{\link[base]{list.files}} (excluding \code{full.names}).
#'
#' @details This function requires the \code{\link[stringr]{stringr}} package.
#'
#' @return The function returns an ordered list of file paths.
#'
#' @examples
#' #### Example (1): Example with the AC algorithm
#' # Define a directory in which to save files
#' root <- paste0(tempdir(), "/ac/")
#' dir.create(root)
#' # Implement the AC algorithm for some example time series
#' acc <- dat_acoustics[dat_acoustics$individual_id == 25, ][1:5, ]
#' out_ac <- ac(acoustics = acc,
#'              step = 120,
#'              bathy = dat_gebco,
#'              detection_centroids = dat_centroids,
#'              mobility = 250,
#'              write_record_spatial_for_pf = list(filename = root))
#' # List the files for pf()
#' files <- pf_setup_record(root, type = "acs", pattern = "*.grd")
#' utils::head(files)
#' # Implement pf() using files (not shown).
#'
#' #### Example (2): Example with the DC algorithm
#' # Define a directory in which to save files
#' root <- paste0(tempdir(), "/dc/")
#' dir.create(root)
#' # Implement the DC algorithm for some example time series
#' depth <- dat_archival[dat_archival$individual_id == 25, ][1:5, ]
#' out_dc <- dc(archival = depth,
#'              bathy = dat_gebco,
#'              write_record_spatial_for_pf = list(filename = root))
#' # List the files for pf()
#' files <- pf_setup_record(root, type = "dc", pattern = "*.grd")
#' utils::head(files)
#' # Implement pf() using files (not shown).
#'
#' @seealso This function is designed to list outputs from \code{\link[flapper]{ac}}, \code{\link[flapper]{dc}} or \code{\link[flapper]{acdc}} (see the \code{write_record_spatial_for_pf} argument) as required by \code{\link[flapper]{pf}} (see the \code{record} argument).
#' @author Edward Lavender
#' @export

pf_setup_record <- function(root, type = c("acs", "dc"), use_absolute_paths = FALSE,...){
  if(!requireNamespace("stringr", quietly = TRUE)){
    stop("This function requires the 'stringr' package. Please install it before continuing with install.packages('stringr').")
  }
  check...("full.names",...)
  type <- match.arg(type)
  check_dir(input = root)
  files <- list.files(root,...)
  msg_unrecognised <- "File naming structure is unrecognised."
  if(type == "acs"){
    if(!grepl("chu", files[1], fixed = TRUE)){
      warning("File naming structure is unrecognised.", immediate. = TRUE)
      if(grepl("arc", files[1], fixed = TRUE)){
        if(utils::askYesNo("...Did you mean type = 'dc'?")){
          type <- "dc"
        } else stop(msg_unrecognised)
      } else stop(msg_unrecognised)
    }
  } else if(type == "dc"){
    if(!grepl("arc", files[1], fixed = TRUE)){
      stop(msg_unrecognised)
    }
  }
  if(length(unique(tools::file_ext(files))) != 1L)
    warning("Multiple file types (extensions) identified in 'root'. Did you forget to pass 'pattern' to list.files()?",
            immediate. = TRUE, call. = FALSE)
  files <- data.frame(index = 1:length(files), name = files)
  if(type == "acs"){
    files[, c("chu_id", "acc_id", "arc_id")] <- stringr::str_split_fixed(files$name, "_", 6)[, c(2, 4, 6)]
    files$chu_id <- as.integer(files$chu_id)
    files$acc_id <- as.integer(files$acc_id)
  } else if(type == "dc"){
    files[, "arc_id"] <- stringr::str_split_fixed(files$name, "_", 2)[, 2]
  }
  ext <- tools::file_ext(files$name)
  n <- nchar(ext) + 1
  files$arc_id <- as.integer(substr(files$arc_id, 1, nchar(files$arc_id)-n))
  if(type == "acs"){
    files <- files %>% dplyr::arrange(.data$chu_id, .data$acc_id, .data$arc_id)
  } else if(type == "dc"){
    files <- files %>% dplyr::arrange(.data$arc_id)
  }
  files <- list.files(root, full.names = TRUE,...)[files$index]
  if(use_absolute_paths) {
    files <- sapply(files, function(f) tools::file_path_as_absolute(f))
    names(files) <- NULL
  }
  return(files)
}


######################################
######################################
#### pf_setup_optimisers()

#' @title Optimisation settings for \code{\link[flapper]{pf}}
#' @description This function defines optimisation settings for \code{\link[flapper]{pf}}. These settings control under-the-hood implementation routines in \code{\link[flapper]{pf}} that may improve computation time if adjusted.
#'
#' @param use_raster_operations (experimental) A logical input that defines whether or not to use \code{\link[raster]{raster}} operations, where applicable (e.g., \code{\link[raster]{calc}}), which are memory-safe, or to extract \code{\link[raster]{raster}} values into a \code{\link[data.table]{data.table}} and perform arithmetic operations on the \code{\link[data.table]{data.table}}. This option is only implemented for the `fast Euclidean distances' method in \code{\link[flapper]{pf}}. Trials suggest that \code{use_raster_operations = FALSE} does not improve computation time.
#' @param use_calc_distance_euclid_backend_grass A logical input that defines whether or not to use GRASS as the backend for Euclidean distances calculations in \code{\link[flapper]{pf}}. The default is \code{FALSE}, in which case \code{\link[raster]{distanceFromPoints}} is used for these calculations. If \code{TRUE}, the \code{\link[fasterRaster]{fasterRaster}} package is required and \code{\link[fasterRaster]{fasterVectToRastDistance}} is used instead.
#' @param use_grass_dir If \code{use_calc_distance_euclid_backend_grass = TRUE}, \code{use_grass_dir} is a character that defines the directory where GRASS is installed on your system and should be supplied.
#'
#' @details \code{\link[flapper]{pf}} is a computationally intensive routine. To reduce computation time, the most effective approaches are to minimise data volume and reduce the size (dimensions and/or resolution) of the grid over which particle filtering is implemented; use the `fast Euclidean distances' method for distance calculations; and minimise the number of particles. For small numbers of particles, it may be faster to specify the \code{mobility} parameter; for large numbers of particles, it is probably faster to set \code{mobility = NULL}. Adjusting \code{\link{raster}{rasterOptions}} such as \code{chunksize} and/or \code{maxmemory} may help in some circumstances too. Following optimisation of these settings, \code{\link[flapper]{pf_setup_optimisers}} facilitates the adjustment of under-the-hood implementation routines which may further reduce computation time in some settings.
#'
#' @return The function returns \code{pf_optimiser} S3 class object, which is simply a named list of optimisation options that can be passed to \code{\link[flapper]{pf}} via the \code{optimisers} argument.
#'
#' @examples
#' #### Example (1): The default implementation
#' pf_setup_optimisers()
#'
#' #### Example (2): Use GRASS for Euclidean distance calculations
#' # Specification for GRASS-7.4.4 on MacOS
#' pf_setup_optimisers(use_calc_distance_euclid_backend_grass = TRUE,
#'                     use_grass_dir = "/Applications/GRASS-7.4.4.app/Contents/Resources")
#' # This list should be passed to the 'optimisers' argument in pf().
#'
#' @seealso \code{\link[flapper]{pf}}
#' @author Edward Lavender
#' @export

pf_setup_optimisers <- function(use_raster_operations = TRUE,
                                use_calc_distance_euclid_backend_grass = FALSE,
                                use_grass_dir = NULL){
  out <- list(use_raster_operations = use_raster_operations,
              use_calc_distance_euclid_backend_grass = use_calc_distance_euclid_backend_grass,
              use_grass_dir = use_grass_dir)
  if(use_calc_distance_euclid_backend_grass & is.null(use_grass_dir)){
    warning("'use_calc_distance_euclid_backend_grass' specified but 'use_grass_dir' is NULL.",
            immediate. = TRUE, call. = FALSE)
  }
  class(out) <- c(class(out), "pf_optimiser")
  return(out)
}
