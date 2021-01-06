#' @title Convert a detection matrix into a dataframe
#' @description This function converts a matrix of detections (0, 1) by time stamp and receiver into a dataframe.
#'
#' @param acoustics A detection matrix (time stamps x receivers) in which the cells define whether (1) or not (0) a detection was made at each time stamp/receiver combination. 'Meaningful' time stamps and receiver IDs can be taken from the row and column names of this matrix, if specified (see \code{set_names}).
#' @param only_keep_detections A logical variable that defines whether or not to retain only observations that correspond to detections. (If \code{only_keep_detections = FALSE} the returned dataframe includes time stamps without detections.)
#' @param set_names A logical variable that defines whether or not to take the row and column names of \code{acoustics} as the time stamps and receiver IDs. (If \code{set_names = FALSE}, time stamps and receiver IDs are simply given as integer vectors of 1 to the number of rows or columns respectively.)
#' @param as_POSIXct If \code{set_names = TRUE}, \code{as_POSIXct} is a function that converts the row names of \code{acoustics} into POSIXct time stamps.
#'
#' @return The function returns a dataframe with time stamps ('timestamp') and receivers ('receiver_id'). If \code{set_names = FALSE}, these are integer vectors that match the dimensions of \code{acoustics}. Otherwise they are are taken from row and column names of \code{acoustics}. In this case, if \code{as_POSIXct} is defined, time stamps are returned in POSIXct format and receivers are returned as a factor. If \code{only_keep_detections = FALSE}, the dataframe also includes a 'detection' column that defines whether (1) or not (0) a detection was made for each observation; otherwise, this column is dropped (mirroring real-world data).
#'
#' @examples
#' #### Define detection matrix
#' # Simulate array
#' array <- sim_array(boundaries = raster::extent(-1000, 1000, -1000, 1000),
#'                    n_receivers = 24, seed = 1)
#' # Simulate movement in this area
#' path <- sim_path_sa(n = 50, area = array$array$area, seed = 1)
#' # Simulate a detection matrix
#' detections <- sim_detections(n = 100,
#'                              path = path$xy_mat,
#'                              xy = sp::coordinates(array$array$xy),
#'                              calc_detection_pr = function(dist) ifelse(dist < 425, 1, 0),
#'                              )
#' # Extract matrix
#' mat <- detections$det_mat
#' # Define row names
#' rownames(mat) <-
#'   as.character(
#'     seq(as.POSIXct("2016-01-01"), by = "2 mins", length.out = nrow(mat))
#'   )
#'
#' #### Examples: convert the matrix to a dataframe
#' utils::str(mat)
#' dat <- make_df_detections(mat)
#' utils::str(dat)
#' dat <- make_df_detections(mat, only_keep_detections = TRUE)
#' utils::str(dat)
#' dat <- make_df_detections(mat, only_keep_detections = TRUE, set_names = TRUE)
#' utils::str(dat)
#'
#' @author Edward Lavender
#' @export
#'

make_df_detections <- function(acoustics, only_keep_detections = FALSE, set_names = FALSE, as_POSIXct = as.POSIXct){
  nrw <- nrow(acoustics)
  dat_ls <- lapply(1:ncol(acoustics), function(i) {
    d <- data.frame(timestamp = 1:nrw, receiver_id = i, detection = acoustics[, i])
    return(d)
  })
  dat <- do.call(rbind, dat_ls)
  rownames(dat) <- NULL
  if(only_keep_detections) {
    dat <- dat[dat$detection == 1, ]
    dat$detection <- NULL
  }
  if(set_names) {
    if(!is.null(rownames(acoustics))) {
      dat$timestamp <- rownames(acoustics)[match(dat$timestamp, 1:nrow(acoustics))]
      if(!is.null(as_POSIXct)) dat$timestamp <- as_POSIXct(dat$timestamp)
    } else message("'set_names' not implemented for rows: 'acoustics' does not contain row names.")
    if(!is.null(colnames(acoustics))) {
      dat$receiver_id <- factor(colnames(acoustics))
    } else message("'set_names' not implemented for columns: 'acoustics' does not contain column names.")
  }
  dat <- dat[order(dat$timestamp, dat$receiver), ]
  return(dat)
}


