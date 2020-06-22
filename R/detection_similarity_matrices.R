#######################################
#######################################
#### compute_det_sim()

#' @title Compute a detection history similarity matrix
#' @description The function computes a detection history similarity matrix. For all combinations of individuals, this shows the total number (or percentage) of detections 'nearby' in space and time, which can help elucidate possible interactions among individuals that affect space use (see Details). To compute this matrix, the function pairs detections for each individual with the detections nearest in time for each other individual. The function computes the time (minutes) between paired detection timeseries, and the distance (m) between the receiver(s) at which paired detections occurred, dropping any detection pairs that are further apart in time or space than user-defined thresholds (which depend on the mobility of the species under investigation). For each combination of individuals, the function returns total number (or percentage) of detections that are closely associated in time and space. For many combinations of individuals, especially those with long, overlapping timeseries, the function may take some time (minutes to hours) to run; therefore, testing the function on a small subset of individuals first is advisable. However, parallelisation can be used to improve computation time. Similarity matrices can be visualised with \code{\link[plot.pretty]{pretty_mat}}.
#'
#' @param acoustics_ls A list of dataframes, with one element for each individual, which contain each individual's detection timeseries. Each dataframe must include the following columns: 'individual_id', a factor which specifies unique individuals; 'timestamp', a \code{\link[base]{DateTimeClasses}} object which specifies the time of each detection; 'receiver_long', the longitude (decimal degrees) of the receiver(s) at the individual was detected; and 'receiver_lat', the latitude (decimal degrees) of the receiver(s) at which individual was detected. Each dataframe should be ordered by 'individual_id' and then by 'timestamp'. Careful ordering of 'individual_id' factor levels (e.g. perhaps by population group, then by the number of detections of each individual) can aid visualisation of similarity matrices, in which the order or rows/columns corresponds directly to the order of individuals in \code{acoustics_ls}. Sequential elements in \code{acoustics_ls} should correspond to sequential factor levels for 'individual_id', which should be the same across all dataframes.
#' @param thresh_time A number which specifies the time, in minutes, after which detections at nearby receivers are excluded.
#' @param thresh_dist A number which specifies the (Euclidean) distance between receivers, in metres, beyond which detections are excluded (see Details).
#' @param cl (optional) A cluster object created by \code{\link[parallel]{makeCluster}}. This is required if you want to run the algorithm in parallel. The default is \code{NULL} (i.e. the algorithm is run on a single processor: see examples). If supplied, the connection to the cluster is stopped within the function.
#' @param varlist (optional) A character vector of names of objects to export, to be passed to the \code{varlist} argument of \code{\link[parallel]{clusterExport}}. This is required if \code{cl} is supplied and you specify some function arguments via objects, rather than directly. These objects must be located in the global environment.
#' @param output A number which specifies the output type. If \code{output = 1}, the function returns a (usually) symmetric similarity matrix in which each number represents the number of detections nearby in space and time for each pair of individuals. Row and column names are assigned from the 'individual_id' column in \code{acoustics_ls} dataframes. This matrix is usually symmetric, but this is not necessarily the case for data collected from tags which transmit at random intervals around a nominal delay: under this scenario, the tag for a given individual (i) may transmit multiple signals in the space of time that the tag for another individual (j) only releases a single transmission. In this case, the pairing i,j will comprise all unique transmissions for individual i, paired with the nearest observations for individual j, some of which will be duplicate observations. Therefore, this pairing will contain more 'shared observations' than the pairing j,i. However, even under random transmission, the matrix will usually be either symmetric or very nearly symmetric, with only small differences between identical pairs. If \code{output = 2}, the function returns a list with the following elements: (1) 'mat_sim', the symmetric similarity matrix (see above); 'mat_nobs', a matrix with the same dimensions as 'mat_sim' which specifies the number of observations for each individual (by row, used to calculate 'mat_pc', see later); 'mat_pc', a non-symmetric matrix in which each cell represent the percent of observations of the individual in the i'th row that are shared with the individual in the j'th column; and 'dat', a nested list, with one element for each individual which comprises a list of dataframes, one for each other individual, each one of which contains the subset of observations that are shared between the two individuals. Each dataframe contains the same columns as in the \code{acoustics_ls} dataframes with the following columns added: 'pos_in_acc2', 'timestamp_acc2', 'receiver_lat_acc2' and 'receiver_long_acc2', which represent the positions, timestamps and locations of corresponding observations in the second individual's dataframe to the first individual's dataframe, and 'difftime_abs' and 'dist_btw_receivers' which represent the duration (minutes) and distances (m) between corresponding observations. When there are no shared observations between a pair of individuals, the element simply contains \code{NULL}. Note that 'mat_pc' is computed by (mat_sim/mat_nobs)*100. The matrix is therefore non-symmetric (if individuals have differing numbers of observations); i.e., mat_pc[i, j] is the percent of individual i's observations that are shared with individual j; while mat_pc[j, i] is the percent of individual j's observations that are shared with individual i. NaN elements are possible in 'mat_pc' for levels of the factor 'individual_id' without observations.
#' @param verbose A logical input which specifies whether or not to print messages to the console which relay function progress. This is ignored if \code{cl} is supplied.
#'
#' @details
#' \subsection{Background}{
#' Passive acoustic telemetry is widely used to study animal space use, and the possible drivers of spatiotemporal patterns in space use, in aquatic environments. Patterns in space use are widely related to environmental conditions, such as temperature, but the role of interactions among individuals is usually unexplored due to a paucity of data, despite their likely importance. However, discrete detections also contain information on interactions among individuals that may influence space use through their similarities and differences among individuals over time and space. Similarities and differences can take different forms with differing ecological implications. For example, for individuals that are frequently detected in similar areas, detections may indicate (a) prolonged associations among individuals, if detections are usually closely associated in time and space (for example, due to parent-offspring relationships, group-living and/or mating); or (b) avoidance and/or territorial-like behaviour if detections, while close in space, are usually at different receivers and/or disjointed in time. Likewise, detection similarities among individuals that are rarely detected, or usually detected at disparate receivers, may reflect important interactions among those individuals at particular times (e.g. mating). To explore similarities and differences in patterns of space use, visualisation of detection histories with abacus plots and maps is beneficial. However, for many individuals and/or dense/large receiver arrays, quantification of the similarities in detections over time and space is challenging. To this end, \code{compute_det_sim()} computes a similarity matrix across all individuals, defining the number (or percentage) of detections for each individual that are nearby in time, or space, to detections for each other individual.}
#' \subsection{Assumptions}{
#' The distances beyond which detections of different individuals nearby in time are considered to demonstrate that those individuals are not closely associated are Euclidean. This may be problematic (e.g. when receivers hug complex coastlines).}
#'
#' @return The function returns a matrix or a list, depending on the input to \code{output} (see above).
#'
#' @examples
#' #### Prepare data
#' # acoustics_ls requires a dataframe with certain columns
#' # individual_id should be a factor
#' dat_acoustics$individual_id <- factor(dat_acoustics$individual_id)
#' # ensure dataframe ordered by individual, then timestamp
#' dat_acoustics <- dat_acoustics[order(dat_acoustics$individual_id, dat_acoustics$timestamp), ]
#' # define list of dataframes
#' acoustics_ls <- split(dat_acoustics, factor(dat_acoustics$individual_id))
#'
#' #### Example (1): Compute detection similarity matrix using default options
#' # mat_sim contains the number of observations shared among individuals
#' mat_sim <- compute_det_sim(acoustics_ls = acoustics_ls,
#'                            thresh_time = 90,
#'                            thresh_dist = 0)
#' plot.pretty::pretty_mat(mat_sim, col_diag = "dimgrey")
#'
#' #### Example (2): Return list of outputs
#' out_ls <- compute_det_sim(acoustics_ls = acoustics_ls,
#'                           thresh_time = 90,
#'                           thresh_dist = 0,
#'                           output = 2)
#' names(out_ls)
#' # Examine number of observations for each individual
#' plot.pretty::pretty_mat(out_ls$mat_nobs)
#' # Examine % shared detections between individuals
#' plot.pretty::pretty_mat(out_ls$mat_pc, col_diag = "dimgrey")
#'
#' #### Example (3): Turn off messages with verbose = FALSE
#' out_ls_non_verb <- compute_det_sim(acoustics_ls = acoustics_ls,
#'                                    thresh_time = 90,
#'                                    thresh_dist = 0,
#'                                    verbose = FALSE)
#'
#' #### Example (4): Implement algorithm in parallel
#' out_ls_pl <- compute_det_sim(acoustics_ls = acoustics_ls,
#'                              thresh_time = 90,
#'                              thresh_dist = 0,
#'                              cl = parallel::makeCluster(2L),
#'                              output = 2)
#' names(out_ls_pl)
#'
#' @author Edward Lavender
#' @export

compute_det_sim <-
  function(acoustics_ls,
           thresh_time,
           thresh_dist,
           cl = NULL,
           varlist = NULL,
           output = 1,
           verbose = TRUE
  ){

    #### Initial checks
    # Check that acoustics dataframes contains required columns
    mapply(acoustics_ls, 1:length(acoustics_ls), FUN = function(acoustics, i){
      if(!is.null(acoustics)){
        if(any(!(c("individual_id", "timestamp", "receiver_long", "receiver_lat") %in% colnames(acoustics)))){
          stop(paste0("acoustic_ls[[", i, "]] does not contain all required column names."))
        }
      }
    })

    #### Set up
    # Check that acoustics_ls is a factor
    first_non_NULL <- min(which(sapply(acoustics_ls, function(acoustics) return(!is.null(acoustics)))))
    if(!inherits(acoustics_ls[[first_non_NULL]]$individual_id, "factor")) stop("acoustics_ls dataframe must contain individual_id column that is a factor.")
    # Define factor levels (i.e. the names of individuals )
    id_names <- levels(acoustics_ls[[first_non_NULL]]$individual_id)
    # id_names <- as.character(sapply(acoustics_ls, function(acoustics) return(acoustics$individual_id[1])))
    nid <- length(id_names)
    # Check that there is one element in acoustics_ls for every individual
    if(length(acoustics_ls) != nid) stop("'acoustics_ls' needs one element for every individual individual_id level. Add NULL elements to 'acoustics_ls' for remaining individuals.")
    # Create a blank similarity matrix which we'll fill in
    mat_sim <- matrix(NA, nrow = nid, ncol = nid, dimnames = list(id_names, id_names))

    #### Loop over all combinations of individuals, pair timeseries and identify
    # ... nearby observations in time and space
    if(all(!is.null(varlist), !is.null(cl))) parallel::clusterExport(cl = cl, varlist = varlist)
    lout <-
      pbapply::pblapply(acoustics_ls, cl = cl, function(acc1){

        #### Testing:
        # acc1 = acoustics_ls[[1]]; acc2 = acoustics_ls[[2]];
        # acc1 = acoustics_ls[[2]]; acc2 = acoustics_ls[[1]];

        #### For each individual, loop over each other individual...
        lint <-
          lapply(acoustics_ls, function(acc2){

            #### Ignore NULL/empty elements
            if(any(is.null(nrow(acc1)), is.null(nrow(acc2)), nrow(acc1) == 0, nrow(acc2) == 0)) return(NULL)

            if(acc1$individual_id[1] != acc2$individual_id[1]){

              #### Print individual
              if(verbose){
                cat("\n===================================================================================\n")
                cat(paste("Individual (", as.character(acc1$individual_id[1]),
                          ") and individual (", as.character(acc2$individual_id[1]), ").\n"))
              }

              #### Match timeseries based on closest observations in time using pair_ts()
              if(verbose) cat("Matching detection timeseries...\n")
              # Remove any observations from the second individual more than some limit outside of the timeseries of the first individual
              # ... and vice versa (for speed when matching).
              acc2 <- acc2[acc2$timestamp >= (min(acc1$timestamp) - thresh_time*60*2) &
                             acc2$timestamp <= (max(acc1$timestamp) + thresh_time*60*2), ]
              if(nrow(acc2) == 0) return(NULL)
              acc1 <- acc1[acc1$timestamp >= (min(acc2$timestamp) - thresh_time*60*2) &
                           acc1$timestamp <= (max(acc2$timestamp) + thresh_time*60*2), ]
              if(nrow(acc1) == 0) return(NULL)

              # Check for duplicated timestamps in each individuals dataframe, and adjust these
              # ... by a small fraction prior to matching, so that all are included. Unless there are 100,000s
              # ... of duplicate timestamps, this approach does not produce any duplicated observations.
              dup1 <- duplicated(acc1$timestamp)
              dup2 <- duplicated(acc2$timestamp)
              adj <- (thresh_time*60)/4
              if(any(dup1)){
                pos_dups1 <- which(dup1)
                lpd1 <- length(pos_dups1)
                adj_dups1 <- stats::runif(lpd1, -adj, adj)
                acc1$timestamp[pos_dups1] <- acc1$timestamp[pos_dups1] + adj_dups1
                if(any(duplicated(acc1$timestamp))) warning(paste("Duplicate timestamps in, ", as.character(acc1$individual_id[1]), "element in acoustic_ls."))
              }
              if(any(dup2)){
                pos_dups2 <- which(dup2)
                lpd2 <- length(pos_dups2)
                adj_dups2 <- stats::runif(lpd2, -adj, adj)
                acc2$timestamp[pos_dups2] <- acc2$timestamp[pos_dups2] + adj_dups2
                if(any(duplicated(acc2$timestamp))) warning(paste("Duplicate timestamps in, ", as.character(acc2$individual_id[1]), "element in acoustic_ls."))
              }
              # Match timeseries, readjusting any adjusted timestamps back to their original values
              # ... before these are added to the dataframe.
              acc1$pos_in_acc2 <- Tools4ETS::match_ts_nearest(acc1$timestamp, acc2$timestamp)
              if(any(dup1)) acc1$timestamp[pos_dups1] <- acc1$timestamp[pos_dups1] - adj_dups1
              if(any(dup2)) acc2$timestamp[pos_dups2] <- acc2$timestamp[pos_dups2] - adj_dups2
              acc1$timestamp_acc2 <- acc2$timestamp[acc1$pos_in_acc2]

              #### Exclude any timestamps more than timestamp beyond each other (could be 0 mins)
              # Implement this now, before a threshold based on distance, below, for speed.
              if(verbose) cat("Processing timeseries by theshold time difference...\n")
              acc1$difftime_abs <- abs(difftime(acc1$timestamp, acc1$timestamp_acc2, units = "mins"))
              acc1 <- acc1[acc1$difftime_abs <= thresh_time, ]
              if(nrow(acc1) == 0) return(NULL)

              #### Distances between pairs of receivers
              if(verbose) cat("Computing differences between pairs of receivers...\n")
              # Add receivers
              acc1$receiver_lat_acc2 <- acc2$receiver_lat[acc1$pos_in_acc2]
              acc1$receiver_long_acc2 <- acc2$receiver_long[acc1$pos_in_acc2]
              # Compute distances
              acc1$dist_btw_rec <- geosphere::distGeo(acc1[, c("receiver_long", "receiver_lat")],
                                                      acc1[, c("receiver_long_acc2", "receiver_lat_acc2")])

              #### Exclude any receivers more than some threshold distance beyond each other (could be 0 m):
              if(verbose) cat("Processing timeseries by threshold distance...\n")
              acc1 <- acc1[acc1$dist_btw_rec <= thresh_dist, ]
              if(nrow(acc1) == 0) return(NULL) else return(acc1)

            }
          })

        return(lint)

      })
    if(!is.null(cl)) parallel::stopCluster(cl = cl)

    #### Populate similarity matrix
    for(i in 1:nid){
      for(j in 1:nid){
        if(i != j){
          d <- lout[[i]][[j]]
          if(!is.null(d)) mat_sim[i, j] <- nrow(d) else mat_sim[i, j] <- 0
        }
      }
    }

    #### Matrix based on % similarity
    mat_nobs <- mat_sim[]
    for(i in 1:nrow(mat_nobs)){
      if(!is.null(acoustics_ls[[i]])) mat_nobs[i, ] <- nrow(acoustics_ls[[i]]) else mat_nobs[i, ] <- 0
    }
    mat_pc <- (mat_sim/mat_nobs)*100

    #### Return outputs
    if(!(output %in% 1:2)){
      warning(paste("'output ", output, " not supported; defaulting to output = 1."))
      output <- 1
    }
    if(output == 1) {
      out <- mat_sim
    } else if(output == 2) {
      out <- list(mat_sim = mat_sim, mat_nobs = mat_nobs, mat_pc = mat_pc, dat = lout)
    }
    return(out)

  }


#### End of code.
#######################################
#######################################
