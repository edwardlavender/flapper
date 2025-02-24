#' @title \code{flapper}: Routines for the analysis of passive acoustic telemetry data
#' @author Edward Lavender
#'
#' @description \code{flapper} is an R package that provides routines for the analysis of passive acoustic telemetry data, including for data processing, description (via detection statistics and movement metrics), modelling and simulation. These routines are supported by a series of helper functions (including distance calculators, spatial tools and parallelisation options). The main contribution is the provision of a family of algorithms for movement modelling in passive acoustic telemetry systems that permits the reconstruction of fine-scale movement paths and emergent patterns of space use. Package development has been motivated by passive acoustic telemetry data from a Critically Endangered benthic elasmobranch (the flapper skate, \emph{Dipturus intermedius}) off the west coast of Scotland.
#'
#' @section Data processing:
#' Some functions facilitate the acquisition, assembly, processing and checking of passive acoustic telemetry time series:
#' \itemize{
#'   \item Data acquisition
#'   \itemize{
#'      \item \link{query_open_topo} queries the Topo Data Application Programming Interface for elevation/bathymetry data;
#'      }
#'   \item Data assembly
#'    \itemize{
#'      \item \link{assemble_sentinel_counts} assembles counts of transmissions/detections from sentinel tags for modelling purposes (i.e., to model detection probability);
#'      \item \link{make_matrix_ids} matricises individual deployment time series;
#'      \item \link{make_matrix_receivers} matricises receiver deployment time series;
#'      \item \link{make_matrix_detections} matricises detection time series;
#'      \item \link{make_df_detections} reverses this process;
#'   }
#'   \item Data processing
#'   \itemize{
#'     \item \link{process_receiver_id} adds unique receiver IDs to a dataframe (e.g., if the same receiver has been deployed more than once);
#'     \item \link{process_false_detections_sf} passes putative false detections through a spatial filter which incorporates ancillary information on receiver locations and animal swimming speeds to interrogate their plausibility;
#'     \item \link{process_quality_check} passes acoustic data through some basic quality checks prior to analysis;
#'     \item \link{process_surface} determines an 'optimum' \code{\link[raster]{raster}} aggregation method and error induced by this process;
#'    }
#' }
#'
#' @section Spatial tools:
#' Some functions facilitate spatial operations that support common tasks and modelling algorithms:
#' \itemize{
#'   \item \link{buffer_and_crop} buffers a spatial object (e.g., receiver locations) and uses this buffered object to crop another (e.g., the local bathymetry);
#'   \item \link{get_intersection} intersects spatial geometries;
#'   \item \link{xy_from_click} gets location coordinates from mouse clicks;
#'   \item \link{crop_from_click} crops a \code{\link[raster]{raster}} to an area defined by mouse clicks;
#'   \item \link{cells_from_val} returns the cells (or a \code{\link[raster]{raster}} of the cells) of a \code{\link[raster]{raster}} that are equal to a specified value or lie within a specified range of values;
#'   \item \link{invert_poly} inverts a polygon (e.g, to define the `sea' from a polygon of the `land');
#'   \item \link{mask_io} masks values in a \code{\link[raster]{raster}} that lie inside or outside of a spatial mask (e.g., to mask the `land' from the `sea');
#'   \item \link{sim_surface} populates a \code{\link[raster]{raster}} with simulated values;
#'   \item \link{split_raster_equally} splits a \code{\link[raster]{raster}} into equal pieces (using code from the greenbrown (http://greenbrown.r-forge.r-project.org) package);
#'   \item \link{update_extent} shrinks or inflates an extent object;
#'   \item \link{segments_cross_barrier} determines if Euclidean transects cross a barrier;
#' }
#'
#' @section Distance calculations:
#' Some functions facilitate distance calculations, including the calculation of distances between receivers, along 3-dimensional movement paths, and of the shortest paths over a surface.
#' \itemize{
#'   \item Euclidean distances
#'   \itemize{
#'     \item \link{dist_btw_clicks} calculates distances and draws segments between sequential mouse clicks on a map;
#'     \item \link{dist_btw_receivers} calculates the Euclidean distances between all combinations of receivers;
#'     \item \link{dist_btw_points_3d} calculates the Euclidean distances between points in three-dimensional space;
#'     \item \link{dist_over_surface} calculates the Euclidean distance along a path over a three-dimensional surface;
#'   }
#'   \item Shortest (least-cost) distances
#'   \itemize{
#'     \item \link{lcp_costs} calculates the distances between connected cells in a \code{\link[raster]{raster}}, accounting for planar (x, y, diagonal) and vertical (z) distances;
#'     \item \link{lcp_graph_surface} constructs connected graphs for least-cost paths analysis;
#'     \item \link{lcp_from_point} calculates least-cost distances from a point on a \code{\link[raster]{raster}} to all of the other cells of a \code{\link[raster]{raster}};
#'     \item \link{lcp_over_surface} calculates the shortest path(s) and/or the distances of the shortest path(s) over a surface between origin and destination coordinates;
#'     \item \link{lcp_interp} interpolates paths between sequential locations using least-cost paths analysis;
#'     \item \link{lcp_comp} compares Euclidean and shortest distance metrics for an area;
#'   }
#' }
#'
#' @section Detection statistics:
#' Some functions facilitate the calculation of detection statistics, including those related to sampling effort and to detections of individuals:
#' \itemize{
#'   \item \link{get_detection_pr} calculates detection probability given a model for detection probability with distance;
#'   \item \link{get_detection_containers} defines detection containers (areas within the maximum detection range) around receivers;
#'   \item \link{get_detection_containers_overlap} identifies receivers with overlapping detection containers in space and time;
#'   \item \link{get_detection_containers_envir} extracts environmental conditions from within receiver detection ranges, accounting for detection probability;
#'   \item \link{get_detection_area_sum} calculates the total area surveyed by receivers;
#'   \item \link{get_detection_area_ts} defines a time series of the area surveyed by receivers;
#'   \item \link{get_n_operational_ts} defines a time series of the number of operational units (e.g., individuals at liberty or active receivers);
#'   \item \link{get_id_rec_overlap} calculates the overlap between the deployment periods of tagged individuals and receivers;
#'   \item \link{get_detection_days} calculates the total number of days on which each individual was detected (termed `detection days;);
#'   \item \link{get_detection_clumps} identifies detection `clumps' and calculates their lengths;
#'   \item \link{get_detection_overlaps} identifies `overlapping' detections;
#'   \item \link{get_residents} identifies `resident' individuals;
#'   \item \link{make_matrix_cooccurence} computes a detection history similarity matrix across individuals;
#' }
#'
#' @section Movement metrics:
#' Building on the analysis of detection time series, some functions provide movement metrics:
#' \itemize{
#'   \item \link{get_mvt_mobility} functions estimate swimming speeds:
#'   \itemize{
#'     \item \link{get_mvt_mobility_from_acoustics} estimates swimming speeds from acoustic detections;
#'     \item \link{get_mvt_mobility_from_archival} estimates swimming speeds from archival time series;
#'   }
#'   \item \link{get_mvt_resting} identifies `resting' behaviour from archival time series;
#'   \itemize{
#'     \item \link{get_hr} functions get animal `home ranges':
#'     \itemize{
#'       \item \link{get_hr_prop} gets a custom range from a utilisation distribution (UD);
#'       \item \link{get_hr_core} gets the `core range' from a UD;
#'       \item \link{get_hr_home} gets the `home range' from a UD;
#'       \item \link{get_hr_full} gets the `full range' from a UD;
#'     }
#'   }
#' }
#'
#' @section Modelling algorithms:
#' The main thrust of flapper is the implementation of existing and new algorithms designed to reconstruct fine-scale movement paths and emergent patterns of space use in passive acoustic telemetry systems:
#' \itemize{
#'   \item The centres of activity (COA) algorithm
#'   \itemize{
#'     \item \link{coa_setup_delta_t} informs decisions as to an appropriate time interval over which to calculate COAs;
#'     \item \link{make_matrix_detections} summarises matrices over time intervals (see above);
#'     \item \link{coa} implements the arithmetic version of the mean-position algorithm to calculate COAs;
#'     \item \link{kud_habitat}, \link{kud_around_coastline} and \link{kud_around_coastline_fast} facilitate the estimation of home ranges (e.g., from estimated COAs) in areas of complex coastline;
#'   }
#'   \item The `flapper` family of algorithms
#'   \itemize{
#'     \item The AC/DC branch
#'     \itemize{
#'       \item The depth-contour (DC) algorithm
#'       \itemize{
#'          \item \link{dc} implements the DC algorithm;
#'          \item \link{dcq} implements the quick DC (DCQ) algorithm;
#'          }
#'      \item The acoustic-container* (AC*) algorithms
#'      \itemize{
#'         \item \link{acs_setup_mobility} examines the assumption of a constant `mobility' parameter;
#'         \item \link{acs_setup_containers} defines the detection containers for the algorithm(s);
#'         \item \link{acs_setup_detection_kernels} defines detection probability kernels for the algorithm(s);
#'         \item \link{ac} and \link{acdc} implement the acoustic-container (AC) and acoustic-container depth-contour (ACDC) algorithms, via \link{.acs_pl} and \link{.acs};
#'          }
#'      \item AC/DC processing
#'      \itemize{
#'         \item \link{acdc_simplify} simplifies \link{acdc_archive-class} objects into \link{acdc_record-class} objects;
#'         \item \link{acdc_access} functions provide short-cuts to different elements of \link{acdc_record-class} objects:
#'         \item \link{acdc_access_dat} accesses stored dataframes in an \link{acdc_record-class} object object;
#'         \item \link{acdc_access_timesteps} accesses the total number of time steps in an \link{acdc_record-class} object;
#'         \item \link{acdc_access_maps} accesses stored maps in an \link{acdc_record-class} object;
#'         \item \link{acdc_plot_trace} plots acoustic container dynamics;
#'         \item \link{acdc_plot_record} plots the results of the algorithm(s);
#'         \item \link{acdc_animate_record} creates html animations of the algorithm(s);
#'         }
#'      }
#'     \item The particle filtering branch
#'     \itemize{
#'       \item \link{pf_setup_movement_pr} provides a simple movement model that defines the probability of movement between locations given the distance between them;
#'       \item \link{pf_setup_record} creates an ordered list of input files;
#'       \item \link{pf_setup_optimisers} controls optimisation settings;
#'       \item \link{pf} implements the PF routine, building on the AC, DC and ACDC algorithms to form the ACPF, DCPF and ACDCPF algorithms;
#'       \item \link{pf_access_history_files} lists particle histories saved to file;
#'       \item \link{pf_access_history} accesses particle histories;
#'       \item \link{pf_access_particles_unique} accesses unique particle samples;
#'       \item \link{pf_plot_history} plots particle histories;
#'       \item \link{pf_animate_history} animates particle histories;
#'       \item \link{pf_simplify} assembles movement paths from particle histories;
#'       \item \link{pf_plot_map} maps the `proportion-of =use' across an area based on sampled particles or reconstructed paths;
#'       \item \link{pf_kud} smooths POU maps using kernel smoothing;
#'       \item \link{pf_kud_1} and \link{pf_kud_2} apply kernel smoothing to sampled particles or reconstructed paths;
#'       \item \link{pf_loglik} calculates the log-probability of reconstructed paths, given the movement model;
#'       \item \link{pf_plot_1d} plots the depth time series from observed and reconstructed paths;
#'       \item \link{pf_plot_2d} maps the reconstructed paths in two-dimensions;
#'       \item \link{pf_plot_3d} maps the reconstructed paths in three-dimensions;
#'       }
#'   }
#' }
#'
#' @section Simulations:
#' A set of \code{sim_*()} functions provide an integrated workflow for the simulation of receiver arrays, movement paths and detections under different array designs, movement models and detection models:
#' \itemize{
#'   \item \link{sim_array} simulates alternative array designs;
#'   \item \link{sim_path_*} functions simulate discrete-time movement paths, including:
#'   \itemize{
#'     \item \link{sim_path_sa}, supported by \link{sim_steps} and \link{sim_angles}, simulates movement paths (possibly in restricted areas) from step lengths and turning angles;
#'     \item \link{sim_path_ou_1} simulates movement paths under Ornstein-Uhlenbeck processes;
#'   }
#'   \item \link{sim_detections} simulates detections at receivers arising from movement paths under a diversity of detection probability models;
#' }
#'
#' Another set of functions facilitate the evaluation of the performance of alternative algorithms for inferring patterns of space use under different array designs, movement models and detections models:
#' \itemize{
#'   \item \link{eval_by_kud} compares patterns of space use inferred from simulated and estimated movement paths using kernel utilisation distributions;
#' }
#'
#' @section Parallelisation:
#' Parallelisation is facilitated by the \code{cl_*()} function family:
#' \itemize{
#'   \item \link{cl_lapply} is a wrapper for \code{\link[pbapply]{pblapply}} that handles cluster checking set up and closure (see \link{flapper-tips-parallel}) using the following functions:
#'   \itemize{
#'     \item \link{cl_check} checks a cluster;
#'     \item \link{cl_cores} identifies the number of cores;
#'     \item \link{cl_chunks} defines chunks for parallelisation;
#'     \item \link{cl_export} exports objects required by a cluster;
#'     \item \link{cl_stop} closes a cluster;
#'   }
#' }
#'
#' @seealso https://github.com/edwardlavender/flapper
#'
#' @docType package
#' @name flapper
NULL
