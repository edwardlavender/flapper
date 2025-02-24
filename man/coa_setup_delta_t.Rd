% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/coa.R
\name{coa_setup_delta_t}
\alias{coa_setup_delta_t}
\title{Suggest time intervals over which to calculate centres of activity}
\usage{
coa_setup_delta_t(
  acoustics,
  delta_t,
  method = 1L,
  implementation = 1L,
  moorings = NULL,
  xlab = NULL,
  ylab = NULL,
  main = NULL,
  xlim = NULL,
  ylim = NULL,
  add_additional = NULL,
  type = "b",
  ...,
  verbose = TRUE
)
}
\arguments{
\item{acoustics}{A dataframe with passive acoustic telemetry detection time series (see \code{\link[flapper]{dat_acoustics}} for an example) for a specific individual. This should contain the following columns: a vector of receiver IDs, named `receiver_id'; and a POSIXct vector of time stamps, named `timestamp'.}

\item{delta_t}{A number or character, passed to the `breaks' argument of \code{\link[base]{cut.POSIXt}}, that defines the time interval to be evaluated.}

\item{method}{An integer (\code{1L}, \code{2L} or \code{1:2L}) that defines the method(s) used to evaluate the time intervals (see Description).}

\item{implementation}{An integer (\code{1L} or \code{2L}) that defines the metric(s) used in the evaluation: \code{1}, the total number of receivers/detections; or \code{2}, the percent of receivers/number of detections per receiver (see Description). If \code{implementation = 2L} and the total number of receivers is constant through time, the function reverts to \code{implementation = 1L}, the results of which are more interpretable.}

\item{moorings}{(optional) If \code{implementation = 2L}, \code{moorings} is a dataframe that defines, for each receiver, the deployment period. This is required to calculate the number of operational receivers over time (via \code{\link[flapper]{get_n_operational_ts}}). This must contain the following columns: a vector of receiver IDs that includes all receivers in \code{acoustics}, named `receiver_id'; and the times of receiver deployment and retrieval, named `receiver_start_date' and `receiver_end_date' respectively. These are coerced to POSIXct vectors, if required, to match \code{acoustics$timestamp}.}

\item{xlab, ylab, main}{Character strings that define the x and y axis labels and the plot title. If \code{method = 1:2L}, a single input will label both plots identically, whereas a list with one element for each method will label the plots differently.}

\item{xlim, ylim}{X and y axis limits. If \code{method = 1:2L}, a single vector will affect both plots identically, whereas a list with one element for each method will affect plots differently.}

\item{add_additional}{(optional) A stand-alone function, to be executed after a plot has been made, to customise the result.}

\item{type, ...}{Additional arguments, passed to \code{\link[prettyGraphics]{pretty_plot}}, that affect all plots.}

\item{verbose}{A logical variable that defines whether or not to print messages to the console to relay function progress.}
}
\value{
The function returns a named list and a plot for each \code{method}. The list contains a `data' element that is a named list of dataframes for each method (`m1', `m2') and implementation (`i1', `i2'). For `m1_i1', the dataframe contains the number of receivers with detections (`n_receiver_with_detections') and the percent of the time series for which that number of receivers made detections (`pc_of_ts'). For `m1_i2', the percentage of receivers with detections (`pc_receiver_with_detections') is given instead. For `m2_i1', the number of detections (`n_detections') is given; for `m2_i2', it is the number of detections per receiver (`n_detections_per_receiver'). The list also contains an `args' element that records the inputs to \code{acoustics} and \code{delta_t} for reference.
}
\description{
This function implements the two approaches introduced by Simpfendorfer et al. (2002) to suggest suitable time interval(s) over which to estimate centres of activity (COAs). Both approaches rely on prior specification of a candidate time interval (\code{delta_t}) and provide a qualitative assessment of its suitability. The manual iteration of this scheme over multiple intervals provides a measure of their relative suitability.

The first method (\code{method = 1L}) is based on the premise that estimated COAs are more accurate when they are estimated from detections at more receivers and thus examines the frequency distribution (over the whole time series) of the number of receivers at which detections were made within a given time interval (e.g., one hour). The second method (\code{method = 2L}) is based on the premise that estimated COAs are more accurate when they are estimated from a larger number of detections and thus examines the frequency distribution (over the whole time series) of the number of detections in a given time interval. In a comparison of multiple candidate intervals, the most appropriate interval is the one which optimises a trade off between sufficiently meeting these criteria while remaining sufficiently small in duration such that COAs remain a meaningful representation of the individual's location or short-term centre of activity.

To implement these approaches, it is necessary to supply a dataframe with detections for a particular individual (\code{acoustics}), the interval to be evaluated (\code{delta_t}) and the \code{method}. For each method, two implementations are possible. The first implementation (\code{implementation = 1L}) follows Simpfendorfer et al.'s (2002) approach and examines the frequency distribution(s) for the total number of receivers and/or detections, which is appropriate if the total number of receivers is constant over the duration of detections. The second implementation (\code{implementation = 2L}) examines the distribution(s) of the percent of receivers with detections or the number of detections per receiver and can be more appropriate if the number of receivers changes substantially through time. Using these inputs, the function plots the specified distributions and returns a named list with the data plotted. This approach can be applied iteratively over multiple \code{delta_t} values to evaluate their relative suitability.
}
\details{
The selection of the most appropriate time interval over which to calculate COAs reflects a trade off between (a) larger intervals, over which time more detections at more receivers are made, which enable more accurate estimates, versus (b) smaller intervals, over which time movement is more restricted and centres of activity are representative of the individual's location or short-term centre of activity (Simpfendorfer et al., 2002). As a starting point, Simpfendorfer et al. (2002) suggest intervals from 5 - 60 minutes may often be suitable, although longer intervals may be required.
}
\examples{
#### Example (1): For a specified delta_t, use a specific method and implementation
pp <- graphics::par(mfrow = c(2, 2))
dat <- coa_setup_delta_t(
  acoustics = dat_acoustics[dat_acoustics$individual_id == 25, ],
  delta_t = "6 hours",
  method = 1L,
  implementation = 1L,
  main = "1(1)"
)
utils::str(dat)
dat <- coa_setup_delta_t(
  acoustics = dat_acoustics[dat_acoustics$individual_id == 25, ],
  delta_t = "6 hours",
  method = 1L,
  implementation = 2L,
  moorings = dat_moorings,
  main = "1(2)"
)
utils::str(dat)
dat <- coa_setup_delta_t(
  acoustics = dat_acoustics[dat_acoustics$individual_id == 25, ],
  delta_t = "6 hours",
  method = 2L,
  implementation = 1L,
  main = "2(1)"
)
utils::str(dat)
dat <- coa_setup_delta_t(
  acoustics = dat_acoustics[dat_acoustics$individual_id == 25, ],
  delta_t = "6 hours",
  method = 2L,
  implementation = 2L,
  moorings = dat_moorings,
  main = "2(2)"
)
utils::str(dat)
graphics::par(pp)

#### Example (2) For a specified delta_t, use both methods
pp <- graphics::par(mfrow = c(1, 2))
dat <- coa_setup_delta_t(
  acoustics = dat_acoustics[dat_acoustics$individual_id == 25, ],
  delta_t = "6 hours",
  method = 1:2L,
  implementation = 1L
)
graphics::par(pp)

#### Example (3) Plot customisation options
# xlim,ylim,xlab,ylab,main accept a vector that affects all plots
# ... or a list that affects each plot (see 'ylim' versus 'xlim' below).
# ... other arguments can be passed to prettyGraphics::pretty_plot() via...
pp <- graphics::par(mfrow = c(1, 2))
dat <- coa_setup_delta_t(
  acoustics = dat_acoustics[dat_acoustics$individual_id == 25, ],
  delta_t = "6 hours",
  xlim = list(c(0, 50), c(0, 450)), ylim = c(0, 100),
  method = 1:2L,
  implementation = 1L,
  main = list("A", "B"),
  col = "royalblue", lwd = 2
)
graphics::par(pp)

#### Example (4) Compare multiple delta_t values
pp <- graphics::par(mfrow = c(3, 2))
delta_t_opts <- c("6 hours", "12 hours", "24 hours")
lapply(delta_t_opts, function(delta_t) {
  coa_setup_delta_t(
    acoustics = dat_acoustics[dat_acoustics$individual_id == 25, ],
    delta_t = delta_t,
    method = 1:2L,
    implementation = 2L,
    moorings = dat_moorings,
    main = delta_t, col = "royalblue"
  )
})
graphics::par(pp)

}
\references{
Simpfendorfer, C. A., M. R. Heupel, and R. E. Hueter. 2002. Estimation of short-term centers of activity from an array of omnidirectional hydrophones and its use in studying animal movements. Canadian Journal of Fisheries and Aquatic Sciences 59:23-32.
}
\seealso{
\code{\link[flapper]{coa}} calculates centres of activity.
}
\author{
Edward Lavender
}
