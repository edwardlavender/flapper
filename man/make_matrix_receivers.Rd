% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/make_matrices.R
\name{make_matrix_receivers}
\alias{make_matrix_receivers}
\title{Matricise receiver deployment time series}
\usage{
make_matrix_receivers(
  moorings,
  services = NULL,
  start = NULL,
  end = NULL,
  delta_t = "120 mins",
  as_POSIXct = as.POSIXct,
  set_names = TRUE
)
}
\arguments{
\item{moorings}{A dataframe that defines receiver IDs and deployment times. This must contain the following columns: an identifier for receivers (named `receiver_id'), the start time of receiver' deployment periods (`receiver_start_date') and the end time of receivers' deployment periods (`receiver_end_date') (see \code{\link[flapper]{dat_moorings}} for an example). Deployment times can be recorded as Date or POSIXct objects.}

\item{services}{(optional) A dataframe that defines receiver IDs and servicing dates (times during the deployment period of a receiver when it was not active due to servicing). If provided, this must contain the following columns: an identifier for serviced receivers (named `receiver_id') and two columns that define the time of the service(s) (`service_start_date' and `service_end_date'). Times can be recorded as Date or POSIXct objects. Before/after service events, receivers are assumed to have been deployed in the same locations; receiver deployments in different locations before/after servicing should be treated as distinct deployments in \code{moorings}.}

\item{start, end}{Date or POSIXct objects that define the start and end time. If unspecified, these are taken from the range of deployment times in \code{moorings}.}

\item{delta_t}{A number or character that defines the time interval between successive time steps. This is passed to the `by' argument of \code{\link[base]{seq.POSIXt}} or \code{\link[base]{seq.Date}} (depending on \code{as_POSIXct}, below).}

\item{as_POSIXct}{(optional) A function that coerces supplied any supplied times (\code{moorings$receiver_start_date}, \code{moorings$receiver_end_date}, \code{services$service_start_date}, \code{services$service_end_date}, \code{start} and \code{end}) that are not POSIXct objects to POSIXct objects. This can be suppressed via \code{as_POSIXct = NULL} if supplied times are Date objects and \code{delta_t} is not less than one day.}

\item{set_names}{A logical variable that defines whether or not to set the row and column names of the matrix to the time steps and the receiver IDs respectively.}
}
\value{
The function returns a matrix with one row for each time step and one column for each receiver. Each cell defines whether (1) or not (0) each receiver was at active during that time step. A `bins' attribute is included, which defines the time steps as a Date or POSIXct vector.
}
\description{
This function creates a matrix that, for each time step (matrix row) in a sequence of user-defined times, defines whether or not each receiver (matrix column) was active. To implement the function, a dataframe with receiver IDs and deployment start and end times must be supplied (via \code{moorings}). Servicing dates can also be accounted for via a dataframe with receiver IDs and servicing times (\code{services}). The times for which to express whether or not each receiver was active are provided by optionally defining a \code{start} and \code{end} date (these can be taken from the range of deployment times in \code{moorings} if unspecified) and the interval (\code{delta_t}) between time steps.
}
\examples{
#### Example (1): Illustration using fake data

## Define some example 'moorings' data
# ... with receiver IDs and deployment times
moorings <- data.frame(
  receiver_id = c(1, 2, 3, 4, 5),
  receiver_start_date = as.Date(c(
    "2016-01-01",
    "2016-01-02",
    "2016-01-03",
    "2016-01-04",
    "2016-01-05"
  )),
  receiver_end_date = as.Date(c(
    "2016-01-06",
    "2016-01-07",
    "2016-01-08",
    "2016-01-09",
    "2016-01-09"
  ))
)

## Define some example 'servicing' data
# ... with receiver IDs and servicing times
# ... Here, receiver 1 was serviced twice
# ... ... from 2016-01-02--3 and 2016-01-04--5
# ... and receiver 5 was serviced
# ... ... on 2016-01-08.
services <- data.frame(
  receiver_id = c(1, 1, 5),
  service_start_date = as.Date(c(
    "2016-01-02",
    "2016-01-04",
    "2016-01-08"
  )),
  service_end_date = as.Date(c(
    "2016-01-03",
    "2016-01-05",
    "2016-01-08"
  ))
)

## Get daily receiver status (0, 1) matrix
make_matrix_receivers(moorings, delta_t = "days", as_POSIXct = NULL)

## Get daily receiver status (0, 1) matrix
# ... accounting for servicing dates
make_matrix_receivers(moorings, services, delta_t = "days", as_POSIXct = NULL)

#### Example (2): Illustration using actual data
# ... for different time windows
mat_days <- make_matrix_receivers(dat_moorings, delta_t = "days", as_POSIXct = NULL)
mat_hours <- make_matrix_receivers(dat_moorings, delta_t = "hours")
utils::str(mat_days)
utils::str(mat_hours)

}
\author{
Edward Lavender
}
