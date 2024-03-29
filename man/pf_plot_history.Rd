% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pf_analyse_archive.R
\name{pf_plot_history}
\alias{pf_plot_history}
\title{Plot particle histories from a PF algorithm}
\usage{
pf_plot_history(
  archive,
  time_steps = 1:length(history),
  add_surface = list(),
  add_particles = list(pch = "."),
  forwards = TRUE,
  prompt = TRUE,
  ...
)
}
\arguments{
\item{archive}{A \code{\link[flapper]{pf_archive-class}} object from \code{\link[flapper]{pf}}, or \code{\link[flapper]{pf}} plus \code{\link[flapper]{pf_simplify}} with the \code{return = "archive"} argument, that contains particle histories.}

\item{time_steps}{An integer vector that defines the time steps for which to plot particle histories.}

\item{add_surface}{A named list, passed to \code{\link[prettyGraphics]{pretty_map}}, to customise the appearance of the surface, which shows the set of possible positions that the individual could have occupied at a given time step (from \code{\link[flapper]{ac}}, \code{\link[flapper]{dc}} and \code{\link[flapper]{acdc}}), on each map.}

\item{add_particles}{A named list, passed to \code{\link[prettyGraphics]{pretty_map}}, to customise the appearance of the particles on each map.}

\item{forwards}{A logical variable that defines whether or not create plots forwards (i.e., from the first to the last \code{time_steps}) or backwards (i.e., from the last to the first \code{time_steps}).}

\item{prompt}{A logical input that defines whether or not to pause between plots (\code{prompt = TRUE}).}

\item{...}{Plot customisation arguments passed to \code{\link[prettyGraphics]{pretty_map}}.}
}
\value{
The function returns a plot, for each time step, of all the possible locations of the individual, with sampled locations overlaid.
}
\description{
This function plots the spatiotemporal particle histories from a particle filtering (PF) algorithm (the acoustic-container PF, the depth-contour PF or the acoustic-container depth-contour PF). This produces, for each time step, a map of the individual's possible locations (from the AC, DC or ACDC algorithm), with sampled locations (derived via the particle filtering routine) overlaid.
}
\examples{
#### Implement pf() algorithm
# Here, we use pre-defined outputs for speed

#### Example (1): The default implementation
pf_plot_history(dat_dcpf_histories, time_steps = 1)

#### Example (2): Plot customisation options, e.g.:
# Customise bathy via add_bathy()
pf_plot_history(dat_dcpf_histories,
  time_steps = 1,
  add_surface = list(col = c(grDevices::topo.colors(2)))
)
# Customise particles via add_particles
pf_plot_history(dat_dcpf_histories,
  time_steps = 1,
  add_particles = list(col = "red")
)
# Pass other arguments to prettyGraphics::pretty_map() via ...
pf_plot_history(dat_dcpf_histories,
  time_steps = 1,
  add_polys = list(x = dat_coast, col = "brown"),
  crop_spatial = TRUE
)

#### Example (3): Plot multiple time steps
pp <- graphics::par(mfrow = c(2, 2))
pf_plot_history(dat_dcpf_histories, time_steps = 1:4, prompt = FALSE)
graphics::par(pp)

#### Example (4): Compare outputs for sampled versus connected particles
dat_dcpf_histories_connected <-
  pf_simplify(dat_dcpf_histories, return = "archive")
pp <- graphics::par(mfcol = c(2, 4))
pf_plot_history(dat_dcpf_histories,
  time_steps = 1:4,
  add_particles = list(pch = 21, bg = "black"),
  prompt = FALSE
)
pf_plot_history(dat_dcpf_histories_connected,
  time_steps = 1:4,
  add_particles = list(pch = 21, bg = "black"),
  prompt = FALSE
)
graphics::par(pp)

}
\seealso{
\code{\link[flapper]{pf}} implements PF. \code{\link[flapper]{pf_simplify}} assembles paths from particle histories. \code{\link[flapper]{pf_plot_map}} creates an overall `probability of use' map from particle histories. \code{\link[flapper]{pf_plot_1d}}, \code{\link[flapper]{pf_plot_2d}} and \code{\link[flapper]{pf_plot_3d}} provide plotting routines for paths. \code{\link[flapper]{pf_loglik}} calculates the log-probability of each path.

\code{\link[flapper]{pf}} implements PF. \code{\link[flapper]{pf_simplify}} assembles paths from particle histories. \code{\link[flapper]{pf_plot_history}} visualises particle histories. \code{\link[flapper]{pf_plot_1d}}, \code{\link[flapper]{pf_plot_2d}} and \code{\link[flapper]{pf_plot_3d}} provide plotting routines for paths. \code{\link[flapper]{pf_loglik}} calculates the log-probability of each path.
}
\author{
Edward Lavender
}
