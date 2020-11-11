
# flapper

***From passive acoustic telemetry to space use: an `R` package for
integrating passive acoustic telemetry and archival tag data to estimate
demersal movement pathways at high resolution***

[![Project Status: Active – The project has reached a stable, usable
state and is being actively
developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)

`flapper` is an `R` package which provides tools for passive acoustic
telemetry (PAT) data. This package has been particularly motivated by
the collection of new acoustic and archival data from a Critically
Endangered elasmobranch, the flapper skate (*Dipturus intermedius*), off
the west coast of Scotland where a static PAT array has been established
to examine the movements of individuals within a Marine Protected Area.
`flapper` has been designed to existing packages for the analyses of
these data (e.g. [Vtrack](https://github.com),
[glatos](https://github.com) and [fishtrack3d](https://github.com)),
with a particular focus on the provision of tools that integrate PAT and
archival data for improved inferrence of patterns of space use,
especially for demersal/benthic species. To this end, `flapper` contains
functions in the following themes:

  - **Data processing tools**, including exploring false detections,
    assembling detection – transmission (i.e., range testing) datasets
    and implementing quality checks;
  - **Spatial tools**, including standard spatial operations, distance
    calculations and the (rapid) computation of least-cost pathways over
    a surface between origin and destination coordinates;
  - **Space use algorithms**, including a more flexible implementation
    of the most widely used approach for inferring space use (i.e. the
    centre of activity – kernel utilisation (COA-KUD) approach) and new
    algorithms designed for demersal/benthic species which integrate PAT
    and archival data to improve estimates of space use;
  - **Simulation tools**, including tools for the comparison of
    simulated and inferred patterns of space use for existing and new
    algorithms;

<img src="vignettes/readme_context.png"/> *flapper: An `R` package
designed to facilitate the integration of acoustic and archival datasets
to improve estimates of space use for demersal and benthic species.
Inserted sample depth and acoustic timeseries were collected by the
Movement Ecology of Flapper Skate project. The insert of the flapper
skate is also courtesy of this project. The bathymetry data are sourced
from the Ireland, Northern Island and Scotland Hydrographic survey (Howe
et al., 2015). Plots were produced using the
[pretttyGraphics](https://github.com/edwardlavender/pretttyGraphics)
package.*

## Installation

You can install the development version of `flapper` from
[GitHub](https://github.com) as shown below:

``` r
devtools::install_github("edwardlavender/flapper", build_vignette = TRUE)
```

To access the vignettes, use `vignette("flapper_intro", package =
"flapper")` for a general introduction to the package.

## Example datasets

`flapper` functions are illustrated using simulated data and the
following sample data collected from flapper skate off the west coast of
Scotland by Marine Scotland Science and Scottish Natural Heritage (now
NatureScot) as part of the Movement Ecology of Flapper Skate project:

  - **dat\_ids** is a dataset containing the characteristics of a sample
    of tagged flapper skate;
  - **dat\_moorings** is a dataset containing some sample passive
    acoustic telemetry receiver locations and associated information;
  - **dat\_acoustics** is a dataset containing some sample detection
    timeseries;
  - **dat\_archival** is a dataset containing some sample depth
    timeseries;
  - **dat\_sentinel** is a dataset containing some sample transmission –
    detection timeseries assembled from sentinel tags;

## Data processing tools

A number of functions facilitate the assembly, processing and checking
of passive acoustic telemetry timeseries.

  - **Data assembly.** Some functions facilitate data assembly. For
    example, `add_receiver_id()` adds unique receiver IDs to a dataframe
    (which is useful if the same receiver has been deployed more than
    once). `assemble_sentinel_counts()` assembles counts of
    transmissions/detections from sentinel tags for modelling purposes
    (i.e., to model detection probability).
  - **Data processing.** Some functions facilitate data processing. The
    identification of false detections is one key component of this
    process. These can be suggested by the ‘short interval’ criterion
    approach using functions in the `glatos`. In `flapper`,
    `false_detections_sf()` passes putative false detections through a
    spatial filter which incorporates ancillary information on receiver
    locations and animal swimming speeds to interrogate their
    plausibility.
  - **Data checking.** Some functions, namely `quality_check()`, pass
    acoustic data through some basic quality checks prior to analysis.

## Spatial tools

A number of other functions faciliate spatial operations, which support
ecological inferences and space use algorithms.

  - **Spatial manipulations.** Some functions provide simple wrappers
    for common spatial operations. For instance, `buffer_and_crop()`
    buffers a spatial object (e.g., receiver locations) and uses this
    buffered object to crop another (e.g., the local bathymetry).
    `update_extent()` shrinks or inflates an extent object.
  - **Euclidean distances.** Some functions facilitate distance
    calculations. For instance, `dist_btw_receivers()` calculates
    Euclidean distances between all receiver combinations and
    `pythagoras_3d()` and `dist_over_surface()` calculate the distances
    between points or along a path over a three-dimensional surface.
  - **Shortest pathways/distances.** Often, Euclidean distances may not
    be a suitable representation of distance. This is especially the
    case for coastal demersal/benthic species in bathymetrically complex
    environments, for which navigation between locations may require
    movement over hilly terrain and around coastline. Thus,
    `lcp_over_surface()` calculates shortest pathways and/or the
    distances of the shortest pathways over a surface between origin and
    destination coordinates.

## Space use algorithms

Some functions facilitate the implementation of existing algorithms
designed to infer space use from PAT data. These include:

  - **The centres of activity – kernel utilisation distribution
    (COA-KUD) algorithm.** This is the most widely used approach for
    inferring animal space use from PAT data. This can be implemented by
    `VTrack`, but `flapper` provides a more flexible implementation
    (`coa()`) which is not so restrictive in terms of data format.

However, the central objective of `flapper` is to implement new methods
designed to provide improved estimates of space use for PAT/archival
data for demersal/benthic species. These include:

  - **The depth contour (DC) algorithm.** This approach is the simplest.
    Under the assumption that individuals are benthic/demersal, this
    algorithm uses observed depths (± some error) to define the subset
    of possible locations of each individual within an defined area.
    This is implemented via `dc()`.
  - **The acoustic centroid – depth contour (ACDC) algorithm.** The ACDC
    algorithm extends the DC algorithm by using passive acoustic
    telemetry data to inform the area within which depth contours are
    most likely to be found. This is implemented via `setup_acdc()` and
    `acdc()`.
  - **The acoustic centroid – depth contour – movement pathway (MP)
    algorithm.** The ACDC-MP algorithm incorporates movement pathways
    into the ACDC process to restrict further the inferred distribution
    of locations within which the individual must have been located at
    each time point. This is implemented with `setup_acdcmp()`,
    `acdcmp()`, and `proc_acdcmp()`.

Some functions provide methods to examine the drivers of patterns of
space use, such as:

  - **Detection similarity matrices.** These are matrices of the total
    number (or percentage) of detections of each individual that are
    ‘nearby’ in space and time to other individuals and can help to
    elucidate possible interactions among individuals that affect space
    use. `compute_det_sim()` computes these matrices.

## Simulation tools

Simulations are a valuable tool in ecology which can elucidate the
relative performance of alternative methods for ecological inferences
(e.g., the `COA` approach versus the `DC` approach for inferring
patterns of space use) and the extent to which new data sources
influence ecological inferences under different circumstances (e.g. the
extent to which sparse or regular PAT detections improve estimates of
space use). To this end, `flapper` provides functions which implement
simulation to evaluate the performance of existing or new approaches for
the processing and analysis of PAT data. This includes the following:

  - **Simulate movement and detections**. Several functions are designed
    to simulate PAT arrays, movement within these arrays and detections
    arising from simulated movements. These include `sim_array()`,
    `sim_beta_surface()` and `sim_det()`. Simulated movement pathways
    and detections can then be used to evaluate the performance of
    processing and space use algorithms, including via the following
    functions.

  - **Simulate false detections.** \[to be written\]

  - **Evaluate the COA (-KUD) approach.** `coa()` and `plot_coa_kud()`
    can be used to calculate COAs/KUDs from simulated movement pathways
    and compare these to simulated patterns.

  - **Evaluate the benefits of PAT data.** \[to be written\]

## Associated packages

  - **[prettyGraphics](https://github.com/edwardlavender/prettyGraphics)**
    facilitates the production of pretty, publication-quality and
    interactive visualisations, with a particular focus on timeseries.
    This makes it easy to create abacus plots, visualise timeseries
    (across factor levels, at different temporal scales and in relation
    to covariates), bathymetric landscapes and movement pathways in
    three-dimensions and detection similarity matrices.
  - **[Tools4ETS](https://github.com/edwardlavender/Tools4ETS)**
    provides a set of general tools for ecological timeseries, including
    the definition of time categories, matching timeseries (e.g.,
    detection observations with environmental covariates), flagging
    independent timeseries and simulating timeseries.
  - **[fvcom.tbx](https://github.com/edwardlavender/fvcom.tbx)**
    provides tools for the integration of hydrodynamic model predictions
    (from the Finite Coastal Ocean Volume Model) with ecological
    datasets (e.g., detection timeseries). This facilitates the
    inclusion of hydrodynamic model predictions as covariates in
    movement models and the validation of hydrodynamic model predictions
    with movement datasets or data collected from static acoustic
    receivers.

## References

Howe, J.A., Anderton, R., Arosio, R., Dove, D., Bradwell, T., Crump, P.,
Cooper, R., Cocuccio, A., 2014. The seabed geomorphology and geological
structure of the Firth of Lorn, western Scotland, UK, as revealed by
multibeam echo-sounder survey. Earth Environ. Sci. Trans. R. Soc.
Edinburgh 105, 273–284. <https://doi.org/10.1017/S1755691015000146>
