
# flapper

***From passive acoustic telemetry to space use: an `R` package for
integrating passive acoustic telemetry and archival tag data to estimate
benthic movement pathways at high resolution***

[![Project Status: WIP – Initial development is in progress, but there
has not yet been a stable, usable release suitable for the
public.](https://www.repostatus.org/badges/latest/wip.svg)](https://www.repostatus.org/#wip)

`flapper` is an `R` package which provides tools for passive acoustic
telemetry (PAT) data. The package has been particularly motivated by the
collection of new acoustic and archival data from a Critically
Endangered elasmobranch, the flapper skate (*Dipturus intermedius*), off
the west coast of Scotland where a static PAT array has been established
to examine the movements of individuals within a Marine Protected Area.
`flapper` has been designed to complement existing packages for the
analyses of these data
(e.g. [VTrack](https://github.com/RossDwyer/VTrack),
[glatos](https://gitlab.oceantrack.org/GreatLakes/glatos) and
[fishtrack3d](https://github.com/aspillaga/fishtrack3d) and
[actel](https://github.com/hugomflavio/actel)), with a particular focus
on the provision of tools that integrate PAT and archival data for
improved inference of patterns of space use, especially for
benthic/demersal species. To this end, `flapper` contains functions in
the following themes:

  - **Data processing tools**, including data assembly (e.g., of
    range-testing datasets), checking false detections and implementing
    quality checks;
  - **Spatial tools**, including common spatial operations for
    manipulation of spatial data, such as polygon inversion;
  - **Distance calculations**, including the calculation of distances
    between receivers, along 3-dimensional movement paths, and of the
    shortest paths over a surface;
  - **Detection statistics**, including metrics of sampling effort, such
    as detection area; and individual detection metrics, such as
    detection days and co-occurrence;
  - **Space use algorithms**, including a straightforward implementation
    of the mean-position algorithm for the estimation of centres of
    activity and new algorithms designed for benthic/demersal species
    which integrate PAT and archival data to improve estimates of space
    use;
  - **Simulation tools**, including tools for the simulation of passive
    acoustic telemetry arrays, movement paths, detections and the
    comparison of simulated and inferred patterns of space use under
    different conditions;

<img src="vignettes/readme_context.png"/> *flapper: An `R` package
designed to facilitate the integration of acoustic and archival datasets
to improve estimates of space use for benthic/demersal species. Inserted
sample depth and acoustic time series were collected as part of the
Movement Ecology of Flapper Skate project by Marine Scotland Science and
NatureScot. The insert of the flapper skate is also courtesy of this
project. The bathymetry data are sourced from the Ireland, Northern
Island and Scotland Hydrographic survey (Howe et al., 2015). Plots were
produced using the
[prettyGraphics](https://github.com/edwardlavender/prettyGraphics)
package.*

## Highlights

The main highlights of the package are the provision of routines for the
rapid calculation of biologically meaningful distances in area with
complex barriers to movement (e.g., coastline) alongside algorithms
(most of which are exclusive to `flapper`) for inferring space use from
discrete detections at receivers, especially:

  - **`lcp_over_surface()`.** This function calculates the shortest
    pathways and/or their distances over a surface using efficient `C++`
    algorithms from the
    [cppRouting](https://github.com/vlarmet/cppRouting) package. This
    makes it easy to use biologically meaningful distances (that account
    for the bathymetric depth over which a benthic animal must move, if
    applicable, and barriers to movement) in movement models.
  - **`coa()`.** This function implements the arithmetic version of the
    mean-position algorithm to estimate centres of activity (COAs) from
    discrete detections at receivers, given only a detection matrix and
    the locations of receivers.
  - **`dc()`.** This function implements the ‘depth-contour’ (DC)
    algorithm to examine patterns of space use of benthic/demersal
    animals. This relates one-dimensional depth time series to a
    two-dimensional bathymetry surface to determine the extent to which
    different parts of an area might have (or have not) been used, or
    effectively represent occupied depths, over time.
  - **`acdc()`.** This function implements the ‘acoustic-centroids
    depth-contour’ (ACDC) algorithm to examine patterns of space use.
    This is a new approach that integrates the locational information
    provided by acoustic detections and concurrent depth observations to
    infer where benthic/demersal animals could have spent more or less
    time over the period of observations.
  - **`acdcpf()`** and **`acdcmp()`** are extensions of the ACDC
    algorithm that will be added in the package in the near future.
  - **`sim_*()`.** These functions provide flexible routines for the
    simulation of receiver arrays, movement paths and detections to
    evaluate alternative algorithms for inferences about patterns of
    space use under different conditions.

## Installation

This package requires `R` version ≥ 4.0. You can check your current
version with `R.version.string`. Subsequent installation steps (may)
require the `devtools` and `pkgbuild` packages, which can be installed
with `install.packages(c("devtools", "pkgbuild"))`. On Windows, package
building requires `Rtools`. You can check whether `Rtools` is installed
with `pkgbuild::has_rtools()`. If `Rtools` is not installed, it is
necessary to download and install the appropriate version of `Rtools`
before proceeding by following the instructions
[here](https://cran.r-project.org/bin/windows/Rtools/).

Three packages
([prettyGraphics](https://github.com/edwardlavender/prettyGraphics),
[Tools4ETS](https://github.com/edwardlavender/Tools4ETS) and
[glatos](https://gitlab.oceantrack.org/GreatLakes/glatos)) required or
suggested from [GitHub](https://github.com) repositories (since they are
not currently available from [CRAN](https://cran.r-project.org)). These
can be installed during the installation process (see below), but it is
safer to install them sequentially as follows:

    devtools::install_github("edwardlavender/prettyGraphics")
    devtools::install_github("edwardlavender/Tools4ETS")
    devtools::install_url("https://gitlab.oceantrack.org/GreatLakes/glatos/repository/master/archive.zip",
                          build_opts = c("--no-resave-data", "--no-manual"))

To install these packages with their vignettes, add `dependencies =
TRUE` and `build_vignettes = TRUE` as arguments to the code above (see
`?devtools::install_github` or `?devtools::install_url` for further
information). Then, you can install the development version of `flapper`
from [GitHub](https://github.com/edwardlavender/flapper) as shown below:

``` r
devtools::install_github("edwardlavender/flapper", dependencies = TRUE, build_vignettes = TRUE)
```

The `dependencies = TRUE` argument will also install any suggested
packages, which may be required by some functions and to build vignettes
(which will be added to the package in due course). To access the
vignettes, use `vignette("flapper_intro", package = "flapper")` for a
general introduction to the package. *Note that vignettes have not yet
been added to the package.*

## Example datasets

A key feature of the `flapper` package is that almost all functions are
designed to be implemented using standard object types (e.g., dataframes
and matrices) rather than package-specific object classes. For
simplicity, `flapper` makes some assumptions about variable names that
follow a consistent and logical structure (e.g., individual IDs are
given as `individual_id` and receiver IDs are given as `receiver_id`)
but, apart from these assumptions, this structure means the functions in
the package are accessible and straightforward to use.

`flapper` functions are illustrated using simulated data and the
following sample data collected from flapper skate off the west coast of
Scotland:

  - `dat_ids` is a dataset containing the characteristics of a sample of
    tagged flapper skate;
  - `dat_moorings` is a dataset containing some sample passive acoustic
    telemetry receiver locations and associated information;
  - `dat_acoustics` is a dataset containing some sample detection time
    series;
  - `dat_archival` is a dataset containing some sample depth time
    series;
  - `dat_sentinel` is a dataset containing some sample transmission –
    detection time series assembled from sentinel tags;

These example datasets were collected by Marine Scotland Science and
NatureScot as part of the Movement Ecology of Flapper Skate project and
belong to these organisations. If you wish to use these data, please
contact Marine Scotland Science and NatureScot for further information.

## Data processing tools

A number of functions facilitate the assembly, processing and checking
of passive acoustic telemetry time series:

  - **Data assembly.**
      - `assemble_sentinel_counts()` assembles counts of
        transmissions/detections from sentinel tags for modelling
        purposes (i.e., to model detection probability);
      - `make_matrix_*()` functions create matrices of individual and
        receiver deployment time series and detection time series:
          - `make_matrix_ids()` matricises individual deployment time
            series;
          - `make_matrix_receivers()` matricises receiver deployment
            time series;
          - `make_matrix_detections()` matricises detection time series;
      - `make_df_*()` functions (i.e., `make_df_detections()`) reverse
        this process;
  - **Data processing.**
      - `process_receiver_id()` adds unique receiver IDs to a dataframe
        (which is useful if the same receiver has been deployed more
        than once);
      - `process_false_detections_sf()` passes putative false detections
        through a spatial filter which incorporates ancillary
        information on receiver locations and animal swimming speeds to
        interrogate their plausibility;
      - `process_quality_check()` passes acoustic data through some
        basic quality checks prior to analysis;

## Spatial tools

A number of functions facilitate spatial operations, which support
ecological investigations and space use algorithms:

  - `buffer_and_crop()` buffers a spatial object (e.g., receiver
    locations) and uses this buffered object to crop another (e.g., the
    local bathymetry);
  - `cells_from_val()` returns the cells (or a raster of the cells) of a
    raster that are equal to a specified value or lie within a specified
    range of values;
  - `invert_poly()` inverts a polygon (e.g, to define the ‘sea’ from a
    polygon of the ‘land’);
  - `mask_io()` masks values in a raster that lie inside or outside of a
    spatial mask (e.g., to mask the ‘land’ from the ‘sea’);
  - `sim_surface()` populates a raster with simulated values;
  - `split_raster_equally()` splits a raster into equal pieces (using
    code from the [greenbrown](http://greenbrown.r-forge.r-project.org)
    package);
  - `update_extent()`shrinks or inflates an extent object;

## Distance calculations

Some functions facilitate standard distance calculations using Euclidean
distances:

  - `dist_btw_receivers()` calculates the Euclidean distances between
    all combinations of receivers;
  - `dist_btw_points_3d()` calculates the Euclidean distances between
    points in three-dimensional space;
  - `dist_over_surface()` calculates the Euclidean distance along a path
    over a three-dimensional surface;

Often, Euclidean distances may not be a suitable representation of
distance. This is especially the case for coastal benthic/demersal
species in bathymetrically complex environments, for which navigation
between locations may require movement over hilly terrain and around
coastline. Thus, `lcp_over_surface()` calculates shortest pathways
and/or the distances of the shortest pathways over a surface between
origin and destination coordinates.

## Detection statistics

A number of functions facilitate the calculation of detection
statistics, including those related to sampling effort and to detections
of individuals:

  - `get_detection_pr()` calculates detection probability given a model
    for detection probability with distance;
  - `get_detection_centroids()` defines detection centroids (areas
    within the maximum detection range) around receivers;
  - `get_detection_centroids_envir()` extracts environmental conditions
    from within receiver detection ranges, accounting for detection
    probability;
  - `get_detection_area_sum()` calculates the total area surveyed by
    receivers;
  - `get_detection_area_ts()` defines a time series of the area surveyed
    by receivers;
  - `get_n_operational_ts()`defines a time series of the number of
    operational units (e.g., individuals at liberty or active receivers)
  - `get_id_rec_overlap()` calculates the overlap between the deployment
    periods of tagged individuals and receivers;
  - `make_matrix_cooccurence()` computes a detection history similarity
    matrix across individuals;

## Space use algorithms

The main thrust of `flapper` is the implementation of existing and new
algorithms designed to infer space use from PAT data and their
evaluation under different circumstances (e.g., array designs, movement
and detection models).

### The centres of activity (COA) algorithm

Centres of activity (COA) are one of the most widely used metrics for
the reconstruction of patterns of space use from PAT data. Several
methods have been developed to calculate COAs, but the mean-position
algorithm is the commonest. To generate estimates of space use, COAs are
usually taken as point estimates from which utilisation distributions
(typically kernel utilisation distributions, KUDs) are estimated.
`flapper` facilitates the implementation of this approach with the
following functions:

  - `coa_setup_delta_t()` informs decisions as to an appropriate time
    interval over which to calculate COAs;
  - `make_matrix_detections()` summarises matrices over time intervals
    (see above);
  - `coa()` implements the arithmetic version of the mean-position
    algorithm to calculate COAs;
  - `kud_around_coastline()` facilitates the estimation of home ranges
    (e.g., from estimated COAs) in areas of complex coastline;

### The depth contour (DC) algorithm

Alongside the COA algorithm, `flapper` introduces a number of new
algorithms for the inferring patterns of space use. The depth-contour
(DC) algorithm is the simplest. Whereas the COA approach only makes use
of detections, the DC approach only uses depth observations.
Specifically, under the assumption that individuals are
benthic/demersal, this algorithm uses observed depths (± some error) to
define the subset of possible locations of each individual within a
defined area. This is implemented via `dc()`.

### The acoustic-centroid depth-contour (ACDC) algorithm

The acoustic-centroid depth-contour (ACDC) algorithm extends the DC
algorithm by using PAT data to inform the area within which depth
contours are most likely to be found. This algorithm is supported by a
number of functions:

  - `acdc_setup_n_centroids()` suggests the number of acoustic centroids
    for the algorithm;
  - `acdc_setup_centroids()` defines the acoustic centroids for the
    algorithm;
  - `acdc()` implements the algorithm, via the back-end function
    `.acdc()`;
  - `acdc_simplify()` simplifies the results of the algorithm;
  - `acdc_plot()` plots the results of the algorithm;
  - `acdc_animate()` creates html animations of the algorithm;

### The ACDC particle filtering and movement pathway (ACDCPF and ACDCMP) algorithms

The ACDCPF and ACDCMP algorithms incorporate movement pathways into the
ACDC process to restrict further the inferred distribution of locations
within which the individual must have been located at each time point.
*These algorithms are not currently available in the public version of
this package.*

## Simulation tools

Simulations are a valuable tool in ecology which can elucidate the
relative performance of alternative methods for ecological inferences
(e.g., the `COA` approach versus the `DC` approach for inferring
patterns of space use) and the extent to which new data sources
influence ecological inferences under different circumstances (e.g. the
extent to which sparse or regular PAT detections improve estimates of
space use). To this end, `flapper` provided joined-up routines for the
evaluation of approaches for the estimation of patterns of space use
under different conditions; namely:

  - `sim_array()` simulates alternative array designs;  
  - `sim_path_*()` functions simulate discrete-time movement paths,
    including:
      - `sim_path_sa()`, supported by `sim_steps()` and `sim_angles()`,
        simulates movement paths (possibly in restricted areas) from
        step lengths and turning angles;
      - `sim_path_ou_1()` simulates movement paths under
        Ornstein-Uhlenbeck processes;
  - `sim_detections()` simulates detections at receivers arising from
    movement paths under a diversity of detection probability models;

To evaluate the performance of alternative algorithms for inferring
patterns of space use under different array designs, movement and
detections models, `eval_by_kud()` compares patterns of space use
inferred from simulated and estimated movement paths using kernel
utilisation distributions.

## Associated packages

  - **[prettyGraphics](https://github.com/edwardlavender/prettyGraphics)**
    facilitates the production of pretty, publication-quality and
    interactive visualisations, with a particular focus on time series.
    This makes it easy to create abacus plots, visualise time series
    (across factor levels, at different temporal scales and in relation
    to covariates), bathymetric landscapes and movement pathways in
    three-dimensions, and detection similarity matrices.
  - **[Tools4ETS](https://github.com/edwardlavender/Tools4ETS)**
    provides a set of general tools for ecological time series,
    including the definition of time categories, matching time series
    (e.g., detection observations with environmental covariates),
    flagging independent time series and simulating time series.
  - **[fvcom.tbx](https://github.com/edwardlavender/fvcom.tbx)**
    provides tools for the integration of hydrodynamic model predictions
    (from the Finite Coastal Ocean Volume Model) with ecological
    datasets (e.g., detection time series). This facilitates the
    inclusion of hydrodynamic model predictions as covariates in
    movement models and the validation of hydrodynamic model predictions
    with movement datasets or data collected from static acoustic
    receivers. This package was particularly motivated by the West
    Scotland Coastal Ocean Modelling System (WeStCOMS).

## References

Howe, J.A., Anderton, R., Arosio, R., Dove, D., Bradwell, T., Crump, P.,
Cooper, R., Cocuccio, A., 2014. The seabed geomorphology and geological
structure of the Firth of Lorn, western Scotland, UK, as revealed by
multibeam echo-sounder survey. Earth Environ. Sci. Trans. R. Soc.
Edinburgh 105, 273–284. <https://doi.org/10.1017/S1755691015000146>
