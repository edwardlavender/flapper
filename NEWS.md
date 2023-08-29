# `flapper` (development version)

* **Datasets**
    * Rebuild (`RasterLayer`) datasets to fix `no slot of name "srs" for this object of class "RasterLayer"` warning for `dat_gebco` (and downstream errors) & error in CRS slot in algorithm outputs
    * Replace non-ASCII characters in the `comment` attributes of Coordinate Reference System (CRS) objects to fix the
    warnings for these during R CMD check (`Warning: found non-ASCII strings PROJCRS...`)
    * Set `flapper_run_slow` and `flapper_run_parallel` to `FALSE` on remote

* **Function fixes**
    * Fix bug in `.acs_pl()` in internal check of receiver overlaps
    * Fix re-sampling issue in `pf_simplify()`
    * Fix and refine calculations in `get_detection_days()`
    * Fixes to functions after dependencies update, including: 
        - `get_detection_*()` functions (e.g., issue with data slot in SpatialPointsDataFrame);
        - `lcp_over_surface()`
        - `process_false_detections_sf()` examples 

* **Function improvements**
    * Update `.acs_pl()` and `.acs()` internals to check for centroid overlap (to prevent unclear error messages)
    * Update `get_mvt_mobility_from_archival()` with `step_check` argument to check for regular time series
    * Update `dist_btw_receivers()` to handle lon/lat and planar coordinates and return dataframes or matrices
    * Update `check_*()` functions with internal improvements
    * Update `update_extent()` to handle Spatial* objects
    * Update selected `sim_*()` functions (namely `sim_array()` and `sim_path_sa()`) to support RasterLayer inputs for the coastline (resolving issues associated with the discrepancy between simulations and gridded algorithm implementations)
        

* **Function tests**
    * Set up testing with `testthat`;
    
* **Documentation**
    * Update `{glatos}` source (thanks @chrisholbrook)
    * Add `flapper_algorithms_faqs` vignette
    * Fix typos
    * Add vignette & startup message

# `flapper` 1.0.0

* The first release of `flapper`

# `flapper` 0.0.0.9000

* The developmental version of `flapper`
