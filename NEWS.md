# flapper (development version)

* **Function fixes**
    * Fix and refine calculations in `get_detection_days()`
    * Fixes to functions after dependencies update, including: 
        - `get_detection_clumps()` 
        - `lcp_over_surface()`

* **Function improvements**
    * Update `get_mvt_mobility_from_archival()` with `step_check` argument to check for regular time series
    * Update `dist_btw_receivers()` to handle lon/lat and planar coordinates and return dataframes or matrices
    * Update `check_*()` functions with internal improvements

* **Function tests**
    * Set up testing with `testthat`;
    
* **Documentation**
    * Add `flapper_algorithms_faqs` vignette

# flapper 1.0.0

* The first release of `flapper`

# flapper 0.0.9000

* The developmental version of `flapper`
