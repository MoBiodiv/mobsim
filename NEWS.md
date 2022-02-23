mobsim 0.2.1
================================================================================

### POTENTIAL CODE BREAK
* In sample_quadrats() default value for `avoid_overlap` argument was changed from FALSE to TRUE.
* sim_sad() now throws an error if s_pool or n_sim is not a whole number.

### NEW FEATURES
In `sim_thomas_community()` and `sim_thomas_coords()` users can now specify:
* a range (rectangular) for each species
* different numbers of mother points (0 (no clustering), 1 or more for each species
* different standard deviations for each species
* coordinates for the mother points

### MINOR IMPROVEMENTS
* spatstat.core updated to spatstat.random. Thanks to @rubak for pointing it out.
* increased speed for sim_thomas_coords() when there is no clustering.
* significant increase of test code coverage.

### DEPENDENCE
testthat version 3.0.0 or above is needed to test and build from source.

### ACKNOWLEDGEMENTS CHANGES
* new contributor to the package: @AlbanSagouis
* added a CITATION file for easier citation of the 2018 article in MEE by May et al.
* added a CONTRIBUTING file to encourage and describe how to warn about bugs and
propose changes.
* added a codemeta.json file to improve mobsim referencing.

mobsim 0.2.0
================================================================================

* Update of `sim_thomas_process`.
* Addition of `dist_decay_quadrats`.
