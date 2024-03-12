## Resubmission

This is a resubmission. In this version we have:

* Fixed the problems that caused mobsim to be archived on CRAN: essentially by
removing the dependency on spatstat.core and using spatstat.random and 
spatstat.geom instead.
* Fixed inst/CITATION and DESCRIPTION
* Added features to major functions
* Added Depends R >= 4.0.0 because mobsim depends on sads which itself depends on VGAM which depends on R 4.0.0.

## Test environments

* local - Darwin, R 4.3.1
* win-builder (release and devel)
* r hub - Ubuntu Linux 20.04.1 LTS, R-release, GCC
* r hub - Fedora Linux, R-devel, clang, gfortran

## R CMD check results

0 errors | 0 warnings | 0 notes

We ran revdepcheck::revdep_check() locally and no packages depend on mobsim
any more.

Best wishes,

Felix may and Alban Sagouis
