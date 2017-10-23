# Test environments
* local ubuntu 16.04 LTS, R 3.4.2
* win-builder (devel and release)

# R CMD check results

test


## ubuntu 16.04

0 errors | 0 warnings | 0 notes

R CMD check succeeded

## win-builder (devel)

Status: 1 WARNING, 1 NOTE

* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Felix May <felix.may@posteo.de>'

* checking top-level files ... WARNING
Conversion of 'README.md' failed:
pandoc.exe: Could not fetch README-unnamed-chunk-4-1.png
README-unnamed-chunk-4-1.png: openBinaryFile: does not exist (No such file or directory)

I use an R markdown file (README.Rmd) to create README.md
In the ode examples in README.Rmd figures are created. This seems to be a problem
for the windows build. If this is a serious problem I can easily remove the plotting
functions from the examples. However, it would be nice to have a graphical example
in the README file.

## win-builder (release)

Status: 1 WARNING, 1 NOTE

* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Felix May <felix.may@posteo.de>'

* checking top-level files ... WARNING
Conversion of 'README.md' failed:
pandoc.exe: Could not fetch README-unnamed-chunk-4-1.png
README-unnamed-chunk-4-1.png: openBinaryFile: does not exist (No such file or directory)

# Downstream dependencies

I ran devtools::revdep_check() locally

0 packages with problems

No ERRORs or WARNINGs found :)
