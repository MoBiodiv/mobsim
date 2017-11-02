## Resubmission

This is a resubmission. In this version I have:

* Formatted the reference according to the suggestions of Uwe Ligges and changed from 
a reference using the bioRXiv labe to the DOI

* Included a reference to the paper describing the package in the DESCRIPTION file, see
[](https://www.biorxiv.org/content/early/2017/10/26/209502), as suggested by the CRAN team
member Swetlana Herbrandt

* Fixed small typos in the vignettes and documentation

## Test environments
* local ubuntu 16.04 LTS, R 3.4.2
* win-builder (devel and release)

## R CMD check results

There were no ERRORs or WARNINGs. 

There was 1 NOTE:

* checking CRAN incoming feasibility ... NOTE
Maintainer: ‘Felix May <felix.may@posteo.de>’

It seems this is normal for the first submission of a package

## Downstream dependencies

I ran devtools::revdep_check() locally

0 packages with problems

No ERRORs or WARNINGs found :)
