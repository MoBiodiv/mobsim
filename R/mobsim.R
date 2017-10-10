#' mobsim: A package for spatial analysis of scale-dependent biodiversity changes.
#'
#' The package includes functions to simulate species distributions in space
#' as well as for the analysis of spatially-explicit data, where each individual
#' is described by its xy-coordinates and a species identity label.
#'
#' @section Functions to simulate species abundances and distributions:
#'
#' \code{\link{sim_sad}}
#'
#' \code{\link{sim_poisson_coords}}
#'
#' \code{\link{sim_thomas_coords}}
#'
#' \code{\link{sim_poisson_community}}
#'
#' \code{\link{sim_thomas_community}}
#'
#'
#' @section Functions to analyse species abundances and distributions:
#'
#' \code{\link{rare_curve}}
#'
#' \code{\link{spec_sample_curve}}
#'
#' \code{\link{divar}}
#'
#' \code{\link{dist_decay}}
#'
#' \code{\link{sample_quadrats}}
#'
#' @author Felix May
#'
#' @docType package
#'
#' @name mobsim
NULL

#' @useDynLib mobsim, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL
