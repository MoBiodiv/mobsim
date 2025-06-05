#' Adds the same fixed movement to all individuals of a species on x
#'
#' @param drift (numeric) A single value added to the x or y coordinates.
#' @inheritParams jitter_species
#' @rdname drift_species
#'
#' @returns A named list of class \code{community} as described in
#' \code{\link[mobsim]{community}}
#'
#' @details The function takes a \code{community} object and adds a fixed
#' drift value to either x or y coordinate per species. In a 1 by 1 landscape, a
#'  drift of 0.5 means moving over half the map.
#' @seealso drift_motherpoints
#' @author Alban Sagouis
#'
#'
#' @examples
#' # Default behaviour
#' simdat <- mobsim::sim_thomas_community(s_pool = 5L, n_sim = 100L,
#'                                        mother_points = 1L)
#' simdatD <- simdat |>
#'              drift_x_species(0.1) |>
#'              drift_y_species(-0.2)
#' par(mfrow = c(1L, 2L))
#' plot(simdat)
#' plot(simdatD)
#'
#' @export
#'

drift_x_species <- function(comm, drift = NULL) {
  comm$census$x <- comm$census$x + drift

  return(comm)
}

#' Adds the same fixed movement to all individuals of a species on y
#' @inherit drift_x_species
#' @rdname drift_species
#' @export
#'

drift_y_species <- function(comm, drift = NULL) {
  comm$census$y <- comm$census$y + drift

  return(comm)
}
