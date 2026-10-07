#' Adds the same random movement to all individuals of a species
#'
#' @param comm (community) List of class \code{community}.
#' @inheritParams jitter_motherpoints
#'
#' @returns A named list of class \code{community} as described in
#' \code{\link[mobsim]{community}}
#'
#' @details The function takes a \code{community} object and adds a random
#' jitter value to each coordinate per species. The jitter follows a normal
#' probability distribution of \code{mean} = 0 and \code{sd} specified by the
#' used with argument \code{sd}.
#' @seealso jitter_motherpoints
#' @author Alban Sagouis
#'
#'
#' @examples
#' # Default behaviour
#' simdat <- mobsim::sim_thomas_community(s_pool = 5L, n_sim = 100L,
#'                                        mother_points = 1L)
#' simdatJ <- jitter_species(simdat)
#' par(mfrow = c(1L, 2L))
#' plot(simdat)
#' plot(simdatJ)
#'
#' # With more variability for the movement of individuals
#' simdatJ <- jitter_species(simdat, sd = 2)
#'
#' # Fixing the seed
#' simdatJ <- jitter_species(simdat, seed = 42L)
#'
#' @export

jitter_species <- function(comm, sd = 0.01, seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }

  sp <- unique(comm$census$species)
  x_shift <- stats::setNames(
    stats::rnorm(n = length(sp), mean = 0L, sd = sd),
    nm = sp
  )
  y_shift <- stats::setNames(
    stats::rnorm(n = length(sp), mean = 0L, sd = sd),
    nm = sp
  )
  comm$census$x <- comm$census$x + x_shift[comm$census$species]
  comm$census$y <- comm$census$y + y_shift[comm$census$species]

  return(comm)
}
