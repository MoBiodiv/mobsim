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
   if (!is.null(seed)) set.seed(seed)

   comm$census$x <- unlist(tapply(X = comm$census$x,
                                  INDEX = comm$census$species,
                                  FUN = function(x) x + stats::rnorm(n = 1L,
                                                                     mean = 0,
                                                                     sd = sd)
   )
   )
   comm$census$y <- unlist(tapply(X = comm$census$y,
                                  INDEX = comm$census$species,
                                  FUN = function(y) y + stats::rnorm(n = 1L,
                                                                     mean = 0,
                                                                     sd = sd)))

   return(comm)
}

