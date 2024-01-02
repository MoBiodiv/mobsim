#' Adds a random movement to mother point coordinates
#'
#' @param mpcoords (mother_map) List of class \code{mother_map}.
#' @param sd (numeric) Standard deviation for the \code{\link[stats]{rnorm}}.
#' @inheritParams create_motherpoints
#'
#' @returns A named list of class \code{mother_map} with three elements:
#' \enumerate{
#'   \item \code{n_mother_points} a vector of integers giving the number of mother points per species
#'   \item \code{xmother} List of length equal to the number of species. Each list element
#' is a vector of x coordinates for every mother points. If one element is NA, the
#' the corresponding species is not clustered.
#'   \item \code{ymother} List of length equal to the number of species. Each list element
#' is a vector of y coordinates for every mother points. If one element is NA, the
#' the corresponding species is not clustered.
#' }
#'
#' @details The function takes a \code{mother_map} object and adds a random
#' jitter value to each coordinate. The jitter follows a normal probability
#' distribution of \code{mean} = 0 and \code{sd} specified by the used with
#' argument \code{sd}.
#'
#' @author Alban Sagouis
#'
#'
#' @examples
#' # Default behaviour
#' mpcoords <- create_motherpoints(s_local = 10L)
#' mpcoordsJ <- jitter_motherpoints(mpcoords)
#' par(mfrow = c(1L, 2L))
#' plot(mpcoords)
#' plot(mpcoordsJ)
#'
#' # With more variability for the movements of motherpoints
#' mpcoordsJ <- jitter_motherpoints(mpcoords, sd = 2)
#'
#' # Fixing the seed
#' mpcoordsJ <- jitter_motherpoints(mpcoords, seed = 42L)
#'
#' # Integration with mobsim::sim_thomas_coords()
#' abund_vec <- 11L:100L # a community with 90 species
#' s_local <- length(abund_vec)
#' mpcoords <- create_motherpoints(s_local, seed = 42L)
#' simdat <- mobsim::sim_thomas_coords(abund_vec,
#'  xmother = mpcoords$xmother,	# list of vectors
#'  ymother = mpcoords$ymother	# list of vectors
#' )
#' ## computing biodiversity metrics
#' div1 <- mobsim::div_rand_rect(prop_area = 0.01, comm = simdat, n_rect = 10L, exclude_zeros = TRUE)
#'
#' ## moving motherpoints slightly
#' mpcoordsJ <- jitter_motherpoints(mpcoords)
#' simdat2 <- mobsim::sim_thomas_coords(abund_vec,
#'  xmother = mpcoords$xmother,	# list of vectors
#'  ymother = mpcoords$ymother	# list of vectors
#' )
#' ## computing biodiversity metrics
#' div2 <- mobsim::div_rand_rect(prop_area = 0.01, comm = simdat2,
#'                               n_rect = 10L, exclude_zeros = TRUE)
#'
#' @export

jitter_motherpoints <- function(mpcoords, sd = 0.01, seed = NULL) {
   if (!is.null(seed)) set.seed(seed)

   mpcoords$xmother <- lapply(mpcoords$xmother, function(x) x + stats::rnorm(
      n = length(x),
      mean = 0,
      sd = sd))

   mpcoords$ymother <- lapply(mpcoords$ymother, function(y) y + stats::rnorm(
      n = length(y),
      mean = 0,
      sd = sd))

   return(mpcoords)
}
