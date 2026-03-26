#' Adds a fixed drift movement to mother points x coordinates
#'
#' @param drift (numeric) Amount of spatial drift
#' @inheritParams jitter_motherpoints
#' @rdname drift_motherpoints
#'
#' @returns A named list of class \code{mother_map} with three elements:
#' \enumerate{
#'   \item \code{n_mother_points} a vector of integers giving the number of
#'    mother points per species
#'   \item \code{xmother} List of length equal to the number of species. Each
#'    list element is a vector of x coordinates for every mother points. If one
#'     element is NA, the the corresponding species is not clustered.
#'   \item \code{ymother} List of length equal to the number of species. Each
#'    list element is a vector of y coordinates for every mother points. If one
#'    element is NA, the the corresponding species is not clustered.
#' }
#'
#' @details The function takes a \code{mother_map} object and adds a fixed
#' drift value to either the x or y coordinates. In a 1 by 1 landscape, a drift
#' of 0.5 means moving over half the map.
#'
#' @author Alban Sagouis
#'
#'
#' @examples
#' # Default behaviour
#' mpcoords <- create_motherpoints(s_local = 10L)
#' mpcoordsD <- drift_x_motherpoints(mpcoords, 0.1)
#' par(mfrow = c(1L, 2L))
#' plot(mpcoords)
#' plot(mpcoordsD)
#'
#' # With more distance for the movements of motherpoints
#' mpcoordsD <- mpcoords |>
#'                drift_x_motherpoints(drift = 0.5) |>
#'                drift_y_motherpoints(drift = 0.5)
#'
#' # Integration with mobsim::sim_thomas_coords()
#' abund_vec <- 11L:100L # a community with 90 species
#' s_local <- length(abund_vec)
#' mpcoords <- create_motherpoints(s_local, seed = 42L)
#' simdat <- mobsim::sim_thomas_coords(abund_vec,
#'   xmother = mpcoords$xmother,	# list of vectors
#'   ymother = mpcoords$ymother	# list of vectors
#' )
#' ## computing biodiversity metrics
#' div1 <- mobsim::div_rand_rect(prop_area = 0.01, comm = simdat,
#'                             n_rect = 10L, exclude_zeros = TRUE)
#'
#' ## moving motherpoints slightly
#' mpcoordsD <- drift_x_motherpoints(mpcoords, 0.01) |>
#'                drift_y_motherpoints(drift = 0.01)
#' simdat2 <- mobsim::sim_thomas_coords(abund_vec,
#'  xmother = mpcoords$xmother,	# list of vectors
#'  ymother = mpcoords$ymother	# list of vectors
#' )
#' ## computing biodiversity metrics
#' div2 <- mobsim::div_rand_rect(prop_area = 0.01, comm = simdat2,
#'                               n_rect = 10L, exclude_zeros = TRUE)
#'
#' @export

drift_x_motherpoints <- function(mpcoords, drift) {
  mpcoords$xmother <- lapply(mpcoords$xmother, function(x) x + drift)

  return(mpcoords)
}

#' Adds a fixed drift movement to mother points y coordinates
#' @inherit drift_x_motherpoints
#' @rdname drift_motherpoints
#' @export

drift_y_motherpoints <- function(mpcoords, drift) {
  mpcoords$ymother <- lapply(mpcoords$ymother, function(x) x + drift)

  return(mpcoords)
}
