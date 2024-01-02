#' Creates a random number of mother points for a given number of species
#'
#' @param s_local (integer) Species number.
#' @param n_motherpoints_range (integer vector) Vector of integer values for the number of motherpoints a species can get.
#' @param xrange (numeric vector) Extent of the community in x-direction.
#' @param yrange (numeric vector) Extent of the community in y-direction.
#' @param seed (integer) Any integer passed to \code{\link[base]{set.seed}} for reproducibility.
#'
#' @returns A named list of class \code{mother_map} with three elements:
#' \enumerate{
#'   \item   \code{n_mother_points} a vector of integers giving the number of mother points per species
#'   \item   \code{xmother} List of length equal to the number of species. Each list element
#' is a vector of x coordinates for every mother points. If one element is NA, the
#' the corresponding species is not clustered.
#'   \item   \code{ymother} List of length equal to the number of species. Each list element
#' is a vector of y coordinates for every mother points. If one element is NA, the
#' the corresponding species is not clustered.
#' }
#'
#' @details For each species, the function first draws a random number of
#' mother points in the given \code{n_motherpoints_range}. Then the function
#' places as many mother points as needed in the window given by \code{xrange}
#' and \code{yrange}.
#' @author Alban Sagouis
#'
#'
#' @examples
#' # Default behaviour
#' mpcoords <- create_motherpoints(s_local = 10L)
#'
#' # With more variability for the number of motherpoints
#' mpcoords <- create_motherpoints(s_local = 10L, n_motherpoints_range = 1L:100L)
#'
#' # Allowing no clustering at all
#' mpcoords <- create_motherpoints(s_local = 10L, n_motherpoints_range = 0L:5L)
#'
#' # Fixing the seed
#' mpcoords <- create_motherpoints(10L, seed = 42L)
#'
#' # Integration with mobsim::sim_thomas_coords()
#' abund_vec <- 11L:100L # a community with 90 species
#' s_local <- length(abund_vec)
#' mpcoords <- create_motherpoints(s_local, seed = 42L)
#' simdat <- mobsim::sim_thomas_coords(abund_vec,
#'  xmother = mpcoords$xmother,	# list of vectors
#'  ymother = mpcoords$ymother	# list of vectors
#' )
#' # Integration with mobsim::sim_thomas_coords() without clustering
#' abund_vec <- 20L:25L
#' s_local <- length(abund_vec)
#' mpcoords <- create_motherpoints(s_local, n_motherpoints_range = 0L:1L, seed = 42L)
#' simdat <- mobsim::sim_thomas_coords(abund_vec,
#'  xmother = mpcoords$xmother,	# list of vectors
#'  ymother = mpcoords$ymother	# list of vectors
#' )
#' plot(simdat)
#'
#' @export


create_motherpoints <- function(s_local, n_motherpoints_range = 1:4,
                                xrange = c(0, 1), yrange = c(0, 1),
                                seed = NULL) {
   if (!is.null(seed)) set.seed(seed)

   n_mother_points <- sample(x = n_motherpoints_range,
                             size = s_local,
                             replace = TRUE)
   xmother <- sapply(n_mother_points, stats::runif,
                     min = xrange[1L], max = xrange[2L])
   ymother <- sapply(n_mother_points, stats::runif,
                     min = yrange[1L], max = yrange[2L])

   mpcoords <- list(
      n_mother_points = n_mother_points,
      xmother = xmother, ymother = ymother,
      x_min_max = xrange, y_min_max = yrange
   )
   class(mpcoords) <- c("mother_map", "list")

   return(mpcoords)
}


#' Plot mother_map object object
#'
#' Plot positions and species identities of all mother points created by
#' \code{\link{create_motherpoints}} or \code{\link{jitter_motherpoints}}.
#'
#' @param x mother_map object
#' @param col Colour vector to mark species identities
#' @param pch Plotting character to mark species identities.
#'     pch 16 is advised for large data sets.
#' @param ... Other parameters to \code{\link[graphics:plot.default]{graphics::plot}}
#'
#' @rdname create_motherpoints
#'
#' @examples
#' mpcoords <- create_motherpoints(10L)
#' plot(mpcoords)
#'
#' @export
#'
plot.mother_map <- function(x, col = NULL, pch = NULL, ...) {
   nspec <- length(x$n_mother_points)
   if (is.null(col))  col <- grDevices::rainbow(nspec, alpha = 1)
   if (is.null(pch))  pch <- 16

   graphics::plot(x = unlist(x$xmother), y = unlist(x$ymother),
                  xlim = x$x_min_max, ylim = x$y_min_max,
                  col = rep(x = seq_along(x$n_mother_points),
                            times = x$n_mother_points),
                  xlab = "X", ylab = "Y",
                  pch = pch, las = 1, asp = 1, ...)
}
