# Utilities

## torusify() ----
#' Rounds coordinates around like on a torus
#'
#' @param mpcoords Can be either a data.frame with columns \code{xmother} and
#' \code{ymother} or an object of class \code{mother_map} or a \code{numeric}
#' vector.
#' @param ... Kept for future uses.
#'
#' @details \code{torusify()} is a generic function that calls functions adapted
#' to the class of argument \code{mapcoords}. Classes \code{data.frame},
#'  \code{numeric} and \code{mother_map} are accepted at the moment.
#'
#' \code{torusify()} takes coordinates and makes sure that all
#' coordinates actually inside the window specified
#' in \code{x_min_max} and \code{y_min_max}.
#'
#' @returns Object returned depends on the object given. For an object of class
#' \code{mother_map} as given by \code{\link{create_motherpoints}}, see
#' \code{\link{torusify.mother_map}}. For a \code{data.frame}, see
#' \code{\link{torusify.data.frame}} and for a \code{numeric}  vector, see
#' \code{\link{torusify.numeric}}.
#'
#' @examples
#' # Default behaviour
#' mpcoords <- create_motherpoints(s_local = 10L)
#' mpcoordsJ <- jitter_motherpoints(mpcoords)
#' mcoordsT <- torusify(mpcoordsJ)
#' @export


torusify <- function(mpcoords, ...) {
   UseMethod("torusify", mpcoords)
}


#' Rounds coordinates around like on a torus
#'
#' @param mpcoords A named list of class \code{mother_map} with three elements:
#' \enumerate{
#'   \item \code{n_mother_points} a vector of integers giving the number of
#'    mother points per species \item \code{xmother} List of length equal to the
#'    number of species. Each list element is a vector of x coordinates for
#'    every mother points. If one element is NA, the the corresponding species
#'    is not clustered.
#'   \item \code{ymother} List of length equal to the number of species. Each
#'   list element is a vector of y coordinates for every mother points. If one
#'   element is NA, the the corresponding species is not clustered.
#' }
#' @param ... Kept for future uses.
#'
#' @returns An object of class \code{mother_map} of the same dimensions as the
#' one provided to the function.
#' @details The function takes a \code{mother_map} object and makes sure that
#' all coordinates that were just jittered are actually inside the window
#' specified in the \code{x_min_max} and \code{y_min_max} elements of the
#' \code{mother_map} object.
#'
#' @author Alban Sagouis
#'
#' @examples
#' # Default behaviour
#' mpcoords <- create_motherpoints(s_local = 10L)
#' mpcoordsJ <- jitter_motherpoints(mpcoords)
#' mcoordsT <- torusify(mpcoordsJ)
#'
#' @rdname torusify
#' @export



torusify.mother_map <- function(mpcoords, ...) {
   mpcoords$xmother <- lapply(
      X = mpcoords$xmother,
      FUN = function(MPs)
         (MPs %% diff(mpcoords$x_min_max)) + mpcoords$x_min_max[1L])
   mpcoords$ymother <- lapply(
      X = mpcoords$ymother,
      FUN = function(MPs)
         (MPs %% diff(mpcoords$y_min_max)) + mpcoords$y_min_max[1L])

   return(mpcoords)
}


#' Rounds coordinates around like on a torus
#'
#' @param mpcoords A list of class \code{\link[mobsim]{community}}
#' @param ... Kept for future uses.
#'
#' @returns A list of class \code{community} of same dimension as
#' \code{mpcoords}.
#'
#' @details The function takes a \code{community} and makes sure that all
#' individual coordinates are actually inside the window specified by the
#' \code{x_min_max} and \code{y_min_max} elements of the \code{community}
#' object.
#'
#' @author Alban Sagouis
#'
#' @examples
#' simdat <- mobsim::sim_thomas_community(10L, 100L, mother_points = 1)
#' simdatJ <- jitter_species(simdat, sd = 1)
#' simdatJ <- torusify(simdatJ)
#'
#' @rdname torusify
#' @export

torusify.community <- function(mpcoords, ...) {
   mpcoords$census$x <- mpcoords$census$x %% diff(mpcoords$x_min_max) + mpcoords$x_min_max[1L]
   mpcoords$census$y <- mpcoords$census$y %% diff(mpcoords$y_min_max) + mpcoords$y_min_max[1L]
   return(mpcoords)
}


#' Rounds coordinates around like on a torus
#'
#' @param mpcoords A data.frame with two columns \code{xmother} and \code{ymother}
#' @param x_min_max vector of 2 values describing the window on the x axis
#' @param y_min_max vector of 2 values describing the window on the y axis
#' @param ... Kept for future uses.
#'
#' @returns A \code{data.frame} of same dimension as \code{mpcoords}.
#'
#' @details The function takes a \code{data.frame} and makes sure that all
#' coordinates are actually inside the window specified
#' by the \code{x_min_max} and \code{y_min_max} arguments.
#'
#' @author Alban Sagouis
#'
#' @examples
#' x_min_max <- c(0, 1)
#' y_min_max <- c(0, 1)
#' mpcoords <- data.frame(xmother = stats::runif(100L),
#'                        ymother = stats::runif(100L))
#' mcoordsT <- torusify(mpcoords, x_min_max, y_min_max)
#'
#' @rdname torusify
#' @export

torusify.data.frame <- function(mpcoords, x_min_max, y_min_max, ...) {
   if (!all(colnames(mpcoords) %in% c("xmother","ymother")))
      stop("column names should be xmother and ymother")
   mpcoords$xmother <- mpcoords$xmother %% diff(x_min_max) + x_min_max[1L]
   mpcoords$ymother <- mpcoords$ymother %% diff(y_min_max) + y_min_max[1L]
   return(mpcoords)
}



#' Rounds coordinates around like on a torus
#'
#' @param mpcoords A vector of class \code{numeric}.
#' @param range vector of 2 values describing the window.
#' @param ... Kept for future uses.
#'
#' @returns A \code{numeric} vector of same length as \code{mpcoords}.
#'
#' @details The function takes a \code{numeric} vector and makes sure that all
#' values are actually inside the window specified
#' by the \code{range} argument.
#'
#' @author Alban Sagouis
#'
#' @examples
#' mpcoords <- stats::runif(10L, min = -1, max = 2)
#' torusify(mpcoords, range = c(0, 1))
#'
#' @rdname torusify
#' @export

torusify.numeric <- function(mpcoords, range, ...) {
   return(mpcoords %% diff(range) + range[1L])
}



## community2ppp() ----
#' Converts an object of class \code{\link[mobsim]{community}} into an object of
#' class \code{\link[spatstat.geom]{ppp}}.
#'
#' @param comm Community object
#'
#' @details For each species, the function first draws a random number of
#' mother points in the given \code{n_motherpoints_range}. Then the function
#' places as many mother points as needed in the window given by \code{xrange}
#' and \code{yrange}.
#' @author Alban Sagouis
#'
#' @examples
#' # Integrated between mobsim::sim_thomas_coords() and spatstat.explore::marktable()
#' simdat <- mobsim::sim_thomas_coords(abund_vec = 11L:100L)
#' simdat <- community2ppp(comm = simdat)
#' spatstat.explore::marktable(X = simdat, N = 10L)
#'
#' @export

community2ppp <- function(comm) {
   spatstat.geom::ppp(
      x = comm$census$x,
      y = comm$census$y,
      marks = comm$census$species,
      window = spatstat.geom::owin(
         xrange = comm$x_min_max,
         yrange = comm$y_min_max
      )
   )
}

## create_random_ID() ----
#' Creates a random string.
#'
#' @param n (integer) number of random strings to create. Default is 1.
#' @inheritParams jitter_motherpoints
#'
#' @export
#'

create_random_ID <- function(n = 1L, seed = NULL) {
   oldseed <- .Random.seed
   on.exit({.Random.seed <<- oldseed})
   if (!is.null(seed)) {
      set.seed(seed)
   } else {
      set.seed(Sys.time())
   }

   a <- do.call(paste0,
                replicate(5L,
                          sample(LETTERS, n, replace = TRUE),
                          simplify = FALSE))
   return(
      paste0(a,
             sprintf("%04d", sample(9999L, n, replace = TRUE)),
             sample(LETTERS, n, replace = TRUE))
   )
}
