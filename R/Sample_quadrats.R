#' Plot-based samples from a spatially-explicit census
#'
#' This function allows to sample quadratic subplots from a spatially-explicit
#' community. The output format are a sites x species abundance table and a
#' sites x xy-coordinates table. The sites x species abundance is
#' a classical data format used in community ecology. The table generated
#' can be for instance be further analysed with the package \code{\link{vegan}}.
#'
#' @param comm Community object from which the samples are generated
#' @param n_quadrats (integer) Number of sampling quadrats
#' @param quadrat_area (numeric) Area of the sampling quadrats
#' @param plot (logical) Should the sampling design be plotted? default to TRUE.
#' @param method (character) Available methods are \code{"random", "transect", "grid"}
#' @param avoid_overlap (logical) For the random sampling try to generate a design
#' without overlap of quadrats . Default is TRUE.
#' @param x0,y0 (numeric value) Lower left corner of the first quadrat in transect and grid sampling
#' @param delta_x (numeric value) Distance between consecutive quadrats in transect and grid sampling
#' in x-direction (the distance between the left sides is measured)
#' @param delta_y (numeric value) Distance between consecutive quadrats in transect and grid sampling
#' in y-direction (the distance between the lower sides is measured)
#' @param seed (integer) Any integer passed to \code{set.seed} for reproducibility.
#'
#' @return  A list with two items, \code{spec_dat} and \code{xy_dat}.
#' \code{spec_dat} is a data.frame with sampling quadrats in rows and species abundances
#' in columns, and \code{xy_dat} is a data.frame with sampling quadrats in rows
#' and the xy-coordinates of the quadrats (lower left corner) in columns.
#'
#' @examples
#' library(vegan)
#' sim_com1 <- sim_poisson_community(100, 10000)
#' comm_mat1 <- sample_quadrats(sim_com1, n_quadrats = 100,
#' quadrat_area = 0.002, method = "grid")
#' specnumber(comm_mat1$spec_dat)
#' diversity(comm_mat1$spec_dat, index = "shannon")
#'
#'
#' @export
#'
sample_quadrats <- function(comm, n_quadrats = 20, quadrat_area = 0.01,
                            plot = TRUE, method = "random",
                            avoid_overlap = TRUE,
                            x0 = 0, y0 = 0, delta_x = 0.1, delta_y = 0.1,
                            seed = NULL)
{
   if (class(comm) != "community")
      stop("comm has to be a community object")
   if (round(n_quadrats, 0) != n_quadrats) stop("n_quadrats has to be an integer")

   method <- match.arg(method, c("random", "transect", "grid"))

   quadrat_size <- sqrt(quadrat_area)

   community_size <- (comm$x_min_max[2] - comm$x_min_max[1]) *
      (comm$y_min_max[2] - comm$y_min_max[1])

   if (quadrat_size < community_size) { # This works as long as we use square quadrats

      if (n_quadrats == 1L) {
         xy_dat <- sampling_one_quadrat(
            xmin = comm$x_min_max[1], xmax = comm$x_min_max[2] - quadrat_size,
            ymin = comm$y_min_max[1], ymax = comm$y_min_max[2] - quadrat_size,
            seed = seed
         )
      } else {

         min_dist <- sqrt(2*quadrat_area)

         if (method == "random") {

            if (isTRUE(avoid_overlap)) {

               if (requireNamespace("spatstat.random", quietly = TRUE)) {
                  xy_dat <- sampling_random_spatstat(
                     n_quadrats = n_quadrats, min_dist = min_dist,
                     xmin = comm$x_min_max[1], xmax = comm$x_min_max[2] - quadrat_size,
                     ymin = comm$y_min_max[1], ymax = comm$y_min_max[2] - quadrat_size,
                     seed = seed
                  )

               } else {
                  xy_dat <- sampling_random_bruteforce(
                     n_quadrats = n_quadrats, min_dist = min_dist,
                     xmin = comm$x_min_max[1], xmax = comm$x_min_max[2] - quadrat_size,
                     ymin = comm$y_min_max[1], ymax = comm$y_min_max[2] - quadrat_size,
                     seed = seed
                  )

               } # end of no spatstat

            } else { # if isFALSE(avoid_overlap)
               xy_dat <- sampling_random_overlap(
                  n_quadrats = n_quadrats, min_dist = min_dist,
                  xmin = comm$x_min_max[1], xmax = comm$x_min_max[2] - quadrat_size,
                  ymin = comm$y_min_max[1], ymax = comm$y_min_max[2] - quadrat_size,
                  seed = seed
               )
            }

         } # end method == random

         if (method == "transect") {
            xy_dat <- sampling_transects(
               n_quadrats = n_quadrats,
               xmin = comm$x_min_max[1], xmax = comm$x_min_max[2],
               ymin = comm$y_min_max[1], ymax = comm$y_min_max[2],
               x0 = x0, y0 = y0, delta_x = delta_x, delta_y = delta_y,
               quadrat_size = quadrat_size
            )
         } # end transect

         if (method == "grid") {
            xy_dat <- sampling_grids(
               n_quadrats = n_quadrats,
               xmin = comm$x_min_max[1], xmax = comm$x_min_max[2],
               ymin = comm$y_min_max[1], ymax = comm$y_min_max[2],
               x0 = x0, y0 = y0, delta_x = delta_x, delta_y = delta_y,
               quadrat_size = quadrat_size
            )
         } # end grid
      }

      comm_tab <- mapply(abund_rect, xy_dat$x, xy_dat$y,
                         MoreArgs = list(xsize = quadrat_size, ysize = quadrat_size,
                                         comm = comm))
      spec_dat <- as.data.frame(t(comm_tab))
      rownames(spec_dat) <- paste("site", 1L:nrow(spec_dat), sep = "")

      rownames(xy_dat) <- rownames(spec_dat)

   } else {
      comm_tab <- table(comm$census$species)

      xy_dat <- data.frame(
         x = comm$x_min_max[1],
         y = comm$y_min_max[1]
      )

      rownames(xy_dat) <- "site1"
      # spec_dat should be a data.frame, 1 row, n_species columns, rowname would be "site1"
      spec_dat <- as.data.frame(matrix(comm_tab, nrow = 1L, byrow = TRUE, dimnames = list("site1", names(comm_tab))))

   }

   # plot sampling design
   if (plot == TRUE) {
      plot(comm)
      graphics::rect(xy_dat$x, xy_dat$y, xy_dat$x + quadrat_size, xy_dat$y + quadrat_size, lwd = 2,
                     col = grDevices::adjustcolor("white", alpha.f = 0.6))
   }

   return(list(spec_dat = spec_dat, xy_dat = xy_dat))

}




#' Creates coordinates (lower left corner of a quadrat) randomly distributed but
#' without overlapping each other
#'
#' Efficient algorithm from package \code{spatstat.random} is used.
#' Produces similar results as \code{\link{sampling_random_bruteforce}}.
#' @inheritParams  sample_quadrats
#' @param n_quadrats Number of sampling quadrats
#' @param min_dist (numeric) minimal distance between two points to avoid overlap.
#' Equal to the length of a quadrat diagonal
#' @param xmin (numeric) minimum possible value on the x axis a quadrat can cover.
#' @param xmax (numeric) maximum possible value on the x axis a quadrat can cover.
#' @param ymin (numeric) minimum possible value on the y axis a quadrat can cover.
#' @param ymax (numeric) maximum possible value on the y axis a quadrat can cover.
#'
#' @return  a data.frame with 2 columns x and y giving the coordinates of the
#' lower left corner of the square quadrats.
#' @export
#'

sampling_random_spatstat <- function(n_quadrats, min_dist, xmin, xmax, ymin, ymax, seed = NULL) {
   if (!is.null(seed)) set.seed(seed)
   hc_mod <- list(cif = "hardcore", par = list(beta = n_quadrats, hc = min_dist),
                  w = spatstat.geom::owin(c(xmin, xmax), c(ymin, ymax)))
   hc_points <- spatstat.random::rmh(model = hc_mod, start = list(n.start = n_quadrats),
                                     control = list(p = 1, nrep = 1e6),
                                     saveinfo = F, verbose = F)
   xpos <- hc_points$x
   ypos <- hc_points$y

   coords <- data.frame(x = xpos, y = ypos)
   if (min(stats::dist(coords)) < min_dist)
      warning("There are overlapping sampling squares in the design. Use less quadrats or smaller quadrat area.")
   return(coords)
}


#' Creates coordinates (lower left corner of a quadrat) randomly distributed but without overlapping each
#' other
#'
#' This function works without having the \code{spatstat.random} package install.
#'
#' @inheritParams sampling_random_spatstat
#'
#' @return  a data.frame with 2 columns x and y giving  the coordinates of the
#' lower left corner of the square quadrats.
#' @export
#'
sampling_random_bruteforce <- function(n_quadrats, min_dist, xmin, xmax, ymin, ymax, seed = NULL) {
   if (!is.null(seed)) set.seed(seed)
   count <- 0L
   maxCount <- 10000L

   xpos <- stats::runif(n_quadrats, min = xmin, max = xmax)
   ypos <- stats::runif(n_quadrats, min = ymin, max = ymax)
   coords <- cbind(xpos, ypos)

   while (min(stats::dist(coords)) < min_dist && count <= maxCount) {
      xpos <- stats::runif(n_quadrats, min = xmin, max = xmax)
      ypos <- stats::runif(n_quadrats, min = ymin, max = ymax)

      coords <- cbind(xpos,ypos)
      count <- count + 1L
   }

   if (count > maxCount) warning("Cannot find a sampling layout with no overlap.
Install the package spatstat for an improved method for non-overlapping squares,
Use less quadrats or smaller quadrat area, or set avoid_overlap to FALSE.")
   colnames(coords) <- c("x", "y")

   return(as.data.frame(coords))
}


#' Creates coordinates (lower left corner of a quadrat) randomly distributed
#' that may overlap each other
#'
#' @inheritParams  sampling_random_spatstat
#'
#' @return  a data.frame with 2 columns x and y giving  the coordinates of the
#' lower left corner of the square quadrats.
#' @export
#'
sampling_random_overlap <- function(n_quadrats, min_dist, xmin, xmax, ymin, ymax, seed = NULL) {
   if (!is.null(seed)) set.seed(seed)
   xpos <- stats::runif(n_quadrats, min = xmin, max = xmax)
   ypos <- stats::runif(n_quadrats, min = ymin, max = ymax)

   coords <- data.frame(x = xpos, y = ypos)
   if (min(stats::dist(coords)) < min_dist)
      warning("There are overlapping sampling squares in the design")

   return(coords)
}

#' Creates square quadrats aligned along a transect
#'
#' @inheritParams  sample_quadrats
#' @inheritParams  sampling_random_spatstat
#' @param quadrat_size (numeric) width of the quadrats.
#'
#' @return  a data.frame with 2 columns x and y giving  the coordinates of the
#' lower left corner of the square quadrats.
#' @export
#'
sampling_transects <- function(n_quadrats, xmin, xmax, ymin, ymax, x0, y0, delta_x, delta_y, quadrat_size) {
   t_xmin <- x0
   t_ymin <- y0

   t_xmax <- x0 + (n_quadrats - 1) * delta_x + quadrat_size
   t_ymax <- y0 + (n_quadrats - 1) * delta_y + quadrat_size

   if (t_xmin < xmin || t_xmax > xmax)
      stop("x-extent of sampling design is larger than landscape")

   if (t_ymin < ymin || t_ymax > ymax)
      stop("y-extent of sampling design is larger than landscape")

   xpos <- seq(from = x0, by = delta_x, len = n_quadrats)
   ypos <- seq(from = y0, by = delta_y, len = n_quadrats)

   coords <- data.frame(x = xpos, y = ypos)
   if (min(stats::dist(coords)) < 0.9999*quadrat_size)
      warning("There are overlapping sampling squares in the design")
   return(coords)
}


#' Creates square quadrats aligned on a regular grid
#'
#' @inheritParams  sample_quadrats
#' @inheritParams  sampling_random_spatstat
#' @param quadrat_size (numeric) width of the quadrats.
#'
#' @return  a data.frame with 2 columns x and y giving  the coordinates of the
#' lower left corner of the square quadrats.
#' @export
#'
sampling_grids <- function(n_quadrats, xmin, xmax, ymin, ymax, x0, y0, delta_x, delta_y, quadrat_size) {
   grid_dim <- sqrt(ceiling(sqrt(n_quadrats))^2)

   x1 <- seq(from = x0, by = delta_x, len = grid_dim)
   y1 <- seq(from = y0, by = delta_y, len = grid_dim)

   if (min(x1) < xmin || max(x1) > xmax)
      stop("x-extent of sampling design is larger than landscape")

   if (min(y1) < ymin || max(y1) > ymax)
      stop("y-extent of sampling design is larger than landscape")

   coords <- expand.grid(x = x1, y = y1)[1L:n_quadrats, ]

   if (min(stats::dist(coords)) < 0.9999*quadrat_size)
      warning("There are overlapping sampling squares in the design")

   return(coords)
}

#' Creates one square quadrat randomly located in the landscape
#'
#' @inheritParams  sampling_random_spatstat
#'
#' @return  a data.frame with 2 columns x and y giving  the coordinates of the
#' lower left corner of the square quadrat.
#' @export
#'
sampling_one_quadrat <- function(xmin, xmax, ymin, ymax, seed = NULL) {
   if (!is.null(seed)) set.seed(seed)
   data.frame(
      x = stats::runif(1L, min = xmin, max = xmax),
      y = stats::runif(1L, min = ymin, max = ymax)
   )
}





