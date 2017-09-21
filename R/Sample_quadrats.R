#' Plot-based samples from a spatially-explicit census
#'
#' This function allows to sample quadratic subplots from a spatially-explicit
#' community. The output format are a sites x species abundance table and a
#' sites x xy-coordinates table. The sites x species abundance is
#' a classical data format used in community ecology. The table generated
#' can be for instance be further analysed with the package \code{\link{vegan}}.
#'
#' @param comm Community object from which the samples are generated
#' @param n_quadrats Number of sampling quadrats
#' @param quadrat_area Area of the sampling quadrats
#' @param plot Should the sampling design be plotted? (logical)
#' @param method Available methods are \code{"random", "transect", "grid"}
#' @param avoid_overlap For the random sampling try to generate a design
#' without overlap of quadrats (logical)
#' @param x0,y0 Lower left corner of the first quadrat in transect and grid sampling
#' @param delta_x Distance between consecutive quadrats in transect and grid sampling
#' in x-direction (the distance between the left sides is measured)
#' @param delta_y Distance between consecutive quadrats in transect and grid sampling
#' in y-direction (the distance between the lower sides is measured)
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
                            avoid_overlap = F,
                            x0 = 0, y0 = 0, delta_x = 0.1, delta_y = 0.1)
{
   if (class(comm) != "community")
      stop ("comm has to be a community object")

   method <- match.arg(method, c("random", "transect","grid"))

   quadrat_size <- sqrt(quadrat_area)

   community_size <- (comm$x_min_max[2] - comm$x_min_max[1]) *
      (comm$y_min_max[2] - comm$y_min_max[1])

   if (quadrat_size < community_size){

      if (n_quadrats > 1){

         min_dist <- sqrt(2*quadrat_area)

         if (method == "random"){

            if (avoid_overlap == TRUE){

               if (requireNamespace("spatstat", quietly = TRUE)){

                  # Hard core process:
                  hc_mod <- list(cif = "hardcore", par = list(beta = n_quadrats, hc = min_dist),
                                 w = spatstat::owin(c(0, 1 - quadrat_size), c(0, 1 - quadrat_size)))
                  hc_points <- spatstat::rmh(model = hc_mod, start=list(n.start = n_quadrats),
                                             control=list(p = 1, nrep = 1e6),
                                             saveinfo = F, verbose = F)
                  xpos <- hc_points$x
                  ypos <- hc_points$y

                  coords <- cbind(xpos,ypos)
                  if (min(stats::dist(coords)) < min_dist)
                     warning("There are overlapping sampling squares in the design. Use less quadrats or smaller quadrat area.")
               } else {

                  count <- 0

                  xpos <- stats::runif(n_quadrats, min = comm$x_min_max[1],
                                       max = comm$x_min_max[2] - quadrat_size)
                  ypos <- stats::runif(n_quadrats, min = comm$y_min_max[1],
                                       max = comm$y_min_max[2] - quadrat_size)
                  coords <- cbind(xpos,ypos)

                  while(min(stats::dist(coords)) < min_dist && count <= 999){
                     xpos <- stats::runif(n_quadrats, min = comm$x_min_max[1],
                                                      max = comm$x_min_max[2] - quadrat_size)
                     ypos <- stats::runif(n_quadrats, min = comm$y_min_max[1],
                                                      max = comm$y_min_max[2] - quadrat_size)

                     coords <- cbind(xpos,ypos)
                     count <- count + 1
                  }

                  if (count > 999) warning("Cannot find a sampling layout with no overlap.
                                            Install the package spatstat for an improved method for non-overlapping squares,
                                            Use less quadrats or smaller quadrat area, or set avoid_overlap to FALSE.")

               } # of no spatstat

            } else {
               xpos <- stats::runif(n_quadrats, min = comm$x_min_max[1],
                                                max = comm$x_min_max[2] - quadrat_size)
               ypos <- stats::runif(n_quadrats, min = comm$y_min_max[1],
                                                max = comm$y_min_max[2] - quadrat_size)

               coords <- cbind(xpos,ypos)
               if (min(stats::dist(coords)) < 0.9999*quadrat_size)
                  warning("There are overlapping sampling squares in the design")
            }

         } # end method == random

         if (method == "transect"){

            xmin <- x0
            ymin <- y0

            xmax <- x0 + (n_quadrats - 1) * delta_x + quadrat_size
            ymax <- y0 + (n_quadrats - 1) * delta_y + quadrat_size

            if (xmin < comm$x_min_max[1] || xmax > comm$x_min_max[2])
               stop ("x-extent of sampling desing is larger than landscape")

            if (ymin < comm$y_min_max[1] || ymax > comm$y_min_max[2])
               stop ("y-extent of sampling desing is larger than landscape")

            xpos <- seq(from = x0, by = delta_x, len = n_quadrats)
            ypos <- seq(from = y0, by = delta_y, len = n_quadrats)

            coords <- cbind(xpos,ypos)
            if (min(stats::dist(coords)) < 0.9999*quadrat_size)
               warning("There are overlapping sampling squares in the design")
         } # end transect

         if (method == "grid"){

            grid_dim <- sqrt(ceiling(sqrt(n_quadrats))^2)

            x1 <- seq(from = x0, by = delta_x, len = grid_dim)
            y1 <- seq(from = y0, by = delta_y, len = grid_dim)

            if (min(x1) < comm$x_min_max[1] || max(x1) > comm$x_min_max[2])
               stop ("x-extent of sampling desing is larger than landscape")

            if (min(y1) < comm$y_min_max[1] || max(y1) > comm$y_min_max[2])
               stop ("y-extent of sampling desing is larger than landscape")

            coords <- expand.grid(xpos = x1, ypos = y1)

            xpos <- coords$xpos[1:n_quadrats]
            ypos <- coords$ypos[1:n_quadrats]

            if (min(stats::dist(coords)) < 0.9999*quadrat_size)
               warning("There are overlapping sampling squares in the design")

         } # end grid

      } else { # if n_quadrats == 1

         xpos <- stats::runif(1, min = comm$x_min_max[1], max = comm$x_min_max[2] - quadrat_size)
         ypos <- stats::runif(1, min = comm$y_min_max[1], max = comm$y_min_max[2] - quadrat_size)

      }

      comm_tab <- mapply(abund_rect, xpos, ypos,
                        MoreArgs=list(xsize = quadrat_size, ysize = quadrat_size,
                                      comm = comm))

   } else {
      comm_tab <- table(comm$census$species)

      xpos <- comm$x_min_max[1]
      ypos <- comm$y_min_max[1]
   }

   spec_dat <- as.data.frame(t(comm_tab))
   rownames(spec_dat) <- paste("site", 1:nrow(spec_dat), sep = "")

   xy_dat <- data.frame(x = xpos, y = ypos)
   rownames(xy_dat) <- rownames(spec_dat)

   # plot sampling design
   if (plot == TRUE){
      plot(comm)
      graphics::rect(xpos, ypos, xpos + quadrat_size, ypos + quadrat_size, lwd = 2,
                     col = grDevices::adjustcolor("white", alpha.f = 0.6))
   }

   return(list(spec_dat = spec_dat, xy_dat = xy_dat))

}


