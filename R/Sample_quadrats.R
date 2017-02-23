#' Plot-based samples from a spatially-explicit census
#'
#' This function allows to sample quadratic subplots from a spatially-explicit
#' community. The output format is a sites x species abundance table which is
#' the classical data format used in community ecology. The table generated
#' can be for instance be further analysed with the package \code{\link[vegan]}.
#'
#' @param comm Community object from which the samples are generated
#' @param n_quadrats Number of sampling quadrats
#' @param quadrat_area Area of the sampling quadrats
#' @param plot Logical: should the sampling design be plotted?
#' @param method Available methods are \code{"random", "transect", "grid"}
#' @param avoid_overlap Logical: for the random sampling try to generate a design
#' without overlap of quadrats
#' @param x0, y0 Lower left corner of the first quadrat in transect and grid sampling
#' @param delta_x Distance between consecutive quadrats in transect and grid sampling
#' in x-direction (the distance between the left sides is measured)
#' @param delta_y Distance between consecutive quadrats in transect and grid sampling
#' in y-direction (the distance between the lower sides is measured)
#'
#' @return  A matrix with species abundances where the rows represents sampling
#' quadrats and the columns represent species
#'
#' @examples
#' library(vegan)
#' sim_com1 <- Sim.Poisson.Community(100, 10000)
#' comm_mat1 <- sample_quadrats(sim_com1, n_quadrats = 100,
#' quadrat_area = 0.002, method = "grid")
#' renyi(comm_mat1, hill = T)
#'
#' @export
#'
sample_quadrats <- function(comm, n_quadrats = 10, quadrat_area = 0.01,
                            plot = T, method = "random",
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

      if (method == "random"){

         xpos <- runif(n_quadrats, min = comm$x_min_max[1], max = comm$x_min_max[2] - quadrat_size)
         ypos <- runif(n_quadrats, min = comm$y_min_max[1], max = comm$y_min_max[2] - quadrat_size)

         coords <- cbind(xpos,ypos)
         min_dist <- sqrt(2*quadrat_size^2)

         if (avoid_overlap == T){

            count <- 0

            while(min(dist(coords)) < min_dist && count <= 9999){
               xpos <- runif(n_quadrats, min = comm$x_min_max[1], max = comm$x_min_max[2] - quadrat_size)
               ypos <- runif(n_quadrats, min = comm$y_min_max[1], max = comm$y_min_max[2] - quadrat_size)

               coords <- cbind(xpos,ypos)
               count <- count + 1
            }

            if (count > 9999) stop("Cannot find a sampling layout with no overlap")

         } else {
            if (min(dist(coords)) < min_dist)
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
         if (min(dist(coords)) < 0.9999*quadrat_size)
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

         if (min(dist(coords)) < 0.9999*quadrat_size)
            warning("There are overlapping sampling squares in the design")

      } # end grid

      comm_tab <- mapply(abund_rect, xpos, ypos,
                        MoreArgs=list(xsize = quadrat_size, ysize = quadrat_size,
                                      comm = comm))

   } else {
      comm_tab <- table(comm$census$species)

      xpos <- comm$x_min_max[1]
      ypos <- comm$y_min_max[1]
   }

   comm_tab <- t(comm_tab)
   rownames(comm_tab) <- paste("site", 1:nrow(comm_tab), sep = "")

   # plot sampling design
   if (plot == TRUE){
      plot(comm)
      rect(xpos, ypos, xpos + quadrat_size, ypos + quadrat_size, lwd = 2,
           col = adjustcolor("white", alpha.f = 0.6))
   }

   return(comm_tab)

}


