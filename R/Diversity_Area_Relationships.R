
#' Local diversity indices
#'
#' Get diversity indices including species richness, no. of endemics,
#' Shannon and Simpson diversity for one rectangle subplot in the community
#'
#' @param x0 x-coordinate of lower left corner
#' @param y0 x-coordinate of lower left corner
#' @param xsize size of the subplot in x-direction
#' @param ysize size of the subplot in y-direction
#' @param comm \code{\link{community}} object
#'
#' @return Named vector with four values
#' \enumerate{
#'    \item n_species - the number of species
#'    \item n_endemics - the number of endemics
#'    \item shannon - Shannon diversity index defined as \eqn{H = - \sum p_{i} * log(p_[i])},
#'    where \eqn{p_i} is the relative abundance of species i:
#'    \item simpson - Simpson diversity index (= probabiliy of interspecific encounter PIE)
#'    defined as \eqn{D =  1 - \sum p_i^2}
#' }
#'
#' @export
#'
div_rect <- function(x0, y0, xsize, ysize, comm)
{
   x <- comm$census$x
   y <- comm$census$y

   # logical vector which trees are in the sampling rectangle
   in_rect <- (x >= x0 & x < (x0 + xsize) & y >= y0 & y < (y0 + ysize))

   spec_in <- unique(comm$census$species[in_rect])
   spec_out <- unique(comm$census$species[!in_rect])

   n_species <- length(spec_in)
   n_endemics <- length(spec_in[!spec_in %in% spec_out])

   abund <- table(comm$census$species[in_rect])
   abund <- abund[abund > 0]
   relabund <- abund/sum(abund)

   shannon <- - sum(relabund * log(relabund))

   n <- sum(abund)
   if (n > 1)
      simpson <- (n/(n-1)) * (1 - sum(relabund^2))
   else
      simpson <- NA

   return(c(n_species = n_species,
            n_endemics = n_endemics,
            shannon = shannon,
            ens_shannon = exp(shannon),
            simpson = simpson,
            ens_simpson = 1/(1 - simpson)
            ))
}

#' -----------------------------------------------------------------------------
#' Distribution of local diversity indices
#'
#' Get mean and sd of diversity indices in several equally sized subplots
#' of a community
#'
#' @param prop_area Size of subplots as proportion of the total area
#' @param comm \code{\link{community}} object
#' @param n_rect Number of randomly located subplots
#' @param exclude_zeros logical - should subplots without individuals be excluded?
#'
#' @return Vector with mean and standard deviation of the following diversity
#' indices: (1) Number of species (2) Number of endemics (3) Shannon-diversity
#' (4) Simpson diversity
#'
#' @seealso \code{\link{div_rect}}
#'
#' @export
#'
div_rand_rect <- function(prop_area = 0.25, comm, n_rect = 100,
                          exclude_zeros = F)
{
   dx_plot <- comm$x_min_max[2] - comm$x_min_max[1]
   dy_plot <- comm$y_min_max[2] - comm$y_min_max[1]

   area <- dx_plot * dy_plot * prop_area
   square_size <- sqrt(area)

   if (square_size <= min(c(dx_plot, dy_plot))){
      dx_rect <- square_size
      dy_rect <- square_size
   } else
   {
      if (dx_plot >= dy_plot){
         dx_rect <- dx_plot*prop_area
         dy_rect <- dy_plot
      } else {
         dx_rect <- dx_plot
         dy_rect <- dy_plot*prop_area
      }
   }

   xpos <- runif(n_rect, min = comm$x_min_max[1], max = comm$x_min_max[2] - dx_rect)
   ypos <- runif(n_rect, min = comm$y_min_max[1], max = comm$y_min_max[2] - dy_rect)

   div_plots <- mapply(div_rect, xpos, ypos,
                       MoreArgs=list(xsize = dx_rect, ysize = dy_rect,
                                     comm = comm))
   if (exclude_zeros == T)
      div_plots <- div_plots[, div_plots["n_species",] > 0]

   return(c(species    = mean(div_plots["n_species",]),
            # se_species      = sd(div_plots["n_species",]),
            endemics     = mean(div_plots["n_endemics",]),
            # se_end       = sd(div_plots["n_endemics",]),
            shannon = mean(div_plots["shannon",]),
            # sd_shannon   = sd(div_plots["shannon",]),
            ens_shannon = mean(div_plots["ens_shannon",]),
            # sd_ens_shannon   = sd(div_plots["ens_shannon",]),
            simpson = mean(div_plots["simpson",], na.rm = T),
            # sd_simpson   = sd(div_plots["simpson",], na.rm = T),
            ens_simpson = mean(div_plots["ens_simpson",], na.rm = T)
            # sd_ens_simpson   = sd(div_plots["ens_simpson",], na.rm = T))
          )
   )
}

#' -----------------------------------------------------------------------------
#' Diversity-area relationships
#'
#' Estimate diversity indices in subplots of different sizes. This includes the
#' well-known species-area and endemics-area relationships.
#'
#' @param prop_area numeric vector with subplot sizes as proportion of the total area
#' @param comm \code{\link{community}} object
#' @param n_samples Number of randomly located subplots per subplot size
#' @param exclude_zeros logical - should subplots without individuals be excluded?
#'
#' @return Dataframe with the proportional area of the subplots and mean and
#' standard deviation of the following diversity indices: (1) Number of species
#' (2) Number of endemics (3) Shannon-diversity (4) Simpson diversity
#'
#' @seealso \code{\link{div_rand_rect}}, \code{\link{div_rect}}
#'
#' @examples
#' sim_com1 <- Sim.Thomas.Community(100, 1000)
#' divar1 <- DivAR(sim_com1, prop_area = seq(0.01, 1.0, length = 20))
#' plot(meanSpec ~ propArea, data = divar1, xlab = "Proportion of area",
#'      ylab = "No. of species", type = "b", ylim = c(0,100))
#' points(meanEnd ~ propArea, data = divar1, type = "b", col = 2)
#'
#' @export
#'
divar <- function(comm, prop_area = seq(0.1, 1, by = 0.1), n_samples = 100,
                  exclude_zeros = T)
{
   if (any(prop_area > 1))
      warning("Subplot areas larger than the community size are ignored!")
   prop_area <- prop_area[prop_area <= 1]

   if (class(comm) != "community")
      stop("DiVAR requires a community object as input. See ?community.")

   n_scales <- length(prop_area)
   dx_plot <- comm$x_min_max[2] - comm$x_min_max[1]
   dy_plot <- comm$y_min_max[2] - comm$y_min_max[1]

   div_area <- sapply(prop_area,
                      div_rand_rect,
                      comm = comm,
                      n_rect = n_samples,
                      exclude_zeros = exclude_zeros)

   div_dat <- as.data.frame(t(div_area))
   div_dat <- cbind(prop_area = prop_area, div_dat)

   class(div_dat) <- c("divar", "data.frame")

   return(div_dat)
}

#' -----------------------------------------------------------------------------
#' Plot diversity-area relationships
#'
#' @param divar
#'
#' @export
plot.divar <- function(divar)
{
   #http://stackoverflow.com/questions/8929663/r-legend-placement-in-a-plot

   #Plot an empty graph and legend to get the size of the legend
   max_spec <- max(divar$species)
   plot(species ~ prop_area, data = divar, ylim = c(0, max_spec),
        type = "n")

   my_legend_size <-legend("topright",
                           legend = c("Species","Endemics", "ENS Shannon ","ENS Simpson"),
                           lwd = 2, plot = FALSE)
   dev.off()

   #custom ylim. Add the height of legend to upper bound of the range
   my_range <- c(0, max_spec)
   my_range[2] <- 1.1*(my_range[2] + my_legend_size$rect$h)

   #draw the plot with custom ylim
   plot(species ~ prop_area, data = divar, ylim = my_range,
        type = "n", las = 1,
        xlab = "Proportion of area sampled", ylab = "No. of species",
        main = "Diversity-area relationships")
   lines(species ~ prop_area, data = divar, type = "b", col = 1)
   lines(endemics ~ prop_area, data = divar, type = "b", col = 2)
   lines(ens_shannon ~ prop_area, data = divar, type = "b", col = 3)
   lines(ens_simpson ~ prop_area, data = divar, type = "b", col = 4)

   legend("topleft",legend = c("Species","Endemics",
                               "ENS Shannon (common species) ",
                               "ENS Simpson (dominant species)"),
          col = 1:4, lwd = 2, cex = 0.9)
}


#' -----------------------------------------------------------------------------
#' Local species abundance distribution
#'
#' Get local abundance distribution in rectangle bounded by x0, y0, x0 + xsize,
#' y0 + ysize
#'
#' @param x0 x-coordinate of lower left corner
#' @param y0 x-coordinate of lower left corner
#' @param xsize size of the subplot in x-direction
#' @param ysize size of the subplot in y-direction
#' @param comm \code{\link{community}} object
#'
#' @return Integer vector with local species abundances
#' @export
#'
abund_rect <- function(x0, y0, xsize, ysize, comm)
{
   x <- comm$census$x
   y <- comm$census$y

   # logical vector which trees are in the sampling rectangle
   in_rect <- (x >= x0 & x < (x0+xsize) & y >= y0 & y < (y0+ysize))

   abund <- table(comm$census$species[in_rect])
   return(abund)
}

#'------------------ -----------------------------------------------------------
#' Distance decay
#'
#' Estimate pairwise similarities of abundance distributions of subplots as
#' function of distance
#'
#' @param comm \code{\link{community}} object
#' @param prop_area Subplot size as proportion of the total area
#' @param n_samples Number of randomly located subplots per subplot size
#' @param method Choise of (dis)similarity index. See \code{\link[vegan]{vegdist}}
#' @param binary Perform presence/absence standardization before analysis.
#'
#' @return Dataframe with distances between subplots and the respective similarity
#' indices
#'
#' @examples
#' sim_com1 <- Sim.Thomas.Community(100, 10000)
#' dd1 <- dist_decay(sim_com1)
#' plot(similarity ~ distance, data = dd1)
#' dd_loess <- loess(similarity ~ distance, data = dd1)
#' new_dist <- data.frame(distance = seq(0.01,1,0.01))
#' pred_dd <- predict(dd_loess, newdata = new_dist$distance)
#' lines(new_dist$distance, pred_dd, lwd=2, col = "red")
#'
#'@export
#'
dist_decay <- function(comm, prop_area = 0.05, n_samples = 30,
                       method = "bray", binary = F, plot = F)
{
   if (any(prop_area > 1))
      warning("Subplot areas larger than the community size are ignored!")
   prop_area <- prop_area[prop_area <= 1]

   if (class(comm) != "community")
      stop("divar requires a community object as input. See ?community.")

   dx_plot <- comm$x_min_max[2] - comm$x_min_max[1]
   dy_plot <- comm$y_min_max[2] - comm$y_min_max[1]

   area <- dx_plot * dy_plot * prop_area
   square_size <- sqrt(area)

   xpos <- runif(n_samples, min = comm$x_min_max[1],
                 max = comm$x_min_max[2] - square_size)
   ypos <- runif(n_samples, min = comm$y_min_max[1],
                 max = comm$y_min_max[2] - square_size)

   d <- dist(cbind(xpos,ypos))

   com_tab <- mapply(abund_rect, xpos, ypos,
                     MoreArgs=list(xsize = square_size, ysize = square_size,
                                   comm = comm ))

   similarity <- 1 - vegdist(t(com_tab), method = method, binary = binary)
   similarity[!is.finite(similarity)] <- NA

   dat_out <- data.frame(distance = as.numeric(d),
                         similarity = as.numeric(similarity))

   # order by increasing distance
   dat_out <- dat_out[order(dat_out$distance), ]

   class(dat_out) <- c("dist_decay", "data.frame")

   return(dat_out)
}

#' Plot distance decay object
#'
#' @param dist_decay
#'
#' @return
#' @export
#'
plot.dist_decay <- function(dist_decay)
{
   plot(similarity ~ distance, data = dist_decay, las = 1,
        xlab = "Distance", ylab = "Similarity", main = "Distance decay")
   dd_loess <- loess(similarity ~ distance, data = dist_decay)
   pred_sim <- predict(dd_loess)
   lines(dist_decay$distance, pred_sim, col = "red")
}

# # -----------------------------------------------------------
# abund.rand.rect <- function(prop_area = 0.25, community,
#                             xext = c(0,1), yext=c(0,1))
# {
#    x <- community[,1]
#    y <- community[,2]
#
#    dx_plot <- xext[2] - xext[1]
#    dy_plot <- yext[2] - yext[1]
#
#    area <- dx_plot*dy_plot*prop_area
#    square_size <- sqrt(area)
#
#    if (square_size <= min(c(dx_plot, dy_plot))){
#       dx_rect <- square_size
#       dy_rect <- square_size
#    } else
#    {
#       if (dx_plot >= dy_plot){
#          dx_rect <- dx_plot*prop_area
#          dy_rect <- dy_plot
#       } else {
#          dx_rect <- dx_plot
#          dy_rect <- dy_plot*prop_area
#       }
#    }
#
#    xpos <- runif(1, min = xext[1], max = xext[2] - dx_rect)
#    ypos <- runif(1, min = yext[1], max = yext[2] - dy_rect)
#
#    abund.plots <- mapply(abund.rect, xpos, ypos,
#                          MoreArgs=list(xsize = dx_rect, ysize = dy_rect,
#                                        community = community))
#
#    return(abund.plots[,1])
# }
#




