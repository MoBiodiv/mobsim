#' Get local diversity indices
#'
#' Get diversity indices including species richness, no. of endemics,
#' Shannon and Simpson diversity for one rectangle subplot in the community.
#'
#' @param x0 x-coordinate of lower left corner
#' @param y0 y-coordinate of lower left corner
#' @param xsize Size of the subplot in x-direction
#' @param ysize Size of the subplot in y-direction
#' @param comm \code{\link{community}} object
#'
#' @return Named vector with six diversity indices
#' \enumerate{
#'    \item n_species: Number of species
#'    \item n_endemics: Number of endemics
#'    \item shannon: Shannon index index defined as \eqn{H = - \sum p_i * log(p_i)},
#'    where \eqn{p_i} is the relative abundance of species i:
#'    \item ens_shannon: Effective number of species (ENS) based on the Shannon index exp(H)
#'    \item simpson: Simpson index index (= probability of interspecific encounter PIE)
#'    defined as \eqn{D =  1 - \sum p_i^2}
#'    \item ens_simpson: Effective number of species (ENS) based on the Simpson index \eqn{1/D}
#'
#' }
#'
#' @details The effective number of species is defined as the number of equally
#' abundant species that produce the same value of a certain diversity index
#' as an observed community (Jost 2006). According to Chao et al. 2014 and
#' Chiu et al. 20 ENS_shannon
#' can be interpreted as the number of common species and ENS_simpson as the
#' number of dominant species in a community.
#'
#' @references
#' Jost 2006. Entropy and diversity. Oikos, 113, 363-375.
#'
#' Chao et al. 2014. Rarefaction and extrapolation with Hill numbers: a framework
#' for sampling and estimation in species diversity studies.
#' Ecological Monographs, 84, 45-67.
#'
#' Hsieh et al. 2016. iNEXT: an R package for rarefaction and extrapolation of
#' species diversity (Hill numbers). Methods Ecol Evol, 7, 1451-1456.

#'
#' @examples
#' sim1 <- sim_poisson_community(100,1000)
#' div_rect(0, 0, 0.3, 0.3, sim1)
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
      #simpson <- (n/(n-1)) * (1 - sum(relabund^2))
      simpson <- 1- sum(relabund^2)
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

#' Distribution of local diversity indices
#'
#' Get mean and standard deviation of diversity indices in several equally
#' sized subplots of a community
#'
#' @param prop_area Size of subplots as proportion of the total area
#' @param comm \code{\link{community}} object
#' @param n_rect Number of randomly located subplots
#' @param exclude_zeros Should subplots without individuals be excluded? (logical)
#'
#' @return Vector with mean and standard deviation of the following diversity
#' indices:
#'
#' \enumerate{
#'    \item Number of species
#'    \item Number of endemics
#'    \item Shannon index
#'    \item Effective number of species (ENS) based on Shannon index
#'    \item Simpson index
#'    \item Effective number of species (ENS) based on Simpson index
#' }
#'
#' See the documentation of  \code{\link{div_rect}} for detailed information on the
#' definition of the diversity indices.
#'
#' @examples
#' sim1 <- sim_poisson_community(100,1000)
#' div_rand_rect(prop_area = 0.1, comm = sim1)
#'
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

   xpos <- stats::runif(n_rect, min = comm$x_min_max[1],
                                max = comm$x_min_max[2] - dx_rect)
   ypos <- stats::runif(n_rect, min = comm$y_min_max[1],
                                max = comm$y_min_max[2] - dy_rect)

   div_plots <- mapply(div_rect, xpos, ypos,
                       MoreArgs=list(xsize = dx_rect, ysize = dy_rect,
                                     comm = comm))
   if (exclude_zeros == T)
      div_plots <- div_plots[, div_plots["n_species",] > 0]

   return(c(m_species      = mean(div_plots["n_species",]),
            sd_species     = stats::sd(div_plots["n_species",]),
            m_endemics     = mean(div_plots["n_endemics",]),
            sd_endemics    = stats::sd(div_plots["n_endemics",]),
            m_shannon      = mean(div_plots["shannon",]),
            sd_shannon     = stats::sd(div_plots["shannon",]),
            m_ens_shannon  = mean(div_plots["ens_shannon",]),
            sd_ens_shannon = stats::sd(div_plots["ens_shannon",]),
            m_simpson      = mean(div_plots["simpson",], na.rm = T),
            sd_simpson     = stats::sd(div_plots["simpson",], na.rm = T),
            m_ens_simpson  = mean(div_plots["ens_simpson",], na.rm = T),
            sd_ens_simpson = stats::sd(div_plots["ens_simpson",], na.rm = T)
            )
   )
}

#' Diversity-area relationships
#'
#' Estimate diversity indices in subplots of different sizes. This includes the
#' well-known species-area and endemics-area relationships.
#'
#' @param prop_area Subplot sizes as proportion of the total area (numeric)
#' @param comm \code{\link{community}} object
#' @param n_samples Number of randomly located subplots per subplot size (single integer)
#' @param exclude_zeros Should subplots without individuals be excluded? (logical)
#'
#' @return Dataframe with the proportional area of the subplots and mean and
#' standard deviation of the following diversity indices:
#' \enumerate{
#'    \item Number of species
#'    \item Number of endemics
#'    \item Shannon index
#'    \item Effective number of species (ENS) based on Shannon index
#'    \item Simpson index
#'    \item Effective number of species (ENS) based on Simpson index
#' }
#'
#' See the documentation of  \code{\link{div_rect}} for detailed information on the
#' definition of the diversity indices.
#'
#' @seealso \code{\link{div_rand_rect}}, \code{\link{div_rect}}
#'
#' @examples
#' sim1 <- sim_thomas_community(100, 1000)
#' divar1 <- divar(sim1, prop_area = seq(0.01, 1.0, length = 20))
#' plot(divar1)
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
#' @param x Dataframe generated function \code{\link{divar}}.
#'
#' @param ... Additional graphical parameters used in \code{\link[graphics]{plot}}.

#'
#' @export
plot.divar <- function(x, ...)
{
   #http://stackoverflow.com/questions/8929663/r-legend-placement-in-a-plot

   # #Plot an empty graph and legend to get the size of the legend
   max_spec <- max(x$m_species)
   # plot(m_species ~ prop_area, data = divar, ylim = c(0, max_spec),
   #      type = "n")
   #
   # my_legend_size <-legend("topright",
   #                         legend = c("Species","Endemics", "ENS Shannon ","ENS Simpson"),
   #                         lwd = 2, plot = FALSE)
   # dev.off()

   #custom ylim. Add the height of legend to upper bound of the range
   my_range <- c(0, max_spec)
   #my_range[2] <- 1.05*(my_range[2] + my_legend_size$rect$h)
   my_range[2] <- 1.3 * my_range[2]


   #draw the plot with custom ylim
   graphics::plot(m_species ~ prop_area, data = x, ylim = my_range,
                  type = "n", las = 1,
                  xlab = "Proportion of area sampled", ylab = "No. of species",
                  main = "Diversity-area relationships", ...)
   graphics::lines(m_species ~ prop_area, data = x, type = "b", col = 1)
   graphics::lines(m_endemics ~ prop_area, data = x, type = "b", col = 2)
   graphics::lines(m_ens_shannon ~ prop_area, data = x, type = "b", col = 3)
   graphics::lines(m_ens_simpson ~ prop_area, data = x, type = "b", col = 4)

   graphics::legend("topleft",legend = c("Species","Endemics",
                                         "ENS Shannon ",
                                         "ENS Simpson"),
   col = 1:4, lwd = 2, ncol = 2, cex = 0.9)
}


#' -----------------------------------------------------------------------------
#' Get local species abundance distribution
#'
#' Get local abundance distribution in rectangle bounded by x0, y0, x0 + xsize,
#' y0 + ysize
#'
#' @param x0 x-coordinate of lower left corner
#' @param y0 y-coordinate of lower left corner
#' @param xsize Size of the subplot in x-direction
#' @param ysize Size of the subplot in y-direction
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

#' Distance decay of similarity
#'
#' Estimate pairwise similarities of communities in subplots as
#' function of distance
#'
#' @param comm \code{\link{community}} object
#' @param prop_area Subplot size as proportion of the total area
#' @param n_samples Number of randomly located subplots
#' @param method Choice of (dis)similarity index. See \code{\link[vegan]{vegdist}}
#' @param binary Perform presence/absence standardization before analysis?
#' See \code{\link[vegan]{vegdist}}
#'
#' @return Dataframe with distances between subplot pairs and the respective
#' similarity indices
#'
#' @examples
#' sim_com1 <- sim_thomas_community(100, 10000, sigma = 0.1, mother_points = 2)
#' dd1 <- dist_decay(sim_com1, prop_area = 0.005, n_samples = 20)
#' plot(dd1)
#'
#'@export
#'
dist_decay <- function(comm, prop_area = 0.005, n_samples = 20,
                       method = "bray", binary = F)
{
   if (any(prop_area > 1))
      warning("Subplot areas larger than the community size are ignored!")
   prop_area <- prop_area[prop_area <= 1]

   if (class(comm) != "community")
      stop("divar requires a community object as input. See ?community.")

   dx_plot <- comm$x_min_max[2] - comm$x_min_max[1]
   dy_plot <- comm$y_min_max[2] - comm$y_min_max[1]
   area <- dx_plot * dy_plot * prop_area

   samples1 <- sample_quadrats(comm, n_quadrats = n_samples, quadrat_area = area,
                               avoid_overlap = T, plot = F)
   com_mat <- samples1$spec_dat[rowSums(samples1$spec_dat) > 0,]
   d <- stats::dist(samples1$xy_dat[rowSums(samples1$spec_dat) > 0,])

   similarity <- 1 - vegan::vegdist(com_mat, method = method,
                                    binary = binary)
   similarity[!is.finite(similarity)] <- NA

   dat_out <- data.frame(distance = as.numeric(d),
                         similarity = as.numeric(similarity))

   # order by increasing distance
   dat_out <- dat_out[order(dat_out$distance), ]

   class(dat_out) <- c("dist_decay", "data.frame")

   return(dat_out)
}

#' Plot distance decay of similarity
#'
#' @param x Dataframe generated by \code{\link{dist_decay}}
#'
#' @param ... Additional graphical parameters used in \code{\link[graphics]{plot}}.
#'
#' @details The function plots the similarity indices between all pairs of
#' subplots as function of distance. To indicate the relationship a
#' \code{\link[stats]{loess}} smoother is added to the plot.
#'
#' @examples
#' sim_com1 <- sim_thomas_community(100, 10000)
#' dd1 <- dist_decay(sim_com1)
#' plot(dd1)
#'
#' @export
#'
plot.dist_decay <- function(x, ...)
{
   graphics::plot(similarity ~ distance, data = x, las = 1,
                  xlab = "Distance", ylab = "Similarity",
                  main = "Distance decay", ...)
   dd_loess <- stats::loess(similarity ~ distance, data = x)
   pred_sim <- stats::predict(dd_loess)
   graphics::lines(x$distance, pred_sim, col = "red", lwd = 2)
}

# -----------------------------------------------------------
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




