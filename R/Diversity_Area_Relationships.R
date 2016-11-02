
# -----------------------------------------------------------
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
#'    \item nSpecies - the number of species
#'    \item nEndemics - the number of endemics
#'    \item shannon - Shannon diversity index defined as \eqn{H = - \sum p_{i} * log(p_[i])},
#'    where \eqn{p_i} is the relative abundance of species i:
#'    \item simpson - Simpson diversity index (= probabiliy of interspecific encounter PIE)
#'    defined as \eqn{D =  1 - \sum p_i^2}
#' }

div_rect <- function(x0, y0, xsize, ysize, comm)
{
   x <- comm$census$X
   y <- comm$census$X

   # logical vector which trees are in the sampling rectangle
   in.rect <- (x >= x0 & x < (x0+xsize) & y >= y0 & y < (y0+ysize))

   spec.in <- unique(comm$census$Species[in.rect])
   spec.out <- unique(comm$census$Species[!in.rect])

   nSpecies <- length(spec.in)
   nEndemics <- length(spec.in[!spec.in %in% spec.out])

   abund <- table(comm$census$Species[in.rect])
   abund <- abund[abund > 0]
   rel.abund <- abund/sum(abund)

   shannon <- - sum(rel.abund * log(rel.abund))

   n <- sum(abund)
   if (n > 1)
      simpson <- (n/(n-1)) * (1 - sum(rel.abund^2))
   else
      simpson <- NA

   return(c(nSpecies = nSpecies,
            nEndemics = nEndemics,
            shannon = shannon,
            simpson = simpson))
}

# ------------------------------------------------------------------------------
#' Distribution of local diversity indices
#'
#' Get mean and sd of diversity indices in several equally sized subplots
#' of a community
#'
#' @param prop.A Size of subplots as proportion of the total area
#' @param comm \code{\link{community}} object
#' @param nrect Number of randomly located subplots
#' @param exclude_zeros logical - should subplots without individuals be excluded?
#'
#' @return Vector with mean and standard deviation of the following diversity
#' indices: (1) Number of species (2) Number of endemics (3) Shannon-diversity
#' (4) Simpson diversity
#'
#' @seealso \code{\link{div_rect}}
#'
div_rand_rect <- function(prop.A = 0.25, comm, nrect = 100, exclude_zeros = F)
{
   dx.plot <- comm$x_min_max[2] - comm$x_min_max[1]
   dy.plot <- comm$y_min_max[2] - comm$y_min_max[1]

   area <- dx.plot * dy.plot * prop.A
   square.size <- sqrt(area)

   if (square.size <= min(c(dx.plot, dy.plot))){
      dx.rect <- square.size
      dy.rect <- square.size
   } else
   {
      if (dx.plot >= dy.plot){
         dx.rect <- dx.plot*prop.A
         dy.rect <- dy.plot
      } else {
         dx.rect <- dx.plot
         dy.rect <- dy.plot*prop.A
      }
   }

   xpos <- runif(nrect, min = comm$x_min_max[1], max = comm$x_min_max[2] - dx.rect)
   ypos <- runif(nrect, min = comm$y_min_max[1], max = comm$y_min_max[2] - dy.rect)

   div_plots <- mapply(div_rect, xpos, ypos,
                       MoreArgs=list(xsize = dx.rect, ysize = dy.rect,
                                     comm = comm))
   if (exclude_zeros == T)
      div_plots <- div_plots[, div_plots["nSpecies",] > 0]

   return(c(meanSpec    = mean(div_plots["nSpecies",]),
            sdSpec      = sd(div_plots["nSpecies",]),
            meanEnd     = mean(div_plots["nEndemics",]),
            sdEnd       = sd(div_plots["nEndemics",]),
            meanShannon = mean(div_plots["shannon",]),
            sdShannon   = sd(div_plots["shannon",]),
            meanSimpson = mean(div_plots["simpson",], na.rm = T),
            sdSimpson   = sd(div_plots["simpson",], na.rm = T))
          )
}

#' Diversity-area relationships
#'
#' Estimate diversity indices in subplots of different sizes. This includes the
#' well-known species-area and endemics-area relationships.
#'
#' @param prop.A numeric vector with subplot sizes as proportion of the total area
#' @param comm \code{\link{community}} object
#' @param nsamples Number of randomly located subplots per subplot size
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
#' divar1 <- DivAR(sim_com1, prop.A = seq(0.01, 1.0, length = 20))
#' plot(meanSpec ~ propArea, data = divar1, xlab = "Proportion of area",
#'      ylab = "No. of species", type = "b", ylim = c(0,100))
#' points(meanEnd ~ propArea, data = divar1, type = "b", col = 2)
DivAR <- function(comm, prop.A = seq(0.1, 1, by = 0.1), nsamples = 100,
                  exclude_zeros = F)
{
   if (any(prop.A > 1))
      warning("Subplot areas larger than the community size are ignored!")
   prop.A <- prop.A[prop.A <= 1]

   if (class(comm) != "community")
      stop("DiVAR requires a community object as input. See ?community.")

   nscales <- length(prop.A)
   dx.plot <- comm$x_min_max[2] - comm$x_min_max[1]
   dy.plot <- comm$y_min_max[2] - comm$y_min_max[1]

   div_area <- sapply(prop.A,
                      div_rand_rect,
                      comm = comm,
                      nrect = nsamples,
                      exclude_zeros = exclude_zeros)

   div_dat <- as.data.frame(t(div_area))
   div_dat <- cbind(propArea = prop.A, div_dat)

   return(div_dat)
}


# -----------------------------------------------------------
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

abund_rect <- function(x0, y0, xsize, ysize, comm)
{
   x <- comm$census$X
   y <- comm$census$Y

   # logical vector which trees are in the sampling rectangle
   in.rect <- (x >= x0 & x < (x0+xsize) & y >= y0 & y < (y0+ysize))

   abund <- table(comm$census$Species[in.rect])
   return(abund)
}

# -----------------------------------------------------------
#' Distance decay
#'
#' Estimate pairwise similarities of abundance distributions of subplots as
#' function of distance
#'
#' @param comm \code{\link{community}} object
#' @param prop.A Subplot size as proportion of the total area
#' @param nsamples Number of randomly located subplots per subplot size
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
dist_decay <- function(comm, prop.A = 0.05, nsamples = 30,
                       method = "bray", binary = F)
{
   require(vegan)

   if (any(prop.A > 1))
      warning("Subplot areas larger than the community size are ignored!")
   prop.A <- prop.A[prop.A <= 1]

   if (class(comm) != "community")
      stop("DiVAR requires a community object as input. See ?community.")

   dx.plot <- comm$x_min_max[2] - comm$x_min_max[1]
   dy.plot <- comm$y_min_max[2] - comm$y_min_max[1]

   area <- dx.plot * dy.plot * prop.A
   square.size <- sqrt(area)

   xpos <- runif(nsamples, min = comm$x_min_max[1], max = comm$x_min_max[2] - square.size)
   ypos <- runif(nsamples, min = comm$y_min_max[1], max = comm$y_min_max[2] - square.size)

   d <- dist(cbind(xpos,ypos))

   com.tab <- mapply(abund_rect, xpos, ypos,
                     MoreArgs=list(xsize = square.size, ysize = square.size,
                                   comm = comm ))

   similarity <- 1 - vegdist(t(com.tab), method = method, binary = binary)
   similarity[!is.finite(similarity)] <- NA

   dat_out <- data.frame(distance = as.numeric(d),
                         similarity = as.numeric(similarity))

   # order by increasing distance
   dat_out <- dat_out[order(dat_out$distance), ]

   return(dat_out)
}

# # -----------------------------------------------------------
# abund.rand.rect <- function(prop.A = 0.25, community,
#                             xext = c(0,1), yext=c(0,1))
# {
#    x <- community[,1]
#    y <- community[,2]
#
#    dx.plot <- xext[2] - xext[1]
#    dy.plot <- yext[2] - yext[1]
#
#    area <- dx.plot*dy.plot*prop.A
#    square.size <- sqrt(area)
#
#    if (square.size <= min(c(dx.plot, dy.plot))){
#       dx.rect <- square.size
#       dy.rect <- square.size
#    } else
#    {
#       if (dx.plot >= dy.plot){
#          dx.rect <- dx.plot*prop.A
#          dy.rect <- dy.plot
#       } else {
#          dx.rect <- dx.plot
#          dy.rect <- dy.plot*prop.A
#       }
#    }
#
#    xpos <- runif(1, min = xext[1], max = xext[2] - dx.rect)
#    ypos <- runif(1, min = yext[1], max = yext[2] - dy.rect)
#
#    abund.plots <- mapply(abund.rect, xpos, ypos,
#                          MoreArgs=list(xsize = dx.rect, ysize = dy.rect,
#                                        community = community))
#
#    return(abund.plots[,1])
# }
#


