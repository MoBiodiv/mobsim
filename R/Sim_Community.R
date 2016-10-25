# ----------------------------------------------------------------------------------
#' Simulate log-normal species abundance distributions.
#'
#' Simulate log-normal abundance data in a local community with fixed number of
#' individuals (N.local), number of species in the pool (S.pool), and coefficient
#' of variation of relative abundances (cv.abund).
#'
#' @param S.pool integer - number of species in the pool
#' @param N.local integer - number of individuals in the local community
#' @param cv.abund numeric - coefficient of variation ( = sd/mean) of relative
#' abundances. The higher \code{cv.abund}, the lower the evenness of the
#' simulated community. This means with increasing \code{cv.abund} there are more
#' rare and more dominant species.
#'
#' @param fixS.local logical - should the simulation constrain the number of
#'   species in the local community? This can result in deviations from mean and
#'   sd of local abundances from the theoretical distributions
#'
#' @details When \code{fixS.local = FALSE} the species number in the local
#'   community might deviate from \code{S.pool} due to sampling. When
#'   \code{fix.S.local = TRUE} the local number of species will equal
#'   \code{S.pool}, but this constraint can result in biases from the
#'   theoretical parameters. Therefore the simulated distribution parameters are
#'   provided as model output.
#'
#' @return List with five named items:
#' \enumerate{
#'     \item abund - simulated abundance vector in the local community
#'     \item mean.theor - theoretical mean of the abundance distribution
#'     \item sd.theor - theoretical standard deviation
#'     \item mean.sim - simulated mean
#'     \item sd.sim - simulated standard deviation
#' }
#'
#' @author Felix May
#'
#' @examples
#' SAD.lognorm(S.pool = 100, N.local = 10000, cv.abund = 1)
#'
#' ## Create preston SAD plot with log2 abundance classes
#'
#' require(untb)
#'
#' abund.vec <- SAD.lognorm(S.pool = 100, N.local = 10000, cv.abund = 1)$abund
#' sad1 <- preston(census(abund.vec))
#' barplot(height = as.numeric(sad1), names.arg = names(sad1),
#'         xlab = "Abundance class", ylab ="No. of species")
#'
SAD.lognorm <- function(S.pool, N.local, cv.abund = 1, fixS.local = F)
{
   if (fixS.local == T)
      mean.abund <- N.local/S.pool
   else
      mean.abund <- 100
   sd.abund <- mean.abund*cv.abund

   sigma1 <- sqrt(log(sd.abund^2/mean.abund^2 +1))
   mu1 <- log(mean.abund) - sigma1^2/2

   if (fixS.local == T){

      n <- 0
      while (n < N.local){
         abund1 <- rlnorm(S.pool, meanlog = mu1, sdlog = sigma1)
         abund.local <- round(abund1)
         abund.local[abund.local == 0] <- 1
         n <- sum(abund.local)
      }

      #randomly remove individuals until target level is reached
      while (n > N.local){
         relabund <- abund.local/sum(abund.local)
         irand <- sample(1:S.pool, size=1, prob=relabund) # draw proportional to relative abundance
         if (abund.local[irand]>1) abund.local[irand] <- abund.local[irand]-1
         n <- sum(abund.local)
      }

      abund.local <- sort(abund.local, decreasing = T)
      names(abund.local) <- paste("species", 1:length(abund.local), sep = "")

   } else {

      abund.pool <- rlnorm(S.pool, meanlog = mu1, sdlog = sigma1)
      relabund.pool <- sort(abund.pool/sum(abund.pool), decreasing = T)
      abund.local <- table(sample(1:S.pool, N.local, replace = T, prob = relabund.pool))
      names(abund.local) <- paste("species", names(abund.local), sep = "")

   }

   return(list(abund = abund.local,
               mean.theor = N.local/S.pool,
               sd.theor   = N.local/S.pool * cv.abund,
               mean.sim   = mean(abund.local),
               sd.sim     = sd(abund.local))
   )
}

# ----------------------------------------------------------------------------------
#' Create spatial community object.
#'
#' Creates a spatial community object with a certain extent and with coordinates
#' and species identities of all individuals in the community.
#'
#' @param x numeric - x coordinates
#' @param y numeric - y coordinates
#' @param specID character vector - species names or IDs. Can be integers, characters or factors
#' @param xrange numeric vector of length 2 - extent of the community in x-direction
#' @param yrange numeric vector of length 2 - extent of the community in y-direction
#'
#' @return community object which includes three items:
#' \enumerate{
#'    \item census - data.frame with three columns: X, Y and species names for each individual
#'    \item x_min_max - extent of the community in x-direction
#'    \item y_min_max - extent of the community in y-direction
#' }
#'
#' @examples
#' x <- runif(100)
#' y <- runif(100)
#' species_names <- rep(paste("species",1:10, sep = ""), each = 10)
#'
#' com1 <- community(x,y, species_names)
#' plot(com1)
#' summary(com1)
#'
community <- function(x, y, specID, xrange = c(0,1), yrange = c(0,1))
{
   if (length(xrange) < 2 | length(yrange) < 2) ("Érror: missing ranges for x or y!")

   if (xrange[1] > min(x)) return("Error: Inappropriate ranges for x!")
   if (xrange[2] < max(x)) return("Error: Inappropriate ranges for x!")

   if (yrange[1] > min(y)) return("Error: Inappropriate ranges for y!")
   if (yrange[2] < max(y)) return("Error: Inappropriate ranges for y!")


   points <- data.frame(X = as.numeric(x), Y = as.numeric(y),
                       Species = as.factor(specID))

   com <- list(census = points,
               x_min_max = as.numeric(xrange[1:2]),
               y_min_max = as.numeric(yrange[1:2])
               )

   class(com) <- "community"

   return(com)
}

#' Print summary of spatial community object
summary.community <- function(com1)
{
   cat("No. of individuals: ", nrow(com1$census), "\n")
   cat("No. of species: ", length(unique(com1$census$Species)), "\n")
   cat("x-extent: ", com1$x_min_max, "\n")
   cat("y-extent: ", com1$y_min_max, "\n\n")
   print(summary(com1$census))
}

#' Plot spatial community object
plot.community <- function(com1, col = NA, pch = NA, ...)
{
   nSpec <- length(table(com1$census$Species))
   if (is.na(col))  col <- rainbow(nSpec)
   if (is.na(pch))  pch <- 19

   plot(Y ~ X, data = com1$census, xlim = com1$x_min_max, ylim = com1$y_min_max,
        col = col[com1$census$Species], pch = pch, ...)
}

# ----------------------------------------------------------------------------------
#' Random spatial coordinates
#'
#' Add random spatial positions to a species abundance distribution.
#'
#' @param abund.vec integer vector with species abundances
#' @param xrange numeric vector of length 2 - extent of the community in x-direction
#' @param yrange numeric vector of length 2 - extent of the community in y-direction
#'
#' @return A community object as defined by \code{\link{community}}.
#'
#' @author Felix May
#' @examples
#' abund <- SAD.lognorm(S.pool = 100, N.local = 1000)$abund
#' sim.com1 <- Sim.Poisson.Coords(abund)
#' plot(sim.com1)
#' summary(sim.com1)

Sim.Poisson.Coords <- function(abund.vec,
                               xrange = c(0,1),
                               yrange = c(0,1)
                               )
{
   abund.vec <- trunc(abund.vec)
   if (length(names(abund.vec)) < length(abund.vec))
      names(abund.vec) <- paste("species", 1:length(abund.vec), sep = "")

   N <- sum(abund.vec)
   x <- runif(N, xrange[1], xrange[2])
   y <- runif(N, yrange[1], yrange[2])

   id.spec <- factor(rep(names(abund.vec), times=abund.vec))

   sim.dat1 <- community(x, y, id.spec, xrange, yrange)
   return(sim.dat1)
}

# ----------------------------------------------------------------------------------
#' Simulate community with random spatial positions.
#'
#' This function simulates a community with a certain abundance distribution and
#' and random spatial coordinates. This function consecutively calls
#' \code{\link{SAD.lognorm}} and \code{\link{Sim.Poisson.Coords}}
#'
#' @inheritParams SAD.lognorm
#' @param xrange numeric vector of length 2 - extent of the community in x-direction
#' @param yrange numeric vector of length 2 - extent of the community in y-direction
#'
#' @return A community object as defined by \code{\link{community}}.

#' @author Felix May
#'
#' @examples
#' com1 <- Sim.Poisson.Community(S = 20, N = 500, cv = 1)
#' plot(com1)
#'
Sim.Poisson.Community <- function(S.pool,
                                  N.local,
                                  cv.abund = 1,
                                  fixS.local = F,
                                  xrange= c(0,1),
                                  yrange = c(0,1)
                                  )
{
   sim1 <- SAD.lognorm(S.pool = S.pool, N.local = N.local, cv.abund = cv.abund,
                       fixS.local = fixS.local)
   abund.vec <- sim1$abund

   sim.dat <- Sim.Poisson.Coords(abund.vec = abund.vec,
                                 xrange = xrange, yrange = yrange)
   return(sim.dat)
}


#-----------------------------------------------------------------------------
#' Clumped spatial coordinates
#'
#' Add clumped (aggregated) positions to a species abundance distribution.
#' Clumping is simulated using a Thomas cluster process, also known as Poisson
#' cluster process (Morlon et al. 2008, Wiegand & Moloney 2014)
#'
#' @param abund.vec integer vector with species abundances
#'
#' @param sigma numeric vector of length = 1 or length = 2. Mean displacement
#' (along each coordinate axes) of a point from its mother point (= cluster centre).
#' Therefore sigma correlates with cluster radius. When length(sigma) == 1
#' all species have the same cluster radius.When length(sigma) == 2 a linear
#' relationship between species-specific sigma values and log(relative abundance)
#' is simulated. Thereby sigma[1] is the cluster parameter for the least abundant
#' species and sigma[2] the cluster parameter for the most abundant species.
#' When sigma of any species is more than twice as large as the largest
#' plot dimension, than a random Poisson distribution is simulated, which is more
#' efficient than a Thomas cluster process. The parameter \code{sigma} corresponds to the
#' \code{scale} parameter of \code{\link[spatstat]{rThomas}}.
#'
#' @param mother.points numeric - number of mother points (= cluster centres).
#' If this is a single value, all species have the same number of clusters.
#' For example \code{mother.points = 1} can be used to simulate only one cluster
#' per species, which then represents the complete species range.
#' If \code{mother.points} is a vector of the same length as \code{abund.vec},
#' each species has a specific number of clusters. If no value is provided, the
#' number of clusters is determined from the abundance and the number of points
#' per cluster (\code{cluster.points}).
#'
#' @param cluster.points numeric - mean number of points per cluster. If this is
#' a single value, species have the same average number of points per cluster.
#' If this is a vector of the same length as \code{abund.vec}, each species has
#' a specific mean number of points per cluster.  If no value is provided, the
#' number of points per cluster is determined from the abundance and from
#' \code{mother.points}.  The parameter \code{cluster.points} corresponds to the
#' \code{mu} parameter of \code{\link[spatstat]{rThomas}}.
#'
#' @param xrange numeric vector of length 2 - extent of the community in x-direction
#' @param yrange numeric vector of length 2 - extent of the community in y-direction
#'
#' @details To generate a Thomas cluster process of a single species this
#' function uses a C++ re-implementation if the function
#' \code{\link[spatstat]{rThomas}} in the package \code{spatstat}.
#'
#' There is an inherent link between the parameters \code{abund.vec},
#' \code{mother.points}, and \code{cluster.points}. For every species the
#' abundance has to be equal to the product of the number of clusters
#' (\code{mother.points}) times the numbers of points per cluster
#' (\code{mother.points}).
#'
#' \deqn{abundance = mother.points * cluster.points}
#'
#' Accordingly, if one of the parameters is provided the other one is directly
#' calculated from the abundance. Values for \code{mother.points} override values
#' for \code{cluster.points}. If none of the parameters is specified, it is assumed
#' that for every species there is a similar number of clusters and of points
#' per cluster.
#'
#' \deqn{mother.points = cluster.points = \sqrt(abundance),}
#'
#' In this case rare species have few clusters with few points per
#' cluster, while abundant species have many clusters with many points per cluster.
#'
#' @return A community object as defined by \code{\link{community}}.
#'
#' @references
#' Morlon et al. 2008. A general framework for the distance--decay of similarity
#' in ecological communities. Ecology Letters 11, 904-917.
#'
#' Wiegand and Moloney 2014. Handbook of Spatial Point-Pattern Analysis in Ecology.
#' CRC Press
#'
#' @author Felix May
#'
#' @seealso \code{spatstat::rThomas()}
#'
#' @examples
#'
#' abund <- c(10,20,50,100)
#' sim1 <- Sim.Thomas.Coords(abund, sigma = 0.02)
#' plot(sim1)
#'
#' # Simulate species "ranges"
#' sim2 <- Sim.Thomas.Coords(abund, sigma = 0.02, mother.points = 1)
#' plot(sim2)
#'
#' # Equal numbers of points per cluster
#' sim3 <- Sim.Thomas.Coords(abund, sigma = 0.02, cluster.points = 5)
#' plot(sim3)

Sim.Thomas.Coords <- function(abund.vec,
                              sigma = 0.02,
                              mother.points = NA,
                              cluster.points = NA,
                              xrange = c(0,1),
                              yrange = c(0,1)
                              )
{
   abund.vec <- trunc(abund.vec)
   if (length(names(abund.vec)) < length(abund.vec))
      names(abund.vec) <- paste("species", 1:length(abund.vec), sep = "")

   xext <- xrange[2] - xrange[1]
   yext <- yrange[2] - yrange[1]

   max.dim <- ifelse(xext >= yext, xext, yext)

   cum.abund <- cumsum(abund.vec)
   S.local <- length(abund.vec)
   N <- sum(abund.vec)

   if (length(sigma) == 2){
      # linear relationship between sigma and log(relabund)
      # sigma = a1 + b1 * log(relabund)
      log.relabund <- log(abund.vec/sum(abund.vec))
      range.abund <- max(log.relabund) - min(log.relabund)

      if (range.abund != 0) b1 <- (sigma[2] - sigma[1])/range.abund
      else b1 <- 0

      a1 <- sigma[2] - b1*max(log.relabund)
      sigma.vec <- a1 + b1*log.relabund
   }
   else {
      sigma.vec <- rep(sigma[1], times = S.local)
   }

   X = numeric(N)
   Y = numeric(N)
   id.spec <- factor(rep(names(abund.vec), times=abund.vec))

   # determine points per cluster and number of mother points
   if (is.numeric(mother.points)){

      if (length(mother.points) == S.local)
         n.mother.points <- mother.points
      else
         n.mother.points <- rep(mother.points[1], S.local)

      points.per.cluster <- abund.vec / n.mother.points

   } else {

      if (is.numeric(cluster.points)){

         if (length(cluster.points) == S.local)
             points.per.cluster <- cluster.points
         else
            points.per.cluster <- rep(cluster.points[1], S.local)

         lambda.mother <- abund.vec / points.per.cluster

      } else {
         lambda.mother <- points.per.cluster <- sqrt(abund.vec)
      }
      #n.mother.points <- rpois(S.local, lambda = lambda.mother)
      n.mother.points <- ceiling(lambda.mother)
   }

   # create map for first species
   if (sigma.vec[1] < 2 * max.dim){
      dat1 <- rThomas_rcpp(abund.vec[1],
                           nMotherPoints = n.mother.points[1],
                           sigma = sigma.vec[1],
                           mu = points.per.cluster[1],
                           xmin = xrange[1], xmax = xrange[2],
                           ymin = yrange[1], ymax = yrange[2])
   } else {
      x <- runif(abund.vec[1], xrange[1], xrange[2])
      y <- runif(abund.vec[1], yrange[1], yrange[2])
      dat1 <- data.frame(X = x, Y = y)
   }

   irange <- 1:cum.abund[1]
   X[irange] <- dat1$X
   Y[irange] <- dat1$Y

   for (ispec in 2:S.local){

      if (sigma.vec[ispec] < 2 * max.dim){
         dat1 <- rThomas_rcpp(abund.vec[ispec],
                              nMotherPoints = n.mother.points[ispec],
                              sigma = sigma.vec[ispec],
                              mu = points.per.cluster[ispec],
                              xmin = xrange[1], xmax = xrange[2],
                              ymin = yrange[1], ymax = yrange[2])
      } else {
         x <- runif(abund.vec[ispec], xrange[1], xrange[2])
         y <- runif(abund.vec[ispec], yrange[1], yrange[2])
         dat1 <- data.frame(X = x, Y = y)
      }

      irange <- (cum.abund[ispec-1] + 1):cum.abund[ispec]
      X[irange] <- dat1$X
      Y[irange] <- dat1$Y
   }

   sim.dat1 <- community(X, Y, id.spec, xrange, yrange)
   return(sim.dat1)
}




# -----------------------------------------------------------------------------------------------
#' Simulate community with clumped spatial positions.
#'
#' This function simulates a community with a certain abundance distribution and
#' and clumped, i.e. aggregated, spatial coordinates. This function consecutively calls
#' \code{\link{SAD.lognorm}} and \code{\link{Sim.Thomas.Coords}}
#'
#' @inheritParams SAD.lognorm
#'
#' @param sigma numeric vector of length = 1 or length = 2. Mean displacement
#' (along each coordinate axes) of a point from its mother point (=cluster centre).
#'
#' @param mother.points numeric - number of mother points (= cluster centres).
#' @param cluster.points numeric - mean number of points per cluster.

#' @param xrange numeric vector of length 2 - extent of the community in x-direction
#' @param yrange numeric vector of length 2 - extent of the community in y-direction
#'
#' @return A community object as defined by \code{\link{community}}.
#'
#' @seealso \code{SAD.lognorm}; \code{Sim.Thomas.Coords}

#' @author Felix May
#'
#' @examples
#' com1 <- Sim.Thomas.Community(S.pool = 20, N.local = 500, cv.abund = 1,
#'                              sigma = 0.01)
#' plot(com1)
#'
Sim.Thomas.Community <- function(S.pool,
                                 N.local,
                                 cv.abund = 1,
                                 fixS.local = F,
                                 sigma = 0.02,
                                 cluster.points = NA,
                                 mother.points = NA,
                                 xrange = c(0,1),
                                 yrange = c(0,1)
                                 )
{
   sim1 <- SAD.lognorm(S.pool = S.pool, N.local = N.local, cv.abund = cv.abund,
                       fixS.local = fixS.local)
   abund.vec <- sim1$abund

   sim.dat <- Sim.Thomas.Coords(abund.vec = abund.vec,
                                sigma = sigma,
                                mother.points = mother.points,
                                cluster.points = cluster.points,
                                xrange = xrange,
                                yrange = yrange)

   return(sim.dat)
}


# old version using spatstat function - should do the same but less efficient
# # ----------------------------------------------------------------------------------
# # Simulate community with log-normal SAD and Thomas process clustering
# # when sigma > 2*(max(xmax,ymax)) a Poisson distribution is simulated which is more efficient
# Sim.Thomas.Community <- function(S,N,
#                                  cv.abund=1,
#                                  sigma=0.02, #single value for all species or
#                                              #min and max for least and most abundant species
#                                  xmax=1,
#                                  ymax=1,
#                                  points.cluster=NULL)
# {
#    require(spatstat)
#
#    max.dim <- ifelse(xmax>=ymax,xmax,ymax)
#
#    sim1 <- SAD.log.normal(S,N,cv.abund)
#    abund.vec <- sim1$abund
#
#    if (length(sigma)==2){
#       # linear relationship between sigma and log(relabund)
#       # sigma = a1 + b1 * log(relabund)
#       log.relabund <- log(abund.vec/sum(abund.vec))
#       range.abund <- max(log.relabund)-min(log.relabund)
#
#       if (range.abund != 0) b1 <- (sigma[2]-sigma[1])/range.abund
#       else b1 <- 0
#       a1 <- sigma[2] - b1*max(log.relabund)
#       sigma.vec <- a1 + b1*log.relabund
#
#    #    plot(sigma.vec~log.relabund)
#    #    abline(a1,b1,col="red")
#    #    abline(h=c(sigma[1],sigma[2]),col=1:2)
#    #    abline(v=c(min(log.relabund),max(log.relabund)))
#    }
#    else {
#
#       if (sigma > 2*max.dim){
#          x <- runif(N,0,xmax)
#          y <- runif(N,0,ymax)
#          id.spec <- rep.int(1:S,times=abund.vec)
#
#          dat1 <- data.frame(X=x,Y=y,SpecID=id.spec)
#          return(dat1)
#       }
#
#       sigma.vec <- rep(sigma[1],times=S)
#    }
#
#    print("Test")
#
#    # create map for first species
#    if (!is.numeric(points.cluster)){
#       n.mother.points <- points.per.cluster <- sqrt(abund.vec) # assumption : similar numbers of cluster and of
#                                                                #individuals per cluster
#    }
#    else {
#       points.per.cluster <- rep(points.cluster,S)
#       n.mother.points <- abund.vec/points.per.cluster
#       n.mother.points <- ifelse(n.mother.points<0.1,0.1,n.mother.points)
#    }
#
#    n <- 0
#    #n.trials <- 0
#    while (n<abund.vec[1]){
#       pp1 <- rThomas(kappa=n.mother.points[1]/(xmax*ymax),
#                      scale=sigma.vec[1],
#                      mu=points.per.cluster[1],
#                      win=owin(c(0,xmax),c(0,ymax)))
#       n <- pp1$n
#       #n.trials <- n.trials +1
#    }
#    dat1 <- as.data.frame(pp1)
#    dat1 <- dat1[sample(1:nrow(dat1),abund.vec[1]),]
#    dat1$SpecID <- rep(1,abund.vec[1])
#
#    #plot(y~x,dat1,pch=19,xlim=c(0,1),ylim=c(0,1))
#    #n.trials
#
#    for (ispec in 2:S){
#       n <- 0
#       while (n<abund.vec[ispec]){
#          pp1 <- rThomas(kappa=n.mother.points[ispec]/(xmax*ymax),
#                         scale=sigma.vec[ispec],
#                         mu=points.per.cluster[ispec],
#                         win=owin(c(0,xmax),c(0,ymax)))
#          n <- pp1$n
#       }
#       dat2 <- as.data.frame(pp1)
#       dat2 <- dat2[sample(1:nrow(dat2),abund.vec[ispec]),]
#       dat2$SpecID <- rep(ispec,abund.vec[ispec])
#       dat1 <- rbind(dat1,dat2)
#    }
#
#    dat3 <- dat1
#    dat3$SpecID <- factor(dat1$SpecID)
#
# #    require(ggplot2)
# #    ggplot(data=dat3,aes(x=x,y=y,color=SpecID)) + geom_point(size=4)
#
#    names(dat3) <- c("X","Y","SpecID")
#    return(dat3)
# }




