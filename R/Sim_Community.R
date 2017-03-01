#' Simulate log-normal species abundance distributions.
#'
#' Simulate log-normal abundance data in a local community with fixed number of
#' individuals (n_sim), number of species in the pool (s_pool), and coefficient
#' of variation of relative abundances (cv_abund).
#'
#' @param s_pool integer - number of species in the pool
#' @param n_sim integer - number of individuals in the simulate local community
#' @param cv_abund numeric - coefficient of variation ( = sd/mean) of relative
#' abundances. The higher \code{cv_abund}, the lower the evenness of the
#' simulated community. This means with increasing \code{cv_abund} there are more
#' rare and more dominant species.
#'
#' @param fix_s_sim logical - should the simulation constrain the number of
#'   species in the local community? This can result in deviations from mean and
#'   sd of local abundances from the theoretical distributions
#'
#' @details When \code{fix_s_sim = FALSE} the species number in the local
#'   community might deviate from \code{s_pool} due to sampling. When
#'   \code{fix_s_sim = TRUE} the local number of species will equal
#'   \code{s_pool}, but this constraint can result in biases from the
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
#' abund1 <- sim_sad(s_pool = 100, n_sim = 10000, cv_abund = 1)
#' plot(abund1, method = "preston")
#' plot(abund1, method = "rank")
#'
#' @export
#'
sim_sad <- function(s_pool, n_sim, cv_abund = 1, fix_s_sim = F)
{
   if (!is.numeric(s_pool) || s_pool <= 0)
      stop("s_pool has to be a positive integer number")
   if (!is.numeric(n_sim) || n_sim <= 0)
      stop("n_sim has to be a positive integer number")

   if (s_pool %% as.integer(s_pool) > 0)
      warning("s_pool is rounded to the nearest integer")
   if (n_sim %% as.integer(n_sim) > 0)
      warning("n_sim is rounded to the nearest integer")

   s_pool <- round(s_pool, digits = 0)
   n_sim <- round(n_sim, digits = 0)

   if (fix_s_sim == T)
      mean_abund <- n_sim/s_pool
   else
      mean_abund <- 100

   sd_abund <- mean_abund * cv_abund

   sigma1 <- sqrt(log(sd_abund^2/mean_abund^2 +1))
   mu1 <- log(mean_abund) - sigma1^2/2

   if (fix_s_sim == T){

      n <- 0
      while (n < n_sim){
         abund1 <- rlnorm(s_pool, meanlog = mu1, sdlog = sigma1)
         abund_local <- round(abund1)
         abund_local[abund_local == 0] <- 1
         n <- sum(abund_local)
      }

      #randomly remove individuals until target level is reached
      while (n > n_sim){
         relabund <- abund_local/sum(abund_local)
         # draw proportional to relative abundance
         irand <- sample(1:s_pool, size = 1, prob = relabund)
         if (abund_local[irand] > 1) abund_local[irand] <- abund_local[irand]-1
         n <- sum(abund_local)
      }

      abund_local <- sort(abund_local, decreasing = T)
      names(abund_local) <- paste("species", 1:length(abund_local), sep = "")

   } else {

      abund_pool <- rlnorm(s_pool, meanlog = mu1, sdlog = sigma1)
      relabund_pool <- sort(abund_pool/sum(abund_pool), decreasing = T)
      abund_local <- table(sample(1:s_pool, n_sim, replace = T,
                                  prob = relabund_pool))
      names(abund_local) <- paste("species", names(abund_local), sep = "")

   }

   # return(list(abund = abund.local,
   #             mean.theor = n_sim/s_pool,
   #             sd.theor   = n_sim/s_pool * cv_abund,
   #             mean.sim   = mean(abund.local),
   #             sd.sim     = sd(abund.local))
   # )

   class(abund_local) <- "sad"
   return(abund_local)
}



#' Plot species abundance distribution
#'
#' @param abund
#'
#' @param method
#'
#' @export
plot.sad <- function(abund, method = c("octave","rank"))
{
   method <- match.arg(method)

   if (method == "rank")
      plot(sort(as.numeric(abund), decreasing = T), type="b", log="y", las = 1,
           xlab="Species rank", ylab="Species abundance",
           main = "Rank-abundance curve", las = 1)

   if (method == "octave"){

      # code adopted from untb:preston()
      max_abund <- max(abund)
      n <- 1 + ceiling(log(max_abund)/log(2)) # number of abundance classes

      if (n < 2) breaks <- c(0, 1)
      else       breaks <- c(0, 2^(0:(n - 2)), max_abund)

      r <- hist(abund, plot = FALSE, breaks = breaks, right = TRUE)
      abund_dist <- r$counts

      if (n <= 2) names(abund_dist) <- c("1","2")[1:n]
      else        names(abund_dist) <- c("1", "2",
                                         paste(breaks[-c(1:2, length(breaks))] + 1,
                                          "-", breaks[-c(1:3)], sep = ""))

      barplot(height = as.numeric(abund_dist), names.arg = names(abund_dist),
              xlab = "Abundance class", ylab ="No. of species",
              main = "Preston octave plot", las = 1)
   }
}

#' Create spatial community object.
#'
#' Creates a spatial community object with a certain extent and with coordinates
#' and species identities of all individuals in the community.
#'
#' @param x numeric - x coordinates
#' @param y numeric - y coordinates
#' @param spec_id character vector - species names or IDs. Can be integers, characters or factors
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
#' @export
#'
community <- function(x, y, spec_id, xrange = c(0,1), yrange = c(0,1))
{
   if (length(xrange) < 2 | length(yrange) < 2) ("Érror: missing ranges for x or y!")

   if (xrange[1] > min(x)) return("Error: Inappropriate ranges for x!")
   if (xrange[2] < max(x)) return("Error: Inappropriate ranges for x!")

   if (yrange[1] > min(y)) return("Error: Inappropriate ranges for y!")
   if (yrange[2] < max(y)) return("Error: Inappropriate ranges for y!")


   points <- data.frame(x = as.numeric(x), y = as.numeric(y),
                        species = as.factor(spec_id))

   com <- list(census = points,
               x_min_max = as.numeric(xrange[1:2]),
               y_min_max = as.numeric(yrange[1:2])
               )

   class(com) <- "community"

   return(com)
}

#' Print summary of spatial community object
#'
#' @param com1
#'
#' @export
#'
summary.community <- function(com1)
{
   cat("No. of individuals: ", nrow(com1$census), "\n")
   cat("No. of species: ", length(unique(com1$census$species)), "\n")
   cat("x-extent: ", com1$x_min_max, "\n")
   cat("y-extent: ", com1$y_min_max, "\n\n")
   print(summary(com1$census))
}

#' Plot spatial community object
#'
#' @param com1
#' @param col
#' @param pch
#' @param pattern
#' @param ...
#'
#' @export
#'
plot.community <- function(com1, col = NA, pch = NA, ...)
{
   nspec <- length(table(com1$census$species))
   if (is.na(col))  col <- rainbow(nspec)
   if (is.na(pch))  pch <- 19

   plot(y ~ x, data = com1$census, xlim = com1$x_min_max, ylim = com1$y_min_max,
        col = col[com1$census$species], pch = pch, las = 1, ...)

}

#' Get species abundance distribution from community object
#'
#' @param com
#'
#' @export
#'
community_to_sad <- function(com)
{
   if (class(com) != "community")
      stop("community_to_sad requires a community object as input. See ?community.")

   abund <- table(com$census$species)
   class(abund) <- "sad"

   return(abund)
}


#' Random spatial coordinates
#'
#' Add random spatial positions to a species abundance distribution.
#'
#' @param abund_vec integer vector with species abundances
#' @param xrange numeric vector of length 2 - extent of the community in x-direction
#' @param yrange numeric vector of length 2 - extent of the community in y-direction
#'
#' @return A community object as defined by \code{\link{community}}.
#'
#' @author Felix May
#' @examples
#' abund <- sim_sad(s_pool = 100, n_sim = 1000)
#' sim_com1 <- sim_poisson_coords(abund)
#' plot(sim_com1)
#' summary(sim_com1)
#'
#' @export
#'
sim_poisson_coords <- function(abund_vec,
                               xrange = c(0,1),
                               yrange = c(0,1)
                               )
{
   abund_vec <- trunc(abund_vec)
   if (length(names(abund_vec)) < length(abund_vec))
      names(abund_vec) <- paste("species", 1:length(abund_vec), sep = "")

   n <- sum(abund_vec)
   x <- runif(n, xrange[1], xrange[2])
   y <- runif(n, yrange[1], yrange[2])

   id_spec <- factor(rep(names(abund_vec), times = abund_vec))

   sim_dat1 <- community(x, y, id_spec, xrange, yrange)
   return(sim_dat1)
}

#' Simulate community with random spatial positions.
#'
#' This function simulates a community with a certain abundance distribution and
#' and random spatial coordinates. This function consecutively calls
#' \code{\link{sim_sad}} and \code{\link{sim_poisson_coords}}
#'
#' @inheritParams sim_sad
#' @param xrange numeric vector of length 2 - extent of the community in x-direction
#' @param yrange numeric vector of length 2 - extent of the community in y-direction
#'
#' @return A community object as defined by \code{\link{community}}.

#' @author Felix May
#'
#' @examples
#' com1 <- sim_poisson_community(S = 20, N = 500, cv = 1)
#' plot(com1)
#'
#' @export
#'
sim_poisson_community <- function(s_pool,
                                  n_sim,
                                  cv_abund = 1,
                                  fix_s_sim = F,
                                  xrange= c(0,1),
                                  yrange = c(0,1)
                                  )
{
   sim1 <- sim_sad(s_pool = s_pool, n_sim = n_sim, cv_abund = cv_abund,
                       fix_s_sim = fix_s_sim)
   abund_vec <- sim1

   sim_dat <- sim_poisson_coords(abund_vec = abund_vec,
                                 xrange = xrange, yrange = yrange)
   return(sim_dat)
}


#' Clumped spatial coordinates
#'
#' Add clumped (aggregated) positions to a species abundance distribution.
#' Clumping is simulated using a Thomas cluster process, also known as Poisson
#' cluster process (Morlon et al. 2008, Wiegand & Moloney 2014)
#'
#' @param abund_vec integer vector with species abundances
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
#' @param mother_points numeric - number of mother points (= cluster centres).
#' If this is a single value, all species have the same number of clusters.
#' For example \code{mother_points = 1} can be used to simulate only one cluster
#' per species, which then represents the complete species range.
#' If \code{mother_points} is a vector of the same length as \code{abund_vec},
#' each species has a specific number of clusters. If no value is provided, the
#' number of clusters is determined from the abundance and the number of points
#' per cluster (\code{cluster_points}).
#'
#' @param cluster_points numeric - mean number of points per cluster. If this is
#' a single value, species have the same average number of points per cluster.
#' If this is a vector of the same length as \code{abund_vec}, each species has
#' a specific mean number of points per cluster.  If no value is provided, the
#' number of points per cluster is determined from the abundance and from
#' \code{mother_points}.  The parameter \code{cluster_points} corresponds to the
#' \code{mu} parameter of \code{\link[spatstat]{rThomas}}.
#'
#' @param xrange numeric vector of length 2 - extent of the community in x-direction
#' @param yrange numeric vector of length 2 - extent of the community in y-direction
#'
#' @details To generate a Thomas cluster process of a single species this
#' function uses a C++ re-implementation if the function
#' \code{\link[spatstat]{rThomas}} in the package \code{spatstat}.
#'
#' There is an inherent link between the parameters \code{abund_vec},
#' \code{mother_points}, and \code{cluster_points}. For every species the
#' abundance has to be equal to the product of the number of clusters
#' (\code{mother_points}) times the numbers of points per cluster
#' (\code{mother_points}).
#'
#' \deqn{abundance = mother_points * cluster_points}
#'
#' Accordingly, if one of the parameters is provided the other one is directly
#' calculated from the abundance. Values for \code{mother_points} override values
#' for \code{cluster_points}. If none of the parameters is specified, it is assumed
#' that for every species there is a similar number of clusters and of points
#' per cluster.
#'
#' \deqn{mother_points = cluster_points = \sqrt(abundance),}
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
#' sim1 <- sim_thomas_coords(abund, sigma = 0.02)
#' plot(sim1)
#'
#' # Simulate species "ranges"
#' sim2 <- sim_thomas_coords(abund, sigma = 0.02, mother_points = 1)
#' plot(sim2)
#'
#' # Equal numbers of points per cluster
#' sim3 <- sim_thomas_coords(abund, sigma = 0.02, cluster_points = 5)
#' plot(sim3)
#'
#' @export
#'
sim_thomas_coords <- function(abund_vec,
                              sigma = 0.02,
                              mother_points = NA,
                              cluster_points = NA,
                              xrange = c(0,1),
                              yrange = c(0,1)
                              )
{
   abund_vec <- trunc(abund_vec)
   if (length(names(abund_vec)) < length(abund_vec))
      names(abund_vec) <- paste("species", 1:length(abund_vec), sep = "")

   xext <- xrange[2] - xrange[1]
   yext <- yrange[2] - yrange[1]

   max_dim <- ifelse(xext >= yext, xext, yext)

   cum_abund <- cumsum(abund_vec)
   s_local <- length(abund_vec)
   n <- sum(abund_vec)

   # if (length(sigma) == 2){
   #    # linear relationship between sigma and log(relabund)
   #    # sigma = a1 + b1 * log(relabund)
   #    log_relabund <- log(abund_vec/sum(abund_vec))
   #    range_abund <- max(log_relabund) - min(log_relabund)
   #
   #    if (range_abund != 0) b1 <- (sigma[2] - sigma[1])/range_abund
   #    else b1 <- 0
   #
   #    a1 <- sigma[2] - b1*max(log_relabund)
   #    sigma_vec <- a1 + b1*log_relabund
   # }

   if (length(sigma) == s_local){
      sigma_vec <- sigma
   } else {
      sigma_vec <- rep(sigma[1], times = s_local)
   }

   x = numeric(n)
   y = numeric(n)
   id_spec <- factor(rep(names(abund_vec), times = abund_vec))

   # determine points per cluster and number of mother points
   if (!is.na(mother_points)){

      if (length(mother_points) == s_local)
         n_mothers <- mother_points
      else
         n_mothers <- rep(mother_points[1], s_local)

      points_per_cluster <- abund_vec / n_mothers

   } else {

      if (!is.na(cluster_points)){

         if (length(cluster_points) == s_local)
             points_per_cluster <- cluster_points
         else
            points_per_cluster <- rep(cluster_points[1], s_local)

         lambda_mother <- abund_vec / points_per_cluster

      } else {
         lambda_mother <- points_per_cluster <- sqrt(abund_vec)
      }
      #n.mother_points <- rpois(s_local, lambda = lambda_mother)
      n_mothers <- ceiling(lambda_mother)
   }

   # create map for first species
   if (sigma_vec[1] < 2 * max_dim){
      dat1 <- rThomas_rcpp(abund_vec[1],
                           n_mother_points = n_mothers[1],
                           sigma = sigma_vec[1],
                           mu = points_per_cluster[1],
                           xmin = xrange[1], xmax = xrange[2],
                           ymin = yrange[1], ymax = yrange[2])
   } else {
      x1 <- runif(abund_vec[1], xrange[1], xrange[2])
      y1 <- runif(abund_vec[1], yrange[1], yrange[2])
      dat1 <- data.frame(x = x1, y = y1)
   }

   irange <- 1:cum_abund[1]
   x[irange] <- dat1$x
   y[irange] <- dat1$y

   if (s_local > 1){
      for (ispec in 2:s_local){

         if (sigma_vec[ispec] < 2 * max_dim){
            dat1 <- rThomas_rcpp(abund_vec[ispec],
                                 n_mother_points = n_mothers[ispec],
                                 sigma = sigma_vec[ispec],
                                 mu = points_per_cluster[ispec],
                                 xmin = xrange[1], xmax = xrange[2],
                                 ymin = yrange[1], ymax = yrange[2])
         } else {
            x1 <- runif(abund_vec[ispec], xrange[1], xrange[2])
            y1 <- runif(abund_vec[ispec], yrange[1], yrange[2])
            dat1 <- data.frame(x = x1, y = y1)
         }

         irange <- (cum_abund[ispec-1] + 1):cum_abund[ispec]
         x[irange] <- dat1$x
         y[irange] <- dat1$y
      }
   }

   sim_dat1 <- community(x, y, id_spec, xrange, yrange)
   return(sim_dat1)
}

#' Simulate community with clumped spatial positions.
#'
#' This function simulates a community with a certain abundance distribution and
#' and clumped, i.e. aggregated, spatial coordinates. This function consecutively calls
#' \code{\link{sim_sad}} and \code{\link{sim_thomas_coords}}
#'
#' @inheritParams sim_sad
#'
#' @param sigma numeric vector of length = 1 or length = 2. Mean displacement
#' (along each coordinate axes) of a point from its mother point (=cluster centre).
#'
#' @param mother_points numeric - number of mother points (= cluster centres).
#' @param cluster_points numeric - mean number of points per cluster.

#' @param xrange numeric vector of length 2 - extent of the community in x-direction
#' @param yrange numeric vector of length 2 - extent of the community in y-direction
#'
#' @return A community object as defined by \code{\link{community}}.
#'
#' @seealso \code{sim_sad}; \code{sim_thomas_coords}

#' @author Felix May
#'
#' @examples
#' com1 <- sim_thomas_community(s_pool = 20, n_sim = 500, cv_abund = 1,
#'                              sigma = 0.01)
#' plot(com1)
#'
#' @export
#'
sim_thomas_community <- function(s_pool,
                                 n_sim,
                                 cv_abund = 1,
                                 fix_s_sim = F,
                                 sigma = 0.02,
                                 cluster_points = NA,
                                 mother_points = NA,
                                 xrange = c(0,1),
                                 yrange = c(0,1)
                                 )
{
   sim1 <- sim_sad(s_pool = s_pool, n_sim = n_sim, cv_abund = cv_abund,
                       fix_s_sim = fix_s_sim)
   abund_vec <- sim1

   sim_dat <- sim_thomas_coords(abund_vec = abund_vec,
                                sigma = sigma,
                                mother_points = mother_points,
                                cluster_points = cluster_points,
                                xrange = xrange,
                                yrange = yrange)

   return(sim_dat)
}


# old version using spatstat function - should do the same but less efficient
# # ----------------------------------------------------------------------------------
# # Simulate community with log-normal SAD and Thomas process clustering
# # when sigma > 2*(max(xmax,ymax)) a Poisson distribution is simulated which is more efficient
# sim_thomas_community <- function(S,N,
#                                  cv_abund=1,
#                                  sigma=0.02, #single value for all species or
#                                              #min and max for least and most abundant species
#                                  xmax=1,
#                                  ymax=1,
#                                  points.cluster=NULL)
# {
#    require(spatstat)
#
#    max_dim <- ifelse(xmax>=ymax,xmax,ymax)
#
#    sim1 <- SAD.log.normal(S,N,cv_abund)
#    abund_vec <- sim1$abund
#
#    if (length(sigma)==2){
#       # linear relationship between sigma and log(relabund)
#       # sigma = a1 + b1 * log(relabund)
#       log_relabund <- log(abund_vec/sum(abund_vec))
#       range_abund <- max(log_relabund)-min(log_relabund)
#
#       if (range_abund != 0) b1 <- (sigma[2]-sigma[1])/range_abund
#       else b1 <- 0
#       a1 <- sigma[2] - b1*max(log_relabund)
#       sigma_vec <- a1 + b1*log_relabund
#
#    #    plot(sigma_vec~log_relabund)
#    #    abline(a1,b1,col="red")
#    #    abline(h=c(sigma[1],sigma[2]),col=1:2)
#    #    abline(v=c(min(log_relabund),max(log_relabund)))
#    }
#    else {
#
#       if (sigma > 2*max_dim){
#          x <- runif(N,0,xmax)
#          y <- runif(N,0,ymax)
#          id_spec <- rep.int(1:S,times=abund_vec)
#
#          dat1 <- data.frame(X=x,Y=y,spec_id=id_spec)
#          return(dat1)
#       }
#
#       sigma_vec <- rep(sigma[1],times=S)
#    }
#
#    print("Test")
#
#    # create map for first species
#    if (!is.numeric(points.cluster)){
#       n.mother_points <- points_per_cluster <- sqrt(abund_vec) # assumption : similar numbers of cluster and of
#                                                                #individuals per cluster
#    }
#    else {
#       points_per_cluster <- rep(points.cluster,S)
#       n.mother_points <- abund_vec/points_per_cluster
#       n.mother_points <- ifelse(n.mother_points<0.1,0.1,n.mother_points)
#    }
#
#    n <- 0
#    #n.trials <- 0
#    while (n<abund_vec[1]){
#       pp1 <- rThomas(kappa=n.mother_points[1]/(xmax*ymax),
#                      scale=sigma_vec[1],
#                      mu=points_per_cluster[1],
#                      win=owin(c(0,xmax),c(0,ymax)))
#       n <- pp1$n
#       #n.trials <- n.trials +1
#    }
#    dat1 <- as.data.frame(pp1)
#    dat1 <- dat1[sample(1:nrow(dat1),abund_vec[1]),]
#    dat1$spec_id <- rep(1,abund_vec[1])
#
#    #plot(y~x,dat1,pch=19,xlim=c(0,1),ylim=c(0,1))
#    #n.trials
#
#    for (ispec in 2:S){
#       n <- 0
#       while (n<abund_vec[ispec]){
#          pp1 <- rThomas(kappa=n.mother_points[ispec]/(xmax*ymax),
#                         scale=sigma_vec[ispec],
#                         mu=points_per_cluster[ispec],
#                         win=owin(c(0,xmax),c(0,ymax)))
#          n <- pp1$n
#       }
#       dat2 <- as.data.frame(pp1)
#       dat2 <- dat2[sample(1:nrow(dat2),abund_vec[ispec]),]
#       dat2$spec_id <- rep(ispec,abund_vec[ispec])
#       dat1 <- rbind(dat1,dat2)
#    }
#
#    dat3 <- dat1
#    dat3$spec_id <- factor(dat1$spec_id)
#
# #    require(ggplot2)
# #    ggplot(data=dat3,aes(x=x,y=y,color=spec_id)) + geom_point(size=4)
#
#    names(dat3) <- c("X","Y","spec_id")
#    return(dat3)
# }




