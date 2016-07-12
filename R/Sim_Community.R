#' Simulate log-normal species abundance distributions
#'
#' Simulate log-normal abundance data in a local community with fixed number of
#' individuals (N.local), number of species in the pool (S.pool) and
#' coefficient of variation in abundance (cv)
#'
#' @param S.pool single integer - number of species in the pool
#' @param N.local single integer - number of individuals in the local community
#' @param cv.abund single numeric - coefficient of variation ( = sd/mean) in abundance
#' @param fix.S.local - logical - should the simulation constrain the number of
#' species in the local community. This can result in deviations from
#' mean and sd of local abundances from the theoretical distributions
#'
#' @details The function can precisely control for the number of species and
#' the number of individuals in the simulated community. These constraints
#' can result in biases from the theoretical parameters. Therefore the simulated
#' distribution parameters are provided as model output.
#'
#' @return List with five named items:
#' \enumerate{
#'    \item abund - simulated abundance vector in the local community
#'    \item mean.theor - theoretical mean of the distribution
#'    \item sd.theor - theoretical standard deviation
#'    \item mean.sim - simulated mean
#'    \item sd.sim - simulated standard deviation
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
#' names(abund.vec) <- paste("spec", 1:length(abund.vec),sep="")
#' sad1 <- preston(census(abund.vec))
#' barplot(height = as.numeric(sad1), names.arg = names(sad1),
#'         xlab = "Abundance class", ylab ="No. of species")

SAD.lognorm <- function(S.pool, N.local, mean.abund = 100, cv.abund = 1, fix.S.local = F)
{
   if (fix.S.local == T)
      mean.abund <- N.local/S.pool
   sd.abund <- mean.abund*cv.abund

   sigma1 <- sqrt(log(sd.abund^2/mean.abund^2 +1))
   mu1 <- log(mean.abund) - sigma1^2/2

   if (fix.S.local == T){

      n <- 0
      while (n < N.local){
         abund1 <- rlnorm(S.pool, meanlog=mu1, sdlog=sigma1)
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

      abund.pool <- rlnorm(S.pool, meanlog=mu1, sdlog=sigma1)
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
#' Simulate community with log-normal SAD and random spatial positions
#'
#' @param S scalar integer - number of species
#' @param N scalar integer - number of individuals
#' @param cv scalar numeric - coefficient of variation ( = sd/mean) in abundance
#' @param xmax scalar numeric - size of the simulated area in x-direction
#' @param ymax scalar numeric - size of the simulated area in y-direction
#'
#' @return A dataframe with three columns
#' \enumerate{
#'    \item X : x-coordinate of individual
#'    \item Y : y-coordinate of individual
#'    \item SpecID : factor with species identity label (from 1 to number of species)
#' }

#' @author Felix May
#'
#' @examples
#' community1 <- Sim.Poisson.Community(S = 20, N = 500, cv = 1)
#'
#' spec.cols <- rainbow(length(levels(community1$SpecID)))
#' plot(Y ~ X, data =community1, col = spec.cols[community1$SpecID], pch = 19)
#'
Sim.Poisson.Community <- function(S, N, cv.abund = 1, xmax = 1, ymax = 1, fix.S = F)
{
   sim1 <- SAD.lognorm(S.pool = S, N.local = N, cv.abund = cv.abund, fix.S.local = fix.S)
   abund.vec <- sim1$abund

   x <- runif(N, 0, xmax)
   y <- runif(N, 0, ymax)

   id.spec <- factor(rep(names(abund.vec), times=abund.vec))

   dat1 <- data.frame(X=x, Y=y, SpecID=id.spec)
   return(dat1)
}

#' Simulate community with log-normal SAD and Thomas process clustering
#'

#' @param S scalar integer - number of species
#' @param N scalar integer - number of individuals
#' @param cv scalar numeric - coefficient of variation ( = sd/mean) in abundance
#' @param sigma numeric vector of length = 1 or length = 2. Standard deviation of
#'        random displacement (along each coordinate axis) of a point from
#'        its cluster centre. Therefore sigma correlates with cluster radius.
#'        When length(sigma) == 1 all species have the same cluster radius.
#'        When length(sigma) == 2 a linear relationship between species-specific sigma values
#'        and log(relative abundance) is simulated. Thereby sigma[1] is the cluster parameter
#'        for the least abundant species and sigma[2] the cluster parameter for the most
#'        abundant species.
#'        When length(sigma) == 1 and sigma is more than twice as large as the largest
#'        plot dimension, than a random Poisson distribution is simulated, which is more
#'        efficient than a Thomas cluster process.
#'
#' @param xmax scalar numeric - size of the simulated area in x-direction
#' @param ymax scalar numeric - size of the simulated area in y-direction

#' @param points.cluster scalar integer - numbers of individuals per cluster.
#'        If it is defined as sqrt(species abundance), i.e. there will be on average the number
#'        of clusters will equal the number of individuals per cluster.
#'
#' @details To generate a Thomas cluster process of a single species this function uses a C++
#'          re-implementation if the function \code{rThomas()} in the package \code{spatstat}.
#'
#' @return A dataframe with three columns
#' \enumerate{
#'    \item X : x-coordinate of individual
#'    \item Y : y-coordinate of individual
#'    \item SpecID : factor with species identity label (from 1 to number of species)
#' }
#'
#' @author Felix May
#'
#' @seealso \code{spatstat::rThomas()}
#'
Sim.Thomas.Community <- function(S, N,
                                 cv.abund = 1,
                                 sigma = 0.02,
                                 xmax = 1,
                                 ymax = 1,
                                 points.cluster = NULL,
                                 fix.S = F)
{
   #require(spatstat)

   max.dim <- ifelse(xmax >= ymax, xmax ,ymax)

   sim1 <- SAD.lognorm(S.pool = S, N.local = N, cv.abund = cv.abund, fix.S.local = fix.S)
   abund.vec <- sim1$abund
   cum.abund <- cumsum(abund.vec)

   S.local <- length(abund.vec)

   if (length(sigma) == 2){
      # linear relationship between sigma and log(relabund)
      # sigma = a1 + b1 * log(relabund)
      log.relabund <- log(abund.vec/sum(abund.vec))
      range.abund <- max(log.relabund) - min(log.relabund)

      if (range.abund != 0) b1 <- (sigma[2]-sigma[1])/range.abund
      else b1 <- 0
      a1 <- sigma[2] - b1*max(log.relabund)
      sigma.vec <- a1 + b1*log.relabund

      #    plot(sigma.vec~log.relabund)
      #    abline(a1,b1,col="red")
      #    abline(h=c(sigma[1],sigma[2]),col=1:2)
      #    abline(v=c(min(log.relabund),max(log.relabund)))
   }
   else {

      if (sigma > 2*max.dim){
         x <- runif(N,0,xmax)
         y <- runif(N,0,ymax)
         id.spec <- rep.int(1:S,times=abund.vec)

         dat1 <- data.frame(X=x, Y=y, SpecID=id.spec)
         return(dat1)
      }

      sigma.vec <- rep(sigma[1], times = S)
   }

   sim.dat <- data.frame(X = numeric(N),
                         Y = numeric(N))

   # create map for first species
   if (!is.numeric(points.cluster)){
      n.mother.points <- points.per.cluster <- sqrt(abund.vec) # assumption : similar numbers of cluster and of
      #individuals per cluster
   }
   else {
      points.per.cluster <- rep(points.cluster, S.local)
      n.mother.points <- abund.vec/points.per.cluster
      n.mother.points <- ifelse(n.mother.points < 0.1, 0.1, n.mother.points)
   }

   dat1 <- rThomas_rcpp(abund.vec[1], sigma=sigma.vec[1], mu=points.per.cluster[1])

   irange <- 1:cum.abund[1]
   sim.dat$X[irange] <- dat1$x
   sim.dat$Y[irange] <- dat1$y

   for (ispec in 2:S.local){
      dat1 <- rThomas_rcpp(abund.vec[ispec], sigma = sigma.vec[ispec], mu = points.per.cluster[ispec])

      irange <- (cum.abund[ispec-1] + 1):cum.abund[ispec]
      sim.dat$X[irange] <- dat1$x
      sim.dat$Y[irange] <- dat1$y
   }

   sim.dat$SpecID <- factor(rep(names(abund.vec), abund.vec))
   names(sim.dat) <- c("X","Y","SpecID")
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



