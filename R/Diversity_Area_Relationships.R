
# -----------------------------------------------------------
#' Get diversity indices for one rectangle subplot in the community
#'
#' @param x0 x-coordinate of lower left corner
#' @param y0 x-coordinate of lower left corner
#' @param xsize size of the subplot in x-direction
#' @param ysize size of the subplot in x-direction
#' @param community Dataframe with three columns: x, y, species identity
#'
#' @return Named vector with four values
#' \enumerate{
#'    \item nSpecies - the number of species
#'    \item nEndemics - the number of endemicsw
#'    \item shannon - Shannon diversity index defined as \eqn{H = - \sum p_{i} * log(p_[i])},
#'    where \eqn{p_i} is the relative abundance of species i:
#'    \item simpson - Simpson diversity index (= probabiliy of interspecific encounter PIE)
#'    defined as \eqn{D =  1 - \sum p_i^2}
#' }

diversity.rect <- function(x0, y0, xsize, ysize, community)
{

   x <- community[,1]
   y <- community[,2]
   id.spec <- community[,3]

   # logical vector which trees are in the sampling rectangle
   in.rect <- (x >= x0 & x < (x0+xsize) & y >= y0 & y < (y0+ysize))

   spec.in <- unique(id.spec[in.rect])
   spec.out <- unique(id.spec[!in.rect])

   nSpecies <- length(spec.in)
   nEndemics <- length(spec.in[!spec.in %in% spec.out])

   abund <- table(id.spec[in.rect])
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

# # -----------------------------------------------------------
# # get the abundance distribution in a rectangle with upper left corner in x0,y0
# # size of the rectangle: xsize*ysize
# abund.rect <- function(x0, y0, xsize, ysize, community)
# {
#    x <- community[,1]
#    y <- community[,2]
#    id.spec <- community[,3]
#
#    # logical vector which trees are in the sampling rectangle
#    in.rect <- (x >= x0 & x < (x0+xsize) & y >= y0 & y < (y0+ysize))
#
#    return(abund)
# }

#' Get mean and sd of diversity indices in several equally sized subplots
#' of a community
#'
#' @param prop.A Proportion of the total area samples
#' @param community Dataframe with three columns: x, y, species identity
#' @param nrect Number of randomly located subplots
#' @param xext Minimum and maximum x-coordinates of the community
#' @param yext Minimum and maximum y-coordinates of the community
#'
#' @return Vector with mean and standard deviation of the following diversity
#' indices: (1) Number of species (2) Number of endemics (3) Shannon-diversity
#' (4) Simpson diversity

diversity.rand.rect <- function(prop.A = 0.25, community, nrect = 100, xext = c(0,1), yext=c(0,1),
                                exclude_zeros = F)
{
   x <- community[,1]
   y <- community[,2]

   dx.plot <- xext[2] - xext[1]
   dy.plot <- yext[2] - yext[1]

   area <- dx.plot*dy.plot*prop.A
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

   xpos <- runif(nrect, min = xext[1], max = xext[2] - dx.rect)
   ypos <- runif(nrect, min = yext[1], max = yext[2] - dy.rect)

   div.plots <- mapply(diversity.rect, xpos, ypos,
                       MoreArgs=list(xsize = dx.rect, ysize = dy.rect, community=community))

   if (exclude_zeros == T)
      div.plots <- div.plots[, div.plots["nSpecies",] > 0]

   return(c(meanSpec    = mean(div.plots["nSpecies",]),
            sdSpec      = sd(div.plots["nSpecies",]),
            meanEnd     = mean(div.plots["nEndemics",]),
            sdEnd       = sd(div.plots["nEndemics",]),
            meanShannon = mean(div.plots["shannon",]),
            sdShannon   = sd(div.plots["shannon",]),
            meanSimpson = mean(div.plots["simpson",], na.rm = T),
            sdSimpson   = sd(div.plots["simpson",], na.rm = T))
          )
}



#' Estimate diversity-area relationships
#'
#' @param prop.A Proportion of area sampled
#' @param community Dataframe with three columns: x, y, species identity
#' @param nrect Number of randomly located subplots
#' @param xext Minimum and maximum x-coordinates of the community
#' @param yext Minimum and maximum y-coordinates of the community
#'
#' @return Vector with mean and standard deviation of the following diversity
#' indices: (1) Number of species (2) Number of endemics (3) Shannon-diversity
#' (4) Simpson diversity
DivAR <- function(community, prop.A = seq(0.1, 1, by=0.1), nsamples=100, xext=c(0,1), yext=c(0,1),
                  exclude_zeros = F)
{
   nscales <- length(prop.A)
   xrange <- xext[2] - xext[1]
   yrange <- yext[2] - yext[1]

   div.area <- sapply(prop.A,
                      diversity.rand.rect,
                      community = community,
                      nrect = nsamples,
                      xext = xext,
                      yext = yext,
                      exclude_zeros = exclude_zeros)
   div.dat <- as.data.frame(t(div.area))
   div.dat <- cbind(propArea = prop.A, div.dat)

   return(div.dat)
}

# -----------------------------------------------------------
# Distance decay
distance.decay <- function(community, prop.A=0.05, nsamples=30, xext=c(0,1), yext=c(0,1), method="jaccard")
{
   require(vegan)

   dx.plot <- xext[2] - xext[1]
   dy.plot <- yext[2] - yext[1]

   area <- dx.plot*dy.plot*prop.A
   square.size <- sqrt(area)

   xpos <- runif(nsamples,min=0, max=xext[2]-square.size)
   ypos <- runif(nsamples,min=0, max=yext[2]-square.size)

   d <- dist(cbind(xpos,ypos))

   com.tab <- mapply(abund.rect, xpos, ypos,
                     MoreArgs=list(xsize=square.size, ysize=square.size, community=community))

   jaccard <- 1 - vegdist(t(com.tab),method=method,binary=T)
   dat.out <- data.frame(distance = as.numeric(d), similarity = as.numeric(jaccard))
   return(dat.out)
}

