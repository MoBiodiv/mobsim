#' Simulate species abundance distributions
#'
#' Simulate species abundance distribution (SAD) of a local community with
#' user-defined number of species and relative abundance distribution in the pool,
#' and user-defined number of individuals in the simulated local community.
#'
#' @param s_pool Number of species in the pool (integer)
#'
#' @param n_sim  Number of individuals in the simulated community (integer)
#'
#' @param sad_type Root name of the species abundance distribution model of the
#'   species pool (character) - e.g., "lnorm" for the lognormal distribution
#'   (\code{\link[stats]{rlnorm}}); "geom" for the geometric distribution
#'   (\code{\link[stats]{rgeom}}), or "ls" for Fisher's log-series distribution
#'   (\code{\link[sads]{rls}}).
#'
#'   See the table in \strong{Details} below, or \code{\link[sads]{rsad}}
#'   for all SAD model options.
#'
#' @param sad_coef List with named arguments to be passed to the distribution
#'   function defined by the argument \code{sad_type}. An overview of parameter
#'   names is given in the table below.
#'
#'   In \code{mobsim} the log-normal and the Poisson log-normal distributions
#'   can alternatively be parameterized by the coefficient of variation (cv)
#'   of the relative abundances in the species pool. Accordingly, \code{cv_abund}
#'   is the standard deviation of abundances divided by the mean abundance
#'   (no. of individuals / no. of species). \code{cv_abund} is thus negatively
#'   correlated with the evenness of the species abundance distribution.
#'
#'   Please note that the parameters \emph{mu} and \emph{sigma} are not equal
#'   to the mean and standard deviation of the log-normal distribution.
#'
#' @param fix_s_sim Should the simulation constrain the number of
#'   species in the simulated local community? (logical)
#'
#' @details The function \code{sim_sad} was built using code of the function
#'   \code{\link[sads]{rsad}} from the R package \code{\link{sads}}. However, in
#'   contrast to \code{\link[sads]{rsad}}, the function \code{sim_sad} allows to
#'   define the number of individuals in the simulated local community. This is
#'   implemented by converting the abundance distribution simulated based on
#'   \code{\link[sads]{rsad}} into a relative abundance distribution. This
#'   relative abundance distribution is considered as the species pool for the
#'   local community. In a second step the required no. of individuals \code{(n_sim)}
#'   is sampled (with replacement) from this relative abundance distribution.
#'
#'   Please note that this might effect the interpretation of the parameters of
#'   the underlying statistical distribution, e.g. the mean abundance will always
#'   be \code{n_sim/n_pool} irrespective of the settings of \code{sad_coef}.
#'
#'   When \code{fix_s_sim = FALSE} the species number in the local
#'   community might deviate from \code{s_pool} due to stochastic sampling. When
#'   \code{fix_s_sim = TRUE} the local number of species will equal
#'   \code{s_pool}, but this constraint can result in systematic biases from the
#'   theoretical distribution parameters. Generally, with \code{fix_s_sim = TRUE}
#'   additional very rare species will be added to the community, while the abundance
#'   of the most common ones is reduced to keep the defined number of individuals.
#'
#'   Here is an overview of all available models (\code{sad_type}) and their
#'   respective coefficients (\code{sad_coef}). Further information is provided
#'   by the documentation of the specific functions that can be accesses by the
#'   links. Please note that the coefficient \code{cv_abund} for the log-normal
#'   and Poisson log-normal model are only available within \code{mobsim}.
#'
#' \tabular{lllll}{
#'    \strong{SAD function} \tab \strong{Distribution name} \tab \strong{coef #1} \tab \strong{coef #2} \tab \strong{coef #3} \cr
#'    \code{\link[sads]{rbs}} \tab Mac-Arthur's brokenstick \tab N \tab S \tab \cr
#'    \code{\link[stats]{rgamma}} \tab Gamma distribution \tab shape \tab rate \tab scale \cr
#'    \code{\link[stats]{rgeom}} \tab Geometric distribution \tab prob \tab \tab \cr
#'    \code{\link[stats]{rlnorm}} \tab	Log-normal distributions \tab	meanlog \tab sdlog \tab cv_abund \cr
#'    \code{\link[sads]{rls}} \tab Fisher's log-series distribution \tab N \tab alpha \tab \cr
#'    \code{\link[sads]{rmzsm}} \tab Metacommunity zero-sum multinomial \tab J \tab theta \tab \cr
#'    \code{\link[stats]{rnbinom}} \tab Negative binomial distribution \tab size \tab	prob \tab mu \cr
#'    \code{\link[sads]{rpareto}} \tab Pareto distribution \tab shape \tab scale \tab \cr
#'    \code{\link[sads]{rpoilog}} \tab Poisson-lognormal distribution \tab	mu	 \tab sigma \tab cv_abund \cr
#'    \code{\link[sads]{rpower}} \tab Power discrete distributions \tab s \tab \tab \cr
#'    \code{\link[sads]{rpowbend}} \tab Puyeo's Power-bend discrete distribution \tab s \tab omega \tab \cr
#'    \code{\link[stats]{rweibull}} \tab Weibull distribution \tab shape \tab scale \tab \cr
#'}
#'
#'
#' @return Object of class \code{sad}, which contains a named integer vector
#'  with species abundances
#'
#' @author Felix May
#'
#' @examples
#' #Simulate log-normal species abundance distribution
#' sad_lnorm1 <- sim_sad(s_pool = 100, n_sim = 10000, sad_type = "lnorm",
#'                       sad_coef = list("meanlog" = 5, "sdlog" = 0.5))
#' plot(sad_lnorm1, method = "octave")
#' plot(sad_lnorm1, method = "rank")
#'
#' # Alternative parameterization of the log-normal distribution
#' sad_lnorm2 <- sim_sad(s_pool = 100, n_sim = 10000, sad_type = "lnorm",
#'                       sad_coef = list("cv_abund" = 0.5))
#' plot(sad_lnorm2, method = "octave")
#'
#' # Fix species richness in the simulation by adding rare species
#' sad_lnorm3a <- sim_sad(s_pool = 500, n_sim = 10000, sad_type = "lnorm",
#'                        sad_coef = list("cv_abund" = 5), fix_s_sim = TRUE)
#' sad_lnorm3b <- sim_sad(s_pool = 500, n_sim = 10000, sad_type = "lnorm",
#'                        sad_coef = list("cv_abund" = 5))
#'
#' plot(sad_lnorm3a, method = "rank")
#' points(1:length(sad_lnorm3b), sad_lnorm3b, type = "b", col = 2)
#' legend("topright", c("fix_s_sim = TRUE","fix_s_sim = FALSE"),
#'        col = 1:2, pch = 1)
#'
#' # Different important SAD models
#'
#' # Fisher's log-series
#' sad_logseries <- sim_sad(s_pool = NULL, n_sim = 10000, sad_type = "ls",
#'                          sad_coef = list("N" = 1e5, "alpha" = 20))
#'
#' # Poisson log-normal
#' sad_poilog <- sim_sad(s_pool = 100, n_sim = 10000, sad_type = "poilog",
#'                       sad_coef = list("mu" = 5, "sig" = 0.5))
#'
#' # Mac-Arthur's broken stick
#' sad_broken_stick <- sim_sad(s_pool = NULL, n_sim = 10000, sad_type = "bs",
#'                             sad_coef = list("N" = 1e5, "S" = 100))
#'
#' # Plot all SADs together as rank-abundance curves
#' plot(sad_logseries, method = "rank")
#' lines(1:length(sad_lnorm2), sad_lnorm2, type = "b", col = 2)
#' lines(1:length(sad_poilog), sad_poilog, type = "b", col = 3)
#' lines(1:length(sad_broken_stick), sad_broken_stick, type = "b", col = 4)
#' legend("topright", c("Log-series","Log-normal","Poisson log-normal","Broken stick"),
#'        col = 1:4, pch = 1)
#'
#' @export

sim_sad <- function(s_pool, n_sim,
                    sad_type = c("lnorm", "bs", "gamma", "geom", "ls",
                                 "mzsm","nbinom", "pareto", "poilog", "power",
                                 "powbend", "weibull"),
                    sad_coef = list("cv_abund" = 1),
                    fix_s_sim = FALSE)
{
   sad_type <- match.arg(sad_type)

   if (!is.numeric(n_sim) || n_sim <= 0)
      stop("n_sim has to be a positive integer number")

   n_sim <- round(n_sim, digits = 0)

   if (class(sad_coef) != "list" | is.null(names(sad_coef))) stop("coef must be a named list!")

   # Handles parameters that give the community size
   if (sad_type %in% c("bs", "ls", "mzsm")) {
      S <- switch(sad_type,
                  bs = sad_coef$S,
                  ls = sad_coef$alpha * log ( 1 + sad_coef$N / sad_coef$alpha ),
                  mzsm = sum(sad_coef$theta / (1:sad_coef$J) *
                             (1 - (1:sad_coef$J)/sad_coef$J)^(sad_coef$theta - 1))
                 )
      S <- round(S)
      if (!is.null(s_pool)){
         warning(paste("For the selected SAD model the value of s_pool is ignored.
  s_pool calculated from the SAD model coefficients is", S, "species."))
      }
      s_pool <- S
   } else {
      if (is.null(s_pool) || is.na(s_pool) || !is.numeric(s_pool) || s_pool <= 0)
         stop("The argument s_pool is mandatory for the selected sad and has to be a positive integer number.")
      s_pool <- round(s_pool, digits = 0)
   }

   if (s_pool > 1){

      #alternative parameterization for lnorm and poilog
      if ((sad_type == "lnorm" || sad_type == "poilog") &&
          names(sad_coef)[1] == "cv_abund"){
         mean_abund <- n_sim/s_pool
         sd_abund <-  mean_abund * sad_coef$cv_abund
         sigma1 <- sqrt(log(sd_abund^2/mean_abund^2 + 1))
         mu1 <- log(mean_abund) - sigma1^2/2

         # mean1 <- exp(mu1 + sigma1^2/2)
         # sd1 <- exp(mu1 + sigma1^2/2) * sqrt(exp(sigma1^2) - 1)
         # cv1 <- sd1/mean1

         if (sad_type == "lnorm")
            sad_coef <- list("meanlog" = mu1, "sdlog" = sigma1)
         if (sad_type == "poilog")
            sad_coef <- list("mu" = mu1, "sig" = sigma1)
      }

      # Generates the "community"
      if (sad_type %in% c("gamma","geom","lnorm","nbinom","weibull")){
         sadr <- utils::getFromNamespace(paste("r", sad_type, sep=""), ns = "stats")
      } else {
         sadr <- utils::getFromNamespace(paste("r", sad_type, sep=""), ns = "sads")
      }
      abund_pool <- do.call(sadr, c(list(n = s_pool), sad_coef))

      abund_pool <- abund_pool[abund_pool > 0]
      rel_abund_pool <- abund_pool/sum(abund_pool)

      sample_vec <- sample(x = length(rel_abund_pool),
                           size = n_sim, replace = TRUE,
                           prob = rel_abund_pool)

      abund2 <- as.numeric(sort(table(sample_vec), decreasing = TRUE))

      if (fix_s_sim == TRUE & length(abund2) < s_pool){
         s_diff <- s_pool - length(abund2)
         abund2 <- c(abund2, rep(1, s_diff))
         n <- sum(abund2)

         #randomly remove individuals until target level is reached
         while (n > n_sim){
            rel_abund <- abund2/sum(abund2)
            # draw proportional to relative abundance
            irand <- sample(1:s_pool, size = 1, prob = rel_abund)
            if (abund2[irand] > 1) abund2[irand] <- abund2[irand] - 1
            n <- sum(abund2)
         }
      }
   } else { # end if(s_pool > 1)
      abund2 <- n_sim
   }

   names(abund2) <- paste("species", 1:length(abund2), sep = "")
   class(abund2) <- c("sad","integer")
   return(abund2)
}

#' Print summary of species abundance distribution object
#'
#' @param object Community object of class \code{sad}
#'
#' @param ... Additional arguments passed to \code{\link{print}}.
#'
#' @seealso \code{\link{sim_sad}}
#'
#' @export
#'
summary.sad <- function(object, ...)
{
   cat("Species abundance distribution\n\n")
   cat("No. of individuals: ", sum(object), "\n")
   cat("No. of species: ", length(object), "\n\n")
   cat("Min. abundance: ", min(object), "\n")
   cat("Mean abundance: ", mean(object), "\n")
   cat("Max. abundance: ", max(object), "\n")
}

#' Plot species abundance distributions
#'
#' @param x Vector with species abundances (integer vector)
#'
#' @param ... Additional graphical parameters used in \code{\link[graphics]{plot}}
#' or \code{\link[graphics]{barplot}}
#'
#' @param method Plotting method, partial match to \code{"octave"} or \code{"rank"}
#'
#' @details With \code{method = "octave"} a histogram showing the number
#' species in several abundance classes is generated. The abundance class
#' are a simplified version of the "octaves" suggested by Preston (1948), which
#' are based on log2-binning. The first abundance class includes species
#' with 1 individual, the second with 2, the third with 3-4, the fourth with 5-8, etc.
#'
#' With \code{method = "rank"} rank-abundance curve is generated with
#' species abundance rank on the x-axis (descending) and species abundance on
#' the y-axis (Hubbell 2001).
#'
#' @references
#' Preston 1948. The Commonness, and rarity, of species. Ecology 29(3):254-283.
#'
#' Hubbell 2001. The unified neutral theory of biodiversity and biogeography.
#' Princeton University Press.
#'
#' @examples
#' abund1 <- sim_sad(s_pool = 100, n_sim = 10000, sad_type = "lnorm",
#'                   sad_coef = list("cv_abund" = 1))
#' plot(abund1, method = "octave")
#' plot(abund1, method = "rank")
#'
#' @export
#'
plot.sad <- function(x, ..., method = c("octave","rank"))
{
   method <- match.arg(method)

   if (method == "rank")
      graphics::plot(sort(as.numeric(x), decreasing = TRUE), type="b", log="y",
                     xlab="Species rank", ylab="Species abundance",
                     main = "Rank-abundance curve", las = 1, ...)

   if (method == "octave"){

      # code adopted from untb:preston()
      max_abund <- max(x)
      n <- 1 + ceiling(log(max_abund)/log(2)) # number of abundance classes

      if (n < 2) breaks <- c(0, 1)
      else       breaks <- c(0, 2^(0:(n - 2)), max_abund)

      r <- graphics::hist(x, plot = FALSE, breaks = breaks, right = TRUE)
                          abund_dist <- r$counts

      if (n <= 2) names(abund_dist) <- c("1","2")[1:n]
      else        names(abund_dist) <- c("1", "2",
                                         paste(breaks[-c(1:2, length(breaks))] + 1,
                                          "-", breaks[-c(1:3)], sep = ""))

      graphics::barplot(height = as.numeric(abund_dist),
                        names.arg = names(abund_dist),
                        xlab = "Abundance class", ylab ="No. of species",
                        main = "Preston octave plot", las = 1, ...)
   }
}

#' Create spatial community object
#'
#' Creates a spatial community object with defined extent and with coordinates
#' and species identities of all individuals in the community.
#'
#' @param x,y Coordinates of individuals (numeric)
#' @param spec_id Species names or IDs; can be integers, characters or factors
#' @param xrange Extent of the community in x-direction (numeric vector of length 2)
#' @param yrange Extent of the community in y-direction (numeric vector of length 2)
#'
#' @return Community object which includes three items:
#' \enumerate{
#'    \item census: data.frame with three columns: x, y, and species names for
#'     each individual
#'    \item x_min_max: extent of the community in x-direction
#'    \item y_min_max: extent of the community in y-direction
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
   if (length(xrange) < 2 | length(yrange) < 2) ("Error: missing ranges for x or y!")

   if (xrange[1] > min(x)) return("Error: Inappropriate ranges for x!")
   if (xrange[2] < max(x)) return("Error: Inappropriate ranges for x!")

   if (yrange[1] > min(y)) return("Error: Inappropriate ranges for y!")
   if (yrange[2] < max(y)) return("Error: Inappropriate ranges for y!")


   points <- data.frame(x = as.numeric(x), y = as.numeric(y),
                        species = as.factor(spec_id))

   comm <- list(census = points,
               x_min_max = as.numeric(xrange[1:2]),
               y_min_max = as.numeric(yrange[1:2])
               )

   class(comm) <- "community"

   return(comm)
}

#' Print summary of spatial community object
#'
#' @param object Community object of class \code{\link{community}}
#'
#' @param ... Additional arguments passed to \code{\link{print}}.
#'
#' @export
#'
summary.community <- function(object, ...)
{
   cat("No. of individuals: ", nrow(object$census), "\n")
   cat("No. of species: ", length(unique(object$census$species)), "\n")
   cat("x-extent: ", object$x_min_max, "\n")
   cat("y-extent: ", object$y_min_max, "\n\n")
   print(summary(object$census))
}

#' Plot spatial community object
#'
#' Plot positions and species identities of all individuals in a community object.
#'
#' @param x Community object
#' @param col Colour vector to mark species identities
#' @param pch Plotting character to mark species identities
#' @param ... Other parameters to \link[graphics]{plot}
#'
#' @examples
#' sim1 <- sim_thomas_community(30, 500)
#' plot(sim1)
#'
#' @export
#'
plot.community <- function(x, ..., col = NULL, pch = NULL)
{
   comm <- x

   nspec <- length(table(comm$census$species))
   if (is.null(col))  col <- grDevices::rainbow(nspec)
   if (is.null(pch))  pch <- 19

   graphics::plot(y ~ x, data = comm$census, xlim = comm$x_min_max,
                  ylim = comm$y_min_max, col = col[comm$census$species],
                  pch = pch, las = 1, asp = 1, ...)
}

#' Get species abundance distribution from community object
#'
#' @param comm Community object
#'
#' @return Object of class \code{sad}, which contains a named integer vector
#'  with species abundances
#'
#' @examples
#' sim1 <- sim_poisson_community(s_pool = 200, n_sim = 20000, sad_type = "lnorm",
#'                               sad_coef = list("cv_abund" = 2))
#' sad1 <- community_to_sad(sim1)
#' plot(sad1, method = "rank")
#' plot(sad1, method = "octave")
#'
#' @export
#'
community_to_sad <- function(comm)
{
   if (class(comm) != "community")
      stop("community_to_sad requires a community object as input. See ?community.")

   abund <- table(comm$census$species)
   class(abund) <- "sad"

   return(abund)
}


#' Simulate random spatial coordinates
#'
#' Add random spatial positions to a species abundance distribution.
#'
#' @param abund_vec Species abundance vector (integer)
#' @param xrange Extent of the community in x-direction (numeric vector of length 2)
#' @param yrange Extent of the community in y-direction (numeric vector of length 2)
#'
#' @return A community object as defined by \code{\link{community}}.
#'
#' @author Felix May
#'
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
   x <- stats::runif(n, xrange[1], xrange[2])
   y <- stats::runif(n, yrange[1], yrange[2])

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
#' @param xrange Extent of the community in x-direction (numeric vector of length 2)
#' @param yrange Extent of the community in y-direction (numeric vector of length 2)
#'
#' @return A community object as defined by \code{\link{community}}.

#' @author Felix May
#'
#' @examples
#' com1 <- sim_poisson_community(s_pool = 20, n_sim = 500, sad_type = "lnorm",
#' sad_coef = list("meanlog" = 2, "sdlog" = 1))
#' plot(com1)
#'
#' @export
#'
sim_poisson_community <- function(s_pool,
                                  n_sim,
                                  sad_type = "lnorm",
                                  sad_coef = list("cv_abund" = 1),
                                  fix_s_sim = FALSE,
                                  xrange= c(0,1),
                                  yrange = c(0,1)
                                  )
{
   sim1 <- sim_sad(s_pool = s_pool, n_sim = n_sim,
                   sad_type = sad_type,
                   sad_coef = sad_coef,
                   fix_s_sim = fix_s_sim)
   abund_vec <- sim1

   sim_dat <- sim_poisson_coords(abund_vec = abund_vec,
                                 xrange = xrange, yrange = yrange)
   return(sim_dat)
}


#' Simulate clumped spatial coordinates
#'
#' Add clumped (aggregated) positions to a species abundance distribution.
#' Clumping is simulated using a Thomas cluster process, also known as Poisson
#' cluster process (Morlon et al. 2008, Wiegand & Moloney 2014)
#'
#' @param abund_vec Species abundance vector (integer)
#'
#' @param sigma Mean displacement (along each coordinate axes) of a point from
#' its mother point (= cluster centre). \code{Sigma} correlates with cluster
#' extent. When \code{length(sigma) == length(abund_vec)}, each species
#' receives a specific cluster extent. Otherwise, the first value of \code{sigma}
#' is recycled and all species share the same cluster extent.
#' When \code{sigma} of any species is more than twice as large as the largest
#' plot dimension, a random Poisson distribution is simulated, which is more
#' efficient than a Thomas cluster process. The parameter \code{sigma} corresponds
#' to the \code{scale} parameter of the function \code{rThomas} in the package
#' \href{https://CRAN.R-project.org/package=spatstat}{spatstat}.
#'
#'
#' @param mother_points Number of mother points (= cluster centres).
#' If this is a single value, all species have the same number of clusters.
#' For example \code{mother_points = 1} can be used to simulate only one cluster
#' per species, which then represents the complete species range.
#' If \code{mother_points} is a vector of the same length as \code{abund_vec},
#' each species has a specific number of clusters. If no value is provided, the
#' number of clusters is determined from the abundance and the number of points
#' per cluster (\code{cluster_points}).
#'
#' @param cluster_points Mean number of points per cluster. If this is
#' a single value, species have the same average number of points per cluster.
#' If this is a vector of the same length as \code{abund_vec}, each species has
#' a specific mean number of points per cluster.  If no value is provided, the
#' number of points per cluster is determined from the abundance and from
#' \code{mother_points}.  The parameter \code{cluster_points} corresponds to the
#' \code{mu} parameter of \code{spatstat::rThomas}.
#'
#' @param xrange Extent of the community in x-direction (numeric vector of length 2)
#' @param yrange Extent of the community in y-direction (numeric vector of length 2)
#'
#' @details To generate a Thomas cluster process of a single species this
#' function uses a C++ re-implementation of the function
#' \code{rThomas} in the package
#' \href{https://CRAN.R-project.org/package=spatstat}{spatstat}.
#'
#' There is an inherent link between the parameters \code{abund_vec},
#' \code{mother_points}, and \code{cluster_points}. For every species the
#' abundance has to be equal to the number of clusters
#' (\code{mother_points}) times the number of points per cluster
#' (\code{cluster_points}).
#'
#' \deqn{abundance = mother_points * cluster_points}
#'
#' Accordingly, if one of the parameters is provided, the other one is directly
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
#' Morlon et al. 2008. A general framework for the distance-decay of similarity
#' in ecological communities. Ecology Letters 11, 904-917.
#'
#' Wiegand and Moloney 2014. Handbook of Spatial Point-Pattern Analysis in Ecology.
#' CRC Press
#'
#' @author Felix May
#'
#' @seealso \code{\link[spatstat]{rThomas}}
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
#' # With large sigma the distribution will be essentially random (see Details)
#' sim4 <- sim_thomas_coords(abund, sigma = 10)
#' plot(sim4)
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
   if (all(!is.na(mother_points))){

      if (length(mother_points) == s_local)
         n_mothers <- mother_points
      else
         n_mothers <- rep(mother_points[1], s_local)

      points_per_cluster <- abund_vec / n_mothers

   } else {

      if (all(!is.na(cluster_points))){

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
      x1 <- stats::runif(abund_vec[1], xrange[1], xrange[2])
      y1 <- stats::runif(abund_vec[1], yrange[1], yrange[2])
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
            x1 <- stats::runif(abund_vec[ispec], xrange[1], xrange[2])
            y1 <- stats::runif(abund_vec[ispec], yrange[1], yrange[2])
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
#' with intraspecific aggregation, i.e. individuals of the same species are
#' distributed in clusters.
#'
#' This function consecutively calls \code{\link{sim_sad}} and
#' \code{\link{sim_thomas_coords}}
#'
#' @inheritParams sim_sad
#'
#' @param sigma Mean displacement (along each coordinate axes) of a point from
#' its mother point (= cluster centre).
#'
#' @param mother_points Number of mother points (= cluster centres).
#'
#' @param cluster_points Mean number of points per cluster.
#'
#' @param xrange Extent of the community in x-direction (numeric vector of length 2)
#' @param yrange Extent of the community in y-direction (numeric vector of length 2)
#'
#' @details See the documentations of \code{\link{sim_sad}} and
#'  \code{\link{sim_thomas_coords}} for details.
#'
#' @return A community object as defined by \code{\link{community}}
#'
#' @author Felix May
#'
#' @examples
#' com1 <- sim_thomas_community(s_pool = 20, n_sim = 500, sad_type = "lnorm",
#'                              sad_coef = list("meanlog" = 2, "sdlog" = 1),
#'                              sigma = 0.01)
#' plot(com1)
#'
#' @export
#'
sim_thomas_community <- function(s_pool, n_sim,
                                 sad_type = "lnorm",
                                 sad_coef = list("cv_abund" = 1),
                                 fix_s_sim = FALSE,
                                 sigma = 0.02,
                                 cluster_points = NA,
                                 mother_points = NA,
                                 xrange = c(0,1),
                                 yrange = c(0,1)
                                 )
{
   sim1 <- sim_sad(s_pool = s_pool, n_sim = n_sim,
                   sad_type = sad_type,
                   sad_coef = sad_coef,
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




