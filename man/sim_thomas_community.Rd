% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Sim_Community.R
\name{sim_thomas_community}
\alias{sim_thomas_community}
\title{Simulate community with clumped spatial positions.}
\usage{
sim_thomas_community(
  s_pool,
  n_sim,
  sad_type = "lnorm",
  sad_coef = list(cv_abund = 1),
  fix_s_sim = FALSE,
  sigma = 0.02,
  cluster_points = NA,
  mother_points = NA,
  xmother = NA,
  ymother = NA,
  xrange = c(0, 1),
  yrange = c(0, 1),
  seed = NULL
)
}
\arguments{
\item{s_pool}{Number of species in the pool (integer)}

\item{n_sim}{Number of individuals in the simulated community (integer)}

\item{sad_type}{Root name of the species abundance distribution model of the
  species pool (character) - e.g., "lnorm" for the lognormal distribution
  (\code{\link[stats]{rlnorm}}); "geom" for the geometric distribution
  (\code{\link[stats]{rgeom}}), or "ls" for Fisher's log-series distribution
  (\code{\link[sads]{rls}}).

  See the table in \strong{Details} below, or \code{\link[sads]{rsad}}
  for all SAD model options.}

\item{sad_coef}{List with named arguments to be passed to the distribution
  function defined by the argument \code{sad_type}. An overview of parameter
  names is given in the table below.

  In \code{mobsim} the log-normal and the Poisson log-normal distributions
  can alternatively be parameterized by the coefficient of variation (cv)
  of the relative abundances in the species pool. Accordingly, \code{cv_abund}
  is the standard deviation of abundances divided by the mean abundance
  (no. of individuals / no. of species). \code{cv_abund} is thus negatively
  correlated with the evenness of the species abundance distribution.

  Please note that the parameters \emph{mu} and \emph{sigma} are not equal
  to the mean and standard deviation of the log-normal distribution.}

\item{fix_s_sim}{Should the simulation constrain the number of
species in the simulated local community? (logical)}

\item{sigma}{Mean displacement (along each coordinate axes) of a point from
its mother point (= cluster centre). \code{Sigma} correlates with cluster
extent. When \code{length(sigma) == length(abund_vec)}, each species
receives a specific cluster extent. Otherwise, the first value of \code{sigma}
is recycled and all species share the same cluster extent.
When \code{sigma} of any species is more than twice as large as the largest
plot dimension, a random Poisson distribution is simulated, which is more
efficient than a Thomas cluster process. The parameter \code{sigma} corresponds
to the \code{scale} parameter of the function \code{\link[spatstat.random]{rThomas}} in the package
\href{https://CRAN.R-project.org/package=spatstat.random}{spatstat.random}.}

\item{cluster_points}{Mean number of points per cluster. If this is
a single value, species have the same average number of points per cluster.
If this is a vector of the same length as \code{abund_vec}, each species has
a specific mean number of points per cluster.  If no value is provided, the
number of points per cluster is determined from the abundance and from
\code{mother_points}. If \code{mother_points} and \code{cluster_points} are given OR
\code{xmother} and \code{ymother}, and cluster points are given, \code{cluster_points} is overridden.
If \code{mother_points}=0, there will be no clustering even if \code{cluster_points}=400 (high clustering) because
\code{cluster_points} is overridden.
The parameter \code{cluster_points} corresponds to the
\code{mu} parameter of \code{spatstat.random::rThomas}.}

\item{mother_points}{Number of mother points (= cluster centres).
If this is a single value, all species have the same number of clusters.
For example \code{mother_points = 1} can be used to simulate only one cluster
per species, which then represents the complete species range.
If \code{mother_points} is a vector of the same length as \code{abund_vec},
each species has a specific number of clusters. If \code{mother_points} equals 0
there is no clustering and the distribution is homogeneous. If no value is provided,
the number of clusters is determined from the abundance and the number of points
per cluster (\code{cluster_points}).}

\item{xmother}{List of length equal to the number of species. Each list element
is a vector of x coordinates for every mother points. If one element is NA, the
the corresponding species is not clustered.}

\item{ymother}{List of length equal to the number of species. Each list element
is a vector of y coordinates for every mother points. If one element is NA, the
the corresponding species is not clustered.}

\item{xrange}{Extent of the community in x-direction. If this a numeric vector
of length 2, all species share the same range. To specify different x ranges for
all species, \code{xrange} should be a data.frame with 2 columns, min and max.}

\item{yrange}{Extent of the community in y-direction. If this a numeric vector
of length 2, all species share the same range. To specify different y ranges for
all species, \code{xrange} should be a data.frame with 2 columns, min and max.}

\item{seed}{Integer. Any integer passed to \code{set.seed} for reproducibility.}
}
\value{
A community object as defined by \code{\link{community}}
}
\description{
This function simulates a community with a certain abundance distribution and
with intraspecific aggregation, i.e. individuals of the same species are
distributed in clusters.
}
\details{
This function consecutively calls \code{\link{sim_sad}} and
\code{\link{sim_thomas_coords}}

See the documentations of \code{\link{sim_sad}} and
 \code{\link{sim_thomas_coords}} for details.
}
\examples{
com1 <- sim_thomas_community(s_pool = 20, n_sim = 500, sad_type = "lnorm",
                             sad_coef = list("meanlog" = 2, "sdlog" = 1),
                             sigma = 0.01)
plot(com1)

}
\author{
Felix May
}
