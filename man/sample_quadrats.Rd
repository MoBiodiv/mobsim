% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Sample_quadrats.R
\name{sample_quadrats}
\alias{sample_quadrats}
\title{Plot-based samples from a spatially-explicit census}
\usage{
sample_quadrats(
  comm,
  n_quadrats = 20,
  quadrat_area = 0.01,
  plot = TRUE,
  method = "random",
  avoid_overlap = TRUE,
  x0 = 0,
  y0 = 0,
  delta_x = 0.1,
  delta_y = 0.1,
  seed = NULL
)
}
\arguments{
\item{comm}{Community object from which the samples are generated}

\item{n_quadrats}{(integer) Number of sampling quadrats}

\item{quadrat_area}{(numeric) Area of the sampling quadrats}

\item{plot}{(logical) Should the sampling design be plotted? default to TRUE.}

\item{method}{(character) Available methods are \code{"random", "transect", "grid"}}

\item{avoid_overlap}{(logical) For the random sampling try to generate a design
without overlap of quadrats . Default is TRUE.}

\item{x0, y0}{(numeric value) Lower left corner of the first quadrat in transect and grid sampling}

\item{delta_x}{(numeric value) Distance between consecutive quadrats in transect and grid sampling
in x-direction (the distance between the left sides is measured)}

\item{delta_y}{(numeric value) Distance between consecutive quadrats in transect and grid sampling
in y-direction (the distance between the lower sides is measured)}

\item{seed}{(integer) Any integer passed to \code{set.seed} for reproducibility.}
}
\value{
A list with two items, \code{spec_dat} and \code{xy_dat}.
\code{spec_dat} is a data.frame with sampling quadrats in rows and species abundances
in columns, and \code{xy_dat} is a data.frame with sampling quadrats in rows
and the xy-coordinates of the quadrats (lower left corner) in columns.
}
\description{
This function allows to sample quadratic subplots from a spatially-explicit
community. The output format are a sites x species abundance table and a
sites x xy-coordinates table. The sites x species abundance is
a classical data format used in community ecology. The table generated
can be for instance be further analysed with the package \code{\link[vegan]{vegan}}.
}
\examples{
library(vegan)
sim_com1 <- sim_poisson_community(100, 10000)
comm_mat1 <- sample_quadrats(sim_com1, n_quadrats = 100,
quadrat_area = 0.002, method = "grid")
specnumber(comm_mat1$spec_dat)
diversity(comm_mat1$spec_dat, index = "shannon")

}
