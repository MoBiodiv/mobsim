% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Sample_quadrats.R
\name{sampling_grids}
\alias{sampling_grids}
\title{Creates square quadrats aligned on a regular grid}
\usage{
sampling_grids(
  n_quadrats,
  xmin,
  xmax,
  ymin,
  ymax,
  x0,
  y0,
  delta_x,
  delta_y,
  quadrat_size
)
}
\arguments{
\item{n_quadrats}{(integer) Number of sampling quadrats}

\item{xmin}{(numeric) minimum possible value on the x axis a quadrat can cover.}

\item{xmax}{(numeric) maximum possible value on the x axis a quadrat can cover.}

\item{ymin}{(numeric) minimum possible value on the y axis a quadrat can cover.}

\item{ymax}{(numeric) maximum possible value on the y axis a quadrat can cover.}

\item{x0, y0}{(numeric value) Lower left corner of the first quadrat in transect and grid sampling}

\item{delta_x}{(numeric value) Distance between consecutive quadrats in transect and grid sampling
in x-direction (the distance between the left sides is measured)}

\item{delta_y}{(numeric value) Distance between consecutive quadrats in transect and grid sampling
in y-direction (the distance between the lower sides is measured)}

\item{quadrat_size}{(numeric) width of the quadrats.}
}
\value{
a data.frame with 2 columns x and y giving  the coordinates of the
lower left corner of the square quadrats.
}
\description{
Creates square quadrats aligned on a regular grid
}
