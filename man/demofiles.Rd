\name{demofiles}
\alias{demofiles}
\title{Returns Rcost Package Demo File Locations}
\description{
This package contains some demo data used for the examples and
vignettes. This function simply returns the location of these
files on the local disc.
}
\usage{
## Basic usage
demofiles( silent = FALSE )
}
\arguments{
   \item{silent}{\code{logical}. If \code{TRUE} console output is supressed.}
}
\value{
Returns a \code{list} object with two elements. \code{sfc} locates
the path to the demo surface data NetCDF file, \code{pl} the one to the
corresponding pressure level NetCDF file.

Both files contain five subsequent ERA ITERIM reanalysis fields from
July 1th, 2000 to July 5th, 2000.
}
\examples{
files <- demofiles( )
print(files)
}
\author{Reto Stauffer}
