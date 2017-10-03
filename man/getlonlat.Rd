\name{getlonlat}
\alias{getlonlat}
\title{Getting Grid Information from a Raster Object}
\description{
Returns the grid specification of a \code{raster} object.
}
\usage{
## Basic usage
getlonlat( x )
}
\arguments{
   \item{x}{\code{RasterLayer} or \code{RasterStack} object.}
}
\value{
Returns a \code{list} object containing the following information:
   - \code{lons}/\code{lats}: vector of unique longitude and latitude points
   - \code{dx}/\code{dy}: coordinate increments in x/y direction (grid spacing)
   - \code{nx}/\code{ny}: number of grid points in x/y direction
}
\examples{
library("raster")
# Create an empty demo RasterLayer object
test <- raster(xmn=9.5,xmx=20.5,ymn=19.5,ymx=30.5,ncols=11,nrows=11)
# Getting grid information
getlonlat( test )
}
\author{Reto Stauffer}
