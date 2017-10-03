\name{weightmask}
\alias{weightmask}
\title{Returns cost733 Style Weighting Mask}
\description{
The cost733 WLK weather type classification algorithm uses a weight mask.
This allows to ``focus'' on specific parts of the area of interest. There
are three different zones called ``core'', ``mid'' and ``margin''. 

By default (cost733) uses a \code{stripsize} of \code{5} (cut's the area
in five rows and five columns). This creates a rectangular core of one
fifth times one fifth of the full area and two surrounding ``bands''.
The inner ``band'' is the ``mid'' zone, the outer one (reaching the
edges of the zone) the ``margin''.

Each of these three zones gets a weight, by default \code{15} for
the ``core'', \code{2} for ``mid'', and \code{1} for ``margin''.
During computation values within these zones are weighted according
to the \code{weights} set.
}
\usage{
## Basic usage
weightmask( lon, lat, weights=c(15,2,1), stripsize=5 )
}
\arguments{
   \item{lon}{\code{numeric} vector of longitude values.}
   \item{lat}{\code{numeric} vector of latitude values.}
   \item{weights}{\code{numeric} vector of length three with the weights
      for the three zones ``core'', ``mid'' and ``margin''.}
   \item{stripsize}{\code{numeric} value which is used to define the zones.}
}
\value{
Returns a matrix of size \code{length(lon)} times \code{length(lat)} with
the corresponding weights.
}
\seealso{
   \code{\link{weightmaskZAMG}},
   \code{\link{getClassification}}
}
\examples{
lons <- seq( 5,30,by=0.25)
lats <- seq(40,60,by=0.25)
w1 <- weightmask( lons, lats )
w2 <- weightmask( lons, lats, weights=c(3,2,1) )
w3 <- weightmask( lons, lats, stripsize = 10 )
w4 <- weightmask( lons, lats, weights=c(3,2,1), stripsize = 4.5 )

# Pretty ugly plots
par(mfrow=c(2,2))
   image(w1,zlim=c(1,15),col=1:15,main="weight mask c(15,2,1) stripsize 5")
   image(w2,zlim=c(1,15),col=1:15,main="weight mask c(3,2,1) stripsize 10")
   image(w3,zlim=c(1,15),col=1:15,main="weight mask c(15,2,1) stripsize 5")
   image(w4,zlim=c(1,15),col=1:15,main="weight mask c(3,2,1) stripsize 4.5")
}
\author{Reto Stauffer}
