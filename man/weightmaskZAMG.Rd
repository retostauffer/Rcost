\name{weightmaskZAMG}
\alias{weightmaskZAMG}
\title{Returns cost733-ZAMG Style Weighting Mask}
\description{
The cost733 WLK weather type classification algorithm uses a weight mask.
This allows to ``focus'' on specific parts of the area of interest. There
are three different zones called ``core'', ``mid'' and ``margin''. 

\code{\link{getClassification}} provides an option which is called
\code{zamg=TRUE}. In case this is set the algorithm used the ZAMG
weighting mask as shown in ``Trendanalyse von hydro-meteorologischen Extremwerten''.
The are using a fixed weight-mask specification for Austria.

This function returns the weights with respect to the WLK classification
as used by the ZAMG in Vienna.
}
\usage{
## Basic usage
weightmaskZAMG( lon, lat, weights=c(15,2,1), ... )
}
\arguments{
   \item{lon}{\code{numeric} vector of longitude values.}
   \item{lat}{\code{numeric} vector of latitude values.}
   \item{weights}{\code{numeric} vector of length three with the weights
      for the three zones ``core'', ``mid'' and ``margin''.}
   \item{...}{Unused additional arguments.}
}
\value{
Returns a matrix of size \code{length(lon)} times \code{length(lat)} with
the corresponding weights.
}
\seealso{
   \code{\link{weightmask}},
   \code{\link{getClassification}}
}
\examples{
lons <- seq( 5,30,by=0.25)
lats <- seq(40,60,by=0.25)
w1 <- weightmaskZAMG( lons, lats )
w2 <- weightmaskZAMG( lons, lats, weights=c(3,2,1) )
w3 <- weightmaskZAMG( lons, lats, stripsize = 10 )
w4 <- weightmaskZAMG( lons, lats, weights=c(3,2,1), stripsize = 4.5 )

# Pretty ugly plots
par(mfrow=c(2,2))
   image(w1,zlim=c(1,15),col=1:15,main="weight mask c(15,2,1) stripsize 5")
   image(w2,zlim=c(1,15),col=1:15,main="weight mask c(3,2,1) stripsize 10")
   image(w3,zlim=c(1,15),col=1:15,main="weight mask c(15,2,1) stripsize 5")
   image(w4,zlim=c(1,15),col=1:15,main="weight mask c(3,2,1) stripsize 4.5")
}
\author{Reto Stauffer}
