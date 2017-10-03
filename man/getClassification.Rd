\name{getClassification}
\alias{getClassification}
\title{Getting the cost733 WLK Classification For a Given Date}
\description{
This function was the reason to set up this small R package.  This package
performs the cost733 WLK classification on gridde data sets such as the ECMWF
INTERIM or ECMWF HIRES.  The classification is performed based on the fortran
routine by Andreas Philipp, UNI Augsburg.

http://cost733class.geo.uni-augsburg.de/moin/cost733wiki/data/pages_bak/WLKC733/attachments/cost733_wlk_1.07.f90

The classification requires the variables ``\code{tcw}'' (total column water),
``\code{u700}'' and ``\code{v700}'' (zonal and meridional wind component in 700
hectopascal), as far as ``\code{z500}'' and ``\code{z925}''. These variables have
to be provided by the \code{nc} files given as input argument. See details on how
the netcdf files should look like.

The defaults of this method use the specificationas as proposed by Thomas Krennert
(ZAMG Vienna). The method here is not 100 percent generic as it e.g, expects a specific
file format (and file names) was I have only planned to use the script on 
ECMWF 00UTC deterministic forecasts at the moment.
}
\usage{
## Basic usage
getClassification(nc, date, steps, subset, pwclim, weights,
                  nsector=8, maindirthreshold=0.4,
                  zamg=FALSE, verbose=0)
}
\arguments{
   \item{nc}{\code{list} containing one or more open netcdf connection opened
      by \code{ncdf4::nc_open}.}
   \item{date}{\code{Date}, \code{POSIXt}, or \code{character} object which can
      be converted into \code{POSIXt}. Date/time when the forecast was initialized.
      Used to find the NetCDF files with a specific (fixed) name at the moment.}
   \item{steps}{Optional, if not set the classification will be performed on all
      steps provided by the NetCDF files.}
   \item{subset}{Optional, default \code{NULL}. See \code{\link{getdata}} for more
      details.} 
   \item{pwclim}{Optional. If not set the dry/wet classification cannot be performed.
      If set it has to be a \code{numeric} vector of length \code{366} containing
      the climatological ``precip. water content'' for each day of the year.}
   \item{weights}{Optional, if not set the ZAMG config with \code{c(15,2,1)}
      will be used for the weighting. If set it has to be a \code{numeric} vector
      of length 3 (for core/mid/margin weight).}
   \item{nsector}{Integer, default \code{8}. Number of wind sectors for the classification.}
   \item{maindirthreshold}{\code{numeric}, default is \code{0.4}. Used for wind sector
      classification. If less than \code{maindirthreshold} can be classified for one main
      wind direction the wind sector classification will be ``variable'' (\code{sector=0}).}
   \item{zamg}{Logical flag, default \code{FALSE}. If set to \code{TRUE} the main
      wind direction will be computed based on the classification by
      ``Trendanalyse Klien Endbericht'' (as suggested by Thomas Kennert, ZAMG).}
   \item{verbose}{Integer. Logical \code{TRUE}/\code{FALSE} will be converted to \code{0}/\code{1}.
      However, verbose levels up to \code{4} (highly verbose) can be set.}
}
\details{
For \code{subset} please see details in \code{\link{getdata}}.

The netcdf files have to provide the variables u and v (wind components) on 700 hectopascal,
z (geopotential height) on 500 and 925 hectopascal, and tcw (total column water).
To be able to extract the necessary information the netcdf file(s) have to provide
the following information:
   - dimension longitude (named longitude)
   - dimension latitude (named latitude)
   - dimension time (named time; units something like ``hours since 1900-01-01 00:00''
     or similar)
   - dimension level (named level) which is used when searching for level variables
I typically use two different files, a surface level netcdf (dimensions: longitude,
latitude, and time; containing at least tcw), and a pressure level file (dimensions:
longitude, latitude, time, and level; containing at least u700, v700, z500, z925) and
use both as input arguments. E.g.,
   - \code{getClassification( list(nc_open('<sfcfile>'),nc_open('<plfile>')), ... )}
}
\value{
Returns a \code{zoo} object with the following items:
   - \code{date}: date/time when the classification is valid
   - \code{yday}: day of the year
   - \code{sector}: classifier: wind sector
   - \code{cyclonic925}/\code{cyclonic500}: classifier: wheter it is
      \code{0=anticyclonic} or \code{1=cyclonic} in 925/500 hPa
   - \code{wet}: classifier: \code{-99=missing} (due to lack of pwclim), 
      \code{0=dry}, and \code{1=wet}.
   - \code{c925}/\code{c500}: curvature values used to create the classification.
   - \code{pwi}: computed pwi values which will be used for classification
}
\references{
Philipp A., C. Beck, R. Huth and J. Jacobeit (2014): Development and comparison
of circulation type classifications using the COST 733 dataset and software.
International Journal of Climatology. DOI: 10.1002/joc.3920
}
\examples{
# Open netcdf connections
require("ncdf4")
files <- demofiles()
nc    <- list(nc_open(files$sfc), nc_open(files$pl)) # Open NetCDF Files

# Extract initial date (hours since 1900-01-01 00:00:0.0)
# Take first date only.
init  <- as.POSIXct( ncvar_get(nc[[1]],"time") * 3600, origin = "1900-01-01" )[1]

# Perform classification
x1 <- getClassification( nc, init)
plot(x1, main = "COST733 Classification Output (x1)")

# Use different areal weights and only 4 wind sectors
x2 <- getClassification( nc, init, weights = c(10,5,2), nsector = 4 )
plot(x2, main = "COST733 Classification Output (x2)")

# Provide pw climatology (in this case simply a constant of
# 24.844 for all 1:365 days of a year): returns additional
# binary classification "wet" (simply pwi > pwiclim).
x3 <- getClassification( nc, init, weights = c(10,5,2), nsector = 8,
                         pwclim = rep(24.844,365) )
plot(x3, main = "COST733 Classification Output (x3)")

# Close NetCDF Files
for ( n in nc ) nc_close(n)
}
\author{Reto Stauffer}
