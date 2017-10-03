\name{datetimeinfo}
\alias{datetimeinfo}
\title{Getting Time/Forecast Time Information from a NetCDF Connection}
\description{
Returns the time or forecast time information from a \code{ncdf4} connection.
Developed for some ECMWF forecasts, time information is stored as hours since
\code{1900-01-01} wherefore \code{init} is required to compute forecast steps.
}
\usage{
## Basic usage
datetimeinfo(nc,init)
}
\arguments{
   \item{nc}{Open \code{ncdf4} connection.}
   \item{init}{\code{Date}, \code{POSIXt}, or \code{character} which can
         be converted into a \code{POSIXt} object.}
}
\value{
Returns a \code{data.frame} with 3 columns. \code{init} contains
the \code{init} datetime, \code{valid} is the time for which the
forecast is valid, and \code{step} is the forecast step or lead time.
}
\examples{
\dontrun{
# Open NetCDF connection
nc <- nc_open("costgrib/ECMWF_cost_sfc_201111180000.nc")
# Get time information
datetimeinfo(nc,"2011-11-18")
# Close NetCDF connection
nc_close(nc)
}
}
\author{Reto Stauffer}
