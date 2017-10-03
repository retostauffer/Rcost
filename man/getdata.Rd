\name{getdata}
\alias{getdata}
\title{Extracting Data from NetCDF File Connection}
\description{
Built to extract certain variables/forecast times from NetCDF files.
}
\usage{
## Basic usage
getdata(nc, init, varname, steps, level = NULL, subset = NULL, silent = FALSE )
}
\arguments{
   \item{nc}{Open \code{ncdf4} connection}
   \item{init}{\code{Date}, \code{POSIXt}, or \code{character} object which can
      be converted into \code{POSIXt}. Date/time when the forecast was initialized.
      Crucial to compute the forecast steps (lead times).}
   \item{varname}{\code{character}, name of the variable to load.}
   \item{steps}{Optional. If not set all forecast steps from the \code{nc} file
      will be returned. Single or multiple steps can be defined (\code{numeric}).}
   \item{level}{Used when loading variables from pressure level files where the
      level has to be specified as a \code{numeric} value.}
   \item{subset}{Optional, default \code{NULL}. If not set the whole grid from the
      \code{nc} file will be returned. Can also be a \code{\link{extent}} 
      object or a \code{\link[=raster]{RasterLayer}} object. See details.}
   \item{silent}{\code{logical}, default is \code{FALSE}. If \code{TRUE} some
      messages will be suppressed.}
}
\details{
\code{subset} allows to only load parts of the data from the \code{nc} connection.
If \code{subset} is of class \code{Extent} (see \code{\link{extent}}) then
only (a spatial) subset of the data will be taken. If \code{subset} is of class \code{RasterLayer}
(see \code{\link{raster}}) the extent of the \code{RasterLayer} will be used to crop
the data, the resolution (or grid specification) of this \code{RasterLayer} will be used to
re-sample the data. This allows to change the resolution when loading the data.
}
\value{
Returns a \code{RasterLayer} (if only one \code{step}) is loaded or a \code{RasterStack}
object where each layer corresponds to a forecast step. Names of the layers contain the
initial initial date plus the forecast step.
}
\examples{
\dontrun{
library('maps'); library("raster")

# Open netcdf connections
init <- as.POSIXct("2011-11-18 00:00")
init_str <- strftime(init,"%Y%m%d%H%M")
ncs <- nc_open(sprintf("costgrib/ECMWF_cost_sfc_%s.nc",init_str))
ncp <- nc_open(sprintf("costgrib/ECMWF_cost_pl_%s.nc", init_str))

# Loading "tcw" from surface file, step 24, full field
x <- getdata(ncs,init,"tcw",24)
plot(x); map(add=T)

# Using an extent subset and loading three different steps
x <- getdata(ncs,init,"tcw",c(12,24,48),subset=extent(c(10,30,40,50)))
plot(x)

# Loading pressure level data (u on 700 hPa) and resample the data
# onto a coarser resolution.
template <- raster(xmn=9.5,xmx=30.5,ymn=39.5,ymx=50.5,nrows=21,ncols=11,
               crs=crs("+proj=longlat +datum=WGS84 +ellps=WGS84 +no_defs +towgs84=0,0,0"))
x <- getdata(ncp,init,"u",144,level=700,subset=template)
plot(x); map(add=T)

# Close connections
nc_close(ncs)
nc_close(ncp)
}
}
\author{Reto Stauffer}
