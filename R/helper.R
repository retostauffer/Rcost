# -------------------------------------------------------------------
# Simply shows available dates/times/steps of the nc file
# -------------------------------------------------------------------
datetimeinfo <- function(nc,init) {
   # Extract steps
   time      <- list('str'=tolower(nc$dim$time$units))
   time$unit <- regmatches(time$str,regexpr("^[a-zA-Z]{1,}",time$str))
   time$orig <- regmatches(time$str,regexpr("[0-9]{4}-[0-9]{2}-[0-9]{2}\ [0-9:.]{1,}",time$str))
   if ( time$unit == "hours" ) {
      time$factor <- 3600
   } else if ( time$unit == "minutes" ) {
      time$factor <- 60
   } else if ( time$unit == "days" ) {
      time$factor <- 86000
   } else if ( time$unit == "seconds" ) {
      time$factor <- 1
   } else {
      stop(sprintf("Sorry, no rule for time unit \"%s\" in datetimeinfo function.",time$unit))
   }
   valid   <- as.POSIXct(ncvar_get(nc,'time')*time$factor,origin=time$orig)
   ncsteps <- as.numeric(valid - as.POSIXct(init),unit="secs") / 3600

   data.frame("init"=rep(init,length(ncsteps)),
              "valid"=valid,
              "step"=ncsteps)
}

# -------------------------------------------------------------------
# Helper function to get proper longitude and latitude vectors
# from a raster::RasterLayer object.
# -------------------------------------------------------------------
# We need longitude and latitude vectors of the layers
getlonlat <- function(x) {
   coord <- data.frame(sp::coordinates(x))
   dx <- unique(round(diff(sort(unique(coord$x))),5))
   dy <- unique(round(diff(sort(unique(coord$y))),5))
   if ( length(dx) != 1 | length(dy) != 1 )
      stop("Error in getlonlat: grid-spacing not regular!")
   nx   <- x@ncols
   ny   <- x@nrows
   #lons <- seq(x@extent@xmin+dx/2,x@extent@xmax-dx/2,length=nx)
   #lats <- seq(x@extent@ymin+dy/2,x@extent@ymax-dy/2,length=ny)
   lons <- unique(coord$x)
   lats <- unique(coord$y)
   return( list(lons=lons,lats=lats,dx=dx,dy=dy,nx=nx,ny=ny) )
}
