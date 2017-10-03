# -------------------------------------------------------------------
# Helper function to extract the data we need from the
# netcdf file.
# param:  nc, list of netcdf connections returned by. Can be more
#         than one, variables will be searched in all of them.
# param:  varname. Character, name of the variable.
# param:  step. Integer, forecast stpe in hours.
# param:  level. Either NULL if 'level' is not required (e.g., for
#         surface variables) or integer indicating the level.
# param:  subset. Either NULL if the full field should be returned,
#         a raster::Extent object, or a raster::RasterLayer object.
# -------------------------------------------------------------------
getdata <- function(nc, init, varname, steps, level = NULL, subset = NULL, silent = FALSE ) {

   found <- NULL
   for ( i in seq_along(nc) ) {
      if ( varname %in% names(nc[[i]]$var) ) { found <- i; break; }
   }
   # If variable cannot be found
   if ( is.null(found) && silent ) {
      warning("Variable not found, return NULL") 
      return(NULL)
   } else if ( is.null(found) ) {
      stop(sprintf("Cannot find variable \"%s\" in the netcdf files.",varname))
   }

   # Else take corresponding netcdf file to retrieve the data set
   nc <- nc[[found]]; rm(found)

   # Extract steps using datetimeinfo helper function
   tmp     <- datetimeinfo(nc,init)
   valid   <- tmp$valid
   ncsteps <- tmp$step

   # Check variable name
   if ( ! varname %in% names(nc$var) ) stop(sprintf("Variable \"%s\" not in netcdf file",varname))
   # Time index
   # If there is only one time step (length of time dimension == 1)
   # we have to change the request. Therefore, set time.idx = FALSE if so.
   time.idx <- which(ncsteps %in% steps)
   if ( length(time.idx) == 1 & nc$dim$time$len == 1 ) time.idx = FALSE
   if ( length(time.idx) != length(steps) ) stop("Problems to find some of the requested steps!")
   # Checking subset range if set
   lons <- ncvar_get(nc,"longitude"); lats <- ncvar_get(nc,"latitude")

   # If the subset is a raster::Extent object: take areal subset
   # with the Extent as bounding box.
   # If subset is a raster::RasterLayer object take extent of the raster.
   if ( class(subset) %in% c("Extent","RasterLayer") ) {
      e <- extent(subset)
      if ( e@xmin < min(lons) | e@xmax > max(lats) |
           e@ymin < min(lats) | e@ymax > max(lats) ) {
         stop("Sorry, subset specification is outside limits.")
      }
      idx.lon <- which(lons >= e@xmin & lons <= e@xmax) 
      idx.lat <- which(lats >= e@ymin & lats <= e@ymax) 
      lons <- lons[idx.lon]; lats <- lats[idx.lat]
   }

   if ( ! "level" %in% names(nc$dim) ) {
      if ( identical(time.idx,FALSE) ) {
         if ( is.null(subset) ) { data <- ncvar_get(nc,varname)[,] }
         else                   { data <- ncvar_get(nc,varname)[idx.lon,idx.lat] }
      } else {
         if ( is.null(subset) ) { data <- ncvar_get(nc,varname)[,,time.idx] }
         else                   { data <- ncvar_get(nc,varname)[idx.lon,idx.lat,time.idx] }
      }
   } else {
      if ( is.null(level) ) stop("Sorry, input \"level\" required for this file")
      lev.idx <- which(ncvar_get(nc,"level")==level)
      if ( length(lev.idx) != 1 ) stop(sprintf("Problems to find level=\"%d\"",level))
      if ( identical(time.idx,FALSE) ) {
         if ( is.null(subset) ) { data <- ncvar_get(nc,varname)[,,lev.idx] }
         else                   { data <- ncvar_get(nc,varname)[idx.lon,idx.lat,lev.idx] }
      } else {
         if ( is.null(subset) ) { data <- ncvar_get(nc,varname)[,,lev.idx,time.idx] }
         else                   { data <- ncvar_get(nc,varname)[idx.lon,idx.lat,lev.idx,time.idx] }
      }
   }
   # Prepare return value
   dx <- abs(mean(diff(lons))/2)
   dy <- abs(mean(diff(lats))/2)
   # If we only have one layer
   if ( length(dim(data)) == 2 ) {
      data <- raster(t(data),xmn=min(lons)-dx,xmx=max(lons)+dx,ymn=min(lats)-dy,ymx=max(lats)+dy,
               crs=crs("+proj=longlat +datum=WGS84 +ellps=WGS84 +no_defs +towgs84=0,0,0"))
   # For multiple layers, multiple steps requested by the user
   } else {
      empty <- raster(ncols=dim(data)[1L],nrows=dim(data)[2L],
               xmn=min(lons)-dx,xmx=max(lons)+dx,ymn=min(lats)-dy,ymx=max(lats)+dy,
               crs=crs("+proj=longlat +datum=WGS84 +ellps=WGS84 +no_defs +towgs84=0,0,0"))
      res <- NULL
      for ( l in 1:dim(data)[3L] ) {
         tmp <- empty; values(tmp) <- as.vector(data[,,l])
         res <- `if`(is.null(res),tmp,stack(res,tmp))
      }
      data <- res; rm(res)
   }
   # If input 'subset' was of class raster::RasterLayer: resample
   if ( class(subset) == "RasterLayer" ) data <- resample(data,subset)
   if ( identical(time.idx,FALSE) ) {
      names(data) <- sprintf("%s_%03d",strftime(init,"%Y%m%d"),ncsteps)
   } else {
      names(data) <- sprintf("%s_%03d",strftime(init,"%Y%m%d"),ncsteps[time.idx])
   }
   data
}
