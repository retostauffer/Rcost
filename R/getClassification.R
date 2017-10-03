# -------------------------------------------------------------------
# This is the main function which does the classification.
# -------------------------------------------------------------------
getClassification <- function( nc, date, steps, subset, pwclim, weights, nsector=8, maindirthreshold=0.4, zamg=FALSE, verbose=0 ) {

   # Store some execution times
   time <- list()

   # Specify forecast initial time
   init <- try(as.Date(date))
   if ( any("try-error" %in% class(init)) )
      stop("Wrong input format for \"date\" when calling getClassification")

   # IF steps missing: take the steps from the ncs file
   if ( missing(steps) ) {
      tmp <- list()
      for ( i in seq_along(nc) ) tmp[[i]] <- datetimeinfo(nc[[i]],init)$step
      # Extracting overlapping time steps (which are available in all nc files)
      steps <- table(unlist(tmp))
      steps <- unique(as.integer(names(steps))[steps == length(nc)])
      if ( length(steps) == 0 ) stop("Sorry, no overlapping steps in the netcdf files")
      rm(tmp)
   }

   # If subset is missing: take the ZAMG Austria definition.
   # No re-gridding, take the ncfiles as they are.
   # examples: subset <- extent(c(1.5,25.5 ,40.5 ,54.5))
   # examples: subset <- raster(nrow=15,ncol=25,xmn=1,xmx=26,ymn=40,ymx=55)
   # If not set: take biggest overlapping subset possible
   if ( missing(subset) ) {
      lon <- lat <- list()
      for ( n in nc ) {
         lon[[length(lon)+1]] <- ncvar_get(n,"longitude")
         lat[[length(lat)+1]] <- ncvar_get(n,"latitude")
      }
      subset <- extent( c(max(unlist(lapply(lon,min))), min(unlist(lapply(lon,max))),
                          max(unlist(lapply(lat,min))), min(unlist(lapply(lat,max)))) )
   }

   time[['getdata']] <- Sys.time()
   raster_tcw   <- getdata(nc,init,"tcw",steps,           subset=subset, silent=T)
   raster_u700  <- getdata(nc,init,"u",  steps,level=700, subset=subset)
   raster_v700  <- getdata(nc,init,"v",  steps,level=700, subset=subset)
   raster_z500  <- getdata(nc,init,"z",  steps,level=500, subset=subset)
   raster_z925  <- getdata(nc,init,"z",  steps,level=925, subset=subset)
   time[['getdata']] <- Sys.time() - time[['getdata']]

   # If tcw not in the netcdf input file:
   if ( is.null(raster_tcw) ) {
      raster_tcw <- raster_u700
      values(raster_tcw) <- NA
   }

   #warning("Multiply geopotential height with 1000.")
   #raster_z500 <- raster_z500 * 1000.
   #raster_z925 <- raster_z925 * 1000.

   warning("Divided z500 by factor of 10 (as in cost733 routine)")
   raster_z500 <- raster_z500 / 10.

   # ----------------------------------------------------------------
   # Classification
   # ----------------------------------------------------------------

   # Check data loaded
   if ( ! all(dim(raster_tcw) == dim(raster_u700)) |
        ! all(dim(raster_tcw) == dim(raster_v700)) |
        ! all(dim(raster_tcw) == dim(raster_z500)) |
        ! all(dim(raster_tcw) == dim(raster_z925)) )
      stop("Dimension mismatch in loaded raster objects.")
   if ( ! extent(raster_tcw) == extent(raster_u700) |
        ! extent(raster_tcw) == extent(raster_v700) |
        ! extent(raster_tcw) == extent(raster_z500) |
        ! extent(raster_tcw) == extent(raster_z925) )
      stop("Extent mismatch in loaded raster objects.")

   # Getting grid information
   coord <- getlonlat( raster_tcw )

   # Prepare fortran inputs
   nj   <- as.integer(coord$nx) # columns: longitudes (x)
   ni   <- as.integer(coord$ny) # rows: latitudes (y)
   lons <- sort( as.numeric(coord$lons) ) # sort, we also transpose the data
   lats <- sort( as.numeric(coord$lats) ) # sort, we also transpose the data

   # The data
   time[['asarray']] <- Sys.time()
   u700 = as.array(raster_u700)
   v700 = as.array(raster_v700)
   z500 = as.array(raster_z500)
   z925 = as.array(raster_z925)
   pw   = as.array(raster_tcw)
   time[['asarray']] <- Sys.time() - time[['asarray']]

   # Number of wind sectors
   nsector = as.integer(nsector)

   # pw climatology, corrently set to -999 (missing)
   if ( missing(pwclim) ) pwclim <- rep(-999.,366)
   pwclim <- as.numeric(pwclim)

   # The weighting factors for the three different areas.
   # If missing, the ZAMG specification will be used.
   if ( missing(weights) ) {
      weights <- as.numeric(c(15,2,1))
   } else {
      if ( ! length(weights) == 3 ) stop("\"weights\" has to be of length 3.")
      weights <- as.numeric(weights)
   }

   # Getting dates and ydays
   # from the raster names (consist of <initdate>_<step>). Note
   # that "dates" are then of format (int)YYYYmmddHH, no minutes.
   # Adding minutes yields integers bigger than allowed.
   inits <- regmatches(names(raster_tcw),regexpr("[0-9]{8}", names(raster_tcw)))
   steps <- regmatches(names(raster_tcw),regexpr("[0-9]{1,}$",names(raster_tcw)))
   tmp   <- strptime(inits,"%Y%m%d")+as.numeric(steps)*3600
   dates <- as.integer(strftime(tmp,"%Y%m%d%H"))
   ydays <- as.integer(as.POSIXlt(tmp)$yday)
   nt <- as.integer(length(dates))

   intresult <- matrix(as.integer(-88),  ncol=6,nrow=nt)
   realresult <- matrix(as.numeric(-88), ncol=3,nrow=nt)

   # Loading weight mask
   if ( ! zamg ) {
      w <- weightmask(lons,lats,weights)
   } else {
      message("\"zamg\" set to true in getClassification: \"nsector\" will not be used.")
      w <- weightmaskZAMG(lons,lats,weights)
   }

   #cat(sprintf("   R weight mask sum: %d\n",sum(w)))
   time[['fortran']] <- Sys.time()
   x <- .Fortran("wlk733",dates,ydays,nt,lats,ni,lons,nj,
        u700,v700,z500,z925,pw,pwclim,nsector,w,
        as.numeric(maindirthreshold[1L]),as.integer(zamg),intresult,realresult,
        as.integer(verbose),
        PACKAGE="Rcost")
   res <- data.frame(cbind(x[[18]],x[[19]]))
   rownames(res) <- NULL
   time[['fortran']] <- Sys.time()-time[['fortran']]

   # Result contains
   # date: YYYYmmddHH: time for which the analysis is valid
   # yday: day of the year for YYYYmmddHH
   # sector: main wind sector, analysis
   # cyclonic925: 925hPa analysis, 0=not cyclonic (or anticyclonic), 1=cyclonic
   # cyclonic500: 500hPa analysis, 0=not cyclonic (or anticyclonic), 1=cyclonic
   # wet: -99 missing as pwclim is missing, 0=dry, 1=wet
   # c925/c500: the curvature values (real) which are used internally
   # pwi: the real value used for the internal computation for "dry/wet"
   names(res) <- c("date","yday","sector","cyclonic925","cyclonic500","wet","c925","c500","pwi")
   res$date <- as.POSIXct(strptime(res$date,"%Y%m%d%H"))
   #print(x)

   z <- zoo::zoo( subset(res,select=-c(date)), res$date ) 
   attr(z,'init')    <- init
   attr(z,'created') <- Sys.time()

   # Some stats
   if ( verbose ) {
      cat(sprintf("   Some execution time summaries:\n"))
      for ( n in names(time) ) { cat(sprintf("   - %-10s  ",n)); print(time[[n]]) }
   }
   
   return(z)

}











