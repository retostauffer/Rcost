# -------------------------------------------------------------------
# - NAME:        weightmask.R
# - AUTHOR:      Reto Stauffer
# - DATE:        2017-02-26
# -------------------------------------------------------------------
# - DESCRIPTION: Weightmask
# -------------------------------------------------------------------
# - EDITORIAL:   2017-02-26, RS: Created file on thinkreto.
# -------------------------------------------------------------------
# - L@ST MODIFIED: 2017-10-03 18:14 on thinkreto
# -------------------------------------------------------------------

# lon/lat: vector of longitudes and latitudes
# weights: weights for core/mid/margin
weightmask <- function(lon,lat,weights=c(15,2,1),stripsize=5) {

   if ( ! is.numeric(weights) | ! length(weights) == 3 )
      stop("Misspecification of \"weights\" in function weightmask.")

   # Compute strip width
   ni <- length(lat)
   nj <- length(lon)
   wi <- floor(ni / stripsize)
   wj <- floor(nj / stripsize)

   # Give 'core' weight to all cells
   res <- matrix(weights[1L],nrow=ni,ncol=nj)

   # mid
   res[1:(wi*2),]       <- weights[2L]
   res[(ni-wi*2+1):ni,] <- weights[2L]
   res[,1:(wj*2)]       <- weights[2L]
   res[,(nj-wj*2+1):nj] <- weights[2L]

   # margin
   res[1:wi,]           <- weights[3L]
   res[(ni-wi+1):ni,]   <- weights[3L]
   res[,1:wj]           <- weights[3L]
   res[,(nj-wj+1):nj]   <- weights[3L]

   return(res)
}



# lon/lat: vector of longitudes and latitudes
# weights: weights for core/mid/margin
# Try to reproduce the ZAMG weight mask as proposes by Thomas Kennert, ZAMG
weightmaskZAMG <- function(lon,lat,weights=c(15,2,1),...) {

   lon <- sort(lon); lat <- sort(lat)
   message("Using ZAMG weight mask specification. Please think about it!")
   if ( ! is.numeric(weights) | ! length(weights) == 3 )
      stop("Misspecification of \"weights\" in function weightmask.")

   res <- matrix(weights[3L],ncol=length(lon),nrow=length(lat))

   core.lon <- which( lon >=  9.0 & lon <= 18.0 )
   core.lat <- which( lat >= 45.0 & lat <= 50.0 )

   mid.lon <- which( lon >=  6.0 & lon <= 21.0 )
   mid.lat <- which( lat >= 42.0 & lat <= 52.5 )

   res[mid.lat,mid.lon]   <- weights[2L]
   res[core.lat,core.lon] <- weights[1L]

   return(res)
}
