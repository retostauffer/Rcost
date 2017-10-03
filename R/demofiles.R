# -------------------------------------------------------------------
# - NAME:        demofiles.R
# - AUTHOR:      Reto Stauffer
# - DATE:        2017-10-03
# -------------------------------------------------------------------
# - DESCRIPTION:
# -------------------------------------------------------------------
# - EDITORIAL:   2017-10-03, RS: Created file on thinkreto.
# -------------------------------------------------------------------
# - L@ST MODIFIED: 2017-10-03 16:04 on thinkreto
# -------------------------------------------------------------------


demofiles <- function( silent = FALSE ) {

   sfcfile <- system.file("extdata", "ITERIM_sfc.nc", package = "Rcost")
   plfile  <- system.file("extdata", "ITERIM_pl.nc",  package = "Rcost")

   if ( ! silent ) {
      cat(sprintf("Surface file:        %s\n", sfcfile))
      cat(sprintf("Pressure level file: %s\n", plfile))
   }

   return( list( sfc = sfcfile, pl = plfile ) )

}

