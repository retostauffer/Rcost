# -------------------------------------------------------------------
# Forces the system to use UTC timezone as soon as the package
# is attached in an R session.
# -------------------------------------------------------------------
.onAttach <- function(libname, pkgname) {
   Sys.setenv("TZ"="UTC")
   info <- sprintf("\n    Warning: R timezone has been set to UTC!\n\n")
   packageStartupMessage(info)
}
