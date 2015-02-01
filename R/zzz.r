##########################################################################################################################################
# PrimSRC
##########################################################################################################################################

.onAttach <- function(libname, pkgname) {
    SSver <- read.dcf(file=system.file("DESCRIPTION", package=pkgname),
                      fields="Version")
    packageStartupMessage("\n")
    packageStartupMessage(paste(pkgname, SSver))
    packageStartupMessage("\n")
    packageStartupMessage("Type PrimSRC.news() to see new features, changes, and bug fixes.")
    packageStartupMessage("\n")
}
