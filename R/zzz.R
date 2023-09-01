.onAttach <- function(libname, pkgname) {
  packageStartupMessage(paste0("This is {flapper} v.", utils::packageVersion("flapper"), ". For an overview, see ?flapper. For support, contact edward.lavender@eawag.ch."))
}
