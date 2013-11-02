## .First.lib <- function(lib, pkg) library.dynam("waveslim", pkg, lib)

.onAttach <- function(lib, pkg) {
  txt <- paste("\n",
               pkg,
               ": Wavelet Method for 1/2/3D Signals (version = ",
               packageDescription(pkg, lib)[["Version"]],
               ")\n",
               sep="")
  packageStartupMessage(txt)
}
