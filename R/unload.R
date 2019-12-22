.onUnload <- function (libpath) {
  library.dynam.unload("SparsePosterior", libpath)
}