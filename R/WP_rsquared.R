# WPR2 <- function(Y, nu, p = 2, method = "exact",...) {UseMethod("WPR2")}

WPR2.matrix <- function(Y, nu, p = 2, method = "exact", ...) {
  
  stopifnot(p >= 1)
  
  n <- nrow(Y)
  d <- ncol(Y)
  
  wp_mod <- limbs::wasserstein(Y, nu, p = p, ground_p = p,
                           method = method, ...)^p
  mu <- matrix(colMeans(Y), n, d, byrow=TRUE)
  # wp_base <- if(method == "exact") {
  #   mean(colSums((Y - mu)^p))
  # } else {
  #   limbs::wasserstein(Y, mu, p = p, 
  #                      ground_p = p,
  #                      method = method, 
  #                      ...)^p
  # }
  wp_base <- limbs::wasserstein(Y, mu, p = p, 
                                ground_p = p,
                                method = method, 
                                ...)^p
  
  r2 <- pmax(1 - wp_mod/wp_base,0)
  
  return(r2)
  
}

WPR2.distcompare <- function(Y=NULL, nu, ...) {
  
  stopifnot(inherits(nu, "distcompare"))
  
  df <- nu$mean
  p <- nu$p
  method <- df$method
  
  stopifnot(p >= 1)
  
  
  if(!is.null(Y)) {
    stopifnot(inherits(Y, "matrix"))
    meth.table <- table(method)
    method.use <- names(meth.table)[which.max(meth.table)]
    wass.args <- list(X = Y, Y = as.matrix(rowMeans(Y)),
                      p = as.numeric(p), method = method.use,
                      ...)
    wass.args <- wass.args[!duplicated(names(wass.args))]
    argn <- lapply(names(wass.args), as.name)
    names(argn) <- names(wass.args)
    
    wass.call <- as.call(c(list(as.name("wasserstein")), argn))
    
    max_vals <- eval(wass.call, envir = wass.args)
    max_vec <- rep(max_vals, length(df$dist))
    base <- "dist.from.expectation"
  } else {
    max_vals <- tapply(df$dist, df$groups, max)
    max_vec <- max_vals[as.numeric(df$groups)]
    base <- "dist.from.null"
  }
  
  
  
  r2 <- pmax(1- df$dist^p/max_vec^p, 0)
  
  df$dist <- r2
  df$p <- p
  df$base <- base
  colnames(df)[1] <- "r2"
  
  
  return(df)
  
}

setGeneric("WPR2", function(Y = NULL, nu, p = 2, method = "exact",...) UseMethod("WPR2"))
setMethod("WPR2", c("Y" = "matrix", "nu" = "matrix"), WPR2.matrix)
setMethod("WPR2", c("nu" = "distcompare"), WPR2.distcompare)

