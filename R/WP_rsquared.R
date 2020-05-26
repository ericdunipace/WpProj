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

WPR2.distcompare <- function(Y, nu, ...) {
  
  stopifnot(inherits(Y, "distcompare"))
  
  stopifnot(p >= 1)
  
  df <- Y$mean
  p <- Y$p
  
  max_vals <- tapply(df$dist, df$groups, max)
  
  max_vec <- max_vals[as.numeric(df$groups)]
  
  r2 <- pmax(1- df$dist^p/max_vec^p, 0)
  
  df$dist <- r2
  df$p <- p
  colnames(df)[1] <- "r2"
  
  return(df)
  
}

setGeneric("WPR2", function(Y, nu, p = 2, method = "exact",...) UseMethod("WPR2"))
setMethod("WPR2", c("Y" = "matrix", "nu" = "matrix"), WPR2.matrix)
setMethod("WPR2", c("Y" = "distcompare"), WPR2.distcompare)

