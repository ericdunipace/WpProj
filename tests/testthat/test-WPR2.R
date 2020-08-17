wpr2.data <- function(n, p, s) {
  x <- matrix( rnorm( p * n ), nrow = n, ncol = p )
  x_ <- t(x)
  beta <- (1:p)/p
  y <- x %*% beta + rnorm(n)
  post_beta <- matrix(beta, nrow=p, ncol=s) + rnorm(p*s, 0, 0.1)
  post_mu <- x %*% post_beta
  transp <- "exact"
  model.size <- c(2,4,8)
  
  test <- W2IP(X = x, Y = post_mu, theta = post_beta, transport.method = transp, 
               infimum.maxit = 10, 
               tol = 1e-7, solution.method = "cone",
               display.progress = FALSE,model.size = model.size)
  
  proj <- W2L1(x, post_mu, post_beta, penalty = "lasso", method = "projection", infimum.maxit = 1)
  sel <- W2L1(x, post_mu, post_beta, penalty = "selection.lasso", method = "selection.variable")
  
  out <- list(test, proj, sel)
  dist <- distCompare(out, list(posterior = post_beta, mean = post_mu), p = 2, ground_p = 2, quantity = c("posterior", "mean"))
  return(dist)
}

wpr2.prep <- function(n, p, s) {
  out <- wpr2.data(n,p,s)
  
  
  r2 <- WPR2(nu = out, p = 2, method = "exact")
  return(r2)
}

testthat::test_that("WPR2 works", {
  set.seed(203402)
  
  n <- 128
  p <- 10
  s <- 100
  
  x <- matrix( rnorm( p * n ), nrow = n, ncol = p )
  x_ <- t(x)
  beta <- (1:p)/p
  y <- x %*% beta + rnorm(n)
  post_beta <- matrix(beta, nrow=p, ncol=s) + rnorm(p*s, 0, 0.1)
  post_mu <- x %*% post_beta
  transp <- "exact"
  model.size <- c(2,4,8)
  
  test <- W2IP(X = x, Y = post_mu, theta = post_beta, transport.method = transp, 
               infimum.maxit = 10, 
               tol = 1e-7, solution.method = "cone",
               display.progress = FALSE,model.size = model.size)
  
  proj <- W2L1(x, post_mu, post_beta, penalty = "lasso", method = "projection", infimum.maxit = 1)
  sel <- W2L1(x, post_mu, post_beta, penalty = "selection.lasso", method = "selection.variable")
  
  out <- list(test, proj, sel)
  
  
  dist <- distCompare(out, list(posterior = post_beta, mean = post_mu), p = 2, ground_p = 2, quantity = c("posterior", "mean"))
  
  r2 <- WPR2(Y = post_mu, nu = dist, p = 2, method = "exact")
  r2 <- WPR2(nu = dist, p = 2, method = "exact")
  
  maxes <- tapply(dist$mean$dist, dist$mean$groups, max)
  r2_check <- 1 - dist$mean$dist^2/maxes[as.numeric(dist$mean$groups)]^2
    
  r2_mat <- WPR2.matrix(post_mu, test$eta[[1]], p = 2, method  ="exact")
  r2_mat_check <- 1 - (limbs::wasserstein(post_mu, test$eta[[1]],
                                         p = 2, ground_p = 2,
                                         method = "exact", 
                                         observation.orientation = "colwise")^2/
    limbs::wasserstein(post_mu, 
                       matrix(colMeans(post_mu), nrow(post_mu),
                          ncol(post_mu), byrow=TRUE),
                       p = 2, ground_p = 2,
                       method = "exact", 
                       observation.orientation = "colwise")^2)
  
  testthat::expect_silent(r2_limbs <- WPR2.list(post_mu, out, p = 2, method  ="exact"))
  testthat::expect_silent(r2_limbs <- WPR2(post_mu, out, p = 2, method  ="exact"))
  
  names(out) <- c("BP", "L2", "relaxed bp")
  r2_limbs <- WPR2(post_mu, out, p = 2, method  ="exact")
  r2_limbs_check <- 1 - (limbs::wasserstein(post_mu, proj$eta[[1]],
                                            p = 2, ground_p = 2,
                                            method = "exact", 
                                            observation.orientation = "colwise")^2/
                           limbs::wasserstein(post_mu, 
                                              matrix(colMeans(post_mu), nrow(post_mu),
                                                     ncol(post_mu), byrow=TRUE),
                                              p = 2, ground_p = 2,
                                              method = "exact", 
                                              observation.orientation = "colwise")^2)
  
  testthat::expect_equivalent(r2$r2, r2_check)
  testthat::expect_equivalent(r2_mat[1,1], r2_mat_check)
  testthat::expect_equivalent(r2_limbs$r2[r2_limbs$groups == "L2"][1], r2_limbs_check, )
})

testthat::test_that("WPR2 combining works", {
  set.seed(203402)
  
  n <- 128
  p <- 10
  s <- 100
  
  out1 <- wpr2.prep(n,p,s)
  out2 <- wpr2.prep(n,p,s)
  # debugonce(distCompare)
  
  comb <- combine.WPR2(out1,out2)
  comb2 <- combine.WPR2(list(out1,out2))
  
  testthat::expect_equal(comb, comb2)
  
})

testthat::test_that("WPR2 plotting works", {
  set.seed(203402)
  
  n <- 128
  p <- 10
  s <- 100
  reps <- 10
  out <- lapply(1:reps, function(i) wpr2.prep(n,p,s))
  # debugonce(combine.WPR2)
  comb <- combine.WPR2(out)
  # debugonce(plot.WPR2)
  p <- plot(comb)
  testthat::expect_true(ggplot2::is.ggplot(p))
})