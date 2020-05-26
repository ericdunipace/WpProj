testthat::test_that("multiplication works", {
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
  # debugonce(distCompare)
  dist <- distCompare(out, list(posterior = post_beta, mean = post_mu), p = 2, ground_p = 2, quantity = c("posterior", "mean"))
  
  r2 <- WPR2(Y=dist, p = 2, method = "exact")
  r2 <- WPR2(dist, p = 2, method = "exact")
  
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
  
  testthat::expect_equivalent(r2$r2, r2_check)
  testthat::expect_equivalent(r2_mat, r2_mat_check)
})
