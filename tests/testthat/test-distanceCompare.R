test_that("distance compare gives correct values for wass", {
  
  set.seed(84370158)
  
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
  compost <- unlist(sapply(out, function(o) sapply(o$theta, function(tt)  SparsePosterior::wasserstein(tt, post_beta, 2, 2, "colwise","exact"))))
  commean <- unlist(sapply(out, function(o) sapply(o$eta, function(tt)  SparsePosterior::wasserstein(tt, post_mu, 2, 2, "colwise","exact"))))
  
  testthat::expect_equal(dist$posterior$dist, compost)
  testthat::expect_equal(dist$mean$dist, commean)
})

test_that("distance compare gives correct values for mse", {
  
  set.seed(84370158)
  
  n <- 128
  p <- 10
  s <- 100
  
  x <- matrix( rnorm( p * n ), nrow = n, ncol = p )
  x_ <- t(x)
  beta <- (1:p)/p
  y <- x %*% beta + rnorm(n)
  post_beta <- matrix(beta, nrow=p, ncol=s) + rnorm(p*s, 0, 0.1)
  post_mu <- x %*% post_beta
  mu <- x %*% beta
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
  mse <- distCompare(out, list(posterior = beta, mean =  mu), p = 2, ground_p = 2, quantity = c("posterior", "mean"), method = "mse")
  compost <- unlist(sapply(out, function(o) sapply(o$theta, function(tt)  SparsePosterior::wasserstein(tt, as.matrix(beta), 2, 2, "colwise","exact"))))^2/p
  commean <- unlist(sapply(out, function(o) sapply(o$eta, function(tt)  SparsePosterior::wasserstein(tt,  as.matrix(mu), 2, 2, "colwise","exact"))))^2/n
  
  testthat::expect_equal(mse$posterior$dist, compost)
  testthat::expect_equal(mse$mean$dist, commean)
})

test_that("distance compare gives correct group names", {

  set.seed(84370158)

  n <- 128
  p <- 10
  s <- 100

  x <- matrix( rnorm( p * n ), nrow = n, ncol = p )
  x_ <- t(x)
  beta <- (1:p)/p
  y <- x %*% beta + rnorm(n)
  post_beta <- matrix(beta, nrow=p, ncol=s) + rnorm(p*s, 0, 0.1)
  post_mu <- x %*% post_beta
  mu <- x %*% beta
  transp <- "exact"
  model.size <- c(2,4,8)

  test <- W2IP(X = x, Y = post_mu, theta = post_beta, transport.method = transp,
               infimum.maxit = 10,
               tol = 1e-7, solution.method = "cone",
               display.progress = FALSE,model.size = model.size)

  proj <- W2L1(x, post_mu, post_beta, penalty = "lasso", method = "projection", infimum.maxit = 1)
  sel <- W2L1(x, post_mu, post_beta, penalty = "selection.lasso", method = "selection.variable")

  out <- list(Test=test, Projection = proj, Selection = sel)
  # debugonce(distCompare)
  mse <- distCompare(out, list(posterior = beta, mean =  mu), p = 2, ground_p = 2, quantity = c("posterior", "mean"), method = "mse")
  compost <- unlist(sapply(out, function(o) sapply(o$theta, function(tt)  SparsePosterior::wasserstein(tt, as.matrix(beta), 2, 2, "colwise","exact"))))^2/p
  commean <- unlist(sapply(out, function(o) sapply(o$eta, function(tt)  SparsePosterior::wasserstein(tt,  as.matrix(mu), 2, 2, "colwise","exact"))))^2/n

  expectnames <- c('Test', 'Test', 'Test', 'Projection', 'Projection', 'Projection', 
                   'Projection', 'Projection', 'Projection', 'Projection', 'Projection', 
                   'Projection', 'Projection', 'Selection', 'Selection', 'Selection', 
                   'Selection', 'Selection', 'Selection', 'Selection', 'Selection', 
                   'Selection', 'Selection')
  testthat::expect_equivalent(mse$posterior$dist, compost)
  testthat::expect_equivalent(mse$mean$dist, commean)
  testthat::expect_equal(as.character(mse$posterior$groups), expectnames)
  testthat::expect_equal(as.character(mse$mean$groups), expectnames)
})
