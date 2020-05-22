test_that("model sizes work for W2IP", {
  set.seed(84370158)
  
  n <- 100
  p <- 10
  s <- 1000
  
  x <- matrix( rnorm( p * n ), nrow = n, ncol = p )
  x_ <- t(x)
  beta <- (1:p)/p
  y <- x %*% beta + rnorm(n)
  post_beta <- matrix(beta, nrow=p, ncol=s) + rnorm(p*s, 0, 0.1)
  post_mu <- x %*% post_beta
  transp <- "hilbert"
  model.size <- c(2,4,8)
    
  test <- W2IP(X = x, Y = post_mu, theta = post_beta, transport.method = transp, 
               infimum.maxit = 10, 
               tol = 1e-7, solution.method = "cone",
               display.progress = FALSE,model.size = model.size)
  testthat::expect_equal(test$nzero, model.size)
  
})


test_that("model sizes work for gurobi", {
  set.seed(84370158)
  
  n <- 100
  p <- 10
  s <- 1000
  
  x <- matrix( rnorm( p * n ), nrow = n, ncol = p )
  x_ <- t(x)
  beta <- (1:p)/p
  y <- x %*% beta + rnorm(n)
  post_beta <- matrix(beta, nrow=p, ncol=s) + rnorm(p*s, 0, 0.1)
  post_mu <- x %*% post_beta
  transp <- "hilbert"
  model.size <- c(2,4,8)
  
  # debugonce(W2IP)
  test <- W2IP(X = x, Y = post_mu, theta = post_beta, transport.method = transp, 
               infimum.maxit = 10, 
               tol = 1e-7, solution.method = "gurobi",
               display.progress = FALSE,model.size = model.size)
  testthat::expect_equal(test$nzero, model.size)
  
})

test_that("model sizes work for mosek", {
  set.seed(84370158)
  
  n <- 100
  p <- 10
  s <- 1000
  
  x <- matrix( rnorm( p * n ), nrow = n, ncol = p )
  x_ <- t(x)
  beta <- (1:p)/p
  y <- x %*% beta + rnorm(n)
  post_beta <- matrix(beta, nrow=p, ncol=s) + rnorm(p*s, 0, 0.1)
  post_mu <- x %*% post_beta
  transp <- "hilbert"
  model.size <- c(2,4,8)
  
  # debugonce(W2IP)
  test <- W2IP(X = x, Y = post_mu, theta = post_beta, transport.method = transp, 
               infimum.maxit = 10, 
               tol = 1e-7, solution.method = "mosek",
               display.progress = FALSE,model.size = model.size)
  testthat::expect_equal(test$nzero, model.size)
  
})


# test_that("times work for W2IP LP", { # not work for lpsolve
#   set.seed(84370158)
#   
#   n <- 1000
#   p <- 500
#   s <- 1000
#   
#   x <- matrix( rnorm( p * n ), nrow = n, ncol = p )
#   x_ <- t(x)
#   beta <- (1:p)/p
#   y <- x %*% beta + rnorm(n)
#   post_beta <- matrix(beta, nrow=p, ncol=s) + rnorm(p*s, 0, 0.1)
#   post_mu <- x %*% post_beta
#   transp <- "exact"
#   
#   time.start <- proc.time()
#   testthat::expect_warning(W2IP(X = x, Y = post_mu, theta = post_beta, transport.method = transp, 
#                infimum.maxit = 10, 
#                tol = 1e-7,model.size = 1, solution.method = "lp",
#                display.progress = FALSE, control = list(tm_limit = 1)))
#   time.end <- proc.time()
#   testthat::expect_lt((time.end - time.start)[3], 100)
#   
# })
