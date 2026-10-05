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
  nvars <- c(2,4,8)
    
  test <- W2IP(X = x, Y = post_mu, theta = post_beta, transport.method = transp, 
               infimum.maxit = 10, 
               tol = 1e-7, solver = "cone",
               display.progress = FALSE,nvars = nvars)
  testthat::expect_equal(test$nzero, nvars)
  
})


# test_that("model sizes work for gurobi", {
#   check_gurobi()
#   
#   set.seed(84370158)
#   
#   n <- 100
#   p <- 10
#   s <- 1000
#   
#   x <- matrix( rnorm( p * n ), nrow = n, ncol = p )
#   x_ <- t(x)
#   beta <- (1:p)/p
#   y <- x %*% beta + rnorm(n)
#   post_beta <- matrix(beta, nrow=p, ncol=s) + rnorm(p*s, 0, 0.1)
#   post_mu <- x %*% post_beta
#   transp <- "hilbert"
#   nvars <- c(2,4,8)
#   
#   # debugonce(W2IP)
#   WpProj:::check_gurobi()
#   test <- W2IP(X = x, Y = post_mu, theta = post_beta, transport.method = transp, 
#                infimum.maxit = 10, 
#                tol = 1e-7, solver = "gurobi",
#                display.progress = FALSE,nvars = nvars)
#   testthat::expect_equal(test$nzero, nvars)
#   
# })

test_that("model sizes work for mosek", {
  check_mosek()
  set.seed(84370158)
  
  n <- 100
  p <- 10
  s <- 1000
  
  x <- matrix( stats::rnorm( p * n ), nrow = n, ncol = p )
  x_ <- t(x)
  beta <- (1:p)/p
  y <- x %*% beta + stats::rnorm(n)
  post_beta <- matrix(beta, nrow=p, ncol=s) + stats::rnorm(p*s, 0, 0.1)
  post_mu <- x %*% post_beta
  transp <- "hilbert"
  nvars <- c(2,4,8)
  
  WpProj:::check_mosek()
  test <- W2IP(X = x, Y = post_mu, theta = post_beta, transport.method = transp, 
               infimum.maxit = 10, 
               tol = 1e-7, solver = "mosek",
               display.progress = FALSE,nvars = nvars)
  testthat::expect_equal(test$nzero, nvars)
  
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
#                tol = 1e-7,nvars = 1, solver = "lp",
#                display.progress = FALSE, control = list(tm_limit = 1)))
#   time.end <- proc.time()
#   testthat::expect_lt((time.end - time.start)[3], 100)
#   
# })

test_that("model sizes work for highs and match lpsolve", {
  check_highs()
  set.seed(84370158)
  
  n <- 100
  p <- 10
  s <- 1000
  
  x <- matrix( stats::rnorm( p * n ), nrow = n, ncol = p )
  beta <- (1:p)/p
  post_beta <- matrix(beta, nrow=p, ncol=s) + stats::rnorm(p*s, 0, 0.1)
  post_mu <- x %*% post_beta
  transp <- "hilbert"
  nvars <- c(2,4,8)
  
  test_highs <- W2IP(X = x, Y = post_mu, theta = post_beta, transport.method = transp, 
                     infimum.maxit = 10, 
                     tol = 1e-7, solver = "highs",
                     display.progress = FALSE, nvars = nvars)
  test_lp <- W2IP(X = x, Y = post_mu, theta = post_beta, transport.method = transp, 
                  infimum.maxit = 10, 
                  tol = 1e-7, solver = "lp",
                  display.progress = FALSE, nvars = nvars)
  testthat::expect_equal(test_highs$nzero, nvars)
  testthat::expect_equal(test_highs$beta, test_lp$beta)
  
})

test_that("model sizes work for scip and match lpsolve", {
  check_scip()
  set.seed(84370158)
  
  n <- 100
  p <- 10
  s <- 1000
  
  x <- matrix( stats::rnorm( p * n ), nrow = n, ncol = p )
  beta <- (1:p)/p
  post_beta <- matrix(beta, nrow=p, ncol=s) + stats::rnorm(p*s, 0, 0.1)
  post_mu <- x %*% post_beta
  transp <- "hilbert"
  nvars <- c(2,4,8)
  
  test_scip <- W2IP(X = x, Y = post_mu, theta = post_beta, transport.method = transp, 
                    infimum.maxit = 10, 
                    tol = 1e-7, solver = "scip",
                    display.progress = FALSE, nvars = nvars)
  test_lp <- W2IP(X = x, Y = post_mu, theta = post_beta, transport.method = transp, 
                  infimum.maxit = 10, 
                  tol = 1e-7, solver = "lp",
                  display.progress = FALSE, nvars = nvars)
  testthat::expect_equal(test_scip$nzero, nvars)
  testthat::expect_equal(test_scip$beta, test_lp$beta)
  
  # user supplied control, including a time limit, is passed through
  test_scip_tl <- W2IP(X = x, Y = post_mu, theta = post_beta, transport.method = transp, 
                       infimum.maxit = 10, 
                       tol = 1e-7, solver = "scip",
                       display.progress = FALSE, nvars = nvars,
                       control = list(verbose = FALSE, time_limit = 30))
  testthat::expect_equal(test_scip_tl$beta, test_lp$beta)
})

test_that("augmented lagrangian algorithm matches the exact algorithm", {
  set.seed(84370158)
  
  n <- 100
  p <- 12
  s <- 200
  
  x <- matrix( stats::rnorm( p * n ), nrow = n, ncol = p ) %*% chol(0.5^abs(outer(1:p, 1:p, "-")))
  post_beta <- matrix(stats::rnorm(p), nrow=p, ncol=s) + stats::rnorm(p*s, 0, 0.3)
  post_mu <- x %*% post_beta
  
  ss <- sufficientStatistics(x, post_mu, post_beta, 
                             list(same = TRUE, method = "selection.variable",
                                  transport.method = "hilbert", epsilon = 0.05, niter = 0L))
  
  # single solves of the penalized problem give the hard constrained optimum
  check_scip()
  for (k in c(1, 3, 6, 11)) {
    QP <- qp_w2(ss$XtX, ss$XtY, 1)
    QP$constraints$rhs[1] <- k
    Q <- as.matrix(QP$objective$Q)
    L <- as.numeric(as.matrix(QP$objective$L))
    obj <- function(a) c(0.5 * crossprod(a, Q %*% a)) + sum(L * a)
    
    pen <- augmented_lagrangian_penalty(QP)
    a_lagr  <- scip_solver(add_cardinality_penalty(QP, pen$nu, pen$rho))[1:p]
    a_hard <- scip_solver(QP)
    testthat::expect_equal(sum(a_lagr), k)
    testthat::expect_equal(obj(a_lagr), obj(a_hard))
  }
  QP$constraints$rhs[1] <- p
  testthat::expect_null(augmented_lagrangian_penalty(QP))
  
  # full algorithm with each free solver
  nvars <- c(2, 4, 8)
  fit <- function(solver, algorithm) {
    W2IP(X = x, Y = post_mu, theta = post_beta, transport.method = "hilbert", 
         infimum.maxit = 10, tol = 1e-7, nvars = nvars,
         solver = solver, algorithm = algorithm)
  }
  test_exact <- fit("lp", "exact")
  testthat::expect_equal(test_exact$nzero, nvars)
  for (solver in c("scip", "lp", "highs")) {
    if (solver == "highs") check_highs()
    test_lagr <- fit(solver, "augmented.lagrangian")
    testthat::expect_equal(test_lagr$beta, test_exact$beta, info = solver)
  }
  # ECOS's branch and bound isn't always optimal so only check the sizes
  testthat::expect_equal(fit("cone", "augmented.lagrangian")$nzero, nvars)
})

test_that("mosek problem includes the augmented lagrangian slack variables", {
  testthat::skip_if_not_installed("Rmosek")
  set.seed(1)
  p <- 6
  X <- matrix(stats::rnorm(20 * p), 20, p)
  QP <- qp_w2(crossprod(X), stats::rnorm(p), 2)
  pen <- augmented_lagrangian_penalty(QP)
  TP <- add_cardinality_penalty(QP, pen$nu, pen$rho)
  
  captured <- NULL
  testthat::local_mocked_bindings(
    mosek = function(prob, opts) {
      captured <<- prob
      list(sol = list(int = list(xx = c(1, 1, rep(0, length(prob$c) - 2)))))
    },
    .package = "Rmosek")
  mosek_solver(TP)
  
  # binaries, then t (free) and |t| (>= 0) before mosek_qptoprob's extra columns
  testthat::expect_equal(captured$intsub, 1:p)
  testthat::expect_equal(unname(captured$bx[, p + 1]), c(-Inf, Inf))
  testthat::expect_equal(unname(captured$bx[, p + 2]), c(0, Inf))
  testthat::expect_equal(captured$c[p + 1:2], c(-pen$nu, pen$rho))
  # cardinality row has the slack and both |t| rows are present
  A <- as.matrix(captured$A)
  testthat::expect_equal(A[, p + 1][A[, 1] != 0 & A[, 2] != 0 & A[, p] != 0][1], 1)
  testthat::expect_equal(sum(A[, p + 2] != 0), 2)
})
