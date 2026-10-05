
# Test for typical use cases
test_that("method_lookup works with typical inputs", {
  expect_equal(
    method_lookup(2.0, "L1", "lasso", NULL)$method, 
    "L1"
  )
  expect_equal(
    method_lookup(2.0, "L1", "lasso", NULL)$solver, 
    "lasso"
  )
  expect_equal(
    method_lookup(1.0, "L1", NULL, NULL)$solver, 
    resolve_default_solver("clarabel")
  )
  
  expect_equal(
    method_lookup(Inf, "L1", NULL, NULL)$solver, 
    resolve_default_solver("clarabel")
  )
  
  expect_equal(
    method_lookup(2.0, "binary program", NULL, NULL)$solver, 
    "lasso"
  )
  
})

# Test for boundary conditions
test_that("method_lookup handles boundary conditions correctly", {
  expect_error(method_lookup(0.5, "L1", "lasso", NULL))  # power < 1.0 should error
})

test_that("method_lookup handles different power values correctly", {
  # Power exactly at boundary (1.0 and 2.0)
  expect_equal(method_lookup(1.0, NULL, NULL, NULL)$power, 1.0)
  expect_equal(method_lookup(2.0, NULL, NULL, NULL)$power, 2.0)
  
  # Power just over boundary (slightly more than 1.0 and 2.0)
  expect_silent(result <- method_lookup(1.0001, NULL, NULL, NULL))
  expect_equal(result$power, 1.0001) # Adjusts to 3.0 in implementation
  expect_silent(result <- method_lookup(2.0001, NULL, NULL, NULL))
  expect_equal(result$power, 2.0001)
  expect_equal(result$fun, "WPL1")
  
  # Power just under boundary (slightly less than 1.0 and 2.0)
  expect_error(method_lookup(0.9999, NULL, NULL, NULL))
  expect_silent(result <- method_lookup(1.9999, NULL, NULL, NULL))
  expect_equal(result$power, 1.9999)
  
  # Infinite power
  expect_equal(method_lookup(Inf, NULL, NULL, NULL)$power, Inf)
})

# Test method and solver inputs for edge cases
test_that("method_lookup handles method and solver edge cases correctly", {
  # Valid method but invalid solver for a power
  expect_warning(result <- method_lookup(1.0, "L1", "nonexistent_solver", NULL))
  expect_equal(result$solver, resolve_default_solver("clarabel")) # Default to first valid solver
  
  # Invalid method and valid solver for a power
  expect_warning(result <- method_lookup(1.0, "nonexistent_method", "ecos", NULL))
  expect_equal(result$method, "L1") # Default to first valid method
  
  # Invalid method and solver
  expect_warning(result <- method_lookup(1.0, "nonexistent_method", "nonexistent_solver", NULL))
  expect_equal(result$method, "L1")
  expect_equal(result$solver, resolve_default_solver("clarabel"))
  
})

# Test for error handling
test_that("method_lookup handles invalid inputs appropriately", {
  expect_error(method_lookup("invalid", "L1", "lasso", NULL))  # Invalid power type
  expect_warning(method_lookup(2.0, "invalid_method", "lasso", NULL),
                 "Using method invalid_method with Wasserstein power 2 is not allowed. Switching to method L1.")  # Invalid method
  # Add more tests for other invalid inputs
})

# Test for typical use cases
test_that("method_lookup returns expected output for typical inputs", {
  # Test with typical valid inputs
  result <- method_lookup(1.0, "L1", "ecos", NULL)
  expect_equal(result$power, 1.0)
  expect_equal(result$method, "L1")
  expect_equal(result$solver, "ecos")
  # Add more assertions for different valid combinations of inputs
})

# Test for handling of power values
test_that("method_lookup handles different power values correctly", {
  # Power value exactly 1.0 or 2.0
  expect_equal(method_lookup(1.0, NULL, NULL, NULL)$power, 1.0)
  expect_equal(method_lookup(2.0, NULL, NULL, NULL)$power, 2.0)
  
  # Power value not 1.0, 2.0, or Inf (should default to 3.0 in this implementation)
  expect_silent(result <- method_lookup(2.5, NULL, NULL, NULL))
  expect_equal(result$power, 2.5)
  
  # Infinite power value
  expect_equal(method_lookup(Inf, NULL, NULL, NULL)$power, Inf)
  
  # Power less than 1.0 should throw an error
  expect_error(method_lookup(0.5, NULL, NULL, NULL))
})

# Test for method and solver defaults and error handling
test_that("method_lookup handles method and solver inputs correctly", {
  # Default method and solver for a given power
  result <- method_lookup(1.0, NULL, NULL, NULL)
  expect_equal(result$method, "L1")
  expect_equal(result$solver, resolve_default_solver("clarabel"))
  
  # Non-existent method for a power should give a warning and default to first method
  expect_warning(result <- method_lookup(1.0, "nonexistent", NULL, NULL))
  expect_equal(result$method, "L1")
  
  # Non-existent solver for a method should give a warning and default to first solver
  expect_warning(result <- method_lookup(1.0, "L1", "nonexistent", NULL))
  expect_equal(result$solver, resolve_default_solver("clarabel"))
  
  result <- method_lookup(3.0, NULL, NULL, NULL)
  expect_equal(result$method, "L1")
  expect_equal(result$solver, "lasso")
  
  testthat::expect_warning(result <- method_lookup(3.0, "L1", "ecos", NULL))
  expect_equal(result$method, "L1")
  expect_equal(result$solver, "lasso")
})


# Test for typical use cases of WpProj
testthat::test_that("WpProj works with typical inputs", {
  set.seed(84370158)
  
  n <- 32
  p <- 10
  s <- 99
  
  x <- matrix( stats::rnorm( p * n ), nrow = n, ncol = p )
  beta <- (1:p)/p
  y <- x %*% beta + stats::rnorm(n)
  post_beta <- matrix(beta, nrow=p, ncol=s) + stats::rnorm(p*s, 0, 0.1)
  post_mu <- x %*% post_beta
  
  slotnames <- c("call",
                 "theta",
                 "fitted.values",
                 "power",
                 "method",
                 "solver",
                 "niter",
                 "nzero")
  
  # Test with typical valid inputs
  result <- WpProj(X = x, eta = post_mu, theta = post_beta, power = 2.0, method = "L1",
                   solver = "lasso")
  
  # Test that output is of class WpProj
  testthat::expect_s3_class(result, "WpProj")
  testthat::expect_equal(result$power, 2.0)
  testthat::expect_equal(names(result),
                         slotnames)
  
  # Test with typical valid inputs
  testthat::expect_warning(result2 <- WpProj(X = x, eta = post_mu, theta = post_beta, 
                                             power = 2.0, method = "L1", solver = "ecos"))
  
  testthat::expect_equal(result[-1], result2[-1])
  
  # Test that output is of class WpProj
  testthat::expect_s3_class(result2, "WpProj")
  testthat::expect_equal(result2$power, 2.0)
  testthat::expect_equal(names(result2),
                         slotnames)
  
  # Test L1, W1
  result <- WpProj(X = x, eta = post_mu, theta = post_beta, power = 1.0, method = "L1",
                   solver = "ecos", options = list(nlambda = 2L))
  testthat::expect_s3_class(result, "WpProj")
  testthat::expect_equal(result$power, 1.0)
  testthat::expect_equal(names(result),
                         slotnames)
  
  # test L1, WInf
  result <- WpProj(X = x, eta = post_mu, theta = post_beta, power = Inf,
                   options = list(nlambda = 2L))
  testthat::expect_s3_class(result, "WpProj")
  testthat::expect_equal(result$power, Inf)
  testthat::expect_equal(names(result),
                         slotnames)
  
  # test binary program, W2
  result <- WpProj(X = x, eta = post_mu, theta = post_beta, power = 2.0, method = "binary program",
                   solver = "ecos")
  testthat::expect_s3_class(result, "WpProj")
  testthat::expect_equal(result$power, 2.0)
  testthat::expect_equal(names(result),
                         slotnames)
  
  testthat::expect_error(WpProj(X = x, eta = post_mu, power = 2.0, method = "binary program",                    solver = "ecos"))
  
  result <- WpProj(X = x, eta = post_mu, theta = post_beta, power = 2.0, method = "binary program",
                   solver = "lasso")
  
  # test simulated annealing, W2
  result <- WpProj(X = x, eta = post_mu, theta = post_beta, power = 2.0, method = "simulated annealing")
  testthat::expect_s3_class(result, "WpProj")
  testthat::expect_equal(result$power, 2.0)
  testthat::expect_equal(names(result),
                         slotnames)
  
  # test stepwise, W2
  result <- WpProj(X = x, eta = post_mu, theta = post_beta, power = 2.0, method = "stepwise")
  testthat::expect_s3_class(result, "WpProj")
  testthat::expect_equal(result$power, 2.0)
  testthat::expect_equal(names(result),
                         slotnames)
  
  
  # test L0, W2
  result <- WpProj(X = x, eta = post_mu, theta = post_beta, power = 2.0, method = "L0")
  testthat::expect_s3_class(result, "WpProj")
  testthat::expect_equal(result$power, 2.0)
  testthat::expect_equal(names(result),
                         slotnames)
  
  
})



test_that("method_lookup routes open-source solvers", {
  expect_silent(result <- method_lookup(1.0, "L1", "clarabel", NULL))
  expect_equal(result$solver, "clarabel")
  expect_equal(result$fun, "W1L1")
  
  expect_silent(result <- method_lookup(Inf, "L1", "clarabel", NULL))
  expect_equal(result$solver, "clarabel")
  expect_equal(result$fun, "WInfL1")
  
  expect_silent(result <- method_lookup(2.0, "binary program", "highs", NULL))
  expect_equal(result$solver, "highs")
  expect_equal(result$fun, "W2IP")
  
  expect_silent(result <- method_lookup(2.0, "binary program", "scip", NULL))
  expect_equal(result$solver, "scip")
  expect_equal(result$fun, "W2IP")
  
  # clarabel can't solve binary programs and highs can't solve the cone programs
  expect_warning(result <- method_lookup(2.0, "binary program", "clarabel", NULL))
  expect_equal(result$solver, "lasso")
  expect_warning(result <- method_lookup(1.0, "L1", "highs", NULL))
  expect_equal(result$solver, resolve_default_solver("clarabel"))
})

test_that("WpProj runs with open-source solvers", {
  set.seed(87897)
  n <- 32
  p <- 10
  s <- 21
  x <- matrix(stats::rnorm(p*n), nrow=n, ncol=p)
  post_beta <- matrix((1:p)/p, nrow=p, ncol=s) + stats::rnorm(p*s, 0, 0.1)
  post_mu <- x %*% post_beta
  
  check_clarabel()
  fit_ecos <- WpProj(X=x, eta=post_mu, power = 1.0, method = "L1", solver = "ecos",
                     options = list(nlambda = 5))
  fit_clarabel <- WpProj(X=x, eta=post_mu, power = 1.0, method = "L1", solver = "clarabel",
                         options = list(nlambda = 5))
  expect_equal(fit_clarabel$solver, "clarabel")
  expect_equal(fit_clarabel$nzero, fit_ecos$nzero)
  
  check_highs()
  fit_highs <- WpProj(X=x, eta=post_mu, theta = post_beta, power = 2.0, 
                      method = "binary program", solver = "highs",
                      options = list(nvars = c(2, 4)))
  expect_equal(fit_highs$solver, "highs")
  expect_equal(fit_highs$nzero, c(2, 4))
  
  check_scip()
  fit_scip <- WpProj(X=x, eta=post_mu, theta = post_beta, power = 2.0, 
                     method = "binary program", solver = "scip",
                     options = list(nvars = c(2, 4),
                                    solver.options = list(control = list(time_limit = 30))))
  expect_equal(fit_scip$solver, "scip")
  expect_equal(fit_scip$nzero, c(2, 4))
})

test_that("default solvers fall back when optional solver packages are missing", {
  local_mocked_bindings(solver_installed = function(pkg) FALSE)
  
  expect_equal(suppressMessages(method_lookup(1.0, "L1", NULL, NULL))$solver, "ecos")
  expect_equal(suppressMessages(method_lookup(Inf, "L1", NULL, NULL))$solver, "ecos")
  expect_equal(suppressMessages(resolve_default_solver("clarabel", internal = TRUE)), "cone")
  expect_equal(suppressMessages(resolve_default_solver("scip", internal = TRUE)), "lp")
  expect_equal(resolve_default_solver("highs"), "highs")
  
  # explicitly requested solvers are not swapped
  expect_equal(method_lookup(1.0, "L1", "clarabel", NULL)$solver, "clarabel")
  
  # the fitting functions fall back too
  set.seed(87897)
  n <- 32
  p <- 5
  s <- 21
  x <- matrix(stats::rnorm(p*n), nrow=n, ncol=p)
  post_beta <- matrix((1:p)/p, nrow=p, ncol=s) + stats::rnorm(p*s, 0, 0.1)
  post_mu <- x %*% post_beta
  
  fit_default <- suppressMessages(W2IP(X = x, Y = post_mu, theta = post_beta, nvars = c(2, 3)))
  fit_lp <- W2IP(X = x, Y = post_mu, theta = post_beta, nvars = c(2, 3), solver = "lp")
  expect_equal(fit_default$beta, fit_lp$beta)
  
  fit_default <- suppressMessages(W1L1(X = x, Y = post_mu, penalty = "lasso", nlambda = 3))
  fit_ecos <- W1L1(X = x, Y = post_mu, penalty = "lasso", nlambda = 3, solver = "cone")
  expect_equal(fit_default$beta, fit_ecos$beta)
})

test_that("WpProj passes the augmented lagrangian algorithm to the binary program", {
  check_scip()
  set.seed(87897)
  n <- 32
  p <- 10
  s <- 21
  x <- matrix(stats::rnorm(p*n), nrow=n, ncol=p)
  post_beta <- matrix((1:p)/p, nrow=p, ncol=s) + stats::rnorm(p*s, 0, 0.1)
  post_mu <- x %*% post_beta
  
  fit_exact <- WpProj(X=x, eta=post_mu, theta = post_beta, power = 2.0, 
                      method = "binary program", solver = "scip",
                      options = list(nvars = c(2, 4)))
  fit_lagr <- WpProj(X=x, eta=post_mu, theta = post_beta, power = 2.0, 
                     method = "binary program", solver = "scip",
                     options = list(nvars = c(2, 4), algorithm = "augmented.lagrangian"))
  expect_equal(fit_lagr$theta, fit_exact$theta)
  
  # ignored by the approximate lasso method
  expect_silent(WpProj(X=x, eta=post_mu, theta = post_beta, power = 2.0, 
                       method = "binary program", solver = "lasso",
                       options = list(algorithm = "augmented.lagrangian")))
})
