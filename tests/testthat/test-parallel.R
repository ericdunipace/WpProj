# Parallel computation is set up by the user with future::plan(). The
# deprecated `parallel` argument is still used as the plan until it is removed.

parallel_test_data <- function() {
  set.seed(2918)
  n <- 32
  p <- 5
  s <- 20
  x <- matrix(stats::rnorm(p * n), nrow = n, ncol = p)
  theta <- matrix((1:p)/p, nrow = p, ncol = s) + stats::rnorm(p * s, 0, 0.1)
  list(x = x, theta = theta)
}

# the workers only load WpProj if a loop was sent to them. Only check that some
# worker was used: future reuses a node once its future has resolved, so on a
# fast machine every chunk can land on the same worker.
workers_used <- function(cl) {
  any(unlist(parallel::clusterEvalQ(cl, "WpProj" %in% loadedNamespaces())))
}

test_that("loops run on the workers of the plan set with future::plan()", {
  testthat::skip_on_cran()
  dat <- parallel_test_data()
  seq_fit <- WpProj(dat$x, theta = dat$theta, method = "L0")

  cl <- parallel::makeCluster(2)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  oplan <- future::plan(future::cluster, workers = cl)
  on.exit(future::plan(oplan), add = TRUE, after = FALSE)

  par_fit <- WpProj(dat$x, theta = dat$theta, method = "L0")
  expect_true(workers_used(cl))
  expect_equal(par_fit$theta, seq_fit$theta)
})

test_that("seeded loops give the same results sequentially and in parallel", {
  testthat::skip_on_cran()
  dat <- parallel_test_data()
  opts <- simulated_annealing_method_options(nvars = 2:4, maxit = 5L, temps = 5L)

  set.seed(1)
  seq_fit <- WpProj(dat$x, theta = dat$theta, method = "simulated annealing",
                    options = opts)

  oplan <- future::plan(future::multisession, workers = 2)
  on.exit(future::plan(oplan), add = TRUE)
  set.seed(1)
  par_fit <- WpProj(dat$x, theta = dat$theta, method = "simulated annealing",
                    options = opts)
  expect_equal(par_fit$theta, seq_fit$theta)
})

test_that("the deprecated parallel argument still works and resets the plan", {
  testthat::skip_on_cran()
  dat <- parallel_test_data()
  seq_fit <- WpProj(dat$x, theta = dat$theta, method = "L0")
  seq_vi <- WPVI(X = dat$x, eta = dat$x %*% dat$theta, theta = dat$theta)

  cl <- parallel::makeCluster(2)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  plan_before <- future::plan()

  lifecycle::expect_deprecated(opts <- L0_method_options(parallel = cl))
  par_fit <- WpProj(dat$x, theta = dat$theta, method = "L0", options = opts)
  expect_true(workers_used(cl))
  expect_equal(par_fit$theta, seq_fit$theta)
  expect_identical(future::plan(), plan_before)

  lifecycle::expect_deprecated(
    par_vi <- WPVI(X = dat$x, eta = dat$x %*% dat$theta, theta = dat$theta,
                   parallel = cl)
  )
  expect_equal(par_vi, seq_vi)
  expect_identical(future::plan(), plan_before)
})
