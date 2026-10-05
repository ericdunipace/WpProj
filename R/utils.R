mround <- function(x, base) {
  base * round(x/base)
}

cround <- function(x, base) {
  base * ceiling(x/base)
}

fround <- function(x, base) {
  base * floor(x/base)
}

getDigits <- function(x) {
  if (abs(x) < 1) {
    digi <- nchar(strsplit(sub('0+$', '', as.character(x)), ".", fixed = TRUE)[[1]][[2]])
    quant <- x + 5/(10 ^(digi + 1))
  } else {
    digi <- (-(nchar(as.character(round(x)))-1))
    quant <- x + 5*10^abs(digi)
  }
  return(mround(quant, digi))
}

check_mosek <- function() {
  skip.fun <- !rlang::is_installed("Rmosek")
  
  if(skip.fun) {
    testthat::skip("Rmosek not found for tests")
  } else {
    mosek.err <- tryCatch(
      !is.character(Rmosek::mosek_version()),
      error = function(e) {TRUE}
    )
    if (mosek.err) testthat::skip("Rmosek installed but mosek optimizer not found for tests.")
  }
}

check_gurobi <- function() {
  skip.fun <- !rlang::is_installed("gurobi")
  if(skip.fun) {
    testthat::skip("gurobi not found for tests")
  }
}

check_highs <- function() {
  if(!rlang::is_installed("ROI.plugin.highs")) {
    testthat::skip("ROI.plugin.highs not found for tests")
  }
}

check_clarabel <- function() {
  if(!rlang::is_installed("ROI.plugin.clarabel")) {
    testthat::skip("ROI.plugin.clarabel not found for tests")
  }
}

check_scip <- function() {
  if(!rlang::is_installed("scip")) {
    testthat::skip("scip not found for tests")
  }
}

# Suggested solver packages. Defaults with a fallback are swapped for an
# always-installed solver when their package is missing (e.g., clarabel needs
# Rust to build from source). `fallback` is the name used by WpProj() and
# `fallback_internal` the name used by the underlying fitting functions.
optional_solvers <- list(
  clarabel = list(pkg = "ROI.plugin.clarabel", fallback = "ecos",    fallback_internal = "cone"),
  scip     = list(pkg = "scip",                fallback = "lpsolve", fallback_internal = "lp"),
  highs    = list(pkg = "ROI.plugin.highs",    fallback = NULL,      fallback_internal = NULL)
)

solver_installed <- function(pkg) rlang::is_installed(pkg)

# only for solvers the user did not explicitly request
resolve_default_solver <- function(solver, internal = FALSE) {
  opt <- optional_solvers[[solver]]
  if (is.null(opt) || is.null(opt$fallback) || solver_installed(opt$pkg)) return(solver)

  rlang::inform(sprintf("Package `%s` is not installed so using solver \"%s\" instead of the recommended \"%s\".",
                        opt$pkg, opt$fallback, solver),
                .frequency = "once", .frequency_id = paste0("WpProj_default_", solver))
  if (internal) opt$fallback_internal else opt$fallback
}

register_solver <- function(solution.method) {
  opt <- optional_solvers[[solution.method]]
  if (!is.null(opt)) {
    rlang::check_installed(opt$pkg, reason = sprintf("to use solver = \"%s\".", solution.method))
  }
  switch(solution.method,
         cone = ROI::ROI_require_solver("ecos"),
         lp =  ROI::ROI_require_solver("lpsolve"),
         cplex = ROI::ROI_require_solver("cplex"),
         clarabel = ROI::ROI_require_solver("clarabel"),
         highs = ROI::ROI_require_solver("highs")
  )
}

# solvers that use the ROI conic formulation built in lp_prob_to_model()
roi_cone_solvers <- function() c("cone", "clarabel")

is_inst <- function(pkg) {
  nzchar(find.package(pkg, quiet=TRUE))
}

# Set the future plan used by `%dofuture%` from a `parallel` argument: a cluster
# from parallel::makeCluster() or a number of workers. Returns the previous plan
# so the caller can restore it with on.exit(), or NULL if `parallel` is NULL, in
# which case the plan the user set with future::plan() is used.
set_parallel_plan <- function(parallel) {
  if (is.null(parallel)) return(NULL)
  if (inherits(parallel, "cluster")) {
    future::plan(future::cluster, workers = parallel)
  } else if (is.numeric(parallel)) {
    future::plan(future::multisession, workers = as.integer(parallel))
  } else {
    stop("parallel must be a cluster, the number of workers desired, or NULL")
  }
}

# Signal that the `parallel` argument of `what` is deprecated in favor of
# future::plan(). `user_env` is the environment of the code calling `what`, so
# lifecycle only warns when the argument is used directly.
deprecate_parallel <- function(what, user_env = rlang::caller_env(2)) {
  lifecycle::deprecate_soft(
    when = "0.3.0",
    what = paste0(what, "(parallel)"),
    details = "Set up parallel workers with `future::plan()` instead, e.g. `future::plan(future::multisession, workers = 4)`.",
    user_env = user_env
  )
}
