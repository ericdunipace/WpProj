# WpProj 

## Version 0.2.9000 (development version)

### Breaking Changes
* Clarabel (`solver = "clarabel"`) replaces ECOS as the default solver for the `power = 1` and `power = Inf` L1 methods, and SCIP replaces ECOS as the default exact solver in the binary program (`W2IP`), where ECOS's branch and bound could return suboptimal subsets. ECOS remains available with `solver = "ecos"`. Both new solvers are suggested packages: if `ROI.plugin.clarabel` or `scip` is not installed (e.g., no Rust toolchain to build Clarabel from source), the defaults fall back to ECOS and lpSolve, respectively, with a once-per-session message.
* Moving from `doRNG` to `doFuture`. There had been some
issues with getting replies from the `doParallel` team so
concerned about ongoing support/updating. Parallel computation is now set up
with `future::plan()`, e.g. `future::plan(future::multisession, workers = 4)`,
and all functions use whatever plan is set.
* The `parallel` argument of `WPVI()`, `distCompare()`, and the
`*_method_options()` functions is deprecated in favor of `future::plan()`.
For now, a cluster from `parallel::makeCluster()` or a number of workers is
still accepted and used as the plan for the duration of the call.

### New Features
* Added free open-source solvers as alternatives to `mosek`: `solver = "clarabel"` for the `power = 1` and `power = Inf` L1 methods and `solver = "scip"` and `solver = "highs"` for the exact binary program method. The `"scip"` solver works on the binary quadratic program directly and supports a time limit. These require the suggested packages `ROI.plugin.clarabel`, `scip`, and `ROI.plugin.highs`.
* Added the augmented Lagrangian binary program of Gu, Ahmed, and Dey (2020) <doi:10.1137/19M1271695>, `binary_program_method_options(algorithm = "augmented.lagrangian")`, for the exact binary program method. It replaces the constraint on the number of coefficients with an exact augmented Lagrangian penalty from the continuous relaxation, giving the same solution as `algorithm = "exact"`, and works with every exact solver.

### Minor Improvements and Bug Fixes
* Making sure that final step of projection method will try to use an OLS regression rather than SVD iterations if possible
* Fixing some bugs (eg Eigen not doing explicit flattening anymroe, etc)
* Adding some more tests
* Merging PR from RcppEigen mainteners so they can update Eigen version!

## Version 0.2.3

### Minor Improvements and Bug Fixes
* Adding support for `approxOT` header files directly
* Fixing some bugs

## Version 0.2.1

### Minor Improvements and Bug Fixes
* Differences in sorting by system caused one test to fail on Fedora.
* Added URLs for Bug Reports
* Adding `lifecycle` badges to main functions


## Version 0.2

* Initial CRAN submission.
* Added a `NEWS.md` file to track changes to the package.
* Please note interface is currently a bit experimental and may change with little warning in the future
