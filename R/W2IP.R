#' 2-Wasserstein distance selection by Integer Programming
#'
#' @param X Covariates
#' @param Y Predictions from arbitrary model
#' @param theta Parameters of original linear model. Required
#' @param transport.method Method for Wasserstein distance calculation. Should be one of the outputs of [transport_options()].
#' @param model.size Maximum number of coefficients in interpretable model
#' @param nvars The number of variables to explore. Should be an integer vector of model sizes. Default is NULL which will explore all models from 1 to `model.size`.
#' @param maxit Maximum number of solver iterations
#' @param infimum.maxit Maximum iterations to alternate binary program and Wasserstein distance calculation
#' @param tol Tolerance for convergence of coefficients
#' @param solver The solver to use. Must be one of "scip", "cone","lp", "highs", "cplex", "gurobi","mosek". 
#' @param algorithm How the cardinality constraint is handled. "exact" (default) enforces it directly. "augmented.lagrangian" is the augmented Lagrangian binary program of Gu, Ahmed, and Dey (2020) <doi:10.1137/19M1271695>: the constraint is softened with an augmented Lagrangian penalty from the continuous relaxation that is large enough to give the same solution. Works with any solver. See details.
#' @param display.progress Should progress be printed?
#' @param parallel `r lifecycle::badge("deprecated")` Use [future::plan()] to run the computations in parallel instead. A cluster from [parallel::makeCluster()] or a number of workers is still accepted for now and is used as the plan for the duration of the call.
#' @param ... Extra args to Wasserstein distance methods
#' 
#' @details
#' For argument `solver`, the default "scip" solves the binary quadratic program directly with the free SCIP solver (requires package `scip`; "lp" is used by default if it is not installed). Pass `control = list(time_limit = <seconds>)` to return the best solution found within a time limit. Options "cone" and "lp" use the free solvers "ECOS" and "lpSolver", respectively. Option "highs" uses the free HiGHS mixed-integer solver on the same linear reformulation as "lp" and requires package `ROI.plugin.highs`. "cplex", "gurobi" and "mosek" require installing the corresponding commercial solvers.
#' 
#' For `algorithm = "augmented.lagrangian"`, the constraint \eqn{\sum_j \alpha_j = k} is replaced by \eqn{\sum_j \alpha_j + t = k} with penalty \eqn{-\nu t + \rho |t|} in the objective, where \eqn{\nu} is the Lagrange multiplier of the continuous relaxation and \eqn{\rho} is the gap between a feasible solution and the relaxation's objective. Following Gu, Ahmed, and Dey (2020) <doi:10.1137/19M1271695>, a finite \eqn{\rho} closes the duality gap: by Lagrangian duality any \eqn{\alpha} with \eqn{t \neq 0} has penalized objective no better than that feasible solution, so the solution is the same as for `algorithm = "exact"`; how fast it is found depends on the solver. If the penalized problem ever returns a solution of the wrong size, the hard constrained problem is solved instead.
#' 
#' @keywords internal
# @examples
# if(rlang::is_installed("stats")) {
# n <- 128
# p <- 10
# s <- 100
# 
# x <- matrix( stats::rnorm( p * n ), nrow = n, ncol = p )
# x_ <- t(x)
# beta <- (1:p)/p
# y <- x %*% beta + stats::rnorm(n)
# post_beta <- matrix(beta, nrow=p, ncol=s) + stats::rnorm(p*s, 0, 0.1)
# post_mu <- x %*% post_beta
# 
# test <- W2IP(X = x, Y = post_mu, theta = post_beta, transport.method = "exact",
#              infimum.maxit = 10,
#              tol = 1e-7, solution.method = "cone",
#              display.progress = FALSE,nvars = c(2,4,8))
#              }
W2IP <- function(X, Y=NULL, theta,
                 transport.method = transport_options(),
                 model.size = NULL,
                 nvars = NULL,
                 maxit = 100L,
                 infimum.maxit = 100L,
                 tol = 1e-7,
                 solver = c("scip", "cone","lp", "highs", "mosek", "cplex", "gurobi"),
                 algorithm = c("exact", "augmented.lagrangian"),
                 display.progress=FALSE, parallel = NULL, ...) 
{
  this.call <- as.list(match.call()[-1])
  
  solution.method <- if (missing(solver)) NULL else solver
  
  dots <- list(...)
  if(!is.matrix(X)) X <- as.matrix(X)
  if(dim(X)[2] == 1) X <- t(X)
  if(!is.matrix(Y)) Y <- as.matrix(Y)  
  if(!is.matrix(theta)) theta <- as.matrix(theta)
  dims <- dim(X)
  p <- dims[2]
  varnames <- colnames(X)
  if (is.null(varnames))
    varnames = paste("V", seq(p), sep = "")
  infm.maxit <- infimum.maxit
  if(is.null(infm.maxit)){
    infm.maxit <- 100
  }
  
  if (is.null(model.size)) {
    model.size <- p
  }
  
  if(is.null(nvars)) nvars <- 1:model.size
  
  p_star <- length(nvars)
  
  if(is.null(transport.method)){
    transport.method <- "exact"
  } else {
    transport.method <- match.arg(transport.method, transport_options())
  }
  
  if(is.null(solution.method)) {
    solution.method <- resolve_default_solver("scip", internal = TRUE)
  } else {
    solution.method <- match.arg(solution.method, choices = c("scip", "cone","lp", "highs", "mosek", "cplex", "gurobi"))
  }
  algorithm <- match.arg(algorithm)
  
  
  # solves the binary program, optionally with the cardinality constraint
  # softened by `penalty` (see augmented_lagrangian_penalty()), and returns the binary vector
  solve_binary_program <- function(QP, control, solution.method, start, penalty = NULL) {
    soften <- function(op) {
      if (is.null(penalty)) return(op)
      add_cardinality_penalty(op, penalty$nu, penalty$rho)
    }
    TP <- switch(solution.method, 
                 # these reformulations need all binary variables so soften after.
                 # for a binary QP "socp" gives a linearization similar to bqp_to_lp
                 cone = soften(ROI::ROI_reformulate(QP, to = "socp")),
                 lp = soften(ROI::ROI_reformulate(QP, "lp", method = "bqp_to_lp")),
                 highs = soften(ROI::ROI_reformulate(QP, "lp", method = "bqp_to_lp")),
                 soften(QP))
    sol <- switch(solution.method, 
                  cone = ROI::ROI_solve(TP, solver = "ecos", control),
                  lp =  ROI::ROI_solve(TP, solver = "lpsolve", control),
                  highs = ROI::ROI_solve(TP, solver = "highs", control),
                  scip = scip_solver(TP, control),
                  cplex = ROI::ROI_solve(TP, solver = "cplex", control),
                  # gurobi = gurobi_solver(TP, control, start),
                  mosek = mosek_solver(TP, control, start)
    )
    switch(solution.method,
           "gurobi" = sol[1:p],
           "mosek" = sol[1:p],
           "scip" = sol[1:p],
           # HiGHS can return binaries off by ~1e-15
           "highs" = round(ROI::solution(sol)[1:p]),
           ROI::solution(sol)[1:p])
  }
  
  register_solver(solution.method) # registers ROI solver if needed
  
  # if ( solution.method %in% c("cone","lpsolve","cplex") ) ROI::ROI_require_solver(solution.method)
  
  if(ncol(theta) == ncol(X)){
    theta_ <- t(theta)
  } else {
    theta_ <- theta
  }
  if(nrow(theta_) != p) stop("dimensions of theta must match X")
  theta_save <- theta_
  
  #transpose X
  X_ <- t(X)
  
  same <- FALSE
  if(is.null(Y)) {
    same <- TRUE
    Y_ <- crossprod(X_,theta_)
  } else{
    if(!any(dim(Y) %in% dim(X_))) stop("dimensions of Y must match X")
    if(!is.matrix(Y)) Y <- as.matrix(Y)
    if(nrow(Y) == ncol(X_)){ 
      # print("Transpose")
      Y_ <- Y
    } else{
      Y_ <- t(Y)
    }
    if(all(Y_==crossprod(X_, theta_))) same <- TRUE
  }
  if(ncol(Y_) != ncol(theta_)) stop("ncol of Y should be same as ncols of theta")
  if(nrow(Y_) != ncol(X_)) stop("The number of observations in Y and X don't line up. Make sure X is input with observations in rows.")
  rmv.idx <- NULL
  if(any(apply(theta_,1, function(x) all(x == 0)))) {
    rmv.idx <- which(apply(theta_,1, function(x) all(x == 0)))
    
    X_ <- X_[-rmv.idx, ]
    theta_ <- theta_[-rmv.idx,]
    penalty.factor <- penalty.factor[-rmv.idx]
    warning("Some dimensions of theta have no variation. These have been removed")
  }
  
  # get control functions
  control <- dots$control
  if(is.null(control)) {
    control <- list()
  } 
  
  epsilon <- dots$epsilon
  if(is.null(epsilon)) epsilon <- 0.05
  OTmaxit <- dots$OTmaxit
  if(is.null(OTmaxit) || missing(OTmaxit)) OTmaxit <- switch(transport.method, "exact" = 0L, 100L)
  # else if (solution.method == "lp") {
  #   if(!is.null(control$verbose)) control$verbose <- as.logical(control$verbose)
  #   if(!is.null(control$presolve)) control$presolve <- as.logical(control$presolve)
  #   if(!is.null(control$tm_limit)) control$tm_limit <- as.integer(control$tm_limit)
  #   if(!is.null(control$canonicalize_status)) control$canonicalize_status <- as.logical(control$canonicalize_status)
  # } else if (solution.method == "cone") {
  #   
  #   control <- ecos.control.better(control)
  # }
  
  #make R types align with c types
  infm.maxit <- as.integer(infm.maxit)
  display.progress <- as.logical(display.progress)
  transport.method <- as.character(transport.method)
  nvars <- as.integer(nvars)
  
  if (infm.maxit <=0) {
    stop("infimum.maxit should be greater than 0")
  }
  
  oplan <- set_parallel_plan(parallel)
  if (!is.null(oplan)) {
    on.exit(future::plan(oplan), add = TRUE)
    display.progress <- FALSE
  }
  
  options <- list(infm_maxit = infm.maxit,
                  display_progress = display.progress, 
                  model_size = nvars)
  OToptions <- list(same = same,
                    method = "selection.variable",
                    transport.method = transport.method,
                    epsilon = epsilon,
                    niter = OTmaxit)
  
  ss <- sufficientStatistics(X, Y_, theta_, OToptions)
  xtx <- ss$XtX
  xty <- xty_init <- ss$XtY
  Ytemp <- Y_
  
  if(display.progress){
    pb <- utils::txtProgressBar(min = 0, max = p_star, style = 3)
    utils::setTxtProgressBar(pb, 0)
  }
  
  QP <- QP_orig <- qp_w2(ss$XtX,ss$XtY,1)
  # LP <- ROI::ROI_reformulate(QP,"lp",method = "bqp_to_lp" )
  alpha <- alpha_save <- rep(0,p)
  obj <- obj_save <- Inf
  beta <- matrix(0, nrow = p, ncol = p_star)
  iter.seq <- rep(0, p_star)
  comb <- function(x, ...) {
    # from https://stackoverflow.com/questions/19791609/saving-multiple-outputs-of-foreach-dopar-loop
    lapply(seq_along(x),
           function(i) c(x[[i]], lapply(list(...), function(y) y[[i]]))
    )
  }
  
  idx <- NULL
  
  output <- foreach::foreach(idx=1:p_star, .combine='comb', .multicombine=TRUE,
                             .init=list(list(), list()),
                             .errorhandling = 'pass', 
                             .inorder = FALSE,
                             .options.future = list(seed = TRUE)) %dofuture% 
    {
       m <- options$model_size[idx]
       QP <- QP_orig
       QP$constraints$rhs[1L] <- m
       results <- list(NULL, NULL)
       obj_save <- Inf
       for(inf in 1:options$infm_maxit) {
         penalty <- if (algorithm == "augmented.lagrangian") augmented_lagrangian_penalty(QP) else NULL
         # sol.meth <- if ( solution.method == "cone" && !("cone" %in% names(TP)) ) {
         #   "lp"
         # } else {
         #   solution.method
         # }
         # browser()
         # sol <- ROI::ROI_solve(LP, "glpk")
         # can use ROI.plugin.glpk:::.onLoad("ROI.plugin.glpk","ROI.plugin.glpk") to use base solver ^
         alpha <- solve_binary_program(QP, control, solution.method, start = alpha, penalty = penalty)
         if (!is.null(penalty) && (anyNA(alpha) || sum(round(alpha)) != m)) {
           alpha <- solve_binary_program(QP, control, solution.method, start = alpha)
         }
         obj <- c(0.5 * t(alpha) %*% (QP$objective$Q) %*% alpha - QP$objective$L %*% alpha)
         if(all(is.na(alpha))) {
           warning("Likely terminated early")
           break
         }
         if(not.converged(alpha, alpha_save, tol) || 
            not.converged(obj, obj_save, tol)){
           alpha_save <- alpha
           obj_save <- obj
           
           Ytemp <- selVarMeanGen(X_, theta_, as.double(alpha))
           xty   <- xtyUpdate(X, Ytemp, theta_, result_ = alpha, 
                                             OToptions)
           QP$objective$L$v <- c(-2*xty)
           
         } else {
           break
         }
       }
       if(display.progress) utils::setTxtProgressBar(pb, idx)
       results[[2]] <- inf
       results[[1]] <- alpha
       return(results)
       # iter.seq[idx] <- inf
       # if ( same ) QP <- QP_orig
       # QP$constraints$rhs[1] <- m + 1
       # beta[,idx] <- alpha
    }
  if (display.progress) close(pb)
  names(output) <- c("beta","niter")
  output$beta <- do.call("cbind", output$beta)
  output$niter <- unlist(output$niter)
  
  # if (!is.null(parallel) ){
  #   parallel::stopCluster(parallel)
  # }
  output[c("xtx", "xty_init","xty_final")] <- list(xtx, xty_init, xty)
  
  output$nvars <- p
  output$varnames <- varnames
  output$call <- formals(W2L1)
  output$call[names(this.call)] <- this.call
  output$remove.idx <- rmv.idx
  output$nonzero_beta <- colSums(output$beta != 0)
  # output$nzero <- nz
  class(output) <- c("WpProj","IP")
  extract <- extractTheta(output, theta_)
  output$nzero <- extract$nzero
  output$eta <- lapply(extract$theta, function(tt) crossprod(X_, tt))
  output$theta <- extract$theta
  if(!is.null(rmv.idx)) {
    for(i in seq_along(output$theta)){
      output$theta[[i]] <- theta_save
      output$theta[[i]][-rmv.idx,] <- extract$theta[[i]]
    }
  }
  
  return(output)
  
}

qp_w2 <- function(xtx, xty, K) {
  d <- NCOL(xtx)
  Q0 <- 2 * xtx # *2 since the Q part is 1/2 a^\top (x^\top x) a in ROI!!!
  L0 <- c(a = c(-2*xty))
  op <- ROI::OP(objective = ROI::Q_objective(Q = Q0, L = L0, names = as.character(1:d)),
                maximum = FALSE)
  ## sum(alpha) = K 
  A1 <- rep(1,d) 
  LC1 <- ROI::L_constraint(A1, ROI::eq(1), K)
  ROI::constraints(op) <- LC1
  ROI::types(op) <- rep.int("B", d)
  
  op$Upper <- chol(Q0)
  
  return(op)
}


# gurobi_solver <- function(problem,opts = NULL, start) {
#   
#   prob <-  list()
#   prob$Q <- Matrix::sparseMatrix(i=problem$objective$Q$i,
#                                  j = problem$objective$Q$j,
#                                  x = problem$objective$Q$v/2)
#   prob$modelsense <- 'min'
#   prob$obj <- as.numeric(problem$objective$L$v)
#   num_param <- length(problem$objective$L$v)
#   
#   prob$A <- Matrix::sparseMatrix(i=problem$constraints$L$i,
#                                  j = problem$constraints$L$j,
#                                  x = problem$constraints$L$v)
#   
#   prob$sense <- ifelse(problem$constraints$dir == "==", "=", NA)
#   # prob$sense <- rep(NA, length(qp$LC$dir))
#   # prob$sense[qp$LC$dir=="E"] <- '='
#   # prob$sense[qp$LC$dir=="L"] <- '<='
#   # prob$sense[qp$LC$dir=="G"] <- '>='
#   prob$rhs <- problem$constraints$rhs
#   prob$vtype <- rep("B", num_param)
#   prob$start <- start
#   
#   if(is.null(opts) | length(opts) == 0) {
#     opts <- list(OutputFlag = 0)
#   }
#   
#   res <- gurobi::gurobi(prob, opts)
#   
#   sol <- as.integer(res$x)
#   
#   return(sol)
# }

# Augmented Lagrangian binary program (Gu, Ahmed, and Dey 2020, SIAM J. Optim.):
# exact penalty for the cardinality constraint sum(a) = k:
# it is softened to sum(a) + t = k with penalty -nu * t + rho * |t|, where nu
# is the multiplier of the continuous relaxation. By duality, any a with
# sum(a) != k has penalized objective >= z_relax + rho, so setting rho to the
# gap between a feasible point and z_relax keeps the solution exact.
# Returns NULL if no penalty is needed or the relaxation can't be solved.
augmented_lagrangian_penalty <- function(problem) {
  Q <- as.matrix(problem$objective$Q)
  L <- as.numeric(as.matrix(problem$objective$L))
  p <- length(L)
  k <- problem$constraints$rhs[1L]
  if (k >= p) return(NULL)
  
  obj_fun <- function(a) c(0.5 * crossprod(a, Q %*% a)) + sum(L * a)
  
  # continuous relaxation min 0.5 a'Qa + L'a s.t. sum(a) = k via its KKT system
  kkt <- rbind(cbind(Q, 1), c(rep(1, p), 0))
  relax <- tryCatch(solve(kkt, c(-L, k)), error = function(e) NULL)
  if (is.null(relax)) return(NULL) # singular Q, keep the hard constraint
  
  a_relax <- relax[1:p]
  z_relax <- obj_fun(a_relax)
  
  # feasible point from the k largest relaxed values
  a_feas <- as.numeric(rank(-a_relax, ties.method = "first") <= k)
  rho    <- max(obj_fun(a_feas) - z_relax, 0)
  rho    <- rho * (1 + 1e-6) + 1e-8 # so infeasible points can't tie
  
  return(list(nu = relax[p + 1L], rho = rho))
}

# adds slack t (free) and |t| bound u to an ROI OP whose first constraint is
# the cardinality constraint
add_cardinality_penalty <- function(op, nu, rho) {
  A <- op$constraints$L
  m <- nrow(A)
  n <- ncol(A)
  
  A <- cbind(A, slam::simple_triplet_matrix(i = 1L, j = 1L, v = 1, nrow = m, ncol = 2L))
  A <- rbind(A, slam::simple_triplet_matrix(i = c(1L, 1L, 2L, 2L), j = n + c(1L, 2L, 1L, 2L),
                                            v = c(1, -1, 1, 1), nrow = 2L, ncol = n + 2L))
  
  L <- c(as.numeric(as.matrix(op$objective$L)), -nu, rho)
  Q <- op$objective$Q
  objective <- if (is.null(Q)) {
    ROI::L_objective(L)
  } else {
    Q <- slam::as.simple_triplet_matrix(Q)
    ROI::Q_objective(Q = slam::simple_triplet_matrix(i = Q$i, j = Q$j, v = Q$v,
                                                     nrow = n + 2L, ncol = n + 2L),
                     L = L)
  }
  
  vtypes <- ROI::types(op)
  if (is.null(vtypes)) vtypes <- rep("C", n)
  bnds <- op_bounds(op)
  
  out <- ROI::OP(objective = objective,
                 constraints = ROI::L_constraint(A, c(op$constraints$dir, "<=", ">="),
                                                 c(op$constraints$rhs, 0, 0)),
                 types = c(vtypes, "C", "C"),
                 bounds = ROI::V_bound(li = seq_len(n + 2L), ui = seq_len(n + 2L),
                                       lb = c(bnds$lb, -Inf, 0), ub = c(bnds$ub, Inf, Inf),
                                       nobj = n + 2L),
                 maximum = FALSE)
  if (!is.null(op$Upper)) out$Upper <- cbind(op$Upper, matrix(0, nrow(op$Upper), 2L))
  return(out)
}

# dense lower and upper variable bounds of an ROI OP
op_bounds <- function(op) {
  n  <- ncol(op$constraints$L)
  lb <- rep(0, n)
  ub <- rep(Inf, n)
  vtypes <- ROI::types(op)
  if (!is.null(vtypes)) ub[vtypes == "B"] <- 1
  b <- ROI::bounds(op)
  if (!is.null(b)) {
    if (length(b$lower$ind)) lb[b$lower$ind] <- b$lower$val
    if (length(b$upper$ind)) ub[b$upper$ind] <- b$upper$val
  }
  return(list(lb = lb, ub = ub))
}

scip_solver <- function(problem, opts = NULL) {
  # mixed binary QP min 0.5 a'Qa + L'a s.t. linear constraints, written as
  # min t s.t. 0.5 a'Qa + L'a - t <= 0 since SCIP needs a linear objective
  Q <- as.matrix(problem$objective$Q)
  L <- as.numeric(as.matrix(problem$objective$L))
  num_param <- length(L)
  
  ctrl <- if (inherits(opts, "scip_control")) {
    opts
  } else {
    do.call(scip::scip_control, utils::modifyList(list(verbose = FALSE), as.list(opts)))
  }
  
  model <- scip::scip_model("W2IP")
  on.exit(scip::scip_model_free(model))
  for (nm in names(ctrl$scip_params)) {
    scip::scip_set_param(model, nm, ctrl$scip_params[[nm]])
  }
  
  vtypes <- ROI::types(problem)
  if (is.null(vtypes)) vtypes <- rep("C", num_param)
  bnds <- op_bounds(problem)
  scip::scip_add_vars(model, obj = rep(0, num_param), lb = bnds$lb, ub = bnds$ub, vtype = vtypes)
  t_idx <- scip::scip_add_var(model, obj = 1, lb = -Inf, ub = Inf, vtype = "C")
  
  # a_i^2 = a_i for binaries so their diagonal of Q moves to the linear term
  bin  <- vtypes == "B"
  keep <- upper.tri(Q)
  diag(keep) <- !bin
  quad <- which(keep & Q != 0, arr.ind = TRUE)
  quadcoefs <- ifelse(quad[, 1] == quad[, 2], 0.5, 1) * Q[quad]
  scip::scip_add_quadratic_cons(model,
                                linvars = c(seq_len(num_param), t_idx),
                                lincoefs = c(L + 0.5 * diag(Q) * bin, -1),
                                quadvars1 = quad[, 1],
                                quadvars2 = quad[, 2],
                                quadcoefs = quadcoefs,
                                rhs = 0)
  
  A   <- as.matrix(problem$constraints$L)
  dir <- problem$constraints$dir
  rhs <- problem$constraints$rhs
  for (k in seq_len(nrow(A))) {
    nz <- which(A[k, ] != 0)
    scip::scip_add_linear_cons(model, vars = nz, coefs = A[k, nz],
                               lhs = if (dir[k] %in% c(">=", "==")) rhs[k] else -Inf,
                               rhs = if (dir[k] %in% c("<=", "==")) rhs[k] else Inf)
  }
  
  scip::scip_set_objective_sense(model, "minimize")
  scip::scip_optimize(model)
  
  if (scip::scip_get_nsols(model) == 0) return(rep(NA_real_, num_param))
  
  sol <- round(scip::scip_get_solution(model)$x[1:num_param])
  
  return(sol)
}

mosek_solver <- function(problem, opts = NULL, start) {
  
  cc        <- as.numeric(as.matrix(problem$objective$L))
  num_param <- length(cc)
  
  # Upper is chol(Q) so the objective is 0.5 ||Upper a||^2 + cc'a
  Upper <- problem$Upper
  if (ncol(Upper) < num_param) Upper <- cbind(Upper, matrix(0, nrow(Upper), num_param - ncol(Upper)))
  
  A   <- problem$constraints$L
  A   <- Matrix::sparseMatrix(i = A$i, j = A$j, x = A$v, dims = c(A$nrow, A$ncol))
  dir <- problem$constraints$dir
  rhs <- problem$constraints$rhs
  eq  <- dir == "=="
  ineq_sign <- ifelse(dir[!eq] == ">=", -1, 1) # mosek_qptoprob wants A a <= b
  
  bnds <- op_bounds(problem)
  prob <- mosek_qptoprob(F = Upper, f = cc, 
                         A = if (any(!eq)) ineq_sign * A[!eq, , drop = FALSE] else NA,
                         b = if (any(!eq)) ineq_sign * rhs[!eq] else NA,
                         Aeq = A[eq, , drop = FALSE],
                         beq = rhs[eq],
                         lb = bnds$lb,
                         ub = bnds$ub)
  
  vtypes <- ROI::types(problem)
  prob$intsub <- if (is.null(vtypes)) seq_len(num_param) else which(vtypes == "B")
  
  if(is.null(opts) | length(opts) == 0) opts <- list(verbose = 0)
  
  res <- Rmosek::mosek(prob, opts)
  
  sol <- round(res$sol$int$xx[1:num_param])
  
  return(sol)
}

mosek_qptoprob <- function (F = NA, f = NA, A = NA, b = NA, Aeq = NA, beq = NA, 
                            lb = NA, ub = NA) 
{
  # code directly from Rmosek::mosek_qptoprob
  stopifnot(all(!is.na(F)))
  stopifnot(all(!is.na(f)))
  stopifnot(length(f) == ncol(F))
  if (all(!is.na(A)) || all(!is.na(b))) {
    stopifnot(all(!is.na(A)) && all(!is.na(b)))
    if (!methods::is(A, "TsparseMatrix")) {
      A <- methods::as(A, "CsparseMatrix")
    }
    stopifnot(nrow(A) == length(b))
    stopifnot(ncol(A) == length(f))
  }
  else {
    A <- Matrix::Matrix(0, nrow = 0, ncol = length(f), sparse = TRUE)
    b <- numeric(0)
  }
  if (all(!is.na(Aeq)) || all(!is.na(beq))) {
    stopifnot(all(!is.na(Aeq)) && all(!is.na(beq)))
    if (!methods::is(Aeq, "TsparseMatrix")) {
      Aeq <- methods::as(Aeq, "CsparseMatrix")
    }
    stopifnot(nrow(Aeq) == length(beq))
    stopifnot(ncol(Aeq) == length(f))
  }
  else {
    Aeq <- Matrix::Matrix(0, nrow = 0, ncol = length(f), sparse = TRUE)
    beq <- numeric(0)
  }
  stopifnot(all(!is.na(lb)))
  stopifnot(length(lb) == length(f))
  stopifnot(all(!is.na(ub)))
  stopifnot(length(ub) == length(f))
  prob <- list(sense = "min")
  nt <- nrow(F)
  nx <- ncol(F)
  nrA <- nrow(A)
  nrEQ <- nrow(Aeq)
  prob$c <- c(f, 1, 0, rep(0, nt))
  prob$A <- rbind(cbind(A, Matrix::Matrix(0, nrA, 1), Matrix::Matrix(0, nrA, 
                                                     1), Matrix::Matrix(0, nrA, nt)), cbind(Aeq, Matrix::Matrix(0, nrEQ, 
                                                                                                1), Matrix::Matrix(0, nrEQ, 1), Matrix::Matrix(0, nrEQ, nt)), cbind(F, 
                                                                                                                                                    Matrix::Matrix(0, nt, 1), Matrix::Matrix(0, nt, 1), -1 * Matrix::Diagonal(nt)))
  prob$bc <- rbind(blc = c(rep(-Inf, nrA), beq, rep(0, nt)), 
                   buc = c(b, beq, rep(0, nt)))
  prob$bx <- rbind(blx = c(lb, 0, 1, rep(-Inf, nt)), bux = c(ub, 
                                                             Inf, 1, rep(Inf, nt)))
  prob$cones <- matrix(nrow = 2, dimnames = list(c("type", 
                                                   "sub"), c()), list("RQUAD", nx + (1:(2 + nt))), )
  return(prob)
}

