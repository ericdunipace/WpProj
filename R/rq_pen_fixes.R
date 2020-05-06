l1.group.fit <- function (x, y, groups, lambda, intercept = TRUE, 
                          penalty = "MCP", alg = "QICD", a = 3.7, penGroups = NULL, 
                          ...) 
{
  tau <- 0.5
  p <- ncol(x)
  n <- nrow(x)
  if(is.null(alg) | missing(alg)) alg <- "QICD"
  if(is.null(penalty) | missing(penalty)) penalty <- "MCP"
  
  if (!penalty %in% c("SCAD", "MCP")) {
    stop("Penalty must be SCAD or MCP")
  }
  if (is.null(dim(x))) {
    stop("x must be matrix with at least 1 column")
  }
  if (length(groups) != ncol(x)) {
    stop("length(groups) must be equal to ncol(x)")
  }
  if (lambda <= 0) {
    stop("lambda must be positive")
  }
  if (penalty == "LASSO") {
    pen_func <- rqPen::lasso
  }
  if (penalty == "SCAD") {
    pen_func <- rqPen::scad
  }
  if (penalty == "MCP") {
    pen_func <- rqPen::mcp
  }
  if (alg == "QICD") {
    if (length(lambda) != 1) 
      stop("QICD Algorithm only allows 1 lambda value")
    coefs <- rqPen::QICD.group(y = y, x = x, groups = groups, 
                               lambda = lambda, intercept = intercept, 
                               penalty = penalty, a = a, ...)
    coefnames <- paste("x", 1:p, sep = "")
    if (intercept) 
      coefnames <- c("(Intercept)", coefnames)
    names(coefs) <- coefnames
    if (intercept) {
      residuals <- c(y - x %*% (coefs[-1]) - coefs[1])
      pen_vars <- coefs[-1]
    }
    else {
      residuals <- c(y - x %*% coefs)
      pen_vars <- coefs
    }
    if (penalty == "LASSO") {
      pen_val <- sum(pen_func(tapply(abs(pen_vars), groups, 
                                     sum), lambda = lambda))
    }
    else {
      pen_val <- sum(pen_func(tapply(abs(pen_vars), groups, 
                                     sum), lambda = lambda, a = a))
    }
    rho <- sum(rqPen::check(residuals))
    PenRho <- rho + pen_val
    return_val <- list(coefficients = coefs, PenRho = PenRho, 
                       residuals = residuals, rho = rho, tau = tau, n = n, 
                       intercept = intercept, penalty = penalty, alg = "QICD")
    class(return_val) <- c("rq.group.pen", "rq.pen")
  }
  else {
    group_num <- length(unique(groups))
    if (length(lambda) == 1) {
      lambda <- rep(lambda, group_num)
    }
    if (length(lambda) != group_num) {
      stop("lambdas do not match with group number")
    }
    if (sum(groups == 0) > 0) {
      stop("0 cannot be used as a group")
    }
    if (dim(x)[2] != length(groups)) {
      stop("length of groups must be equal to number of columns in x")
    }
    if (penalty == "LASSO") {
      stop("doesn't handle lasso penalties")
    }
    else {
      return_val <- rqPen::rq.group.lin.prog(x, y, groups, tau, 
                                             lambda, intercept = intercept, penalty = penalty, 
                                             penGroups = penGroups, a = a, ...)
      return_val$alg <- "lin.prog"
      class(return_val) <- c("rq.group.pen", "rq.pen")
    }
  }
  return_val
}


rqGroupLambda <- function(x, y, groups, lambda, intercept = FALSE, tau = 0.5,
                          penalty = "MCP", alg = "QICD_warm", penGroups = NULL, a = 3.7,
                          model.size = NULL,
                          ...) 
{
  return_val <- vector("list", length(lambda))
  pos <- 1
  tau <- 0.5
  intercept <- FALSE
  if(is.null(model.size) | length(model.size) == 0) model.size <- ncol(x)
  if(is.null(alg) | missing(alg)) alg <- "QICD_warm"
  alg <- match.arg(alg, choices = c("QICD_warm","lin.prog"))
  if(is.null(penalty) | missing(penalty)) penalty <- "MCP"
  
  
  if(penalty == "LASSO") {
    # dots <- list(...)
    # for(lam in lambda) {
    #   #X, Y, lambda, groups, solver, opts = NULL, init = NULL,...
    #   return_val[[pos]] <- l1_norm(X = x, Y = y, lambda = lam,
    #                                groups = groups, solver = dots$solver, 
    #                                opts = dots$opts, init = dots$init)
    #   pos <- pos + 1L
    # }
    # pen.lp <- switch(penalty,
    #                  "MCP" = "mcp",
    #                  "SCAD" = "scad",
    #                  "LASSO" = "lasso")
    temp_beta <- GroupLambda(X = x, Y = y, power = 1, groups = groups, lambda = lambda,
                             penalty = "lasso",
                             gamma = a, solver = list(...)$solver,
                             model.size = model.size, 
                             options = list(...)$options, 
                             ...)
    return_val <- lapply(1:ncol(temp_beta), function(tt) list(coefficients = temp_beta[,tt]))
    return(return_val)
  }
  if (alg != "QICD_warm") {
    pos <- 1
    for (lam in lambda) {
      return_val[[pos]] <- l1.group.fit(x = x, y = y, groups = groups,
                                        tau = tau, lambda = lam, intercept = intercept, a = a,
                                        penalty = penalty, alg = alg, penGroups = penGroups, 
                                        ...)
      if(sum(return_val[[pos]]$coefficients != 0) > model.size) {
        if(pos != length(lambda)) return_val[(pos):length(return_val)] <- NULL
        break
      }
      pos <- pos + 1
    }
  }
  else {
    p <- dim(x)[2]
    pos <- 1
    alg = "QICD"
    if (intercept) {
      initial_beta <- list(c(quantile(y, tau), rep(0, p)))
    }
    else {
      initial_beta <- list(rep(0, p))
    }
    for (lam in lambda) {
      return_val[[pos]] <- rqPen::QICD.group(y = y, x = x, groups = groups, tau = tau, 
                                             lambda = lam, intercept = intercept, initial_beta = initial_beta[[1]],
                                             penalty = "LASSO", ...)
      # return_val[[pos]] <- rqPen::rq.group.fit(x = x, y = y, groups = groups,
      #                                   tau = tau, lambda = lam, intercept = intercept,
      #                                   penalty = "LASSO", alg = alg, initial_beta = initial_beta[[1]],
      #                                   penGroups = penGroups, ...)
      initial_beta[[1]] <- return_val[[pos]]
      pos <- pos + 1
    }
    # temp_beta <- GroupLambda(X = x, Y = y, power = 1, groups = groups, lambda = lambda,
    #               penalty = "lasso",
    #               gamma = a, solver = list(...)$solver,
    #               options = list(...)$options, ...)
    # return_val <- lapply(1:ncol(temp_beta), function(tt) temp_beta[,tt])
    # if (penalty != "LASSO") {
    pos <- 1
    for (lam in lambda) {
      initial_beta[[1]] <- return_val[[pos]]
      return_val[[pos]] <- l1.group.fit(x = x, y = y,
                                        groups = groups, tau = tau, lambda = lam, intercept = intercept,
                                        a = a,
                                        penalty = penalty, alg = alg, initial_beta = initial_beta[[1]], 
                                        penGroups = penGroups, ...)
      if(sum(return_val[[pos]]$coefficients != 0) > model.size) {
        if(pos != length(lambda)) return_val[(pos):length(return_val)] <- NULL
        break
      }
      pos <- pos + 1
      
    }
    # }
  }
  return_val
}

# l1_norm <- function(X, Y, lambda, groups, solver, opts = NULL, init = NULL,...) {
#   
#   d <- ncol(X)
#   
#   group_length <- length(groups)
#   if(group_length > 0 ) {
#     group_idx <- lapply(unique(groups), function(i) which(groups == i))
#     ngroups <- length(group_idx)
#     group_length <- sapply(group_idx, length)
#   }
#   
#   if(is.null(init)) init <- rep(0, d)
#   
#   lambda_update <- list(rep(lambda, ngroups))
#   
#   problem <- lp_prob_w1(X = X, Y = Y, lambda = lambda_update[[1]], groups = groups)
#   
#   if (solver == "mosek") {
#     prob <-  list(sense = problem$sense,
#                   c = problem$C,
#                   A = problem$Const,
#                   bc = rbind(problem$Const_lower, problem$Const_upper),
#                   bx = rbind(problem$LB, problem$UB))
#     if(!is.null(init)) prob$sol <- list(bas = list(xx = init))
#     if(is.null(opts)) opts <- list(verbose = 0, )
#   } else if (solver == "gurobi") {
#     prob <-  list()
#   }
#   
#   beta <- lp_solve(prob, opts, solver)
#   return(beta)
# }
