WPL1 <- function(X, Y=NULL, theta = NULL, power = 2.0,
                 penalty =  c("lasso", "ols", "mcp", "elastic.net",
                              "selection.lasso",
                              "scad", "mcp.net",
                              "scad.net",
                              "grp.lasso",
                              "grp.lasso.net", "grp.mcp",
                              "grp.scad", "grp.mcp.net",
                              "grp.scad.net",
                              "sparse.grp.lasso"),
                 model.size = NULL,
                 lambda = numeric(0),
                 nlambda = 100L,
                 lambda.min.ratio = NULL, alpha = 1,
                 gamma = 1, maxit = 500L,
                 tol = 1e-07, ...)
{
  this.call <- as.list(match.call()[-1])
  stopifnot(power >= 1)
  if(power == 2.0) {
    pot.args <- c(as.list(environment()), list(...))
    w2l1.arg.sel <- which(names(pot.args) %in% formalArgs("W2L1"))
    w2l1.args <- pot.args[w2l1.arg.sel]
    argn <- lapply(names(w2l1.args), as.name)
    names(argn) <- names(w2l1.args)
    f.call <- as.call(c(list(as.name("W2L1")), argn))
    output <- eval(f.call, envir = w2l1.args)
    # output <- W2L1(X = X, Y = Y, theta = theta, penalty = penalty,
    #                )
  } else if (power == 1) {
    w1l1.args <- c(as.list(environment()), list(...))
    w1l1.args$alpha <- w1l1.args$tau <- NULL
    # w2l1.arg.sel <- which(names(pot.args) %in% formalArgs("W2L1"))
    # w2l1.args <- pot.args[w2l1.arg.sel]
    argn <- lapply(names(w1l1.args), as.name)
    names(argn) <- names(w1l1.args)
    f.call <- as.call(c(list(as.name("W1L1")), argn))
    output <- eval(f.call, envir = w1l1.args)
  } else if (power == Inf) {
    winfl1.args <- c(as.list(environment()), list(...))
    winfl1.args$alpha <- winfl1.args$tau <- NULL
    # w2l1.arg.sel <- which(names(pot.args) %in% formalArgs("W2L1"))
    # w2l1.args <- pot.args[w2l1.arg.sel]
    argn <- lapply(names(winfl1.args), as.name)
    names(argn) <- names(winfl1.args)
    f.call <- as.call(c(list(as.name("WInfL1")), argn))
    output <- eval(f.call, envir = winfl1.args)
  } else {
    output <- lp_reg(x = X, y = Y, theta = theta, power = power, gamma = gamma,
                     alpha = alpha,
                     penalty = penalty,  lambda = lambda,
                     nlambda = nlambda,
                     model.size = model.size,
                     iter = maxit,
                     tol = tol)
  }
  
  return(output)
  
}


lp_reg <- function(x, y, theta = NULL, power, gamma, alpha, tau, penalty, penalty.factor = numeric(0),
                    lambda = numeric(0),
                   nlambda = 100,
                   model.size = NULL,
                   iter = 100,
                   tol = 1e-7) {
  this.call <- as.list(match.call()[-1])
  
  log_sum_exp <- function(x) {
    # if(is.vector(x)) {
    if(all(is.infinite(x))) return(x[1])
    mx <- max(x)
    x_temp <- x - mx
    return(log(sum(exp(x_temp)))+ mx)
    # } else if (is.matrix(x)) {
    #   mx <- apply(x, 1, max)
    #   x_temp <- x - mx
    #   return(log(rowSums(exp(x_temp)))+ mx)
    # }
  }
  n <- nrow(x)
  d <- ncol(x)
  
  if(is.null(y) | missing(y)) {
    if(!(is.null(theta) | missing(theta))) {
      if(nrow(theta) != ncol(X)) theta <- t(theta)
      y <- x %*% theta
    }
  }
  
  s <- ncol(y)
  cols <- lapply(1:s, function(ss) Matrix::sparseMatrix(i = n*(ss-1) + rep(1:n,d),
                                                        j = rep(1:d,each = n),
                                                        x = c(x),
                                                        dims = c(n*s, d)))
  Xmat <- do.call(cbind, cols)
  rm(cols)
  
  if(is.null(model.size) | length(model.size) == 0) {
    model.size <- ncol(Xmat)
  } else {
    model.size <- model.size * s
  }
  
  Y <- as.matrix(c(y))
  beta <- beta_old <- lapply(1:length(lambda), function(l) rep(Inf, d * s))
  obs.weights <- list(rep(1/(n*s),n*s))
  ow <- list()
  Xw <- list()
  XtX <- list()
  XtY <- list()
  oem_holder <- list()
  
  if(penalty != "ols" & !grepl("grp.", penalty)) penalty <- paste0("grp.",penalty)
  
  if(length(penalty.factor) == 0) penalty.factor <- rep(1, d)
  groups <- rep(1:d, s)
  penalty.factor <- rep(penalty.factor, s)
  groups[penalty.factor == 0] <- 0
  
  if(length(lambda) == 0 | is.null(lambda)) {
    max.lambda <- sqrt(rowSums(crossprod(x,y)^2))
    lambda <- max.lambda * exp( seq(0, log(lambda.min.ratio), length.out = nlambda))
  }
  
  for(l in seq_along(lambda)) {
    lam <- lambda[lam]
    for(i in 1:iter) {
      ow[[1]]  <- Matrix::sparseMatrix(i = 1:(n*s), j = 1:(n*s), x = obs.weights[[1]],
                                      dims = c(n*s, n*s))
      Xw[[1]]  <- ow[[1]] * Xmat
      XtX[[1]] <- Matrix::crossprod( Xw[[1]], Xmat)
      XtY[[1]] <- Matrix::crosprod(Xw[[1]], Y)
      oem_holder[[1]] <- oem::oem.xtx(xtx = XtX[[1]], xty = XtY[[1]], family = "gaussian",
                                      penalty = penalty, lambda = lam,
                                      gamma = gamma, alpha = alpha,
                                      tol = tol, maxit = maxit, groups = groups,
                                      group.weights = penalty.factor)
      beta[[l]] <-  c(oem_holder[[1]]$beta)
      if(not.converged(beta[[l]], beta_old[[l]], 1e-7)){
        beta_old[[l]] <- beta[[l]]
        obs.weights[[1]] <- log(abs(Y - Xmat %*% beta[[1]])) * (power - 2)
        obs.weights[[1]] <- pmin(obs.weights[[1]], log(1e4))
        obs.weights[[1]] <- exp(obs.weights[[1]] - log_sum_exp(obs.weights[[1]]) )
      } else {
        break
      }
      if(sum(beta[[l]] != 0) >= model.size) break
      
    }
  }
  
  output <- list()
  output$beta <- do.call("cbind", beta)
  output$penalty <- penalty
  output$lambda <- lambda
  output$nvars <- d
  output$call <- formals(lp_reg)
  output$call[names(this.call)] <- this.call
  output$nonzero_beta <- colSums(output$beta != 0)
  output$method <- "projection"
  output$power <- power
  class(output) <- c("limbs", "optimization")
  
  extract <- extractTheta(output, matrix(0, d,s))
  output$nzero <- extract$nzero
  output$eta <- lapply(extract$theta, function(tt) X %*%tt )
  output$theta <- extract$theta
  output$model <- res
  
  return(output)
}
