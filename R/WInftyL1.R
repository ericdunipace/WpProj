WInfL1 <- function(X, Y, theta = NULL, penalty = c("none","lasso", "mcp","scad"), 
                 lambda = numeric(0), 
                 lambda.min.ratio = 1e-4, 
                 gamma = 1.5,
                 nlambda = 10, 
                 solver = c("mosek","gurobi"),
                 options = list(solver_opts = NULL,
                                init = NULL,
                                tol = 1e-7,
                                iter = 100),
                 model.size = NULL,
                 ...) {
  
  this.call <- as.list(match.call()[-1])
  
  solver <- match.arg(solver)
  
  if(penalty == "ols") penalty <- "none"
  penalty <- match.arg(penalty)
  
  n <- nrow(X)
  d <- ncol(X)
  
  if(is.null(Y) | missing(Y)) {
    if(!(is.null(theta) | missing(theta))) {
      if(nrow(theta) != ncol(X)) theta <- t(theta)
      Y <- X %*% theta
    }
  }
  
  s <- ncol(Y)
  cols <- lapply(1:s, function(ss) Matrix::sparseMatrix(i = n*(ss-1) + rep(1:n,d), 
                                                        j = rep(1:d,each = n), 
                                                        x = c(x),
                                                        dims = c(n*s, d)))
  Xmat <- do.call(cbind, cols)
  
  
  if(penalty != "none" & length(lambda) == 0) {
    if(!is.null(theta)) {
      lambda.max <- max(sqrt(rowSums(theta^2)))
      if(lambda.max == 0) lambda.max <- max(crossprod(x,Y))/(n)
    } else {
      lambda.max <- max(crossprod(x,Y))/(n)
    }
    lambda <-  exp(log(lambda.max) + seq(0, log(lambda.min.ratio), length.out = nlambda))
  } else if (penalty == "none"){
    lambda <- 0
    nlambda <- 0
  }
  
  if(is.null(options$tol)) options$tol <- 1e-7
  if(is.null(options$iter)) options$iter <- 100
  
  if(is.null(model.size) | length(model.size) == 0) {
    model.size <- ncol(Xmat)
  } else {
    model.size <- model.size * s
  }
  
  beta <- GroupLambda(X = Xmat, Y = Y, power = Inf, groups = rep(1:d,s), lambda = lambda,
                           penalty = penalty,
                     gamma = gamma, solver = solver,
                     model.size = model.size,
                     options = options, ...)
    # beta <- linf_norm(X = Xmat, Y = Y, deriv_func = deriv_func, thresholder = thresh_fun,
    #                  lambda = lambda, groups=rep(1:d, s), solver = solver, 
    #                  gamma = gamma, opts = options$solver_opts, init = options$init, iter = options$iter, tol = options$tol)
    
  if(solver == "mosek") Rmosek::mosek_clean()
  
  output <- list()
  output$beta <- as.matrix(beta)
  output$penalty <- penalty
  output$lambda <- lambda
  output$nvars <- p
  output$maxit <- NULL
  output$call <- formals(WInfL1)
  output$call[names(this.call)] <- this.call
  output$nonzero_beta <- colSums(output$beta != 0)
  output$method <- "projection"
  output$power <- Inf
  class(output) <- c("limbs", "optimization")
  
  extract <- extractTheta(output, matrix(0, d,s))
  output$nzero <- extract$nzero
  output$eta <- lapply(extract$theta, function(tt) X %*%tt )
  output$theta <- extract$theta
  return(output)
} 
