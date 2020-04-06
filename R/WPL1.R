WPL1 <- function(X, Y=NULL, theta = NULL, power = 1,
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
                 gamma = 1, tau = 0.5, 
                 groups = numeric(0), 
                 scale.factor = numeric(0), 
                 penalty.factor = NULL, 
                 group.weights = NULL, maxit = 500L, 
                 tol = 1e-07,
                 display.progress=FALSE, ...) 
{
  this.call <- as.list(match.call()[-1])
  if(power == 2.0) {
    w2l1.arg.sel <- names(this.call) %in% formalArgs("W2L1")
    w2l1.args <- this.call[w2l1.arg.sel]
    output <- do.call("W2L1", w2l1.args)
  } else {
    if ("penalty" %in% names(this.call)) {
      penalty <- match.arg(penalty, several.ok = TRUE)
    }
    else {
      penalty <- match.arg(penalty, several.ok = FALSE)
    }
    if(!is.matrix(X)) X <- as.matrix(X)
    dims <- dim(X)
    
    p <- dims[2]
    power <- as.double(power)
    power <- ifelse(power < 1 | is.null(power) | missing(power), 1.0, power)
    
    if (inherits(X, "sparseMatrix")) {
      stop("Sparse matrices not allowed")
    }
    if (is.null(penalty.factor)) {
      penalty.factor <- rep(1, p)
    }
    varnames <- colnames(X)
    if (is.null(varnames))
      varnames = paste("V", seq(p), sep = "")
    if (length(penalty.factor) != p) {
      stop("penalty.factor must have same length as number of columns in x")
    }
    if ( any(penalty.factor < 0) ) {
      penalty.factor <- abs( penalty.factor )
      warning("penalty.factor had negative values. These were turned to positive values via abs().")
    }
    penalty.factor <- penalty.factor * p/sum(penalty.factor)
    penalty.factor <- drop(penalty.factor)
    if (any(grep("grp", penalty) > 0) ) {
      if (length(groups) != p) {
        stop("groups must have same length as number of columns in x")
      }
      unique.groups <- sort(unique(groups))
      zero.idx <- unique.groups[which(unique.groups == 0)]
      groups <- drop(groups)
      if (length(group.weights)!=0) {
        if (length(zero.idx) > 0) {
          group.weights[zero.idx] <- 0
        }
        group.weights <- drop(group.weights)
        if (length(group.weights) != length(unique.groups)) {
          stop("group.weights must have same length as the number of groups")
        }
        group.weights <- as.numeric(group.weights)
      } else {
        group.weights <- numeric(0)
      }
    } else {
      unique.groups <- numeric(0)
      group.weights <- numeric(0)
    }
    if (is.null(lambda.min.ratio)) {
      lambda.min.ratio <- 1e-04
    } else {
      lambda.min.ratio <- as.numeric(lambda.min.ratio)
    }
    
    if (lambda.min.ratio >= 1 | lambda.min.ratio <= 0) {
      stop("lambda.min.ratio must be between 0 and 1")
    }
    
    if (nlambda[1] <= 0) {
      stop("nlambda must be a positive integer")
    }
    
    if (!is.list(lambda)) {
      lambda <- sort(as.numeric(lambda), decreasing = TRUE)
      if (length(lambda) > 0) {
        lambda <- as.double(lambda)
      }
      lambda <- rep(list(lambda), length(penalty))
    } 
    else {
      if (length(lambda) != length(penalty)) {
        stop("If list of lambda vectors is provided, it must be the same length as the number of penalties fit")
      }
      nlambda.tmp <- length(lambda[[1]])
      for (l in 1:length(lambda)) {
        if (is.null(lambda[[l]]) || length(lambda[[l]]) <
            1) {
          stop("Provided lambda vector must have at least one value")
        }
        if (length(lambda[[l]]) != nlambda.tmp) {
          stop("All provided lambda vectors must have same length")
        }
        lambda[[l]] <- as.double(sort(as.numeric(lambda[[l]]),
                                      decreasing = TRUE))
      }
    }
    if (is.null(model.size)) {
      model.size <- p
    }
    
    #transpose X
    X_ <- t(X)
    
    if(is.null(Y)) {
      same <- TRUE
      Y_ <- crossprod(X_,theta_)
    } 
    else {
      if(!any(dim(Y) %in% dim(X_))) stop("Number of observations of Y must match X")
      if(!is.matrix(Y)) Y <- as.matrix(Y)
      if(nrow(Y) == ncol(X_)){ 
        # print("Transpose")
        Y_ <- Y
      } else{
        Y_ <- t(Y)
      }
    }
    if(nrow(Y_) != ncol(X_)) stop("The number of observations in Y and X don't line up. Make sure X is input with observations in rows.")
    nS <- dim(Y_)[2]
    
    #replicate groups to do multivariate regression
    if(length(groups) == 0 ) {
      
      groups <- 1:p
      groups[penalty.factor == 0] <- 0
      group.weights <- rep(penalty.factor, nS)
      groups <- rep(groups, nS)
      penalty.factor <- rep(1, p * nS)
      if(penalty != "ols") penalty <- paste0("grp.",penalty)
      
    } else {
      groups <- rep(groups, nS)
      group.weights <- rep(group.weights, nS)
    }
    
    #make R types align with c types
    groups <- as.integer(groups)
    unique.groups <- as.integer(unique.groups)
    nlambda <- as.integer(nlambda)
    alpha <- as.double(alpha)
    gamma <- as.double(gamma)
    tau <- as.double(tau)
    tol <- as.double(tol)
    maxit <- as.integer(maxit)
    display.progress <- as.logical(display.progress)
    model.size <- as.integer(model.size)
    
    if (length(scale.factor) > 0) {
      if (length(scale.factor) != p)
        stop("scale.factor must be same length as xty (nvars)")
      scale.factor <- as.double(scale.factor)
    }
    
    if (maxit <= 0) {
      stop("maxit should be greater than 0")
    }
    
    if (tol < 0) {
      stop("tol and irls.tol should be nonnegative")
    }
    options <- list(maxit = maxit, tol = tol,
                    display_progress = display.progress, 
                    model_size = model.size)
    # function (X_, Y_, theta_, power_, penalty_, groups_, 
    #           unique_groups_, group_weights_, lambda_, nlambda_, lmin_ratio_, 
    #           alpha_, gamma_, tau_, scale_factor_, penalty_factor_, opts_) 
    output <- WPpenalized(X_,Y_, power,
                          penalty, groups, unique.groups, group.weights, 
                          lambda, nlambda, lambda.min.ratio, alpha, gamma, tau, 
                          scale.factor, penalty.factor, options)
  
    if ( penalty != "ols" ) {
      options_ols <- options
      options_ols$display_progress <- FALSE
      penalty_ols <- "ols"
      options_ols$model_size <- p
      ols.out <- WPpenalized(X_,Y_, power,
                             penalty_ols, groups, unique.groups, group.weights, 
                             lambda, nlambda, lambda.min.ratio, alpha, gamma, tau, 
                             scale.factor, penalty.factor, options_ols)[c("beta","niter")]
        output$niter <- c(output$niter, 0)
        output$niter[ncol(output$niter)] <- ols.out$niter[1]
        output$lambda <- c(output$lambda, 0)
        output$beta <- cbind(output$beta, ols.out$beta)
        rm(ols.out)
    }
    
    output$nvars <- p
    output$power <- power
    output$penalty <- penalty
    output$method <- "projection"
    output$varnames <- varnames
    output$call <- formals(WPL1)
    output$call[names(this.call)] <- this.call
    output$nonzero_beta <- colSums(output$beta != 0)
    class(output) <- c("limbs", "optimization")
    extract <- extractTheta(output, matrix(output$beta[,ncol(output$beta)], p, nS))
    output$nzero <- extract$nzero
    output$eta <- lapply(extract$theta, function(tt) crossprod(X_, tt))
    output$theta <- extract$theta
  }
  
  return(output)
  
}
