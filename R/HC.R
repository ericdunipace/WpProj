HC <- function(X, Y=NULL, theta, family="gaussian", 
               penalty =  c("elastic.net", "selection.lasso",
                            "lasso", "ols", "mcp",
                            "scad", "mcp.net", 
                            "scad.net", 
                            "grp.lasso", 
                            "grp.lasso.net", "grp.mcp",
                            "grp.scad", "grp.mcp.net",
                            "grp.scad.net", 
                            "sparse.grp.lasso"), 
               method = c("selection.variable", "projection"),
               lambda = numeric(0), 
               nlambda = 100L, 
               lambda.min.ratio = NULL, alpha = 1,  
               gamma = 1, tau = 0.5, 
               groups = numeric(0), 
               scale.factor = numeric(0), 
               penalty.factor = NULL, 
               group.weights = NULL, maxit = 500L, 
               tol = 1e-07, irls.maxit = 100L, 
               irls.tol = 0.001,
               intercept = TRUE)
{
  if( family == "cox" | family == "exponential" | family == "survival") {
    if(!(penalty == "lasso")) {
      warning("penalty must be lasso for survival method")
      penalty <- "lasso"
    }
  }
  method <- match.arg(method)
  
  p <- ncol(X)
  n <- nrow(X)
  
  if(ncol(X) == ncol(theta)) {
    theta <- t(theta)
  }
  s <- ncol(theta)
  
  if(is.null(Y)) {
    Y <- X %*% theta
    if ( family == "binomial") {
      Y_ <- c(rep(1,n), rep(0,n))
    } else if (family == "gaussian" ) {
      Y_ <- rowMeans(Y)
    } else if (family == "exponential" | family == "survival" | family == "cox") {
      Y_ <- NULL # need to implement survival times
    }
  } else {
    stopifnot(nrow(Y) == nrow(X))
    stopifnot(ncol(Y) == ncol(theta))
      Y_ <- rowMeans(Y)
  }
  if(family == "binomial") {
    prob <- rowMeans(plogis(X %*% theta))
    weights <- c(prob, 1-prob)
  } else {
    weights <- numeric(0)
  }
  
  hessian.type <- "upper.bound"
  if(family == "binomial"){
    if(n > 100*p){
      hessian.type <- "upper.bound"
    } else {
      hessian.type <- "full"
    }
  }
  if(all(X[,1]==1) & intercept == TRUE) {
    x <- X[,-1]
    if(length(penalty.factor) == ncol(X)) penalty.factor <- penalty.factor[-1]
  } else {
    x <- X
  }
  
  if( family %in% c("gaussian")) {
    
   output <-  oem::oem(x, Y_, family = family, penalty = penalty,
                  weights = weights, lambda = lambda, nlambda = nlambda,
                  lambda.min.ratio = lambda.min.ratio, alpha=alpha, gamma = gamma, tau = tau,
                  groups = groups, penalty.factor = penalty.factor, group.weights = group.weights,
                  standardize = TRUE, intercept = intercept, 
                  maxit = maxit, tol = tol, irls.maxit = irls.maxit,
                  irls.tol = irls.tol, accelerate = FALSE, ncores = -1, 
                  compute.loss = FALSE, hessian.type = hessian.type)
  
   } else if (family == "exponential" | family == "survival" | family == "cox") {
    
    output <-  glmnet::glmnet(x, Y_, family = "cox", penalty = penalty,
             weights = weights, lambda = lambda, nlambda = nlambda,
             lambda.min.ratio = lambda.min.ratio, alpha=alpha, 
             thresh = 1e-07, dfmax = ncol(X) + 1,
             penalty.factor = penalty.factor,
             standardize.response = TRUE, intercept = intercept)
   } else if (family == "binomial") {
     
     output <-  glmnet::glmnet(x, Y_, family = "binomial", penalty = penalty,
             weights = weights, lambda = lambda, nlambda = nlambda,
             lambda.min.ratio = lambda.min.ratio, alpha=alpha, 
             thresh = 1e-07, dfmax = ncol(X) + 1,
             penalty.factor = penalty.factor,
             standardize.response = TRUE, intercept = intercept)
  }
  if(!intercept) {
    output$beta[[1]] <- output$beta[[1]][-1,]
  }
  extract <- extractCoef(output)
  output$nzero <- c(extract$nzero,p)
  output$E_theta <- extract$coefs
  output$E_eta <- lapply( 1:ncol(output$E_theta), function(et) X %*% output$E_theta[,et] )
  xtx <- crossprod(X)
  xty <- crossprod(X,Y)
  output$theta <- lapply(1:ncol(output$E_theta), function(i){
    et <- output$E_theta[,i]
    coefs <- matrix(0, p, s) 
    idx <- which(et != 0)
    if(length(idx) == 0) return(coefs)
    if ( method == "selection.variable" ) {
      coefs[idx,] <- theta[idx,]
    } else {
      coefs <- calc.beta(xtx = xtx, xty = xty, 
                         active.idx = idx, method = "projection",
                          OToptions = NULL)
    }
    return(coefs)
  })
  output$theta[[length(output$theta)+1]] <- theta
  output$eta <- lapply(output$theta, function(tt) X %*% tt)
  class(output) <- c("limbs", "HC")
  
  return(output)
}