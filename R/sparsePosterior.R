SparsePosterior <- function(X, Y= NULL, theta, p = 2, 
                            estimation.method = c("L1", "SA","SW"), ...) {
  
  estimation.method <- match.arg(estimation.method)
  dots <- list(...)
  out <- NULL
  
  if ( estimation.method == "L1" ) {
    
    if ( p == 2) {
      penalty <- dots$penalty
      method <- dots$method
      transport.method <- dots$transport.method
      lambda <- dots$lambda
      nlambda <- dots$nlambda
      lambda.min.ratio = dots$lambda.min.ratio
      alpha <- dots$alpha 
      gamma <- dots$gamma 
      tau <- dots$tau
      groups <- dots$groups 
      scale.factor <- dots$scale.factor 
      penalty.factor <- dots$penalty.factor
      group.weights <- dots$group.weights
      maxit <- dots$maxit
      tol <- dots$tol
      irls.maxit <- dots$irls.maxit 
      irls.tol <- dots$irls.tol
      # pseudo_observations <- dots$pseudo_observations
      infimum.maxit <- dots$infimum.maxit
      display.progress <- dots$display.progress
      
      if(is.null(penalty)) penalty <- c("elastic.net", "selection.lasso",
                   "lasso", "ols", "mcp",
                   "scad", "mcp.net", 
                   "scad.net", 
                   "grp.lasso", 
                   "grp.lasso.net", "grp.mcp",
                   "grp.scad", "grp.mcp.net",
                   "grp.scad.net", 
                   "sparse.grp.lasso")
      if(is.null(method)){
        method = c("selection.variable","projection","location.scale","scale")
      }
      if(is.null(transport.method)){
        transport.method <- c("shortsimplex")
      }
      
      if(is.null(lambda)) lambda = numeric(0)
      if(is.null(nlambda)) nlambda = 100L
      if(is.null(lambda.min.ratio)) lambda.min.ratio = NULL
      if(is.null(alpha)) alpha = 1
      if(is.null(gamma)) gamma = 1
      if(is.null(tau)) tau = 0.5
      if(is.null(groups)) groups = numeric(0)
      if(is.null(scale.factor)) scale.factor = numeric(0)
      if(is.null(penalty.factor)) penalty.factor = NULL
      if(is.null(group.weights)) group.weights = NULL
      if(is.null(maxit)) maxit = 500L
      if(is.null(tol)) tol = 1e-07
      if(is.null(irls.maxit)) irls.maxit = 100L
      if(is.null(irls.tol)) irls.tol = 0.001
      # if(is.null(pseudo_observations)) pseudo_observations = 0.0
      if(is.null(infimum.maxit)) infimum.maxit=NULL
      if(is.null(display.progress)) display.progress=FALSE
      
      
      out <- W2L1(X = X, Y = Y, theta = theta, penalty = penalty, method = method,
                  transport.method = transport.method, lambda = lambda, nlambda = nlambda,
                  lambda.min.ratio = lambda.min.ratio, alpha = alpha, 
                  gamma = gamma, tau = tau, 
                  groups = groups, 
                  scale.factor = scale.factor, 
                  penalty.factor = penalty.factor, 
                  group.weights = group.weights, maxit = maxit, 
                  tol = tol, irls.maxit = irls.maxit, 
                  irls.tol = irls.tol,
                  # pseudo_observations = pseudo_observations,
                  infimum.maxit=infimum.maxit,
                  display.progress=display.progress
                  )
      
    } else {
      
      stop("Other values of 'p' not yet supported for estimation method 'L1'")
    }
    
    
  } else if (estimation.method == "SA") {
    
    force <- dots$force
    model.size <- dots$model.size
    iter <- dots$iter
    temps <- dots$temps
    # pseudo.obs <- dots$pseudo.obs
    proposal.fun <- dots$proposal
    SA_opt <- dots$options
    
    
    out <- WPSA(X = X, Y=Y, theta = theta, 
                force = force, p = p, model.size = model.size, 
                iter = iter, temps = temps, 
                # pseudo.obs = pseudo.obs, 
                proposal = proposal.fun, 
                options = SA_opt)
      
  } else if (estimation.method == "SW") {
    force <- dots$force
    direction <- dots$direction
    method <- dots$method
    
    out <- wPSW(X=X, Y=Y, theta=theta, force = force, p = p,
                direction = direction, 
                method = method)
    
  } else {
    stop("'estimation.method' not found")
  }
  
  return(out)
}
