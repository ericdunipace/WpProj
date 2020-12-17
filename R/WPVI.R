#' p-Wasserstein variable importance
#'
#' @param X covariates
#' @param Y predictions
#' @param theta parameters from prediction
#' @param pred.fun prediction function. must take variables x, theta
#' @param p Power of Wasserstein power
#' @param ground_p Power of distance metric
#' @param transport.method Transport methods. One of "exact", "sinkhorn", "greenkhorn","randkhorn", "gandkhorn","hilbert"
#' @param epsilon Hyperparameter for sinkhorn iterations
#' @param OTmaxit maximum number of iterations for the Wasserstein method
#' @param display.progress Display intermediate progress
#' @param parallel foreach backend
#'
#' @return
#' @export
WPVI <- function(X, Y, theta, pred.fun = NULL, p = 2, ground_p = 2,
                 transport.method = transport_options(),
                 epsilon = 0.05,
                 OTmaxit = 100,
                 display.progress = FALSE,
                 parallel = NULL) {
  this.call <- as.list(match.call()[-1])
  
  d <- ncol(X)
  n <- nrow(X)
  
  if(is.null(pred.fun) ) {
    pred.fun <- function(x,theta) {
      return(x %*% theta)
    }
    stopifnot(is.matrix(theta))
    if(ncol(theta) == ncol(X)) {
      theta <- t(theta)
    }
  }
  stopifnot(is.function(pred.fun))
  
  S <- ncol(theta)
  X_ <- t(X)
  if(is.null(Y)) {
    Y_ <- pred.fun(X, theta)
    same <- TRUE
  } else {
    if(nrow(Y) != n){
      Y_ <- t(Y)
    } else {
      Y_ <- Y
    }
    same <- FALSE
    if(all(Y_==pred.fun(X, theta))) same <- TRUE
  }
  transport.method <- match.arg(transport.method)
  
  if(!is.null(parallel)){
    if(!inherits(parallel, "cluster")) {
      stop("parallel must be a registered cluster backend or NULL")
    }
    doParallel::registerDoParallel(parallel)
    display.progress <- FALSE
  } else{
    foreach::registerDoSEQ()
  }
  
  wp <- foreach::foreach(i  = 1:d) %dorng% {
    x_temp <- X
    x_temp[,i] <- 0
    eta <- pred.fun(x_temp, theta)
    return(
      WpProj::wasserstein(eta, Y_, 
                  p = p, ground_p = ground_p, 
                  observation.orientation = "colwise",
                transport.method = transport.method, 
                epsilon = epsilon, niter = OTmaxit)
    )
  }
  orders <- order(unlist(wp), decreasing = TRUE)
  names(orders) <- colnames(X)
  return(orders)
}
