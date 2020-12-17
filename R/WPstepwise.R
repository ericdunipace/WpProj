#' p-Wasserstein distance linear projections using a stepwise method
#'
#' @param X matrix of covariates
#' @param Y matrix of predictions
#' @param theta optional parameter matrix for selection methods.
#' @param force Any covariates to force into the model?
#' @param p Power of the wasserstein distance
#' @param ground_p Distance metric power. Typically the same as `p`
#' @param direction forward or backward selection
#' @param method "selection.variable" or "projection
#' @param transport.method Method for calculating the wasserstein distance. One of "exact", "sinkhorn", "greenkhorn","randkhorn", "gandkhorn","hilbert"
#' @param epsilon hyperparameter if using sinkhorn iterations to approximate OT
#' @param OTmaxit maximum number of iterations for the opt?imal transport methods
#' @param calc.theta should we get the linear coefficients
#' @param model.size Maximum model size
#' @param parallel foreach backend
#' @param display.progress Display intermediate progress
#'
#' @return
#' @export
WPSW <- function(X, Y, theta, force = NULL, p = 2, ground_p = 2,
                        direction = c("backward","forward"), 
                        method=c("selection.variable","scale","projection"),
                        transport.method = transport_options(),
                        epsilon = 0.05,
                        OTmaxit = 100,
                        calc.theta = TRUE,
                        model.size = NULL,
                        parallel = NULL,
                        display.progress = FALSE) {
  this.call <- as.list(match.call()[-1])
  
  d <- ncol(X)
  n <- nrow(X)
  if(ncol(theta) == ncol(X)) {
    theta <- t(theta)
  } 
  
  S <- ncol(theta)
  X_ <- t(X)
  if(is.null(Y)) {
    Y_ <- crossprod(X_,theta)
    same <- TRUE
  } else {
    if(nrow(Y) != n){
      Y_ <- t(Y)
    } else {
      Y_ <- Y
    }
    same <- FALSE
    if(all(Y_==crossprod(X_, theta))) same <- TRUE
  }
  method <- match.arg(method)
  transport.method <- match.arg(transport.method, transport_options())
  if(!is.null(force)) stopifnot(is.numeric(force))
  if(is.null(epsilon)) epsilon <- 0.05
  if(is.null(OTmaxit)) OTmaxit <- 100
  
  if(!is.null(parallel)){
    if(!inherits(parallel, "cluster")) {
      stop("parallel must be a registered cluster backend")
    }
    doParallel::registerDoParallel(parallel)
    # display.progress <- FALSE
  } else{
    foreach::registerDoSEQ()
  }
  # stopifnot(is.character(diretction))
  
  # if (grepl("univariate", transport.method) ) {
  #   obs.direction <- "rowwise"
  #   # Y_ <- apply(Y,1,sort)
  # } else {
  #   obs.direction <- "colwise"
  #   # Y_ <- Y
  # }
  
  direction <- match.arg(direction)
  not.force.logical <- !(1:d %in% force)
  l_force <- length(force)
  max_iter <- (d - max(l_force, 1))
  sel.idx <- rep(NA, max_iter)
  wP_traj <- rep(NA, max_iter+1)
 
  wP_traj[max_iter + 1] <- 0
  
  xtx <- xty <- NULL
  
  OToptions <- list(
    same = FALSE,
    method = method,
    transport.method = transport.method,
    epsilon = as.double(epsilon),
    niter =  as.integer(OTmaxit)
  )
  
  # theta_norm <- colMeans(theta^2)
  # wt <- n /(n + pseudo.obs)
  
  if(method == "selection.variable") {
    add.idx <- function(j, in.idx = NULL, X = NULL, sort_mu = NULL, p, ground_p,
                        OToptions, obs.direction, ...) {
      idx <- c(which(in.idx),j)
      temp_mu <- crossprod(X[idx,, drop=FALSE], theta[idx,, drop=FALSE])
      wp <- WpProj::wasserstein(X = sort_mu, Y = temp_mu, 
                                         p = p, ground_p = ground_p, 
                                         observation.orientation = obs.direction, 
                                         method = OToptions$transport.method, 
                                         epsilon = OToptions$epsilon, niter = OToptions$niter)
      beta <- rep(0, nrow(X))
      beta[idx] <- 1
      return(list(wp = wp, beta = beta))
    }
    minus.idx <- function(j, in.idx = NULL, X = NULL, sort_mu = NULL, p, ground_p,
                          OToptions, obs.direction,...) {
      temp.in.idx <- in.idx
      temp.in.idx[ j ] <- FALSE
      idx <- which( temp.in.idx )
      temp_mu <- crossprod(X[idx,, drop=FALSE], theta[idx,, drop=FALSE])
      wp <- WpProj::wasserstein(X = sort_mu, Y = temp_mu, 
                                         p = p, ground_p = ground_p, 
                                         observation.orientation = obs.direction, 
                                         method = OToptions$transport.method, 
                                         epsilon = OToptions$epsilon, niter = OToptions$niter)
      beta <- rep(0, nrow(X))
      beta[idx] <- 1
      return(list(wp = wp, beta = beta))
    }
  } else {
    add.idx <- function(j, in.idx = NULL, X = NULL, sort_mu = NULL, p, ground_p, OToptions, obs.direction,...) {
      idx <- c(which(in.idx),j)
      beta <- calc.beta(xtx, xty, idx, method, OToptions, X, theta)
      d <- length(idx)
      if(method != "projection"){
        # beta <- theta %*% diag(beta)
        temp_mu <- selVarMeanGen(X, theta, beta)
      } else if (method == "projection") {
        temp_mu <- crossprod(X, beta)
      } else {
        stop("Error in calculating mu. method not found")
      }
      # tsortmu <- t(sort_mu)
      # if(method == "projection") {
        # transp <- transport_plan(sortmu, temp_mu, p, p, "colwise", "exact")
      wp <- WpProj::wasserstein(X = sort_mu, Y = temp_mu, 
                                         p = p, ground_p = ground_p, 
                                         observation.orientation = obs.direction, 
                                         method = OToptions$transport.method, 
                                         epsilon = OToptions$epsilon, niter = OToptions$niter)
      # } else {
      #   wp <- WpProj::wasserstein(sort_mu, temp_mu, p = p, ground_p = p, "colwise", 
      #                                      method=transport.method)
      # }
      return(list(wp = wp, beta = beta))
    }
    minus.idx <- function(j, in.idx = NULL, X = NULL, sort_mu = NULL, 
                          p, ground_p, OToptions, obs.direction,...) {
      temp.in.idx <- in.idx
      temp.in.idx[ j ] <- FALSE
      idx <- which( temp.in.idx )
      beta <- calc.beta(xtx, xty,idx, method, OToptions, X_, theta)
      d <- length(idx)
      # tsortmu <- t(sort_mu)
      if(method != "projection"){
        # beta <- theta %*% diag(beta)
        temp_mu <- selVarMeanGen(X, theta, beta)
      } else {
        temp_mu <- crossprod(X, beta)
      }
      # tsortmu <- t(sort_mu)
      # if(method == "projection") {
        # transp <- transport_plan(tsortmu, temp_mu, p, p, "colwise", "exact")
      wp <- WpProj::wasserstein(X = sort_mu, Y = temp_mu, 
                                         p = p, ground_p = ground_p, 
                                         observation.orientation = obs.direction, 
                                         method = OToptions$transport.method, 
                                         epsilon = OToptions$epsilon, niter = OToptions$niter)
      # } else {
      #   wp <- WpProj::wasserstein(sort_mu, temp_mu, p = p, ground_p = p, "colwise", 
      #                                      method=transport.method)
      # }
      return (list(wp = wp, beta = beta))
    }
    suffstat <- sufficientStatistics(X, Y_, t(theta), OToptions)
    xtx <- suffstat$XtX #* wt + diag(theta_norm) * (1-wt)
    xty <- suffstat$XtY #* wt + theta_norm * (1-wt)
  }
  
  if(method == "projection") {
    beta_store <- matrix(NA, nrow=S*d, ncol=max_iter)
  } else {
    beta_store <- matrix(NA, nrow=d, ncol=max_iter)
  }
  # if (grepl("univariate", transport.method ) ) {
  #   Y_ <- t(Y_)
  # }
  
  # if (grepl("univariate", transport.method) ) {
  #   obs.direction <- "rowwise"
  #   # X_ <- t(X_)
  #   # Y_ <- t(Y_)
  #   } else {
    obs.direction <- "colwise"
  # }
  
  if(display.progress){
    pb <- txtProgressBar(min = 0, max = max_iter, style = 3)
  }
  if (direction == "forward") {
    in.idx <- rep(FALSE,d)
    in.idx[force] <- TRUE
    wP <- rep(Inf,d)
    temp_idx <- which(in.idx)
    temp_mu <- crossprod(X_[temp_idx, , drop=FALSE], theta[temp_idx, ,drop=FALSE])
    wP_traj[1] <- WpProj::wasserstein(temp_mu, Y_, p, ground_p, obs.direction, transport.method, epsilon = epsilon, niter = OTmaxit)
    # wP_traj[1] <- mean((Y_ - temp_mu)^2)
     
    for(i in 1:max_iter){
      candidates <- which(!in.idx & not.force.logical )
      wP_list <- foreach::foreach(cand = candidates) %dopar% {
        return(add.idx(cand, in.idx = in.idx, X= X_, sort_mu = Y_, p = p, 
                       ground_p = ground_p, OToptions = OToptions, 
                       obs.direction = obs.direction,
                xtx = xtx, xty = xty, theta = theta))
      } #function(j, in.idx = NULL, X = NULL, sort_mu = NULL, p, ground_p, OToptions, obs.direction,...)
      wP <- sapply(wP_list, function(f) f$wp)
      min_cand <- which.min(wP)
      add <- candidates[min_cand]
      in.idx[add] <- TRUE
      sel.idx[i] <- add
      beta_store[,i] <- c(wP_list[[min_cand]]$beta)
      wP_traj[i+1] <- wP[min_cand] 
      if(!is.null(model.size)) if((l_force + i) == model.size) break
      if(display.progress) setTxtProgressBar(pb, i)
    }
  }
  if(direction == "backward") {
    
    in.idx <- rep(TRUE,d)
    wP <- rep(0,d)
    
    for(i in 1:max_iter){
      candidates <- which( in.idx & not.force.logical )
      wP_list <- foreach::foreach(cand = candidates) %dopar% {
        return(minus.idx(cand, in.idx = in.idx, X=X_, sort_mu = Y_, p = p, 
                         ground_p = ground_p,
                         OToptions = OToptions, obs.direction = obs.direction,
                         xtx = xtx, xty = xty, theta = theta))
      } #function(j, in.idx = NULL, X = NULL, sort_mu = NULL, p, ground_p, OToptions, obs.direction,...)
      # wP_list <- lapply(candidates, minus.idx, in.idx, X_, Y_, p = p, ground_p = ground_p,
      #                   OToptions = OToptions, obs.direction = obs.direction,
      #                   xtx = xtx, xty = xty, theta = theta)
      wP <- sapply(wP_list, function(f) f$wp)
      min_cand <- which.min( wP )
      remove <- candidates[min_cand]
      in.idx[remove] <- FALSE
      sel.idx[i] <- remove
      beta_store[,max_iter - i + 1] <- c(wP_list[[min_cand]]$beta)
      wP_traj[max_iter - i + 1] <- wP[min_cand]
      if(!is.null(model.size)) if(d-i == model.size) break
      if(display.progress) setTxtProgressBar(pb, i)
    }
  }
  if(display.progress) close(pb)
  wP_traj[max_iter + 1] <- 0
  sel.idx <- sel.idx[!is.na(sel.idx)]
  num_coef <- (0 + max(l_force,1)):max_iter
  indices <- if(direction=="forward") {
    lapply(seq_along(sel.idx), function(i) sort(c(force, sel.idx[1:i])))
  } else {
    lapply(seq_along(sel.idx), function(i) sort(c(force, sel.idx[(length(sel.idx)-i+1):length(sel.idx)])))
  }
  if(l_force != 0) {
    indices <- c(list(force), indices)
  }
  
  
  # if(direction == "backward") {
  #   beta_store <- beta_store[,rev(1:ncol(beta_store))]
  # }
  if(direction == "backward") {
    beta <- calc.beta(xtx, xty,1:ncol(X), method, OToptions, X_, theta)
    beta_store <- cbind(beta_store, c(beta))
    num_coef <- c(num_coef, ncol(X))
  } else if( direction == "forward" ) {
    if(!is.null(force)) {
      beta <- calc.beta(xtx, xty,force, method, OToptions, X_, theta)
      beta_store <- cbind(c(beta), beta_store)
      num_coef <- c(num_coef, max_iter + l_force)
    }
    # beta <- calc.beta(xtx, xty,1:ncol(X), method, OToptions, X_, theta)
    # beta_store <- cbind(beta_store, c(beta))
    
  }
  num_coef <- num_coef[apply(beta_store,2,function(x) all(!is.na(x)))]
  beta_store <- beta_store[,apply(beta_store,2,function(x) all(!is.na(x)))]
  output <- list(index = indices, 
                 path = sel.idx, wP = wP_traj, p = p, 
                 nzero=num_coef, force = force, 
                 beta= beta_store, call = formals(WPSW), 
                 method=method, direction = direction)
  output$call[names(this.call)] <- this.call
  class(output) <- c("WpProj","stepwise")
  output$method <- method
  if(calc.theta) {
    extract <- extractTheta(output, theta)
    output$theta <- extract$theta
    output$eta <- lapply(output$theta, function(tt) X %*% tt)
  }
  
  return(output)
}
