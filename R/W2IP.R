W2IP <- function(X, Y=NULL, theta,
                 transport.method = transport_options(),
                 model.size = NULL,
                 infimum.maxit = 100,
                 tol = 1e-7,
                 solution.method = c("cone","lp", "cplex", "gurobi","mosek"),
                 display.progress=FALSE, parallel = NULL, ...) 
{
  this.call <- as.list(match.call()[-1])
  
  # `%doRNG%`` <- doRNG::`%dorng%`
  
  dots <- list(...)
  if(!is.matrix(X)) X <- as.matrix(X)
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
    model.size <- 1:p
  }
  p_star <- length(model.size)
  
  if(is.null(transport.method)){
    transport.method <- "exact"
  } else {
    transport.method <- match.arg(transport.method, transport_options())
  }
  
  if(is.null(solution.method)) {
    solution.method <- "cone"
  } else {
    solution.method <- match.arg(solution.method)
  }
  
  
  translate <- function(QP, solution.method) {
    switch(solution.method, 
           cone = ROI::ROI_reformulate(QP,to = "socp"),
           lp = ROI::ROI_reformulate(QP,"lp",method = "bqp_to_lp" ),
           cplex = QP,
           gurobi = QP,
           mosek = QP
    )
  }
  
  solver <- function(obj, control, solution.method) {
    switch(solution.method, 
           cone = ROI.plugin.ecos:::solve_OP(obj, control),
           lp =  ROI.plugin.lpsolve:::solve_OP(obj, control),
           cplex = ROI.plugin.cplex:::solve_OP(obj, control),
           gurobi = gurobi_solver(obj, control),
           mosek = mosek_solver(obj, control)
    )
  }
  
  if(ncol(theta) == ncol(X)){
    theta_ <- t(theta)
  } else {
    theta_ <- theta
  }
  if(nrow(theta) != p) stop("dimensions of theta must match X")
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
  if(is.null(OTmaxit)) OTmaxit <- 100
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
  model.size <- as.integer(model.size)
  
  if (infm.maxit <=0) {
    stop("infimum.maxit should be greater than 0")
  }
  
  if(!is.null(parallel)){
    if(!inherits(parallel, "cluster")) {
      stop("parallel must be a registered cluster backend")
    }
    doParallel::registerDoParallel(parallel)
    display.progress <- FALSE
  } else{
    foreach::registerDoSEQ()
  }
  
  options <- list(infm_maxit = infm.maxit,
                  display_progress = display.progress, 
                  model_size = model.size)
  OToptions <- list(same = same,
                    method = "selection.variable",
                    transport.method = transport.method,
                    epsilon = epsilon,
                    niter = OTmaxit)
  
  ss <- limbs:::sufficientStatistics(X, Y_, theta_, OToptions)
  xtx <- ss$XtX
  xty <- xty_init <- ss$XtY
  Ytemp <- Y_
  
  if(display.progress){
    pb <- txtProgressBar(min = 0, max = p_star, style = 3)
    setTxtProgressBar(pb, 0)
  }
  
  QP <- QP_orig <- limbs:::qp_w2(ss$XtX,ss$XtY,1)
  # LP <- ROI::ROI_reformulate(QP,"lp",method = "bqp_to_lp" )
  alpha <- alpha_save <- rep(0,p)
  beta <- matrix(0, nrow = p, ncol = p_star)
  iter.seq <- rep(0, p_star)
  comb <- function(x, ...) {
    # from https://stackoverflow.com/questions/19791609/saving-multiple-outputs-of-foreach-dopar-loop
    lapply(seq_along(x),
           function(i) c(x[[i]], lapply(list(...), function(y) y[[i]]))
    )
  }
  
  output <- foreach::foreach(idx=1:p_star, .combine='comb', .multicombine=TRUE,
                             .init=list(list(), list()),
                             .errorhandling = 'pass', 
                             .inorder = FALSE) %dorng% 
    {
       m <- options$model_size[idx]
       QP <- QP_orig
       QP$constraints$rhs[1] <- m
       results <- list(NULL, NULL)
       for(inf in 1:options$infm_maxit) {
         TP <- translate(QP, solution.method)
         # browser()
         # sol <- ROI::ROI_solve(LP, "glpk")
         # can use ROI.plugin.glpk:::.onLoad("ROI.plugin.glpk","ROI.plugin.glpk") to use base solver ^
         sol <- solver(TP, control, solution.method)
         # print(ROI::solution(sol))
         alpha <- ROI::solution(sol)[1:p]
         if(all(is.na(alpha))) {
           warning("Likely terminated early")
           break
         }
         if(limbs:::not.converged(alpha, alpha_save, tol)){
           alpha_save <- alpha
           Ytemp <- limbs:::selVarMeanGen(X_, theta_, alpha)
           xty <- limbs:::xtyUpdate(X, Ytemp, theta_, result_ = alpha, 
                                             OToptions)
           QP <- limbs:::qp_w2(xtx,xty,m)
         } else {
           break
         }
       }
       if(display.progress) setTxtProgressBar(pb, idx)
       results[[2]] <- inf
       results[[1]] <- alpha
       return(results)
       # iter.seq[idx] <- inf
       # if ( same ) QP <- QP_orig
       # QP$constraints$rhs[1] <- m + 1
       # beta[,idx] <- alpha
    }
  if (display.progress) close(pb)
  names(output) <- c("beta","iter")
  output$beta <- do.call("cbind", output$beta)
  output$iter <- unlist(output$iter)
  
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
  class(output) <- c("limbs","IP")
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
  return(op)
}


