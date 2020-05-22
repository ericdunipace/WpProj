W1L1 <- function(X, Y, theta = NULL, penalty = c("none", "lasso","scad","mcp"), 
                 model.size = NULL,
                 lambda = numeric(0), 
                 lambda.min.ratio = 1e-4, 
                 nlambda = 10, 
                 gamma = 1, ...) {
  
  this.call <- as.list(match.call()[-1])
  
  # if(penalty == "lasso") stop("Lasso group penalty is currently incorrect in rqPen package!")
  penalty <- match.arg(penalty, choices = c("none","lasso","scad","mcp"))
  if(penalty == "ols") penalty <- "none"
  
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
                                                        x = c(X),
                                                        dims = c(n*s, d)))
  Xmat <- do.call(cbind, cols)
  rm(cols)
  
  if(length(lambda) == 0) {
    if(penalty != "lasso") {
      lambda.max <- max(colSums(abs(X))/n)
    } else {
      lambda.max <- max(sqrt(colSums(X^2)))
    }
    lambda <-  exp(log(lambda.max) + seq(0, log(lambda.min.ratio), length.out = nlambda))
  }
  if(length(lambda) == 1) if(lambda == 0) penalty <- "none"
  
  
  if(is.null(model.size) | length(model.size) == 0) {
    model.size <- ncol(Xmat)
  } else {
    model.size <- model.size * s
  }
  
  
  if(penalty != "none") {
    rqpenargs <-  list(
      x = as.matrix(Xmat),
      y = c(Y),
      groups = rep(1:d, s),
      tau = 0.5,
      intercept = FALSE,
      a = gamma,
      model.size = model.size,
      penalty = switch(penalty, "lasso" = "LASSO", 
                       "mcp" = "MCP",
                       "scad" = "SCAD"),
      lambda = lambda,
      ...)
    rqpenargs <- rqpenargs[!duplicated(names(rqpenargs))]
    
    if(is.null(rqpenargs$method)) {
      rqpenargs$method <- if(nrow(Xmat) > 1000 & ncol(Xmat) > 100) { #sfn gives error
        #   "sfn"
        #   rqargs$x <- as(Xmat, "matrix.csr")
        # } else if(nrow(Xmat) > 1000 & ncol(Xmat) > 100) {
        "pfn"
      } else if(nrow(Xmat) < 1000 & ncol(Xmat) > 100) {
        "fn"
      } else if(nrow(Xmat)  < 1000 & ncol(Xmat) < 100) {
        "br"
      }
    }
    
    
    if(length(lambda) ==1 ) {
      rqpenargs <- rqpenargs[names(rqpenargs) %in% names(c(formals(l1.group.fit), 
                                                             formals(quantreg::rq.fit.pfn), 
                                                             formals(quantreg::rq.fit.br), 
                                                             formals(quantreg::rq.fit.fnb)))]
      if(is.null(rqpenargs$alg)) rqpenargs$alg <- "QICD"
      argn <- lapply(names(rqpenargs), as.name) 
      names(argn) <- names(rqpenargs)
      
      f.call <- as.call(c(list(call("::", as.name("rqPen"), 
                                    as.name("l1.group.fit"))), argn))
      # res <- do.call(l1.group.fit, rqpenargs)
      # res <- l1.group.fit(x = rqpenargs$x, y = rqpenargs$y, 
      #                            groups = rqpenargs$groups, tau = 0.5, lambda = rqpenargs$lambda,
      #                            intercept = FALSE, 
      #                            penalty = rqpenargs$penalty, 
      #                            alg = rqpenargs$alg, penGroups = NULL, ...)
    } else {
      rqpenargs <- rqpenargs[names(rqpenargs) %in% c(names(c(formals(l1.group.fit), 
                                                             formals(quantreg::rq.fit.pfn), 
                                                             formals(quantreg::rq.fit.br), 
                                                             formals(quantreg::rq.fit.fnb))), 
                                                     "model.size")]
      # if(is.null(rqpenargs$alg)) rqpenargs$alg <- "QICD_warm"
      argn <- lapply(names(rqpenargs), as.name) 
      names(argn) <- names(rqpenargs) 
      
      f.call <- as.call(c(list(as.name("rqGroupLambda")), argn))
      # res1 <- rqGroupLambda(x = rqpenargs$x, y = rqpenargs$y,
      #                      groups = rqpenargs$groups, tau = 0.5, lambda = rqpenargs$lambda,
      #                      intercept = FALSE,
      #                      penalty = rqpenargs$penalty,
      #                      a = rqpenargs$a,
      #                      # alg = rqpenargs$alg,
      #                      ...)
      
    }
    res <- eval(f.call, envir = rqpenargs)
    beta <- sapply(res, function(r) r$coefficients)
  } else {
    lambda <- 0
    nlambda <- 0
    rqargs <- list(
      x = as.matrix(Xmat),
      y = c(Y),
      tau = 0.5,
      ...
    )
    # rqargs <- list(formula = formula("Y ~ . + 0"),
    #   data = data.frame(Y = c(Y), X = as.matrix(Xmat)),
    #   tau = 0.5,
    #   ...
    # )
    rqargs <- rqargs[!duplicated(names(rqargs))]
    if(is.null(rqargs$method)) {
      rqargs$method <- if(nrow(Xmat) > 1000 & ncol(Xmat) > 100) { #sfn gives error
        #   "sfn"
        #   rqargs$x <- as(Xmat, "matrix.csr")
        # } else if(nrow(Xmat) > 1000 & ncol(Xmat) > 100) {
        "pfn"
      } else if(nrow(Xmat) < 1000 & ncol(Xmat) > 100) {
        "fn"
      } else if(nrow(Xmat)  < 1000 & ncol(Xmat) < 100) {
        "br"
      }
    }
    arg.names <- switch(rqargs$method,
                        "pfn" = names(c(formals(quantreg::rq.fit.pfn))), 
                        "br" = names(formals(quantreg::rq.fit.br)), 
                        "fn" = names(formals(quantreg::rq.fit.fnb)))
    rqargs <- rqargs[names(rqargs) %in% arg.names ]
    
    argn <- lapply(names(rqargs), as.name)
    names(argn) <- names(rqargs)
    f.call <- as.call(c(list(call("::", as.name("quantreg"), 
                                as.name("rq.fit"))), argn))
    res <- eval(f.call, envir = rqargs)
    # res <- do.call(quantreg::rq, rqargs)
    
    beta <- as.matrix(res$coefficients)
  }
  output <- list()
  output$beta <- beta
  output$penalty <- penalty
  output$lambda <- lambda
  output$nvars <- d
  output$maxit <- NULL
  output$call <- formals(W1L1)
  output$call[names(this.call)] <- this.call
  output$nonzero_beta <- colSums(output$beta != 0)
  output$method <- "projection"
  output$power <- 1
  class(output) <- c("limbs", "optimization")
  
  extract <- extractTheta(output, matrix(0, d,s))
  output$nzero <- extract$nzero
  output$eta <- lapply(extract$theta, function(tt) X %*%tt )
  output$theta <- extract$theta
  output$model <- res
  return(output)
} 