#retired R files
# .W2L1R <- function(X, Y=NULL, theta, family="gaussian",
#                     penalty = c("elastic.net", "selection.lasso",
#                                  "lasso", "ols", "mcp",
#                                  "scad", "mcp.net",
#                                  "scad.net",
#                                  "grp.lasso",
#                                  "grp.lasso.net", "grp.mcp",
#                                  "grp.scad", "grp.mcp.net",
#                                  "grp.scad.net",
#                                  "sparse.grp.lasso"),
#                      lambda = numeric(0),
#                      nlambda = 100L,
#                      lambda.min.ratio = NULL, alpha = 1,
#                      gamma = 1, tau = 0.5,
#                      groups = numeric(0),
#                      scale.factor = numeric(0),
#                      penalty.factor = NULL,
#                      group.weights = NULL, maxit = 500L,
#                      tol = 1e-07, irls.maxit = 100L,
#                      irls.tol = 0.001,
#                      infimum.maxit = NULL)
# {
# 
#   this.call <- match.call()
# 
#   if (is.null(lambda.min.ratio)) {
#     lambda.min.ratio <- 1e-04
#   }
#   else {
#     lambda.min.ratio <- as.numeric(lambda.min.ratio)
#   }
#   if (lambda.min.ratio >= 1 | lambda.min.ratio <= 0) {
#     stop("lambda.min.ratio must be between 0 and 1")
#   }
# 
#   if(is.null(Y)) {
#     same <- TRUE
#     Y <- tcrossprod(X, theta)
#   }
#   data <- augmentMatrix(X,Y,theta,same)
#   xtx <- data$XtX
#   xty <- data$XtY
# 
#   if(! is.list(lambda)){
#     lambda.max <- max(abs(xty))
#     lambda.min <- lambda.max * lambda.min.ratio
#     lambda <- list(exp(seq(log(lambda.max), log(lambda.min), length.out = nlambda)))
#   }
# 
#   if(is.null(infm.maxit)) infm.maxit <- maxit
# 
#   output <- list(xtx=xtx,
#                  xty=xty,
#                  gamma = NULL,
#                  gamma_old=rep(1e6, ncol(xtx)),
#                  beta = matrix(NA, ncol=nlambda,nrow=ncol(xtx)))
#   alt.iter <- rep(NA,nlambda)
# 
#   for(l in seq_along(lambda[[1]])){
#     notconv <- TRUE
#     i<-0
#     while((i < infm.maxit) & notconv) {
#       oemfit <- oem.xtx(output$xtx,output$xty, family,
#                         penalty,
#                         list(lambda[[1]][l]),
#                         nlambda,
#                         lambda.min.ratio , alpha,
#                         gamma, tau,
#                         groups,
#                         scale.factor,
#                         penalty.factor,
#                         group.weights, maxit,
#                         tol, irls.maxit,
#                         irls.tol)
#       output$gamma <- drop(oemfit$beta[[1]])
#       if(any(output$gamma != 0)) {
#         output$xty <- xtyUpdate(X, Y,
#                               theta, output$gamma)
#       }
#       notconv <- not.converged(output$gamma, output$gamma_old, tol)
# 
#       if(notconv) {
#         output$gamma_old <- output$gamma
#       } else {
#         output$gamma_old <- rep(1e6, ncol(xtx))
#         output$beta[,l] <- output$gamma
#       }
#       i <- i + 1
#     }
#     alt.iter[l] <- i
# 
#   }
# 
#   return(list(beta=output$beta, iter=alt.iter))
# 
# }
