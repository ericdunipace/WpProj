# library(slam)
# Sys.setenv(ROI_LOAD_PLUGINS = FALSE)
# library(ROI)
# library(ROI.plugin.glpk)
library(SparsePosterior)

set.seed(84370158)

n <- 100
p <- 10
s <- 1000

x <- matrix( rnorm( p * n ), nrow = n, ncol = p )
x_ <- t(x)
beta <- (1:10)/10
y <- x %*% beta + rnorm(n)
post_beta <- matrix(beta, nrow=p, ncol=s) + rnorm(p*s, 0, 0.1)
post_mu <- x %*% post_beta

xtx <- crossprod(x)/n #* wt + diag(1,p,p) * (1 - wt)
xty <- crossprod(x, post_mu)/n #* wt + post_beta * (1 - wt)


#exact
transp <- "exact"
# suffStat <- sufficientStatistics(x, post_mu, t(post_beta), TRUE, "selection.variable",transp)
# xtx <- suffStat$XtX #* wt + diag(post_beta_norm) * (1-wt)
# xty <- suffStat$XtY #* wt + post_beta_norm * (1-wt)

poss_means <- lapply(1:p,function(i) crossprod(x_[i,,drop=FALSE], post_beta[i,,drop=FALSE]))
possible <- sapply(poss_means, function(pm) wasserstein(post_mu, pm,
                                                2,2,"colwise","exact"))
which.min(possible)

size2 <- combn(1:p,2)
poss_means2 <- lapply(1:ncol(size2),function(i) crossprod(x_[size2[,i],,drop=FALSE], post_beta[size2[,i],,drop=FALSE]))
possible2 <- sapply(poss_means2, function(pm) wasserstein(post_mu, pm,
                                                        2,2,"colwise","exact"))
size2[,which.min(possible2)]
# dbind <- function(...) {
#   .dbind <- function(x, y) {
#     A <- simple_triplet_zero_matrix(NROW(x), NCOL(y))
#     B <- simple_triplet_zero_matrix(NROW(y), NCOL(x))
#     rbind(cbind(x, A), cbind(B, y))
#   }
#   Reduce(.dbind, list(...))
# }

# qp_w2 <- function(xtx, xty, K) {
#   stzm <- simple_triplet_zero_matrix
#   stdm <- simple_triplet_diag_matrix
#   d <- NCOL(xtx)
#   Q0 <- xtx
#   # Q0 <- dbind(stzm(n), stdm(1, m), stzm(n))
#   # a0 <- c(a = double(n), g = double(m), t = lambda * rep(1, n))
#   L0 <- c(a = c(-2*xty))
#   op <- OP(objective = Q_objective(Q = Q0, L = L0, names = as.character(1:d)))
#   ## sum(alpha) = K 
#   A1 <- rep(1,d) # beta portion, gamma porition, t portion is 0. multiplies all vars each time!!!
#   LC1 <- L_constraint(A1, eq(1), K)
#   # ##  -t <= beta  <=>  0 <= beta + t
#   # A2 <- cbind(stdm(1, n), stzm(n, m), stdm(1, n)) #beta gamma t
#   # LC2 <- L_constraint(A2, geq(n), double(n))  #reframes in terms of constant!!!
#   # ##   beta <= t  <=>  beta - t <= 0
#   # A3 <- cbind(stdm(1, n), stzm(n, m), stdm(-1, n))
#   # LC3 <- L_constraint(A3, leq(n), double(n))
#   constraints(op) <- LC1 #rbind(LC1, LC2, LC3)
#   # bounds(op) <- V_bound(li = 1:d, lb = rep.int(0, d), ui = 1:d, ub = rep.int(1, d))
#   types(op) <- rep.int("B", d)
#   op
# }
# 
# op1 <- qp_w2(xtx, xty, 1)
# ROI:::get_solver_methods(OP_signature(op1))
# lp1 <- ROI_reformulate(op1, "lp")
# ROI:::get_solver_methods(OP_signature(lp1))
# (qp1 <- ROI_solve(op1, "glpk"))
# (qp1 <- ROI_solve(lp1, "glpk"))
# solution(qp1)[1:10]
# 
# op3 <- qp_w2(xtx, xty, 3)
# lp3 <- ROI_reformulate(op3, "lp")
# (qp3 <- ROI_solve(lp3, "glpk"))
# 
# op9 <- qp_w2(xtx, xty, 9)
# lp9 <- ROI_reformulate(op9, "lp")
# (qp9 <- ROI_solve(lp9, "glpk"))
# 
# debugonce(W2IP)
test <- W2IP(X = x, Y = post_mu, theta = post_beta, transport.method = transp, 
     infimum.maxit = 10, 
     tol = 1e-7,
     display.progress = TRUE,model.size = 2)

test0 <- W2L1(X = x, Y = post_mu, theta = post_beta, transport.method = transp, 
              method = "selection.variable",
              penalty = "selection.lasso",
              infimum.maxit = 10, 
              tol = 1e-7, 
              display.progress = TRUE)

test1 <- W2L1(X = x, Y = post_mu, theta = post_beta, transport.method = transp, 
              method = "selection.variable",
              penalty = "selection.lasso",
              infimum.maxit = 10, 
              tol = 1e-7, penalty.factor = 1/rowMeans(post_beta^2),
              display.progress = TRUE)
test2 <- W2L1(X = x, Y = post_mu, theta = post_beta, transport.method = "univariate.approximation.pwr", 
              method = "selection.variable",
              penalty = "selection.lasso",
              infimum.maxit = 10, 
              tol = 1e-7, penalty.factor = 1/rowMeans(post_beta^2),
              display.progress = TRUE)

cbind(IP=lapply(test$theta, function(tt) which(rowSums(tt) > 0))[[1]],
      exact = lapply(test0$theta, function(tt) which(rowSums(tt) > 0))[[2]],
      exact.wt = lapply(test1$theta, function(tt) which(rowSums(tt) > 0))[[2]],
      uap.wt = lapply(test2$theta, function(tt) which(rowSums(tt) > 0))[[2]],
      true = which.min(possible)
      )

cbind(IP=lapply(test$theta, function(tt) which(rowSums(tt) > 0))[[2]],
      exact = lapply(test0$theta, function(tt) which(rowSums(tt) > 0))[[3]],
      exact.wt = lapply(test1$theta, function(tt) which(rowSums(tt) > 0))[[3]],
      uap.wt = lapply(test2$theta, function(tt) which(rowSums(tt) > 0))[[3]],
      true = size2[,which.min(possible2)]
)
