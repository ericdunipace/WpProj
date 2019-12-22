library(glmnet)
data(QuickStartExample)
lambda <- 1
mo <- glmnet(x, y, family = "gaussian", alpha = 1, lambda = lambda, 
             intercept = FALSE, standardize = FALSE)
glmnet_beta <- setNames(as.vector(coef(mo)), rownames(coef(mo)))
round(glmnet_beta, 4)

library(slam)
Sys.setenv(ROI_LOAD_PLUGINS = FALSE)
library(ROI)
library(ROI.plugin.qpoases)
library(ROI.plugin.ecos)

dbind <- function(...) {
  .dbind <- function(x, y) {
    A <- simple_triplet_zero_matrix(NROW(x), NCOL(y))
    B <- simple_triplet_zero_matrix(NROW(y), NCOL(x))
    rbind(cbind(x, A), cbind(B, y))
  }
  Reduce(.dbind, list(...))
}

qp_lasso <- function(x, y, lambda) {
  stzm <- simple_triplet_zero_matrix
  stdm <- simple_triplet_diag_matrix
  m <- NROW(x); n <- NCOL(x)
  Q0 <- dbind(stzm(n), stdm(1, m), stzm(n))
  a0 <- c(b = double(n), g = double(m), t = lambda * rep(1, n))
  op <- OP(objective = Q_objective(Q = Q0, L = a0))
  ## y - X %*% beta = gamma  <=>  X %*% beta + gamma = y
  A1 <- cbind(x, stdm(1, m), stzm(m, n)) # beta portion, gamma porition, t portion is 0. multiplies all vars each time!!!
  LC1 <- L_constraint(A1, eq(m), y)
  ##  -t <= beta  <=>  0 <= beta + t
  A2 <- cbind(stdm(1, n), stzm(n, m), stdm(1, n)) #beta gamma t
  LC2 <- L_constraint(A2, geq(n), double(n))  #reframes in terms of constant!!!
  ##   beta <= t  <=>  beta - t <= 0
  A3 <- cbind(stdm(1, n), stzm(n, m), stdm(-1, n))
  LC3 <- L_constraint(A3, leq(n), double(n))
  constraints(op) <- rbind(LC1, LC2, LC3)
  bounds(op) <- V_bound(ld = -Inf, nobj = ncol(Q0))
  op
}

op <- qp_lasso(x, y, 0)
(qp0 <- ROI_solve(op, "qpoases"))

op <- qp_lasso(x, y, lambda * NROW(x))
(qp1 <- ROI_solve(op, "qpoases"))

n <- ncol(x)
cbind(lm = coef(lm.fit(x, y)), qp = head(solution(qp0), n))

cbind(lm = round(glmnet_beta, 4), qp = round(c(0, head(solution(qp1), n)), 4))
