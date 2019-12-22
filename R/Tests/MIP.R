library(dplyr)
library(ROI)
library(ROI.plugin.glpk)
library(ompr)
library(ompr.roi)
library(SparsePosterior)

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
suffStat_star <- sufficientStatistics(x, post_mu, t(post_beta), TRUE,"selection.variable",transp)
xtx_star <- suffStat_star$XtX #* wt + diag(post_beta_norm) * (1-wt)
xty_star <- suffStat_star$XtY #* wt + post_beta_norm * (1-wt)

xstar <- t(t(x)*post_beta[1,])
ystar <- post_mu[,1,drop=FALSE]
result <- MIPModel() %>%
  add_variable(alpha[i], i = 1:p, type = "binary") %>%
  # set_bounds(alpha, lb = 0, ub=1) %>%
  # set_objective(-2*sum_expr(xty_star[i] * alpha[i],i=1:p) + 2*sum_expr(alpha[i] * xtx_star[i,j] * alpha[j], i = 1:p, j = 1:p), "min") %>%
  set_objective(sum_expr((ystar[i] -  sum_expr(xstar[i,j] * alpha[j], j = 1:p)), i = n) , "min") %>%
  add_constraint( sum_expr(alpha[i],i =1:p)  <= 2) %>%
  solve_model(with_ROI(solver = "glpk")) 
get_solution(result, alpha[i])
#can't do non-linear
