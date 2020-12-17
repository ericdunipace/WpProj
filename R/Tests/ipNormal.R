# library(slam)
# Sys.setenv(ROI_LOAD_PLUGINS = FALSE)
# library(ROI)
# library(ROI.plugin.glpk)
library(SparsePosterior)
library(CoarsePosteriorSummary)

set.seed(84370158)

n <- 1024
p <- 101
s <- 100

target <- get_normal_linear_model()

target$X$corr <- 0.9
x <- target$X$rX(n, target$X$corr , p)
theta_star <- target$rparam()
data <- target$rdata(n,x,theta_star$theta,theta_star$sigma2)
y <- data$Y

hyperparameters <- list(mu = rep(0.0,p), sigma = diag(p), alpha = as.double(1.0), beta = as.double(1.0))
posterior <- target$rpost(n.samp = as.integer(s), x=x, y=y, hyperparameters = hyperparameters, method = "conjugate")

theta <- (posterior$theta)
eta <- posterior$eta

#exact
transp <- "gandkhorn"
# suffStat <- sufficientStatistics(x, post_mu, t(post_beta), TRUE, "selection.variable",transp)
# xtx <- suffStat$XtX #* wt + diag(post_beta_norm) * (1-wt)
# xty <- suffStat$XtY #* wt + post_beta_norm * (1-wt)

poss_means <- lapply(1:p,function(i) crossprod(t(x[,i,drop=FALSE]), theta[i,,drop=FALSE]))
possible <- sapply(poss_means, function(pm) wasserstein(eta, pm,
                                                        2,2,"colwise",transp))
which.min(possible)

size2 <- combn(1:p,2)
possible2 <- sapply(1:ncol(size2), function(i) wasserstein(
  crossprod(t(x[, size2[,i],drop=FALSE]), theta[size2[,i],,drop=FALSE]), 
  eta,
  2,2,"colwise",transp))
size2[,which.min(possible2)]

# debugonce(W2IP)
time <- proc.time()
test <- W2IP(X = x, Y = NULL, theta = theta, transport.method = transp, 
             infimum.maxit = 5, 
             tol = 1e-7, model.size = 2,
             solution.method = "cone",
             display.progress = FALSE,
             control = list(VERBOSE = 1L))
print(proc.time() - time)


cl <- parallel::makeCluster(3)
# debugonce(W2IP)
test <- W2IP(X = x, Y = NULL, theta = theta, transport.method = transp, 
             infimum.maxit = 5, 
             tol = 1e-7, model.size = 1:10,
             solution.method = "cone",
             display.progress = TRUE,
             parallel = cl)
parallel::stopCluster(cl)

# test0 <- W2L1(X = x, Y = post_mu, theta = post_beta, transport.method = transp, 
#               method = "selection.variable",
#               penalty = "selection.lasso",
#               infimum.maxit = 10, 
#               tol = 1e-7, 
#               display.progress = TRUE)

test1 <- W2L1(X = x, Y = eta, theta = theta, transport.method = transp, 
              method = "selection.variable",
              penalty = "selection.lasso",
              infimum.maxit = 10, nlambda = 1000,
              tol = 1e-7, penalty.factor = 1/rowMeans(theta^2),
              display.progress = TRUE)
test2 <- W2L1(X = x, Y = eta, theta = theta, transport.method = "univariate.approximation.pwr", 
              method = "selection.variable",
              penalty = "selection.lasso",
              infimum.maxit = 10, nlambda = 1000,
              tol = 1e-7, penalty.factor = 1/rowMeans(theta^2),
              display.progress = TRUE)

cbind(IP=lapply(test$theta, function(tt) which(rowSums(tt) !=0 ))[[1]],
      # exact = lapply(test0$theta, function(tt) which(rowSums(tt) > 0))[[2]],
      exact.wt = lapply(test1$theta, function(tt) which(rowSums(tt) !=0 ))[[2]],
      uap.wt = lapply(test2$theta, function(tt) which(rowSums(tt) != 0))[[2]]
      # true = which.min(possible)
)

cbind(IP=lapply(test$theta, function(tt) which(rowSums(tt) != 0))[[2]],
      # exact = lapply(test0$theta, function(tt) which(rowSums(tt) > 0))[[3]],
      exact.wt = lapply(test1$theta, function(tt) which(rowSums(tt) != 0))[[3]],
      uap.wt = lapply(test2$theta, function(tt) which(rowSums(tt) != 0))[[3]]
      # true = size2[,which.min(possible2)]
)

cl <- parallel::makeCluster(3)
debugonce(WPSA)
test <- WPSA(X = x, Y = NULL, theta = theta,
             force = NULL, p=2, model.size = 1:2, iter = 5, temps = 5,
             options = list(method = "selection.variable",
                            energy.distribution = "boltzman",
                            transport.method = "exact",
                            cooling.schedule="exponential",
                            proposal.method = "random"),
             display.progress=FALSE, max.time=Inf,
             parallel = cl, const = 34)
parallel::stopCluster(cl)

cl <- parallel::makeCluster(2)
# debugonce(WPSA)
test <- WPSA(X = x, Y = NULL, theta = theta,
             force = NULL, p=2, model.size = 1:2, iter = 5, temps = 5,
             options = list(method = "selection.variable",
                            energy.distribution = "boltzman",
                            transport.method = "exact",
                            cooling.schedule="exponential",
                            proposal.method = "random"),
             display.progress=FALSE, max.time=0.1,
             parallel = cl, const = 34)
parallel::stopCluster(cl)

debugonce(WPSA)
test <- WPSA(X = x, Y = NULL, theta = theta,
             force = NULL, p=2, model.size = 1, iter = 5, temps = 5,
             options = list(method = "selection.variable",
                            energy.distribution = "boltzman",
                            transport.method = "exact",
                            cooling.schedule="exponential",
                            proposal.method = "random"),
             display.progress=FALSE, max.time=Inf,
             const = 34)

test <- WPSA(X = x, Y = NULL, theta = theta,
             force = NULL, p=2, model.size = 1:2, iter = 5, temps = 5,
             options = list(method = "selection.variable",
                            energy.distribution = "boltzman",
                            transport.method = "exact",
                            cooling.schedule="exponential",
                            proposal.method = "random"),
             display.progress=FALSE, max.time=1,
             const = 34)

time <- proc.time()
test <- W2IP(X = x, Y = NULL, theta = theta, transport.method = transp, 
             infimum.maxit = 5, 
             tol = 1e-7, model.size = 51,
             solution.method = "cone",
             display.progress = FALSE,
             control = list(VERBOSE = 1L,
                            MI_MAX_ITERS=1e1L))
print(proc.time() - time)

time <- proc.time()
test <- W2IP(X = x, Y = NULL, theta = theta, transport.method = transp, 
             infimum.maxit = 5, 
             tol = 1e-7, model.size = 51,
             solution.method = "lp",
             display.progress = FALSE)
print(proc.time() - time)
