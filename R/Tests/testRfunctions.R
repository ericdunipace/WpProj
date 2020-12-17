rm(list=ls())
require(SparsePosterior)
set.seed(111)

##### Testing R Functions ####
n <- 100
p <- 10
s <- 1000

x <- matrix( rnorm( p * n ), nrow = n, ncol = p )
x_ <- t(x)
beta <- (1:10)/10
y <- x %*% beta + rnorm(n)
post_beta <- matrix(beta, nrow=p, ncol=s) + rnorm(p*s, 0, 0.1)
post_mu <- x %*% post_beta
post_beta_norm <- rowMeans(post_beta^2)
pseudo.obs <- 1
wt <- n/(pseudo.obs + n)

xtx <- crossprod(x)/n #* wt + diag(1,p,p) * (1 - wt)
xty <- crossprod(x, post_mu)/n #* wt + post_beta * (1 - wt)

suffStat_star <- sufficientStatistics(x, post_mu, t(post_beta), TRUE, 0.0, "selection.variable")
xtx_star <- suffStat_star$XtX #* wt + diag(post_beta_norm) * (1-wt)
xty_star <- suffStat_star$XtY #* wt + post_beta_norm * (1-wt)

active.idx <- seq(2,10,2)

#### Test convergence function ####
fake_coefs <- matrix(1:100, ncol=10, nrow=10)
fake_coefs2 <- matrix(100:1, ncol=10, nrow=10)

not.converged(fake_coefs, fake_coefs2, 1e-10) # TRUE
not.converged(fake_coefs, fake_coefs2, 100) # FALSE
not.converged(fake_coefs, fake_coefs, 1e-10) # FALSE

#### Test calculating beta ####
# selection.variable selection (should give 5 ones) 
  # shouldn't give any errors since it just takes number of active and spits out that number of ones

  #specify Y
  out <- calc.beta(xtx=xtx_star, xty=xty_star, active.idx = active.idx, 
                    method="selection.variable",
             x=x_, theta=post_beta,
             Y = post_mu, niter=5e2)
  print(out)
  sum(out)
  
  # dont' specify Y
  out <- calc.beta(xtx=xtx_star, xty=xty_star, active.idx = active.idx, 
                    method="selection.variable",
             x=x_, theta=post_beta, Y = NULL, niter=5e2)
  print(out)
  sum(out)
  
  # theta in right direction
  out <- calc.beta(xtx=xtx_star, xty=xty_star, active.idx = active.idx,
                    method="selection.variable",
                    x=x_, theta=post_beta, Y = NULL, niter=5e2)
  print(out)
  sum(out)
  
  # drop xtx or xty
  calc.beta(xtx=NULL, xty=xty_star, active.idx = active.idx, method="selection.variable",
             x=x_, theta=post_beta,
             Y = NULL, niter=5e2)
  calc.beta(xtx=xtx, xty=NULL, active.idx = active.idx, method="selection.variable",
             x=x_, theta=post_beta,
             Y = NULL, niter=5e2)

# scale selection should give 0.0000000 0.9351839 0.0000000 0.8069530 0.0000000 0.5354074 0.0000000 0.9519285 0.0000000 1.0024417
  #specify Y
  out <- calc.beta(xtx=xtx_star, xty=xty_star, 
                    active.idx = active.idx, method="scale",
                    x=x_, theta=post_beta,
                    Y = post_mu, niter=5e2)
  print(out)
  
  # dont' specify Y
  out <- calc.beta(xtx=xtx_star, xty=xty_star, active.idx = active.idx, 
                    method="scale",
                    x=x_, theta=post_beta, Y = NULL, niter=5e2)
  print(out)
  
  # theta in right direction
  out <- calc.beta(xtx=xtx_star, xty=xty_star, active.idx = active.idx, 
                    method="scale",
                    x=x_, theta=post_beta, Y = NULL, niter=5e2)
  print(out)
  
  # drop xtx or xty
  calc.beta(xtx=NULL, xty=xty_star, active.idx = active.idx, method="scale",
             x=x_, theta=post_beta, Y = NULL, niter=5e2)
  calc.beta(xtx=xtx_star, xty=NULL, active.idx = active.idx, method="scale",
             x=x_, theta=post_beta, Y = NULL, niter=5e2)
  
# calculate projection operator
  #specify Y
  out1 <- calc.beta(xtx=xtx, xty=xty, active.idx = active.idx, method="projection",
                    x=x_, theta=post_beta, Y = post_mu, niter=5e2)
  out2 <- calc.beta(xtx=crossprod(x)/n, xty=crossprod(x, post_mu)/n, 
                     active.idx = active.idx, method="projection",
                     x=x_, theta=post_beta, Y = post_mu, niter=5e2)
  out5 <- calc.beta(xtx=crossprod(x)/n, xty=crossprod(x, post_mu)/n, 
                     active.idx = active.idx, method="projection",
                     x=x_, theta=post_beta, Y = post_mu, niter=5e2)
  print(out1[1,1])
  print(out2[1,1])
  print(out5[1,1])
  print(out1[2,1])
  print(out2[2,1])
  print(out5[2,1])
  
  all.equal(out1,out2)
  all.equal(out2,out5)
  
  # Wasserstein_iid(post_mu, x[,active.idx] %*% out1,2, "rowwise" ) #full wass
  # Wasserstein_iid(post_mu, x[,active.idx] %*% out,2, "rowwise" ) #univariate approx
  # Wasserstein_iid(post_mu, x[,active.idx] %*% out2,2, "rowwise" ) #full wass, no pseudo
  # Wasserstein_iid(post_mu[1,, drop=FALSE], 
  #                 x[1,active.idx, drop=FALSE] %*% out3,
  #                 2, "rowwise" ) #univariate wass, heavy weighting toward post, single obs
  # Wasserstein_iid(post_mu[1,, drop=FALSE], 
  #                 x[1,active.idx, drop=FALSE] %*% out4,
  #                 2, "rowwise" ) #full wass, heavy weighting toward post, single obs
  # Wasserstein_iid(post_mu[,, drop=FALSE], 
  #                 x[,active.idx, drop=FALSE] %*% out5,
  #                 2, "rowwise" ) #KL proj
  # Wasserstein_iid(post_mu[1,, drop=FALSE], 
  #                 x[,active.idx, drop=FALSE] %*% out6,
  #                 2, "rowwise" ) #univariate wass, heavy weighting toward post
  # Wasserstein_iid(post_mu[,, drop=FALSE], 
  #                 x[,active.idx, drop=FALSE] %*% out7,
  #                 2, "rowwise" ) #full wass, heavy weighting toward post
  
  # transport::wasserstein(transport::pp(t(post_mu)), 
  #                        transport::pp(t(x[,active.idx] %*% out1)),2, 
  #                        method="shortsimplex") #full wass
  # transport::wasserstein(transport::pp(t(post_mu)), 
  #                        transport::pp(t(x[,active.idx] %*% out)),2,
  #                        method="shortsimplex") #univariate approx
  # transport::wasserstein(transport::pp(t(post_mu)), 
  #                        transport::pp(t(x[,active.idx] %*% out2)),2,
  #                        method="shortsimplex") #full wass no pseudo
  # transport::wasserstein(transport::pp(t(post_mu)), 
  #                        transport::pp(t(x[,active.idx] %*% out3)),2,
  #                        method="shortsimplex") #univariate wass, heavy weighting, single obs, toward post
  # transport::wasserstein(transport::pp(t(post_mu)), 
  #                        transport::pp(t(x[,active.idx] %*% out4)),2,
  #                        method="shortsimplex") #full wass, heavy weighting toward post, single obs
  # transport::wasserstein(transport::pp(t(post_mu)), 
  #                        transport::pp(t(x[,active.idx] %*% out6)),2,
  #                        method="shortsimplex") #univariate wass, heavy weighting toward post
  # transport::wasserstein(transport::pp(t(post_mu)), 
  #                        transport::pp(t(x[,active.idx] %*% out7)),2,
  #                        method="shortsimplex") #full wass, heavy weighting toward post
  # transport::wasserstein(transport::pp(t(post_mu)), 
  #                        transport::pp(t(x[,active.idx] %*% out5)),2,
  #                        method="shortsimplex") # KL
  
  
  # dont' specify Y
  out <- calc.beta(xtx=xtx, xty=NULL, active.idx = active.idx, method="projection",
                    x=x_, theta=post_beta,
                    Y = NULL, niter=5e2)
  print(out)
  
  # theta in right direction
  out <- calc.beta(xtx=xtx, xty=NULL, active.idx = active.idx, method="projection",
                    x=x_, theta=post_beta,
                    Y = NULL, niter=5e2)
  print(out)
  
  # drop xtx or xty
  calc.beta(xtx=NULL, xty=xty, active.idx = active.idx, method="projection",
             x=x_, theta=post_beta,
             Y = NULL, niter=5e2)
  calc.beta(xtx=xtx, xty=NULL, active.idx = active.idx, method="projection",
             x=x_, theta=post_beta, 
             Y = NULL, niter=5e2)
  
  # calculate beta adaptively
  
  #specify Y
  out <- calc.beta(xtx=xtx, xty=NULL, active.idx = active.idx, method="projection",
                    x=x_, theta=post_beta,
                    Y = post_mu, niter=5e2)
  print(out)
  
  # dont' specify Y
  out <- calc.beta(xtx=xtx, xty=NULL, active.idx = active.idx, method="projection",
                    x=x_, theta=post_beta,
                    Y = NULL, niter=5e2)
  print(out)
  
  # theta in right direction
  out <- calc.beta(xtx=xtx, xty=NULL, active.idx = active.idx, method="projection",
                    x=x_, theta=post_beta, 
                    Y = NULL, niter=5e2)
  print(out)
  
  # drop xtx or xty
  calc.beta(xtx=NULL, xty=crossprod(x, post_mu)/n, active.idx = active.idx, method="projection",
             x=x_, theta=post_beta, 
             Y = NULL, niter=5e2)
  calc.beta(xtx=xtx, xty=NULL, active.idx = active.idx, method="projection",
             x=x_, theta=post_beta, 
             Y = NULL, niter=5e2)
  
  # check that match arg works
  calc.beta(xtx=xtx_star, xty=NULL, active.idx = active.idx, method="sel",
             x=x_, theta=post_beta, 
             Y = NULL, niter=5e2)
  calc.beta(xtx=xtx, xty=NULL, active.idx = active.idx, method="pr",
            x=x_, theta=t(post_beta), 
            Y = NULL, niter=5e2)

  
#### Stepwise Checks ####
  # debugonce(WPSW)
  out <- WPSW(x, post_mu, t(post_beta), force = 1, p = 2, ground_p = 2,
           direction = c("backward"), 
           method=c("selection.variable"))
  out1 <- WPSW(x, post_mu, t(post_beta), force = 1, p = 2, ground_p = 2,
                     direction = c("forward"), 
                     method=c("selection.variable"))
  #check model size truncation
  out <- WPSW(x, t(post_mu), t(post_beta), force = 1, p = 2, ground_p = 2,
              direction = c("backward"), 
              method=c("selection.variable"), model.size = 8)
  which(apply(out$beta, 2, function(i) all(!is.na(i))))
  out1 <- WPSW(x, t(post_mu), t(post_beta), force = 1, p = 2,
               direction = c("forward"), 
               method=c("selection.variable"), model.size = 3)
  which(apply(out1$beta, 2, function(i) all(!is.na(i))))
  
  #check projection
  out2 <- WPSW(x, t(post_mu), t(post_beta), force = 1, p = 2,
                     direction = c("backward"), 
                     method=c("projection"))
  out3 <- WPSW(x, t(post_mu), t(post_beta), force = 1, p = 2,
                     direction = c("forward"), 
                     method=c("projection"))
  #check modelsize truncation
  out2 <- WPSW(x, t(post_mu), t(post_beta), force = 1, p = 2,
               direction = c("backward"), 
               method=c("projection"), model.size = 8)
  which(apply(out2$beta,2,function(i) all(!is.na(i))))
  out3 <- WPSW(x, t(post_mu), t(post_beta), force = 1, p = 2,
               direction = c("forward"), 
               method=c("projection"), model.size = 3)
  which(apply(out3$beta,2,function(i) all(!is.na(i))))
  
#### Annealing ####
  # Rprof("sa.out", line.profiling=TRUE)
  # debugonce(WPSA)
 sv <-  WPSA(X=x, Y=t(post_mu), t(post_beta), 
                   force = 1, p = 2, model.size = 2:p,
                   # groups = NULL,
                   iter=10, temps = 1000,
                   max.time = 100,
                   proposal = proposal.fun,
                   options = list(method = c("selection.variable"),
                                  energy.distribution = "boltzman",
                                  transport.method = "univariate.approximation",
                                  cooling.schedule = "Geman-Geman"),
             display.progress = TRUE,
  )
 # Rprof()
 # summaryRprof("sa.out", lines ="show")
  sv <-  WPSA(X=x, Y=t(post_mu), t(post_beta), 
              force = 1, p = 2, model.size = 2:p,
              # groups = NULL,
              iter=10, temps = 1000,
              max.time = 6,
              proposal = proposal.fun,
              options = list(method = c("selection.variable"),
                             energy.distribution = "boltzman",
                             transport.method = "univariate.approximation",
                             cooling.schedule = "Geman-Geman"),
              display.progress = TRUE,
  )
 sv <-  WPSA(X=x, Y=t(post_mu), t(post_beta), 
             force = 1, p = 2, model.size = 2:p,
             # groups = NULL,
             iter=50, temps = 50,
             max.time = 30,
             proposal = proposal.fun,
             options = list(method = c("selection.variable"),
                            energy.distribution = "boltzman",
                            transport.method = "univariate.approximation",
                            cooling.schedule = "exponential"),
             display.progress = TRUE
 )
 
   #optimal idx should be 1,9,10
  sv <-  WPSA(X=x, Y=t(post_mu), t(post_beta), 
              force = 1, p = 2, model.size = 3,
              # groups = NULL,
              iter=30, temps = 1000,
              max.time = 30,
              proposal = proposal.fun,
              options = list(method = c("selection.variable"),
                             energy.distribution = "boltzman",
                             transport.method = "univariate.approximation.pwr",
                             cooling.schedule = "Geman-Geman"),
              display.progress = TRUE
  )
  sv$optimal$index
  sv <-  WPSA(X=x, Y=t(post_mu), t(post_beta), 
              force = 1, p = 2, model.size = 3,
              # groups = NULL,
              iter=30, temps = 1000,
              max.time = 30,
              proposal = proposal.fun,
              options = list(method = c("selection.variable"),
                             energy.distribution = "boltzman",
                             transport.method = "exact",
                             cooling.schedule = "Geman-Geman"),
              display.progress = TRUE
  )
  sv$optimal$index
  sc <-  WPSA(X=x, Y=t(post_mu), t(post_beta), 
              force = 1, p = 2, model.size = 3,
              # groups = NULL,
              iter=3, temps = 10,
              max.time = 30,
              proposal = proposal.fun,
              options = list(method = c("scale"),
                             energy.distribution = "boltzman",
                             transport.method = "univariate.approximation.pwr",
                             cooling.schedule = "Geman-Geman")
  )
  sv$optimal$index
  # Rprof("sapr.out", line.profiling=TRUE)
  pr <-  WPSA(X=x, Y=t(post_mu), t(post_beta), 
                force = 1, p = 2, model.size = 5,
                # groups = NULL,
                iter=3, temps = 10,
              max.time = 30,
                proposal = proposal.fun,
                options = list(method = c("projection"),
                               energy.distribution = "boltzman",
                               transport.method = "univariate.approximation.pwr",
                               cooling.schedule = "Geman-Geman")
  )
  pr$optimal$index #1,7,8,9,10
  # Rprof()
  # summaryRprof("sapr.out", lines ="show")
  
  #test multiple model code
  sv.mult <-  WPSA(X=x, Y=t(post_mu), t(post_beta), 
              force = 1, p = 2, model.size = c(3,4,5,6),
              iter=3, temps = 10,
              max.time = 30,
              proposal = proposal.fun,
              options = list(method = c("selection.variable"),
                             energy.distribution = "boltzman",
                             transport.method = "univariate.approximation.pwr",
                             cooling.schedule = "Geman-Geman"),
              display.progress=FALSE
  )
  sv.mult <-  WPSA(X=x, Y=t(post_mu), t(post_beta), 
                   force = 1, p = 2, model.size = c(3,4,5,6),
                   iter=3, temps = 10,
                   max.time = 30,
                   proposal = proposal.fun,
                   options = list(method = c("selection.variable"),
                                  energy.distribution = "boltzman",
                                  transport.method = "univariate.approximation.pwr",
                                  cooling.schedule = "Geman-Geman"),
                   display.progress=TRUE
  )
  
  sv <-  WPSA(X=x, Y=t(post_mu), t(post_beta), 
              force = 1, p = 2, model.size = p,
              # groups = NULL,
              iter=30, temps = 1000,
              max.time = 30,
              proposal = proposal.fun,
              options = list(method = c("selection.variable"),
                             energy.distribution = "boltzman",
                             transport.method = "univariate.approximation.pwr",
                             cooling.schedule = "Geman-Geman"),
              display.progress = TRUE
  )
  sv <-  WPSA(X=x, Y=t(post_mu), t(post_beta), 
              force = 1, p = 2, model.size = c(5,p),
              # groups = NULL,
              iter=3, temps = 3,
              max.time = 30,
              proposal = proposal.fun,
              options = list(method = c("selection.variable"),
                             energy.distribution = "boltzman",
                             transport.method = "univariate.approximation.pwr",
                             cooling.schedule = "Geman-Geman"),
              display.progress = TRUE
  )


#### W2L1 ####
  out <- W2L1(x, post_mu, post_beta,
              penalty = "lasso", method = "selection.variable")
  out <- W2L1(x, post_mu, post_beta,
              penalty = "lasso", method = "projection")
  
#### transport functions ####
  A <- matrix(rnorm(1000*1024),nrow=1024,ncol=1000)
  B <- matrix(rnorm(1000*1024),nrow=1024,ncol=1000)
  at <- t(A)
  bt <- t(B)
  cost <- cost_calc(at,bt,2)
  mass_x <- rep(1/ncol(at),ncol(at))
  mass_y <- rep(1/ncol(bt), ncol(bt))
  test <- transport_plan_given_C(mass_x, mass_y, 2, cost, "exact")
  ttest <- transport::transport(mass_x, mass_y, costm = cost^2, method = "shortsimplex")
  sink <- sinkhorn_distance(mass_x, mass_y, cost, 2, 0.05, 100)
  
  
  # microbenchmark::microbenchmark(transport::transport(mass_x, mass_y, costm = cost^2, method = "shortsimplex"), times = 10)
  # microbenchmark::microbenchmark(transport_plan_given_C(mass_x, mass_y, 2, cost, "exact"), times = 10)
  # microbenchmark::microbenchmark(transport_C_(mass_x, mass_y, cost^2, "shortsimplex"), times = 10)
  

#### wasserstein distances ####
  A <- matrix(rnorm(1000*1024),nrow=1024,ncol=1000)
  B <- matrix(rnorm(1000*1024),nrow=1024,ncol=1000)
  at <- t(A)
  bt <- t(B)
  wass.compare.fun <- function(a,b){
    
    cost <- cost_calc(a,b,2)
    mass <- rep(1/ncol(a), ncol(a))
    return(transport::wasserstein(mass,mass,costm=cost, p = 2, method = "shortsimplex"))
  }
  ttest <- wass.compare.fun(at,bt)
  test <- SparsePosterior::wasserstein(at, bt, 2, 2, "colwise", "exact")
  stest <- SparsePosterior::wasserstein(at,bt, 2, 2, "colwise", "sinkhorn")
  all.equal(ttest,test)
  all.equal(test,stest)
  # microbenchmark::microbenchmark(wass.compare.fun(at,bt), times = 10)
  # microbenchmark::microbenchmark(wasserstein(at, bt, 2, 2, "colwise", "exact"), times = 10)
  # microbenchmark::microbenchmark(wasserstein(at, bt, 2, 2, "colwise", "sinkhorn"), times = 100)
  
  # Rprof("wass.out", line.profiling=TRUE, interval=0.001)
  # for(i in 1:10) wasserstein(at, bt, 2, 2, "colwise", "exact")
  # Rprof(NULL)
  # summaryRprof("wass.out", lines = "both")
#### distance comparing functions ####
