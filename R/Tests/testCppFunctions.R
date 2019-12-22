# Test Cpp functions
rm(list=ls())
require(SparsePosterior)
require(transport)
require(oem)
set.seed(222)

#### Setup Data ####
n <- 100
p <- 10
s <- 1000

x <- matrix(rnorm(p*n), nrow=n, ncol=p)
beta <- (1:10)/10
y <- x %*% beta + rnorm(n)

#posterior
prec <- crossprod(x) + diag(1,p,p)*1
mu_post <- solve(prec, crossprod(x,y))
alpha <- 1 + n/2
beta <- 1 + 0.5 * (crossprod(y) + t(mu_post) %*% prec %*% mu_post )
sigma_post <- 1/rgamma(s, alpha, 1/beta)
theta <- sapply(sigma_post, function(ss) mu_post + t(chol(ss * solve(prec))) %*% matrix(rnorm(p, 0, 1),p,1))

post_mu <- x %*% theta
post_diff <- matrix(c(y),nrow=n,ncol=s) + matrix(rnorm(s*n,0,0.01),nrow=n,ncol=s)
post_vdiff <- matrix(rnorm(n*s),nrow=n,ncol=s)
theta_norm <- rowMeans(theta^2)
pseudo.obs <- 2.7
wt <- n/(pseudo.obs + n)

xtx <- crossprod(x)/n * wt + diag(1,p,p) * (1 - wt)
xty <- crossprod(x, post_mu)/n * wt + theta * (1 - wt)

active.idx <- seq(2,10,2)

#### Check sorting functions ####
Rcpp::sourceCpp("src/testing_src/sort_test.cpp")
X <- matrix(rnorm(1000*100), nrow=100,ncol=1000)
sorts <- sort_check_eig(X)  +1
orders <- apply(X,2,order)
all.equal(c(sorts),c(orders))

Y <- matrix(rnorm(1000*100), nrow=100,ncol=1000)
sortY <- sort_check_rel(X,Y)
orders <- apply(X,2,order)
Y_ord <- Y
for(i in 1:ncol(Y)) {
  sY <- sort(Y[,i])
  for(j in 1:100) Y_ord[orders[j,i],i] <- sY[j]
}
all.equal(c(sortY),c(Y_ord))

sortY2 <- sort_check_rel2(X,Y)
all.equal(c(sortY2),c(Y_ord))

#### Augmented Matrix code ####
  ### new data
  dat <- list(temp=matrix(0, n, p), xtx = matrix(0,p,p), xty = rep_len(0, p),
              mu = rep(0, n), idx_mu = rep(0, n),
              sort_y = rep(0, n))

  for(i in 1:n) {
    dat$temp <- t(theta) * matrix(x[i,,drop=FALSE], s,p, byrow = TRUE)
    dat$mu <- rowSums(dat$temp)
    dat$idx_mu <- order(dat$mu)
    dat$sort_y <- sort(post_mu[i,])
    dat$xtx = dat$xtx + crossprod(dat$temp)
    dat$xty = dat$xty + crossprod(dat$temp[dat$idx_mu,,drop=FALSE], dat$sort_y)
  }
  dat$xtx = dat$xtx/(n*s)
  dat$xty = dat$xty/(n*s)
  out <- sufficientStatistics(x,post_mu, t(theta), FALSE, 0.0, "scale")

  all.equal(out$XtX, dat$xtx)
  all.equal(out$XtY, dat$xty)

  dat$xtx <-  matrix(0,p,p)
  dat$xty <- rep_len(0,p)
  theta_norm <- rowMeans(theta^2)

  for(i in 1:n) {
    dat$temp <- t(theta) * matrix(x[i,,drop=FALSE], s,p, byrow = TRUE)
    dat$mu <- rowSums(dat$temp)
    dat$idx_mu <- order(dat$mu)
    dat$sort_y <- sort(post_mu[i,])
    dat$xtx = dat$xtx + crossprod(dat$temp)/(n*s + pseudo.obs*s)
    dat$xty = dat$xty + crossprod(dat$temp[dat$idx_mu,,drop=FALSE], dat$sort_y)/(n*s + pseudo.obs*s)
  }
  dat$xtx = dat$xtx + diag(theta_norm) * pseudo.obs/(n+pseudo.obs)
  dat$xty = dat$xty + theta_norm * pseudo.obs/(n+pseudo.obs)
  outA <- sufficientStatistics(x,post_mu, t(theta), FALSE, pseudo.obs, "scale")

  all.equal(c(outA$XtX), c(dat$xtx))
  all.equal(c(outA$XtY), c(dat$xty))

  #projection
  #full
  keep.idx <- 1:p
  out1 <- sufficientStatistics(x[,keep.idx], post_mu, t(theta), FALSE, 0 , "projection")
  dat$xty <- crossprod(x[,keep.idx,drop=FALSE], post_mu)/n
  dat$xtx <- crossprod(x[,keep.idx,drop=FALSE])/n
  all.equal(c(out1$XtY), c(dat$xty))
  all.equal(c(out1$XtX), c(dat$xtx))

  out1A <- sufficientStatistics(x[,keep.idx], post_mu, t(theta), FALSE, pseudo.obs , "projection")
  dat$xty <- crossprod(x[,keep.idx,drop=FALSE], post_mu)/(n + pseudo.obs) + theta * pseudo.obs/(s*n + s*pseudo.obs)
  dat$xtx <- crossprod(x[,keep.idx,drop=FALSE])/(n + pseudo.obs) + diag(1,p,p) * pseudo.obs/(s*n + s*pseudo.obs)
  all.equal(c(out1A$XtY), c(dat$xty))
  all.equal(c(out1A$XtX), c(dat$xtx))

  #subset
  keep.idx <- 1:2
  out2 <- sufficientStatistics(x[,keep.idx], post_mu, t(theta), FALSE, 0 , "projection")
  dat$xty <- crossprod(x[,keep.idx,drop=FALSE], post_mu)/n
  dat$xtx <- crossprod(x[,keep.idx,drop=FALSE])/n
  all.equal(c(out2$XtY), c(dat$xty))
  all.equal(c(out2$XtX), c(dat$xtx))

  out2A <- sufficientStatistics(x, post_mu, t(theta), FALSE, pseudo.obs , "projection")
  dat$xty <- crossprod(x, post_mu)/(n + pseudo.obs) + theta * pseudo.obs/(s*n + s*pseudo.obs)
  dat$xtx <- crossprod(x)/(n + pseudo.obs) + diag(1,p,p) * pseudo.obs/(s*n + s*pseudo.obs)
  all.equal(c(out2A$XtY), c(dat$xty))
  all.equal(c(out2A$XtX), c(dat$xtx))

  #location.scale
  # out3 <- sufficientStatistics(x, post_mu, t(theta), FALSE, 0 , "location.scale")
  # m_theta <- matrix(rowMeans(theta),p,s)
  # c_theta <- theta - m_theta
  # dat$xty <- rep(0,2*p)
  # dat$xtx <- matrix(0,2*p,2*p)
  # 
  # for(i in 1:n) {
  #   dat$temp <- t(theta) * matrix(x[i,,drop=FALSE], s,p, byrow = TRUE)
  #   dat$temp2 <- t(theta) * matrix(x[i,,drop=FALSE], s,p, byrow = TRUE)
  #   dat$temp3 <-  cbind(t(c_theta)   * matrix(x[i,,drop=FALSE], s,p, byrow = TRUE),
  #                       t(m_theta)   * matrix(x[i,,drop=FALSE], s,p, byrow = TRUE))
  #   dat$mu <- rowSums(dat$temp2)
  #   dat$idx_mu <- order(dat$mu)
  #   dat$sort_y <- sort(post_mu[i,])
  #   dat$xty = dat$xty + crossprod(dat$temp3[dat$idx_mu,,drop=FALSE], dat$sort_y)/(n*s)
  #   dat$xtx = dat$xtx + crossprod(dat$temp3)/(n*s)
  # }
  # all.equal(c(out3$XtY), c(dat$xty))
  # all.equal(c(out3$XtX), c(dat$xtx))
  # 
  # out3A <- sufficientStatistics(x, post_mu, t(theta), FALSE, pseudo.obs , "location.scale")
  # m_theta <- matrix(rowMeans(theta),p,s)
  # c_theta <- theta - m_theta
  # dat$xty <- rep(0,2*p)
  # dat$xtx <- matrix(0,2*p,2*p)
  # full_theta <- rbind(theta,theta)
  # demeaned_theta <- rbind(c_theta, m_theta)
  # theta_norm_locscale <- rowMeans(full_theta * demeaned_theta)

  # testC <- as.matrix(read.table("c_theta.txt"))
  # all.equal(c(c_theta),c(testC))
  # testM <- as.matrix(read.table("m_theta.txt"))
  # all.equal(c(m_theta),c(testM))
  # testT <- as.matrix(read.table("temp.txt"))
  # all.equal(c(testT),c(t(cbind(t(c_theta) * matrix(x[1,,drop=FALSE], s,p, byrow = TRUE),
  #                              t(m_theta) * matrix(x[1,,drop=FALSE], s,p, byrow = TRUE)))))
  # testTheta <- as.matrix(read.table("theta.txt"))f
  # all.equal(c(theta),c(testTheta))
  # testtnorm <- as.matrix(read.table("theta_norm.txt"))
  # all.equal(c(theta_norm* pseudo.obs/(n + pseudo.obs)),c(testtnorm))

  # for(i in 1:n) {
  #   dat$temp <- t(theta) * matrix(x[i,,drop=FALSE], s,p, byrow = TRUE)
  #   dat$temp2 <- t(theta) * matrix(x[i,,drop=FALSE], s,p, byrow = TRUE)
  #   dat$temp3 <-  cbind(t(c_theta) * matrix(x[i,,drop=FALSE], s,p, byrow = TRUE),
  #                       t(m_theta) * matrix(x[i,,drop=FALSE], s,p, byrow = TRUE))
  #   dat$mu <- rowSums(dat$temp2)
  #   dat$idx_mu <- order(dat$mu)
  #   dat$sort_y <- sort(post_mu[i,])
  #   dat$xty = dat$xty + crossprod(dat$temp3[dat$idx_mu,,drop=FALSE], dat$sort_y)/(n*s + pseudo.obs * s)
  #   temp <- svd(dat$temp3)
  #   dat$xtx = dat$xtx + crossprod(dat$temp3)/(n*s + pseudo.obs * s)
  # }
  # # testxtx1 <- as.matrix(read.table("xtx1.txt"))
  # # all.equal(c(testxtx1),c(dat$xtx))
  # # testxtx2 <- as.matrix(read.table("xtx2.txt"))
  # dat$xtx <- dat$xtx + diag(theta_norm_locscale) * pseudo.obs/(n + pseudo.obs)
  # # all.equal(c(testxtx2), c(dat$xtx + diag(theta_norm) * pseudo.obs/(n + pseudo.obs)))
  # # all.equal(c(testxtx1 + diag(theta_norm) * pseudo.obs/(n + pseudo.obs)), c(dat$xtx + diag(theta_norm) * pseudo.obs/(n + pseudo.obs)))
  # dat$xty <- dat$xty + theta_norm_locscale * pseudo.obs/(n + pseudo.obs)
  # all.equal( c(out3A$XtY), c(dat$xty) )
  # all.equal( c(out3A$XtX), c(dat$xtx ))
  # all.equal(c(testxtx2), c(out3A$XtX))

  ## XtY Update
  #scale method
  out4 <- xtyUpdate(x, post_mu,t(theta),c(1, rep(0,p-1)), 0 , "scale")
  dat$xty <- rep(0,p)
  for(i in 1:n) {
    dat$temp <- t(theta) * matrix(x[i,,drop=FALSE], s,p, byrow = TRUE)
    dat$temp2 <- t(theta) * matrix(c(1, rep(0,p-1)), s,p, byrow=TRUE)  * matrix(x[i,,drop=FALSE], s,p, byrow = TRUE)
    dat$mu <- rowSums(dat$temp2)
    dat$idx_mu <- order(dat$mu)
    dat$sort_y <- sort(post_mu[i,])
    dat$xty = dat$xty + crossprod(dat$temp[dat$idx_mu,,drop=FALSE], dat$sort_y)
  }
  dat$xty = dat$xty/(n*s)
  all.equal(c(out4), c(dat$xty))

  out4A <- xtyUpdate(x, post_mu,t(theta),c(1, rep(0,p-1)), pseudo.obs , "scale")
  dat$xty <- rep(0,p)
  for(i in 1:n) {
    dat$temp <- t(theta) * matrix(x[i,,drop=FALSE], s,p, byrow = TRUE)
    dat$temp2 <- t(theta) * matrix(c(1, rep(0,p-1)), s,p, byrow=TRUE)  * matrix(x[i,,drop=FALSE], s,p, byrow = TRUE)
    dat$mu <- rowSums(dat$temp2)
    dat$idx_mu <- order(dat$mu)
    dat$sort_y <- sort(post_mu[i,])
    dat$xty = dat$xty + crossprod(dat$temp[dat$idx_mu,,drop=FALSE], dat$sort_y)/(s*pseudo.obs + n*s)
  }
  orders <- transport_(t(theta), t(theta) * matrix(c(1, rep(0,p-1)), s,p, byrow=TRUE))
  temp_theta_norm <- rowMeans( theta * theta[,orders])
  # testnorm <- as.matrix(read.table("theta_norm.txt"))
  # all.equal(c(testnorm), c(temp_theta_norm *  pseudo.obs/(n + pseudo.obs)))
  dat$xty = dat$xty + temp_theta_norm * pseudo.obs/(n + pseudo.obs)
  all.equal(c(out4A), c(dat$xty))

  out4B <- xtyUpdate(x, post_mu,t(theta),rep(1,p), pseudo.obs , "scale")
  dat$xty <- rep(0,p)
  for(i in 1:n) {
    dat$temp <- t(theta) * matrix(x[i,,drop=FALSE], s,p, byrow = TRUE)
    dat$temp2 <- t(theta) * matrix(rep(1,p), s,p, byrow=TRUE)  * matrix(x[i,,drop=FALSE], s,p, byrow = TRUE)
    dat$mu <- rowSums(dat$temp2)
    dat$idx_mu <- order(dat$mu)
    dat$sort_y <- sort(post_mu[i,])
    dat$xty = dat$xty + crossprod(dat$temp[dat$idx_mu,,drop=FALSE], dat$sort_y)/(pseudo.obs*s +n*s)
  }

  dat$xty = dat$xty +theta_norm * pseudo.obs/(n + pseudo.obs )
  all.equal(c(out4B), c(dat$xty))

#projection method
  #check full data first
  keep.idx <- 1:p
  perp <- coef(lm(post_mu ~ x[,keep.idx] + 0))
  out5 <- xtyUpdate(x[,keep.idx], post_mu, t(theta)[,keep.idx], perp, 0 , "projection")
  dat$xty <- crossprod(x[,keep.idx,drop=FALSE], post_mu)/n
  all.equal(c(out5), c(dat$xty))

  keep.idx <- 1:p
  perp <- coef(lm(post_mu ~ x[,keep.idx] + 0))
  cond_mu <- x[,keep.idx] %*% perp
  orders <- transport_(t(cond_mu),t(post_mu))
  out5A <- xtyUpdate(x[,keep.idx], post_mu, t(theta)[,keep.idx], perp, pseudo.obs, "projection")
  dat$xty <- crossprod(x[,keep.idx,drop=FALSE], post_mu[,orders,drop=FALSE])/(n + pseudo.obs) + theta * pseudo.obs/(pseudo.obs + n)/s
  all.equal(c(out5A), c(dat$xty))

  #subset
  keep.idx <- 1:2
  perp <- coef(lm(post_mu ~ x[,keep.idx] + 0))
  out6 <- xtyUpdate(x[,keep.idx], post_mu, t(theta), perp, 0.0 , "projection")
  dat$xty <- crossprod(x[,keep.idx,drop=FALSE], post_mu)/n
  all.equal(c(out6), c(dat$xty))

  perp <- coef(lm(post_mu ~ x[,keep.idx] + 0))
  cond_mu <- x[,keep.idx] %*% perp
  orders <- transport_(t(cond_mu),t(post_mu))
  perp_dist <- matrix(0,p,s)
  perp_dist[keep.idx,] <- perp
  order_theta <- transport_(t(theta),t(perp_dist))
  out6A <- xtyUpdate(x, post_mu, t(theta), perp_dist, pseudo.obs, "projection")
  dat$xty <- crossprod(x[,,drop=FALSE], post_mu[,orders,drop=FALSE])/(n + pseudo.obs) + theta[,order_theta] * pseudo.obs/(pseudo.obs + n)/s
  all.equal(c(out6A), c(dat$xty))

# location scale
  # m_theta <- matrix(rowMeans(theta),p,s)
  # c_theta <- theta - m_theta
  # out7 <- xtyUpdate(x, post_mu,t(theta),rep(c(1, rep(0,p-1)),2), 0 , "location.scale")
  # dat$xty <- rep(0,2*p)
  # 
  # for(i in 1:n) {
  #   dat$temp <- t(theta) * matrix(x[i,,drop=FALSE], s,p, byrow = TRUE)
  #   dat$temp2 <- t(theta) * matrix(c(1, rep(0,p-1)), s,p, byrow=TRUE)  * matrix(x[i,,drop=FALSE], s,p, byrow = TRUE)
  #   dat$temp3 <-  cbind(t(c_theta)   * matrix(x[i,,drop=FALSE], s,p, byrow = TRUE),
  #                       t(m_theta)   * matrix(x[i,,drop=FALSE], s,p, byrow = TRUE))
  #   dat$mu <- rowSums(dat$temp2)
  #   dat$idx_mu <- order(dat$mu)
  #   dat$sort_y <- sort(post_mu[i,])
  #   dat$xty = dat$xty + crossprod(dat$temp3[dat$idx_mu,,drop=FALSE], dat$sort_y)/(n*s)
  # }
  # all.equal(c(out7), c(dat$xty))
  # 
  # out7A <- xtyUpdate(x, post_mu,t(theta),rep(c(1, rep(0,p-1)),2), pseudo.obs , "location.scale")
  # # res_thet <- as.matrix(read.table("result.txt"))
  # # all.equal(c(res_thet), c(t(cbind(t(c_theta),t(m_theta)) %*% diag(rep(c(1, rep(0,p-1)),2)))))
  # dat$xty <- rep(0,2*p)
  # demeaned_theta <- t(cbind(t(c_theta),t(m_theta)))
  # 
  # for(i in 1:n) {
  #   dat$temp <- t(theta) * matrix(x[i,,drop=FALSE], s,p, byrow = TRUE)
  #   dat$temp2 <- t(theta) * matrix(c(1, rep(0,p-1)), s,p, byrow=TRUE)  * matrix(x[i,,drop=FALSE], s,p, byrow = TRUE)
  #   dat$temp3 <-  cbind(t(c_theta)   * matrix(x[i,,drop=FALSE], s,p, byrow = TRUE),
  #                       t(m_theta)   * matrix(x[i,,drop=FALSE], s,p, byrow = TRUE))
  #   dat$mu <- rowSums(dat$temp2)
  #   dat$idx_mu <- order(dat$mu)
  #   dat$sort_y <- sort(post_mu[i,])
  #   dat$xty = dat$xty + crossprod(dat$temp3[dat$idx_mu,,drop=FALSE], dat$sort_y)/(n*s + pseudo.obs*s)
  # }
  # orders <- transport_(t(theta), t(diag(c(1, rep(0,p-1))) %*% theta))
  # theta_norm_temp <- rowMeans(demeaned_theta[,orders] * demeaned_theta)
  # dat$xty <- dat$xty + theta_norm_temp * pseudo.obs/(n + pseudo.obs)
  # all.equal(c(out7A), c(dat$xty))
  
  
  ### different mu's!
  #scale method
  out8 <- xtyUpdate(x, post_diff,t(theta),c(1, rep(0,p-1)), 0 , "scale")
  dat$xty <- rep(0,p)
  for(i in 1:n) {
    dat$temp <- t(theta) * matrix(x[i,,drop=FALSE], s,p, byrow = TRUE)
    dat$temp2 <- t(theta) * matrix(c(1, rep(0,p-1)), s,p, byrow=TRUE)  * matrix(x[i,,drop=FALSE], s,p, byrow = TRUE)
    dat$mu <- rowSums(dat$temp2)
    dat$idx_mu <- order(dat$mu)
    dat$sort_y <- sort(post_diff[i,])
    dat$xty = dat$xty + crossprod(dat$temp[dat$idx_mu,,drop=FALSE], dat$sort_y)
  }
  dat$xty = dat$xty/(n*s)
  all.equal(c(out8), c(dat$xty))
  
  beta <- c(1, 1,1,1,rep(0,p-4))
  out9 <- xtyUpdate(x, post_diff,t(theta),beta, 0 , "scale")
  dat$xty <- rep(0,p)
  for(i in 1:n) {
    dat$temp <- t(theta) * matrix(x[i,,drop=FALSE], s,p, byrow = TRUE)
    dat$temp2 <- t(theta) * matrix(beta, s,p, byrow=TRUE)  * matrix(x[i,,drop=FALSE], s,p, byrow = TRUE)
    dat$mu <- rowSums(dat$temp2)
    dat$idx_mu <- order(dat$mu)
    dat$sort_y <- sort(post_diff[i,])
    dat$xty = dat$xty + crossprod(dat$temp[dat$idx_mu,,drop=FALSE], dat$sort_y)/(n*s)
  }
  # dat$xty = dat$xty/(n*s)
  all.equal(c(out9), c(dat$xty))
  
  beta <- rep(1,p)
  out10 <- xtyUpdate(x, post_diff,t(theta), beta, 0 , "scale")
  dat$xty <- rep(0,p)
  for(i in 1:n) {
    dat$temp <- t(theta) * matrix(x[i,,drop=FALSE], s,p, byrow = TRUE)
    dat$temp2 <- t(theta) * matrix(beta, s,p, byrow=TRUE)  * matrix(x[i,,drop=FALSE], s,p, byrow = TRUE)
    dat$mu <- rowSums(dat$temp2)
    dat$idx_mu <- order(dat$mu)
    dat$sort_y <- sort(post_diff[i,])
    dat$xty = dat$xty + crossprod(dat$temp[dat$idx_mu,,drop=FALSE], dat$sort_y)
  }
  dat$xty = dat$xty/(n*s)
  all.equal(c(out10), c(dat$xty))
  
  beta <- rep(1,p)
  out11 <- xtyUpdate(x, post_diff,t(theta), beta, 0 , "selection.variable")
  dat$xty <- rep(0,p)
  for(i in 1:n) {
    dat$temp <- t(theta) * matrix(x[i,,drop=FALSE], s,p, byrow = TRUE)
    dat$temp2 <- t(theta) * matrix(beta, s,p, byrow=TRUE)  * matrix(x[i,,drop=FALSE], s,p, byrow = TRUE)
    dat$mu <- rowSums(dat$temp2)
    dat$idx_mu <- order(dat$mu)
    dat$sort_y <- sort(post_diff[i,])
    dat$xty = dat$xty + crossprod(dat$temp[dat$idx_mu,,drop=FALSE], dat$sort_y)
  }
  dat$xty = dat$xty/(n*s)
  all.equal(c(out11), c(dat$xty))
  
 ## XtXY weighted method
  dat <- list(temp=matrix(0, n, p), xtx = matrix(0,p,p), xty = rep_len(0, p),
              mu = rep(0, n), idx_mu = rep(0, n),
              sort_y = rep(0, n))
  
  result <- c(1,rep(0,p-1))
  # weight <- lapply(1:n, function(i) chol(diag(1,s,s)))
  weight <- chol(diag(1,s,s))
  pwr <- 3
  
  for(i in 1:n) {
    dat$temp <- t(theta) * matrix(x[i,,drop=FALSE], s,p, byrow = TRUE)
    dat$mu <- rowSums(dat$temp[,which(result==1),drop=FALSE])
    dat$idx_mu <- order(dat$mu)
    dat$sort_y <- sort(post_mu[i,])
    dat$weight <- sqrt(abs(dat$sort_y - dat$mu[dat$idx_mu,drop=FALSE])^(pwr - 2))
    dat$xtx = dat$xtx + crossprod(dat$weight * dat$temp[dat$idx_mu,,drop=FALSE])
    dat$xty = dat$xty + crossprod(dat$weight * dat$temp[dat$idx_mu,,drop=FALSE], dat$weight * dat$sort_y)
  }
  dat$xtx = dat$xtx/(n*s)
  dat$xty = dat$xty/(n*s)
  out12 <- xtxyUpdate(x,post_mu, t(theta), result, 0.0, 0, pwr, "selection.variable")
  all.equal(out12$XtY, dat$xty)
  all.equal(out12$XtX, dat$xtx)

#### sorting LM code ####
suffStat_star <- sufficientStatistics(x, post_mu, t(theta), TRUE,pseudo.obs,"location.scale")

#sorting projection operator
out <- sorting_lm_proj_approx(x, post_mu, theta, 5, pseudo.obs)
all.equal(out, theta)
transport::wasserstein(transport::pp(t(out)),
                       transport::pp(t(theta)), 2,
                       method="shortsimplex")

out2 <- sorting_lm_proj_approx( x[,active.idx,drop=FALSE], post_mu, theta[active.idx,,drop=FALSE], 5, pseudo.obs )
out2_compare <- calc.gamma(xtx, xty, active.idx, "approx.projection",
                           x, theta, theta_norm, pseudo.obs, post_mu, 5)
transport::wasserstein(transport::pp(t(out2)),
                       transport::pp(t(out2_compare)), 2,
                       method="shortsimplex") #fairly close but cpp doesn't have sorting on theta yet

#### transport test ####
  A <- matrix(rnorm(1000*1024),nrow=1024,ncol=1000)
  B <- matrix(rnorm(1000*1024),nrow=1024,ncol=1000)
  # dist_mat <- as.matrix(dist(rbind(A,B)))[1:1024, 1025:2048]
  # dist_mat <- dist_mat^2
  # dist_check <- matrix(0,1024,1024)
  at <- t(A)
  bt <- t(B)
  # for(i in 1:1024) for(j in 1:1024) dist_check[i,j] <- sum((at[,i] - bt[,j])^2)
  # all.equal(c(dist_mat), c(dist_check))
  indexes <- transport_(at, bt, 2.0, 2.0, "shortsimplex")
  # debugonce(transport::transport.pp)
  index_trans <- transport::transport(transport::pp(A),transport::pp(B),p=2, method = "shortsimplex")
  all.equal(indexes$from, index_trans[["from"]])
  all.equal(indexes$to, index_trans[["to"]])
  all.equal(indexes$mass, index_trans[["mass"]]/1024)
  
  mass_a <- rep(1/ncol(at), ncol(at))
  mass_b <- rep(1/ncol(bt), ncol(bt))
  costm <- cost_calc(at,bt,2)
  indexes2 <- transport_C_(mass_a, mass_b, costm^2, "shortsimplex")
  index_trans2 <- transport::transport.default(mass_a, mass_b, costm^2, method = "shortsimplex")
  check_sink <- sinkhorn_(mass_a, mass_b, costm^2, 0.05*median(costm^2), 100)
  sum(check_sink$transportmatrix * costm^2)
  all.equal(indexes2$from, index_trans2[["from"]])
  all.equal(indexes2$to, index_trans2[["to"]])
  all.equal(indexes2$mass, index_trans2[["mass"]])
  all.equal(indexes2$from, index_trans[["from"]])
  all.equal(indexes2$to, index_trans[["to"]])
  all.equal(indexes2$mass, index_trans[["mass"]]/1024)
  
  C <- t(A[1:100,,drop = FALSE])
  D <- t(B[1:2,,drop = FALSE])
  
  costm <- cost_calc(C,D,2.0)
  mass_c <- rep(1/ncol(C), ncol(C))
  mass_d <- rep(1/ncol(D), ncol(D))
  
  trans_sp <- transport_C_(mass_c, mass_d, costm^2, method = "shortsimplex")
  # debugonce(transport::transport.default)
  trans_t <- transport::transport.default(a=mass_c, b=mass_d, costm=costm^2, method = "shortsimplex")
  all.equal(trans_sp$from, trans_t$from)
  all.equal(trans_sp$to, trans_t$to)
  all.equal(trans_sp$mass, trans_t$mass)
  # microbenchmark::microbenchmark(transport::transport.default(a=mass_c, b=mass_d, costm=costm^2, method = "shortsimplex"), unit="us")
  # microbenchmark::microbenchmark(transport_C_(mass_c, mass_d, costm^2, method = "shortsimplex"), unit = "us")
  
  trans_t <- transport::transport.default(a=mass_d, b=mass_c, costm=t(costm^2), method = "shortsimplex")
  trans_sp <- transport_C_(mass_d, mass_c, t(costm^2), method = "shortsimplex")
  all.equal(trans_sp$from, trans_t$from)
  all.equal(trans_sp$to, trans_t$to)
  all.equal(trans_sp$mass, trans_t$mass)
  
  
#### Wasserstein test ####
  A <- matrix(rnorm(1000*1024),nrow=1024,ncol=1000)
  B <- matrix(rnorm(1000*1024),nrow=1024,ncol=1000)
  at <- t(A)
  bt <- t(B)
  cost <- cost_calc(at,bt,2)
  mass_a <- rep(1/ncol(at),ncol(at))
  mass_b <- rep(1/ncol(bt),ncol(bt))
  
  tplan <- transport_plan_given_C(mass_a, mass_b, 2, cost, "exact")

  loss <- wasserstein_(tplan$mass, cost, p = 2, tplan$from, tplan$to)
  sink <- sinkhorn_(mass_a, mass_b, cost^2, 0.05*median(cost^2), 100)
  sink$distances^(1/2)
  sum(sink$transportmatrix * cost^2)^(1/2)
  # debugonce(transport::wasserstein)
  loss_t <- transport::wasserstein(transport::pp(A),transport::pp(B),p=2, method = "shortsimplex")
  loss_t_def <- transport::wasserstein(mass_a,mass_b, tplan = data.frame(tplan), costm = cost, p=2, method = "shortsimplex")
  all.equal(loss, loss_t)
  all.equal(loss, loss_t_def)
  # microbenchmark::microbenchmark(transport::wasserstein(mass_a,mass_b, tplan = data.frame(tplan), costm = cost, p=2, method = "shortsimplex"), unit="us")
  # microbenchmark::microbenchmark(wasserstein_(tplan$mass, cost, p = 2, tplan$from, tplan$to), unit = "us")
  # microbenchmark::microbenchmark(sinkhorn_(mass_a, mass_b, cost^2, 0.05*median(cost^2), 100), unit="ms")
  
  C <- t(A[1:100,,drop = FALSE])
  D <- t(B[1:2,,drop = FALSE])
  
  cost2 <- cost_calc(C,D,2)
  mass_c <- rep(1/ncol(C),ncol(C))
  mass_d <- rep(1/ncol(D),ncol(D))
  tplan2 <- transport_plan_given_C(mass_c, mass_d, 2, cost2, "exact")
  loss <- wasserstein_(tplan2$mass, cost2, p = 2, tplan2$from, tplan2$to)
  loss_t_def <- transport::wasserstein(mass_c,mass_d, tplan = data.frame(tplan2), costm = cost2, p=2, method = "shortsimplex")
  all.equal(loss, loss_t_def)
  
#### W2 Test ####
  #good luck
  lambda <- 0
  nlambda <- 100
  lambda.min.ratio <- 1e-10
  gamma <- 1
  penalty.factor <- 1/rowMeans(theta^2)
  penalty.factor.null <- rep(1,p)
  suffstat <- sufficientStatistics(x,post_mu, t(theta), TRUE, 0.0, "scale")
  suffstatd <- sufficientStatistics(x,post_diff, t(theta), FALSE, 0.0, "scale")
  suffstatdd <- sufficientStatistics(x,post_vdiff, t(theta), FALSE, 0.0, "scale")
  
  #compare to oem
  check <- oem.xtx(xtx=suffstat$XtX, xty=suffstat$XtY, family="gaussian",penalty="ols", lambda=0,maxit=10000, scale.factor = sqrt(diag(suffstat$XtX)))
  w2 <- W2L1(X=x, Y=NULL, 
             theta=theta, penalty="ols",
             nlambda = 1, lambda.min.ratio = lambda.min.ratio,
             infimum.maxit=1, maxit = 1e3, gamma = gamma, 
             display.progress = FALSE, lambda=0,
             method="scale")
  cbind(check$beta[[1]],w2$beta)
  all.equal(c(check$beta[[1]]), c(w2$beta)) # not quite equal but ok
  
  check2 <- oem.xtx(xtx=suffstatd$XtX, xty=suffstatd$XtY, family="gaussian",penalty="ols", lambda=0,maxit=10000, scale.factor = sqrt(diag(suffstatd$XtX)))
  w22 <- W2L1(X=x, Y=post_diff, 
             theta=theta, penalty="ols",
             nlambda = 1, lambda.min.ratio = lambda.min.ratio,
             infimum.maxit=1, maxit = 1e3, gamma = gamma, 
             display.progress = FALSE, lambda=0,
             method="scale")
  cbind(check2$beta[[1]],w22$beta)
  all.equal(c(check2$beta[[1]]), c(w22$beta)) #ok but worse
  
  check3 <- oem.xtx(xtx=suffstatdd$XtX, xty=suffstatdd$XtY, family="gaussian",penalty="ols", lambda=0,maxit=10000, scale.factor = sqrt(diag(suffstatdd$XtX)))
  w23 <- W2L1(X=x, Y=post_vdiff, 
              theta=theta, penalty="ols",
              nlambda = 1, lambda.min.ratio = lambda.min.ratio,
              infimum.maxit=1, maxit = 1e3, gamma = gamma, 
              display.progress = FALSE, lambda=0,
              method="scale")
  cbind(check3$beta[[1]],w23$beta)
  all.equal(c(check3$beta[[1]]), c(w23$beta)) #worse
  all.equal(c(check3$d), c(w23$d)) #eigenvals good
  c(check3$niter[[1]], w23$niter[1,1])
  all.equal(c(w23$xtx), c(suffstatdd$XtX))
  all.equal(c(w23$xty), c(suffstatdd$XtY)) #xty same!
  
  
  selection <- W2L1(X=x, Y=NULL, 
       theta=theta, penalty="selection.lasso",
       nlambda = nlambda, lambda.min.ratio = lambda.min.ratio,
       infimum.maxit=1e4, maxit = 1e3, gamma = gamma,
       display.progress = FALSE, lambda=lambda,
       penalty.factor = penalty.factor.null, method="selection.variable")
  all(selection$beta==1)
  selection <- W2L1(X=x, Y=NULL, 
                    theta=theta, penalty="selection.lasso",
                    nlambda = 10, lambda.min.ratio = lambda.min.ratio,
                    infimum.maxit=1e4, maxit = 1e3, gamma = gamma,
                    display.progress = FALSE,
                    penalty.factor = penalty.factor.null, method="selection.variable")
  all(selection$beta %in% c(0,1))
  
  selection <- W2L1(X=x, Y=NULL, 
                    theta=theta, penalty="selection.lasso",
                    nlambda = 10, lambda.min.ratio = lambda.min.ratio,
                    infimum.maxit=1e4, maxit = 1e3, gamma = gamma,
                    display.progress = FALSE,
                    penalty.factor = penalty.factor.null, method="selection.variable")
  all(selection$beta %in% c(0,1))
  
  selection <- W2L1(X=x, Y=NULL, 
                    theta=theta, penalty="selection.lasso",
                    nlambda = 10, lambda.min.ratio = lambda.min.ratio,
                    infimum.maxit=1e4, maxit = 1e3, gamma = gamma,
                    display.progress = FALSE,
                    penalty.factor = penalty.factor.null, method="selection.variable", model.size = 3)
  all(selection$beta %in% c(0,1))
  selection <- W2L1(X=x, Y=NULL, 
                    theta=theta, penalty="selection.lasso",
                    nlambda = 100, lambda.min.ratio = 0.1,
                    infimum.maxit=1e4, maxit = 1e3, gamma = gamma,
                    display.progress = FALSE,
                    penalty.factor = penalty.factor.null, method="selection.variable", model.size = 3)
  all(selection$beta %in% c(0,1))
  
  selection <- W2L1(X=x, Y=NULL, 
                    theta=theta, penalty="selection.lasso",
                    nlambda = nlambda, lambda.min.ratio = lambda.min.ratio,
                    infimum.maxit=1e4, maxit = 1e3, gamma = gamma,
                    display.progress = FALSE,
                    penalty.factor = rep(1,p), method="selection.variable")
  all(selection$beta %in% c(0,1))
  any(selection$beta %in% 0)
  any(selection$beta %in% 1)
  
  scale <- W2L1(X=x, Y=NULL, 
             theta=theta, penalty="mcp",
             nlambda = nlambda, lambda.min.ratio = lambda.min.ratio,
             infimum.maxit=1e4, maxit = 1e3, gamma = gamma,
             display.progress = FALSE, lambda=lambda,
             penalty.factor = penalty.factor.null, method="scale")
  all.equal(c(scale$beta), rep(1,p)) #should be pretty close
  scale <- W2L1(X=x, Y=NULL, 
                    theta=theta, penalty="mcp",
                    nlambda = nlambda, lambda.min.ratio = lambda.min.ratio,
                    infimum.maxit=1e4, maxit = 1e3, gamma = gamma,
                    display.progress = FALSE,
                    penalty.factor = penalty.factor, method="scale")
  print(scale$beta)
  
  projection <- W2L1(X=x, Y=NULL, 
                     theta=theta, penalty="ols",
                     nlambda = nlambda, lambda.min.ratio = lambda.min.ratio,
                     infimum.maxit=3, maxit = 1e3, gamma = gamma,
                     display.progress = FALSE, lambda=lambda,
                     penalty.factor = penalty.factor.null, method="projection")
  all.equal(c(projection$beta), c(theta)) #should be pretty close
  
  projection <- W2L1(X=x, Y=NULL, 
                     theta=theta, penalty="mcp",
                     nlambda = nlambda, lambda.min.ratio = lambda.min.ratio,
                     infimum.maxit=3, maxit = 1e3, gamma = gamma,
                     display.progress = FALSE, lambda=lambda,
                     penalty.factor = penalty.factor.null, method="projection")
  all.equal(c(projection$beta), c(theta)) #should be pretty close
  
  projection <- W2L1(X=x, Y=NULL, 
                     theta=theta, penalty="mcp",
                     nlambda = nlambda, lambda.min.ratio = lambda.min.ratio,
                     infimum.maxit=3, maxit = 1e3, gamma = gamma, 
                     display.progress = FALSE, 
                     penalty.factor = penalty.factor, method="projection")
  projection$beta[1:10,]
  
  projection <- W2L1(X=x, Y=NULL, 
                     theta=theta, penalty="scad",
                     nlambda = nlambda, lambda.min.ratio = lambda.min.ratio,
                     infimum.maxit=3, maxit = 1e3, gamma = gamma,
                     display.progress = FALSE, lambda=lambda,
                     penalty.factor = penalty.factor, method="projection")
  all.equal(c(projection$beta), c(theta)) #should be pretty close
  
  projection <- W2L1(X=x, Y=NULL, 
                     theta=theta, penalty="scad",
                     nlambda = nlambda, lambda.min.ratio = lambda.min.ratio,
                     infimum.maxit=3, maxit = 1e3, gamma = gamma, 
                     display.progress = FALSE, 
                     penalty.factor = penalty.factor, method="projection")
  projection$beta[1:10,]
  
  projection <- W2L1(X=x, Y=NULL, 
                     theta=theta, penalty="lasso",
                     nlambda = nlambda, lambda.min.ratio = lambda.min.ratio,
                     infimum.maxit=3, maxit = 1e3, gamma = gamma,
                     display.progress = FALSE, lambda=lambda,
                     penalty.factor = penalty.factor, method="projection")
  all.equal(c(projection$beta), c(theta)) #should be pretty close
  
  projection <- W2L1(X=x, Y=NULL, 
                     theta=theta, penalty="lasso",
                     nlambda = nlambda, lambda.min.ratio = lambda.min.ratio,
                     infimum.maxit=2, maxit = 1e3, gamma = gamma, 
                     display.progress = FALSE, 
                     penalty.factor = penalty.factor, method="projection")
  projection$beta[1:10,]
  
  
  locscale <- W2L1(X=x, Y=NULL, 
                     theta=theta, penalty="mcp",
                     nlambda = nlambda, lambda.min.ratio = lambda.min.ratio,
                     infimum.maxit=2, maxit = 1e3, gamma = gamma,
                     display.progress = FALSE, lambda=lambda,
                     penalty.factor = penalty.factor, method="location.scale")
  all.equal(c(locscale$beta), c(rep(1,2*p))) #should be pretty close
  
  locscale <- W2L1(X=x, Y=NULL, 
                     theta=theta, penalty="mcp",
                     nlambda = nlambda, lambda.min.ratio = lambda.min.ratio,
                     infimum.maxit=2, maxit = 1e3, gamma = gamma, 
                     display.progress = FALSE, 
                     penalty.factor = penalty.factor, method="location.scale")
  locscale$beta[1:(2*p),]
  
  projection <- W2L1(X=x[1,,drop=FALSE], Y=NULL, 
                     theta=theta, penalty="ols",
                     nlambda = nlambda, lambda.min.ratio = lambda.min.ratio,
                     infimum.maxit=1e4, maxit = 1e3, gamma = gamma,
                     display.progress = FALSE, lambda=lambda,
                     penalty.factor = penalty.factor.null, method="projection")
  all.equal(c(projection$beta), c(theta)) #should be pretty close
  
  projection <- W2L1(X=x[,,drop=FALSE], Y=NULL, 
                     theta=theta, penalty="elastic.net",
                     nlambda = nlambda, lambda.min.ratio = lambda.min.ratio,
                     infimum.maxit=2, maxit = 1e3, gamma = gamma, alpha=0.5,
                     display.progress = FALSE, lambda=lambda,
                     penalty.factor = penalty.factor.null, method="projection")
  all.equal(c(projection$beta), c(theta)) #should be pretty close
  
  projection <- W2L1(X=x[1:10,,drop=FALSE], Y=NULL, 
                     theta=theta, penalty="mcp.net",
                     nlambda = nlambda, lambda.min.ratio = lambda.min.ratio,
                     infimum.maxit=2, maxit = 1e3, gamma = gamma, alpha=0.5,
                     display.progress = FALSE, lambda=lambda,
                     penalty.factor = penalty.factor.null, method="projection")
  all.equal(c(projection$beta), c(theta)) #should be pretty close
  projection <- W2L1(X=x[1:10,,drop=FALSE], Y=NULL, 
                     theta=theta, penalty="scad.net",
                     nlambda = nlambda, lambda.min.ratio = lambda.min.ratio,
                     infimum.maxit=2, maxit = 1e3, gamma = gamma, alpha=0.5,
                     display.progress = FALSE, lambda=lambda,
                     penalty.factor = penalty.factor.null, method="projection")
  all.equal(c(projection$beta), c(theta)) #should be pretty close
  