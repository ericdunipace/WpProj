#compare distances
  set.seed(12308947)
  n <- 512
  d <- 51
  p <- 2.0
  set.seed(293897)
  A <- matrix(rnorm(n*d),nrow=d,ncol=n)
  B <- matrix(rnorm(n*d),nrow=d,ncol=n)
  transp.meth <- c("exact","sinkhorn", "greenkhorn","randkhorn","gandkhorn", "hilbert","rank")
  niter = 5e4
  dist <- trans <- vector("list",length(transp.meth))
  names(dist) <- transp.meth
  mass_a <- rep(1/n,n)
  mass_b <- rep(1/n,n)
  cost <- cost_calc(A,B,p)
  
  for(i in seq_along(trans)) {
    trans[[i]] <- transport_plan(A, B,p,p, "colwise", method = transp.meth[[i]], niter = niter)
    dist[[i]] <- wasserstein_(trans[[i]]$tplan$mass, cost, p, trans[[i]]$tplan$from, trans[[i]]$tplan$to)
  }
print(unlist(dist))
