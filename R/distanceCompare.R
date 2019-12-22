distCompare <- function(models, target = list(posterior = NULL, mean = NULL), p = 2, ground_p = 2, 
                         method = "exact", 
                         quantity = c("posterior","mean"),
                         parallel=NULL, transform = function(x){return(x)}, ...) 
{
  method <- match.arg(method, c("mse", transport_options()), 
                      several.ok = TRUE)
  quantity <- match.arg(quantity, several.ok = TRUE)
  
  if ( !is.null(parallel) ) {
      if ( !inherits(parallel, "cluster") ) {
        stop("parallel must be a registered cluster backend or null")
      }
      doParallel::registerDoParallel(parallel)
    
    # display.progress <- FALSE
  } else {
    foreach::registerDoSEQ()
  }
  
  if (inherits(models, "sparse-posterior") ) {
    models <- list(models)
  } else {
    stopifnot(all(sapply(models, inherits, "sparse-posterior")))
  }
  
  dist_df <- dist_mu_df <- nactive <- groups <- plot <- plot_mu <- NULL
  
  nzero <- lapply(models, function(mm) mm$nzero)
  
  data <- set_dist_data(target, models, quantity, method, transform)
  group_names <- data$group_names
  n <- data$n
  d <- data$d
  s <- data$s
  
  if ("posterior" %in% quantity){
    posterior <- data$posterior
    dist_list <- mapply(function(mc, proj){
        dist_fun(mc, mu = posterior$target, p = p, ground_p = ground_p, 
                 method = posterior$method, observation.orientation = posterior$obs.direction,
                 projection = proj, ...) }, mc = posterior$source, proj = posterior$projection)
    
    dist <- unlist(dist_list)
    nactive <- unlist(nzero)
    groups <- mapply(function(x,z){return(rep(x, each=z))}, 
                     x = group_names, z = sapply(dist_list, length))
    
    if(posterior$isMSE)
    {
      dist <- (dist^2)/d
      posterior$method <- "mse"
    }    
    dist_df <- data.frame(dist = dist,
                          nactive = unlist(nactive),
                          groups=factor(unlist(groups)),
                          method = posterior$method)
  }
  
  if ("mean" %in% quantity){
    
    mu <- data$mean
    dist_list_mu <- mapply(function(mc, proj){
      dist_fun(mc, mu = mu$target, p = p, ground_p = ground_p,
               method = mu$method, observation.orientation = mu$obs.direction,
               projection = proj, ...)}, mc = mu$source, proj = mu$projection)
    
    dist_mu <- unlist(dist_list_mu)
    if (is.null(nactive)) nactive <- unlist(nzero)
    
    # if (is.null(groups)) groups <- sapply(group_names, function(g) rep(g, nrow(dist_list_mu)))
    if (is.null(groups)) groups <- mapply(function(x,z){return(rep(x, each=z))}, 
                                          x = group_names, z = sapply(dist_list_mu, length))
    
    if( mu$isMSE )
    {
      dist_mu <- (dist_mu^2)/n
      mu$method <- "mse"
    }   
    
    dist_mu_df <- data.frame(dist = dist_mu,
                             nactive = unlist(nactive),
                             groups=factor(unlist(groups)),
                             method = mu$method)
    
  }
  
  # if (parallel) parallel::stopCluster(cl)
  output <- list(posterior = dist_df, mean = dist_mu_df, p = p)
  class(output) <- c("distcompare","sparse-posterior")
  
  return(output)
}

is.distcompare <- function(x) inherits(x, "distcompare")

dist_fun <- function(mulist, mu, p, ground_p, method, observation.orientation, projection, ...) {
  dist <-
    foreach::foreach(m=mulist, .combine = c) %dopar% {
        wp <- if(projection & method == "exact" & p == 2) {
          denom <- max(ncol(mu), ncol(m))
          sqrt(sum((c(m) - c(mu))^2)/denom)
        } else {
          SparsePosterior::wasserstein(X = m, Y = mu, p = p, 
                                       ground_p = ground_p, 
                                       observation.orientation = observation.orientation, 
                                       method = method, ...)
        }
        return(wp)
      }
  return(dist)
}

set_dist_data <- function(target, models, quantity, method, transform) {
  
  post_targ <- mean_targ <- NULL
  if (length(quantity) > 2) {
    quantity <- sort(unique(quantity), decreasing=TRUE)
    warning("Arbitrarily shortening quantity length. Length should be less than or equal to 2.")
  }
  isMSE <- rep(FALSE, 2)
  
  ### method checks ###
  if (length(quantity) == 1 &  1 < length(method)) {
    if ( !is.null(names(method)) ) {
      method <- method[order(names(method), decreasing=TRUE)]
    } 
    if("posterior" %in% quantity) {
      meth.rm <- 2L
    } else if ("mean" %in% quantity) {
      meth.rm <- 1L
    }
    method[meth.rm] <- NA
    warning("1 == length(quantity) < length(method). Arbitrarily shortening method length to agree with quantity length.")
  } else if (length(method) < length(quantity)) {
    method <- rep(method, length(quantity))
  }
  if (length(method) != 2) {
    method <- rep(unique(method), 2)
  }
  if ( length(method ) > 2 ) {
    method <- method[1:2]
  }
  ### Target checks ###
  if (!is.null(target) ) {
    if (length(quantity) > 1 & (!is.list(target) | (is.list(target) & length(target) ==1) ) ) {
      stop("For more than one quantity, target should be a list of length 2. The list should either be named with slots 'posterior' and 'mean' OR the posterior target should be in the first slot and the mean target in the second.")
    }
    if (!is.list(target)) {
      target <- list(target)
    }
    if(is.null(names(target))) names(target) <- sort(quantity, decreasing = TRUE) 
    if ( any (names(target) != quantity) ) {
      target <- target[sort(quantity, decreasing = TRUE) ]
    }
    if ("posterior" %in% quantity) {
      post_targ <- target$posterior
    }
    if ("mean" %in% quantity) {
      mean_targ <- target$mean
    }
  } else {
    target <- list(posterior = NULL, mean = NULL)
  }
  
  if(any(sapply(target, is.null)) ){
    avail.models <- sapply(models, function(x) x$method)
    if (!("selection.variable" %in% avail.models)) {
      stop("No way of setting target(s). Must either provide a named list or one of the methods should be the 'selection.variable' method since it preserves the original posterior.")
    } else {
      get_idx <- which(avail.models == "selection.variable")
    }
    selMod <- models[[get_idx]]
    last <- length(selMod$nzero)
    if (is.null(target$posterior)  & "posterior" %in% quantity) {
      post_targ <- selMod$theta[[last]]
    }
    if (is.null(target$mean) & "mean" %in% quantity) {
      mean_targ <- transform(selMod$eta[[last]])
    }
  }
  
  ### mse checks ###
  if (any(method == "mse") ) {
    mse.sel <- which(method == "mse")
    method[mse.sel] <- "exact"
    isMSE[mse.sel] <- TRUE
  }
  
  ### observation direction checks ###
  # obs.direction <- ifelse(grepl("univariate", method),"rowwise","colwise")
  obs.direction <- rep("colwise", length(method))
  
  
  ### model checks ###
  if (inherits(models, "sparse-posterior")) {
    models <- list(models)
  } else {
    if(!is.list(models)) stop("models must be a SparsePosterior fit or a list of fits")
  }
  
  ### set group names ###
  group_names <- names(models)
  if (is.null(group_names)) group_names <- seq.int(length(models))
  
  ### get dimensions ###
  n <- dim(models[[1]]$eta[[1]])[1]
  d <- dim(models[[1]]$theta[[1]])[1]
  s <- dim(models[[1]]$theta[[1]])[2]
  
  ### set posterior data if present ###
  if ("posterior" %in% quantity) {
    theta <- lapply(models, function(mm) mm$theta)
    if(!is.matrix(post_targ)) post_targ <- as.matrix(post_targ)
    if(nrow(post_targ) != d) post_targ <- t(post_targ)
    if(nrow(post_targ) != d) stop("Number of parameters in posterior target isn't equal to the number of parameters in theta")
    
    if(ncol(post_targ) == 1 & method[1] != "exact") {
      method[1] <- "exact"
      obs.direction[1] <- "colwise"
      warning("posterior target only has 1 observation. Changing to exact method which will be fast in this case.")
    }
    
    posterior <- list(source = theta,
                      target = post_targ,
                      method = method[1],
                      obs.direction = obs.direction[1],
                      isMSE = isMSE[1],
                      projection = rep(FALSE, length(theta))
    )
  } else {
    posterior <- NULL
  }
  
  ### set mean data if present ###
  if ("mean" %in% quantity) {
    eta <- lapply(models, function(mm) lapply(mm$eta, transform))
    
    if(!is.matrix(mean_targ)) mean_targ <- as.matrix(mean_targ)
    if(any(dim(mean_targ) %in% c(n,s))) {
      if(nrow(mean_targ) == s) mean_targ <- t(mean_targ)
    }
    if( nrow (mean_targ) != n) stop("Number of obsersvations of mean target not equal to number of observations in sample")
    if(ncol(mean_targ) == 1 & method[2] != "exact") {
      method[2] <- "exact"
      obs.direction[2] <- "colwise"
      warning("Mean target only has 1 observation. Changing to exact method which will be fast in this case")
    }
    
    mean <- list( source = eta,
                  target = mean_targ,
                  method = method[2],
                  obs.direction = obs.direction[2],
                  isMSE = isMSE[2],
                  projection = ifelse(sapply(models, function(mm) mm$method) == "projection",TRUE, FALSE)
    )
  } else {
    mean <- NULL
  }
  
  ### set output matrix ###
  output <- list( posterior = posterior,
                  mean = mean,
                  models = models,
                  quantity = quantity,
                  n = n, d = d, s = s,
                  group_names = group_names
                )
  
  return(output)
}

set_equal_y_limits.distcompare <- function(distance_data){
  dist <- ylim <- list(posterior = NULL, mean = NULL)
  for(i in c("posterior", "mean")){
    dist[[i]] <- list(dist = unlist(sapply(distance_data, function(x) x[[i]]$dist)))
    ylim[[i]] <- set_y_limits(dist, ylim[[i]], i)
  }
  return(ylim)
}

rank.distCompare <- function(distances) {
  if(!is.distcompare(distances)) stop("Must be distcompare object")
  rank.fun <- function(distance, quant) {
    dist <- distance[[quant]]$dist
    idx <- distance[[quant]]$nactive
    names(dist) <- distance[[quant]]$groups
    ranks <- tapply(dist, idx, rank, ties.method = "min")
    nactive <- unlist(mapply(FUN = function(n,r) {return(rep(n, r))}, 
                             n = names(ranks), r = sapply(ranks, length)))
    groups <- unlist(sapply(ranks, names))
    if(is.null(names(groups)) & length(levels(distance[[quant]]$groups)) == 1) groups <- levels(distance[[quant]]$groups)
    return( data.frame(ranks = unlist(ranks), nactive=nactive, groups = groups) )
  }
  rp.df <- rm.df <- NULL
  
  if ( !is.null(distances$posterior) ) {
    rp.df <- rank.fun(distances, "posterior")
  }
  
  if ( !is.null(distances$mean) ) {
    rm.df <- rank.fun(distances, "mean")
  }
  
  return(list(posterior = rp.df, mean = rm.df))
}