combine.dist.compare <- function(distances) {
  stopifnot(is.list(distances))
  if (!all(sapply(distances, is.distcompare))) {
    stop("All members must be  distcompare object")
  }
  niter <- length(distances)
  length.each <- sapply(distances, function(i) nrow(i$mean))
  cmb <- list(posterior = NULL, mean = NULL, p = NULL)
  
  ps <- sapply(distances, function(d) d$p)
  stopifnot(all(diff(ps)==0))
  
  ranks.list <- lapply(distances, rank.distCompare)
  ranks.post <- unlist(sapply(ranks.list, function(r) r$posterior$ranks))
  ranks.mean <- unlist(sapply(ranks.list, function(r) r$mean$ranks))
  iter <- unlist(sapply(1:niter, function(i) rep(i, length.each[i])))
  
  cmb$p <- p[1]
  post <- do.call("rbind", lapply(distances, function(d) d$posterior))
  means <- do.call("rbind", lapply(distances, function(d) d$mean))
  # methods <- do.call("rbind", lapply(distances, function(d) d$method))
  
  if(! is.null(post)) {
    cmb$posterior <- post
    cmb$posterior <- cbind(cmb$posterior, ranks = ranks.post, iter = iter)
  }
  
  if (!is.null(means)) {
    cmb$mean <- means
    if(!is.null(cmb$mean)) cmb$mean <- cbind(cmb$mean, ranks = ranks.mean, iter = iter)
  }
  
  class(cmb) <- c("limbs","combine.dist.compare")
  
  return(cmb)
}


plot.combine.dist.compare <- function (distances, ylim = NULL, ylabs = c(NULL,NULL), facet.group = NULL, ...) {
  stopifnot(inherits(distances, "combine.dist.compare"))
  dots <- list(...)
  alpha <- dots$alpha
  base_size <- dots$base_size
  ribbon <- dots$ribbon
  xlab <- dots$xlab
  leg.pos <- dots$legend.position
  if(is.null(ribbon)) ribbon <- FALSE
  if(is.null(alpha)) alpha <- 0.3
  if(is.null(base_size)) base_size <- 11
  if(is.null(xlab)) xlab <- "Number of active coefficients"
  if(is.null(leg.pos)) leg.pos <- NULL
  
  # methods <- levels(distances$mean$groups)
  d <- max(distances$mean$nactive)
  numactive <- 1:d
  
  # df <- data.frame(numactive = numactive, method = rep(methods, each = d))
  
  ppost <- pmean <- NULL
  
  if ( !is.null(distances$posterior) ) {
    dd <- distances$posterior
    
    dd$groups <- factor(dd$groups)
    if(!is.null(facet.group)) {
      grouping <- c("groups", facet.group, "nactive")
    } else {
      grouping <- c("groups", "nactive")
    }
    df <- dd %>% dplyr::group_by(.dots = grouping) %>% dplyr::summarise(
                                                    low = quantile(dist, 0.025),
                                                    hi = quantile(dist, 0.975),
                                                    dist = mean(dist)
                                                    )
    # E <- tapply(dd$dist, INDEX = grouping, mean)
    # sigma <- tapply(dd$dist, INDEX = grouping, sd)
    # 
    # M <- tapply(dd$dist, INDEX = grouping, median)
    # hi <- tapply(dd$dist, INDEX = grouping, quantile, 0.975)
    # low <- tapply(dd$dist, INDEX = grouping, quantile, 0.025)
    
    # df <- data.frame(dist = c(M), low = c(low), hi = c(hi), 
    #                  groups = rep(colnames(M), each = nrow(M)),
    #                  nactive = as.numeric(rep(rownames(M), ncol(M))),
    #                  row.names = NULL)
    # df <- df[complete.cases(df),]
    
    ylim_post <- set_y_limits(df, ylim, "posterior")
    ppost <- ggplot2::ggplot( df, 
                              ggplot2::aes(x=nactive, y=dist, 
                                           color = groups, fill = groups,
                                           group=groups ))
    if(ribbon) {
      ppost <- ppost + ggplot2::geom_ribbon(ggplot2::aes(ymin = low, ymax = hi), alpha = alpha, linetype=0)
    } else {
      ppost <- ppost + ggplot2::geom_errorbar(ggplot2::aes(ymin = low, ymax = hi), alpha = alpha, position = ggplot2::position_dodge(width=0.25))
    }
    ppost <- ppost + ggplot2::geom_line(position = ggplot2::position_dodge(width=0.25)) +
      ggplot2::geom_point(position = ggplot2::position_dodge(width=0.25)) +
      ggsci::scale_color_jama() + 
      ggsci::scale_fill_jama() +
      ggplot2::labs(fill ="Method", color="Method") +
      ggplot2::xlab(xlab) + 
      ggplot2::ylab(ylabs[1]) + ggplot2::theme_bw(base_size) +
      ggplot2::scale_x_continuous(expand = c(0, 0)) +
      ggplot2::scale_y_continuous(expand = c(0, 0), limits = ylim_post ) + 
      ggplot2::theme(legend.position = leg.pos)
    if(!is.null(facet.group)) {
      ppost <- ppost + ggplot2::facet_grid(stats::reformulate(facet.group))
    }
  }
  
  if (!is.null(distances$mean)){
    dd <- distances$mean
    dd$groups <- factor(dd$groups)
    if(!is.null(facet.group)) {
      grouping <- c("groups", facet.group, "nactive")
    } else {
      grouping <- c("groups", "nactive")
    }
    df <- dd %>% dplyr::group_by(.dots = grouping) %>% dplyr::summarise(
      low = quantile(dist, 0.025),
      hi = quantile(dist, 0.975),
      dist = mean(dist)
    )
    # E <- tapply(dd$dist, INDEX = list(dd$nactive, dd$groups), mean)
    # sigma <- tapply(dd$dist, INDEX = list(dd$nactive, dd$groups), sd)
    # 
    # M <- tapply(dd$dist, INDEX = list(dd$nactive, dd$groups), median)
    # hi <- tapply(dd$dist, INDEX = list(dd$nactive, dd$groups), quantile, 0.975)
    # low <- tapply(dd$dist, INDEX = list(dd$nactive, dd$groups), quantile, 0.025)
    # 
    # df <- data.frame(dist = c(M), low = c(low), hi = c(hi), 
    #                  groups = rep(colnames(M), each = nrow(M)),
    #                  nactive = as.numeric(rep(rownames(M), ncol(M))),
    #                  row.names = NULL)
    
    ylim_mean <- set_y_limits(df, ylim, "mean")
    
    pmean <- ggplot2::ggplot( df, 
                              ggplot2::aes(x=nactive, y=dist, 
                                           color = groups, fill = groups,
                                           group=groups ))
    if(ribbon) {
      pmean <- pmean + ggplot2::geom_ribbon(ggplot2::aes(ymin = low, ymax = hi), alpha = alpha, linetype=0)
    } else {
      pmean <- pmean + ggplot2::geom_errorbar(ggplot2::aes(ymin = low, ymax = hi), alpha = alpha, position = ggplot2::position_dodge(width=0.25))
    }
    pmean <- pmean + ggplot2::geom_line(position = ggplot2::position_dodge(width=0.25)) +
      ggplot2::geom_point(position = ggplot2::position_dodge(width=0.25)) +
      ggsci::scale_color_jama() + 
      ggsci::scale_fill_jama() +
      ggplot2::labs(fill ="Method", color="Method") +
      ggplot2::xlab(xlab) + 
      ggplot2::ylab(ylabs[1]) + ggplot2::theme_bw(base_size) +
      ggplot2::scale_x_continuous(expand = c(0, 0)) +
      ggplot2::scale_y_continuous(expand = c(0, 0), limits = ylim_mean ) + 
      ggplot2::theme(legend.position = leg.pos)
    if(!is.null(facet.group)) {
      pmean <- pmean + ggplot2::facet_grid(stats::reformulate(facet.group))
    }
  }
  
  plots <- list(posterior = ppost, mean = pmean)
  class(plots) <- c("plotcombine","limbs")
  return(plots)
}

print.plotcombine <- function(x) {
  for(i in 1:length(x)) {
    if(is.null(x[[i]])) next
    print(x[[i]])
  }
}

plot_ranks <- function(distances, ylim = NULL, ylabs = c(NULL,NULL), ...) {
  stopifnot(inherits(distances, "combine.dist.compare"))
  dots <- list(...)
  alpha <- dots$alpha
  if(is.null(alpha)) alpha <- 0.1
  ppost <- pmean <- NULL
  
  countfun <- function(x) {
    tab <- table(x)
    n <- sum(tab)
    return(tab/n)
  }
  
  if ( !is.null(distances$posterior) ) {
    dd <- distances$posterior
    index <- list(dd$nactive, dd$groups)
    
    M <- tapply(dd$ranks, INDEX = index, mean)
    low <- tapply(dd$ranks, INDEX = index, quantile, 0.025)
    hi <- tapply(dd$ranks, INDEX = index, quantile, 0.975)
    
    df <- data.frame(dist = c(M), low = c(low), hi = c(hi), 
                     groups = rep(colnames(M), each = nrow(M)),
                     nactive = as.numeric(rep(rownames(M), ncol(M))),
                     row.names = NULL)
    
    ylim <- set_y_limits(distances, ylim, "posterior")
    ppost <- ggplot2::ggplot( df, 
                              ggplot2::aes(x=nactive, y=dist, 
                                           color = groups, fill = groups,
                                           group=groups )) +
      ggplot2::geom_line() + 
      ggplot2::geom_ribbon(ggplot2::aes(ymin = low, ymax = hi), alpha = alpha, linetype=0) + 
      ggsci::scale_color_jama() + 
      ggsci::scale_fill_jama() +
      ggplot2::labs(fill ="Method", color="Method") +
      ggplot2::xlab(xlab) + 
      ggplot2::ylab(ylabs[1]) + ggplot2::theme_bw() +
      ggplot2::scale_x_continuous(expand = c(0, 0)) +
      ggplot2::scale_y_continuous(expand = c(0, 0), limits = ylim )
  }
  
  if (!is.null(distances$mean)){
    dd <- distances$mean
    index <- list(dd$nactive, dd$groups)
    
    M <- tapply(dd$ranks, INDEX = index, mean)
    low <- tapply(dd$ranks, INDEX = index, quantile, 0.025)
    hi <- tapply(dd$ranks, INDEX = index, quantile, 0.975)
    
    df <- data.frame(dist = c(M), low = c(low), hi = c(hi), 
                     groups = rep(colnames(M), each = nrow(M)),
                     nactive = as.numeric(rep(rownames(M), ncol(M))),
                     row.names = NULL)
    
    ylim <- set_y_limits(distances, ylim, "mean")
    pmean <- ggplot2::ggplot( df, 
                              ggplot2::aes(x=nactive, y=dist, 
                                           color = groups, fill = groups,
                                           group=groups )) +
      ggplot2::geom_line() + 
      ggplot2::geom_ribbon(ggplot2::aes(ymin = low, ymax = hi), alpha = alpha, linetype=0) + 
      ggsci::scale_color_jama() + 
      ggsci::scale_fill_jama() +
      ggplot2::labs(fill ="Method", color="Method") +
      ggplot2::xlab(xlab) + 
      ggplot2::ylab(ylabs[1]) + ggplot2::theme_bw() +
      ggplot2::scale_x_continuous(expand = c(0, 0)) +
      ggplot2::scale_y_continuous(expand = c(0, 0), limits = ylim )
  }
  
  plots <- list(posterior = ppost, mean = pmean)
  class(plots) <- c("plotrank","limbs")
  return(plots)
}

print.plotrank <- function(x) {
  for(i in 1:length(x)) {
    if(is.null(x[[i]])) next
    print(x[[i]])
  }
}
