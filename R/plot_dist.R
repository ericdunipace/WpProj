plot.distcompare <- function(distance = NULL, models = NULL, ylim = NULL, ylabs = c(NULL,NULL),
                             xlab = NULL, xlim = NULL,
                             linesize = 1, pointsize = 1, ...) {
  
  mc <- match.call(expand.dots = TRUE)
  if(is.null(distance)) {
    distance <- distCompare(models, ...)
  }
  if (!inherits(distance, "distcompare")) stop("distance must be output of distCompare function")
  
  ppost <- pmean <- NULL
  
  if(is.null(xlab)) xlab <- "Number of active coefficients"
  
  if ( !is.null(distance$posterior) ) {
    ylim_post <- set_y_limits(distance, ylim, "posterior")
    xlim_post <- set_x_limits(distance, xlim, "posterior")
    
    ppost <- ggplot2::ggplot( distance$posterior, 
                              ggplot2::aes(x=nactive, y=dist, color = groups, group=groups )) +
      ggplot2::geom_line(size = linesize, position = ggplot2::position_dodge(width = 0.25)) + 
      ggplot2::geom_point(size = pointsize, position = ggplot2::position_dodge(width = 0.25)) + 
      ggsci::scale_color_jama() + 
      ggplot2::labs(color ="Method") +
      ggplot2::xlab(xlab) + 
      ggplot2::ylab(ylabs[1]) + ggplot2::theme_bw() +
      ggplot2::scale_x_continuous(expand = c(0, 0), limits = xlim_post) +
      ggplot2::scale_y_continuous(expand = c(0, 0), limits = ylim_post )
  }
  
  if (!is.null(distance$mean)){
    ylim_mean <- set_y_limits(distance, ylim, "mean")
    xlim_mean <- set_x_limits(distance, xlim, "mean")
    pmean <- ggplot2::ggplot( distance$mean, 
                              ggplot2::aes(x=nactive, y=dist, color = groups, group=groups )) +
      ggplot2::geom_line(size = linesize, position = ggplot2::position_dodge(width = 0.25)) + 
      ggplot2::geom_point(size = pointsize, position = ggplot2::position_dodge(width = 0.25)) +
      ggsci::scale_color_jama() + 
      ggplot2::labs(color ="Method") +
      ggplot2::xlab(xlab) + 
      ggplot2::ylab(ylabs[length(ylabs)]) + ggplot2::theme_bw() +
      ggplot2::scale_x_continuous(expand = c(0, 0), limits = xlim_mean) +
      ggplot2::scale_y_continuous(expand = c(0, 0), limits = ylim_mean )
  }
  
  plots <- list(posterior = ppost, mean = pmean)
  class(plots) <- c("plotcompare","sparse-posterior")
  return(plots)
}

print.plotcompare <- function(x) {
  for(i in 1:length(x)) {
    if(is.null(x[[i]])) next
    print(x[[i]])
  }
}

set_y_limits <- function(distance_data, ylim, quantity){
  idx <- if (quantity == "posterior"){
    1L
  } else if (quantity == "mean") {
    2L
  }
  
  if (!is.null(ylim)) {
    if (is.numeric(ylim)){
      if (length(ylim) == 4){
        return(ylim[(idx-1)*2 + 1:2])
      } else {
        return(ylim)
      }
    } 
    if (is.list(ylim) & !is.null(ylim[[idx]])) return(ylim[[idx]])
  }
  
  df <- distance_data[[idx]]
  if (is.null(df)) return(NULL)
  if (is.null(df$dist)) return(NULL)
  range.size <- diff(range(df$dist))
  add.factor <- range.size * 1.2 - range.size
  min_y <- max(0, min(df$dist) - add.factor)
  max_y <- max(df$dist) + add.factor
  ylim <- c(min_y, max_y)
  return(ylim)
}


set_x_limits <- function(distance_data, xlim, quantity){
  idx <- if (quantity == "posterior"){
    1L
  } else if (quantity == "mean") {
    2L
  }
  
  if (!is.null(xlim)) {
    if (is.numeric(xlim)){
      if (length(xlim) == 4){
        return(xlim[(idx-1)*2 + 1:2])
      } else {
        return(xlim)
      }
    } 
    if (is.list(xlim) & !is.null(xlim[[idx]])) return(xlim[[idx]])
  }
  
  df <- distance_data[[idx]]
  if (is.null(df)) return(NULL)
  if (is.null(df$nzero)) return(NULL)
  min_x <- min(df$nzero)
  max_x <- max(df$nzero)
  xlim <- c(min_y, max_y)
  return(xlim)
}

set_equal_y_limits.plotcompare <- function(distance_data){
  # dist.list <- list(dist = unlist(sapply(distance_data, function(x) x[[quantity]]$data$dist )))
  dist <- ylim <- list(posterior = NULL, mean = NULL)
  for(i in c("posterior", "mean")){
    dist[[i]] <- list(dist = unlist(sapply(distance_data, function(x) x[[i]]$data$dist )))
    ylim[[i]] <- set_y_limits(dist, ylim[[i]], i)
  }
  for(j in seq_along(distance_data)) {
    for(i in c("posterior", "mean")) {
      distance_data[[j]][[i]] <- distance_data[[j]][[i]] + ggplot2::scale_y_continuous(expand = c(0, 0), limits = ylim[[i]] )
    }
  }
  return(distance_data)
}