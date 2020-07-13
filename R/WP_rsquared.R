# WPR2 <- function(Y, nu, p = 2, method = "exact",...) {UseMethod("WPR2")}
setClass("WPR2",
         slots = c(r2 = "numeric", 
              nactive = "integer", 
              groups = "factor",
              method = "factor",
              p = "numeric",
              base = "factor"),
         contains = "data.frame")

WPR2.matrix <- function(Y, nu, p = 2, method = "exact", ...) {
  
  stopifnot(p >= 1)
  
  n <- nrow(Y)
  d <- ncol(Y)
  
  wp_mod <- limbs::wasserstein(Y, nu, p = p, ground_p = p,
                           method = method, ...)^p
  mu <- matrix(colMeans(Y), n, d, byrow=TRUE)
  # wp_base <- if(method == "exact") {
  #   mean(colSums((Y - mu)^p))
  # } else {
  #   limbs::wasserstein(Y, mu, p = p, 
  #                      ground_p = p,
  #                      method = method, 
  #                      ...)^p
  # }
  wp_base <- limbs::wasserstein(Y, mu, p = p, 
                                ground_p = p,
                                method = method, 
                                ...)^p
  
  r2 <- pmax(1 - wp_mod/wp_base,0)
  output <- data.frame(r2 = r2, method = method)
  class(output) <- c("WPR2", class(output))
  return(output)
  
}

WPR2.distcompare <- function(Y=NULL, nu, ...) {
  
  stopifnot(inherits(nu, "distcompare"))
  
  df <- nu$mean
  p <- nu$p
  method <- df$method
  
  stopifnot(p >= 1)
  
  
  if(!is.null(Y)) {
    stopifnot(inherits(Y, "matrix"))
    meth.table <- table(method)
    method.use <- names(meth.table)[which.max(meth.table)]
    wass.args <- list(X = Y, Y = as.matrix(rowMeans(Y)),
                      p = as.numeric(p), method = method.use,
                      ...)
    wass.args <- wass.args[!duplicated(names(wass.args))]
    argn <- lapply(names(wass.args), as.name)
    names(argn) <- names(wass.args)
    
    wass.call <- as.call(c(list(as.name("wasserstein")), argn))
    
    max_vals <- eval(wass.call, envir = wass.args)
    max_vec <- rep(max_vals, length(df$dist))
    base <- "dist.from.expectation"
  } else {
    max_vals <- tapply(df$dist, df$groups, max)
    max_vec <- max_vals[as.numeric(df$groups)]
    base <- "dist.from.null"
  }
  
  
  
  r2 <- pmax(1- df$dist^p/max_vec^p, 0)
  
  df$dist <- r2
  df$p <- p
  df$base <- base
  colnames(df)[1] <- "r2"
  
  class(df) <- c("WPR2", class(df))
  return(df)
  
}

setGeneric("WPR2", function(Y = NULL, nu, p = 2, method = "exact",...) UseMethod("WPR2"))
setMethod("WPR2", c("Y" = "matrix", "nu" = "matrix"), WPR2.matrix)
setMethod("WPR2", c("nu" = "distcompare"), WPR2.distcompare)


combine.WPR2 <- function(...) {
  if(...length()>1){
    objects <- list(...)
  } else {
    objects <- c(...)
  }
  stopifnot(is.list(objects))
  if (!all(sapply(objects, inherits, "WPR2"))) {
    stop("All objects must be WPR2 object")
  }
  niter <- length(objects)
  length.each <- sapply(objects, nrow)
  
  ps <- unlist(sapply(objects, function(d) d$p))
  if(!all(diff(ps)==0)) {
    stop("Some of the wasserstein powers in the WPR2 objects are different")
  }
  
  base <- unlist(sapply(objects, function(d) d$base))
  if((any(is.na(base)) & !all(is.na(base))) | !all(base == base[1])) {
    stop("Some of the objects are using different comparison points (i.e. null model vs expectation of full model)")
  }
  
  wpmeth <- unlist(sapply(objects, function(d) d$method))
  if(!all(wpmeth == wpmeth[1])) {
    warning("Some of the Wasserstein distances were calculated with different methods. 
            Advisable not to compare some of these objects")
  }
  
  cmb <- do.call("rbind", objects)
  
  # cmb$iter <- unlist(sapply(1:niter, function(i) rep(i, length.each[i])))
  
  # class(cmb) <- c("WPR2", class(cmb))
  
  return(cmb)
}

plot.WPR2 <- function(object, xlim = NULL, ylim = NULL, linesize =1, pointsize = 1, ...) {
  obj <- object
  stopifnot(inherits(obj, "WPR2"))
  dots <- list(...)
  alpha <- dots$alpha
  base_size <- dots$base_size
  ribbon <- dots$ribbon
  xlab <- dots$xlab
  ylab <- dots$ylab
  leg.pos <- dots$legend.position
  if(is.null(ribbon)) ribbon <- FALSE
  if(is.null(alpha)) alpha <- 0.3
  if(is.null(base_size)) base_size <- 11
  if(is.null(xlab)) xlab <- "Number of active coefficients"
  if(is.null(leg.pos)) leg.pos <- NULL
  
  wp_calc_method <- obj$method
  opt_method <- obj$groups
  p <- unique(obj$p)[1]
  base <- obj$base
  nactive <- obj$nactive
  
  
  ylim <- set_y_limits_gen(obj$r2, ylim)
  if(is.null(ylab)) {
    ylab <- bquote(W[.(p)]~R^2)
  }
  
  if(all(!is.na(nactive))) {
    if(is.null(xlab)) {
     xlab <- "Number of active coefficients" 
    }
    xlim <- set_x_limits_gen(obj$nactive, xlim)
    
    E <- tapply(obj$r2, INDEX = list(obj$nactive, obj$groups), mean)
    # sigma <- tapply(obj$dist, INDEX = list(obj$nactive, obj$groups), sd)
    
    M <- tapply(obj$r2, INDEX = list(obj$nactive, obj$groups), median)
    hi <- tapply(obj$r2, INDEX = list(obj$nactive, obj$groups), quantile, 0.975)
    low <- tapply(obj$r2, INDEX = list(obj$nactive, obj$groups), quantile, 0.025)
    
    df <- data.frame(dist = c(M), low = c(low), hi = c(hi), 
                     groups = rep(colnames(M), each = nrow(M)),
                     nactive = as.numeric(rep(rownames(M), ncol(M))),
                     row.names = NULL)
    plot <- ggplot2::ggplot( df, 
                     ggplot2::aes(x=nactive, y=dist, 
                                  color = groups, fill = groups,
                                  group=groups ))
    if(ribbon) {
      plot <- plot + ggplot2::geom_ribbon(ggplot2::aes(ymin = low, ymax = hi), alpha = alpha, linetype=0)
    } else {
      plot <- plot + ggplot2::geom_errorbar(ggplot2::aes(ymin = low, ymax = hi), alpha = alpha, position = ggplot2::position_dodge(width=0.25))
    }
    plot <- plot + ggplot2::geom_line(position = ggplot2::position_dodge(width=0.25), size = linesize) +
      ggplot2::geom_point(position = ggplot2::position_dodge(width=0.25), size = pointsize) +
      ggsci::scale_color_jama() + 
      ggsci::scale_fill_jama() +
      ggplot2::labs(fill ="Method", color="Method") +
      ggplot2::xlab(xlab) + 
      ggplot2::ylab(ylab) + ggplot2::theme_bw(base_size) +
      ggplot2::scale_x_continuous(expand = c(0, 0), limits = xlim) +
      ggplot2::scale_y_continuous(expand = c(0, 0), limits = ylim ) + 
      ggplot2::theme(legend.position = leg.pos)
  } else {
    plot <- ggplot2::ggplot(data = obj, mapping = ggplot2::aes(y = r2, fill = groups)) +
      geom_bar() + ggsci::scale_fill_jama() +
      ggplot2::labs(fill ="Method", color="Method") +
      ggplot2::xlab(xlab) + 
      ggplot2::ylab(ylab) + ggplot2::theme_bw(base_size) +
      ggplot2::scale_x_continuous(expand = c(0, 0), limits = xlim) +
      ggplot2::scale_y_continuous(expand = c(0, 0), limits = ylim ) + 
      ggplot2::theme(legend.position = leg.pos)
    if(nrow(obj) == 1) {
      plot <- plot + ggplot2::theme(legend.position = "none")
    }
  }

  return(plot)
}

set_y_limits_gen <- function(vals, ylim){
  if (!is.null(ylim)) {
    if (is.numeric(ylim)) {
      if(length(ylim) == 2) {
        return(ylim)
      } else {
        stop("if provided, ylim must be of length two")
      }
    } else {
      stop("ylim must be a numeric vector")
    }
  }
  if(is.null(vals)) return(c(0,1))
  # range.size <- diff(range(vals))
  # add.factor <- range.size * 1.2 - range.size
  # min_y <- max(0, min(vals) - add.factor)
  # max_y <- max(vals) + add.factor
  # ylim <- c(min_y, max_y)
  ylim <- c(0,1)
  return(ylim)
}


set_x_limits_gen <- function(vals, xlim){
  
  if (!is.null(xlim)) {
    if (is.numeric(xlim)) {
        if(length(xlim) == 2) {
          return(xlim)
        } else {
          stop("if provided, xlim must be of length two")
        }
    } else {
      stop("xlim must be a numeric vector")
    }
  }
  
  if (is.null(vals)) return(NULL)
  min_x <- min(vals)
  max_x <- max(vals)
  xlim <- c(min_x, max_x)
  return(xlim)
}
