ecos.control.dots <- function(maxit = 100L, feastol = 1e-08, reltol = 1e-08, abstol = 1e-08, 
                              feastol_inacc = 1e-04, abstol_inacc = 5e-05, reltol_inacc = 5e-05, 
                              verbose = 0L, mi_max_iters = 1000L, mi_int_tol = 1e-04, mi_abs_eps = 1e-06, 
                              mi_rel_eps = 1e-06, ...) {
  ECOSolveR::ecos.control(maxit = maxit, feastol = feastol, reltol = reltol, abstol = abstol,
                          feastol_inacc = feastol_inacc, abstol_inacc = abstol_inacc, reltol_inacc,
                          verbose = verbose, mi_max_iters = mi_max_iters, mi_int_tol, mi_abs_eps, mi_rel_eps)
}

ecos.control.better <-  function(control = NULL) {
  pot.names <- names(formals(ECOSolveR::ecos.control))
  list.names <- pot.names[pmatch(names(control),  pot.names)]
  if(is.na(list.names)) list.names <- NULL
  names(control) <- list.names
  return(do.call(ecos.control.better, control))
}

