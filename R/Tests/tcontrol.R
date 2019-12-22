trcontrol <- function(method = c("revsimplex", "shortsimplex", "primaldual", "aha", "shielding", "auction", "auctionbf"),
                      para=list(), start = c("auto", "modrowmin", "nwcorner", "russell"), nscales = 1, scmult = 2, returncoarse = FALSE,
                      a=NULL, b=NULL, M=NULL, N=NULL) {
  # a,b,M,N serve for computing parameters or start solutions automatically
  # by default M=a$N, N=b$N, which are overridden if M and/or N are specified
  method <- match.arg(method)
  start <- match.arg(start)
  if (is.null(M) && !is.null(a)) {
    M <- ifelse(class(a) %in% c("pgrid","pp","wpp"), a$N, length(a))
    # length(a) important if called from transport.default
  }
  if (is.null(N) && !is.null(b)) {
    N <- ifelse(class(b) %in% c("pgrid","pp","wpp"), b$N, length(b))
    # length(b) important if called from transport.default
  }
  if (is.null(M) && !is.null(N)) {M <- N}
  if (is.null(N) && !is.null(M)) {N <- M}
  if (is.null(N) && method %in% c("shortsimplex","auction")) {
    stop("At least M or N must be specified for method ", method)
  }
  
  if (method == "shortsimplex") {
    
    if (is.numeric(para) && length(para) == 3) {
      para = list(slength=para[1], kfound=para[2], psearched=para[3])
      message("Parameter vector interpreted by trcontrol as: \n slength = ",para[1],"; kfound = ",para[2],"; psearched = ",para[3],".")
    }
    
    propnames <- c("slength","kfound","psearched")
    done <- pmatch(names(para), propnames)
    if (!is.list(para) || (length(para) > 0 && length(done) == 0) || any(is.na(done))) {
      stop("expected for 'para' with method='shortlist' either an empty list or a named list with 1-3 components out of 'slength', 'cfound', or 'psearched' (after partial matching)")
    }
    
    newpara <- list(slength=0, kfound=0, psearched=0)
    
    if (1 %in% done) {
      newpara$slength <- para[[which(done == 1)]]
      #      if (newpara$slength > b$N) {
      #      	warning("parameter 'slength' for shortlist algorithm was larger then no. of target points... Fixed.")
      #      	newpara$slength <- b$N
      if (newpara$slength < 1) {
        stop("parameter 'slength' for shortlist algorithm has to be >= 1")
      }
    } else {
      newpara$slength <- min(N,15 + max(0, floor(15 * log(N/400)/log(2)))) 
      
      # /400 since after reduction about half of the N producers disappear
      #newpara$slength <- min(b$N, newpara$slength)
      # Can't choose the shortlist longer then the number of consumers
      # This is caught now in transport functions
    }
    
    if (2 %in% done) {
      newpara$kfound <- para[[which(done == 2)]]
      if (newpara$kfound < 1) {
        stop("parameter 'kfound' for shortlist algorithm has to be >= 1")
      }
    } else {
      newpara$kfound <- newpara$slength
    }
    
    if (3 %in% done) {
      newpara$psearched <- para[[which(done == 3)]]
      if (newpara$psearched > 1 || newpara$psearched <= 0) {
        stop("parameter 'psearched' for shortlist algorithm has to be > 0 and <= 1")
      }
    } else {
      newpara$psearched <- 0.05
    }
    # Note that at least one row is searched anyway (even if psearch was 0 or negative)
    
    para <- newpara
  }
  # fi (method == "shortsimplex")
  
  
  res <- list(method=method, para=para, start=start, nscales=nscales, scmult=scmult, returncoarse=returncoarse)
  class(res) <- "trc"
  
  return(res)
}

