mround <- function(x, base) {
  base * round(x/base)
}

cround <- function(x, base) {
  base * ceiling(x/base)
}

fround <- function(x, base) {
  base * floor(x/base)
}

getDigits <- function(x) {
  if (abs(x) < 1) {
    digi <- nchar(strsplit(sub('0+$', '', as.character(x)), ".", fixed = TRUE)[[1]][[2]])
    quant <- x + 5/(10 ^(digi + 1))
  } else {
    digi <- (-(nchar(as.character(round(x)))-1))
    quant <- x + 5*10^abs(digi)
  }
  return(mround(quant, digi))
}