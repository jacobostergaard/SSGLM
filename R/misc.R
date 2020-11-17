add.alpha <- function(col, alpha=1){
  ## Add an alpha value to a color to make it transparent
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2,
        function(x)
          rgb(x[1], x[2], x[3], alpha=alpha))
}

ks_crit <- function(a){
  # Get critical values of Kolmogorov-Smirnov test
  # see: https://en.wikipedia.org/wiki/Kolmogorov–Smirnov_test or
  # "The Art of Computer Programming, Volume 2" by DE Knuth, Eq (15) in section 3.3.1

  sqrt(-.5*log(a/2))
}




