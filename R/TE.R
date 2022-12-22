TE <- function(P, x){
  if(!is.numeric(x))
    stop("'x' should be a numerical value or a numerical vector!")
  x[x < -1] <- -1
  x[x >  1] <- 1
  alpha0 <- P[1]
  alpha1 <- P[2]
  alpha2 <- P[3] 
  exp( alpha0 + alpha1*x + alpha2*x^2 ) * sqrt(1 - x^2)
}
