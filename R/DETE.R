DETE <- function(P, x){

  if(P[1] <= 0)
    stop("a should be positive real numbers!")
  p <- length(P)
  if(p != 4) 
    stop("The number of parameters should equal 4!")

  x[x < -P[1]] <- -P[1]
  x[x >  P[1]] <- P[1]
  a      <- P[1]
  alpha0 <- P[2]
  alpha1 <- P[3]
  alpha2 <- P[4]
  (alpha1+2*alpha2/a*x-x/a*(1-(x/a)^2)^(-1)) * ETE(P=P, x=x)/a

}



