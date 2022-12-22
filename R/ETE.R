ETE <- function(P, x){
  if(!is.numeric(x))
    stop("'x' should be a numerical value or a numerical vector!")
  a      <- P[1]
  alpha0 <- P[2]
  alpha1 <- P[3]
  alpha2 <- P[4] 
  x[x < -a] <- -a
  x[x >  a] <- a
  a * exp( alpha0 + alpha1*(x/a) + alpha2*(x/a)^2 ) * sqrt(1 - (x/a)^2)
}
