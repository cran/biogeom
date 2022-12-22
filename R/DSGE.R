DSGE <- function(P, phi){

  p <- length(P)
  if(p != 3) 
    stop("The number of parameters should equal 3!")

  a  <- P[1]
  n1 <- P[2]
  n2 <- P[3]
  a/4*n2/n1*(cos(phi/4)^n2 + sin(phi/4)^n2)^(-1/n1-1)*
  ( (cos(phi/4))^(n2-1)*sin(phi/4) - (sin(phi/4))^(n2-1)*cos(phi/4) )
}

