DNRGE <- function(P, x){

  if(P[1] <= 0)
    stop("A should be positive real numbers!")
  p <- length(P)
  if(p != 4) 
    stop("The number of parameters should equal 4!")

  x[x < -P[1]/2] <- -P[1]/2
  x[x >  P[1]/2] <- P[1]/2
  A  <- P[1]
  B  <- P[2]
  C  <- P[3]
  D  <- P[4]
  f1 <- B/2*sqrt((A^2-4*x^2)/(A^2+8*C*x+4*C^2))
  f2 <- sqrt((A*(A^2+8*C*x+4*C^2))/(2*(A-2*C)*
         x^2+(A^2+8*A*C-4*C^2)*x+2*A*C^2+A^2*C+A^3))
  f3 <- A^2 - 4*x^2
  f4 <- A^2 + 8*C*x + 4*C^2
  E  <- (sqrt(5.5*A^2+11*B*C+4*C^2)*
          (sqrt(3)*B*A-2*D*sqrt(A^2+2*C*A+4*C^2)))/
          (sqrt(3)*B*A*(sqrt(5.5*A^2+11*A*C+4*C^2)-2*
           sqrt(A^2+2*C*A+4*C^2)))
  F <- 2*(A-2*C)
  G <- A^2 + 8*A*C - 4*C^2
  H <- 2*A*C^2 + A^2*C + A^3

  -4*f1*(C*f3+x*f4)/(f3*f4)*(1-E*(1-f2))-A*E/2*f1/f2*(f4*(2*F*x+G)/(F*x^2+G*x+H)^2)
}



