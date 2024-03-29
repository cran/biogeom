SarabiaE <- function(P, x){
    lambda <- P[1]
    eta    <- P[2]
    a1     <- P[3]
    a2     <- P[4]
    (1-lambda+eta)*x + lambda*x^(a1+1) - eta*(1-(1-x)^(a2+1))
}