SHE <- function(P, x){
    delta <- P[1]
    rho   <- P[2]
    omega <- P[3]
    P     <- P[4]
    myfun <- function(z){
      z[z < delta] <- delta
      return(z)
    }
    x <- myfun(x)
    (1-rho) * ( (2/(P+1)) * ((x-delta)/(1-delta)) ) + rho*(
    (1-omega)*((x-delta)/(1-delta))^P + omega * (1-(1-((x-delta)/(1-delta)))^(1/P)) )
}
